"""한 모드 식 시험 (docs/GATE_JACOBIAN_DAMPING.md): 근사 Jacobian 의 감쇠를 식 1 / (1 + h² λ) 로 맞힐 수 있는가.

사용:  python tests/gate_jacobian_damping.py <브랜치 cfmsDrape.dll>     (리포 루트에서, run_gate 밖 — 약 10 분)
       (내부) python tests/gate_jacobian_damping.py --child <인자 json>
A  행렬의 성질 — 엔진의 디버그 고리 debug_stiffness_product 로 좌변 행렬 ∂f/∂x 를 읽는다(두 팔: off = Oh, gn_state = GN).
B0 한 모드의 감쇠 · B1 강체 회전의 감쇠(무중력) · B2 중력 아래의 종단 속도 · C 판이 굽는가(보고).
예측은 A 의 행렬에서 계산해 B 를 돌리기 전에 찍는다. 결과 results/gate_jacobian_damping.json. 통과 0, 미달 1, 준비 안 됨 2.
"""
import sys
sys.dont_write_bytecode = True

import concurrent.futures  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
import uuid  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import link_l0 as L0  # noqa: E402

ROOT = L0.ROOT
OUT = os.path.join(ROOT, 'results', 'gate_jacobian_damping.json')
PLATE = 30.0                    # cm
GRAV = 980.0
ARMS = ('off', 'gn_state')      # Oh, GN
MASS0, DX0 = 0.1, 0.375         # 장면 P 의 노드 질량(g) · 간격(cm)
PCG_TOL, PCG_CAP = 1e-6, 2000
PARALLEL = max(1, min(16, (os.cpu_count() or 2) // 2))
BASE = (81, 1e4)
MATRIX_CASES = [(81, 1e3), (81, 1e4), (81, 1e5), (41, 1e4), (21, 1e4)]
B0_SPS, B0_MODES = (200, 50), (1, 2, 3)
B1_SPS, B1_RATE = 200, 0.05     # rad/s
B2_CASES = ([(81, b, s) for b in (1e3, 1e4, 1e5) for s in (800, 400, 200, 100, 50)]
            + [(n, 1e4, s) for n in (41, 21) for s in (200, 100, 50)])
B2_CLEAN = [(81, 1e4, s) for s in (200, 100, 50)]
B2_T = 1.0
TERMINAL = 20.0                 # h λ1 ≥ 20 / s
C_CASE = (81, 1e4, 800)


def scene(n):
    xs = np.linspace(-PLATE / 2, PLATE / 2, n)
    gx, gz = np.meshgrid(xs, xs)
    V = np.column_stack([gx.ravel(), np.full(gx.size, 20.0), gz.ravel()])
    F = []
    for j in range(n - 1):
        for i in range(n - 1):
            a = j * n + i
            F += [[a, a + 1, (j + 1) * n + i + 1], [a, (j + 1) * n + i + 1, (j + 1) * n + i]]
    uv = np.column_stack([gx.ravel(), gz.ravel()])
    s = gz.ravel() + PLATE / 2                       # 고정 변에서 잰 거리
    row = np.repeat(np.arange(n), n)                 # 줄 번호(0 = 고정 줄)
    return V, np.array(F, np.int32), uv, list(range(n)), s, row


def node_mass(n):
    dx = PLATE / (n - 1)
    return MASS0 * (dx / DX0) ** 2


def make(dll, n, bend, h, clean, gravity, extra=None):
    Drape = L0.load_drape()
    V, F, uv, pins, s, row = scene(n)
    d = Drape(dll)
    d.set_mesh(V, F)
    d.set_mesh_2d(uv)
    mat = dict(BE_U=float(bend), mass=node_mass(n))
    if clean:
        mat.update(BV=0.0, air_damping=0.0, internal_damping=0.0)
    if extra:
        mat.update(extra)
    d.set_material(**mat)
    d.set_gravity((0.0, -float(gravity), 0.0))
    d.set_timestep(float(h))
    d.set_dynamic_pcg(True, PCG_TOL, PCG_CAP)
    skipped = []
    if clean:
        for name, call in (('set_damping', lambda: d.set_damping(0.0)), ('set_air_drag', lambda: d.set_air_drag(0.0)),
                           ('set_bend_damping', lambda: d.set_bend_damping(0))):
            try:
                call()
            except Exception as e:   # 지원하지 않는 호출은 건너뛰고 기록한다
                skipped.append('%s: %s' % (name, str(e)[:80]))
    for p in pins:
        d.pin_vertex(p)
    return d, V, s, row, skipped


def row_mean(n, a):
    return a.reshape(n, n).mean(axis=1)


def pcg_flags(d):
    hist = d.get_pcg_solve_history()
    return (sum(int(r['iterations']) for r in hist), any(int(r['attempt']) > 0 for r in hist),
            any(r['exit_reason'] == 'numerical_failure' for r in hist))


def child_matrix(a, out):
    n, bend = a['n'], a['bend']
    d, V, s, row, _ = make(a['dll'], n, bend, 1.0 / 200, False, 0.0)
    m = node_mass(n)
    res = {'build_info': d.build_info, 'threads': int(d.openmp_threads)}
    try:
        down = np.array([0.0, -1.0, 0.0])
        K = np.zeros((n - 1, n - 1))
        for j in range(1, n):
            e = np.zeros((n * n, 3))
            e[row == j] = down
            p = np.asarray(d.debug_stiffness_product(V, e), float)
            K[:, j - 1] = (p[:, 1].reshape(n, n).sum(axis=1))[1:]        # −eᵢᵀ J eⱼ = Σ_줄 i p_y
        fields = {'linear': (s, PLATE), 'quadratic': (s * s / PLATE, 4.0 * PLATE / 3.0),
                  'sine': (PLATE * np.sin(np.pi * s / (2 * PLATE)), np.pi ** 2 * PLATE / 8.0)}
        cbar, quot = {}, {}
        for name, (u, grad2) in fields.items():
            f = u[:, None] * down[None, :]
            p = np.asarray(d.debug_stiffness_product(V, f), float)
            q = -float((f * p).sum())
            quot[name] = q
            cbar[name] = q / (PLATE * grad2)
        res['rayleigh_flat'] = quot['linear'] / (m * float((s * s).sum()))
        ph = math.radians(45.0)
        X = V.copy()
        X[:, 1] = 20.0 - s * math.sin(ph)
        X[:, 2] = -PLATE / 2 + s * math.cos(ph)
        f = np.column_stack([np.zeros_like(s), -s * math.cos(ph), -s * math.sin(ph)])
        p = np.asarray(d.debug_stiffness_product(X, f), float)
        res['rayleigh_rot45'] = -float((f * p).sum()) / (m * float((s * s).sum()))
        res['cbar'] = cbar
    finally:
        d.close()
    d, V, s, row, _ = make(a['dll'], n, bend, 1.0 / 200, False, 0.0, extra=dict(E=0.0, SE=0.0))   # 굽힘만
    try:
        f = s[:, None] * np.array([0.0, -1.0, 0.0])[None, :]
        p = np.asarray(d.debug_stiffness_product(V, f), float)
        res['rayleigh_flat_bend_only'] = -float((f * p).sum()) / (m * float((s * s).sum()))
    finally:
        d.close()
    np.save(out + '.K.npy', K)
    with open(out, 'w', encoding='utf-8') as fh:
        json.dump(res, fh)


def child_decay(a, out):
    n, sps, N = a['n'], a['sps'], a['steps']
    w0 = np.load(a['w0'])                                # 줄(1 … n−1)마다의 법선 속도
    d, V, s, row, skipped = make(a['dll'], n, a['bend'], 1.0 / sps, a.get('clean', True), 0.0)   # clean: 2 차에서 더함
    res = {'skipped': skipped, 'build_info': d.build_info}
    try:
        d.step(1)
        v0 = np.zeros((n * n, 3))
        v0[:, 1] = -np.concatenate([[0.0], w0])[row]
        d.set_velocities(v0)
        W = [row_mean(n, -np.asarray(d.velocities(), float)[:, 1])[1:]]
        iters, retry, fail, inplane = [], False, False, 0.0
        for _ in range(N):
            d.step(1)
            v = np.asarray(d.velocities(), float)
            W.append(row_mean(n, -v[:, 1])[1:])
            it, r, f = pcg_flags(d)
            iters.append(it)
            retry, fail = retry or r, fail or f
            inplane = max(inplane, float(np.abs(v[:, [0, 2]]).max()))
        res.update({'W': np.array(W).tolist(), 'pcg_iters_max': int(max(iters)), 'pcg_retry': retry, 'pcg_failure': fail,
                    'inplane_speed_max': inplane})
    finally:
        d.close()
    with open(out, 'w', encoding='utf-8') as fh:
        json.dump(res, fh)


def child_drop(a, out):
    n, sps = a['n'], a['sps']
    d, V, s, row, skipped = make(a['dll'], n, a['bend'], 1.0 / sps, a['clean'], GRAV)
    if a.get('air_drag') is not None:               # 2 차에서 더함 — 공기 저항만 바꾼 회
        d.set_air_drag(float(a['air_drag']))
    hinge = np.array([0.0, 20.0, -PLATE / 2])
    tip = row == n - 1
    res = {'skipped': skipped, 'build_info': d.build_info}
    stop, steps, vtip, vmax, chord, iters = '', 0, [], [], [], []
    try:
        for k in range(1, int(round(B2_T * sps)) + 1):
            try:
                d.step(1)
            except Exception as e:   # StepDiverged 등
                stop = 'diverged: %s' % str(e)[:100]
                break
            X = np.asarray(d.positions(), float)
            v = np.asarray(d.velocities(), float)
            if not (np.isfinite(X).all() and np.isfinite(v).all()):
                stop = 'nonfinite'
                break
            it, r, f = pcg_flags(d)
            iters.append(it)
            if f:
                stop = 'pcg_numerical_failure'
                break
            if r:
                stop = 'pcg_retry'
                break
            steps = k
            sp = np.linalg.norm(v, axis=1)
            vtip.append(float(sp[tip].mean()))
            vmax.append(float(sp.max()))
            chord.append(float(np.linalg.norm(X[tip].mean(axis=0) - hinge) / PLATE))
    finally:
        d.close()
    kp = int(np.argmax(vtip)) if vtip else 0
    res.update({'usable': stop == '', 'stop': stop, 'steps': steps, 'v_tip_peak': max(vtip) if vtip else 0.0,
                't_peak': (kp + 1) / sps, 'v_node_max': max(vmax) if vmax else 0.0,
                'chord_min': min(chord) if chord else None, 'chord_at_peak': chord[kp] if chord else None,
                'pcg_iters_mean': float(np.mean(iters)) if iters else 0.0, 'pcg_iters_max': int(max(iters)) if iters else 0,
                'v_tip_series': vtip[::max(1, sps // 100)]})
    with open(out, 'w', encoding='utf-8') as fh:
        json.dump(res, fh)


def spawn(td, arm, kind, args):
    tag = '%s_%s_%s' % (kind, arm, uuid.uuid4().hex[:12])
    out = os.path.join(td, tag + '.json')
    arg = os.path.join(td, tag + '.arg.json')
    with open(arg, 'w', encoding='utf-8') as fh:
        json.dump(dict(args, kind=kind, out=out), fh)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1')
    env.pop('CFMS_MGCG_ADAPTIVE', None)
    env['CFMS_CPU_SXJ_BEND'] = arm           # 끔도 'off' 로 명시한다(브랜치의 기본은 gn_state)
    env['CFMS_CPU_SXJ_MEMBRANE'] = 'off'     # 인장 · 전단은 엔진 손 커널(L4 진단과 같다)
    r = subprocess.run([sys.executable, os.path.abspath(__file__), '--child', arg], env=env, cwd=td, capture_output=True,
                       text=True, encoding='utf-8', errors='replace')
    if r.returncode != 0 or not os.path.isfile(out):
        raise RuntimeError('회 실패 (%s, %s, %s): %s' % (arm, kind, args, r.stderr[-1500:]))
    with open(out, encoding='utf-8') as fh:
        res = json.load(fh)
    if kind == 'matrix':
        res['K'] = np.load(out + '.K.npy')
    return res



def run_all(td, jobs):
    """jobs: {키: (팔, 종류, 인자)} → {키: 결과}"""
    got = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=PARALLEL) as ex:
        futs = {ex.submit(spawn, td, arm, kind, args): key for key, (arm, kind, args) in jobs.items()}
        for fu in concurrent.futures.as_completed(futs):
            got[futs[fu]] = fu.result()
    return got


def eig(K, n):
    """K̃_rows u = λ M_rows u. M-정규 직교 모드(열)와 오름차순 λ."""
    Ks = 0.5 * (K + K.T)
    mrow = n * node_mass(n)
    lam, U = np.linalg.eigh(Ks / mrow)
    return lam, U / math.sqrt(mrow), mrow


def gate(fail, ok, label, text):
    print('  %-4s %s  %s' % ('PASS' if ok else 'FAIL', label, text))
    if not ok:
        fail.append(label)
    return bool(ok)


def main():
    dll = os.path.abspath(sys.argv[1])
    if not os.path.isfile(dll):
        print('DLL 없음:', dll)
        return 2
    td = tempfile.mkdtemp(prefix='cfmsSXJ_jdamp_')
    t0 = time.perf_counter()
    fail, out = [], {'dll': dll, 'plate_cm': PLATE, 'gravity': GRAV, 'pcg': [PCG_TOL, PCG_CAP]}
    try:
        # ---- A. 행렬
        mats = run_all(td, {(arm, n, b): (arm, 'matrix', dict(dll=dll, n=n, bend=b)) for arm in ARMS for n, b in MATRIX_CASES})
        info = mats[('off',) + BASE]['build_info']
        out['engine'] = {'git_commit': info.get('git_commit'), 'dirty': info.get('dirty')}
        A = {}
        for (arm, n, b), r in mats.items():
            lam, U, mrow = eig(r['K'], n)
            A['%s|%d|%g' % (arm, n, b)] = {'rayleigh_flat': r['rayleigh_flat'], 'rayleigh_rot45': r['rayleigh_rot45'],
                                           'rayleigh_flat_bend_only': r['rayleigh_flat_bend_only'], 'cbar': r['cbar'],
                                           'lambda_1_3': [float(x) for x in lam[:3]],
                                           'K_asymmetry': float(np.abs(r['K'] - r['K'].T).max() / max(np.abs(r['K']).max(), 1e-300))}
        out['A'] = A
        oh, gn = mats[('off',) + BASE], mats[('gn_state',) + BASE]
        print('A. 행렬의 성질 (엔진 %s, dirty %s)' % (info.get('git_commit'), info.get('dirty')))
        print('  Rayleigh 몫 R(δ) [1/s²]: Oh %.6g (굽힘만 %.6g, 45° %.6g) · GN %.3g' % (
            oh['rayleigh_flat'], oh['rayleigh_flat_bend_only'], oh['rayleigh_rot45'], gn['rayleigh_flat']))
        zero = not (oh['rayleigh_flat'] > 0.0)
        gate(fail, (not zero) and abs(gn['rayleigh_flat']) <= 1e-6 * oh['rayleigh_flat'], 'A1',
             'GN / Oh = %.2e (≤ 1e-6)' % (abs(gn['rayleigh_flat']) / oh['rayleigh_flat'] if not zero else float('nan')))
        cb = oh['cbar']
        spread = (max(cb.values()) - min(cb.values())) / max(cb.values()) if not zero else float('nan')
        gate(fail, (not zero) and spread <= 0.05, 'A2', 'c̄ = %s dyn/cm, 폭 %.2f %% (≤ 5 %%)' % (
            ' · '.join('%s %.5g' % kv for kv in cb.items()), 100 * spread))
        cl = {b: mats[('off', 81, b)]['cbar']['linear'] for b in (1e3, 1e4, 1e5)}
        prop = max(abs(cl[b] / b - cl[1e4] / 1e4) / (cl[1e4] / 1e4) for b in cl) if not zero else float('nan')
        kap = {n: mats[('off', n, 1e4)]['cbar']['linear'] * (PLATE / (n - 1)) ** 2 / 1e4 for n in (21, 41, 81)}
        ksp = (max(kap.values()) - min(kap.values())) / max(kap.values()) if not zero else float('nan')
        gate(fail, (not zero) and prop <= 1e-6 and ksp <= 0.10, 'A3', 'BE_U 비례 편차 %.1e (≤ 1e-6), κ = %s, 폭 %.1f %% (≤ 10 %%)' % (
            prop, ' · '.join('n%d %.4g' % kv for kv in sorted(kap.items())), 100 * ksp))
        rot = abs(oh['rayleigh_rot45'] - oh['rayleigh_flat']) / oh['rayleigh_flat'] if not zero else float('nan')
        gate(fail, (not zero) and rot <= 1e-6, 'A4', '45° 회전 뒤 상대 차 %.1e (≤ 1e-6)' % rot)
        out['A_summary'] = {'cbar_spread': spread, 'be_proportionality': prop, 'kappa': {str(k): v for k, v in kap.items()},
                            'kappa_spread': ksp, 'rot45_rel': rot}

        # ---- 예측 (B 를 돌리기 전에 찍는다)
        n0, b0 = BASE
        lam, U, mrow = eig(oh['K'], n0)
        srow = np.linspace(0.0, PLATE, n0)[1:]
        pred = {'lambda': [float(x) for x in lam[:6]], 'B0': {}, 'B1': {}, 'B2': {}}
        b0_jobs = {}
        for sps in B0_SPS:
            for k in B0_MODES:
                fac = 1.0 / (1.0 + lam[k - 1] / sps ** 2)
                N = int(min(200, max(5, math.ceil(math.log(0.5) / math.log(fac))))) if fac < 1.0 else 200
                pred['B0']['%d|%d' % (sps, k)] = {'factor': fac, 'steps': N}
                u = U[:, k - 1]
                w0 = os.path.join(td, 'w0_%d_%d.npy' % (sps, k))
                np.save(w0, u / np.abs(u).max())
                for arm in ARMS:
                    b0_jobs[(arm, sps, k)] = (arm, 'decay', dict(dll=dll, n=n0, bend=b0, sps=sps, steps=N, w0=w0))
        w0r = B1_RATE * srow
        bk = U.T @ (mrow * w0r)
        dec = 1.0 / (1.0 + lam / B1_SPS ** 2)
        tipser = np.array([float((bk * U[-1, :] * dec ** k).sum()) for k in range(0, 2001)]) / float((bk * U[-1, :]).sum())
        n_half = int(np.argmax(tipser <= 0.5)) if (tipser <= 0.5).any() else 2000
        n_fifth = int(np.argmax(tipser <= 0.2)) if (tipser <= 0.2).any() else 2000
        pred['B1'] = {'n_half': n_half, 'ratio_half': float(tipser[n_half]), 'n_fifth': n_fifth, 'ratio_fifth': float(tipser[n_fifth])}
        w0f = os.path.join(td, 'w0_rigid.npy')
        np.save(w0f, w0r)
        b1_jobs = {arm: (arm, 'decay', dict(dll=dll, n=n0, bend=b0, sps=B1_SPS, steps=min(n_fifth, 400), w0=w0f)) for arm in ARMS}
        for n, b, sps in B2_CASES:
            l2, _, mr = eig(mats[('off', n, b)]['K'], n)
            Ks = 0.5 * (mats[('off', n, b)]['K'] + mats[('off', n, b)]['K'].T)
            vss = float(np.linalg.solve(Ks, np.full(n - 1, mr * GRAV))[-1]) * sps
            pred['B2']['%d|%g|%d' % (n, b, sps)] = {'v_ss': vss, 'h_lambda1': float(l2[0]) / sps, 'terminal': float(l2[0]) / sps >= TERMINAL}
        out['predictions'] = pred
        print('예측 (Oh 팔, n = %d · BE_U %g): λ1 – λ3 = %.5g · %.5g · %.5g 1/s²' % (n0, b0, lam[0], lam[1], lam[2]))
        for key, p in pred['B0'].items():
            print('  B0 %s: 단계당 인자 %.5f, N = %d' % (key, p['factor'], p['steps']))
        print('  B1: 비 0.5 → %d 단계(%.4f), 0.2 → %d 단계(%.4f)' % (n_half, tipser[n_half], n_fifth, tipser[n_fifth]))
        for key, p in pred['B2'].items():
            print('  B2 %s: V_ss %.4g cm/s, h λ1 %.4g 1/s%s' % (key, p['v_ss'], p['h_lambda1'], ' (종단)' if p['terminal'] else ''))
        with open(OUT + '.pred', 'w', encoding='utf-8') as fh:   # 실행 전 예측의 사본(끝나면 결과 json 에 든다)
            json.dump(pred, fh, indent=1)

        # ---- B0 · B1 · B2 · C
        b2_jobs = {}
        for n, b, sps in B2_CASES:
            for arm in ARMS:
                b2_jobs[(arm, n, b, sps, False)] = (arm, 'drop', dict(dll=dll, n=n, bend=b, sps=sps, clean=False))
        for n, b, sps in B2_CLEAN:
            for arm in ARMS:
                b2_jobs[(arm, n, b, sps, True)] = (arm, 'drop', dict(dll=dll, n=n, bend=b, sps=sps, clean=True))
        jobs = {('B0',) + k: v for k, v in b0_jobs.items()}
        jobs.update({('B1', k): v for k, v in b1_jobs.items()})
        jobs.update({('B2',) + k: v for k, v in b2_jobs.items()})
        got = run_all(td, jobs)

        print('B0. 한 모드의 감쇠 (무중력, 깨끗한 설정)')
        B0, ok0 = {}, True
        for sps in B0_SPS:
            for k in B0_MODES:
                p = pred['B0']['%d|%d' % (sps, k)]
                u, N = U[:, k - 1], p['steps']
                row = {}
                for arm in ARMS:
                    r = got[('B0', arm, sps, k)]
                    W = np.array(r['W'])
                    c = (W * (mrow * u)[None, :]).sum(axis=1)
                    row[arm] = {'factor': float((c[N] / c[0]) ** (1.0 / N)) if c[N] / c[0] > 0 else float('nan'),
                                'pcg_retry': r['pcg_retry'], 'pcg_failure': r['pcg_failure'], 'inplane_speed_max': r['inplane_speed_max'],
                                'skipped': r['skipped']}
                ratio = math.log(row['off']['factor']) / math.log(p['factor'])
                gnr = abs(math.log(row['gn_state']['factor'])) / abs(math.log(p['factor']))
                good = 0.95 <= ratio <= 1.05 and gnr <= 0.02 and not (row['off']['pcg_retry'] or row['off']['pcg_failure'])
                ok0 = ok0 and good
                B0['%d|%d' % (sps, k)] = {'predicted': p['factor'], 'steps': N, 'oh': row['off'], 'gn': row['gn_state'],
                                          'log_ratio': ratio, 'gn_exponent_over_oh_predicted': gnr}
                print('  1/%-3d 모드 %d: 예측 %.5f · Oh %.5f (ln 비 %.4f) · GN %.6f (지수 비 %.2e)%s' % (
                    sps, k, p['factor'], row['off']['factor'], ratio, row['gn_state']['factor'], gnr, '' if good else '  ← 벗어남'))
        out['B0'] = B0
        gate(fail, ok0, 'B0', 'Oh 의 ln 비 ∈ [0.95, 1.05] 이고 GN 지수 ≤ Oh 예측의 2 %')

        print('B1. 강체 회전의 감쇠 (무중력, 깨끗한 설정, 1/%d)' % B1_SPS)
        B1, ok1 = {}, True
        for arm in ARMS:
            W = np.array(got[('B1', arm)]['W'])
            ser = W[:, -1] / W[0, -1]
            B1[arm] = {'ratio_half': float(ser[min(n_half, len(ser) - 1)]), 'ratio_fifth': float(ser[min(n_fifth, len(ser) - 1)]),
                       'steps_run': int(len(ser) - 1)}
        for name, npt in (('half', n_half), ('fifth', n_fifth)):
            rel = B1['off']['ratio_' + name] / pred['B1']['ratio_' + name]
            ok1 = ok1 and 0.9 <= rel <= 1.1 and npt <= B1['off']['steps_run']
            print('  %d 단계: 예측 비 %.4f · Oh %.4f (측정 / 예측 %.4f) · GN %.4f' % (
                npt, pred['B1']['ratio_' + name], B1['off']['ratio_' + name], rel, B1['gn_state']['ratio_' + name]))
        out['B1'] = B1
        gate(fail, ok1, 'B1', '두 단계에서 측정 비 / 예측 비 ∈ [0.9, 1.1]')

        print('B2. 중력 아래의 종단 속도 (끝 줄 평균 속력의 1 초 안 최댓값, cm/s)')
        B2, ok2, nterm = {}, True, 0
        for clean in (False, True):
            for n, b, sps in (B2_CLEAN if clean else B2_CASES):
                p = pred['B2']['%d|%g|%d' % (n, b, sps)]
                ro, rg = got[('B2', 'off', n, b, sps, clean)], got[('B2', 'gn_state', n, b, sps, clean)]
                ratio = ro['v_tip_peak'] / p['v_ss'] if ro['usable'] else float('nan')
                judged = p['terminal'] and ro['usable'] and not clean
                good = (0.8 <= ratio <= 1.25) if judged else True
                nterm += int(judged)
                ok2 = ok2 and good
                B2['%d|%g|%d|%s' % (n, b, sps, 'clean' if clean else 'orig')] = {
                    'v_ss': p['v_ss'], 'h_lambda1': p['h_lambda1'], 'terminal': p['terminal'], 'judged': judged, 'ratio': ratio,
                    'oh': {k: ro[k] for k in ('usable', 'stop', 'v_tip_peak', 't_peak', 'v_node_max', 'chord_min', 'chord_at_peak',
                                              'pcg_iters_mean', 'pcg_iters_max', 'skipped')},
                    'gn': {k: rg[k] for k in ('usable', 'stop', 'v_tip_peak', 't_peak', 'v_node_max', 'chord_min', 'chord_at_peak',
                                              'pcg_iters_mean', 'pcg_iters_max', 'skipped')}}
                print('  n %-2d BE_U %-6g 1/%-3d %-5s hλ1 %8.4g %s 예측 %8.4g · Oh %8.4g%s (비 %.3f%s) · GN %8.4g%s 노드 최대 Oh %.4g GN %.4g' % (
                    n, b, sps, 'clean' if clean else 'orig', p['h_lambda1'], '종단' if p['terminal'] else '    ', p['v_ss'],
                    ro['v_tip_peak'], '' if ro['usable'] else '[' + ro['stop'][:12] + ']', ratio, '' if good else ' ← 벗어남',
                    rg['v_tip_peak'], '' if rg['usable'] else '[' + rg['stop'][:12] + ']', ro['v_node_max'], rg['v_node_max']))
        out['B2'] = B2
        gate(fail, ok2 and nterm > 0, 'B2', '종단 영역의 Oh 회 %d 개 모두 0.8 ≤ V_peak / V_ss ≤ 1.25' % nterm)

        rc = got[('B2', 'gn_state') + C_CASE + (False,)]
        out['C'] = {'chord_min': rc['chord_min'], 'chord_at_peak': rc['chord_at_peak'], 'v_tip_peak': rc['v_tip_peak'],
                    'v_node_max': rc['v_node_max'], 't_peak': rc['t_peak'], 'rigid_plate_tip_speed': math.sqrt(3.0 * GRAV * PLATE),
                    'usable': rc['usable'], 'v_tip_series': rc['v_tip_series']}
        print('C. 판이 굽는가 (GN, n 81 · BE_U 1e4 · 1/800, 보고): 끝 줄까지의 직선 거리 / L 최소 %.3f (최대 속력 때 %.3f), '
              '끝 줄 평균 속력 최대 %.4g · 노드 최대 %.4g cm/s (t = %.3f s), √(3gL) = %.4g' % (
                  rc['chord_min'], rc['chord_at_peak'], rc['v_tip_peak'], rc['v_node_max'], rc['t_peak'], math.sqrt(3.0 * GRAV * PLATE)))
    finally:
        shutil.rmtree(td, ignore_errors=True)
    out['fail'] = fail
    out['verdict'] = 'PASS' if not fail else 'FAIL'
    out['seconds'] = time.perf_counter() - t0
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
        fh.write('\n')
    if os.path.isfile(OUT + '.pred'):
        os.remove(OUT + '.pred')
    print('종합: %s%s  (%.0f s, 결과 %s)' % (out['verdict'], '' if not fail else ' — ' + ' '.join(fail), out['seconds'],
                                        os.path.relpath(OUT, ROOT)))
    return 0 if not fail else 1


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--child':
        with open(sys.argv[2], encoding='utf-8') as fh:
            a = json.load(fh)
        {'matrix': child_matrix, 'decay': child_decay, 'drop': child_drop}[a['kind']](a, a['out'])
        sys.exit(0)
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main())
