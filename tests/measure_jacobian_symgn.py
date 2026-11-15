"""논문용 측정 4 (docs/GATE_JACOBIAN_DAMPING.md 「논문용 측정 4」, 판정 아님): 블록 대칭으로 만든 Gauss–Newton 행렬(symgn_state).

사용:  python tests/measure_jacobian_symgn.py <실험 cfmsDrape.dll>     (리포 루트에서, run_gate 밖 — 약 10 분)
실험 DLL = 엔진 브랜치에 symgn_state 를 넣은 빌드(%TEMP%\\cfmsDrape_build_sxjsymgn). 먼저 tests/engine_bend_symgn.py 가 통과해야 한다.
팔 셋: off(블록 대칭 근사, 식 (10)) · gn_state · symgn_state. 깨끗한 설정(tests/gate_jacobian_damping.py 의 회).
  1. 평평한 천의 행렬 — 강체 회전의 Rayleigh 몫 · 장력 계수 (적분 없음)
  2. 접힌 천의 행렬 — 주름진 판(굽힘만)에서 강체 회전의 Rayleigh 몫과 −J 의 가장 작은 / 큰 고윳값 (적분 없음)
  3. 낙하 — 끝 줄 평균 속력의 최댓값
전제: off · gn_state 낙하 회가 results/measure_jacobian_drop.json 의 같은 회와 결정적 값이 끝자리까지 같다.
결과 results/measure_jacobian_symgn.json. 전제가 서면 0, 깨지면 1, 준비 안 됨 2.
"""
import sys
sys.dont_write_bytecode = True

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
import gate_jacobian_damping as G  # noqa: E402
import link_l0 as L0  # noqa: E402

OUT = os.path.join(G.ROOT, 'results', 'measure_jacobian_symgn.json')
REC = os.path.join(G.ROOT, 'results', 'measure_jacobian_drop.json')
GATE = os.path.join(G.ROOT, 'results', 'gate_engine_bend_symgn.json')
ARMS = ('off', 'gn_state', 'symgn_state')
N, BEND = 81, 1e4
SPS_ALL = (50, 100, 200, 400, 800, 1600, 3200)
EXTRA = [(81, 1e3, 200), (81, 1e5, 200), (41, 1e4, 200), (21, 1e4, 200)]      # symgn_state 만 — gn_state 는 기록과 견준다
FOLD_N, FOLD_AMP = 21, (0.0, 0.25, 0.5, 1.0, 2.0)                              # 주름의 진폭(cm), 파장 15 cm
FP = ('usable', 'v_tip_peak', 'pcg_iters_mean', 'pcg_iters_max')
KEEP = ('usable', 'stop', 'steps', 'v_tip_peak', 't_peak', 'v_node_max', 'chord_min', 'chord_at_peak', 'pcg_iters_mean', 'pcg_iters_max', 'skipped')


def fold_angles(X, F):
    """힌지마다의 접힌 각(도). 명세의 꼭짓점 순서(h, i, j, k): nA = (xi − xh) × (xj − xh), nB = (xj − xk) × (xi − xk)."""
    out = []
    for h, i, j, k in L0.hinges(F):
        e = X[j] - X[i]
        nA, nB = np.cross(X[i] - X[h], X[j] - X[h]), np.cross(X[j] - X[k], X[i] - X[k])
        nA, nB, e = nA / np.linalg.norm(nA), nB / np.linalg.norm(nB), e / np.linalg.norm(e)
        out.append(math.degrees(math.atan2(float(np.dot(np.cross(nA, nB), e)), float(np.dot(nA, nB)))))
    return np.array(out)


def child_folded(a, out):
    """주름진 판(굽힘만)에서 강체 회전의 Rayleigh 몫과 굽힘 행렬의 고윳값."""
    n = a['n']
    d, V, s, row, _ = G.make(a['dll'], n, a['bend'], 1.0 / 200, False, 0.0, extra=dict(E=0.0, SE=0.0))
    m = G.node_mass(n)
    _, F, _, _, _, _ = G.scene(n)
    res = {'build_info': d.build_info, 'cases': []}
    wave = V[:, 0] if a.get('orient', 'x') == 'x' else V[:, 2]      # 'x' = 주름 선이 고정 변과 수직(측정 4), 'z' = 나란히(판정 뒤 진단)
    try:
        for amp in a['amps']:
            X = V.copy()
            X[:, 1] = 20.0 + amp * np.sin(2.0 * math.pi * 2.0 * (wave + G.PLATE / 2) / G.PLATE)
            # 고정 변(z = −L/2, y = 20 을 지나는 x 축)을 축으로 한 강체 회전: δ = x̂ × r = (0, −(z + L/2), y − 20)
            delta = np.column_stack([np.zeros(len(X)), -s, X[:, 1] - 20.0])
            p = np.asarray(d.debug_stiffness_product(X, delta), float)
            rot = -float((delta * p).sum()) / (m * float((delta * delta).sum()))
            n3 = 3 * len(X)
            J = np.zeros((n3, n3))
            for k in range(n3):
                e = np.zeros((len(X), 3))
                e[k // 3, k % 3] = 1.0
                J[:, k] = np.asarray(d.debug_stiffness_product(X, e), float).ravel()
            lam = np.linalg.eigvalsh(-0.5 * (J + J.T))
            th = fold_angles(X, F)
            res['cases'].append({'amp_cm': amp, 'rayleigh_rot': rot, 'eig_min': float(lam.min()), 'eig_max': float(lam.max()),
                                 'asymmetry': float(np.abs(J - J.T).max() / max(np.abs(J).max(), 1e-300)),
                                 'theta_rms_deg': float(np.sqrt((th ** 2).mean())), 'theta_max_deg': float(np.abs(th).max())})
    finally:
        d.close()
    with open(out, 'w', encoding='utf-8') as fh:
        json.dump(res, fh)


def spawn_folded(td, arm, args):
    tag = 'folded_%s_%s' % (arm, uuid.uuid4().hex[:12])
    out, arg = os.path.join(td, tag + '.json'), os.path.join(td, tag + '.arg.json')
    with open(arg, 'w', encoding='utf-8') as fh:
        json.dump(dict(args, out=out), fh)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1', CFMS_CPU_SXJ_BEND=arm, CFMS_CPU_SXJ_MEMBRANE='off')
    env.pop('CFMS_MGCG_ADAPTIVE', None)
    r = subprocess.run([sys.executable, os.path.abspath(__file__), '--child-folded', arg], env=env, cwd=td, capture_output=True, text=True,
                       encoding='utf-8', errors='replace')
    if r.returncode != 0 or not os.path.isfile(out):
        raise RuntimeError('접힌 천 회 실패 (%s): %s' % (arm, r.stderr[-1500:]))
    with open(out, encoding='utf-8') as fh:
        return json.load(fh)


def main():
    dll = os.path.abspath(sys.argv[1])
    if not all(os.path.isfile(p) for p in (dll, REC, GATE)):
        print('없는 것:', [p for p in (dll, REC, GATE) if not os.path.isfile(p)])
        return 2
    with open(GATE, encoding='utf-8') as fh:
        gate = json.load(fh)
    if gate.get('gates_failed') != 0 or os.path.normcase(gate['dll']['new']) != os.path.normcase(dll):
        print('커널 확인(tests/engine_bend_symgn.py)이 이 DLL 로 통과하지 않았다')
        return 2
    with open(REC, encoding='utf-8') as fh:
        rec = {(r['n'], r['BE_U'], r['sps']): r for r in json.load(fh)['runs']}
    t0 = time.perf_counter()
    td = tempfile.mkdtemp(prefix='cfmsSXJ_symgn_')
    try:
        print('논문용 측정 4 — 블록 대칭 Gauss–Newton (symgn_state)')
        mats = G.run_all(td, {arm: (arm, 'matrix', dict(dll=dll, n=N, bend=BEND)) for arm in ARMS})
        folded = {arm: spawn_folded(td, arm, dict(dll=dll, n=FOLD_N, bend=BEND, amps=list(FOLD_AMP))) for arm in ARMS}
        jobs = {(arm, N, BEND, s): (arm, 'drop', dict(dll=dll, n=N, bend=BEND, sps=s, clean=True)) for arm in ARMS for s in SPS_ALL}
        for n, b, s in EXTRA:
            jobs[('symgn_state', n, b, s)] = ('symgn_state', 'drop', dict(dll=dll, n=n, bend=b, sps=s, clean=True))
        print(' 낙하 %d 회 (동시 프로세스 %d)' % (len(jobs), G.PARALLEL))
        got = G.run_all(td, jobs)
    finally:
        shutil.rmtree(td, ignore_errors=True)
    info = got[('off', N, BEND, 200)]['build_info']

    # 전제 (b): 기존 두 팔이 기록과 같다
    name, bad = {'off': 'oh', 'gn_state': 'gn'}, []
    for arm in name:
        for s in SPS_ALL:
            want, r = rec[(N, BEND, s)][name[arm]], got[(arm, N, BEND, s)]
            if any(r[k] != want[k] for k in FP):
                bad.append('%s|1/%d' % (arm, s))
    valid = not bad
    print(' 전제 %s  off · gn_state %d 회의 결정적 값이 기록과 %s (엔진 %s, dirty %s)' % (
        'PASS' if valid else 'FAIL', 2 * len(SPS_ALL), '모두 같음' if valid else '다름 %s' % bad, info.get('git_commit'), info.get('dirty')))

    print(' 1. 평평한 천 (n %d, 굽힘 %g erg) — 강체 회전의 Rayleigh 몫 [1/s²] · 장력 계수 c̄ [dyn/cm]' % (N, BEND))
    flat = {}
    for arm in ARMS:
        r = mats[arm]
        flat[arm] = {'rayleigh_flat': r['rayleigh_flat'], 'rayleigh_flat_bend_only': r['rayleigh_flat_bend_only'],
                     'rayleigh_rot45': r['rayleigh_rot45'], 'cbar': r['cbar']}
        print('  %-12s R %.6g (굽힘만 %.6g, 45° 기운 판 %.6g) · c̄ 선형 %.6g' % (
            arm, r['rayleigh_flat'], r['rayleigh_flat_bend_only'], r['rayleigh_rot45'], r['cbar']['linear']))
    kdiff = float(np.abs(mats['symgn_state']['K'] - mats['gn_state']['K']).max() / np.abs(mats['gn_state']['K']).max())
    print('  줄 단위 행렬: symgn_state 대 gn_state 의 상대 차 %.2e' % kdiff)

    print(' 2. 접힌 천 (n %d, 굽힘만) — 진폭 · 접힌 각 rms / 최대 · 강체 회전의 R [1/s²] · −J 의 최소 / 최대 고윳값' % FOLD_N)
    for k, amp in enumerate(FOLD_AMP):
        c0 = folded['off']['cases'][k]
        print('  A %.2f cm (접힌 각 rms %.1f° · 최대 %.1f°)' % (amp, c0['theta_rms_deg'], c0['theta_max_deg']))
        for arm in ARMS:
            c = folded[arm]['cases'][k]
            print('     %-12s R %.6g   고윳값 비 %.3e   (블록 비대칭 %.1e)' % (arm, c['rayleigh_rot'], c['eig_min'] / c['eig_max'], c['asymmetry']))

    peak = {arm: [(got[(arm, N, BEND, s)]['v_tip_peak'] if got[(arm, N, BEND, s)]['usable'] else None) for s in SPS_ALL] for arm in ARMS}
    its = {arm: [got[(arm, N, BEND, s)]['pcg_iters_mean'] for s in SPS_ALL] for arm in ARMS}
    rel = [(a / g if (a is not None and g is not None) else None) for a, g in zip(peak['symgn_state'], peak['gn_state'])]
    print(' 3. 낙하 — 끝 줄 평균 속력의 1 초 안 최댓값 (cm/s), n %d · 굽힘 %g erg' % (N, BEND))
    print('  %-12s' % '1/h' + ''.join('%11d' % s for s in SPS_ALL))
    for arm in ARMS:
        cells = []
        for s in SPS_ALL:
            r = got[(arm, N, BEND, s)]
            cells.append('%11s' % ('%.2f' % r['v_tip_peak'] if r['usable'] else '(%s@%d)' % (r['stop'][:6], r['steps'])))
        print('  %-12s' % arm + ''.join(cells))
        print('  %-12s' % '  (PCG 평균)' + ''.join('%11.0f' % v for v in its[arm]))
    print('  %-12s' % 'symgn / gn' + ''.join('%11s' % ('-' if v is None else '%.4f' % v) for v in rel))
    extra = {}
    print('  다른 강성 · 메시 (1/200, symgn_state 대 기록의 gn_state)')
    for n, b, s in EXTRA:
        r, g = got[('symgn_state', n, b, s)], rec[(n, b, s)]['gn']
        extra['%d|%g|%d' % (n, b, s)] = {'symgn': {k: r.get(k) for k in KEEP}, 'gn_record': {k: g.get(k) for k in FP},
                                         'ratio': (r['v_tip_peak'] / g['v_tip_peak'] if r['usable'] else None)}
        print('   n %d · 굽힘 %g: %s (gn_state %.2f, 비 %s, PCG %.0f 대 %.0f)' % (
            n, b, '%.2f' % r['v_tip_peak'] if r['usable'] else '(%s@%d)' % (r['stop'][:8], r['steps']), g['v_tip_peak'],
            '%.4f' % (r['v_tip_peak'] / g['v_tip_peak']) if r['usable'] else '-', r['pcg_iters_mean'], g['pcg_iters_mean']))

    res = {'dll': dll, 'engine': {'git_commit': info.get('git_commit'), 'dirty': info.get('dirty')}, 'arms': list(ARMS),
           'premise': {'old_arms_match_record': valid, 'different': bad}, 'valid': valid,
           'flat': {'n': N, 'BE_U': BEND, 'arms': flat, 'reduced_matrix_symgn_vs_gn_rel': kdiff},
           'folded': {'n': FOLD_N, 'BE_U': BEND, 'wavelength_cm': G.PLATE / 2, 'arms': {arm: folded[arm]['cases'] for arm in ARMS}},
           'drop': {'n': N, 'BE_U': BEND, 'sps': list(SPS_ALL), 'peak': peak, 'pcg_iters_mean': its, 'symgn_over_gn': rel,
                    'runs': {'%s|%d|%g|%d' % k: {kk: r.get(kk) for kk in KEEP} for k, r in sorted(got.items())}, 'extra': extra},
           'seconds': time.perf_counter() - t0}
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
        fh.write('\n')
    print('전제 %s → %s (%.0f s)' % ('섬' if valid else '깨짐 — 값을 쓰지 않는다', os.path.relpath(OUT, G.ROOT), res['seconds']))
    return 0 if valid else 1


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--child-folded':
        with open(sys.argv[2], encoding='utf-8') as fh:
            arg = json.load(fh)
        child_folded(arg, arg['out'])
        sys.exit(0)
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main())
