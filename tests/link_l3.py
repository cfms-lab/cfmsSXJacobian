"""엔진 연결 L3 (docs/GATE_LINK_L3.md): 엔진 브랜치 sxjacobian-link 의 생성 인장 · Provot 전단 갈래를 엔진 손 커널과 견준다.

사용:  python tests/link_l3.py <main cfmsDrape.dll> <브랜치 cfmsDrape.dll>
       (내부) python tests/link_l3.py --measure <dll> <out.npz>
              — 한 DLL 을 다른 프로세스에서 잰다. 모드는 그 프로세스의 환경 변수 CFMS_CPU_SXJ_BEND · CFMS_CPU_SXJ_MEMBRANE 로 준다.
엔진 리포(main)는 CFMS_DRAPE_DIR(없으면 ..\\cfmsDrape_dev), 브랜치 워크트리는 CFMS_DRAPE_LINK_DIR
(없으면 ..\\cfmsDrape_dev_sxjacobian-link). 장면은 tests/link_l0.py 와 같고 속도 탐침만 더한다.
L3-0 (f) 는 results/gate_link_l1.json · gate_link_l2.json 을 읽는다 — 이 스크립트보다 먼저 L1 · L2 를 돌려야 한다.
리포 루트에서 실행. 결과는 results/gate_link_l3.json. 통과 0, 실패 1, 준비 안 됨 2.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import link_l0 as L0  # noqa: E402
import link_l2 as L2  # noqa: E402

ROOT = L0.ROOT
ENGINE_DIR = L0.ENGINE_DIR
LINK_DIR = L2.LINK_DIR
BASE = L0.ENGINE_COMMIT
HEADERS = ('BridsonBend_Engine.h', 'BaraffStretch_Engine.h', 'ProvotShear_Engine.h')
LIM = 1e-12                    # L3-2 · L3-3 · L3-7
LIM_STEP = 1e-9                # L3-4
MIN_BRANCH = 5                 # L3-5
VEL_SEED = L0.SEED + 3         # 속도 탐침(시드 고정)
VEL_MAX = 20.0                 # cm/s
MEMBRANE_MATERIALS = ('stretch', 'shear', 'all')
PROBES = ('v0', 'vel')
B_SHEAR = float(np.float32(1.41421356237))   # 엔진 b_su · b_sv (float 리터럴)
# 살아 있는 매단 천(탐침 상태 X 에서 시작, 중력 −Y — 엔진은 중력을 −Y 로만 당긴다, 윗줄 두 모서리 고정)
HANGX_PINS = (20, 24)
HANGX_GRAVITY = 980.0
HANG_STEPS = L2.HANG_STEPS


def velocity_probe(shape):
    return np.random.default_rng(VEL_SEED).uniform(-VEL_MAX, VEL_MAX, shape)


def new_drape(Drape, dll, V0, F, uv, mat, gravity=(0.0, 0.0, 0.0)):
    d = Drape(dll)
    d.set_mesh(V0, F)
    d.set_mesh_2d(uv)
    d.set_material(mass=0.1, **mat)
    d.set_gravity(gravity)
    d.set_timestep(0.005)
    return d


def probe_measure(Drape, dll, V, F, uv, X, R, mat, vel):
    """디버그 고리로 f · J(단위 벡터 곱) · 무작위 곱을 잰다. vel 이 있으면 set_velocities 뒤에."""
    d = new_drape(Drape, dll, V, F, uv, mat)
    try:
        if vel is not None:
            d.set_velocities(vel)
        out = {'threads': int(d.openmp_threads), 'build_info': d.build_info}
        n3 = X.size
        out['f'] = d.debug_internal_force(X).ravel()
        J = np.empty((n3, n3))
        e = np.zeros(n3)
        for k in range(n3):
            e[:] = 0.0
            e[k] = 1.0
            J[:, k] = d.debug_stiffness_product(X, e.reshape(-1, 3)).ravel()
        out['J'] = J
        out['P'] = np.array([d.debug_stiffness_product(X, r).ravel() for r in R])
        out['f2'] = d.debug_internal_force(X).ravel()
        out['P2'] = np.array([d.debug_stiffness_product(X, r).ravel() for r in R])
        return out
    finally:
        d.close()


def one_step(Drape, dll, V, F, uv, X, mat, vel):
    """L3-4: 속도 탐침에서 20 회 고정 PCG 로 한 단계. 설정 뒤 첫 단계는 PCG 가 1 회에 끝나므로(판정문 「측정 메모」),
    쉬는 평면 천에서 한 단계(아무것도 움직이지 않는다)를 먼저 밟고 탐침 위치 · 속도를 넣는다."""
    d = new_drape(Drape, dll, V, F, uv, mat)
    try:
        d.set_dynamic_pcg(False)
        d.step(1)
        warm = float(max(np.abs(np.array(d.positions()) - V).max(), np.abs(np.array(d.velocities())).max()))
        d.set_positions(X)
        d.set_velocities(vel)
        x0, v0 = np.array(d.positions()), np.array(d.velocities())
        d.step(1)
        h = d.get_pcg_solve_history()[-1]
        return {'v': np.array(d.velocities()).ravel(), 'x': np.array(d.positions()).ravel(),
                'meta': {'warmup_motion': warm, 'start_x_err': float(np.abs(x0 - X).max()),
                         'start_v_err': float(np.abs(v0 - vel).max()), 'iterations': int(h['iterations']),
                         'exit_reason': str(h['exit_reason']), 'iteration_cap': int(h['iteration_cap'])}}
    finally:
        d.close()


def hang_x(Drape, dll, F, uv, X, mat, fixed20):
    d = new_drape(Drape, dll, X, F, uv, mat, gravity=(0.0, -HANGX_GRAVITY, 0.0))
    try:
        if fixed20:
            d.set_dynamic_pcg(False)
        for p in HANGX_PINS:
            d.pin_vertex(p)
        traj, err = [], ''
        for _ in range(HANG_STEPS):
            try:
                d.step(1)
            except Exception as e:   # StepDiverged 등
                err = '%s: %s' % (type(e).__name__, e)
                break
            traj.append(np.array(d.positions(), dtype=float).copy())
        return np.array(traj), err
    finally:
        d.close()


def measure(dll, out):
    Drape = L0.load_drape()
    V, F, uv, X, R = L0.scene()
    vel = velocity_probe(X.shape)
    arrs, info = {}, {'one_step': {}, 'hang_error': {}}
    mats = dict(L0.MATERIALS)
    for name, mat in L0.MATERIALS:
        for probe in PROBES:
            m = probe_measure(Drape, dll, V, F, uv, X, R, mat, vel if probe == 'vel' else None)
            for k in ('f', 'J', 'P', 'f2', 'P2'):
                arrs['%s_%s_%s' % (k, name, probe)] = m[k]
            info['threads'] = m['threads']
            info['build_info'] = m['build_info']
        s = one_step(Drape, dll, V, F, uv, X, mat, vel)
        arrs['stepv_' + name], arrs['stepx_' + name] = s['v'], s['x']
        info['one_step'][name] = s['meta']
    traj, err = L2.hang_run(Drape, dll, V, F, uv, mats['all'])   # L2 와 같은 매단 천
    arrs['hang'], info['hang_error']['l2'] = traj, err
    for key, fixed in (('hangx', False), ('hangx20', True)):
        traj, err = hang_x(Drape, dll, F, uv, X, mats['all'], fixed)
        arrs[key], info['hang_error'][key] = traj, err
    np.savez(out, **arrs)
    with open(out + '.json', 'w', encoding='utf-8') as fh:
        json.dump(info, fh, ensure_ascii=False)
    return 0


def run_measure(dll, bend, membrane, threads1=False):
    """다른 프로세스에서 잰다. bend · membrane None = 그 환경 변수 없음."""
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    for k, v in (('CFMS_CPU_SXJ_BEND', bend), ('CFMS_CPU_SXJ_MEMBRANE', membrane)):
        env.pop(k, None)
        if v is not None:
            env[k] = v
    env['CFMS_CPU_SXJ_MEMBRANE_REPORT'] = '1'   # 브랜치 DLL 이 끝날 때 물러선 요소 수를 stderr 에 적는다(main 은 모른다)
    if threads1:
        env['OMP_NUM_THREADS'] = '1'
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, 'm.npz')
        r = subprocess.run([sys.executable, os.path.abspath(__file__), '--measure', dll, out],
                           env=env, cwd=ROOT, capture_output=True, text=True, encoding='utf-8', errors='replace')
        if r.returncode != 0:
            raise RuntimeError('measure 실패 (bend %s, membrane %s): %s' % (bend, membrane, r.stderr[-2000:]))
        arrs = dict(np.load(out))
        with open(out + '.json', encoding='utf-8') as fh:
            info = json.load(fh)
        info['stderr'] = r.stderr.strip()
        m = re.search(r'SXJ membrane: tension elements (\d+) \(fallback (\d+)\), shear elements (\d+) \(fallback (\d+)\)', r.stderr)
        info['membrane_report'] = ({'tension_elements': int(m.group(1)), 'tension_fallback': int(m.group(2)),
                                    'shear_elements': int(m.group(3)), 'shear_fallback': int(m.group(4))} if m else None)
        return arrs, info


def spec_conditions(X, vel, F, uv):
    """명세 식으로 탐침의 C · Ċ 를 낸다(L3-5 갈래 수, 보고). 성분 순서: 인장 (u, v), 전단 (+, −)."""
    out = {'stretch': {'C': [], 'Cdot': []}, 'shear': {'C': [], 'Cdot': []}}
    frame = []
    for t in F:
        (u0, v0), (u1, v1), (u2, v2) = uv[t[0]], uv[t[1]], uv[t[2]]
        du1, dv1, du2, dv2 = u1 - u0, v1 - v0, u2 - u0, v2 - v0
        det = du1 * dv2 - du2 * dv1
        sa = np.sqrt(abs(det) / 2.0)

        def w(Q):
            d1, d2 = Q[t[1]] - Q[t[0]], Q[t[2]] - Q[t[0]]
            return (d1 * dv2 - d2 * dv1) / det, (-d1 * du2 + d2 * du1) / det
        wu, wv = w(X)
        wud, wvd = w(vel)
        for key, vecs, b in (('stretch', ((wu, wud), (wv, wvd)), (1.0, 1.0)),
                             ('shear', ((wu + wv, wud + wvd), (wu - wv, wud - wvd)), (B_SHEAR, B_SHEAR))):
            for (a, ad), bb in zip(vecs, b):
                n = np.linalg.norm(a)
                out[key]['C'].append(sa * (n - bb))
                out[key]['Cdot'].append(sa * float(a @ ad) / n)
        # 엔진의 쉬는 틀(_Invariant::Prepare): α = ½|(p1 − p0) × (p2 − p0)| (3D), 명세 a = |det| / 2
        p = np.concatenate([uv[t], np.zeros((3, 1))], axis=1)
        alpha = 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
        frame.append(abs(alpha - abs(det) / 2.0) / (abs(det) / 2.0))
    return {k: {q: np.array(v[q]) for q in v} for k, v in out.items()}, float(max(frame))


def gate(fail, ok, label, fmt, value, limit):
    print('  %-12s %s  %s  (기준 %s)' % (label, 'PASS' if ok else 'FAIL', fmt % value, limit))
    return fail + (0 if ok else 1)


def read_json(name):
    p = os.path.join(ROOT, 'results', name)
    if not os.path.isfile(p):
        return None
    with open(p, encoding='utf-8') as fh:
        return json.load(fh)


def main():
    if len(sys.argv) == 4 and sys.argv[1] == '--measure':
        return measure(sys.argv[2], sys.argv[3])
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    dll_main, dll_br = (os.path.abspath(a) for a in sys.argv[1:3])
    gen = [os.path.join(ROOT, 'generated', h) for h in HEADERS]
    missing = [p for p in [dll_main, dll_br] + gen if not os.path.isfile(p)]
    missing += [p for p in (os.path.join(ENGINE_DIR, 'drape'), os.path.join(LINK_DIR, 'drape')) if not os.path.isdir(p)]
    if missing:
        print('없는 것: %s' % ', '.join(missing))
        return 2

    V, F, uv, X, R = L0.scene()
    vel = velocity_probe(X.shape)
    print('엔진 연결 L3 — 노드 %d, 요소 %d\n  main DLL %s\n  브랜치 DLL %s (워크트리 %s)' % (len(V), len(F), dll_main, dll_br, LINK_DIR))

    # 굽힘은 끈다('off'). 인장 · 전단: 'off' = 손 커널, 'on' = 생성 커널. 없음은 브랜치의 기본(on, 사용자 결정 2026-11-13)
    # 이라 끔은 'off' 로 명시한다 — L3-8 이 없음을 'on' 과 견준다(판정문 끝 절).
    meas = {'main': run_measure(dll_main, None, None),
            'unset': run_measure(dll_br, 'off', None),
            'off': run_measure(dll_br, 'off', 'off'),
            'on': run_measure(dll_br, 'off', 'on'),
            'on_t1': run_measure(dll_br, 'off', 'on', threads1=True),
            'gn_off': run_measure(dll_br, None, 'off'),
            'gn_on': run_measure(dll_br, None, 'on')}
    A = {k: v[0] for k, v in meas.items()}
    I = {k: v[1] for k, v in meas.items()}
    names = [n for n, _ in L0.MATERIALS]

    bi_main, bi_br = I['main']['build_info'], I['off']['build_info']
    st = L2.branch_state()
    main_status = L2.git('status', '--porcelain', cwd=ENGINE_DIR).stdout
    sha = {h: {'generated': L2.sha256_file(os.path.join(ROOT, 'generated', h)),
               'engine_copy': (L2.sha256_file(os.path.join(LINK_DIR, 'drape', 'src', 'sde_core', 'sxgen', h))
                               if os.path.isfile(os.path.join(LINK_DIR, 'drape', 'src', 'sde_core', 'sxgen', h)) else '')}
           for h in HEADERS}

    fail = 0
    print(' L3-0 전제')
    a_ok = str(bi_main.get('git_commit', '')).startswith(BASE) and bi_main.get('dirty') is False
    fail = gate(fail, a_ok, '(a)', '%s', 'main DLL git_commit %s · dirty %s' % (bi_main.get('git_commit'), bi_main.get('dirty')), '%s · dirty false' % BASE)
    b_ok = (st['head'] != '' and str(bi_br.get('git_commit', '')).startswith(st['head'][:7]) and st['descends_from_base']
            and not st['outside_allowed'] and bi_br.get('dirty') is False and not st['uncommitted'])
    fail = gate(fail, b_ok, '(b)', '%s', '브랜치 DLL git_commit %s · 워크트리 HEAD %s · %s 자손 %s · 허용 밖 %s · dirty %s · 미커밋 %s' % (
        bi_br.get('git_commit'), st['head'], BASE, st['descends_from_base'], st['outside_allowed'] or '없음', bi_br.get('dirty'),
        st['uncommitted']), 'HEAD 일치 · 자손 · 허용 목록 안 · dirty false')
    fail = gate(fail, main_status.strip() == '', '(c)', '%s', '엔진 main git status %s' % ('비어 있음' if not main_status.strip() else '줄 %d 개' % len(main_status.splitlines())), '비어 있음')
    d_ok = all(s['generated'] == s['engine_copy'] for s in sha.values())
    fail = gate(fail, d_ok, '(d)', '%s', '헤더 사본 셋 %s' % ('같음' if d_ok else '다름 ' + ', '.join(h for h, s in sha.items() if s['generated'] != s['engine_copy'])), 'SHA-256 같음')

    l30e = {}
    for name in names:
        df = L0.rel(A['off']['f_%s_vel' % name], A['off']['f_%s_v0' % name])
        dP = L0.rel(A['off']['P_%s_vel' % name], A['off']['P_%s_v0' % name])
        l30e[name] = {'f': df, 'P': dP}
    e_ok = all(l30e[n]['f'] > 0 and l30e[n]['P'] > 0 for n in MEMBRANE_MATERIALS)
    fail = gate(fail, e_ok, '(e)', '%s', '손 커널(끔) 속도 탐침 대 속도 0: ' + ' · '.join(
        '%s f %.2e 곱 %.2e' % (n, l30e[n]['f'], l30e[n]['P']) for n in MEMBRANE_MATERIALS), '다름')

    l1, l2 = read_json('gate_link_l1.json'), read_json('gate_link_l2.json')
    l1_ok = l1 is not None and l1.get('gates_failed') == 0 and 'L1_5b' in l1
    l2_commit = (l2 or {}).get('build_info_branch', {}).get('git_commit', '')
    l2_ok = (l2 is not None and l2.get('gates_failed') == 0 and 'L2_6' in l2 and st['head'] != ''
             and str(l2_commit).startswith(st['head'][:7]) and l2.get('build_info_branch', {}).get('dirty') is False)
    fail = gate(fail, l1_ok and l2_ok, '(f)', '%s', 'L1 %s · L2 %s (L2 의 브랜치 DLL 커밋 %s)' % (
        'PASS' if l1_ok else 'FAIL/없음', 'PASS' if l2_ok else 'FAIL/없음/다른 커밋', l2_commit or '-'), 'PASS · 같은 브랜치 커밋')

    print(' L3-1 끄면 그대로 (브랜치 굽힘 off · 인장 · 전단 off 대 main, 비트)')
    l31 = {}
    for mode in ('off',):
        for name in names:
            crit = all(L0.same_bits(A['main']['%s_%s_v0' % (k, name)], A[mode]['%s_%s_v0' % (k, name)]) for k in ('f', 'J', 'P'))
            extra = (all(L0.same_bits(A['main']['%s_%s_vel' % (k, name)], A[mode]['%s_%s_vel' % (k, name)]) for k in ('f', 'J', 'P'))
                     and L0.same_bits(A['main']['stepv_' + name], A[mode]['stepv_' + name]))
            l31['%s_%s' % (mode, name)] = {'f_J_P': crit, 'velocity_probe_and_step': extra}
            fail = gate(fail, crit and extra, '%s/%s' % (mode, name), '%s', 'f · J · 무작위 곱 %s · 속도 탐침 · 한 단계 %s' % (
                '같음' if crit else '다름', '같음' if extra else '다름'), '비트까지 같음')
        hang = L0.same_bits(A['main']['hang'], A[mode]['hang']) and A[mode]['hang'].shape[0] == HANG_STEPS
        hangx = all(L0.same_bits(A['main'][k], A[mode][k]) and A[mode][k].shape[0] == HANG_STEPS for k in ('hangx', 'hangx20'))
        l31['%s_hang' % mode] = {'l2_hang': hang, 'hang_from_probe': hangx}
        fail = gate(fail, hang and hangx, '%s/매단천' % mode, '%s', '%d 단계 위치 %s · 탐침에서 시작한 매단 천(수렴 · 20 회) %s' % (
            HANG_STEPS, '같음' if hang else '다름', '같음' if hangx else '다름'), '비트까지 같음')

    print(' L3-2 힘 · L3-3 J (켬 대 손 커널, 같은 브랜치 DLL, 굽힘 off)')
    l32, l33 = {}, {}
    for name in MEMBRANE_MATERIALS:
        for probe in PROBES:
            key = '%s_%s' % (name, probe)
            rf = L0.rel(A['on']['f_' + key], A['off']['f_' + key])
            rJ = L0.rel(A['on']['J_' + key], A['off']['J_' + key])
            rP = L0.rel(A['on']['P_' + key], A['off']['P_' + key])
            bits = L0.same_bits(A['on']['J_' + key], A['off']['J_' + key])
            l32[key] = rf
            l33[key] = {'J': rJ, 'P': rP, 'J_same_bits': bits}
            fail = gate(fail, rf <= LIM, 'L3-2 ' + key, 'f rel %.2e', rf, '%.0e' % LIM)
            fail = gate(fail, rJ <= LIM and rP <= LIM, 'L3-3 ' + key, 'J rel %.2e · 무작위 곱 rel %.2e', (rJ, rP), '%.0e' % LIM)

    print(' L3-4 한 단계 (속도 탐침, 20 회 고정 PCG, 켬 대 끔)')
    l34 = {}
    for name in MEMBRANE_MATERIALS:
        von, voff = A['on']['stepv_' + name], A['off']['stepv_' + name]
        r = float(np.abs(von - voff).max()) / float(np.abs(voff).max())
        mon, moff = I['on']['one_step'][name], I['off']['one_step'][name]
        full = all(m['iterations'] == 20 and m['exit_reason'] == 'iteration_cap' and m['start_x_err'] == 0.0
                   and m['start_v_err'] == 0.0 and m['warmup_motion'] == 0.0 for m in (mon, moff))
        l34[name] = {'dv_over_v': r, 'vmax': float(np.abs(voff).max()), 'on': mon, 'off': moff, 'twenty_iterations_from_probe': full}
        fail = gate(fail, r <= LIM_STEP and full, name, 'max|Δv| / max|v| %.2e · PCG 켬 %d 회(%s) · 끔 %d 회(%s) · 탐침에서 시작 %s', (
            r, mon['iterations'], mon['exit_reason'], moff['iterations'], moff['exit_reason'], full), '%.0e · 20 회' % LIM_STEP)

    print(' L3-5 갈래 (속도 탐침의 요소 성분)')
    cond, frame_diff = spec_conditions(X, vel, F, uv)
    l35 = {}
    for key in ('stretch', 'shear'):
        C, Cd = cond[key]['C'], cond[key]['Cdot']
        n = {'C>0': int((C > 0).sum()), 'C<=0': int((C <= 0).sum()), 'Cdot>0': int((Cd > 0).sum()), 'Cdot<=0': int((Cd <= 0).sum())}
        l35[key] = n
        fail = gate(fail, min(n.values()) >= MIN_BRANCH, key, 'C > 0 %d · C ≤ 0 %d · Ċ > 0 %d · Ċ ≤ 0 %d', (
            n['C>0'], n['C<=0'], n['Cdot>0'], n['Cdot<=0']), '각 %d 이상' % MIN_BRANCH)

    print(' L3-6 결정성 (켠 브랜치: 같은 프로세스 두 번 · 단일 스레드 프로세스 대 기본)')
    l36 = {}
    for name in MEMBRANE_MATERIALS:
        for probe in PROBES:
            key = '%s_%s' % (name, probe)
            a, a1 = A['on'], A['on_t1']
            same_proc = L0.same_bits(a['f_' + key], a['f2_' + key]) and L0.same_bits(a['P_' + key], a['P2_' + key])
            across = L0.same_bits(a['f_' + key], a1['f_' + key]) and L0.same_bits(a['J_' + key], a1['J_' + key])
            l36[key] = {'same_process': same_proc, 'threads_1_vs_default': across}
            fail = gate(fail, same_proc and across, key, '%s', '같은 프로세스 %s · 스레드 %d 대 %d %s' % (
                '같음' if same_proc else '다름', I['on_t1']['threads'], I['on']['threads'], '같음' if across else '다름'), '비트까지 같음')

    print(' L3-7 굽힘과 함께 (all 재질 · 굽힘 gn_state: 인장 · 전단 켬 대 끔)')
    l37 = {}
    for probe in PROBES:
        key = 'all_%s' % probe
        rf = L0.rel(A['gn_on']['f_' + key], A['gn_off']['f_' + key])
        rJ = L0.rel(A['gn_on']['J_' + key], A['gn_off']['J_' + key])
        l37[key] = {'f': rf, 'J': rJ}
        fail = gate(fail, rf <= LIM and rJ <= LIM, key, 'f rel %.2e · J rel %.2e', (rf, rJ), '%.0e' % LIM)

    print(' L3-8 브랜치의 기본 (인장 · 전단 환경 변수 없음 대 on, 굽힘 off, 비트)')
    keys = sorted(A['on'].keys())
    differ = [k for k in keys if k not in A['unset'] or not L0.same_bits(A['unset'][k], A['on'][k])]
    l38 = {'arrays': len(keys), 'differ': differ}
    fail = gate(fail, not differ and set(A['unset'].keys()) == set(keys), '없음/on', '%s', '배열 %d 개(f · J · 곱 · 한 단계 · 매단 천) 중 다른 것 %d' % (
        len(keys), len(differ)), '비트까지 같음')

    # 보고 (판정 밖)
    rep = {'membrane_report': {m: I[m]['membrane_report'] for m in ('unset', 'off', 'on', 'on_t1', 'gn_off', 'gn_on')},
           'rest_frame_alpha_vs_spec_a_max_rel': frame_diff, 'threads': {m: I[m]['threads'] for m in I},
           'on_vs_off_bits_differ': {k: not v['J_same_bits'] for k, v in l33.items()},
           'hang_error': {m: I[m]['hang_error'] for m in ('off', 'on')}, 'hang': {}}
    for k in ('hang', 'hangx', 'hangx20'):
        on, off = A['on'][k], A['off'][k]
        r = {'steps_on': int(on.shape[0]), 'steps_off': int(off.shape[0]), 'finite': bool(np.isfinite(on).all() and np.isfinite(off).all())}
        if on.shape == off.shape and on.shape[0] > 0:
            start = V if k == 'hang' else X
            r['end_max_diff_cm'] = float(np.abs(on[-1] - off[-1]).max())
            r['motion_cm'] = float(np.abs(off[-1] - start).max())
        rep['hang'][k] = r
    print(' 보고')
    mr = rep['membrane_report']['on']
    print('  켠 프로세스의 생성 커널 호출: %s' % (mr if mr else '보고 줄 없음'))
    print('  끈 프로세스: %s' % (rep['membrane_report']['off'],))
    print('  엔진 쉬는 틀 α 대 명세 a (uv 는 2 칸으로 넘긴다) 최대 상대 차이 %.2e' % frame_diff)
    for k, label in (('hang', 'L2 매단 천'), ('hangx', '탐침에서 시작 · 수렴 PCG'), ('hangx20', '탐침에서 시작 · 20 회 고정')):
        r = rep['hang'][k]
        print('  %-24s 켬 대 끔 끝 위치 차이 %s cm (움직인 거리 %s cm, 단계 %d · %d)' % (
            label, '%.3g' % r['end_max_diff_cm'] if 'end_max_diff_cm' in r else '-',
            '%.3g' % r['motion_cm'] if 'motion_cm' in r else '-', r['steps_on'], r['steps_off']))
    print('  스레드: %s' % rep['threads'])

    res = {
        'engine_dir': ENGINE_DIR.replace('\\', '/'), 'link_dir': LINK_DIR.replace('\\', '/'), 'base': BASE,
        'build_info_main': bi_main, 'build_info_branch': bi_br, 'branch': st, 'header_sha256': sha,
        'dll_sha256': {'main': L2.sha256_file(dll_main), 'branch': L2.sha256_file(dll_br)},
        'scene': {'nodes': int(len(V)), 'elements': int(len(F)), 'velocity_seed': VEL_SEED, 'velocity_max': VEL_MAX,
                  'hang_steps': HANG_STEPS, 'hang_from_probe_pins': list(HANGX_PINS), 'hang_from_probe_gravity_minus_y': HANGX_GRAVITY,
                  'materials': {n: m for n, m in L0.MATERIALS}},
        'limit': LIM, 'limit_step': LIM_STEP, 'min_branch': MIN_BRANCH,
        'L3_0': {'main_ok': a_ok, 'branch_ok': b_ok, 'engine_main_clean': main_status.strip() == '', 'headers_same': d_ok,
                 'e_velocity_effect': l30e, 'e_ok': e_ok, 'f_l1_ok': l1_ok, 'f_l2_ok': l2_ok, 'f_l2_branch_commit': l2_commit},
        'L3_1': l31, 'L3_2': l32, 'L3_3': l33, 'L3_4': l34, 'L3_5': l35, 'L3_6': l36, 'L3_7': l37, 'L3_8': l38,
        'report': rep,
        'gates_failed': fail,
    }
    with open(os.path.join(ROOT, 'results', 'gate_link_l3.json'), 'w', encoding='utf-8') as fh:
        json.dump(res, fh, indent=2, ensure_ascii=False)
    print('\n엔진 연결 L3 %s (실패 %d) → results/gate_link_l3.json' % ('PASS' if fail == 0 else 'FAIL', fail))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
