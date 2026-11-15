"""엔진 연결 L2 (docs/GATE_LINK_L2.md): 엔진 브랜치 sxjacobian-link 의 생성 Bridson GN 갈래를 디버그 고리로 잰다.

사용:  python tests/link_l2.py <main cfmsDrape.dll> <브랜치 cfmsDrape.dll> <link_l0_ref.exe>
       (내부) python tests/link_l2.py --measure <dll> <full|light> <out.npz>
              — 한 DLL 을 다른 프로세스에서 잰다. 모드는 그 프로세스의 환경 변수 CFMS_CPU_SXJ_BEND 로 준다.
엔진 리포(main)는 CFMS_DRAPE_DIR(없으면 ..\\cfmsDrape_dev), 브랜치 워크트리는 CFMS_DRAPE_LINK_DIR
(없으면 ..\\cfmsDrape_dev_sxjacobian-link). Python 포장 cfmsdrape 는 main 의 것을 읽는다(브랜치는 포장을 바꾸지 않는다).
장면 · 기준 계산은 tests/link_l0.py 와 같다. 리포 루트에서 실행. 결과는 results/gate_link_l2.json. 통과 0, 실패 1, 준비 안 됨 2.
"""
import sys
sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import link_l0 as L0  # noqa: E402

ROOT = L0.ROOT
ENGINE_DIR = L0.ENGINE_DIR
LINK_DIR = os.path.abspath(os.environ.get('CFMS_DRAPE_LINK_DIR')
                           or os.path.join(ROOT, '..', 'cfmsDrape_dev_sxjacobian-link'))
BASE = L0.ENGINE_COMMIT                                  # 378d20c
HEADER_GEN = os.path.join(ROOT, 'generated', 'BridsonBend_Engine.h')
HEADER_ENGINE = os.path.join(LINK_DIR, 'drape', 'src', 'sde_core', 'sxgen', 'BridsonBend_Engine.h')
ALLOWED = {'drape/CMakeLists.txt', 'drape/src/sde_core/SDECore.h', 'drape/src/sde_core/SDECore.cpp',
           'drape/src/sde_core/SDECore_SXJBend.cpp', 'drape/src/sde_core/sxgen/BridsonBend_Engine.h',
           # 2026-11-13 더함 — docs/GATE_LINK_L4.md 「2 차 (다)」(켬 경로의 전치 처리 곱)
           'drape/src/sde_core/SDESparseMatrix.h', 'drape/src/sde_core/SDESparseMatrix.cpp',
           # 2026-11-13 더함 — docs/GATE_LINK_L3.md (생성 인장 · Provot 전단 갈래)
           'drape/src/sde_core/SDECore_SXJMembrane.cpp', 'drape/src/sde_core/sxgen/BaraffStretch_Engine.h',
           'drape/src/sde_core/sxgen/ProvotShear_Engine.h',
           # 2026-11-13 더함 — 엔진 리포 규칙(기본값 전환은 같은 커밋에 CLAUDE.md 에 적는다)에 따른 기록. 코드가 아니다.
           'CLAUDE.md',
           # 2026-11-14 더함 — docs/GATE_ENGINE_SUITE.md (막 비선형 시험 기준 2 를 상대 1e-12 로). 시험 스크립트지 엔진 코드가 아니다.
           'drape/scripts/py_membrane_nonlinear_test.py'}
MODES = ('gn_rest', 'gn_state')
LIM = 1e-12
HANG_STEPS = 40
HANG_PINS = (0, 4)
HANG_GRAVITY = 980.0
LIGHT_MATERIALS = ('bend', 'all')


def hang_run(Drape, dll, V, F, uv, mat):
    d = Drape(dll)
    try:
        d.set_mesh(V, F)
        d.set_mesh_2d(uv)
        d.set_material(mass=0.1, **mat)
        d.set_gravity((0.0, 0.0, -HANG_GRAVITY))
        d.set_timestep(0.005)
        for p in HANG_PINS:
            d.pin_vertex(p)
        traj, err = [], ''
        for _ in range(HANG_STEPS):
            try:
                d.step(1)
            except Exception as e:   # StepDiverged 등 — 판정은 유한성으로
                err = '%s: %s' % (type(e).__name__, e)
                break
            traj.append(np.array(d.positions(), dtype=float).copy())
        return np.array(traj), err
    finally:
        d.close()


def measure(dll, kind, out):
    Drape = L0.load_drape()
    V, F, uv, X, R = L0.scene()
    arrs, info = {}, {}
    names = [n for n, _ in L0.MATERIALS] if kind == 'full' else list(LIGHT_MATERIALS)
    for name, mat in L0.MATERIALS:
        if name not in names:
            continue
        m = L0.engine_measure(Drape, dll, V, F, uv, X, R, mat, full=(kind == 'full'))
        for k in ('f', 'J', 'f2', 'P', 'P2', 'Jfd'):
            if k in m:
                arrs['%s_%s' % (k, name)] = m[k]
        info['threads'] = m['threads']
        info['build_info'] = m['build_info']
    if kind == 'full':
        traj, err = hang_run(Drape, dll, V, F, uv, dict(L0.MATERIALS)['all'])
        arrs['hang'] = traj
        info['hang_error'] = err
    np.savez(out, **arrs)
    with open(out + '.json', 'w', encoding='utf-8') as fh:
        json.dump(info, fh, ensure_ascii=False)
    return 0


def run_measure(dll, kind, mode, threads1=False):
    """다른 프로세스에서 잰다. mode None = 환경 변수 없음(main DLL, 그리고 브랜치의 기본 = gn_state), 'off' = 끔."""
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    env.pop('CFMS_CPU_SXJ_BEND', None)
    if mode is not None:
        env['CFMS_CPU_SXJ_BEND'] = mode
    # 인장 · 전단은 엔진 손 커널로 고정한다 — 브랜치의 기본이 on 이 된 뒤(2026-11-13, docs/GATE_LINK_L3.md 끝 절)로는
    # 없음 ≠ 끔이고, L2 는 굽힘만 본다. main DLL 은 이 변수를 모른다.
    env['CFMS_CPU_SXJ_MEMBRANE'] = 'off'
    if threads1:
        env['OMP_NUM_THREADS'] = '1'
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, 'm.npz')
        r = subprocess.run([sys.executable, os.path.abspath(__file__), '--measure', dll, kind, out],
                           env=env, cwd=ROOT, capture_output=True, text=True, encoding='utf-8', errors='replace')
        if r.returncode != 0:
            raise RuntimeError('measure 실패 (%s, %s): %s' % (mode, kind, r.stderr[-2000:]))
        arrs = dict(np.load(out))
        with open(out + '.json', encoding='utf-8') as fh:
            info = json.load(fh)
        info['stderr'] = r.stderr.strip()
        return arrs, info


def git(*args, cwd):
    return subprocess.run(['git', '-C', cwd] + list(args), capture_output=True, text=True, encoding='utf-8')


def sha256_file(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def branch_state():
    head = git('rev-parse', '--short', 'HEAD', cwd=LINK_DIR).stdout.strip()
    is_desc = git('merge-base', '--is-ancestor', BASE, 'HEAD', cwd=LINK_DIR).returncode == 0
    changed = set(git('diff', '--name-only', BASE, cwd=LINK_DIR).stdout.split())
    untracked = sorted(git('ls-files', '--others', '--exclude-standard', cwd=LINK_DIR).stdout.split())
    changed |= set(untracked)
    status = git('status', '--porcelain', cwd=LINK_DIR).stdout
    h = hashlib.sha256(git('diff', '--binary', BASE, cwd=LINK_DIR).stdout.encode('utf-8'))
    for p in untracked:
        h.update(p.encode('utf-8'))
        h.update(open(os.path.join(LINK_DIR, p), 'rb').read())
    return {'head': head, 'descends_from_base': is_desc, 'changed_files': sorted(changed),
            'outside_allowed': sorted(changed - ALLOWED), 'uncommitted': status.strip() != '',
            'diff_sha256': h.hexdigest()}


def gate(fail, ok, label, fmt, value, limit):
    print('  %-6s %s  %s  (기준 %s)' % (label, 'PASS' if ok else 'FAIL', fmt % value, limit))
    return fail + (0 if ok else 1)


def main():
    if len(sys.argv) == 5 and sys.argv[1] == '--measure':
        return measure(sys.argv[2], sys.argv[3], sys.argv[4])
    if len(sys.argv) != 4:
        print(__doc__)
        return 2
    dll_main, dll_br, ref_exe = (os.path.abspath(a) for a in sys.argv[1:4])
    missing = [p for p in (dll_main, dll_br, ref_exe, HEADER_GEN) if not os.path.isfile(p)]
    missing += [p for p in (os.path.join(ENGINE_DIR, 'drape'), os.path.join(LINK_DIR, 'drape')) if not os.path.isdir(p)]
    if missing:
        print('없는 것: %s' % ', '.join(missing))
        return 2

    V, F, uv, X, R = L0.scene()
    H = L0.hinges(F)
    print('엔진 연결 L2 — 노드 %d, 힌지 %d\n  main DLL %s\n  브랜치 DLL %s (워크트리 %s)' % (len(V), len(H), dll_main, dll_br, LINK_DIR))

    # 끔은 'off' 로 명시한다. 환경 변수 없음은 브랜치의 기본(gn_state, 사용자 결정 2026-11-13)이다 — L2-6 이 견준다.
    meas = {'main': run_measure(dll_main, 'full', None), 'off': run_measure(dll_br, 'full', 'off'),
            'default': run_measure(dll_br, 'full', None)}
    for mode in MODES:
        meas[mode] = run_measure(dll_br, 'full', mode)
        meas[mode + '_t1'] = run_measure(dll_br, 'light', mode, threads1=True)

    # 이 리포 기준 (k = BE_U = 1)
    U3 = np.concatenate([uv, np.zeros((len(uv), 1))], axis=1)
    w_state = np.array([L0.bridson_w(X, uv, h) for h in H])
    w_rest = np.array([L0.bridson_w(U3, uv, h) for h in H])
    ref = {}
    for mode, w in (('gn_state', w_state), ('gn_rest', w_rest)):
        f_ref, Jgg, _, _, _ = L0.reference(ref_exe, 'bridson', H, w, X, L0.K_BEND)
        ref[mode] = (f_ref, Jgg)

    bi_main, bi_br = meas['main'][1]['build_info'], meas['off'][1]['build_info']
    st = branch_state()
    main_status = git('status', '--porcelain', cwd=ENGINE_DIR).stdout
    sha_gen, sha_eng = sha256_file(HEADER_GEN), (sha256_file(HEADER_ENGINE) if os.path.isfile(HEADER_ENGINE) else '')

    fail = 0
    print(' L2-0 전제')
    a_ok = str(bi_main.get('git_commit', '')).startswith(BASE) and bi_main.get('dirty') is False
    fail = gate(fail, a_ok, '(a)', '%s', 'main DLL git_commit %s · dirty %s' % (bi_main.get('git_commit'), bi_main.get('dirty')), '%s · dirty false' % BASE)
    b_ok = (st['head'] != '' and str(bi_br.get('git_commit', '')).startswith(st['head'][:7]) and st['descends_from_base']
            and not st['outside_allowed'])
    fail = gate(fail, b_ok, '(b)', '%s', '브랜치 DLL git_commit %s · 워크트리 HEAD %s · %s 자손 %s · 허용 밖 %s · dirty %s' % (
        bi_br.get('git_commit'), st['head'], BASE, st['descends_from_base'], st['outside_allowed'] or '없음', bi_br.get('dirty')),
        'HEAD 일치 · 자손 · 허용 목록 안')
    fail = gate(fail, main_status.strip() == '', '(c)', '%s', '엔진 main git status %s' % ('비어 있음' if not main_status.strip() else '줄 %d 개' % len(main_status.splitlines())), '비어 있음')
    fail = gate(fail, sha_gen == sha_eng, '(d)', '%s', '헤더 사본 %s' % ('같음' if sha_gen == sha_eng else '다름 %s / %s' % (sha_eng[:16], sha_gen[:16])), 'SHA-256 같음')

    print(' L2-1 끄면 main 그대로 (브랜치 CFMS_CPU_SXJ_BEND=off 대 main, 비트)')
    am, ao = meas['main'][0], meas['off'][0]
    l21 = {}
    for name, _ in L0.MATERIALS:
        same = all(L0.same_bits(am['%s_%s' % (k, name)], ao['%s_%s' % (k, name)]) for k in ('f', 'J', 'P'))
        l21[name] = same
        fail = gate(fail, same, name, '%s', 'f · J · 무작위 곱 %s' % ('같음' if same else '다름'), '비트까지 같음')
    hang_same = L0.same_bits(am['hang'], ao['hang']) and am['hang'].shape[0] == HANG_STEPS
    l21['hang'] = hang_same
    fail = gate(fail, hang_same, '매단천', '%s', '%d 단계 위치 %s' % (HANG_STEPS, '같음' if hang_same else '다름'), '비트까지 같음')

    print(' L2-2 힘 (bend 재질, k = BE_U = 1)')
    gs, gr = meas['gn_state'][0], meas['gn_rest'][0]
    l22 = {'gn_state_vs_engine_oh': L0.rel(gs['f_bend'], am['f_bend']),
           'gn_state_vs_ref_state_w': L0.rel(gs['f_bend'], ref['gn_state'][0]),
           'gn_rest_vs_ref_rest_w': L0.rel(gr['f_bend'], ref['gn_rest'][0])}
    fail = gate(fail, l22['gn_state_vs_engine_oh'] <= LIM, 'state', 'f 대 엔진 Oh f rel %.2e', l22['gn_state_vs_engine_oh'], '%.0e' % LIM)
    fail = gate(fail, l22['gn_state_vs_ref_state_w'] <= LIM, 'state', 'f 대 기준(w = 그 상태) rel %.2e', l22['gn_state_vs_ref_state_w'], '%.0e' % LIM)
    fail = gate(fail, l22['gn_rest_vs_ref_rest_w'] <= LIM, 'rest', 'f 대 기준(w̄ = 쉬는 값) rel %.2e', l22['gn_rest_vs_ref_rest_w'], '%.0e' % LIM)

    print(' L2-3 J (bend 재질)')
    l23 = {}
    for mode in MODES:
        a = meas[mode][0]
        jr = L0.rel(a['J_bend'], ref[mode][1])
        pc = max(L0.rel(a['J_bend'] @ r.ravel(), p) for r, p in zip(R, a['P_bend']))
        l23[mode] = {'J_vs_ref_gn': jr, 'product_consistency': pc}
        fail = gate(fail, jr <= LIM, mode[3:], 'J 대 GN 기준 rel %.2e', jr, '%.0e' % LIM)
        fail = gate(fail, pc <= LIM, mode[3:], '곱 일관성 rel max %.2e', pc, '%.0e' % LIM)

    print(' L2-4 합성 (all 재질: J_all − J_bend 가 스위치 끔과 같은가)')
    l24 = {}
    Doff_J, Doff_f = ao['J_all'] - ao['J_bend'], ao['f_all'] - ao['f_bend']
    sJ, sf = float(np.abs(ao['J_all']).max()), float(np.abs(ao['f_all']).max())
    for mode in MODES:
        a = meas[mode][0]
        eJ = float(np.abs((a['J_all'] - a['J_bend']) - Doff_J).max()) / sJ
        ef = float(np.abs((a['f_all'] - a['f_bend']) - Doff_f).max()) / sf
        l24[mode] = {'J': eJ, 'f': ef}
        fail = gate(fail, eJ <= LIM and ef <= LIM, mode[3:], 'J %.2e · f %.2e', (eJ, ef), '%.0e' % LIM)

    print(' L2-5 결정성 (브랜치, 같은 프로세스 두 번 · 단일 스레드 프로세스 대 기본)')
    l25 = {}
    for mode in MODES:
        a, a1 = meas[mode][0], meas[mode + '_t1'][0]
        for name in LIGHT_MATERIALS:
            same_proc = L0.same_bits(a['f_' + name], a['f2_' + name]) and L0.same_bits(a['P_' + name], a['P2_' + name])
            across = L0.same_bits(a['f_' + name], a1['f_' + name]) and L0.same_bits(a['J_' + name], a1['J_' + name])
            l25['%s_%s' % (mode, name)] = {'same_process': same_proc, 'threads_1_vs_default': across}
            fail = gate(fail, same_proc and across, '%s/%s' % (mode[3:], name), '%s', '같은 프로세스 %s · 스레드 %d 대 %d %s' % (
                '같음' if same_proc else '다름', meas[mode + '_t1'][1]['threads'], meas[mode][1]['threads'],
                '같음' if across else '다름'), '비트까지 같음')

    print(' L2-6 브랜치의 기본 (환경 변수 없음 대 gn_state, 비트)')
    ad, ag = meas['default'][0], meas['gn_state'][0]
    l26 = {}
    for name, _ in L0.MATERIALS:
        same = all(L0.same_bits(ad['%s_%s' % (k, name)], ag['%s_%s' % (k, name)]) for k in ('f', 'J', 'P'))
        l26[name] = same
        fail = gate(fail, same, name, '%s', 'f · J · 무작위 곱 %s' % ('같음' if same else '다름'), '비트까지 같음')
    hang26 = L0.same_bits(ad['hang'], ag['hang']) and ad['hang'].shape[0] == HANG_STEPS
    l26['hang'] = hang26
    fail = gate(fail, hang26, '매단천', '%s', '%d 단계 위치 %s' % (HANG_STEPS, '같음' if hang26 else '다름'), '비트까지 같음')

    # 보고 (판정 밖)
    rep = {'gn_rest_f_vs_engine_oh': L0.rel(gr['f_bend'], am['f_bend']),
           'w_rest_over_state': {'min': float((w_rest / w_state).min()), 'median': float(np.median(w_rest / w_state)),
                                 'max': float((w_rest / w_state).max())}}
    for mode in MODES:
        a = meas[mode][0]
        rep[mode] = {'J_vs_engine_oh_J': L0.rel(a['J_bend'], am['J_bend']), 'J_vs_fd_of_own_force': L0.rel(a['J_bend'], a['Jfd_bend']),
                     'hang_steps_done': int(a['hang'].shape[0]), 'hang_error': meas[mode][1].get('hang_error', ''),
                     'hang_finite': bool(np.isfinite(a['hang']).all()),
                     'stderr': meas[mode][1].get('stderr', '')}
        if a['hang'].shape[0] == ao['hang'].shape[0]:
            d_end = float(np.abs(a['hang'][-1] - ao['hang'][-1]).max())
            drop = float(np.abs(ao['hang'][-1] - ao['hang'][0]).max())
            rep[mode]['hang_end_max_diff_cm'] = d_end
            rep[mode]['hang_end_diff_over_motion'] = d_end / drop if drop > 0 else None
    rep['engine_oh_J_vs_fd'] = L0.rel(am['J_bend'], am['Jfd_bend'])
    print(' 보고')
    print('  gn_rest f 대 엔진 Oh f %.2e · w̄ / w 최소 %.3f · 중앙 %.3f · 최대 %.3f' % (
        rep['gn_rest_f_vs_engine_oh'], rep['w_rest_over_state']['min'], rep['w_rest_over_state']['median'], rep['w_rest_over_state']['max']))
    for mode in MODES:
        r = rep[mode]
        print('  %-8s J 대 엔진 Oh J %.2e · J 대 자기 힘의 유한차분 %.2e · 매단 천 %d 단계 유한 %s · 끝 위치 차이(대 끔) %s cm' % (
            mode, r['J_vs_engine_oh_J'], r['J_vs_fd_of_own_force'], r['hang_steps_done'], r['hang_finite'],
            '%.3g' % r['hang_end_max_diff_cm'] if 'hang_end_max_diff_cm' in r else '-'))
    print('  엔진 Oh J 대 유한차분 %.2e' % rep['engine_oh_J_vs_fd'])

    res = {
        'engine_dir': ENGINE_DIR.replace('\\', '/'), 'link_dir': LINK_DIR.replace('\\', '/'), 'base': BASE,
        'build_info_main': bi_main, 'build_info_branch': bi_br, 'branch': st,
        'header_sha256': {'generated': sha_gen, 'engine_copy': sha_eng},
        'dll_sha256': {'main': sha256_file(dll_main), 'branch': sha256_file(dll_br)},
        'threads': {m: meas[m][1]['threads'] for m in meas},
        'scene': {'nodes': int(len(V)), 'hinges': int(len(H)), 'hang_steps': HANG_STEPS, 'hang_pins': list(HANG_PINS),
                  'hang_gravity': HANG_GRAVITY, 'materials': {n: m for n, m in L0.MATERIALS}},
        'limit': LIM,
        'L2_0': {'main_ok': a_ok, 'branch_ok': b_ok, 'engine_main_clean': main_status.strip() == '', 'header_same': sha_gen == sha_eng},
        'L2_1': l21, 'L2_2': l22, 'L2_3': l23, 'L2_4': l24, 'L2_5': l25, 'L2_6': l26,
        'report': rep,
        'gates_failed': fail,
    }
    with open(os.path.join(ROOT, 'results', 'gate_link_l2.json'), 'w', encoding='utf-8') as fh:
        json.dump(res, fh, indent=2, ensure_ascii=False)
    print('\n엔진 연결 L2 %s (실패 %d) → results/gate_link_l2.json' % ('PASS' if fail == 0 else 'FAIL', fail))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
