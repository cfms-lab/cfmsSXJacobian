"""엔진 브랜치의 실험용 굽힘 모드 symgn_state 확인 (docs/GATE_JACOBIAN_DAMPING.md 「논문용 측정 4 — 커널 확인」).

사용:  python tests/engine_bend_symgn.py <기준 DLL> <실험 DLL> <link_l0_ref.exe>     (리포 루트에서, run_gate 밖 — 약 2 분)
기준 DLL = 엔진 브랜치 277425b 의 빌드(%TEMP%\\cfmsDrape_build_sxjexact), 실험 DLL = symgn_state 를 넣은 빌드(%TEMP%\\cfmsDrape_build_sxjsymgn).
symgn_state = gn_state 의 노드쌍 블록을 블록마다 대칭화한 것(블록 대칭 Gauss–Newton). 장면 · 기준 계산은 tests/link_l0.py, 회는 tests/link_l2.py.
결과 results/gate_engine_bend_symgn.json. 통과 0, 실패 1, 준비 안 됨 2.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import link_l0 as L0  # noqa: E402
import link_l2 as L2  # noqa: E402

ROOT = L0.ROOT
OUT = os.path.join(ROOT, 'results', 'gate_engine_bend_symgn.json')
COMMIT = '277425b'
OLD = ('off', 'gn_rest', 'gn_state', 'exact_state', 'psd_state')
NEW = 'symgn_state'
ALLOWED = {'drape/src/sde_core/SDECore_SXJBend.cpp', 'drape/src/sde_core/SDECore.h', 'CLAUDE.md'}
LIM = 1e-12


def symgn_ref(ref_exe, H, w, X, k):
    """힌지마다 GN 블록 −k ∇_aC ∇_bCᵀ 를 블록마다 대칭화하고 대각을 −Σ 로 두어 전역으로 조립한다. GN 그대로의 조립도 함께 준다."""
    rows = [np.concatenate([X[s].ravel(), [ev]]) for s, ev in zip(H, w)]
    out = L0.run_ref(ref_exe, 'bridson', rows)
    n3 = X.size
    Jsym, Jgn = np.zeros((n3, n3)), np.zeros((n3, n3))
    for s, row in zip(H, out):
        idx = np.array([3 * v + c for v in s for c in range(3)])
        g = row[1:13]
        Je = -k * np.outer(g, g)
        Js = np.zeros((12, 12))
        for a in range(4):
            for b in range(4):
                if a != b:
                    B = Je[3 * a:3 * a + 3, 3 * b:3 * b + 3]
                    S = 0.5 * (B + B.T)
                    Js[3 * a:3 * a + 3, 3 * b:3 * b + 3] = S
                    Js[3 * a:3 * a + 3, 3 * a:3 * a + 3] -= S
        ix = np.ix_(idx, idx)
        Jsym[ix] += Js
        Jgn[ix] += Je
    return Jsym, Jgn


def eig_ratio(J):
    """−J 의 가장 작은 고윳값 / 가장 큰 고윳값."""
    lam = np.linalg.eigvalsh(-0.5 * (J + J.T))
    return float(lam.min() / lam.max())


def gate(fail, ok, label, text, limit):
    print('  %-12s %s  %s  (기준 %s)' % (label, 'PASS' if ok else 'FAIL', text, limit))
    return fail + (0 if ok else 1)


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        return 2
    dll_base, dll_new, ref_exe = (os.path.abspath(a) for a in sys.argv[1:4])
    missing = [p for p in (dll_base, dll_new, ref_exe) if not os.path.isfile(p)]
    if missing:
        print('없는 것: %s' % ', '.join(missing))
        return 2
    V, F, uv, X, R = L0.scene()
    H = L0.hinges(F)
    print('실험용 굽힘 모드 symgn_state 확인 — 노드 %d, 힌지 %d\n  기준 DLL %s\n  실험 DLL %s' % (len(V), len(H), dll_base, dll_new))
    base = {m: L2.run_measure(dll_base, 'full', m) for m in OLD}
    new = {m: L2.run_measure(dll_new, 'full', m) for m in OLD + (NEW,)}
    names = [n for n, _ in L0.MATERIALS]

    fail = 0
    print(' G0 전제')
    bi_b, bi_n = base['off'][1]['build_info'], new['off'][1]['build_info']
    st = L2.branch_state()
    status = L2.git('status', '--porcelain', cwd=L2.LINK_DIR).stdout.split('\n')
    dirty_files = {ln[3:].strip() for ln in status if ln.strip()}
    fail = gate(fail, str(bi_b.get('git_commit', '')).startswith(COMMIT) and bi_b.get('dirty') is False, '(a)',
                '기준 DLL git_commit %s · dirty %s' % (bi_b.get('git_commit'), bi_b.get('dirty')), '%s · dirty false' % COMMIT)
    fail = gate(fail, str(bi_n.get('git_commit', '')).startswith(COMMIT) and st['head'].startswith(COMMIT[:7]), '(b)',
                '실험 DLL git_commit %s · dirty %s · 워크트리 HEAD %s' % (bi_n.get('git_commit'), bi_n.get('dirty'), st['head']), COMMIT)
    fail = gate(fail, dirty_files <= ALLOWED and 'drape/src/sde_core/SDECore_SXJBend.cpp' in dirty_files, '(c)',
                '워크트리에서 바뀐 파일 %s' % sorted(dirty_files), '%s 안' % sorted(ALLOWED))

    print(' G1 기존 모드 그대로 (두 DLL, 비트)')
    g1 = {}
    for m in OLD:
        a, b = base[m][0], new[m][0]
        same = all(L0.same_bits(a['%s_%s' % (k, n)], b['%s_%s' % (k, n)]) for n in names for k in ('f', 'J', 'P'))
        hang = L0.same_bits(a['hang'], b['hang']) and a['hang'].shape[0] == L2.HANG_STEPS
        g1[m] = {'f_J_P': bool(same), 'hang': bool(hang)}
        fail = gate(fail, same and hang, m, 'f · J · 무작위 곱(재질 %d) %s · 매단 천 %d 단계 %s' % (
            len(names), '같음' if same else '다름', L2.HANG_STEPS, '같음' if hang else '다름'), '비트까지 같음')

    print(' G2 힘 (symgn_state 대 같은 DLL 의 gn_state, 비트)')
    gs, an = new['gn_state'][0], new[NEW][0]
    g2 = all(L0.same_bits(an['f_' + n], gs['f_' + n]) for n in ('bend', 'all'))
    fail = gate(fail, g2, NEW, '굽힘 재질 · 모든 재질의 힘 %s' % ('같음' if g2 else '다름'), '비트까지 같음')

    w_state = np.array([L0.bridson_w(X, uv, h) for h in H])
    Jsym, Jgn = symgn_ref(ref_exe, H, w_state, X, L0.K_BEND)

    print(' G3 symgn_state 의 J (굽힘 재질)')
    g3 = {'J_vs_ref': L0.rel(an['J_bend'], Jsym),
          'product_consistency': max(L0.rel(an['J_bend'] @ r.ravel(), p) for r, p in zip(R, an['P_bend'])),
          'gn_state_vs_gn_ref': L0.rel(gs['J_bend'], Jgn)}
    fail = gate(fail, g3['J_vs_ref'] <= LIM, 'J', '대 블록마다 대칭화한 Σ −k ∇C∇Cᵀ rel %.2e' % g3['J_vs_ref'], '%.0e' % LIM)
    fail = gate(fail, g3['product_consistency'] <= LIM, '곱', '일관성 rel max %.2e' % g3['product_consistency'], '%.0e' % LIM)
    fail = gate(fail, g3['gn_state_vs_gn_ref'] <= LIM, '(gn_state)', '같은 기준 계산의 GN 조립과 gn_state rel %.2e' % g3['gn_state_vs_gn_ref'], '%.0e' % LIM)

    print(' G4 합성 (모든 재질: J_all − J_bend · f_all − f_bend 가 off 와 같은가)')
    ao = new['off'][0]
    Doff_J, Doff_f = ao['J_all'] - ao['J_bend'], ao['f_all'] - ao['f_bend']
    sJ, sf = float(np.abs(ao['J_all']).max()), float(np.abs(ao['f_all']).max())
    g4 = {'J': float(np.abs((an['J_all'] - an['J_bend']) - Doff_J).max()) / sJ, 'f': float(np.abs((an['f_all'] - an['f_bend']) - Doff_f).max()) / sf}
    fail = gate(fail, g4['J'] <= LIM and g4['f'] <= LIM, NEW, 'J %.2e · f %.2e' % (g4['J'], g4['f']), '%.0e' % LIM)

    asym = max(float(np.abs(an['J_bend'][3 * a:3 * a + 3, 3 * b:3 * b + 3] - an['J_bend'][3 * a:3 * a + 3, 3 * b:3 * b + 3].T).max())
               for a in range(len(V)) for b in range(len(V)) if a != b) / float(np.abs(an['J_bend']).max())
    rep = {'hinges': int(len(H)), 'symgn_vs_gn_rel': L0.rel(an['J_bend'], gs['J_bend']), 'offdiag_block_asymmetry': asym,
           'min_eig_ratio': {m: eig_ratio(new[m][0]['J_bend']) for m in ('off', 'gn_state', NEW)},
           'J_vs_fd_of_own_force': {m: L0.rel(new[m][0]['J_bend'], new[m][0]['Jfd_bend']) for m in ('off', 'gn_state', NEW)},
           'hang': {'steps_done': int(an['hang'].shape[0]), 'error': new[NEW][1].get('hang_error', ''), 'finite': bool(np.isfinite(an['hang']).all())}}
    print(' 보고 (변형된 시험 메시)')
    print('  symgn_state 대 gn_state 의 J 상대 차 %.3g · 비대각 블록의 비대칭 %.2e' % (rep['symgn_vs_gn_rel'], asym))
    print('  전역 −J 의 최소 / 최대 고윳값: ' + ' · '.join('%s %.3g' % (m, v) for m, v in rep['min_eig_ratio'].items()))
    print('  J 대 제 힘의 유한차분: ' + ' · '.join('%s %.3g' % (m, v) for m, v in rep['J_vs_fd_of_own_force'].items()))
    print('  매단 천: %d 단계 유한 %s %s' % (rep['hang']['steps_done'], rep['hang']['finite'], rep['hang']['error'][:60]))

    res = {'dll': {'base': dll_base, 'new': dll_new}, 'build_info': {'base': bi_b, 'new': bi_n}, 'branch': st,
           'dirty_files': sorted(dirty_files), 'limits': {'J': LIM}, 'G1': g1, 'G2': bool(g2), 'G3': g3, 'G4': g4, 'report': rep,
           'gates_failed': fail}
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
        fh.write('\n')
    print('\nsymgn_state 확인 %s (실패 %d) → %s' % ('PASS' if fail == 0 else 'FAIL', fail, os.path.relpath(OUT, ROOT)))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
