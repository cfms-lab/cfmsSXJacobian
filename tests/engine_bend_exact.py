"""엔진 브랜치의 실험용 굽힘 모드(exact_state · psd_state) 확인 (docs/GATE_JACOBIAN_DAMPING.md 「논문용 측정 3 — 커널 확인」).

사용:  python tests/engine_bend_exact.py <기준 브랜치 DLL> <실험 DLL> <link_l0_ref.exe>     (리포 루트에서, run_gate 밖 — 약 2 분)
기준 DLL = 실험 모드를 넣기 전의 브랜치 빌드(%TEMP%\\cfmsDrape_build_sxjlink), 실험 DLL = 넣은 뒤의 빌드(%TEMP%\\cfmsDrape_build_sxjexact).
장면 · 기준 계산은 tests/link_l0.py, 회는 tests/link_l2.py 의 --measure 를 그대로 부른다(모드는 환경 변수 CFMS_CPU_SXJ_BEND).
결과 results/gate_engine_bend_exact.json. 통과 0, 실패 1, 준비 안 됨 2.
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
OUT = os.path.join(ROOT, 'results', 'gate_engine_bend_exact.json')
COMMIT = '441b685'
OLD = ('off', 'gn_rest', 'gn_state')
NEW = ('exact_state', 'psd_state')
CHANGED = {'drape/src/sde_core/SDECore_SXJBend.cpp'}
# 커밋 뒤에는 모드 번호 주석(SDECore.h)과 그 리포의 기록(CLAUDE.md)이 함께 들어간다 — 코드 변화는 CHANGED 뿐이다
CHANGED_COMMITTED = CHANGED | {'drape/src/sde_core/SDECore.h', 'CLAUDE.md'}
LIM, LIM_PSD, LIM_EIG = 1e-12, 1e-10, -1e-12


def element_refs(ref_exe, H, w, X, k):
    """힌지마다 −k(∇C∇Cᵀ + C∇²C) 와 그 준정부호 투영을 전역으로 조립한다. 부정치 힌지 수도 센다."""
    rows = [np.concatenate([X[s].ravel(), [ev]]) for s, ev in zip(H, w)]
    out = L0.run_ref(ref_exe, 'bridson', rows)
    n3 = X.size
    Jex, Jpsd, indef = np.zeros((n3, n3)), np.zeros((n3, n3)), 0
    for s, row in zip(H, out):
        idx = np.array([3 * v + c for v in s for c in range(3)])
        C, g, Hf = row[0], row[1:13], L0.upper_to_full(row[13:], 12)
        Je = -k * (np.outer(g, g) + C * Hf)
        lam, U = np.linalg.eigh(0.5 * (Je + Je.T))
        indef += int(lam.max() > 1e-12 * np.abs(lam).max())
        Jp = (U * np.minimum(lam, 0.0)) @ U.T
        ix = np.ix_(idx, idx)
        Jex[ix] += Je
        Jpsd[ix] += Jp
    return Jex, Jpsd, indef


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
    print('실험용 굽힘 모드 확인 — 노드 %d, 힌지 %d\n  기준 DLL %s\n  실험 DLL %s' % (len(V), len(H), dll_base, dll_new))
    base = {m: L2.run_measure(dll_base, 'full', m) for m in OLD}
    new = {m: L2.run_measure(dll_new, 'full', m) for m in OLD + NEW}
    names = [n for n, _ in L0.MATERIALS]

    fail = 0
    print(' E0 전제')
    bi_b, bi_n = base['off'][1]['build_info'], new['off'][1]['build_info']
    st = L2.branch_state()
    status = L2.git('status', '--porcelain', cwd=L2.LINK_DIR).stdout.split('\n')
    dirty_files = {ln[3:].strip() for ln in status if ln.strip()}
    fail = gate(fail, str(bi_b.get('git_commit', '')).startswith(COMMIT) and bi_b.get('dirty') is False, '(a)',
                '기준 DLL git_commit %s · dirty %s' % (bi_b.get('git_commit'), bi_b.get('dirty')), '%s · dirty false' % COMMIT)
    committed = not st['head'].startswith(COMMIT[:7])
    if committed:
        # 커밋 뒤(판정문 「커밋 뒤」): 실험 DLL = 워크트리 HEAD · dirty 아님, HEAD 가 기준 커밋의 자손, 기준 커밋에서 바뀐 파일이 허용 목록 안
        desc = L2.git('merge-base', '--is-ancestor', COMMIT, 'HEAD', cwd=L2.LINK_DIR).returncode == 0
        changed = set(L2.git('diff', '--name-only', COMMIT, 'HEAD', cwd=L2.LINK_DIR).stdout.split())
        fail = gate(fail, str(bi_n.get('git_commit', '')).startswith(st['head'][:7]) and bi_n.get('dirty') is False and desc, '(b)',
                    '실험 DLL git_commit %s · dirty %s · 워크트리 HEAD %s · %s 의 자손 %s' % (
                        bi_n.get('git_commit'), bi_n.get('dirty'), st['head'], COMMIT, desc), 'HEAD 일치 · dirty false · 자손')
        fail = gate(fail, changed <= CHANGED_COMMITTED and CHANGED <= changed and not dirty_files, '(c)',
                    '%s 에서 바뀐 파일 %s · 미커밋 %s' % (COMMIT, sorted(changed), sorted(dirty_files) or '없음'),
                    '%s 안 · 미커밋 없음' % sorted(CHANGED_COMMITTED))
        dirty_files = changed
    else:
        fail = gate(fail, str(bi_n.get('git_commit', '')).startswith(COMMIT) and st['head'].startswith(COMMIT[:7]), '(b)',
                    '실험 DLL git_commit %s · dirty %s · 워크트리 HEAD %s' % (bi_n.get('git_commit'), bi_n.get('dirty'), st['head']), COMMIT)
        fail = gate(fail, dirty_files == CHANGED, '(c)', '워크트리에서 바뀐 파일 %s' % sorted(dirty_files), sorted(CHANGED))

    print(' E1 기존 모드 그대로 (두 DLL, 비트)')
    e1 = {}
    for m in OLD:
        a, b = base[m][0], new[m][0]
        same = all(L0.same_bits(a['%s_%s' % (k, n)], b['%s_%s' % (k, n)]) for n in names for k in ('f', 'J', 'P'))
        hang = L0.same_bits(a['hang'], b['hang']) and a['hang'].shape[0] == L2.HANG_STEPS
        e1[m] = {'f_J_P': bool(same), 'hang': bool(hang)}
        fail = gate(fail, same and hang, m, 'f · J · 무작위 곱(재질 %d) %s · 매단 천 %d 단계 %s' % (
            len(names), '같음' if same else '다름', L2.HANG_STEPS, '같음' if hang else '다름'), '비트까지 같음')

    print(' E2 힘 (실험 모드 대 같은 DLL 의 gn_state, 비트)')
    gs = new['gn_state'][0]
    e2 = {}
    for m in NEW:
        a = new[m][0]
        same = all(L0.same_bits(a['f_' + n], gs['f_' + n]) for n in ('bend', 'all'))
        e2[m] = bool(same)
        fail = gate(fail, same, m, '굽힘 재질 · 모든 재질의 힘 %s' % ('같음' if same else '다름'), '비트까지 같음')

    w_state = np.array([L0.bridson_w(X, uv, h) for h in H])
    Jex, Jpsd, indef = element_refs(ref_exe, H, w_state, X, L0.K_BEND)

    print(' E3 exact_state 의 J (굽힘 재질)')
    ae = new['exact_state'][0]
    e3 = {'J_vs_ref': L0.rel(ae['J_bend'], Jex),
          'product_consistency': max(L0.rel(ae['J_bend'] @ r.ravel(), p) for r, p in zip(R, ae['P_bend']))}
    fail = gate(fail, e3['J_vs_ref'] <= LIM, 'J', '대 Σ −k(∇C∇Cᵀ + C∇²C) rel %.2e' % e3['J_vs_ref'], '%.0e' % LIM)
    fail = gate(fail, e3['product_consistency'] <= LIM, '곱', '일관성 rel max %.2e' % e3['product_consistency'], '%.0e' % LIM)

    print(' E4 psd_state 의 J (굽힘 재질)')
    ap = new['psd_state'][0]
    e4 = {'J_vs_ref': L0.rel(ap['J_bend'], Jpsd), 'min_eig_ratio': eig_ratio(ap['J_bend']),
          'product_consistency': max(L0.rel(ap['J_bend'] @ r.ravel(), p) for r, p in zip(R, ap['P_bend']))}
    fail = gate(fail, e4['J_vs_ref'] <= LIM_PSD, 'J', '대 힌지마다 자른 기준 rel %.2e' % e4['J_vs_ref'], '%.0e' % LIM_PSD)
    fail = gate(fail, e4['min_eig_ratio'] >= LIM_EIG, '고윳값', '−J 의 최소 / 최대 %.2e' % e4['min_eig_ratio'], '≥ %.0e' % LIM_EIG)
    fail = gate(fail, e4['product_consistency'] <= LIM, '곱', '일관성 rel max %.2e' % e4['product_consistency'], '%.0e' % LIM)

    print(' E5 합성 (모든 재질: J_all − J_bend · f_all − f_bend 가 off 와 같은가)')
    ao = new['off'][0]
    Doff_J, Doff_f = ao['J_all'] - ao['J_bend'], ao['f_all'] - ao['f_bend']
    sJ, sf = float(np.abs(ao['J_all']).max()), float(np.abs(ao['f_all']).max())
    e5 = {}
    for m in NEW:
        a = new[m][0]
        eJ = float(np.abs((a['J_all'] - a['J_bend']) - Doff_J).max()) / sJ
        ef = float(np.abs((a['f_all'] - a['f_bend']) - Doff_f).max()) / sf
        e5[m] = {'J': eJ, 'f': ef}
        fail = gate(fail, eJ <= LIM and ef <= LIM, m, 'J %.2e · f %.2e' % (eJ, ef), '%.0e' % LIM)

    rep = {'hinges': int(len(H)), 'indefinite_exact_element': int(indef), 'exact_global_min_eig_ratio': eig_ratio(ae['J_bend']),
           'gn_global_min_eig_ratio': eig_ratio(gs['J_bend']), 'psd_vs_exact_rel': L0.rel(ap['J_bend'], ae['J_bend']),
           'psd_vs_gn_rel': L0.rel(ap['J_bend'], gs['J_bend']),
           'J_vs_fd_of_own_force': {m: L0.rel(new[m][0]['J_bend'], new[m][0]['Jfd_bend']) for m in ('off', 'gn_state') + NEW},
           'hang': {m: {'steps_done': int(new[m][0]['hang'].shape[0]), 'error': new[m][1].get('hang_error', ''),
                        'finite': bool(np.isfinite(new[m][0]['hang']).all())} for m in NEW}}
    print(' 보고')
    print('  정확한 요소 행렬이 부정치인 힌지 %d / %d · 전역 −J 의 최소 / 최대 고윳값: exact %.3g · psd %.3g · GN %.3g' % (
        indef, len(H), rep['exact_global_min_eig_ratio'], e4['min_eig_ratio'], rep['gn_global_min_eig_ratio']))
    print('  psd 대 exact %.3g · psd 대 GN %.3g' % (rep['psd_vs_exact_rel'], rep['psd_vs_gn_rel']))
    print('  J 대 제 힘의 유한차분: ' + ' · '.join('%s %.3g' % (m, v) for m, v in rep['J_vs_fd_of_own_force'].items()))
    print('  매단 천: ' + ' · '.join('%s %d 단계 유한 %s %s' % (m, r['steps_done'], r['finite'], r['error'][:60]) for m, r in rep['hang'].items()))

    res = {'dll': {'base': dll_base, 'new': dll_new}, 'build_info': {'base': bi_b, 'new': bi_n}, 'branch': st,
           'dirty_files': sorted(dirty_files), 'limits': {'J': LIM, 'J_psd': LIM_PSD, 'min_eig_ratio': LIM_EIG},
           'E1': e1, 'E2': e2, 'E3': e3, 'E4': e4, 'E5': e5, 'report': rep, 'gates_failed': fail}
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
        fh.write('\n')
    print('\n실험용 굽힘 모드 확인 %s (실패 %d) → %s' % ('PASS' if fail == 0 else 'FAIL', fail, os.path.relpath(OUT, ROOT)))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
