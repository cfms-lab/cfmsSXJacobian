"""논문용 계산 (docs/GATE_JACOBIAN_DAMPING.md 「모드별 기여」, 판정 아님): 자유 변의 종단 속도를 모드별 기여로 나눈다.

사용:  python tests/measure_mode_sum.py <브랜치 cfmsDrape.dll>     (리포 루트에서, run_gate 밖 — 약 20 초, 적분 없음)
엔진의 좌변 행렬을 tests/gate_jacobian_damping.py 의 「matrix」 회와 같은 방법으로 읽어(줄마다 같은 법선 변위 → (n−1)×(n−1) 축약 행렬 K_r,
줄 질량 M_r = n m I) 일반화 고유 문제 K_r u_k = λ_k M_r u_k (u_kᵀ M_r u_l = δ_kl)를 풀고, 깨끗한 설정(k_d = 0)의 종단 속도

    V_∞ = e_끝ᵀ (h K_r)⁻¹ M_r g1 = Σ_k  u_k(끝) · u_kᵀ M_r g1 / (h λ_k)

의 항을 낸다. 합이 results/measure_jacobian_drop.json 의 v_ss 와 같아야 한다(같은 행렬 · 같은 식 — 다르면 멈춘다).
g / (h λ_1) 은 중력을 모드의 구동값으로 그대로 넣은 값이다 — 첫 모드의 기여(투영 포함)와 다르다는 것을 함께 적는다.
결과 results/measure_mode_sum.json.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gate_jacobian_damping as G  # noqa: E402

OUT = os.path.join(G.ROOT, 'results', 'measure_mode_sum.json')
DROP = os.path.join(G.ROOT, 'results', 'measure_jacobian_drop.json')
N_TERMS = 6            # 적어 두는 모드 수
REL_TOL = 1e-9         # 합과 기록의 v_ss


def main():
    dll = os.path.abspath(sys.argv[1])
    if not os.path.isfile(dll):
        print('DLL 없음:', dll)
        return 2
    with open(DROP, encoding='utf-8') as fh:
        drop = json.load(fh)
    cases = [(r['n'], r['BE_U'], r['sps'], r['v_ss']) for r in drop['runs']]
    keys = sorted({(n, b) for n, b, _, _ in cases})
    td = tempfile.mkdtemp(prefix='cfmsSXJ_modesum_')
    t0 = time.perf_counter()
    try:
        mats = G.run_all(td, {k: ('off', 'matrix', dict(dll=dll, n=k[0], bend=k[1])) for k in keys})
    finally:
        shutil.rmtree(td, ignore_errors=True)
    info = mats[keys[0]]['build_info']
    rows, worst = [], 0.0
    print('모드별 기여 — 깨끗한 설정의 종단 속도 (엔진 %s, dirty %s), cm/s' % (info.get('git_commit'), info.get('dirty')))
    for n, b, sps, v_ref in cases:
        lam, U, mrow = G.eig(mats[(n, b)]['K'], n)             # U 는 M_r-정규: uᵀ M_r u = 1
        drive = U.T @ (mrow * G.GRAV * np.ones(n - 1))         # u_kᵀ M_r g1
        terms = U[-1, :] * drive / (lam / sps)                  # u_k(끝) · 구동 / (h λ_k)
        total = float(terms.sum())
        rel = abs(total - v_ref) / v_ref
        worst = max(worst, rel)
        rows.append({'n': n, 'BE_U': b, 'sps': sps, 'lambda': [float(x) for x in lam[:N_TERMS]],
                     'terms': [float(x) for x in terms[:N_TERMS]], 'sum_all_modes': total, 'v_ss_recorded': v_ref, 'rel_diff': rel,
                     'g_over_h_lambda1': float(G.GRAV / (lam[0] / sps)), 'first_term_over_sum': float(terms[0] / total),
                     'shape_factor_1': float(U[-1, 0] * (U[:, 0] @ (mrow * np.ones(n - 1))))})
        print('  n %-2d BE_U %-6g 1/%-4d  합 %10.6f (기록 %10.6f, 상대 차 %.1e)  첫 모드 %10.6f · 다음 %s  |  g/(hλ1) %10.6f' % (
            n, b, sps, total, v_ref, rel, terms[0], ' '.join('%+.6f' % x for x in terms[1:4]), rows[-1]['g_over_h_lambda1']))
    if worst > REL_TOL:
        print('합이 기록의 v_ss 와 다르다(최대 상대 차 %.2e) — 결과를 쓰지 않는다' % worst)
        return 1
    res = {'dll': dll, 'engine': {'git_commit': info.get('git_commit'), 'dirty': info.get('dirty')}, 'gravity': G.GRAV, 'plate_cm': G.PLATE,
           'formula': 'V = sum_k u_k(tip) * (u_k^T M_r g 1) / (h lambda_k), k_d = 0; K_r u = lambda M_r u, u^T M_r u = 1',
           'n_terms_listed': N_TERMS, 'max_rel_diff_to_recorded_v_ss': worst, 'runs': rows, 'seconds': time.perf_counter() - t0}
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
        fh.write('\n')
    print('결과 %s (%.0f s)' % (os.path.relpath(OUT, G.ROOT), res['seconds']))
    return 0


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main())
