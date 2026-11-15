"""판정 밖 계산 (docs/GATE_JACOBIAN_DAMPING.md 「Oh 외(2008) 식 (14) 와의 대조」): 블록 대칭 굽힘 행렬의 장력 계수 κ = c̄ Δ² / k_b 를 식에서 바로 계산한다.

사용:  python tests/diag_block_symmetric_kappa.py     (엔진 불필요 — numpy 만)
평평한 30 cm 판(한 변 노드 n, 장면 P 의 삼각화)에서 힌지마다 ∇_mθ = α_m n̂A + β_m n̂B 의 α · β 를 구하고(tests/test_bridson.cpp 의 BridsonU 와 같은 식 —
Oh · Noh · Wohn 2008 부록 식 (13) 의 q), 같은 법선끼리의 블록 −k w̄ (α_mα_n + β_mβ_n) n̂n̂ᵀ 로 법선 변위 w = s 의 2 차 형식을 낸다. 대각 블록은 행 합이다.
  engine  : 변의 두 끝 (i, j) 사이의 블록을 두지 않는다 — 엔진의 변형(GATE_BRIDSON B7). 엔진에서 읽은 κ = 17.70 · 17.85 · 17.93 과 맞아야 한다.
  eq14    : (i, j) 블록을 둔다 — Oh 외 식 (14) 그대로.
  gn      : 교차 항까지 둔 ∇θ∇θᵀ(Gauss–Newton 의 꼴) — 강체 회전에서 0 이어야 한다.
가중치는 엔진과 같은 w̄ = 3|e|² / (A_A + A_B), k_b = 1e4. 결과는 화면에만 찍는다.
"""
import sys
sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import link_l0 as L0  # noqa: E402

PLATE, KB = 30.0, 1.0e4


def scene(n):
    xs = np.linspace(-PLATE / 2, PLATE / 2, n)
    gx, gz = np.meshgrid(xs, xs)
    V = np.column_stack([gx.ravel(), np.full(gx.size, 20.0), gz.ravel()])
    F = []
    for j in range(n - 1):
        for i in range(n - 1):
            a = j * n + i
            F += [[a, a + 1, (j + 1) * n + i + 1], [a, (j + 1) * n + i + 1, (j + 1) * n + i]]
    return V, np.array(F, np.int32), gz.ravel() + PLATE / 2


def quad(n, kind):
    V, F, s = scene(n)
    total = 0.0
    for h, i, j, k in L0.hinges(F):
        xh, xi, xj, xk = V[h], V[i], V[j], V[k]
        e = xj - xi
        nA, nB = np.cross(xi - xh, xj - xh), np.cross(xj - xk, xi - xk)
        la, lb, le = np.linalg.norm(nA), np.linalg.norm(nB), np.linalg.norm(e)
        alpha = np.array([le / la, e.dot(xh - xj) / (le * la), e.dot(xi - xh) / (le * la), 0.0])
        beta = np.array([0.0, -e.dot(xj - xk) / (le * lb), -e.dot(xk - xi) / (le * lb), le / lb])
        kw = KB * 3.0 * le * le / (0.5 * la + 0.5 * lb)
        w = s[[h, i, j, k]]                                   # 법선 변위(강체 회전 δ = s n̂), 두 법선이 같은 방향인 평평한 판
        if kind == 'gn':
            g = alpha + beta
            total += kw * float(g.dot(w)) ** 2
            continue
        for m in range(4):
            for q in range(m + 1, 4):
                if kind == 'engine' and (m, q) == (1, 2):
                    continue
                c = kw * (alpha[m] * alpha[q] + beta[m] * beta[q])   # 블록 J_mq = −c n̂n̂ᵀ, 대각 = 행 합
                total += -c * (w[m] - w[q]) ** 2
    return total / (PLATE * PLATE)                              # c̄ = −wᵀJw / (W ∫ w′² ds), W = L, ∫ = L


def main():
    print('블록 대칭 굽힘 행렬의 장력 계수 (평평한 %g cm 판, k_b = %g erg, w̄ = 3|e|²/(A_A + A_B))' % (PLATE, KB))
    for n in (21, 41, 81):
        dx = PLATE / (n - 1)
        row = {k: quad(n, k) for k in ('engine', 'eq14', 'gn')}
        print('  n %-2d  엔진의 변형 c̄ %.5g (κ %.3f)  ·  식 (14) c̄ %.5g (κ %.3f)  ·  Gauss–Newton 꼴 %.2e' % (
            n, row['engine'], row['engine'] * dx * dx / KB, row['eq14'], row['eq14'] * dx * dx / KB, row['gn']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
