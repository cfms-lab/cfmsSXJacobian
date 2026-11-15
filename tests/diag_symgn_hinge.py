"""힌지 하나에서 네 굽힘 행렬이 강체 회전을 막는지 식으로 견준다 (docs/GATE_JACOBIAN_DAMPING.md 「논문용 측정 4」, 엔진 없음 — 몇 초).

  GN      Gauss–Newton: g gᵀ (g = 접힘 각의 기울기를 네 꼭짓점에 늘어놓은 것)
  symGN   블록 대칭 Gauss–Newton: 비대각 블록 ½(g_a g_bᵀ + g_b g_aᵀ), 대각 = −Σ
  Eq14    Oh · Noh · Wohn(2008) 식 (14) 의 꼴: 같은 법선끼리의 곱만(α_aα_b n̂A n̂Aᵀ + β_aβ_b n̂B n̂Bᵀ), 모든 꼭짓점 쌍, 대각 = −Σ
  Eq10    원고 식 (10): Eq14 에서 함께 쓰는 변의 두 끝 사이 블록을 뺀 것
힌지 강성은 1 로 둔다(네 행렬에 같은 계수 — 비만 본다). 격자의 세 가지 힌지(변이 격자 축을 따르는 것 둘 · 대각선을 따르는 것)를
접힌 각 0 · 10 · 30 · 60 · 90° 로 접는다. 재는 것:
  rot   강체 회전에 대한 Rayleigh 몫의 최댓값(회전 중심 = 힌지의 무게중심, 세 축에 대한 일반화 고윳값) / GN 의 가장 큰 고윳값
  emin  가장 작은 고윳값 / GN 의 가장 큰 고윳값 (음이면 부정치)
접힘 각의 기울기 식은 먼저 유한차분과 견준다(1e-8).

사용:  python tests/diag_symgn_hinge.py        (리포 루트에서) → results/diag_symgn_hinge.json
"""
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'results', 'diag_symgn_hinge.json')
ANGLES = (0.0, 10.0, 30.0, 60.0, 90.0)
# 대각선으로 나눈 정사각 격자(간격 1)의 힌지: (날개 h, 변의 시작 i, 변의 끝 j, 날개 k) — 평평할 때의 좌표
HINGES = {
    'edge along a grid line': np.array([[1.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]),
    'edge along the diagonal': np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
}


def fold(X0, theta):
    """날개 k 를 변 i–j 를 축으로 theta 만큼 돌린다."""
    X = X0.copy()
    axis = (X0[2] - X0[1]) / np.linalg.norm(X0[2] - X0[1])
    r = X0[3] - X0[1]
    c, s = np.cos(theta), np.sin(theta)
    X[3] = X0[1] + r * c + np.cross(axis, r) * s + axis * np.dot(axis, r) * (1.0 - c)
    return X


def theta_of(X):
    h, i, j, k = X
    e = (j - i) / np.linalg.norm(j - i)
    nA, nB = np.cross(i - h, j - h), np.cross(j - k, i - k)
    nA, nB = nA / np.linalg.norm(nA), nB / np.linalg.norm(nB)
    return np.arctan2(np.dot(np.cross(nA, nB), e), np.dot(nA, nB))


def gradient(X):
    """접힘 각의 정확한 기울기 g_m = α_m n̂A + β_m n̂B (m = h, i, j, k). (g, n̂A, n̂B, α, β)"""
    h, i, j, k = X
    e = j - i
    el = np.linalg.norm(e)
    nA, nB = np.cross(i - h, j - h), np.cross(j - k, i - k)
    la, lb = np.linalg.norm(nA), np.linalg.norm(nB)
    nA, nB = nA / la, nB / lb
    th, tk = np.dot(h - i, e) / el ** 2, np.dot(k - i, e) / el ** 2          # 날개 꼭짓점을 변에 내린 자리(0 = i, 1 = j)
    ah, bk = -el / la, -el / lb
    alpha = np.array([ah, -(1.0 - th) * ah, -th * ah, 0.0])
    beta = np.array([0.0, -(1.0 - tk) * bk, -tk * bk, bk])
    return alpha[:, None] * nA[None, :] + beta[:, None] * nB[None, :], nA, nB, alpha, beta


def fd_gradient(X, d=1e-6):
    g = np.zeros((4, 3))
    for a in range(4):
        for c in range(3):
            Xp, Xm = X.copy(), X.copy()
            Xp[a, c] += d
            Xm[a, c] -= d
            g[a, c] = (theta_of(Xp) - theta_of(Xm)) / (Xp[a, c] - Xm[a, c])
    return g


def assemble(block):
    K = np.zeros((12, 12))
    for a in range(4):
        for b in range(4):
            if a != b:
                B = block(a, b)
                K[3 * a:3 * a + 3, 3 * b:3 * b + 3] += B
                K[3 * a:3 * a + 3, 3 * a:3 * a + 3] -= B
    return K


def matrices(X):
    g, nA, nB, al, be = gradient(X)
    same = lambda a, b: al[a] * al[b] * np.outer(nA, nA) + be[a] * be[b] * np.outer(nB, nB)     # noqa: E731
    return {'GN': np.outer(g.ravel(), g.ravel()),
            'symGN': assemble(lambda a, b: 0.5 * (np.outer(g[a], g[b]) + np.outer(g[b], g[a]))),
            'Eq14': assemble(same),
            'Eq10': assemble(lambda a, b: np.zeros((3, 3)) if {a, b} == {1, 2} else same(a, b))}


def rotation_quotient(K, X):
    """max_ω δ(ω)ᵀ K δ(ω) / δ(ω)ᵀ δ(ω), δ = ω × (x − 무게중심)."""
    c = X.mean(axis=0)
    D = np.column_stack([np.cross(np.eye(3)[a], X - c).ravel() for a in range(3)])
    A, B = D.T @ (0.5 * (K + K.T)) @ D, D.T @ D
    w = np.linalg.eigvals(np.linalg.solve(B, A)).real
    return float(w.max()), float(w.min())


def main():
    res = {'angles_deg': list(ANGLES), 'hinges': {}, 'gradient_vs_fd': 0.0}
    worst = 0.0
    for name, X0 in HINGES.items():
        rows = []
        for deg in ANGLES:
            X = fold(X0, np.radians(deg))
            g = gradient(X)[0]
            fd = fd_gradient(X)
            worst = max(worst, float(np.abs(g - fd).max() / np.abs(fd).max()))
            Ks = matrices(X)
            scale = float(np.linalg.eigvalsh(Ks['GN']).max())
            row = {'angle_deg': deg, 'theta_check_deg': float(np.degrees(theta_of(X)))}
            for k, K in Ks.items():
                hi, lo = rotation_quotient(K, X)
                row[k] = {'rot_max': hi / scale, 'rot_min': lo / scale, 'eig_min': float(np.linalg.eigvalsh(0.5 * (K + K.T)).min() / scale),
                          'translation': float(np.abs(K @ np.tile([0.3, -0.2, 0.5], 4)).max() / scale),
                          'block_asymmetry': float(max(np.abs(K[3 * a:3 * a + 3, 3 * b:3 * b + 3] - K[3 * a:3 * a + 3, 3 * b:3 * b + 3].T).max()
                                                       for a in range(4) for b in range(4) if a != b) / scale)}
            rows.append(row)
        res['hinges'][name] = rows
    res['gradient_vs_fd'] = worst
    ok = worst <= 1e-8
    print('접힘 각의 기울기 식 대 유한차분: 상대 %.1e (%s, 기준 1e-8)' % (worst, 'PASS' if ok else 'FAIL'))
    for name, rows in res['hinges'].items():
        print('\n%s' % name)
        print('  %5s  %s' % ('각', '   '.join('%-22s' % k for k in ('GN', 'symGN', 'Eq14', 'Eq10'))))
        for row in rows:
            print('  %4.0f°  ' % row['angle_deg'] + '   '.join('rot %.3f emin %+.3f' % (row[k]['rot_max'], row[k]['eig_min']) for k in ('GN', 'symGN', 'Eq14', 'Eq10')))
    res['valid'] = bool(ok)
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
        fh.write('\n')
    print('\n→ %s' % os.path.relpath(OUT, ROOT))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
