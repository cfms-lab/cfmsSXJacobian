"""엔진 cfmsDrape 의 비선형 막 법칙(Provot 전단 성분)을 조건 함수로 — 엔진 비교용(docs/GATE_ENGINE_FORCES.md).

엔진은 전단의 두 대각선 성분에 굳음만 준다(눌림 무름은 없다). 변형률은 ε = |w_u ± w_v| / b − 1 이고
굳음 계수는 k_hard · SE / E 라서 κ = k_hard / E 로 인장과 같다:
  C = √a b g(|w_u ± w_v| / b − 1),
  g = ε (ε ≤ ε_c),   √(ε² + κ (⅔(ε³ − ε_c³) − ε_c(ε² − ε_c²))) (그 위).
ε ≤ ε_c 에서는 선형 Provot 명세(provot_shear.py)의 C = √a(|w_u ± w_v| − b) 와 같다.
스텐실 값 = 쉬는 uv 6 칸 + 쉬는 대각선 길이 2 칸 + 법칙 2 칸(ε_c, κ).
Mathematica 게이트를 거치지 않았다 — 유한차분과 엔진 값으로만 확인했다.
"""
import sympy as sp

NAME = 'ProvotShearNonlinear'
DESCRIPTION = 'nonlinear Provot shear: C = sqrt(a) b g(|w_u +- w_v| / b - 1), hardening kappa beyond eps_c'
N_NODES = 3
N_COMP = 2

POSITION_VAR = 'ppNodes'
VELOCITY_VAR = 'ppVelocities'
STENCIL_VAR = 'ppTriangles'
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pUV', 6), ('pDiag', 2), ('pLaw', 2)]   # 법칙: ε_c, κ
IMPLICIT_TENSION = True


def g(eps, ec, kap):
    hard = sp.sqrt(eps**2 + kap * (sp.Rational(2, 3) * (eps**3 - ec**3) - ec * (eps**2 - ec**2)))
    return sp.Piecewise((eps, eps <= ec), (hard, True))


def condition(ctx, X, E):
    x0, x1, x2 = X
    u0, v0, u1, v1, u2, v2 = E['pUV']
    bp, bm = E['pDiag']
    ec, kap = E['pLaw']
    du1, dv1, du2, dv2 = u1 - u0, v1 - v0, u2 - u0, v2 - v0
    det = du1 * dv2 - du2 * dv1
    dx1, dx2 = x1 - x0, x2 - x0
    wu = (dx1 * dv2 - dx2 * dv1) / det
    wv = (-dx1 * du2 + dx2 * du1) / det
    dp = ctx.stage(wu + wv)
    dm = ctx.stage(wu - wv)
    sa = sp.sqrt(sp.Abs(det) / 2)
    return [sa * bp * g(ctx.norm(dp) / bp - 1, ec, kap), sa * bm * g(ctx.norm(dm) / bm - 1, ec, kap)]
