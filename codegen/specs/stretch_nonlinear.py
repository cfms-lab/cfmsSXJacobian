"""엔진 cfmsDrape 의 비선형 막 법칙(인장 성분)을 조건 함수로 — 엔진 비교용(docs/GATE_ENGINE_FORCES.md).

엔진(`set_membrane_nonlinear(cb, eps_crit, k_hard)`)은 성분마다 변형률 ε = |w| − 1 에서 힘 몫을 k_sec(ε) ε 로 둔다:
  k_sec = cb k (ε < 0, 눌림 무름),  k (0 ≤ ε ≤ ε_c),  k + k_hard (ε − ε_c) (그 위, 늘어남 굳음).
½ k C² 이 아니라 넓이당 에너지 W(ε)(W′ = k_sec ε)다. 그래도 ½ k g(ε)² = W(ε) 가 되게 g 를 잡으면 조건 한 줄이 된다:
  C = √a g(|w| − 1),
  g = √cb ε (ε < 0),   ε (0 ≤ ε ≤ ε_c),   √(ε² + κ (⅔(ε³ − ε_c³) − ε_c(ε² − ε_c²))) (그 위),   κ = k_hard / k.
그러면 정확판 J 가 W 의 참 Hessian 이고, C > 0 판(늘어난 성분만 C∇²C 를 남김)이 엔진 J 와 같다.
스텐실 값 = 쉬는 uv 6 칸 + 법칙 3 칸(cb, ε_c, κ). Mathematica 게이트를 거치지 않았다 — 유한차분과 엔진 값으로만 확인했다.
"""
import sympy as sp

NAME = 'StretchNonlinear'
DESCRIPTION = 'nonlinear stretch: C = sqrt(a) g(|w| - 1), piecewise g (compression softening cb, hardening kappa beyond eps_c)'
N_NODES = 3
N_COMP = 2

POSITION_VAR = 'ppNodes'
VELOCITY_VAR = 'ppVelocities'
STENCIL_VAR = 'ppTriangles'
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pUV', 6), ('pLaw', 3)]   # 법칙: cb, ε_c, κ
IMPLICIT_TENSION = True


def g(eps, cb, ec, kap):
    hard = sp.sqrt(eps**2 + kap * (sp.Rational(2, 3) * (eps**3 - ec**3) - ec * (eps**2 - ec**2)))
    return sp.Piecewise((sp.sqrt(cb) * eps, eps < 0), (eps, eps <= ec), (hard, True))


def condition(ctx, X, E):
    x0, x1, x2 = X
    u0, v0, u1, v1, u2, v2 = E['pUV']
    cb, ec, kap = E['pLaw']
    du1, dv1, du2, dv2 = u1 - u0, v1 - v0, u2 - u0, v2 - v0
    det = du1 * dv2 - du2 * dv1
    dx1, dx2 = x1 - x0, x2 - x0
    wu = ctx.stage((dx1 * dv2 - dx2 * dv1) / det)
    wv = ctx.stage((-dx1 * du2 + dx2 * du1) / det)
    sa = sp.sqrt(sp.Abs(det) / 2)
    return [sa * g(ctx.norm(wu) - 1, cb, ec, kap), sa * g(ctx.norm(wv) - 1, cb, ec, kap)]
