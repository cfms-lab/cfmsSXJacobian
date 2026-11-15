"""엔진 cfmsDrape 의 전단 접선 다항식 법칙을 조건 함수로 — 엔진 비교용(docs/GATE_ENGINE_FORCES.md).

엔진(`set_shear_tangent_polynomial(coefficients, max_angle)`)은 전단각을
  q = 2 w_u·w_v / (|w_u|² + |w_v|²),   γ = asin(q)
로 재고, 접선 K(γ) = Σ aᵢ|γ|ⁱ, 합력 N(γ) = ∫K, 넓이당 에너지 U(γ) = Σ aᵢ|γ|^(i+2) / ((i+1)(i+2)) 를 쓴다(SE 는 안 쓴다).
½ C² = a U 가 되게
  C = √a γ √(Σ 2 aᵢ|γ|ⁱ / ((i+1)(i+2))),   k = 1
로 두면 힘 −C∇C = −a N ∇γ 가 엔진의 힘이고 정확판 J 가 그 참 Hessian 이다. |γ| 는 √(γ²) 로 쓴다 — `Abs` 로 쓰면
미분이 `Derivative` 로 남는다. 최대 각 너머(접선을 고정하는 갈래)와 q 의 자르기는 넣지 않았다.
스텐실 값 = 쉬는 uv 6 칸 + 계수 3 칸(a₀ a₁ a₂). Mathematica 게이트를 거치지 않았다 — 유한차분과 엔진 값으로만 확인했다.
"""
import sympy as sp

NAME = 'ShearPolynomial'
DESCRIPTION = 'shear tangent polynomial: C = sqrt(a) gamma sqrt(sum 2 a_i |gamma|^i / ((i+1)(i+2))), gamma = asin(2 wu.wv / (wu.wu + wv.wv))'
N_NODES = 3
N_COMP = 1

POSITION_VAR = 'ppNodes'
VELOCITY_VAR = 'ppVelocities'
STENCIL_VAR = 'ppTriangles'
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pUV', 6), ('pPoly', 3)]   # 접선 다항식의 계수 a₀ a₁ a₂ (dyn/cm)


def condition(ctx, X, E):
    x0, x1, x2 = X
    u0, v0, u1, v1, u2, v2 = E['pUV']
    a = E['pPoly']
    du1, dv1, du2, dv2 = u1 - u0, v1 - v0, u2 - u0, v2 - v0
    det = du1 * dv2 - du2 * dv1
    dx1, dx2 = x1 - x0, x2 - x0
    wu = ctx.stage((dx1 * dv2 - dx2 * dv1) / det)
    wv = ctx.stage((-dx1 * du2 + dx2 * du1) / det)
    q = 2 * wu.dot(wv) / (wu.dot(wu) + wv.dot(wv))
    gam = sp.asin(q)
    ag = sp.sqrt(gam**2)                              # |γ|
    s = sum(2 * a[i] * ag**i / ((i + 1) * (i + 2)) for i in range(3))
    return [sp.sqrt(sp.Abs(det) / 2) * gam * sp.sqrt(s)]
