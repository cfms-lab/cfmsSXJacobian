"""Sul, Kim, Park 2010 (Textile Research Journal 80(8): 760–766) 의 Poisson 조건 — 삼각형 스텐실 · 성분 1 의 명세.

삼각형 꼭짓점 (0, 1, 2), 쉬는 상태(패턴) 좌표 (uᵢ, vᵢ) 와 Poisson 계수 β 를 스텐실 값으로 받는다.
  (Δx₁ Δx₂) = (w_u w_v) [[Δu₁, Δu₂], [Δv₁, Δv₂]],   Δx₁ = x₁ − x₀ …
  C = √a (log|w_u| + β log|w_v|),   a = uv 넓이            (원문 식 (6), 원문은 a = 1 로 둔다)
힘은 f = −k C ∇C (원문 식 (9), k = k_P). β 는 재질 값이다 — 원문은 β = 1 · 0 · −1 이 넓이 보존 · 가로 변형 없음 ·
가로로 같이 늘어남에 해당한다고 둔다. √a 를 곱하는 것은 인장 · 전단 명세와 같은 처방이다(에너지가 넓이에 비례).

∇²(log|w|) = (I − 2ŵŵᵀ)/|w|² 는 부정치다. 그래서 인장의 C > 0 판(IMPLICIT_TENSION)은 근거가 없어 내지 않는다 —
강성 판은 정확판(EvalImplicit)과 Gauss–Newton(EvalImplicitGN) 둘이다. docs/GATE_POISSON.md.
"""
import sympy as sp

NAME = 'SulPoisson'
DESCRIPTION = 'Sul 2010 Poisson condition C = sqrt(a) (log|w_u| + beta log|w_v|) on triangle (0, 1, 2)'
N_NODES = 3
N_COMP = 1

POSITION_VAR = 'ppNodes'          # 노드 좌표 (nNode × 3)
VELOCITY_VAR = 'ppVelocities'     # 노드 속도 (nNode × 3), 감쇠에 쓴다
STENCIL_VAR = 'ppTriangles'       # 삼각형마다 꼭짓점 번호 (nTri × 3, 정수)
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pUV', 6),   # 삼각형마다 쉬는 uv: u0 v0 u1 v1 u2 v2
                    ('pBeta', 1)]  # 삼각형마다 Poisson 계수 β (원문 식 (6))


def condition(ctx, X, E):
    x0, x1, x2 = X
    u0, v0, u1, v1, u2, v2 = E['pUV']
    beta = E['pBeta'][0]
    du1, dv1, du2, dv2 = u1 - u0, v1 - v0, u2 - u0, v2 - v0
    det = du1 * dv2 - du2 * dv1
    dx1, dx2 = x1 - x0, x2 - x0
    wu = ctx.stage((dx1 * dv2 - dx2 * dv1) / det)      # (Δx₁ Δx₂) [[Δu₁ Δu₂], [Δv₁ Δv₂]]⁻¹ 의 첫 열
    wv = ctx.stage((-dx1 * du2 + dx2 * du1) / det)     # 둘째 열
    sa = sp.sqrt(sp.Abs(det) / 2)                      # √a
    return [sa * (sp.log(ctx.norm(wu)) + beta * sp.log(ctx.norm(wv)))]
