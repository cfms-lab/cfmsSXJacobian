"""Baraff–Witkin 1998 §4.2 인장 조건 — 삼각형 스텐실 · 벡터 조건(성분 2)의 명세.

삼각형 꼭짓점 (0, 1, 2), 쉬는 상태(패턴) 좌표 (uᵢ, vᵢ) 를 스텐실 값으로 받는다.
  (Δx₁ Δx₂) = (w_u w_v) [[Δu₁, Δu₂], [Δv₁, Δv₂]],   Δx₁ = x₁ − x₀ …
  C = √a (|w_u| − b_u,  |w_v| − b_v),   a = uv 넓이,  b_u = b_v = 1
원문은 a 를 곱한다. √a 를 곱하면 에너지 ½k|C|² 가 넓이에 비례한다(사용자 결정 2026-11-13).
역행렬과 넓이는 여기서 기호로 쓴다 — 생성 코드가 계산한다(미리 계산해 두는 손 코드가 없다).
docs/GATE_STRETCH_SHEAR.md.
"""
import sympy as sp

NAME = 'BaraffStretch'
DESCRIPTION = 'Baraff-Witkin 1998 stretch condition C = sqrt(a) (|w_u| - 1, |w_v| - 1) on triangle (0, 1, 2)'
N_NODES = 3
N_COMP = 2

POSITION_VAR = 'ppNodes'          # 노드 좌표 (nNode × 3)
VELOCITY_VAR = 'ppVelocities'     # 노드 속도 (nNode × 3), 감쇠에 쓴다
STENCIL_VAR = 'ppTriangles'       # 삼각형마다 꼭짓점 번호 (nTri × 3, 정수)
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pUV', 6)]   # 삼각형마다 쉬는 uv: u0 v0 u1 v1 u2 v2

B_U = B_V = 1                     # 원문의 b_u, b_v(늘이지 않은 상태)

# 인장 C > 0 판(EvalImplicitTension)도 낸다. 성분마다 C = √a(|w| − 1) 은 좌표의 1 차식 w 의 길이라 ∇²C 가
# 준정부호다 — 그래서 늘어난 성분의 C∇²C 만 남긴 A 가 언제나 준정부호다. docs/GATE_TENSION.md.
IMPLICIT_TENSION = True


def uv_frame(X, E):
    """원문 식: w_u, w_v(좌표의 1 차식)와 uv 넓이 a."""
    x0, x1, x2 = X
    u0, v0, u1, v1, u2, v2 = E['pUV']
    du1, dv1, du2, dv2 = u1 - u0, v1 - v0, u2 - u0, v2 - v0
    det = du1 * dv2 - du2 * dv1
    dx1, dx2 = x1 - x0, x2 - x0
    wu = (dx1 * dv2 - dx2 * dv1) / det          # (Δx₁ Δx₂) [[Δu₁ Δu₂], [Δv₁ Δv₂]]⁻¹ 의 첫 열
    wv = (-dx1 * du2 + dx2 * du1) / det         # 둘째 열
    a = sp.Abs(det) / 2
    return wu, wv, a


def condition(ctx, X, E):
    wu, wv, a = uv_frame(X, E)
    wu, wv = ctx.stage(wu), ctx.stage(wv)       # 중간 변수 표시(선택, docs/GATE_SPEED.md)
    sa = sp.sqrt(a)                             # √a (원문은 a)
    return [sa * (ctx.norm(wu) - B_U),
            sa * (ctx.norm(wv) - B_V)]
