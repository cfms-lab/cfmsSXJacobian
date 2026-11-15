"""Provot 식 전단 조건(대각선 늘림) — 삼각형 스텐실 · 벡터 조건(성분 2)의 명세.

엔진 cfmsDrape 의 기본 전단(SDEisfProvot_Invariant, `_GetImplicitShear_Provot_WithInvariant`)과 같은 모형이다.
Baraff–Witkin 전단(w_u · w_v, baraff_shear.py)과 다르다 — 재질 틀의 두 대각선 방향 w_u ± w_v 의 길이를 본다.
  C = √a (|w_u + w_v| − b₊,  |w_u − w_v| − b₋),   a = uv 넓이
w_u, w_v 와 a 는 인장 명세(baraff_stretch.py)와 같은 식이다. 쉬는 상태에서 |w_u| = |w_v| = 1, w_u ⊥ w_v 이면
|w_u ± w_v| = √2 다. 그 쉬는 대각선 길이 b₊, b₋ 를 삼각형마다 스텐실 값으로 받는다 — 엔진은 √2 를 float 리터럴
(1.41421356237f)로 두고 있어, 값을 맞춰 견주려면 엔진의 값을 그대로 넘길 수 있어야 한다.
엔진은 C = a(…) 에 강성 SE / a 를 쓴다. √a 를 곱하고 강성을 k = SE 로 두면 힘 · 에너지가 같다(인장 명세와 같은 처방).
docs/GATE_PROVOT_SHEAR.md.
"""
import sympy as sp

NAME = 'ProvotShear'
DESCRIPTION = 'Provot-type shear (diagonal stretch) C = sqrt(a) (|w_u + w_v| - b_plus, |w_u - w_v| - b_minus) on triangle (0, 1, 2)'
N_NODES = 3
N_COMP = 2

POSITION_VAR = 'ppNodes'
VELOCITY_VAR = 'ppVelocities'
STENCIL_VAR = 'ppTriangles'
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pUV', 6),       # 삼각형마다 쉬는 uv: u0 v0 u1 v1 u2 v2
                    ('pDiag', 2)]     # 삼각형마다 쉬는 대각선 길이 b₊ b₋ (직교 · 단위 틀이면 √2 √2)

# 늘어난 성분(C > 0)의 기하 강성만 남기는 판(EvalImplicitTension)도 낸다. 성분마다 C = √a(|d| − b) 이고 d = w_u ± w_v 는
# 좌표의 1 차식이라 ∇²C 가 준정부호다(인장 명세와 같은 근거). 엔진 손 코드의 기본 처방이 이것이다.
IMPLICIT_TENSION = True


def condition(ctx, X, E):
    x0, x1, x2 = X
    u0, v0, u1, v1, u2, v2 = E['pUV']
    bp, bm = E['pDiag']
    du1, dv1, du2, dv2 = u1 - u0, v1 - v0, u2 - u0, v2 - v0
    det = du1 * dv2 - du2 * dv1
    dx1, dx2 = x1 - x0, x2 - x0
    wu = (dx1 * dv2 - dx2 * dv1) / det          # 인장 명세와 같은 w_u, w_v
    wv = (-dx1 * du2 + dx2 * du1) / det
    dp = ctx.stage(wu + wv)                     # 두 대각선(좌표의 1 차식) — 중간 변수 표시
    dm = ctx.stage(wu - wv)
    sa = sp.sqrt(sp.Abs(det) / 2)
    return [sa * (ctx.norm(dp) - bp),
            sa * (ctx.norm(dm) - bm)]
