"""엔진 cfmsDrape 의 봉제 스프링 — 2 노드 스텐실 · 벡터 조건(성분 3)의 명세. 엔진 비교용(docs/GATE_ENGINE_FORCES.md).

  C = x₁ − x₀            (길이 0 스프링을 벡터 조건으로)
힘 f = −k C·∇C 는 f₀ = k(x₁ − x₀) = −f₁, 감쇠 d = −kd Ċ·∇C 는 d₀ = kd(v₁ − v₀) — 엔진 `_GetTensionSpring2` 가 쓰는 꼴이다.
∇²C = 0 이라 정확판과 Gauss–Newton 판이 같고 J 는 −k I 블록이다. 스칼라 조건 |x₁ − x₀| 로 써도 탄성 힘과 J 는 같지만
감쇠가 변 방향 성분만 남아 엔진과 달라진다. 스텐실 값은 없다.
Mathematica 게이트를 거치지 않았다 — 유한차분과 엔진 값으로만 확인했다.
"""
import sympy as sp

NAME = 'SeamSpring'
DESCRIPTION = 'zero-rest-length seam spring as a vector condition C = x1 - x0 on stitch (0, 1)'
N_NODES = 2
N_COMP = 3

POSITION_VAR = 'ppNodes'
VELOCITY_VAR = 'ppVelocities'
STENCIL_VAR = 'ppStitches'        # 땀마다 두 노드 번호 (nStitch × 2, 정수)
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = []


def condition(ctx, X, E):
    x0, x1 = X
    d = ctx.stage(x1 - x0)
    return [d[0], d[1], d[2]]
