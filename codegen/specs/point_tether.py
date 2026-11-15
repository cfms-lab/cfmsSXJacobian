"""엔진 cfmsDrape 의 꼭짓점 테더 — 1 노드 스텐실 · 벡터 조건(성분 3)의 명세. 엔진 비교용(docs/GATE_ENGINE_FORCES.md).

  C = x − target,   target = 스텐실 값 3 칸
힘 f = −k(x − target), 감쇠 d = −kd v. ∇C = I, ∇²C = 0. 병진 불변이 아니라(Σ∇C ≠ 0) 엔진 꼴 변형(`EvalEngine*` —
대각을 노드쌍 블록의 행 합으로 채운다)은 이 힘에 맞지 않는다. 엔진도 이 힘은 노드쌍 블록이 아니라 계의 대각에 넣는다.
Mathematica 게이트를 거치지 않았다 — 유한차분과 엔진 값으로만 확인했다.
"""
import sympy as sp

NAME = 'PointTether'
DESCRIPTION = 'vertex tether as a vector condition C = x - target on node (0)'
N_NODES = 1
N_COMP = 3

POSITION_VAR = 'ppNodes'
VELOCITY_VAR = 'ppVelocities'
STENCIL_VAR = 'ppTethers'         # 테더마다 노드 번호 (nTether × 1, 정수)
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pTarget', 3)]


def condition(ctx, X, E):
    x0, = X
    c = ctx.stage(x0 - sp.Matrix(E['pTarget']))
    return [c[0], c[1], c[2]]
