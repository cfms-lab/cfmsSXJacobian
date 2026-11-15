"""Baraff–Witkin 1998 §4.3 굽힘 조건 — 사용자가 쓰는 명세의 예.

이 파일에 적는 것은 (1) 조건 함수 C(x) 와 (2) 그것이 읽는 SXMesh 변수 이름뿐이다.
∇C · ∇²C 의 유도와 C++ 코드는 생성기가, 힘 · Jacobian · 감쇠의 조립은 src/sx_condition_force.h 가 한다.

힌지 = 변 (i, j) 를 공유하는 두 삼각형. 꼭짓점 순서 (h, i, j, k).
  삼각형 A = (h, i, j), 삼각형 B = (k, j, i) — 둘 다 자기 감김 순서대로 적은 것이라 평평하면 θ = 0.
  sin θ = (n̂A × n̂B) · ê,  cos θ = n̂A · n̂B,  C = θ − θ0
"""
import sympy as sp

NAME = 'BaraffBend'
DESCRIPTION = 'Baraff-Witkin 1998 bend condition C = theta - theta0 on hinge (h, i, j, k)'
N_NODES = 4
N_COMP = 1

POSITION_VAR = 'ppNodes'          # 노드 좌표 (nNode × 3)
VELOCITY_VAR = 'ppVelocities'     # 노드 속도 (nNode × 3), 감쇠에 쓴다
STENCIL_VAR = 'ppHinges'          # 힌지마다 꼭짓점 번호 (nHinge × 4, 정수)
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pTheta0', 1)]   # 힌지마다 쉬는 각


def condition(ctx, X, E):
    xh, xi, xj, xk = X
    # ctx.stage = 중간 변수 표시(선택). 생성기가 C 를 이 변수들로 먼저 미분하고 연쇄 법칙을 코드에서 잇는다
    # (EvalStaged). 표시하지 않아도 결과는 같고 느릴 뿐이다 — docs/GATE_SPEED.md.
    nA = ctx.stage((xi - xh).cross(xj - xh))
    nB = ctx.stage((xj - xk).cross(xi - xk))
    e = ctx.stage(xj - xi)
    a, b, eh = ctx.unit(nA), ctx.unit(nB), ctx.unit(e)
    s = a.cross(b).dot(eh)
    c = a.dot(b)
    return [sp.atan2(s, c) - E['pTheta0'][0]]
