"""Baraff 굽힘의 압축 명세 — 같은 θ 를 다르게 쓴다(속도 게이트의 보고 항목, docs/GATE_SPEED.md).

atan2 는 두 인자에 같은 양수를 곱해도 값이 같다. 원 명세의
  sin θ = (n̂A × n̂B)·ê,  cos θ = n̂A·n̂B
에 |nA||nB||e| > 0 을 곱하면
  θ = atan2((nA × nB)·e, (nA·nB)|e|)
가 되어 노멀 정규화가 사라진다. 값은 원 명세와 같아야 한다(시험이 대조한다).
"""
import sympy as sp

NAME = 'BaraffBendCompact'
DESCRIPTION = 'Baraff-Witkin 1998 bend condition, compact form theta = atan2((nA x nB).e, (nA.nB)|e|) - theta0'
N_NODES = 4
N_COMP = 1

POSITION_VAR = 'ppNodes'
VELOCITY_VAR = 'ppVelocities'
STENCIL_VAR = 'ppHinges'
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pTheta0', 1)]


def condition(ctx, X, E):
    xh, xi, xj, xk = X
    nA = ctx.stage((xi - xh).cross(xj - xh))
    nB = ctx.stage((xj - xk).cross(xi - xk))
    e = ctx.stage(xj - xi)
    s = nA.cross(nB).dot(e)
    c = nA.dot(nB) * ctx.norm(e)
    return [sp.atan2(s, c) - E['pTheta0'][0]]
