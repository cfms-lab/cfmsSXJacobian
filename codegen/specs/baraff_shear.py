"""Baraff–Witkin 1998 §4.2 전단 조건 — 삼각형 스텐실의 명세.

  C = √a w_u · w_v,   a = uv 넓이(원문은 a 를 곱한다. √a 는 사용자 결정 2026-11-13, 인장 명세와 같다)
w_u, w_v 와 a 는 인장 명세(baraff_stretch.py)와 같은 식이다. 원문은 인장 · 전단의 강성을 따로 두므로
EvalImplicit 의 계수 {k, kd, h} 도 따로 받게 명세를 나눈다. docs/GATE_STRETCH_SHEAR.md.
"""
import sympy as sp

NAME = 'BaraffShear'
DESCRIPTION = 'Baraff-Witkin 1998 shear condition C = sqrt(a) w_u . w_v on triangle (0, 1, 2)'
N_NODES = 3
N_COMP = 1

POSITION_VAR = 'ppNodes'
VELOCITY_VAR = 'ppVelocities'
STENCIL_VAR = 'ppTriangles'
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pUV', 6)]   # 삼각형마다 쉬는 uv: u0 v0 u1 v1 u2 v2


def condition(ctx, X, E):
    x0, x1, x2 = X
    u0, v0, u1, v1, u2, v2 = E['pUV']
    du1, dv1, du2, dv2 = u1 - u0, v1 - v0, u2 - u0, v2 - v0
    det = du1 * dv2 - du2 * dv1
    dx1, dx2 = x1 - x0, x2 - x0
    wu = ctx.stage((dx1 * dv2 - dx2 * dv1) / det)
    wv = ctx.stage((-dx1 * du2 + dx2 * du1) / det)
    a = sp.Abs(det) / 2
    return [sp.sqrt(a) * wu.dot(wv)]
