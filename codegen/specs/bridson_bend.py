"""Bridson, Marino, Fedkiw 2003 의 굽힘을 조건 함수로 — 힌지 스텐실의 명세(docs/GATE_BRIDSON.md).

Bridson 의 탄성 굽힘 힘 F_m = −k w (sin(θ/2) − sin(θ0/2)) u_m 에서 u_m = ∇_m θ 이다.
w 를 상수로 보고 θ0 = 0 이면 이 힘은 에너지 E = 2 k w (1 − cos(θ/2)) = 4 k w sin²(θ/4) 의 −∇E 라서
  C = 2 √(2 w̄) sin(θ/4)
로 두면 ½ k C² = E, f = −k C ∇C = −k w̄ sin(θ/2) ∇θ 가 Bridson 의 힘이고, K = −k(∇C∇Cᵀ + C∇²C) 가 그 정확한 강성이다.
w̄ = 힌지 가중치(엔진은 discrete shells 의 3|ē|²/(A1 + A2)), 힌지마다 상수인 스텐실 값으로 받는다.
θ 는 압축 명세(baraff_bend_compact.py)와 같이 atan2 의 배율 불변성으로 정규화 없이 쓴다.
θ0 ≠ 0 은 다루지 않는다(엔진 기본값 Theta0 = 0, 판정문 「조건 함수로 쓴 Bridson 굽힘」).
"""
import sympy as sp

NAME = 'BridsonBend'
DESCRIPTION = 'Bridson 2003 bending as a condition C = 2 sqrt(2 w) sin(theta/4), theta = atan2((nA x nB).e, (nA.nB)|e|)'
N_NODES = 4
N_COMP = 1

POSITION_VAR = 'ppNodes'
VELOCITY_VAR = 'ppVelocities'
STENCIL_VAR = 'ppHinges'
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pBendW', 1)]   # 힌지 가중치 w̄


def condition(ctx, X, E):
    xh, xi, xj, xk = X
    nA = ctx.stage((xi - xh).cross(xj - xh))
    nB = ctx.stage((xj - xk).cross(xi - xk))
    e = ctx.stage(xj - xi)
    theta = sp.atan2(nA.cross(nB).dot(e), nA.dot(nB) * ctx.norm(e))
    w = E['pBendW'][0]
    return [2 * sp.sqrt(2 * w) * sp.sin(theta / 4)]
