"""진단용 명세(판정 밖, docs/GATE_ENGINE.md 「판정 뒤 진단」) — 곱의 법칙의 항을 하나씩 떼어 미분한다.

엔진 _dcosth_dxm · _dsinth_dxm 이 ∇P · ∇Q 와 다르므로(diag_engine_terms 의 (가)), 어느 항이 어긋나는지 보려고
다른 벡터를 스텐실 값(고정값)으로 넣고 한 벡터만 x 의 함수로 둔다. 고정값 = 같은 x 에서 잰 nA · nB · e 를 넣으면
  ∇P = ∇PA + ∇PB,  ∇Q = ∇QA + ∇QB + ∇QE
엔진 벡터 그대로(diag_engine_terms_spec.py 와 같다): nA = (xJ − xI) × (xI − xH), nB = (xJ − xK) × (xK − xI), e = xI − xJ.
"""
import sympy as sp

NAME = 'EngineTermsSplit'
DESCRIPTION = 'diagnostic: product-rule pieces of P = nA.nB and Q = (nA x nB).e with fixed partner vectors'
N_NODES = 4
N_COMP = 5

POSITION_VAR = 'ppNodes'
VELOCITY_VAR = 'ppVelocities'
STENCIL_VAR = 'ppHinges'
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('nAfix', 3), ('nBfix', 3), ('efix', 3)]


def condition(ctx, X, E):
    xh, xi, xj, xk = X
    nA = (xj - xi).cross(xi - xh)
    nB = (xj - xk).cross(xk - xi)
    e = xi - xj
    nAf = sp.Matrix(E['nAfix'])
    nBf = sp.Matrix(E['nBfix'])
    ef = sp.Matrix(E['efix'])
    return [nA.dot(nBf),             # PA
            nAf.dot(nB),             # PB
            nA.cross(nBf).dot(ef),   # QA
            nAf.cross(nB).dot(ef),   # QB
            nAf.cross(nBf).dot(e)]   # QE
