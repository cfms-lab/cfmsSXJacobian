"""진단용 명세(판정 밖, docs/GATE_ENGINE.md 「판정 뒤 진단」) — 엔진 Macri 코드가 미분한다고 읽은 두 양.

조건 함수가 아니라 엔진 중간 배열을 견줄 기호 미분을 얻으려고 생성기를 빌려 쓴다. 사용자 명세(codegen/specs/)가 아니다.
엔진 꼭짓점 (H, I, J, K) = (h, i, j, k), 엔진의 벡터 그대로(정규화하지 않는다, SDECore.cpp _Input_b012):
  nA = qA[0] × qA[2] = (xJ − xI) × (xI − xH)
  nB = qB[1] × qB[2] = (xJ − xK) × (xK − xI)
  e  = qB[3]         = xI − xJ
  P  = nA · nB         엔진 _dcosth_dxm 이 이것의 미분이라고 읽었다
  Q  = (nA × nB) · e   엔진 _dsinth_dxm 이 이것의 미분이라고 읽었다
"""

NAME = 'EngineTerms'
DESCRIPTION = 'diagnostic: P = nA.nB and Q = (nA x nB).e with the engine Macri unnormalized vectors'
N_NODES = 4
N_COMP = 2

POSITION_VAR = 'ppNodes'
VELOCITY_VAR = 'ppVelocities'
STENCIL_VAR = 'ppHinges'
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = []


def condition(ctx, X, E):
    xh, xi, xj, xk = X
    nA = (xj - xi).cross(xi - xh)
    nB = (xj - xk).cross(xk - xi)
    e = xi - xj
    return [nA.dot(nB), nA.cross(nB).dot(e)]
