# cfmsSXJacobian — 조건 함수 한 줄에서 의복 시뮬레이션의 힘 · Jacobian · SoA 를

Baraff & Witkin (1998) 식의 암시적 의복 시뮬레이션은 힘마다 조건 함수 `C(x)` 를 두고,
힘 `f = −k ∇C C`, 강성 `K = −k(∇C∇Cᵀ + C ∇²C)`, 감쇠 `d = −kd ∇C Ċ` 와 그 Jacobian 을 선형계에 넣는다.
예전에는 이 미분을 손으로 유도하고, 속도를 위해 데이터를 structure of arrays(SoA)로 바꾸는 데
시간 대부분을 썼다. 특히 굽힘은 식이 길다.

이 프로젝트는 그 두 가지를 없앤다.

1. 사용자가 **조건 함수만** 쓰면 생성기(SymPy)가 `C · ∇C · ∇²C` 를 기호적으로 미분하고
   공통 부분식을 줄여 C++ 함수로 낸다.
2. 명세의 **변수 선언**에서 SXMesh 변수(이름으로 찾는 변수 단위 SoA 배열)를 만들고,
   스텐실(힌지 · 삼각형) 루프가 그 배열에서 값을 모은다.
3. 힘 · 강성 · 감쇠 Jacobian 조립은 조건 함수와 무관한 공통 코드 한 벌이다.

## 명세 예 — Baraff 굽힘 (`codegen/specs/baraff_bend.py`)

```python
NAME = 'BaraffBend'; N_NODES = 4; N_COMP = 1
POSITION_VAR = 'ppNodes'; VELOCITY_VAR = 'ppVelocities'; STENCIL_VAR = 'ppHinges'
NODE_VARS = [(POSITION_VAR, 3), (VELOCITY_VAR, 3)]
ELEM_DOUBLE_VARS = [('pTheta0', 1)]          # 힌지마다 쉬는 각

def condition(ctx, X, E):
    xh, xi, xj, xk = X
    nA = ctx.stage((xi - xh).cross(xj - xh))       # ctx.stage = 중간 변수 표시(선택, 빨라진다)
    nB = ctx.stage((xj - xk).cross(xi - xk))
    e = ctx.stage(xj - xi)
    a, b, eh = ctx.unit(nA), ctx.unit(nB), ctx.unit(e)
    return [sp.atan2(a.cross(b).dot(eh), a.dot(b)) - E['pTheta0'][0]]
```

```
python codegen/sxgen.py codegen/specs/baraff_bend.py     # → generated/BaraffBend.h
```

`ctx.unit` 은 정확한 정규화(`Eval`)와 원문 §4.3 의 「길이를 상수로 본다」 근사(`EvalFrozen`)를 같은 명세에서
둘 다 만든다. `ctx.stage` 를 두면 같은 정확한 미분을 단계판(`EvalStaged`)으로도 내고, 평가와 원문 식 (6) 의 조립을
한 함수로 합친 `EvalImplicit` 도 낸다. 기본 사용 쪽은 다음과 같다.

```cpp
sx::CreateVariables<sxgen::BaraffBend>(mesh, nNode, nHinge);   // 명세 선언대로 SXMesh 변수 생성
sx::StencilView<sxgen::BaraffBend> view(mesh);
const double prm[3] = { k, kd, h };
for (int s = 0; s < view.size(); ++s) {
    double x[12], v[12], ev[1], C, f[12], A[78], b[12];
    view.Gather(s, x, v, ev);
    sxgen::BaraffBend::EvalImplicit(x, ev, v, prm, &C, f, A, b);
    // f = f + d, A = −h ∂f/∂v − h² ∂f/∂x (위 삼각), b = h(f + d) + h² ∂f/∂x v → 전역 선형계로 흩뿌린다
}
```

같은 입출력의 `EvalImplicitGN` 은 ∇²C 를 버린 Gauss–Newton 판이다 — 힘은 같고 `A` 는 `(h kd + h² k) Σ∇C∇Cᵀ` 라 언제나
양의 준정부호이며, 굽힘에서는 정확판보다 4.6 – 6.7 배 빠르다(`docs/GATE_GAUSS_NEWTON.md`). 지금 기본은 Bridson 굽힘 ·
전단 `EvalImplicitGN`, 인장 `EvalImplicitTension`, 나머지 힘(Baraff 굽힘 · 압축 명세)은 `EvalImplicit` 다(2026-11-13 확정).
시간 적분 비교(`docs/GATE_INTEGRATION.md`)에서 인장 · 전단의 GN 판은 정확판이 버티는 곳에서 발산했다. 인장의
`EvalImplicitTension` 은 늘어난 쪽의 기하 강성만 남긴 판으로, 정확판이 버티던 곳을 잃지 않고 더 많이 버틴다
(`docs/GATE_TENSION.md`). 인장을 이 판으로 두면 전단은 GN 판이 정확판보다 더 버틴다(`docs/GATE_SHEAR_GN.md`).

각 행렬(`K` · `∂d/∂x` · `∂d/∂v`)을 따로 보거나 원 논문이 버린 비대칭 감쇠 항을 남기려면 `EvalStaged` 의
`C · ∇C · ∇²C` 를 `sx::AssembleConditionForce`(`src/sx_condition_force.h`)에 넘긴다.

## 현재 상태

**착수 게이트 GO (2026-11-12)** — `docs/MODEL_PROBLEM.md`. 굽힘 하나로 다음을 확인했다.

- 생성된 `C · ∇C · ∇²C` 가 Mathematica(40 자리)와 1e-14 안, 유한차분과 1e-9 안에서 같다.
- 감쇠 Jacobian(원 논문이 버린 비대칭 항 포함)이 유한차분과 같다.
- SXMesh 변수를 거친 스텐실 루프가 직접 모은 값과 비트까지 같다.
- 원문의 길이 고정 근사는 기울기(=힘)에 변 방향 가짜 항 `cosθ·sinθ·ê/|e|` 를 만든다.

**속도 게이트 GO (2026-11-12)** — `docs/GATE_SPEED.md`. 명세가 중간 변수(두 노멀 · 변)를 `ctx.stage(...)` 로 표시하면
생성기가 연쇄 법칙을 코드에서 잇는다. 같은 명세가 5 배 빨라져(힌지당 약 0.39 µs) 길이 고정 근사판의 1.75 배가 됐고,
`atan2` 의 배율 불변성으로 정규화를 없앤 압축 명세(`baraff_bend_compact.py`)는 근사판보다 빠르다(0.15 µs).

**조립 게이트 NO-GO (2026-11-13 닫음)** — `docs/GATE_ASSEMBLY.md`. 원문 식 (6) 이 쓰는 요소 기여(대칭 행렬 위 삼각 · 우변)만
바로 모으는 손 코드 `AssembleImplicit`(1 차)는 맞고 기존 조립보다 2.2 배 싸다. 생성기가 평가와 조립을 한 함수로 내는
`EvalImplicit`(2 차)는 2.6 배(압축 명세는 약 4 배) 싸지만 힘의 반올림 오차가 미리 정한 기준을 1.7 배 넘었다.
둘 다 미리 정한 기준(조립 비용 1/3 이하)에는 못 미쳤다. 기본 조립 경로로는 `EvalImplicit` 를 쓰기로 했다(2026-11-13).

**엔진 손 코드 비교 게이트 GO (2026-11-13)** — `docs/GATE_ENGINE.md`. 기존 의복 엔진의 손으로 짠 SoA 굽힘 코드
(Baraff 경로)를 바이트 그대로 옮겨 견줬다. 같은 일(모든 힌지의 `C · ∇C · ∇²C`)을 생성 코드가 3.7 배 빨리 한다
(힌지당 0.38 µs 대 1.39 µs). 손 코드의 θ 는 맞지만 ∇θ · ∇²θ 는 정확한 미분도 원문의 길이 고정 근사도 아니었다.

**인장 · 전단 게이트 GO (2026-11-13)** — `docs/GATE_STRETCH_SHEAR.md`. 원문 §4.2 의 인장(성분 2) · 전단 조건을
삼각형 명세 두 개(`baraff_stretch.py` · `baraff_shear.py`)로 썼다. 원문은 조건에 넓이 `a` 를 곱하는데, 에너지가 넓이에
비례하도록 `√a` 를 곱한다. 스텐실 값은 삼각형마다 쉬는 uv 좌표 6 칸 그대로이고,
2 × 2 역행렬과 넓이는 명세 안에서 기호로 쓴다. 생성 코드는 Mathematica 와 1e-15 안, 유한차분과 1e-9 안에서 같다.
성분이 둘인 `EvalImplicit` 도 검산용 조립과 같다. 생성기는 절댓값 출력 하나만 더했다.

**Bridson 굽힘 게이트 — 명세 GO, 엔진 비교 NO-GO (2026-11-13)** — `docs/GATE_BRIDSON.md`. Bridson et al. (2003)의
굽힘 힘은 힌지 가중치 w 를 상수로 보면 `C = 2√(2w) sin(θ/4)` 의 `−kC∇C` 다(`bridson_bend.py`). 생성 코드는
Mathematica · 유한차분과 맞는다. 기존 엔진의 기본 굽힘 손 코드를 바이트 그대로 옮겨 견주니, 힘은 이 조건의 힘과
반올림 수준에서 같고 Jacobian 은 정확한 강성도 Gauss–Newton 꼴도 아니었다(계수 · 교차 항 · 한 블록 · ∇²θ 항이 다르다).
정확한 강성까지 내는 생성 코드는 엔진의 근사 경로보다 3.3 배 느려 속도 기준(엔진 이하)에는 못 미쳤다.

**Gauss–Newton 게이트 GO (2026-11-13)** — `docs/GATE_GAUSS_NEWTON.md`. 생성기가 같은 명세에서 ∇²C 를 쓰지 않는 조립
`EvalImplicitGN` 을 함께 낸다. 다섯 명세 모두에서 「∇²C = 0 인 식 (6) 조립」과 반올림 안에서 같고, Bridson 굽힘에서는
엔진 손 코드보다 빠르다(시간이 엔진의 0.65 배). 버리는 `C∇²C` 는 쉬는 모양 가까이에서 탄성 강성의 5 – 7 % 다. 그래서
Bridson 굽힘의 기본 강성을 이 판으로 정했다.

**시간 적분 비교 (2026-11-13)** — `docs/GATE_INTEGRATION.md`. 두 점에 매단 0.5 m 천(노드 441)을 식 (6) 으로 2 초
적분해, 힘마다 정확판 · Gauss–Newton 판 · 엔진 Bridson 을 120 가지 조합으로 견줬다. 미리 정한 규칙의 권고는 네 힘
모두 「정확판 유지」다. 인장 · 전단은 GN 으로 바꾸면 정확판이 버티는 곳에서 발산했다. 굽힘은 GN 이 안정도 궤적도
바꾸지 않았지만(차이가 시간 간격 오차의 천분의 1), 시간을 견줄 기준 짝이 양쪽 다 불안정해 규칙상 판정할 수 없었다.
뻣뻣한 굽힘으로 따로 재면 GN 은 전체 시간(조립 + 풀이)의 0.92(Baraff) · 0.94(Bridson) 배다 — 시간의 90 % 가 풀이라
조립을 반으로 줄여도 이득이 작다. 굽힘이 부드러운 천은 두 인장 판 어느 쪽으로도 h ≥ 1/120 에서 안정하지 않았다.
정확판 인장은 눌린 곳에서, GN 인장은 늘어난 곳에서 실패한다.

**인장 C > 0 판 게이트 (2026-11-13)** — `docs/GATE_TENSION.md`. 위 두 실패를 함께 풀려고, 늘어난 성분의 기하 강성
`k C ∇²C` 만 남기고 눌린 성분의 것은 버리는 조립 `EvalImplicitTension` 을 생성기가 내게 했다(명세가 `IMPLICIT_TENSION = True`
로 고른다 — 지금은 인장만). `A = (h kd + h² k) Σ∇C∇Cᵀ + h² k Σ max(C, 0) ∇²C` 라 언제나 양의 준정부호다. 정의대로이고
(기준 조립과 1e-15), 같은 천 적분에서 정확판이 버티던 곳을 하나도 잃지 않으면서 쓸 수 있는 회가 60 중 31 에서 52 로 늘었다.
부드러운 천도 h = 1/120 에서는 모두, 1/60 에서도 굽힘 · 전단 판에 따라 끝까지 적분된다. 시간은 정확판의 1.003 배다.
미리 정한 규칙의 권고대로 인장의 기본 강성을 이 판으로 바꿨다(사용자 결정 2026-11-13).

**전단 GN 재비교 (2026-11-13)** — `docs/GATE_SHEAR_GN.md`. 인장을 `EvalImplicitTension` 으로 둔 채 전단의 정확판과 GN 판을
같은 천으로 다시 견줬다. 전단 GN 은 정확판이 버티던 곳을 하나도 잃지 않으면서 쓸 수 있는 회가 30 중 24 에서 28 로 늘었다.
늘어난 넷은 모두 부드러운 천에서 정확판 전단의 선형계가 깨지거나 발산한 회다. 결과 차이는 시간 간격을 반으로 줄인 차이의
9 % 안이고, 시간은 정확판의 0.99 – 1.02 배다. 미리 정한 규칙의 권고대로 전단의 기본 강성을 GN 판으로 바꿨다(사용자 결정 2026-11-13).

**시간 측정** — 같은 기계어가 실행 파일의 로드 주소에 따라 함수 하나만 1.6 – 3.3 배 느려진다. 그래서 시간 판정은
`tests/timing_bases.py` 가 고정 로드 주소 8 곳에 다시 링크해 잰 값의 최솟값으로 한다.

## 재현

```
run_gate.bat           # 생성 → 표본 → Mathematica 기준 → MSVC 빌드 → 게이트 → 시간 판정 (약 8.5 분)
run_gate.bat nomma     # Mathematica 없이(커밋된 results/mma_ref.txt 사용)
```

필요한 것: Python + SymPy 1.14, MSVC 2022, (선택) Mathematica `wolframscript`. 엔진 손 코드 비교 게이트는 엔진 리포가
있을 때만 돈다(엔진 코드는 이 리포에 싣지 않고 실행 때 뽑는다).

## 폴더

| 경로 | 무엇 |
|---|---|
| `codegen/` | 생성기 `sxgen.py` 와 명세 `specs/` |
| `src/` | 조건 → 힘 조립(`sx_condition_force.h`), 명세 → SXMesh 변수 · 스텐실 모으기(`sx_stencil.h`) |
| `third_party/sxmesh/` | SXMesh(이름으로 찾는 변수 컨테이너) 사본 |
| `model/` | 게이트 표본과 Mathematica 독립 기준 |
| `tests/` | 게이트 시험 |
| `results/` | 판정 숫자(`*.json`)와 기준 자료 |
| `docs/` | 판정문 |

## 관련 공개 도구

- **SymX** (RWTH Aachen, Apache-2.0) — 기호 에너지에서 미분 · 코드 생성 · 병렬 조립까지 하는 C++ 도구.
  이 프로젝트와는 **별개**이고 그 코드를 쓰지 않는다.
- **TinyAD** (MIT) — 작은 스텐실의 2차 전진 자동미분.

둘 다 **에너지 기반**이다. 이 프로젝트가 덧붙이려는 것은 (1) Baraff 식의 **조건 함수 · 감쇠** 명세,
(2) 근사(길이 고정 등)와 정확판을 같은 명세에서 만들어 **차이를 재는** 장치, (3) 명세 선언에서 바로 나오는
SoA 변수, (4) 의복 계층(봉제 · 패널 · 인체 충돌)이다.

## 라이선스

MIT (`LICENSE`). 엔진 cfmsDrape 에서 온 것은 범위 밖이다 — `third_party/sxmesh/`(SXMesh 사본) ·
`tests/engine_ref/SDEPreHeaders.h`(엔진 헤더의 typedef · enum · 매크로), 그리고 리포에 싣지 않는
`generated/engine_ref/`(빌드 때 뽑는 엔진 코드) · 엔진 DLL.
