# 엔진 손 코드 비교 게이트 — Macri 굽힘(`CheckBendingChange_Macri_SOA`)과 값 · 속도

엔진(cfmsDrape_dev, 커밋 378d20c 2026-10-26)의 Baraff 굽힘 경로는 손으로 짠 SoA 코드다. `__CalculateBending` 이
Baraff 모드에서 `CheckBendingChange_Macri_SOA` 를 부르면, 이 함수가 모든 힌지에 대해 θ · ∇θ · ∇²θ 를 SoA 배열로
만들고 끝에서 `_Input_F_and_J` 가 힘 · Jacobian 을 모은다. Bridson 모드(Oh 경로)도 같은 함수를 부르지만 θ 만 쓴다
(`bUseC2 = 0`). 사용자 기억(2026-11-12): Macri 의 소스를 고친 것이라 오류가 적을 것 — 그래서 값을 견줄 상대로
쓴다. **기준은 아니다**(CLAUDE.md 규약: 정확도 기준은 Mathematica 와 유한차분뿐).

판정 기준은 실행 전에 고정한다(볼트 규칙 §1, 2026-11-13). **판정: GO (2026-11-13)** — 아래 「결과」와 「판정 뒤 진단」.

## 무엇을 견주나

| | 엔진 손 코드 | 생성 코드 |
|---|---|---|
| 하는 일 | `_Input_b012`(모으기 · 변 · 법선) → θ(`atan2`) → `dnA/dnB`(틸드) → `dcos · dsin` → ∇C(가지별) → `_Input_CosD2 · _Input_SinD2` → ∇²C | SXMesh 에서 모으기 → `EvalStaged` |
| 내는 것 | 힌지마다 `_C`, `_dC_dxm[4][3]`, `_d2C_dxmdxn[4][4][3][3]`(12 × 12 전체) | 힌지마다 `C`, `∇C`(12), `∇²C`(위 삼각 78) |
| 빼는 것 | `_Input_F_and_J`(조립) — 엔진 풀이기의 희소 구조에 쓰는 것이라 옮기지 않고 빈 함수로 둔다 | 조립 |

`_Input_F_and_J` 는 옮기지 않으므로 이 게이트는 **미분까지**만 견준다. 조립은 조립 게이트(`docs/GATE_ASSEMBLY.md`)가
다뤘다.

## 옮기는 방법

- `tests/engine_ref/extract_engine.py` 가 엔진 리포에서 `git show 378d20c:<경로>` 로 원본을 읽어(엔진 리포는 읽기만),
  필요한 함수 정의를 **이름으로 찾아 바이트 그대로** `generated/engine_ref/` 에 쓴다. 생성물처럼 git 에 넣지 않는다
  (엔진 코드 사본을 이 리포에 싣지 않기 위해서). 엔진 리포가 없으면 이 게이트는 건너뛴다.
  - `SDECore.cpp`: `CheckBendingChange_Macri_SOA`, `_Input_b012`, `_Input_CosD2`, `_Input_SinD2`, 보조 함수
    (`_Tilde_Operation`, `_Dot31_M33_by_V31` 셋, `_TripleProduct31_V31_M33_W31`, `_Dot33_M33_by_N33`, `_Outer33_V31_W31`,
    `_Add33_M33_N33` 둘, `_AXPY33_M33_N33` 둘, `_Dot33_T333_M33_V3`, `_Dot33_T333_V3`, `_Factor33` 둘),
    `__StartCheckTEMPmem` · `__EndCheckTEMPmem`(빈 함수), `my_sse_crossproduct`, `BEST_SSE_*` 함수 포인터 선언,
    상수 블록(`MINUS` · `PLUS` · `X` · `_H` …), 실행 때 포인터를 넣는 대입문(`_DO_NOT_USE_SSE` 경로)
  - `SDECore.h`: 위 함수들의 선언(기본 인수가 여기 있다)
  - `SDEImplicitVariables.h`: `_SerialBendingVariables` 클래스 선언
  - `SDEmysse.cpp` · `SDEmysse.h`: 파일 전체(`my_SSE_*` 일반 반복문)
- 하네스(`tests/engine_ref/`)가 채우는 것: 최소한의 `SDECore` 껍데기(같은 멤버 이름), `SDEVector3D` · `SDEHalfEdge`
  최소형(꺼진 `if(0)` 블록이 컴파일되도록), `SDEPreHeaders.h` 대역(`dFloat = double`, 엔진과 같은 매크로
  `_DO_NOT_USE_SSE` · `_GLOBAL_IMPLICIT_BENDING_2ND_DERIVATIVE`), `_Input_F_and_J` 빈 함수.
- **`sbv` 할당은 엔진 `SetMem` 을 옮기지 않고 같은 배치로 새로 짠다.** 배열마다 따로 `new`, `_TEMP333[s][t][q] =
  _TEMP[10 + …]` 별칭은 엔진 그대로, `_dnA_dxm` · `_dnB_dxm` 은 엔진처럼 0 으로 시작(엔진은 쓰지 않는 블록이
  0 이라고 가정한다). 나머지도 0 으로 채운다.
- 엔진 빌드와 같게: `/O2`, `/fp:precise`(기본), 함수 포인터로 `my_SSE_*` 를 부름(인라인 안 됨). **단일 스레드** —
  엔진 빌드는 `OMP` 를 정의해 `_Input_b012` 의 모으기와 `atan2` 반복 둘을 병렬로 돌리지만(`SDEmysse.cpp` 는
  `#undef OMP` 라 나머지는 직렬), 생성 코드도 단일 스레드로 재므로 맞춘다. 엔진 코드 번역 단위에는 이 리포의
  `/W4 /WX` 를 걸지 않는다(엔진 코드를 고치지 않으므로).
- 좌표는 `DimensionFactor = 1`, m 단위. `iImplicitBendingForceType = SDEibfImplicitBaraff`, `bUseC2 = 1`.

### 꼭짓점 · 부호 대응 (읽어서 유도, E2 가 확인)

엔진 꼭짓점 `(H, I, J, K)` = 명세 `(h, i, j, k)`. 엔진은 `nA = (xJ − xI) × (xI − xH)`, `nB = (xJ − xK) × (xK − xI)`,
`e = xI − xJ` 라 명세의 `nA` · `nB` · `e` 와 셋 다 부호가 반대다. 그러면 `cos θ` 는 같고 `sin θ` 는 반대라
**θ_엔진 = −θ**. 그래서 ∇C · ∇²C 도 부호를 뒤집어 견준다(힘 `−kC∇C` · 강성 `−k(∇C∇Cᵀ + C∇²C)` 는 부호에 무관).
엔진 θ 에는 쉬는 각이 없다(Baraff 경로가 `_C = atan2(…)` 를 그대로 쓴다) — 견줄 때 명세 쪽은 `C + θ0` 를 쓴다.
이 대응이 E2 에서 틀리면 대응을 고쳐 다시 돌리고 여기에 적는다. 분류 기준(E3)은 바꾸지 않는다.

### 행렬 배치

엔진의 소비자 `_Input_F_and_J` 가 읽는 대로 `_dC_dxm[m][s] = ∂C/∂x_{m,s}`,
`_d2C_dxmdxn[m][n][s][t] = ∂²C/∂x_{m,s}∂x_{n,t}`(행 `3m + s`, 열 `3n + t`)로 읽어 12 × 12 전체를 견준다.
만드는 쪽의 배치가 이와 다르면 그것도 엔진의 오차로 센다 — 소비자가 그렇게 읽기 때문이다.

## 명제

| | 명제 |
|---|---|
| E1 | 하네스가 컴파일하는 엔진 코드가 378d20c 의 원본과 바이트 단위로 같다 |
| E2 | 엔진 θ 가 −θ(Mathematica)와 반올림 수준에서 같다 — 옮김 · 부호 대응이 맞다는 확인 |
| E3 | 엔진 ∇C · ∇²C 가 정확한 미분인지, 원문 §4.3 의 길이 고정 근사(`EvalFrozen`)인지, 어느 쪽도 아닌지 가른다 |
| E4 | 생성 코드가 같은 일(모든 힌지의 C · ∇C · ∇²C 를 내 놓기)을 엔진 손 코드보다 느리지 않게 한다 |

## 판정 기준 (실행 전 고정 · 2026-11-13)

상대오차 = 힌지마다 `max|Δ| / max|기준|`, 표본 전체의 최댓값(착수 게이트 G1 과 같은 잣대).

| 게이트 | 기준 | 표본 |
|---|---|---|
| E1 옮김 | `generated/engine_ref/` 의 파일마다 SHA-256 이 `tests/engine_ref/manifest.json`(378d20c 에서 뽑아 고정한 값)과 같고, 하네스가 그 파일을 고치지 않고 `#include` · 컴파일한다 | — |
| E2 θ | `max |θ_엔진 − (−(C_mma + θ0))| ≤ 1e-12` rad | `results/mma_ref.txt` 200 개 |
| E3 분류 (게이트 아님) | **정확** = ∇C · ∇²C 모두 Mathematica 대비 ≤ 1e-10. **길이 고정 근사** = 둘 다 `EvalFrozen` 대비 ≤ 1e-10. 둘 다 아니면 **어느 쪽도 아님** | 200 개 |
| E4 시간 | `min ns(생성) ≤ min ns(엔진)`. 생성 = SXMesh 에서 모으기 + `EvalStaged` + C · ∇C · ∇²C(91 개)를 힌지별 출력 배열에 저장. 엔진 = `CheckBendingChange_Macri_SOA` 한 번(`_Input_b012` 부터 ∇²C 까지)을 힌지 수로 나눈 값 | 격자 힌지 118,405 개(조립 게이트와 같은 G = 200), 단일 스레드, exe 안 5 회 중 최소, 고정 로드 주소 8 곳의 최솟값(`tests/timing_bases.py`) |
| 보고만 | E3 의 오차(Mathematica 대비 중앙값 · 최댓값, `|sin| > |cos|` 가지별 · 힌지 수), 엔진 ∇²C 의 비대칭, 1000 개에서 `Eval` 대비; E4 의 비, 저장하지 않는 `EvalStaged`, `EvalImplicit`(조립까지), 압축 명세, `EvalFrozen`, 작은 격자(G = 60, 약 1 만 힌지)의 비, 엔진 `sbv` 메모리 | — |

**E1 · E2 가 통과해야 비교가 성립한다.** 둘 중 하나가 실패하면 판정하지 않고 하네스를 고친다(NO-GO 가 아니다).
**그 위에서 E4 를 통과하면 GO** — 생성 코드가 엔진 Macri 손 코드를 대신할 수 있다(값은 착수 게이트에서
Mathematica 와 맞췄고, 속도는 손 코드 이상). E3 는 엔진 코드가 무엇을 계산하는지에 대한 분류라 GO 를 막지 않는다.

E4 를 「저장 포함」으로 정한 이유: 엔진은 결과를 SoA 배열에 남기고 조립이 다시 읽는다. 생성 코드는 그 자리에서
조립까지 이어갈 수 있어 저장이 필요 없지만, 같은 일을 하는지 다툼이 없도록 생성 쪽에 저장을 얹는다(생성 쪽에
불리한 쪽). 엔진은 `_Input_b012` 에서 속도 `_v0` 도 모은다(뒤의 감쇠용) — 엔진 설계의 비용이라 그대로 둔다.

## 읽어서 세운 예측 (판정에 쓰지 않는다)

- θ 는 단위 법선으로 계산하므로 E2 는 통과할 것이다.
- 미분은 **정규화하지 않은** `nA` · `nB` · `e` 로 한다. `_dnA_dxm = tilde(qA)` 이고, 여기에 `1/|nA|` 를 곱하는 블록은
  `if(0)` 로 꺼져 있다(주석 「이게 맞는 거 아닌가?」). 그러면 `dcos` 는 `nA · nB` 의 미분, `dsin` 은 `(nA × nB) · e` 의
  미분이라, 길이 고정 근사와도 코사인 가지에서 `|nA||nB|` 배, 사인 가지에서 `|nA||nB||e|` 배 크기가 다르고
  (변 1 cm 에서 약 1e-8 · 1e-10), 두 가지 사이에서 크기가 끊긴다. 가지 공식 `−dcos/sin` · `dsin/cos` 도 길이 고정
  근사의 `cos·dsin − sin·dcos` 와 같지 않다(정확한 미분에서만 같다). 예측은 **어느 쪽도 아님**이다.
- E3 가 「어느 쪽도 아님」이면 판정 뒤에 중간 배열(`dcos` · `dsin` · `d2cos` · `d2sin`)을 기호 미분과 견주는 진단을
  판정 밖으로 따로 한다.

## 읽기만 한 관찰 (재지 않았다)

`_Input_F_and_J` 는 옮기지 않으므로 아래는 코드를 읽은 것뿐이다.

- 힘이 `_C = atan2(…)` 를 그대로 써 쉬는 각 `Theta0` 를 빼지 않는다(`I_Bend_0 = I_Bend` 라 모든 힌지가 이 경로를 탄다).
- `bUseUpperJacobianMethod` 경로에서 블록을 고르는 조건이 `m == 0 && nn == 2` 처럼 국소 번호 `n` 이 아니라
  전역 노드 번호 `nn` 을 쓴다. 그 경로가 켜져 있으면 (0, 2) · (0, 3) · (1, 2) · (1, 3) · (2, 3) 블록은 노드 번호가
  우연히 2 · 3 일 때만 쓰인다. 첫 조건은 `m==0 & n== 1`(비트 연산 `&`, 결과는 같다).
- 감쇠는 `∂d/∂v = −kd ∇C∇Cᵀ`, `∂d/∂x` 는 `−kd Ċ ∇²C` 만(원문의 대칭 꼴).

## 결과 (2026-11-13) — GO

| 게이트 | 값 | 기준 | |
|---|---|---|---|
| E1 옮김 | 조각 32 개(파일 7 개) SHA-256 이 고정값과 같다. 엔진 코드 번역 단위의 경고는 쓰지 않는 지역 변수(C4101)뿐 | 일치 | PASS |
| E2 θ | `4.4e-16` rad (1000 개 `Eval` 대비 `8.9e-16`) — 부호 대응 θ_엔진 = −θ 가 맞았다 | ≤ 1e-12 | PASS |
| E3 분류 | **어느 쪽도 아님** — Mathematica · `EvalFrozen` 대비 ∇C · ∇²C 상대오차가 중앙값 · 최댓값 모두 1.00 | — | 분류 |
| E4 시간 | **0.272** (생성 + 저장 376 ns / 엔진 1386 ns), 중앙값 0.274, 로드 주소별 0.270 – 0.284 | ≤ 1 | PASS |

E1 · E2 위에서 E4 를 통과했으므로 **GO**: 생성 코드(`EvalStaged`)는 엔진 Macri 손 코드보다 3.7 배 빠르고, 값은
Mathematica 와 맞는다(착수 · 속도 게이트). 엔진 손 코드의 값은 기준과 맞지 않는다(E3, 아래 진단).

**E3 의 내용.** 상대오차 1.00 은 엔진 ∇C 가 정확값보다 훨씬 작다는 뜻이다. 크기 비 `max|∇C_엔진| / max|∇C|` 의
중앙값은 `|sin| > |cos|` 가지(엔진이 `dcos` 를 쓰는 가지, 1000 개 중 549)에서 4.2e-8, `|cos| ≥ |sin|` 가지(`dsin` 을
쓰는 가지, 451 개)에서 1.0e-10 이다. 예측의 「코사인 가지」는 `dcos` 를 쓰는 가지(`|sin| > |cos|`)를 뜻했고, 크기
자릿수(1e-8 · 1e-10)는 예측과 같다. 엔진 ∇²C 는 비대칭이다(`max|H − Hᵀ| / max|H|` = 2.0).

시간(고정 로드 주소 8 곳의 최솟값, 힌지 118,405 개, 단일 스레드, `results/gate_timing_engine.json`):

| | ns / 힌지 | 엔진 대비 |
|---|---|---|
| 엔진 `CheckBendingChange_Macri_SOA` | 1386 | 1 |
| 생성 `EvalStaged` + 저장 (E4) | 376 | 0.272 |
| 생성 `EvalStaged`, 저장 없음 | 381 | 0.275 |
| 생성 `EvalFrozen` | 262 | 0.189 |
| 압축 명세 `EvalStaged` + 저장 | 147 | 0.106 |
| 생성 `EvalImplicit`(조립 · 흩뿌리기까지) | 495 | 0.357 |
| 작은 격자(힌지 10,325): 엔진 1131, 생성 + 저장 380 | | 0.336 |

저장의 비용은 재는 잡음 안이다(376 · 381). 엔진 쪽 `sbv` 는 847 MB(힌지당 7.5 KB)다. 엔진이 할당하지만 이 함수들이
건드리지 않는 `_d2nA_dxmdxn` · `_d2nB_dxmdxn`(배열 864 개)은 할당하지 않았다. 생성 코드는 조립까지 해도(`EvalImplicit`)
엔진이 미분만 하는 시간의 0.36 배다.

## 판정 뒤 진단 (판정 밖, `tests/diag_engine_terms.cpp`)

판정문에 미리 적은 대로, E3 가 「어느 쪽도 아님」이라 엔진 중간 배열을 기호 미분과 견줬다. 엔진의 벡터
(`nA = (xJ − xI) × (xI − xH)`, `nB = (xJ − xK) × (xK − xI)`, `e = xI − xJ`, 정규화하지 않음)로
`P = nA · nB`, `Q = (nA × nB) · e` 를 생성기로 미분한다(`tests/diag_engine_terms_spec.py`). 곱의 법칙의 항은 짝 벡터를
스텐실 값(고정값)으로 넣어 하나씩 뗀다(`tests/diag_engine_split_spec.py`): `∇P = ∇PA + ∇PB`,
`∇Q = ∇QA + ∇QB + ∇QE`(A = `nA` 만, B = `nB` 만, E = `e` 만 x 의 함수). 1000 힌지, `results/diag_engine_terms.json`.

| | 결과 |
|---|---|
| 가지 공식(#b.9 · #b.10) | 엔진 자신의 `dcos` · `dsin` · `d2cos` · `d2sin` 에 읽은 공식을 걸면 엔진 ∇C · ∇²C 와 비트까지 같다(차 0) — 공식은 바르게 읽었다 |
| `_dcosth_dxm` | `∇P` 와 다르다(상대 2.7). 꼭짓점 블록마다 부호 조합을 맞추면 **`dcos = ∇PA − ∇PB`**(오차 0): `∂nB` 항의 부호가 반대다. H 블록은 `∂nB` 항이 없어 맞다 |
| `_dsinth_dxm` | `∇Q` 와 다르다(1.4). **날개 꼭짓점 H · K 에서 0, 힌지 꼭짓점 I · J 에서 `∓(nA × nB)` 뿐**(오차 1e-16 자릿수). `nA × nB` 는 변 `e` 와 나란하다. 정확한 `∂Q/∂x` 의 법선 항(`∇QA` · `∇QB`)이 없다 |
| `_d2costh` · `_d2sinth` | `∇²P` · `∇²Q` 와 그대로도, 블록 안 전치로도 다르다(2.0 · 2.35) |

읽어서 맞춘 원인(측정은 위 표뿐이다):

- `dcos` 의 B 항은 `_Dot31_M33_by_V31(…, sbv._nA, sbv._dnB_dxm[m])`(열 합 = 행렬의 전치)로 만든다. `_dnB_dxm = tilde(qB)` 는
  반대칭이라 전치가 곧 부호 반대다. A 항은 `_Dot31_M33_by_V31(…, sbv._dnA_dxm[m], sbv._nB)`(행 합)라 맞다.
- `dsin` 은 `(−ñB · dnA + ñA · dnB) · e` 로 법선 미분 행렬을 `e` 에 곱한다. 날개 꼭짓점에서는 `qA[H] = −e`,
  `qB[K] = e` 라 `tilde(q) · e = 0` 이 되고, 힌지 꼭짓점에서는 `e` 방향 성분만 남는다. 기울기라면 행렬의 전치를 짝 벡터
  (`nB × e` 등)에 곱해야 한다.

뜻(위 측정에서 따라 나오는 것):

- `|cos| ≥ |sin|`, 곧 `|θ| ≤ 45°`(거의 평평한 천, 가장 흔한 경우)와 `|θ| ≥ 135°` 에서 엔진 Baraff 경로의 굽힘 힘
  `−kC∇C` 는 **힌지 꼭짓점 I · J 에만, 변 방향으로만** 걸리고 크기는 정확값의 약 1e-10 배다. 날개 꼭짓점에는 힘이 없어
  굽힘을 되돌리지 못한다.
- `45° < |θ| < 135°` 에서는 크기가 약 4e-8 배이고 방향도 `∂nB` 항의 부호 탓에 틀린다. 두 가지의 크기 비 중앙값이
  약 400 배 다르므로 가지가 바뀌는 θ = ±45° · ±135° 에서 힘의 크기가 불연속으로 뛴다(한 힌지를 가로질러 잰 것은
  아니다).
- **엔진 기본값은 이 미분을 쓰지 않는다.** 기본 굽힘은 Bridson 이고(`SDECore.cpp:1047`
  `iImplicitBendingForceType = SDEibfImplicitBridson`), 그 경로는 이 함수에서 θ 만 쓴다(`bUseC2 = 0`). θ 는 E2 에서 맞다.
  Baraff 경로는 `SDEVC.cpp:1627` 로 따로 골라야 탄다.
- 사용자 기억(「Macri 의 소스를 고친 것이라 오류가 적을 것」, 2026-11-12)과 달리, 378d20c 의 Baraff 경로에서
  맞는 것은 θ 뿐이고 ∇θ · ∇²θ 는 기준과 맞지 않는다.

## 재현

```
run_gate.bat
```

엔진 리포 위치는 환경 변수 `CFMS_DRAPE_DIR`(없으면 `..\cfmsDrape_dev`). 숫자는 `results/gate_engine.json` 과
`results/gate_timing_engine.json`. 진단은 판정 밖이라 `run_gate.bat` 이 돌리지 않는다 — 생성기로
`tests/diag_engine_terms_spec.py` · `tests/diag_engine_split_spec.py` 를 `generated/` 에 낸 뒤
`tests/diag_engine_terms.cpp` 를 `engine_macri.obj` · `SDEmysse.obj` · `SXMesh.obj` 와 링크해 돌린다.
