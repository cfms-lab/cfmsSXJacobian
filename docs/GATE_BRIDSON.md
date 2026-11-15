# Bridson 굽힘 게이트 — 조건 함수 명세, 그리고 엔진 손 코드(`GetImpBending_Oh_SOA`)의 힘 · Jacobian · 속도

엔진(cfmsDrape_dev, 커밋 378d20c)의 **기본** 굽힘은 Bridson 이다(`SDECore.cpp:1047`
`iImplicitBendingForceType = SDEibfImplicitBridson`). `__CalculateBending` 이 그 갈래에서 `bUseC2 = 0` 으로
`CheckBendingChange_Macri_SOA` 를 불러 θ 만 만들고(엔진 비교 게이트 E2 에서 맞음), 이어서 `GetImpBending_Oh_SOA` 가
힌지마다 힘과 Jacobian 블록을 만든다. 사용자(2026-11-12): 이 함수는 손으로 짠 코드라 **정확도의 기준이 아니다.**
정확도 기준은 Mathematica 와 유한차분뿐이다(CLAUDE.md 규약). 그래서 이 게이트는

1. Bridson 굽힘을 **조건 함수 한 줄**로 쓰고 생성 코드가 기준과 맞는지 판정한다(B1 – B4),
2. 그 정확한 값에 엔진 손 코드의 힘 · Jacobian 을 대어 무엇을 계산하는지 가른다(B6 · B7),
3. 같은 일을 하는 시간을 견준다(B8).

판정 기준은 실행 전에 고정한다(볼트 규칙 §1, 2026-11-13). **판정: 명세 GO · 엔진 비교 NO-GO (2026-11-13)** — 아래 「결과」.

## 조건 함수로 쓴 Bridson 굽힘

Bridson, Marino, Fedkiw (2003) §4 의 탄성 굽힘 힘은 힌지 꼭짓점마다

  F_m = −k w (sin(θ/2) − sin(θ0/2)) u_m

이고 네 벡터 u_m 은 이면각의 기울기 ∇_m θ 다(부호는 θ 의 방향 규약에 따른다). w 는 힌지 가중치로, 엔진은
discrete shells 의 `3|ē|²/(A1 + A2)` 를 쓴다(`SDECore.cpp` Oh 함수 안 주석 「handoff #3」).
w 를 상수로 보고 **θ0 = 0** 이면 이 힘은 에너지

  E = 2 k w (1 − cos(θ/2)) = 4 k w sin²(θ/4)

의 −∇E 다. 그래서 조건 함수를

  **C = 2 √(2 w̄) sin(θ/4)**

로 두면 ½ k C² = E 이고, 생성기가 내는 힘과 강성이 그대로 Bridson 의 힘과 그 정확한 Jacobian 이 된다.

  f = −k C ∇C = −k w̄ sin(θ/2) ∇θ
  K = −k (∇C∇Cᵀ + C ∇²C) = −k w̄ (½ cos(θ/2) ∇θ∇θᵀ + sin(θ/2) ∇²θ)

- **w̄ 는 힌지마다 상수(스텐실 값 `pBendW`)** 로 둔다. 시뮬레이션에서는 쉬는 모양에서 정한다 — 그러면 에너지가 보존적이고
  K 가 그 정확한 Hessian 이다. 엔진은 매 단계 **지금의** 넓이로 w 를 다시 재고(`Length0` 은 쉬는 변 길이) 미분하지 않는다.
  그러면 엔진의 힘 `−k w(x) sin(θ/2) ∇θ` 는 어떤 에너지의 기울기도 아니다. 견줄 때는 w̄ 를 그 상태에서 엔진이 쓴 w 와 같게
  둔다(B6 · B7).
- θ 는 압축 명세(`baraff_bend_compact.py`)와 같은 `atan2((nA × nB)·e, (nA·nB)|e|)` 로 쓴다(정규화가 필요 없다).
- **θ0 ≠ 0 은 다루지 않는다.** 엔진 기본값이 `Theta0 = 0` 이다. θ0 ≠ 0 이면 E 가 조건의 제곱 꼴이 아니다. 엔진은 쉬는 각 항을
  `2·sin(θ0/2)` 로 써서 원문(`sin(θ0/2)`)과 다르다 — 읽기만 한 것이고 재지 않는다.

명세는 `codegen/specs/bridson_bend.py`, 독립 기준은 `model/bridson_bend_ref.wl`(Mathematica).

## 엔진 코드 — 읽은 것

- 꼭짓점: 엔진 힌지 `_vid = (H, I, J, K)` = 명세 `(h, i, j, k)`. Oh 함수의 지역 번호는 `nn = (vid0, vid3, vid2, vid1)`
  = `(h, k, j, i)`.
- 가중치 `_cb = 3·Length0² / (0.5/_nAim + 0.5/_nBim)`, `_nAim = 1/|nA|`(그래서 분모 = A1 + A2, 지금 넓이).
  강성 `_kb = BE_U × GlobalBendingFactor`.
- 힘 `f_m = α u_m`, `α = −_cb·_kb·(sin(θ/2) − 2 sin(θ0/2))`. `u_h = |e| n̂A/|nA|`, `u_k = |e| n̂B/|nB|`,
  `u_j = (dot0332 n̂A/|nA| + dot1332 n̂B/|nB|)/|e|`, `u_i = (dot2032 n̂A/|nA| + dot2132 n̂B/|nB|)/|e|` — Bridson 의 u 와 같은 꼴.
  `pCGSolver->pBendingF[12i + 3·(지역 번호) + 성분]` 에 쓴다.
- Jacobian(기본값 `bUseUpperJacobianMethod = bUseUDJSerialJacobian = true`): 힌지 소유 슬롯 `pBendingJ[6i + 1 … 4]` 에
  블록 넷 — `(h, j) = β·dot0332·n̂An̂Aᵀ`, `(h, i) = β·dot2032·n̂An̂Aᵀ`, `(k, j) = γ·dot1332·n̂Bn̂Bᵀ`, `(k, i) = γ·dot2132·n̂Bn̂Bᵀ`,
  `β = −_cb·_kb/|nA|²`, `γ = −_cb·_kb/|nB|²`. `(h, k)` · `(i, j)` 블록은 쓰지 않는다(0). 대각 블록은 풀이기 규약
  `∂f_m/∂x_m = −Σ_{n≠m} J_mn` 으로 정해진다(`SDECore.cpp:1673` 주석. `SDESparseMatrix.cpp`
  `_Jacobian_VectorProduct` 가 쌍 블록마다 `J (v_j − v_i)` 를 행 i 에 더하고 행 j 에서 뺀다 — 3 × 3 블록 전체).
  루프 뒤 직렬 환원이 넷을 노드쌍 슬롯 `pSerialJacobianJ` 에 더한다.
- 감쇠는 기본값에서 꺼져 있다(`bBendDamping = false`). 퇴화 힌지 가드(넓이 · 변 길이)는 이 표본에서 걸리지 않는다.

### 견주는 방법

- 엔진 힘 12 개를 명세 순서 `(h, i, j, k)` 로 바꿔 생성 코드의 `−k C ∇C` 와 견준다(w̄ = 엔진의 `_cb`, k = 1).
- 엔진 Jacobian 은 블록 넷과 그 전치, 위 규약의 대각으로 12 × 12 전체를 만들어 정확한 K 와 견준다. 풀이기가 이 행렬로
  선형계를 세우므로 이것이 엔진이 쓰는 ∂f/∂x 다.

## 읽어서 세운 예측 (판정에 쓰지 않는다)

- **B6 힘은 정확**할 것이다. 엔진의 u 는 Bridson 의 u 그대로이고, Bridson 의 u 는 ∇θ 이며, `−k C∇C` 는 θ 의 부호 규약에
  무관하다(C∇C 가 θ 에 대해 짝함수).
- **B7 Jacobian 은 어느 쪽도 아닐** 것이다. 엔진 J = `−k w̄ M` 이고 M 은 ∇θ∇θᵀ 에서 (1) 같은 법선끼리의 항만 남기고
  (n̂An̂Bᵀ 교차 항을 버림), (2) `(i, j)` 블록을 쓰지 않아 그 몫이 행 합 규약으로 대각에 옮겨 간 것이다. 정확한 K 와 견주면
  (3) 계수가 `½ cos(θ/2)` 대신 1 — 평평할 때 **두 배**, (4) `sin(θ/2) ∇²θ` 항이 없다.
  이 예측 행렬을 시험이 기하에서 따로 만들어 엔진 J 와 견준다(읽은 대로인지 확인, 보고).
- B8 은 가늠이 안 된다. 엔진은 θ 단계(`_Input_b012` 와 `atan2`)와 짧은 루프뿐이라 싸고, 생성 코드는 정확한 ∇²C 와
  식 (6) 조립까지 한다.

## 옮기는 방법

엔진 비교 게이트(`docs/GATE_ENGINE.md`)의 하네스를 넓힌다.

- `tests/engine_ref/extract_engine.py` 가 더 뽑는 것(378d20c, 바이트 그대로, `generated/engine_ref/` — git 밖):
  `SDECore.cpp` 의 `s_BendSkip` 선언과 `GetImpBending_Oh_SOA`, 전역 `OuterProduct`, `SDEMatrix3f.cpp` 의
  `Set(double*, double)` · `AXPY(double, double[3][3], double)`. 기존 조각 32 개의 SHA-256 은 그대로여야 한다.
- 하네스(`tests/engine_ref/engine_macri.cpp`)가 채우는 것: `SDEVector3D`(엔진처럼 `x, y, z` 와 `Data[3]` 공용체,
  `dFloat × 벡터`, `+=`, `Set`), `SDEMatrix3f`(`Data[3][3]`, `_USE_DYNAMIC_SMatrix3f` 꺼짐 — 엔진 기본), 재질
  `BE_U` · `BV` · `Theta0`, 반변 `Length0`, 천 `p2DNodeMass` · `pSewing`, 풀이기 껍데기(플래그와 배열),
  SDECore 멤버 몇 개. `_USE_MATH_DEFINES`(엔진 `SDEPreHeaders.h:142`).
- 설정: `iImplicitBendingForceType = SDEibfImplicitBridson`, `bUseC2 = 0`(`__CalculateBending` 그대로 두 함수를 차례로
  부른다), `Theta0 = 0`, `BE_U = 1`, `GlobalBendingFactor = 1`, `bBendDamping = false`, `bCalculateSewing = false`,
  `bUseUpperJacobianMethod = true`, `bUseUDJSerialJacobian = true`(기본). `piBendingJ` 는 엔진처럼 힌지의 노드쌍을
  공유 쌍 슬롯에 잇는다(노드 번호 작은 쪽이 행).
- 표본에는 쉬는 모양이 없으므로 `Length0` = 표본의 |e| 로 둔다(B6 · B7 은 w̄ 를 엔진의 w 에 맞추므로 값에 무관하다).

## 명제

| | 명제 |
|---|---|
| B1 – B4 | 조건 함수 명세 하나로 Bridson 굽힘의 정확한 C · ∇C · ∇²C 와 식 (6) 조립이 나온다 |
| B5 | 하네스가 컴파일하는 엔진 코드가 378d20c 원본과 바이트 단위로 같고, 엔진 θ 가 맞다(전제) |
| B6 | 엔진 힘이 정확한 Bridson 힘 `−k C∇C` 인지 가른다 |
| B7 | 엔진 Jacobian 이 정확한 K 인지, Gauss–Newton 꼴 `−k ∇C∇Cᵀ` 인지, 어느 쪽도 아닌지 가른다 |
| B8 | 생성 코드가 같은 일(힌지마다 힘과 강성 내 놓기)을 엔진 손 코드보다 느리지 않게 한다 |

## 판정 기준 (실행 전 고정 · 2026-11-13)

상대오차 = 힌지마다 `max|Δ| / max|기준|`, 표본 전체의 최댓값. 표본은 착수 게이트의 `results/hinges.txt`(1000 개,
처음 넷은 θ = 0 · ±90° · 150°). 표본 w̄ = `3|e|² / (½|nA| + ½|nB|)`(지금 모양 — 표본은 어느 값이든 된다).

| 게이트 | 기준 | 표본 |
|---|---|---|
| B1 Mathematica | `Eval` · `EvalStaged` 둘 다 `|ΔC|/√w̄ ≤ 1e-12`, ∇C · ∇²C rel ≤ 1e-10 | `results/mma_ref_bridson.txt` 200 개 |
| B2 유한차분 | `Eval` · `EvalStaged` 둘 다 ∇C(C 의 중심차분) · ∇²C(∇C 의 중심차분) rel ≤ 1e-6. 간격 = |e| × 1e-6, 분모 = 실제 간격 `xp − xm` | 1000 개 |
| B3 단계판 대 평탄판 | `|ΔC|/√w̄`, ∇C · ∇²C rel 모두 ≤ 1e-12 | 1000 개 |
| B4 조립 | `EvalImplicit` 대 `AssembleImplicit`(검산용 조립): A · b · f + d rel 모두 ≤ 1e-12, k = 1, kd = 0.1, h = 0.01 | 1000 개 |
| B5 옮김 (전제) | 조각마다 SHA-256 이 `tests/engine_ref/manifest.json` 고정값과 같고(새 조각 포함, 기존 32 개 그대로), `max|θ_엔진 − (−θ_mma)| ≤ 1e-12` | 200 개 |
| B6 엔진 힘 (분류) | **정확** = 모든 표본에서 `max|f_엔진 − f| / max|f| ≤ 1e-10`(f = `Eval` 의 `−kC∇C`, w̄ = 엔진 `_cb`, k = 1), 아니면 **다름**. 분모는 첫 실행 전에 `k w̄ max|∇θ|` 로 바꿨다(아래 「측정 방법 변경」) | 1000 개 |
| B7 엔진 Jacobian (분류) | **정확** = 모든 표본에서 엔진 J 대 K rel ≤ 1e-10, **Gauss–Newton** = 대 `−k∇C∇Cᵀ` rel ≤ 1e-10, 둘 다 아니면 **어느 쪽도 아님** | 1000 개 |
| B8 시간 | `min ns(생성) ≤ min ns(엔진)`. 생성 = SXMesh 에서 모으기 + `EvalImplicit` + 출력(C · f + d · A · b, 103 개)을 힌지별 배열에 저장. 엔진 = `CheckBendingChange_Macri_SOA`(`bUseC2 = 0`) + `GetImpBending_Oh_SOA` 를 `bUseUDJSerialJacobian = false` 로(힌지 소유 슬롯에 쓰기까지, 쌍 슬롯 환원 없음) | 격자 힌지 118,405 개(G = 200), 단일 스레드, exe 안 5 회 중 최소, 고정 로드 주소 8 곳의 최솟값(`tests/timing_bases.py bridson`) |
| 보고만 | B6 의 오차 중앙값 · 최댓값, 운동량 `|Σf|`; B7 의 K 대비 중앙값 · 최댓값(전체 · `|θ| < 10°`), 예측 행렬 대비, 오차 나눔(K → `−k w̄ ½cos(θ/2)∇θ∇θᵀ` → 계수 1 → 교차 항 버림 → `(i, j)` 옮김), ∇θ 방향 강성 비 `∇θᵀJ∇θ / ∇θᵀK∇θ`, 엔진 힘을 유한차분한 Jacobian(지금 넓이로 다시 잰 w 의 미분 포함) 대 엔진 J · 대 K; B8 의 엔진 기본값(쌍 슬롯 환원 포함), 생성 `EvalStaged` + 저장, 생성 `EvalImplicit` + 노드 배열로 흩뿌리기, 작은 격자(G = 60) | — |

**B5 가 통과해야 비교가 성립한다** — 실패하면 판정하지 않고 하네스를 고친다(NO-GO 가 아니다).
**B1 – B4 를 모두 통과하면 명세 GO** — Bridson 굽힘을 조건 함수 한 줄로 쓸 수 있다.
**그 위에서 B8 을 통과하면 GO** — 생성 코드가 엔진 Bridson 손 코드를 대신할 수 있다(값은 정확하고 속도는 손 코드 이상).
B6 · B7 은 엔진 코드가 무엇을 계산하는지에 대한 분류라 GO 를 막지 않는다.

B8 을 이렇게 정한 이유:

- 엔진 쪽은 기본값(쌍 슬롯 환원 포함) 대신 **환원을 뺀** 설정으로 잰다 — 생성 쪽도 힌지별 배열에 쓰고 끝나므로 같은 일이고,
  엔진에 유리한 쪽이다. 기본값은 보고한다.
- 생성 쪽은 엔진보다 많이 낸다(감쇠 · 우변 b · 위 삼각 78 칸 전체, 엔진은 힘과 3 × 3 블록 넷). 생성 쪽에 불리한 쪽이다.
- 엔진은 `_Input_b012` 에서 속도도 모으고, 앞 단계(`CheckBendingChange_Macri_SOA`)의 SoA 배열을 거친다 — 엔진 설계의
  비용이라 그대로 둔다. `sbv` 는 엔진처럼 모든 배열을 할당한다(약 850 MB).

## 측정 방법 변경 — B6 의 척도 (첫 실행 전, 2026-11-13)

판정문을 고정한 뒤 시험을 짜면서 B6 의 상대오차 분모를 `max|f|` 에서 **`k w̄ max|∇θ|`**(θ = π 에서의 힘 크기)로 바꿨다.
평평한 표본 0 은 θ = 0 이라 f = 0 이고, 두 쪽의 θ 가 반올림(1e-16 rad)만큼 달라도 `max|f|` 로 나누면 상대오차가 1
자릿수가 된다. 척도는 θ 에 무관한 힘의 크기라 다른 표본의 판정은 사실상 바뀌지 않는다. `max|f|` 대비 값도 보고한다.
실행 결과 `max|f|` 대비 값은 표본 0 에서 1.9e-2, 나머지 999 개의 최댓값은 9.3e-15 였다 — 옛 분모였다면 표본 0 하나 때문에
「다름」이 됐을 것이다.
분류 기준(1e-10)은 바꾸지 않았다.

## 결과 (2026-11-13) — 명세 GO, 엔진 비교 NO-GO (B8)

| 게이트 | 값 | 기준 | |
|---|---|---|---|
| B1 Mathematica | `Eval`: `|ΔC|/√w̄` 5.0e-16, ∇C 7.6e-16, ∇²C 3.4e-15 · `EvalStaged`: 5.0e-16, 9.5e-16, 2.5e-15 | 1e-12 · 1e-10 | PASS |
| B2 유한차분 | `Eval` ∇C 7.2e-10 · ∇²C 5.4e-10, `EvalStaged` 7.2e-10 · 4.1e-10 | 1e-6 | PASS |
| B3 단계판 대 평탄판 | 5.2e-16 · 9.2e-16 · 3.1e-15 | 1e-12 | PASS |
| B4 조립 | A 1.9e-15, b 7.8e-14, f + d 4.6e-14 (중앙값 4.6e-16 · 5.2e-16 · 3.7e-16) | 1e-12 | PASS |
| B5 옮김 (전제) | 조각 37 개(파일 9 개) SHA-256 일치 — 기존 32 개 그대로, 새 5 개. θ 4.4e-16 rad. 퇴화 가드에 걸린 힌지 0 | 일치 · 1e-12 | PASS |
| B6 엔진 힘 | **정확** — 척도 대비 중앙 2.2e-16 · 최댓값 1.1e-15 (`max|f|` 대비 중앙 3.8e-16, 표본 0 을 뺀 최댓값 9.3e-15), 운동량 `|Σf|` 5.9e-16 | 1e-10 | 분류 |
| B7 엔진 Jacobian | **어느 쪽도 아님** — K 대비 중앙 1.15 · 최댓값 2.41, `|θ| < 10°`(56 개) 중앙 0.75 · 최댓값 1.48, `−k∇C∇Cᵀ` 대비 중앙 1.25 | 1e-10 | 분류 |
| B8 시간 | **3.270** (생성 333 ns / 엔진 102 ns), 중앙값 3.208, 로드 주소별 3.14 – 3.29 | ≤ 1 | **FAIL** |

B1 – B4 를 모두 통과했으므로 **명세 GO** — Bridson 굽힘은 조건 함수 한 줄 `C = 2√(2w̄) sin(θ/4)` 로 쓸 수 있고,
생성기는 고치지 않았다(`sin` · `cos` 를 `std::` 로 내는 출력만 더했다 — 굽힘 착수판 해시 그대로, S3 통과).
B5 위에서 B8 을 넘지 못했으므로 **엔진 비교는 NO-GO** — 정확한 강성까지 내는 생성 코드는 엔진 Bridson 손 코드보다
3.3 배 느리다.

**B6 — 엔진 힘은 정확하다.** 엔진의 u 는 −∇θ(명세의 θ, 1000 개 모두, 1.1e-15 안)이고 엔진 θ 는 −θ 라, 힘은
`−k w̄ sin(θ/2) ∇θ` 그대로다. 정확하지 않을 것이라는 사용자 기억과 달리 **힘은 맞는다**(w 를 그 상태의 값으로 둔 Bridson 힘).

**B7 — 엔진 Jacobian 은 읽은 대로다.** 시험이 기하에서 따로 만든 예측 행렬(같은 법선 항만, `(i, j)` 블록을 대각으로 옮김,
계수 1)과 2.1e-15 안에서 같다. 정확한 K 와의 차이를 차례로 나누면(`|K|` 대비 중앙값 · 최댓값):

| 단계 | 중앙값 | 최댓값 |
|---|---|---|
| ∇²θ 항 `sin(θ/2) ∇²θ` 를 버림 | 0.97 | 1.66 |
| 계수 `½cos(θ/2)` → 1 | 0.79 | 2.32 |
| 교차 항 `n̂An̂Bᵀ` 를 버림 | 0.69 | 2.25 |
| `(i, j)` 블록을 대각으로 옮김 | 0.15 | 1.00 |

∇θ 방향(굽힘 모드)의 강성 비 `∇θᵀJ∇θ / ∇θᵀK∇θ` 는 전체 중앙값 1.47, `|θ| < 10°` 에서 중앙값 1.15(0.96 – 1.69),
평평한 표본 0 에서 1.06 이다. 평평할 때 계수 두 배와 교차 항 버림이 굽힘 모드에서는 거의 상쇄된다 — 행렬 전체로는
O(1) 로 다르다.

**엔진 힘의 참 Jacobian.** 엔진 힘을 유한차분한 Jacobian 은 `K + f ∇wᵀ / w` 와 4.5e-10 안에서 같다(w = 지금 넓이로 매 단계
다시 잰 가중치). w 미분 항의 크기는 `|K|` 대비 중앙값 0.47 · 최댓값 1.03 — 평평할 때는 f = 0 이라 사라지지만 굽은 힌지에서는
작지 않다. 엔진 J 는 이 항도, 정확한 K 도 아니다(유한차분 대비 2.97).

**B8 — 시간.** 고정 로드 주소 8 곳의 최솟값(힌지 118,405 개, 단일 스레드, `results/gate_timing_bridson.json`):

| | ns / 힌지 | 엔진 대비 |
|---|---|---|
| 엔진 Bridson(Macri θ + Oh, 쌍 슬롯 환원 없음) | 102 | 1 |
| 엔진 기본값(쌍 슬롯 237,606 개로 환원) | 118 | 1.16 |
| 생성 `EvalImplicit` + 저장 (B8) | 333 | 3.27 |
| 생성 `EvalImplicit` + 노드 배열로 흩뿌리기 | 335 | 3.29 |
| 생성 `EvalStaged` + 저장(C · ∇C · ∇²C 만) | 168 | 1.65 |
| 작은 격자(힌지 10,325): 엔진 82, 생성 329 | | 4.04 |

엔진은 θ 단계에서 `_Input_b012` 의 SoA 배열 몇 개와 `atan2` 만 쓰고(Baraff 경로의 ∇²C 는 건너뜀), Oh 루프는 Bridson 의
닫힌 꼴 u 와 3 × 3 블록 넷이라 싸다. 생성 코드는 ∇²θ 까지 든 정확한 강성과 감쇠 · 우변을 낸다. 미분만(`EvalStaged`) 해도
엔진의 1.6 배다. 같은 코드를 일곱 번 재니(전체 `run_gate.bat` 다섯 번 포함) B8 은 3.19 – 3.32, `EvalStaged` 는 168 – 199 ns 였다 —
판정은 어느 쪽이든 같다. 표는 판정한 실행(4f0da39 에 커밋한 `results/gate_timing_bridson.json`)의 값이고, 지금 커밋된
json 은 그 뒤 마지막 전체 실행의 값이다. 엔진 메모리는 926 MB(Macri 의 `sbv` 전체 + 힌지 소유 슬롯).

## 재현

```
run_gate.bat
```

명세 쪽(B1 – B4)은 엔진 리포가 없어도 돈다. 엔진 쪽(B5 – B8)은 엔진 리포(`CFMS_DRAPE_DIR`, 없으면 `..\cfmsDrape_dev`)가
있을 때만 돈다. 숫자는 `results/gate_bridson.json` 과 `results/gate_timing_bridson.json`.
