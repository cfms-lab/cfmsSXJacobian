# Provot 전단 게이트 — 엔진 기본 전단(대각선 늘림)을 명세 한 줄로 쓴다

상태: **GO (2026-11-13) — V0 – V7 모두 통과.** 기준은 실행 전에 고정했고(사용자가 이 기준 · b 는 스텐실 값으로 진행을 정함),
아래 「결과」는 그 뒤에 적었다. 엔진 연결 L3
(`docs/PLAN_ENGINE_LINK.md`)의 첫 단계다 — 계획대로 Provot 전단 명세는 엔진에 넣기 전에 이 리포에서 Mathematica · 유한차분
게이트를 먼저 거친다.

## 무엇을 보이나

엔진 cfmsDrape 의 기본 전단은 Baraff–Witkin 전단(`w_u · w_v`, `baraff_shear.py`)이 아니라 Provot 식 대각선 늘림이다
(사용자 결정 2026-11-13: 엔진 모형을 명세로 새로 쓴다). 이 게이트는 그 모형을 조건 함수 한 줄로 쓴 명세
`codegen/specs/provot_shear.py` 에서 맞는 `C · ∇C · ∇²C` 와 조립이 나오는지를 본다. 엔진 값과의 비교는 여기서 하지 않는다(L3).

## 엔진에서 읽은 것 (378d20c `SDECore.cpp` `_GetImplicitShear_Provot_WithInvariant`, 읽기만 · 판정에 쓰지 않는다)

- 조건 `C₊ = α(|w_u + w_v| − b_su)`, `C₋ = α(|w_u − w_v| − b_sv)`. α = 쉬는 넓이(`Invariant._alpha`), `w_u` · `w_v` 는 인장과 같다.
- `b_su = b_sv = 1.41421356237f` — √2 를 **float 리터럴**로 둔다(`SDECore.cpp:1302`). double 로는 1.4142135381698608 이라 √2 와
  2.4e-8 다르다.
- 강성 `k = SE × GlobalShearFactor / α`, 감쇠 `kd = SV × GlobalShearFactor / α`. 그래서 요소 에너지가 넓이에 비례한다.
- 힘 `f_m = −k(C₊∇C₊ + C₋∇C₋) − kd(Ċ₊∇C₊ + Ċ₋∇C₋)`. 감쇠는 Baraff 식 그대로다.
- `J`(기본 플래그 `bUseC2InPositiveDefinite` · `bPositiveDefiniteOnly`): `−k Σ ∇C∇Cᵀ` 에 성분마다 `C > 0` 일 때만 `−k C∇²C`,
  `Ċ > 0` 일 때만 `−kd Ċ∇²C` 를 더한다. `V = −kd Σ ∇C∇Cᵀ`. 곧 이 리포의 인장 C > 0 판(`EvalImplicitTension`)과 같은 처방이고,
  감쇠의 `Ċ∇²C` 몫만 더 있다 — 인장과 같은 차이다. 그 몫을 어떻게 다룰지는 L3 판정문에서 정한다.
- 비선형 재질 갈래(`bMembraneNonlinear`) · 전단 다항식 갈래(`bUseShearTangentPolynomial`)는 기본 꺼짐이고 이 명세 밖이다.
- 퇴화 가드: α ≤ 1e-8 이면 요소를 건너뛰고, 길이 ≤ 1e-12 이면 방향을 0 으로 둔다. 명세에는 가드가 없다(엔진 갈래가 맡는다).

## 명세

```
C = √a (|w_u + w_v| − b₊,  |w_u − w_v| − b₋),   a = uv 넓이
```

- 삼각형 스텐실 · 성분 2. 스텐실 값 = 쉬는 uv 6 칸(`pUV`) + **쉬는 대각선 길이 2 칸(`pDiag` = b₊ b₋)**.
- `√a` 를 곱하고 강성을 `k = SE` 로 두면 엔진의 `C = α(…)` · `k = SE / α` 와 힘 · 에너지가 같다 — 인장 명세와 같은 처방이다.
- `IMPLICIT_TENSION = True` — `EvalImplicitTension` 도 낸다. 성분마다 `C = √a(|d| − b)` 이고 `d = w_u ± w_v` 가 좌표의 1 차식이라
  `∇²C` 가 준정부호다(인장과 같은 근거).
- 두 대각선 `w_u ± w_v` 를 `ctx.stage` 로 표시한다.

**정할 것 하나 — b 를 스텐실 값으로 받는다(권고).** 상수 √2 로 명세에 박으면 엔진의 float 리터럴과 2.4e-8 달라, L3 에서 엔진
값과 1e-12 로 견줄 수 없다(쉬는 모양 가까이에서는 C 가 작아 힘의 상대 차이가 훨씬 커진다). 스텐실 값으로 받으면 엔진 갈래가
엔진의 값을 그대로 넘기고, 쓰는 쪽은 정확한 √2 나 요소마다 다른 쉬는 길이를 줄 수 있다. 값은 삼각형마다 double 둘이 는다.

**실행 전에 한 것(판정 밖).** 생성기를 고치지 않고 이 명세에서 모든 변형이 생성되는 것만 봤다(리포 밖 임시 폴더) — 평탄판
연산 399 · 단계판 333 · `EvalImplicit` 689 · GN 406 · T 판 692, 엔진 꼴 넷도 생성됐다. 값은 아무것도 견주지 않았다.

## 명제 · 판정 기준 (실행 전 고정)

인장 · 전단 게이트(`docs/GATE_STRETCH_SHEAR.md`)의 T1 – T5 와 인장 판 게이트(`docs/GATE_TENSION.md`)의 P0 – P2 를 이 명세에
그대로 옮긴다. 문턱 · 표본 수 · 유한차분 방법은 그 두 게이트와 같다.

| 기준 | 무엇 | 통과 |
|---|---|---|
| V0 전제 | 생성기를 고치지 않는다. 고쳤다면 기존 명세 다섯의 변형 해시가 그대로다(`check_flat_unchanged.py` · `check_tension_p0.py` · `check_link_l1_p0.py`) | 모두 통과 |
| V1 Mathematica | `Eval` · `EvalStaged` 대 Mathematica(`model/provot_shear_ref.wl`, uv · b 를 40 자리로 넣고 `D` 로 미분). 처음 200 개 | `|ΔC|/√a` ≤ 1e-12, `∇C` · `∇²C` 상대 ≤ 1e-10 |
| V1b 엔진 꼴과 같은 모형 | Mathematica 가 엔진 꼴 `C_e = a(|w_u ± w_v| − b)` · `k_e = 1 / a` 의 힘 `−k_e Σ C_e ∇C_e` 를 따로 낸다. 생성 `Eval` 의 `−Σ C∇C`(k = 1)와 견준다. 처음 200 개. 척도 = `√a × max|∇C|`(쉬는 표본에서는 힘이 0 이라 힘 자체로 나누지 않는다 — 실행 전에 보탬) | `|Δf|` / 척도 ≤ 1e-10 |
| V2 유한차분 | `Eval` · `EvalStaged` 의 `∇C` · `∇²C`. 중심차분, 간격 = `|x₁ − x₀|` × 1e-6, 분모 = 실제 간격. 1000 개 | 상대 ≤ 1e-6 |
| V3 단계판 대 평탄판 | 1000 개 | `|ΔC|/√a` · `∇C` · `∇²C` ≤ 1e-12 |
| V4 `EvalImplicit` 대 `AssembleImplicit` | `A` · `b` · `f + d`. k = 1, kd = 0.1, h = 0.01. 1000 개 | 상대 ≤ 1e-12 |
| V5 SXMesh 삼각형 루프 | 명세 선언대로 만든 변수(uv 6 칸 + 대각선 2 칸 — 요소 값 변수가 둘인 첫 명세)를 거친 루프 대 직접 모은 값. 격자 200 × 200 노드 | 비트 동일(불일치 0) |
| V6 T 판 정의 | `EvalImplicitTension` 대 기준 조립(`AssembleImplicit(∇²C = 0)` 에 `A += h² k Σ C⁺∇²C`, `b −= h² k Σ C⁺(∇²C v)`). 1000 개 | `|ΔC|/√a` · A · b · f + d ≤ 1e-12 |
| V6 갈래 | 표본의 성분 중 C > 0 과 C ≤ 0 | 각 100 개 이상 |
| V6b 교차 | 두 성분 모두 C > 0 · v = 0 이면 `EvalImplicit` 과, 두 성분 모두 C ≤ 0 이면 `EvalImplicitGN` 과 A · b · f + d 가 같다 | 각 1e-12 이하(해당 표본이 있어야 한다) |
| V7 준정부호 | 모든 표본에서 `A_T` 의 최소 고윳값 / 최대 고윳값 크기 | −1e-12 이상 |

**GO** = V0 – V7 모두 통과. 하나라도 미달이면 NO-GO 이고, 원인을 찾아 명세 · 생성기를 고친 뒤에는 새 절에 까닭을 적고 다시 잰다.

**보고만**: 연산 수(평탄 · 단계 · `EvalImplicit` · GN · T), `∇²C` 를 버린 몫(GN 대 정확판의 탄성 A 차이), 정확판 요소 K 가
부정치인 표본 수, 표본의 `|w_u ± w_v|` 최소 · 최대.

## 표본

- 삼각형 · 속도: 인장 · 전단 게이트의 `results/triangles.txt` 그대로(1000 개, 시드 20261113). 변형이 `w = R F`,
  `F = [[s_u, γ], [0, s_v]]`, s ∈ [0.8, 1.3], γ ∈ [−0.4, 0.4] 라 `|w_u ± w_v| ≥ s_v ≥ 0.8` 이다 — 길이가 0 에 가까운 표본이 없다.
- 대각선 길이: `model/make_provot_diag.py`(새로, 시드 고정)가 `results/provot_diag.txt` 에 줄마다 b₊ b₋ 를 쓴다. 짝수 줄은 √2 √2,
  홀수 줄은 √2 × [0.9, 1.1] 의 무작위 값 둘(b 가 실제 변수로 들어가는지 본다). 처음 4 개(고정 표본)는 √2 √2 다.
- V5 의 격자: 인장 · 전단 게이트 T5 와 같은 격자에 대각선 값 √2 √2.

## 결과 (2026-11-13) — GO

`results/gate_provot_shear.json`. 첫 실행에서 모두 통과했다. 생성기는 고치지 않았다.

| 기준 | 값 | 문턱 |
|---|---|---|
| V0 기존 변형 불변 | `check_flat_unchanged.py` · `check_tension_p0.py` · `check_link_l1_p0.py` 모두 PASS(명세 다섯의 변형 해시 그대로) | 통과 |
| V1 `Eval` `C/√a` · `∇C` · `∇²C` (200 개) | 5.9e-16 · 5.4e-16 · 6.5e-16 | 1e-12 · 1e-10 · 1e-10 |
| V1 `EvalStaged` | 5.9e-16 · 5.4e-16 · 7.0e-16 | 같음 |
| V1b 엔진 꼴 힘 (200 개) | 6.4e-16 | 1e-10 |
| V2 `Eval` `∇C` · `∇²C` (1000 개) | 2.3e-10 · 2.2e-10 | 1e-6 |
| V2 `EvalStaged` | 2.3e-10 · 2.2e-10 | 1e-6 |
| V3 단계판 대 평탄판 `C/√a` · `∇C` · `∇²C` | 0 · 4.4e-16 · 4.9e-16 | 1e-12 |
| V4 `A` · `b` · `f + d` | 9.5e-16 · 8.8e-15 · 1.2e-14 | 1e-12 |
| V5 삼각형 루프 (노드 40,000 · 삼각형 79,202 · 요소 값 변수 둘) | 불일치 0 | 0 |
| V6 T 판 정의 `C/√a` · A · b · f + d (1000 개) | 0 · 8.8e-16 · 9.4e-15 · 1.2e-14 | 1e-12 |
| V6 갈래 | C > 0 1,261 개 · C ≤ 0 739 개 | 각 100 |
| V6b (i) 모두 C > 0 · v = 0 (360 개) | `EvalImplicit` 과 A · b · f + d 모두 0(비트까지 같음) | 1e-12 |
| V6b (ii) 모두 C ≤ 0 (99 개) | `EvalImplicitGN` 과 A 6.4e-16 · b 5.3e-15 · f + d 7.3e-16 | 1e-12 |
| V7 준정부호 (1000 개) | 최소 고윳값 / 최대 크기의 최솟값 −2.9e-16 | −1e-12 |

- **V1b**: 엔진 꼴(`C = a(…)`, 강성 `1/a`)의 힘과 명세(`√a` 곱, k = 1)의 힘이 반올림 수준으로 같다. 곧 강성을 `k = SE` 로 두면
  엔진과 같은 모형이다. 엔진 DLL 의 값과 견주는 것은 L3 몫이다.
- **V5**: 요소 값 변수가 둘(`pUV` 6 칸 · `pDiag` 2 칸)인 첫 명세다. `sx_stencil.h` 는 고치지 않았다 — 선언 순서대로 모은다.

**보고만**

- 연산 수(`results/codegen_ProvotShear.json`): 평탄 399 · 단계 333 · `EvalImplicit` 689 · GN 406 · T 692. 엔진 꼴은 `EvalEngine` 557 ·
  GN 351 · T 561 · 감쇠 없는 GN 285.
- 병진 `|Σ∇C|/max|∇C|` 8.0e-16, 회전 6.3e-16, 쉬는 모양 `|C|/√a` 2.7e-15, `EvalFrozen` 은 `Eval` 과 비트까지 같다.
- 표본의 `|w_u ± w_v|` 는 0.942 – 2.081 이다(0 에 가까운 표본 없음).
- 정확판 A 는 1000 개 중 712 개가 부정치다(`|최소|/최대` 중앙 5.6e-2). T 판은 0 개, GN 은 최솟값 −3.2e-16 이다. 인장과 같은 그림이다 —
  두 대각선 가운데 하나는 거의 언제나 눌린다.
- `|A_T − A_E|/|A_E|` 중앙 0.123 · 최대 1.14, `|A_GN − A_E|/|A_E|` 중앙 0.126 · 최대 1.14.

**이 게이트가 말하지 않는 것.** 엔진 DLL 의 값과 같은지(L3), 이 전단으로 천을 적분했을 때의 안정 · 시간(L4 꼴의 판정), 엔진의
비선형 재질 · 전단 다항식 갈래. 퇴화 가드(넓이 · 길이 0)는 명세 밖이다.

## 하는 일 (확인 뒤)

1. `model/make_provot_diag.py` · `model/provot_shear_ref.wl`(→ `results/mma_ref_provot.txt`, Mathematica 가 없으면 커밋된 것을 쓴다).
2. `tests/test_provot_shear.cpp` — V1 – V7. 인장 · 전단 하네스(`test_stretch_shear.cpp`)와 인장 판 하네스(`test_tension.cpp`)의 꼴을 따른다.
3. `run_gate.bat` 에 생성 한 줄 · 빌드 · 실행 단계를 더한다.
4. 통과하면 L3 판정문(`docs/GATE_LINK_L3.md`)을 쓴다 — 엔진 브랜치에 인장 T 판 · Provot 전단 갈래를 넣고 엔진 손 코드와 견준다.
   감쇠의 `Ċ∇²C` 몫을 어떻게 다룰지를 거기서 먼저 정한다.

## 재현

`run_gate.bat` 22 단계(생성 · 표본 · Mathematica · 빌드는 앞 단계에 섞여 있고, `test_provot_shear.exe` 는 인장 · 전단 게이트 바로
뒤에 돈다). 따로 돌리려면 리포 루트에서:

```
python codegen/sxgen.py codegen/specs/provot_shear.py --out generated --results results
python model/make_provot_diag.py results/provot_diag.txt
wolframscript -file model/provot_shear_ref.wl 200        (없으면 커밋된 results/mma_ref_provot.txt 를 쓴다)
cl … tests/test_provot_shear.cpp … && test_provot_shear.exe
```

22 단계를 넣은 `run_gate.bat` 전체 실행은 아직 안 했다(이 게이트에 필요한 단계만 따로 돌렸다).
