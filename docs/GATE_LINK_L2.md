# 엔진 연결 L2 — 엔진 안 굽힘 값 일치 (생성 Bridson GN 갈래)

상태: **PASS (2026-11-13)** — 기준은 실행 전에 고정했고 아래 「결과」는 그 뒤에 적었다. 계획은 `docs/PLAN_ENGINE_LINK.md` 의 L2 줄이다.

## 왜

L0 는 엔진을 고치지 않고 규약을 쟀고, L1 은 엔진이 받는 꼴(`EvalEngine*`)을 이 리포 안에서 맞췄다. L2 는 처음으로 **엔진
안에** 생성 코드를 넣는다. 엔진 기본 굽힘(Bridson, Oh 손 코드)의 자리에 `BridsonBendEngine::EvalEngineGN` 을 꽂고, 엔진의
디버그 고리로 꺼낸 힘 · `J` 가 이 리포의 정의와 맞는지, 그리고 **스위치를 끄면 엔진이 바이트 하나 다르지 않게 그대로인지**를 본다.
시간 · 안정은 보지 않는다(L4 · L6).

## 엔진 갈래

- 리포 · 브랜치: `cfmsDrape_dev` 의 브랜치 `sxjacobian-link`(main 378d20c 에서 분기), 워크트리 `..\cfmsDrape_dev_sxjacobian-link`
  (사용자 결정 2026-11-13). 엔진 main · 엔진 원격은 건드리지 않는다.
- 바꾸는 파일(허용 목록, 이 밖이 바뀌면 L2-0 실패):
  `drape/CMakeLists.txt`(소스 한 줄), `drape/src/sde_core/SDECore.h`(선언 · 멤버), `drape/src/sde_core/SDECore.cpp`(초기화 한 줄 ·
  분기 한 줄), `drape/src/sde_core/SDECore_SXJBend.cpp`(새 커널), `drape/src/sde_core/sxgen/BridsonBend_Engine.h`(생성 헤더 사본).
- 스위치: 환경 변수 `CFMS_CPU_SXJ_BEND` 를 `SDECore::Init` 에서 읽는다(엔진의 `CFMS_CPU_INTEGRATOR` · `CFMS_CPU_MEMBRANE_NONLINEAR`
  와 같은 관례 — C API · 기능표 · Python 포장은 바꾸지 않는다). 값: 없음 · 빈 값 · `0` · `off` = **끔(엔진 그대로)**, `gn_rest`,
  `gn_state`. 그 밖의 값은 끔 + stderr 한 줄.
- 굽힘 형식은 그대로 `SDEibfImplicitBridson` 이다. `CheckBendingChange_Macri_SOA`(θ · 힌지 기하)는 그대로 부르고, 그 뒤
  `GetImpBending_Oh_SOA` 대신 `GetImpBending_SXJ_GN` 을 부른다. 그래서 Bridson 형식에 걸린 엔진의 다른 준비 단계는 모두 그대로 돈다.
- 커널(힌지마다, OpenMP, 힌지 소유 슬롯 + 힌지 순서 직렬 환원 — 엔진 Oh 와 같은 결정성 구조):
  1. 퇴화 힌지 판별 — Oh 와 같은 식 · 같은 값(`sbv._nAim` · `_nBim` · `_cb` · `x32_size`). 걸리면 힘 슬롯을 0 으로 두고 환원에서 뺀다.
  2. 좌표 `sbv._x0[m]`(m = 엔진 `_vid` 순서 = 명세 `(h, i, j, k)`, Bridson 게이트에서 확인한 대응).
  3. `EvalEngineGN(x, w, v = 0, {k = BE_U · GlobalBendingFactor, kd = 0})`.
  4. 힘 → `pBendingF`(엔진 순서 nn = `(h, k, j, i)`, 엔진 `AssembleForce` 의 모으기에서 확인).
  5. 노드쌍 블록 6 개 → 힌지 슬롯 `pBendingJ[6i + s]`. 생성 쌍 p(명세 사전순 (h,i)(h,j)(h,k)(i,j)(i,k)(j,k)) → 엔진 슬롯
     s = {2, 1, 0, 5, 4, 3}(엔진 슬롯은 nn 쌍 (0,1)(0,2)(0,3)(1,2)(1,3)(2,3), `SDESparseMatrix.cpp` 의 슬롯 만들기에서 확인).
     엔진 슬롯은 (작은 번호, 큰 번호)라 노드 번호가 뒤집힌 쌍은 전치해 넣는다(요소 GN 행렬은 대칭, L1-3 과 같은 규칙).
  6. 직렬 환원은 **슬롯 6 개 모두**를 더한다(Oh 는 1 – 4 만 쓴다). 전치 처리 플래그(`bUseTransposeAwareUpperJacobian`)를 켠다
     — 생성 블록은 비대칭이다(L0 · L1).
- `w`: `gn_state` = Oh 와 같은 식 `3 Length0² / (½/nAim + ½/nBim)`(지금 넓이, 엔진 힘과 같아야 한다). `gn_rest` =
  `3 Length0² / (½|nA⁰| + ½|nB⁰|)`, 넓이는 쉬는 2D 패턴(`ppNodes2D`) — 이 리포 Bridson · GN · 적분 게이트의 `w̄` 식.
  쉬는 넓이가 퇴화해 `w̄` 가 유한한 `[0, 1e8)` 밖이면 그 힌지는 퇴화 힌지처럼 뺀다(`gn_rest` 만, 엔진 판별에 하나 더한 것).
- 하지 않는 것:
  - **굽힘 감쇠**(BV, 엔진 기본 꺼짐). 감쇠가 켜졌거나 upper · serial 조립이 아니면 이 갈래는 Oh 커널로 물러서고 stderr 에
    한 번 알린다. 엔진 감쇠 `−BV·w·θ̇·∇θ` 와 명세 감쇠 `−kd Ċ ∇C`(= `−kd·w̄cos²(θ/4)/2·θ̇∇θ`)가 모형부터 달라 L4 전에 따로 정한다.
  - θ0 ≠ 0. 명세는 θ0 = 0 이다. C API 는 `Theta0` 를 설정하지 않는다(`SDEMaterialProperty` 기본 0, 2026-11-13 읽음).
- 생성 헤더 사본은 이 리포 `generated/BridsonBend_Engine.h` 와 바이트까지 같아야 한다(L2-0 (d)). 엔진 소스는 `/utf-8` 로 빌드된다.

## w 를 어느 쪽으로 — 여기서 정하지 않고 둘 다 잰다

`gn_state` 는 힘이 엔진 Oh 와 같고 `J` 만 다르다 — 엔진 손 Jacobian 을 GN 으로 바꾼 효과만 떼어 본다. 다만 그 힘은 상수 w̄
에너지의 기울기가 아니어서 GN 은 `f∇wᵀ/w` 몫을 버린다(Bridson 게이트 B7). `gn_rest` 는 이 리포 명세 그대로라 `J` 가 참 에너지의
GN 이지만 힘이 엔진과 `w̄ / w` 만큼 다르다. 어느 쪽을 쓸지는 L4 의 적분 규칙으로 정한다(계획 「위험 · 굽힘 가중치」).

## 장면 — L0 와 같다

4 × 4 칸, 간격 2 cm(노드 25 · 삼각형 32 · 힌지 40), uv = 쉬는 (x, y), 같은 탐침 상태(`tests/link_l0.py` 의 `scene()`, 시드
20261113) · 무작위 벡터 5 개 · 질량 0.1 · 중력 0 · dt 0.005. 재질은 L0 의 넷(stretch · shear · bend · all) 그대로다.

매단 천(L2-1 의 시간 적분 확인만): 같은 격자를 쉬는 평면 모양에서 시작, 꼭짓점 0 · 4(한 변의 두 끝)를 고정, 중력 980,
dt 0.005, 재질 all, 40 단계.

DLL 둘:
- **main DLL** — L0 와 같은 빌드(`%TEMP%\cfmsDrape_build_378d20c`, 엔진 main 378d20c).
- **브랜치 DLL** — 워크트리를 `%TEMP%\cfmsDrape_build_sxjlink` 에 빌드(`--target cfmsDrape` 만, 리포 밖). 모드마다 다른
  프로세스에서 환경 변수를 주고 잰다.

## 판정 기준 (실행 전 고정)

상대 오차는 L0 와 같은 식 `max|Δ| / max|기준|` 이다.

| 기준 | 무엇 | 통과 |
|---|---|---|
| L2-0 (a) | main DLL `build_info` | `git_commit` 378d20c · `dirty` false |
| L2-0 (b) | 브랜치 DLL `build_info` · 워크트리 | `git_commit` = 워크트리 HEAD, HEAD 가 378d20c 의 자손, 378d20c 대비 바뀐 파일(커밋 + 미커밋)이 허용 목록 안. `dirty` 는 보고(미커밋이면 true — 그때는 `git diff 378d20c` 의 SHA-256 을 결과에 남긴다) |
| L2-0 (c) | 엔진 main 워크트리 | `git status --porcelain` 비어 있음 |
| L2-0 (d) | 생성 헤더 사본 | SHA-256 = 이 리포 `generated/BridsonBend_Engine.h` |
| L2-1 | 기본 경로 그대로 (브랜치 DLL 스위치 끔 대 main DLL) | 네 재질 모두 f · 단위 벡터 75 개의 곱 · 무작위 벡터 5 개의 곱이 **비트까지 같다**, 매단 천 40 단계의 위치가 비트까지 같다 |
| L2-2 | 힘 (bend 재질) | `gn_state` f 대 main(Oh) f ≤ 1e-12, `gn_state` f 대 이 리포 기준(w = 그 상태) ≤ 1e-12, `gn_rest` f 대 이 리포 기준(w̄ = 쉬는 값) ≤ 1e-12 |
| L2-3 | J (bend 재질) | 두 모드 모두 엔진 `J`(단위 벡터 곱으로 되살림) 대 이 리포 GN 전역 `−k Σ ∇C∇Cᵀ`(그 모드의 w) ≤ 1e-12, 곱 일관성(되살린 J 로 곱한 값 대 엔진 곱, 무작위 5 개) ≤ 1e-12 |
| L2-4 | 합성 (all 재질) | 두 모드 모두 `J_all − J_bend` 가 스위치 끔의 `J_all − J_bend` 와 ≤ 1e-12(척도 = 끔의 max\|J_all\|), f 도 같은 식 — 인장 · 전단과 같이 쓰는 노드쌍 슬롯에서 굽힘 몫이 덮이거나 두 번 더해지지 않는다 |
| L2-5 | 결정성 (브랜치 DLL, 두 모드 · bend · all) | 같은 프로세스 두 번 · 단일 스레드 프로세스 대 기본 스레드에서 f · J 가 비트까지 같다 |

1e-12 의 근거: L1 과 같다. 엔진 환원(힌지 순서)과 기준 합(numpy, 힌지 순서)의 반올림 차이는 한 칸에 기여 수 개 × 1e-16 이다.

보고(판정 밖): `gn_rest` f 대 Oh f(w 차이 몫), 두 모드 J 대 Oh J, 두 모드 J 대 엔진 힘의 유한차분, `w̄ / w` 분포,
매단 천 40 단계에서 두 모드의 유한성 · 끝 위치 차이(대 끔).

## 결과

2026-11-13, `python tests/link_l2.py <main DLL> <브랜치 DLL> <link_l0_ref.exe>`(`run_gate.bat` 19 단계와 같은 호출). 숫자는
`results/gate_link_l2.json`. **실패 0 — PASS.**

| 기준 | 결과 |
|---|---|
| L2-0 (a) | main DLL `git_commit` 378d20c · `dirty` false |
| L2-0 (b) | 브랜치 DLL `git_commit` 54242be = 워크트리 HEAD(브랜치 첫 커밋, 378d20c 의 자식), 바뀐 파일 = 허용 목록 5 개 그대로, `dirty` false. 커밋 전에도 한 번 쟀다 — 그때는 `git_commit` 378d20c · `dirty` true 였고, 나머지 판정 · 보고 값은 지금과 모두 같았다 |
| L2-0 (c) | 엔진 main 워크트리 깨끗 |
| L2-0 (d) | 헤더 사본 SHA-256 = 생성 헤더(`df662d61…`) |
| L2-1 | 네 재질 모두 f · J · 무작위 곱이 비트까지 같다, 매단 천 40 단계 위치도 비트까지 같다 |
| L2-2 | `gn_state` f 대 엔진 Oh 2.8e-16 · 대 기준(w = 그 상태) 4.3e-16, `gn_rest` f 대 기준(w̄ = 쉬는) 3.6e-16 |
| L2-3 | J 대 GN 기준 `gn_rest` 3.5e-16 · `gn_state` 3.3e-16, 곱 일관성 3.5e-16 · 4.5e-16 |
| L2-4 | `J_all − J_bend` 대 끔 `gn_rest` 3.5e-16 · `gn_state` 5.2e-16, f 2.2e-16 · 6.9e-18 |
| L2-5 | 두 모드 · bend · all 모두 같은 프로세스 두 번 · 스레드 1 대 32 가 비트까지 같다 |

보고(판정 밖):

| 무엇 | 값 | 읽기 |
|---|---|---|
| `gn_rest` f 대 엔진 Oh f | 4.3e-2 | w̄ 와 지금 넓이 w 의 차이 몫. 이 탐침 상태에서 `w̄ / w` = 0.867 – 1.195(중앙 1.016) |
| J 대 엔진 Oh J | `gn_rest` 0.52 · `gn_state` 0.50 | Bridson 게이트 B7 대로 엔진 Oh J 는 GN 이 아니다 |
| J 대 자기 힘의 유한차분 | `gn_rest` 0.089 · `gn_state` 0.12 · 엔진 Oh 1.06 | GN 이 버린 몫(`gn_rest` = `C∇²C`, `gn_state` = 그것 + `f∇wᵀ/w`). 엔진 Oh J 는 자기 힘의 참 Jacobian 과 크기만큼 다르다 |
| 매단 천 40 단계 | 두 모드 모두 유한, 끝 위치 차이(대 끔) 3.6e-15 · 4.4e-15 cm | **이 장면은 굽힘을 거의 쓰지 않는다**(재질 all 의 BE_U 1 이 인장 · 전단 1e5 와 중력에 묻힌다) — 단계 경로가 새 갈래를 탄다는 것만 보이고, 적분 품질은 L4 몫이다. (L4 에서 알게 됨: 엔진은 중력을 **−Y 로만** 당긴다 — 이 장면은 천이 z = 0 평면이고 중력을 `(0, 0, −980)` 으로 줬으므로 실제로는 판 안쪽으로 당겼다. 판정과는 무관하다) |

stderr 에 물러섬 알림은 없었다(두 모드 모두 생성 커널이 돌았다). 브랜치 첫 빌드 0.5 분, 새 파일에서 나온 컴파일 경고 없음
(엔진의 기존 경고 C4828 등만).

**L4 에 넘길 것.** w 를 쉬는 값 · 지금 넓이 중 무엇으로 둘지, 굽힘 감쇠를 어떻게 넣을지(엔진 `−BV·w·θ̇·∇θ` 대 명세
`−kd Ċ∇C`), PCG 20 회 고정 아래에서 GN 이 엔진 Oh 보다 나은지. 엔진 브랜치는 54242be 로 커밋했다(로컬만 — 엔진 원격 push · main 합치기는 하지 않았다).

## 브랜치의 기본 굽힘을 `gn_state` 로 (2026-11-13, 사용자 결정)

L4R(`docs/GATE_LINK_L4R.md`)의 권고에 따라 사용자가 엔진 굽힘을 `gn_state` 로 바꾸기로 정했다. 바꾼 곳은 엔진 브랜치
`sxjacobian-link` 뿐이다. 엔진 main 과 엔진 원격은 그대로다(합치기 · push 는 요청이 있을 때).

- `CFMS_CPU_SXJ_BEND` 가 없거나 빈 값이면 `gn_state`, `off`(또는 `0`)면 엔진 Oh, `gn_rest` · `gn_state` 는 그대로다.
  모르는 값은 알리고 기본(`gn_state`)으로 둔다.
- 기본으로 켜진 때는 굽힘 감쇠 등으로 Oh 로 물러서도 알리지 않는다(환경 변수로 직접 골랐을 때만 알린다).
- 바꾼 파일은 `SDECore_SXJBend.cpp`(모드 읽기) 와 `SDECore.cpp` · `SDECore.h` 의 주석 세 줄이다. 계산 경로는 손대지 않았다.
- 이 리포의 하네스(`link_l2.py` · `link_l4.py` · `link_l4r.py` · `diag_link_l4_*.py`)는 이제 끔을 `off` 로 명시한다. 「환경 변수
  없음」이 더는 끔이 아니기 때문이다. `off` 는 옛 브랜치 DLL(e9d8880)도 받는 값이다.

기준(바꾸기 전에 정함): L2-1 은 「`off` 로 끈 브랜치 = main, 비트까지」로 읽는다. 새 **L2-6** = 「환경 변수 없는 브랜치 =
`gn_state`, 비트까지」(네 재질의 f · J · 무작위 곱과 매단 천 40 단계).

결과(`results/gate_link_l2.json`): **L2-0 – L2-6 모두 PASS.** 브랜치를 커밋하기 전 작업 트리(DLL `dirty` true)와 커밋한 뒤
(d784a77, `dirty` false — 로컬만) 두 번 돌렸고 값이 같다. 커밋된 결과는 커밋 뒤의 것이다. L2-1 과 L2-6 은 다섯 항목이 모두 비트까지
같다. L2-2 – L2-5 의 값은 앞 실행과 같다(힘 · J 4.5e-16 이하). 브랜치 DLL 의 커밋이 e9d8880 에서 d784a77 로 바뀌었으므로, L4R 의
전제 (a)(L4 2 차와 같은 커밋)는 `run_gate.bat` 전체 실행에서 20 단계가 먼저 돈 뒤에만 선다. L4 · L4R 은 이 DLL 로 아직 다시 돌리지 않았다.

## 브랜치의 인장 · 전단 기본이 `on` 이 되고 엔진 main 에 합친 뒤 (2026-11-13)

L3(`docs/GATE_LINK_L3.md`) 뒤 사용자 결정으로 브랜치의 인장 · 전단 기본이 생성 커널이 됐고(`CFMS_CPU_SXJ_MEMBRANE` 없음 = `on`),
브랜치를 엔진 main 에 합쳐 push 했다(엔진 92d8838). L2 는 굽힘만 보는 판정이라 `link_l2.py` 가 모든 실행에
`CFMS_CPU_SXJ_MEMBRANE=off` 를 준다(main DLL 은 이 변수를 모른다). 허용 목록에 `SDECore_SXJMembrane.cpp` · 헤더 사본 둘 ·
엔진 `CLAUDE.md`(기본값 전환 기록)를 더했다. 「main DLL」은 378d20c 에 고정한 워크트리 `..\cfmsDrape_dev_base378d20c` 에서 빌드한다.
이 상태로 **L2-0 – L2-6 모두 PASS**(`results/gate_link_l2.json`, 브랜치 DLL 92d8838 · `dirty` false).
