(* sul_poisson_uniaxial_ref.wl — 독립 기준: 단축 인장에서 균일 변형의 가로 늘음률 s_v 를 40 자리로 구한다.
   docs/GATE_POISSON_INTEGRATION.md. 생성 코드와 무관하게, 넓이당 에너지를 손으로 쓴 식이다:
     e(s_v) = 1/2 k_T ((s_u - 1)^2 + (s_v - 1)^2) + 1/2 k_P (Log[s_u] + beta Log[s_v])^2,   s_u = 1 + eps,  r = k_P / k_T
     de/ds_v / k_T = (s_v - 1) + r beta (Log[s_u] + beta Log[s_v]) / s_v = 0
   인장은 C = Sqrt[a] (|w| - 1), Poisson 은 C = Sqrt[a] (Log|w_u| + beta Log|w_v|), 전단은 균일 변형에서 0 이다.
   출력  results/mma_ref_poisson_uniaxial.txt  한 줄 = r beta eps s_v e2 nroots nulin, 유효숫자 25
         e2 = d^2e/ds_v^2 / k_T (근에서 — 양이면 극소), nroots = s_v 가 0.05 부터 5 까지 0.001 간격일 때 de/ds_v 의 부호가 바뀌는 횟수
         (격자 위에서 정확히 0 인 점도 센다), nulin = 선형화한 Poisson 비 r beta / (1 + r beta^2).
   줄의 순서 = r (0.3, 1) 바깥, beta (-0.6, 0, 0.6, 1), eps (0.02, 0.2, 1) 안쪽 — tests/test_poisson_integration.cpp 와 같다.
   사용  wolframscript -file model/sul_poisson_uniaxial_ref.wl *)

rs = {3/10, 1};
betas = {-6/10, 0, 6/10, 1};
epss = {2/100, 2/10, 1};

fmt[v_] := Module[{r = N[v, 30], m, ex},
  If[r == 0, "0",
   {m, ex} = MantissaExponent[r];
   ToString[NumberForm[m, 25, NumberPoint -> ".", ExponentFunction -> (Null &)]] <> "e" <> ToString[ex]]];

f[s_, r_, b_, e_] := (s - 1) + r b (Log[1 + e] + b Log[s])/s;

out = OpenWrite["results/mma_ref_poisson_uniaxial.txt"];
nRow = 0;
Do[
  sol = If[b == 0, 1, s /. FindRoot[f[s, r, b, e] == 0, {s, 1}, WorkingPrecision -> 40]];
  e2 = D[f[s, r, b, e], s] /. s -> sol;
  grid = Table[N[f[s, r, b, e], 30], {s, 5/100, 5, 1/1000}];
  sg = Sign[grid];
  nroots = Count[Partition[sg, 2, 1], {p_, q_} /; p q < 0] + Count[sg, 0];
  nulin = r b/(1 + r b^2);
  WriteString[out, StringRiffle[fmt /@ {r, b, e, sol, e2, nroots, nulin}, " "] <> "\n"];
  nRow++,
  {r, rs}, {b, betas}, {e, epss}];
Close[out];
WriteString[$Output, "wrote " <> ToString[nRow] <> " reference rows to results/mma_ref_poisson_uniaxial.txt\n"];
