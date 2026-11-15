(* sul_poisson_ref.wl — 독립 기준: Mathematica 로 Poisson 조건(Sul 외 2010, 식 (6))의 C, grad C, Hessian C 를 40 자리로 구한다.
   생성기(SymPy)와 다른 CAS 로 같은 값을 내는지 보려는 것이다. docs/GATE_POISSON.md.
   입력  results/triangles.txt     (model/make_triangles.py; 한 줄 = x(9) · uv(6) · v(9))
         results/poisson_beta.txt  (model/make_poisson_beta.py; 한 줄 = β)
   출력  results/mma_ref_poisson.txt  한 줄 = C(1) · g(9) · H(45) · 원문 닫힌 꼴 g(9) · 원문 닫힌 꼴 H(45), 유효숫자 25
         H 는 위 삼각 행 우선이다(생성 코드와 같은 배치).
         C = Sqrt[a] (Log|w_u| + β Log|w_v|) 를 D 로 미분한 것이 앞의 g · H 이고,
         뒤의 둘은 원문 식 (7) · (8) 을 D 없이 그대로 쓴 것에 Sqrt[a] 를 곱한 것이다(Q1b):
           grad_m  = c^u_m w_u/|w_u|^2 + β c^v_m w_v/|w_v|^2
           hess_mn = c^u_m c^u_n (I - 2 ŵ_u ŵ_u^T)/|w_u|^2 + β c^v_m c^v_n (I - 2 ŵ_v ŵ_v^T)/|w_v|^2
         c^u_m, c^v_m = ∂w_u/∂x_m, ∂w_v/∂x_m 의 계수(쉬는 uv 의 역행렬에서 나온다).
   사용  wolframscript -file model/sul_poisson_ref.wl [표본 수, 기본 200]
   (Δx1 Δx2) = (w_u w_v) {{Δu1, Δu2}, {Δv1, Δv2}} 를 그대로 풀어 쓴다(Inverse). uv · β 는 표본마다 40 자리 수로 넣고 x 는 기호로 둔 채 D 로 미분한다. *)

nUse = If[Length[$ScriptCommandLine] >= 2, ToExpression[$ScriptCommandLine[[2]]], 200];

vars = Array[x, 9];
x0 = vars[[1 ;; 3]]; x1 = vars[[4 ;; 6]]; x2 = vars[[7 ;; 9]];

fmt[v_] := Module[{r = N[v, 30], m, ex},
  If[r == 0, "0",
   {m, ex} = MantissaExponent[r];
   ToString[NumberForm[m, 25, NumberPoint -> ".", ExponentFunction -> (Null &)]] <> "e" <> ToString[ex]]];

upper[h_] := Flatten[Table[h[[p, q]], {p, 1, 9}, {q, p, 9}]];

lines = Select[ReadList["results/triangles.txt", String], StringLength[#] > 0 &];
betas = Select[ReadList["results/poisson_beta.txt", String], StringLength[#] > 0 &];
nRow = Min[nUse, Length[lines], Length[betas]];
out = OpenWrite["results/mma_ref_poisson.txt"];
Do[
  nums = ToExpression /@ StringSplit[StringReplace[lines[[n]], {"e" -> "*^"}]];
  beta = SetPrecision[ToExpression[StringReplace[StringTrim[betas[[n]]], {"e" -> "*^"}]], 40];
  xv = SetPrecision[nums[[1 ;; 9]], 40];
  uv = SetPrecision[nums[[10 ;; 15]], 40];
  {u0, v0, u1, v1, u2, v2} = uv;
  M = {{u1 - u0, u2 - u0}, {v1 - v0, v2 - v0}};
  Mi = Inverse[M];
  a = Abs[Det[M]]/2;
  W = Transpose[{x1 - x0, x2 - x0}] . Mi;   (* 3 × 2, 열 = w_u, w_v *)
  wu = W[[All, 1]]; wv = W[[All, 2]];
  cs = Sqrt[a] (Log[Sqrt[wu . wu]] + beta Log[Sqrt[wv . wv]]);   (* 명세: √a 를 곱한다 *)
  rule = Thread[vars -> xv];
  row = {cs /. rule};
  row = Join[row, D[cs, {vars}] /. rule];
  row = Join[row, upper[D[cs, {vars, 2}] /. rule]];
  (* 원문 식 (7) · (8) 의 닫힌 꼴 — D 를 쓰지 않는다 *)
  cu = {-(Mi[[1, 1]] + Mi[[2, 1]]), Mi[[1, 1]], Mi[[2, 1]]};
  cv = {-(Mi[[1, 2]] + Mi[[2, 2]]), Mi[[1, 2]], Mi[[2, 2]]};
  wuN = wu /. rule; wvN = wv /. rule;
  lu2 = wuN . wuN; lv2 = wvN . wvN;
  gP = Sqrt[a] Flatten[Table[cu[[m]] wuN/lu2 + beta cv[[m]] wvN/lv2, {m, 1, 3}]];
  Bu = (IdentityMatrix[3] - 2 Outer[Times, wuN, wuN]/lu2)/lu2;
  Bv = (IdentityMatrix[3] - 2 Outer[Times, wvN, wvN]/lv2)/lv2;
  HP = Sqrt[a] (KroneckerProduct[Outer[Times, cu, cu], Bu] + beta KroneckerProduct[Outer[Times, cv, cv], Bv]);
  row = Join[row, gP, upper[HP]];
  WriteString[out, StringRiffle[fmt /@ row, " "] <> "\n"],
  {n, 1, nRow}];
Close[out];
WriteString[$Output, "wrote " <> ToString[nRow] <> " reference rows to results/mma_ref_poisson.txt\n"];
