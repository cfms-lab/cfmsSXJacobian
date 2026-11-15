(* baraff_stretch_shear_ref.wl — 독립 기준: Mathematica 로 원문 §4.2 인장 · 전단 조건의 C, grad C, Hessian C 를 40 자리로 구한다.
   생성기(SymPy)와 다른 CAS 로 같은 값을 내는지 보려는 것이다. docs/GATE_STRETCH_SHEAR.md.
   입력  results/triangles.txt  (model/make_triangles.py; 한 줄 = x(9) · uv(6) · v(9))
   출력  results/mma_ref_tri.txt  한 줄 = 인장 C(2) · g(2 × 9) · H(2 × 45) · 전단 C(1) · g(9) · H(45), 유효숫자 25
         g · H 는 성분 순서대로 이어 붙이고, H 는 위 삼각 행 우선이다(생성 코드와 같은 배치).
   사용  wolframscript -file model/baraff_stretch_shear_ref.wl [표본 수, 기본 200]
   (Δx1 Δx2) = (w_u w_v) {{Δu1, Δu2}, {Δv1, Δv2}} 를 그대로 풀어 쓴다(Inverse). uv 는 표본마다 40 자리 수로 넣고 x 는 기호로 둔 채 D 로 미분한다. *)

nUse = If[Length[$ScriptCommandLine] >= 2, ToExpression[$ScriptCommandLine[[2]]], 200];

vars = Array[x, 9];
x0 = vars[[1 ;; 3]]; x1 = vars[[4 ;; 6]]; x2 = vars[[7 ;; 9]];

fmt[v_] := Module[{r = N[v, 30], m, ex},
  If[r == 0, "0",
   {m, ex} = MantissaExponent[r];
   ToString[NumberForm[m, 25, NumberPoint -> ".", ExponentFunction -> (Null &)]] <> "e" <> ToString[ex]]];

upper[h_] := Flatten[Table[h[[p, q]], {p, 1, 9}, {q, p, 9}]];

lines = Select[ReadList["results/triangles.txt", String], StringLength[#] > 0 &];
lines = Take[lines, Min[nUse, Length[lines]]];
out = OpenWrite["results/mma_ref_tri.txt"];
Do[
  nums = ToExpression /@ StringSplit[StringReplace[lines[[n]], {"e" -> "*^"}]];
  xv = SetPrecision[nums[[1 ;; 9]], 40];
  uv = SetPrecision[nums[[10 ;; 15]], 40];
  {u0, v0, u1, v1, u2, v2} = uv;
  M = {{u1 - u0, u2 - u0}, {v1 - v0, v2 - v0}};
  a = Abs[Det[M]]/2;
  W = Transpose[{x1 - x0, x2 - x0}] . Inverse[M];   (* 3 × 2, 열 = w_u, w_v *)
  wu = W[[All, 1]]; wv = W[[All, 2]];
  conds = Sqrt[a] {Sqrt[wu . wu] - 1, Sqrt[wv . wv] - 1, wu . wv};   (* √a: 사용자 결정 2026-11-13, 원문은 a *)
  rule = Thread[vars -> xv];
  row = {};
  (* 인장 두 성분: C 둘, g 둘, H 둘 *)
  cs = conds[[1 ;; 2]];
  row = Join[row, cs /. rule];
  row = Join[row, Flatten[Table[D[c, {vars}] /. rule, {c, cs}]]];
  row = Join[row, Flatten[Table[upper[D[c, {vars, 2}] /. rule], {c, cs}]]];
  (* 전단 *)
  c = conds[[3]];
  row = Join[row, {c /. rule}, D[c, {vars}] /. rule, upper[D[c, {vars, 2}] /. rule]];
  WriteString[out, StringRiffle[fmt /@ row, " "] <> "\n"],
  {n, 1, Length[lines]}];
Close[out];
WriteString[$Output, "wrote " <> ToString[Length[lines]] <> " reference rows to results/mma_ref_tri.txt\n"];
