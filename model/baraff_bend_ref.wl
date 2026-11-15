(* baraff_bend_ref.wl — 독립 기준: Mathematica 로 같은 굽힘 조건의 C, grad C, Hessian C 를 40 자리로 구한다.
   생성기(SymPy)와 다른 CAS 로 같은 값을 내는지 보려는 것이다.
   입력  results/hinges.txt  (model/make_hinges.py)
   출력  results/mma_ref.txt  한 줄 = C, g(12), H 위 삼각 행 우선(78), 유효숫자 25
   사용  wolframscript -file model/baraff_bend_ref.wl [표본 수, 기본 200]
   미분은 Mathematica 의 D 로 기호적으로 한다. 결과 식은 잎 9,400만 개라 그대로 대입하면 한 줄에 16 초가 걸려서,
   Mathematica 자신의 공통 부분식 제거(Experimental`OptimizeExpression)를 거친 뒤 대입한다(한 줄 0.15 초).
   두 방식의 처음 27 줄이 글자 그대로 같음을 확인했다(2026-11-12). *)

nUse = If[Length[$ScriptCommandLine] >= 2, ToExpression[$ScriptCommandLine[[2]]], 200];

vars = Array[x, 12];
xh = vars[[1 ;; 3]]; xi = vars[[4 ;; 6]]; xj = vars[[7 ;; 9]]; xk = vars[[10 ;; 12]];
nA = Cross[xi - xh, xj - xh];
nB = Cross[xj - xk, xi - xk];
e = xj - xi;
a = nA/Sqrt[nA . nA]; b = nB/Sqrt[nB . nB]; eh = e/Sqrt[e . e];
s = Cross[a, b] . eh;
c = a . b;
th = ArcTan[c, s];  (* Mathematica ArcTan[x, y] = atan2(y, x) *)
g = D[th, {vars}];
H = D[th, {vars, 2}];
Hup = Flatten[Table[H[[p, q]], {p, 1, 12}, {q, p, 12}]];
opt = Experimental`OptimizeExpression[Join[{th}, g, Hup]];

fmt[v_] := Module[{r = N[v, 30], m, ex},
  If[r == 0, "0",
   {m, ex} = MantissaExponent[r];
   ToString[NumberForm[m, 25, NumberPoint -> ".", ExponentFunction -> (Null &)]] <> "e" <> ToString[ex]]];

lines = Select[ReadList["results/hinges.txt", String], StringLength[#] > 0 &];
lines = Take[lines, Min[nUse, Length[lines]]];
out = OpenWrite["results/mma_ref.txt"];
Do[
  nums = ToExpression /@ StringSplit[StringReplace[lines[[n]], {"e" -> "*^"}]];
  vals = SetPrecision[nums[[1 ;; 12]], 40];
  th0 = SetPrecision[nums[[13]], 40];
  rule = Thread[vars -> vals];
  row = First[opt /. rule];
  row[[1]] = row[[1]] - th0;
  WriteString[out, StringRiffle[fmt /@ row, " "] <> "\n"],
  {n, 1, Length[lines]}];
Close[out];
WriteString[$Output, "wrote " <> ToString[Length[lines]] <> " reference rows to results/mma_ref.txt\n"];
