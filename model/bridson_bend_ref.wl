(* bridson_bend_ref.wl — 독립 기준: Mathematica 로 Bridson 굽힘 조건 C = 2 Sqrt[2 w] Sin[θ/4] 의 C, grad C, Hessian C 를 40 자리로 구한다.
   생성기(SymPy)와 다른 CAS 로 같은 값을 내는지 보려는 것이다. docs/GATE_BRIDSON.md.
   입력  results/hinges.txt  (model/make_hinges.py; 한 줄 = x(12) · θ0 · v(12), θ0 는 쓰지 않는다)
   출력  results/mma_ref_bridson.txt  한 줄 = C, g(12), H 위 삼각 행 우선(78), 유효숫자 25
   사용  wolframscript -file model/bridson_bend_ref.wl [표본 수, 기본 200]
   θ 는 착수 게이트 기준(baraff_bend_ref.wl)과 같이 단위 법선 · 단위 변으로 쓴다(명세는 압축형 — 값이 같아야 한다).
   θ · ∇θ · ∇²θ 는 D 로 기호적으로 구해 OptimizeExpression 을 거친 뒤 대입한다(baraff_bend_ref.wl 과 같은 이유).
   바깥 함수 c(t) = 2 Sqrt[2 w] Sin[t/4] 의 c' · c'' 도 D 로 구하고, 2 계 연쇄 법칙
   ∇C = c' ∇θ,  ∇²C = c'' ∇θ ∇θᵀ + c' ∇²θ  로 잇는다.
   w = 3 e·e / (|nA|/2 + |nB|/2) 를 표본마다 40 자리로 정한다(시험도 같은 식을 double 로 쓴다 — 차이는 1e-16 자릿수). *)

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
wexpr = 3 e . e/(Sqrt[nA . nA]/2 + Sqrt[nB . nB]/2);

cf[t_, w_] := 2 Sqrt[2 w] Sin[t/4];
cp[t_, w_] = D[cf[t, w], t];
cpp[t_, w_] = D[cf[t, w], {t, 2}];

fmt[v_] := Module[{r = N[v, 30], m, ex},
  If[r == 0, "0",
   {m, ex} = MantissaExponent[r];
   ToString[NumberForm[m, 25, NumberPoint -> ".", ExponentFunction -> (Null &)]] <> "e" <> ToString[ex]]];

lines = Select[ReadList["results/hinges.txt", String], StringLength[#] > 0 &];
lines = Take[lines, Min[nUse, Length[lines]]];
out = OpenWrite["results/mma_ref_bridson.txt"];
Do[
  nums = ToExpression /@ StringSplit[StringReplace[lines[[n]], {"e" -> "*^"}]];
  vals = SetPrecision[nums[[1 ;; 12]], 40];
  rule = Thread[vars -> vals];
  row = First[opt /. rule];
  t0 = row[[1]]; gt = row[[2 ;; 13]]; Ht = row[[14 ;; 91]];
  w0 = wexpr /. rule;
  ggT = Flatten[Table[gt[[p]] gt[[q]], {p, 1, 12}, {q, p, 12}]];
  rowC = Join[{cf[t0, w0]}, cp[t0, w0] gt, cpp[t0, w0] ggT + cp[t0, w0] Ht];
  WriteString[out, StringRiffle[fmt /@ rowC, " "] <> "\n"],
  {n, 1, Length[lines]}];
Close[out];
WriteString[$Output, "wrote " <> ToString[Length[lines]] <> " reference rows to results/mma_ref_bridson.txt\n"];
