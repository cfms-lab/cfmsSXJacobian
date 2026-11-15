"""cfmsSXJacobian 생성기 — 조건 함수 C(x) 명세에서 C · ∇C · ∇²C 를 기호 미분해 C++ 로 낸다.

사용자는 codegen/specs/<이름>.py 에 **조건 함수 하나**와 그것이 읽는 변수 이름만 적는다.
힘 · Jacobian · 감쇠(Baraff–Witkin 1998 식 7 · 8 · 11 · 12)는 src/sx_condition_force.h 가
C · ∇C · ∇²C 로부터 조립한다. 이 파일은 그 세 가지를 기호 미분하고 공통 부분식을 묶어
generated/<이름>.h 로 쓴다.

변형을 셋 낸다.
  Eval       — 정확한 미분. 좌표에 대해 통째로 펼친다(평탄판).
  EvalFrozen — Baraff–Witkin §4.3 의 가정("n1, n2, e 의 길이를 상수로 본다")을 그대로 적용한 미분(평탄판).
               옛 손 유도가 기대던 근사를 같은 명세에서 기계적으로 만들어 정확한 미분과 비교하려고 둔다.
  EvalStaged — 정확한 미분, 단계판. 명세가 ctx.stage(...) 로 표시한 중간 변수 u 에 대해 C(u) 를 미분하고
               연쇄 법칙을 코드에서 잇는다(명세에 stage 가 있을 때만):
                 ∇ₓC = Jᵀ ∇ᵤC,  ∇²ₓC = Jᵀ ∇²ᵤC J + Σₖ ∂C/∂uₖ ∇²ₓuₖ,  J = ∂u/∂x.
               docs/GATE_SPEED.md.
  EvalImplicit — 단계판에 원문 식 (6) 의 조립을 합친 것. C · ∇C · ∇²C 를 배열로 쓰지 않고
               요소 기여 A = −h ∂f/∂v − h² ∂f/∂x(위 삼각) · b · 전체 힘 f + d 를 바로 낸다. docs/GATE_ASSEMBLY.md.
  EvalImplicitGN — EvalImplicit 의 Gauss–Newton 판. ∇²C 를 0 으로 둔다(힘은 같고 A 는 α Σ ∇C∇Cᵀ).
               docs/GATE_GAUSS_NEWTON.md.
  EvalImplicitTension — 명세가 IMPLICIT_TENSION = True 일 때만. 늘어난 성분(C > 0)의 기하 강성 k C ∇²C 만 남기고
               ∂d/∂x 는 버린다(A = α Σ ∇C∇Cᵀ + h² k Σ max(C, 0) ∇²C). ∇²C 가 준정부호인 조건에서 A 가 언제나 준정부호.
               docs/GATE_TENSION.md.
  EvalEngine · EvalEngineGN · EvalEngineTension — 엔진(cfmsDrape CPU)용. 따로 된 헤더 generated/<이름>_Engine.h 에 쓴다.
               f + d 와 노드쌍 블록 J = ∂f/∂x · V = ∂f/∂v 만 내고(대각은 엔진의 행 합 규약), 강성 판은 위 셋과 같다.
               docs/GATE_LINK_L1.md.
  EvalEngineGNUndamped — EvalEngineGN 의 감쇠 없는 판. v · kd · V 를 받지 않고 C · f · J 만 낸다(엔진 굽힘 갈래가 감쇠
               없이 돌 때 쓴다). EvalEngineGN(v = 0, kd = 0) 과 C · f · J 가 같다. docs/GATE_LINK_L4.md 「2 차 (다)」.

사용:  python codegen/sxgen.py codegen/specs/baraff_bend.py [--out generated] [--results results]
"""
import argparse
import hashlib
import importlib.util
import json
import os
import sys
import time

import sympy as sp
from sympy.printing.c import C99CodePrinter


# ---------------------------------------------------------------------------
# 명세가 쓰는 도구
# ---------------------------------------------------------------------------
class Ctx:
    """명세의 condition() 에 넘기는 도구 상자. mode = 'exact' | 'frozen', staged = 단계판 여부."""

    def __init__(self, mode, staged=False):
        self.mode = mode
        self.staged = staged
        self.frozen = []          # (길이 기호, 원래 길이 식)
        self.stages = []          # (중간 변수 기호, 좌표의 식)

    def norm(self, v):
        return sp.sqrt((v.T * v)[0])

    def unit(self, v):
        """v / |v|. frozen 모드에서는 |v| 를 미분하지 않는 기호로 둔다(미분 뒤 원래 값으로 바꾼다)."""
        if self.mode == 'exact':
            return v / self.norm(v)
        L = sp.Symbol('_L%d' % len(self.frozen), positive=True)
        self.frozen.append((L, self.norm(v)))
        return v / L

    def stage(self, v):
        """중간 변수 표시. 평탄판에서는 v 를 그대로 돌려준다.
        단계판에서는 성분마다 기호 u_k 를 두고 그 기호들을 돌려준다. v 는 좌표만의 식이어야 한다."""
        if not self.staged:
            return v
        isMat = isinstance(v, sp.MatrixBase)
        M = v if isMat else sp.Matrix([v])
        used = set(u for u, _ in self.stages)
        syms = []
        for e in M:
            if e.free_symbols & used:
                raise ValueError('ctx.stage: 중간 변수는 좌표만의 식이어야 한다(다른 stage 의 기호를 쓰면 안 된다)')
            u = sp.Symbol('u%d' % len(self.stages), real=True)
            self.stages.append((u, e))
            syms.append(u)
        out = sp.Matrix(M.shape[0], M.shape[1], syms)
        return out if isMat else out[0]


# ---------------------------------------------------------------------------
# C++ 출력
# ---------------------------------------------------------------------------
class Printer(C99CodePrinter):
    """작은 정수 거듭제곱은 곱셈으로, 반정수 거듭제곱은 sqrt 로 낸다(pow 호출을 피한다)."""

    def _base(self, b):
        return '(' + self._print(b) + ')'

    def _print_Pow(self, expr):
        b, e = expr.base, expr.exp
        if e.is_Integer:
            n = int(e)
            if n == 1:
                return self._print(b)
            if 2 <= abs(n) <= 4:
                prod = '*'.join([self._base(b)] * abs(n))
                return '(' + prod + ')' if n > 0 else '(1.0/(' + prod + '))'
            if n == -1:
                return '(1.0/' + self._base(b) + ')'
        if e == sp.Rational(1, 2):
            return 'std::sqrt(' + self._print(b) + ')'
        if e == -sp.Rational(1, 2):
            return '(1.0/std::sqrt(' + self._print(b) + '))'
        if e == sp.Rational(3, 2):
            return '(' + self._base(b) + '*std::sqrt(' + self._print(b) + '))'
        if e == -sp.Rational(3, 2):
            return '(1.0/(' + self._base(b) + '*std::sqrt(' + self._print(b) + ')))'
        return 'std::pow(' + self._print(b) + ', ' + self._print(e) + ')'

    def _print_atan2(self, expr):
        return 'std::atan2(' + self._print(expr.args[0]) + ', ' + self._print(expr.args[1]) + ')'

    def _print_Abs(self, expr):   # 인장 · 전단의 uv 넓이 (docs/GATE_STRETCH_SHEAR.md)
        return 'std::fabs(' + self._print(expr.args[0]) + ')'

    def _print_sin(self, expr):   # Bridson 굽힘의 sin(θ/4) (docs/GATE_BRIDSON.md)
        return 'std::sin(' + self._print(expr.args[0]) + ')'

    def _print_cos(self, expr):
        return 'std::cos(' + self._print(expr.args[0]) + ')'

    def _print_log(self, expr):   # Poisson 조건의 log|w| (docs/GATE_POISSON.md)
        return 'std::log(' + self._print(expr.args[0]) + ')'


def _signature(fname, inputs, outputs):
    args = ', '.join('const double* %s' % n for n, _ in inputs) + ', ' + \
           ', '.join('double* %s' % n for n, _ in outputs)
    return '\tstatic inline void %s(%s)' % (fname, args)


def emit_function(fname, inputs, outputs, printer):
    """inputs: [(C 이름, [기호 ...])], outputs: [(C 이름, [식 ...])] → (C++ 소스, 통계)."""
    flat = [e for _, es in outputs for e in es]
    t0 = time.perf_counter()
    repl, reduced = sp.cse(flat, symbols=sp.numbered_symbols('t'), order='none')
    t_cse = time.perf_counter() - t0

    lines = []
    lines.append(_signature(fname, inputs, outputs))
    lines.append('\t{')
    for n, syms in inputs:
        for i, s in enumerate(syms):
            lines.append('\t\tconst double %s = %s[%d];' % (printer._print(s), n, i))
    for s, e in repl:
        lines.append('\t\tconst double %s = %s;' % (printer._print(s), printer.doprint(e)))
    k = 0
    for n, es in outputs:
        for i in range(len(es)):
            lines.append('\t\t%s[%d] = %s;' % (n, i, printer.doprint(reduced[k])))
            k += 1
    lines.append('\t}')
    ops = sum(sp.count_ops(e) for _, e in repl) + sum(sp.count_ops(e) for e in reduced)
    return '\n'.join(lines), {'cse_seconds': round(t_cse, 2), 'temporaries': len(repl), 'ops': int(ops)}


# ---------------------------------------------------------------------------
# 유도
# ---------------------------------------------------------------------------
def _setup(spec):
    nDof = 3 * spec.N_NODES
    xs = sp.symbols('x0:%d' % nDof, real=True)
    X = [sp.Matrix(xs[3 * a:3 * a + 3]) for a in range(spec.N_NODES)]
    ev = []
    E = {}
    for name, ncol in spec.ELEM_DOUBLE_VARS:
        syms = sp.symbols('%s_0:%d' % (name, ncol), real=True)
        E[name] = syms
        ev.extend(syms)
    return nDof, xs, X, ev, E


def derive(spec, mode):
    nDof, xs, X, ev, E = _setup(spec)
    ctx = Ctx(mode)
    Cs = spec.condition(ctx, X, E)

    t0 = time.perf_counter()
    C_out, g_out, H_out = [], [], []
    for C in Cs:
        g = [sp.diff(C, x) for x in xs]
        H = [sp.diff(g[a], xs[b]) for a in range(nDof) for b in range(a, nDof)]   # 위 삼각, 행 우선
        C_out.append(C)
        g_out.extend(g)
        H_out.extend(H)
    if ctx.frozen:   # 미분이 끝난 뒤 길이 기호를 원래 길이로 되돌린다
        sub = dict(ctx.frozen)
        C_out = [e.xreplace(sub) for e in C_out]
        g_out = [e.xreplace(sub) for e in g_out]
        H_out = [e.xreplace(sub) for e in H_out]
    t_diff = time.perf_counter() - t0
    return xs, ev, C_out, g_out, H_out, t_diff


def derive_staged(spec, mode):
    """단계판 유도. 중간 변수가 없으면 None."""
    nDof, xs, X, ev, E = _setup(spec)
    ctx = Ctx(mode, staged=True)
    Cs = spec.condition(ctx, X, E)
    if not ctx.stages:
        return None

    t0 = time.perf_counter()
    us = [u for u, _ in ctx.stages]
    uex = [e for _, e in ctx.stages]
    xset = set(xs)
    direct = sorted({s for C in Cs for s in C.free_symbols if s in xset}, key=xs.index)
    ys = us + direct                    # 조건이 좌표를 직접 쓰면 그 좌표를 항등 중간 변수로 더한다
    yex = uex + direct
    m = len(ys)
    comps = []
    for C in Cs:
        G = [sp.diff(C, y) for y in ys]
        K = [[sp.diff(G[k], ys[l]) if l >= k else None for l in range(m)] for k in range(m)]
        for k in range(m):
            for l in range(k):
                K[k][l] = K[l][k]
        comps.append((C, G, K))
    if ctx.frozen:
        sub = dict(ctx.frozen)
        comps = [(C.xreplace(sub), [e.xreplace(sub) for e in G], [[e.xreplace(sub) for e in row] for row in K])
                 for C, G, K in comps]
    J = [[sp.diff(yex[k], x) for x in xs] for k in range(m)]
    U2 = [[[sp.diff(J[k][a], xs[b]) for b in range(nDof)] for a in range(nDof)] for k in range(m)]
    t_diff = time.perf_counter() - t0
    return dict(nDof=nDof, xs=xs, ev=ev, us=us, uex=uex, ys=ys, m=m, comps=comps, J=J, U2=U2, t_diff=t_diff)


def emit_staged(fname, d, printer, fused=False, gn=False, tension=False):
    """단계판 C++ 함수. 순서: 중간 변수와 J 의 0 아닌 성분 → 성분마다 C(u) · ∇ᵤC · ∇²ᵤC → T = ∇²ᵤC J → g · H.

    fused = True 이면 g · H 대신 원문 식 (6) 의 요소 기여를 바로 낸다(EvalImplicit, docs/GATE_ASSEMBLY.md 2 차):
      K' = α ∇ᵤC∇ᵤCᵀ + β ∇²ᵤC,  G' = β ∇ᵤC  (중간 변수 공간, α = h·kd + h²k, β = h²(kC + kd·Ċ))
      A  = Jᵀ K' J + Σₖ G'ₖ ∇²ₓuₖ  = α ∇C∇Cᵀ + β ∇²C   (위 삼각)
      f  = −(kC + kd·Ċ) ∇C                              (f + d)
      b  = −h k C ∇C − A v                              (= h(f + d) + h²(K + ∂d/∂x) v 를 정리한 꼴)

    gn = True(fused 와 함께)이면 Gauss–Newton 판(EvalImplicitGN, docs/GATE_GAUSS_NEWTON.md): ∇²C 를 0 으로 둔다.
      ∇²ᵤC · ∇²ₓu 를 쓰지 않고 g = Jᵀ ∇ᵤC 만 만든 뒤 좌표 공간에서 바로
      A = α Σ g gᵀ (위 삼각),  f = −Σ s g,  b = −Σ (h k C + α Ċ) g

    tension = True(fused 와 함께)이면 인장 C > 0 판(EvalImplicitTension, docs/GATE_TENSION.md): 정확판과 같은 길로 내되
      β = h² k max(C, 0) 로 둔다(∂d/∂x 를 버리고 눌린 성분의 C∇²C 를 버린다). b = −h k C ∇C − A v 는 그대로 맞다.
    """
    nDof, xs, ev, us, uex, m, J, U2 = d['nDof'], d['xs'], d['ev'], d['us'], d['uex'], d['m'], d['J'], d['U2']
    nHess = nDof * (nDof + 1) // 2
    inputs = [('x', list(xs))] + ([('ev', ev)] if ev else [])
    outputs = [('C', None), ('g', None), ('H', None)]
    if fused:
        vs = list(sp.symbols('v0:%d' % nDof, real=True))
        pk, pkd, ph = sp.symbols('p_k p_kd p_h', real=True)
        inputs = inputs + [('v', vs), ('prm', [pk, pkd, ph])]
        outputs = [('C', None), ('f', None), ('A', None), ('b', None)]
    P = printer.doprint
    lines = [_signature(fname, inputs, outputs), '\t{']
    for n, syms in inputs:
        for i, s in enumerate(syms):
            lines.append('\t\tconst double %s = %s[%d];' % (printer._print(s), n, i))
    ops = 0
    nTemp = 0
    t_cse = 0.

    def assign(sym, expr):
        nonlocal ops, nTemp
        lines.append('\t\tconst double %s = %s;' % (printer._print(sym), P(expr)))
        ops += int(sp.count_ops(expr))
        nTemp += 1

    def local(name, expr):
        """식이 기호 · 수 하나면 그대로, 아니면 지역 값 하나로 둔다."""
        expr = sp.sympify(expr)
        if expr.is_Symbol or expr.is_Number:
            return expr
        s = sp.Symbol(name)
        assign(s, expr)
        return s

    # 1) 중간 변수 u 와 J = ∂u/∂x 의 상수 아닌 성분(공통 부분식을 함께 묶는다)
    Jv = [[None] * nDof for _ in range(m)]
    jlist = []
    for k in range(m):
        for b in range(nDof):
            e = J[k][b]
            if e.is_Number:
                Jv[k][b] = e
            else:
                s = sp.Symbol('J%d_%d' % (k, b))
                Jv[k][b] = s
                jlist.append((s, e))
    block = list(uex) + [e for _, e in jlist]
    t0 = time.perf_counter()
    repl, red = sp.cse(block, symbols=sp.numbered_symbols('a'), order='none')
    t_cse += time.perf_counter() - t0
    lines.append('\t\t// 중간 변수 u 와 J = du/dx')
    for s, e in repl:
        assign(s, e)
    for u, e in zip(us, red[:len(us)]):
        assign(u, e)
    for (s, _), e in zip(jlist, red[len(us):]):
        assign(s, e)

    if fused:
        lines.append('\t\t// 계수: alpha = h kd + h^2 k')
        p_alpha = local('p_alpha', ph * pkd + ph * ph * pk)
        if tension:
            p_hhk = local('p_hhk', ph * ph * pk)
        elif not gn:
            p_hh = local('p_hh', ph * ph)
        p_mhk = local('p_mhk', -ph * pk)
        Aterms = [[[] for _ in range(nDof)] for _ in range(nDof)]   # 위 삼각 (a <= b), 성분마다 더한다
        fterms = [[] for _ in range(nDof)]
        bterms = [[] for _ in range(nDof)]

    for c, (C, G, K) in enumerate(d['comps']):
        # 2) C(u), ∇ᵤC, ∇²ᵤC (위 삼각)
        Gv = [None] * m
        Kv = [[None] * m for _ in range(m)]
        exprs, slots = [C], [('C', None)]
        for k in range(m):
            if G[k].is_Number:
                Gv[k] = G[k]
            else:
                exprs.append(G[k])
                slots.append(('G', k))
        for k in range(m if not gn else 0):   # Gauss–Newton 판은 ∇²ᵤC 를 쓰지 않는다
            for l in range(k, m):
                if K[k][l].is_Number:
                    Kv[k][l] = Kv[l][k] = K[k][l]
                else:
                    exprs.append(K[k][l])
                    slots.append(('K', (k, l)))
        t0 = time.perf_counter()
        repl, red = sp.cse(exprs, symbols=sp.numbered_symbols('b%d_' % c), order='none')
        t_cse += time.perf_counter() - t0
        lines.append('\t\t// 성분 %d: C(u), dC/du' % c if gn else '\t\t// 성분 %d: C(u), dC/du, d2C/du2' % c)
        for s, e in repl:
            assign(s, e)
        for (kind, idx), e in zip(slots, red):
            if kind == 'C':
                if fused:
                    Cc = sp.Symbol('Cv%d' % c)
                    assign(Cc, e)
                    lines.append('\t\tC[%d] = %s;' % (c, printer._print(Cc)))
                else:
                    lines.append('\t\tC[%d] = %s;' % (c, P(e)))
                    ops += int(sp.count_ops(e))
            elif kind == 'G':
                s = sp.Symbol('G%d_%d' % (c, idx))
                assign(s, e)
                Gv[idx] = s
            else:
                k, l = idx
                s = sp.Symbol('K%d_%d_%d' % (c, k, l))
                assign(s, e)
                Kv[k][l] = Kv[l][k] = s
        if gn:
            # 3'') Gauss–Newton: g = Jᵀ ∇ᵤC, Ċ = g · v, A += α g gᵀ, f += −s g, b += −(h k C + α Ċ) g
            lines.append('\t\t// g = J^T dC/du, Cdot = g.v, s = k C + kd Cdot, q = -(h k C + alpha Cdot)')
            gv = [local('g%d_%d' % (c, a), sp.Add(*[Jv[k][a] * Gv[k] for k in range(m) if Jv[k][a] != 0 and Gv[k] != 0]))
                  for a in range(nDof)]
            Cd = local('Cd%d' % c, sp.Add(*[gv[a] * vs[a] for a in range(nDof) if gv[a] != 0]))
            sc = local('s%d' % c, pk * Cc + pkd * Cd)
            msc = local('ms%d' % c, -sc)
            qc = local('q%d' % c, p_mhk * Cc - p_alpha * Cd)
            for a in range(nDof):
                fterms[a].append(msc * gv[a])
                bterms[a].append(qc * gv[a])
            lines.append('\t\t// alpha g (A = alpha g g^T)')
            ag = [local('ag%d_%d' % (c, a), p_alpha * gv[a]) for a in range(nDof)]
            for a in range(nDof):
                for b in range(a, nDof):
                    if ag[a] != 0 and gv[b] != 0:
                        Aterms[a][b].append(ag[a] * gv[b])
            continue
        if fused:
            # 3') g = Jᵀ ∇ᵤC, Ċ = g · v, 계수, K' = α G Gᵀ + β K, G' = β G
            lines.append('\t\t// g = J^T dC/du, Cdot = g.v, s = k C + kd Cdot, beta = h^2 s')
            gv = [local('g%d_%d' % (c, a), sp.Add(*[Jv[k][a] * Gv[k] for k in range(m) if Jv[k][a] != 0 and Gv[k] != 0]))
                  for a in range(nDof)]
            Cd = local('Cd%d' % c, sp.Add(*[gv[a] * vs[a] for a in range(nDof) if gv[a] != 0]))
            sc = local('s%d' % c, pk * Cc + pkd * Cd)
            if tension:
                # 인장 C > 0 판: beta = h^2 k max(C, 0)
                Cp = sp.Symbol('Cp%d' % c)
                lines.append('\t\tconst double %s = %s > 0.0 ? %s : 0.0;' % (printer._print(Cp), printer._print(Cc), printer._print(Cc)))
                ops += 1
                nTemp += 1
                bc = local('beta%d' % c, p_hhk * Cp)
            else:
                bc = local('beta%d' % c, p_hh * sc)
            msc = local('ms%d' % c, -sc)
            mhkC = local('mhkC%d' % c, p_mhk * Cc)
            for a in range(nDof):
                fterms[a].append(msc * gv[a])
                bterms[a].append(mhkC * gv[a])
            lines.append("\t\t// K' = alpha G G^T + beta K,  G' = beta G  (중간 변수 공간)")
            aG = [local('aG%d_%d' % (c, k), p_alpha * Gv[k]) for k in range(m)]
            Gp = [local('Gp%d_%d' % (c, k), bc * Gv[k])
                  if any(U2[k][a][b] != 0 for a in range(nDof) for b in range(nDof)) else sp.S.Zero
                  for k in range(m)]   # ∇²ₓuₖ = 0 이면 쓰이지 않는다
            Kp = [[None] * m for _ in range(m)]
            for k in range(m):
                for l in range(k, m):
                    Kp[k][l] = Kp[l][k] = local('Kp%d_%d_%d' % (c, k, l), aG[k] * Gv[l] + bc * Kv[k][l])
            Kv, Gv = Kp, Gp     # 아래 T · H 단계를 K' · G' 로 그대로 돈다 → A
        # 3) T = ∇²ᵤC J  (m × nDof)
        lines.append('\t\t// T = (d2C/du2) J')
        T = [[0] * nDof for _ in range(m)]
        for k in range(m):
            for b in range(nDof):
                if not any(Jv[k][a] != 0 for a in range(b + 1)):   # H 위 삼각(a <= b)에서 안 쓰인다
                    continue
                e = sp.Add(*[Kv[k][l] * Jv[l][b] for l in range(m) if Jv[l][b] != 0 and Kv[k][l] != 0])
                if e == 0:
                    continue
                if e.is_Symbol or e.is_Number:
                    T[k][b] = e
                else:
                    s = sp.Symbol('T%d_%d_%d' % (c, k, b))
                    assign(s, e)
                    T[k][b] = s
        if fused:
            # 4') A 의 이 성분 몫 = Jᵀ T' + Σₖ G'ₖ ∇²ₓuₖ
            for a in range(nDof):
                for b in range(a, nDof):
                    Aterms[a][b] += [Jv[k][a] * T[k][b] for k in range(m) if Jv[k][a] != 0 and T[k][b] != 0]
                    Aterms[a][b] += [Gv[k] * U2[k][a][b] for k in range(m) if U2[k][a][b] != 0 and Gv[k] != 0]
            continue
        # 4) g = Jᵀ ∇ᵤC,  H = Jᵀ T + Σₖ ∂C/∂uₖ ∇²ₓuₖ
        lines.append('\t\t// g = J^T dC/du,  H = J^T T + sum_k dC/du_k d2u_k/dx2 (위 삼각, 행 우선)')
        for a in range(nDof):
            e = sp.Add(*[Jv[k][a] * Gv[k] for k in range(m) if Jv[k][a] != 0 and Gv[k] != 0])
            lines.append('\t\tg[%d] = %s;' % (c * nDof + a, P(e)))
            ops += int(sp.count_ops(e))
        h = 0
        for a in range(nDof):
            for b in range(a, nDof):
                terms = [Jv[k][a] * T[k][b] for k in range(m) if Jv[k][a] != 0 and T[k][b] != 0]
                terms += [Gv[k] * U2[k][a][b] for k in range(m) if U2[k][a][b] != 0 and Gv[k] != 0]
                e = sp.Add(*terms)
                lines.append('\t\tH[%d] = %s;' % (c * nHess + h, P(e)))
                ops += int(sp.count_ops(e))
                h += 1
    if gn:
        # 5'') A = α Σ g gᵀ (위 삼각) → 출력, f = Σ −s g,  b = Σ q g  (A v 를 쓰지 않는다)
        lines.append('\t\t// A = alpha sum g g^T (위 삼각, 행 우선),  f = sum -s g,  b = sum q g')
        h = 0
        for a in range(nDof):
            for b in range(a, nDof):
                e = sp.Add(*Aterms[a][b])
                lines.append('\t\tA[%d] = %s;' % (h, P(e)))
                ops += int(sp.count_ops(e))
                h += 1
        for a in range(nDof):
            e = sp.Add(*fterms[a])
            lines.append('\t\tf[%d] = %s;' % (a, P(e)))
            ops += int(sp.count_ops(e))
        for a in range(nDof):
            e = sp.Add(*bterms[a])
            lines.append('\t\tb[%d] = %s;' % (a, P(e)))
            ops += int(sp.count_ops(e))
    elif fused:
        # 5') A (위 삼각) → 출력, f = Σ −s g,  b = Σ −h k C g − A v
        lines.append('\t\t// A = alpha g g^T + beta H (위 삼각, 행 우선),  f = -s g,  b = -h k C g - A v')
        Al = [[None] * nDof for _ in range(nDof)]
        h = 0
        for a in range(nDof):
            for b in range(a, nDof):
                Al[a][b] = Al[b][a] = local('A%d_%d' % (a, b), sp.Add(*Aterms[a][b]))
                lines.append('\t\tA[%d] = %s;' % (h, P(Al[a][b])))
                h += 1
        for a in range(nDof):
            e = sp.Add(*fterms[a])
            lines.append('\t\tf[%d] = %s;' % (a, P(e)))
            ops += int(sp.count_ops(e))
        for a in range(nDof):
            e = sp.Add(*bterms[a]) - sp.Add(*[Al[a][j] * vs[j] for j in range(nDof) if Al[a][j] != 0])
            lines.append('\t\tb[%d] = %s;' % (a, P(e)))
            ops += int(sp.count_ops(e))
    lines.append('\t}')
    return '\n'.join(lines), {'cse_seconds': round(t_cse, 2), 'temporaries': nTemp, 'ops': ops, 'stages': m}


def emit_engine(fname, d, printer, kind, damped=True):
    """엔진(cfmsDrape CPU)용 변형 — docs/GATE_LINK_L1.md. kind = 'exact' | 'gn' | 'tension' | 'tension_rate'.
    damped=False 는 감쇠 없는 판(GN 만): 입력 v · prm[1] = kd 와 출력 V 를 빼고 s = k C 로 둔다.

    출력: C[kComp], f[kDof] = f + d, 노드쌍 블록 J[9 kPairs] · V[9 kPairs]
    (지역 꼭짓점 a < b, 사전순, 블록 행 우선: J[9p + 3i + j] = J_{3a+i, 3b+j}).
      J = −k Σ ∇C∇Cᵀ − Σ β'_c ∇²C_c   (β'_c = kC + kdĊ 정확 · 0 GN · k max(C, 0) T ·
                                       k max(C, 0) + kd max(Ċ, 0) 엔진 처방 T — 'tension_rate', L1-5)
      V = −kd Σ ∇C∇Cᵀ
    대각 블록은 내지 않는다(엔진은 행 합 규약 −Σ 로 정한다). prm = {k, kd} — h 는 엔진이 곱한다.
    엔진의 A = M − hV − h²J, b = h(f + d) + h²Jv 가 짝이 되는 EvalImplicit* 의 A · b 와 같다.
    emit_staged 와 같은 단계판 유도를 읽지만 따로 쓴 함수다 — 기존 변형의 출력을 바꾸지 않으려고 고치지 않았다.
    ∇C∇Cᵀ 몫은 좌표 공간의 g 로 한 번 만들어 J · V 가 함께 쓰고, ∇²C 몫은 중간 변수 공간에서 잇는다.
    """
    nDof, xs, ev, us, uex, m, J, U2 = d['nDof'], d['xs'], d['ev'], d['us'], d['uex'], d['m'], d['J'], d['U2']
    nNodes = nDof // 3
    pairs = [(a, b) for a in range(nNodes) for b in range(a + 1, nNodes)]

    def need(a, b):   # 좌표 a, b 가 노드쌍 블록 안인가(행 노드 < 열 노드)
        return a // 3 < b // 3

    vs = list(sp.symbols('v0:%d' % nDof, real=True))
    pk, pkd = sp.symbols('p_k p_kd', real=True)
    if damped:
        inputs = [('x', list(xs))] + ([('ev', ev)] if ev else []) + [('v', vs), ('prm', [pk, pkd])]
        outputs = [('C', None), ('f', None), ('J', None), ('V', None)]
    else:
        assert kind == 'gn', '감쇠 없는 판은 GN 만'
        inputs = [('x', list(xs))] + ([('ev', ev)] if ev else []) + [('prm', [pk])]
        outputs = [('C', None), ('f', None), ('J', None)]
    P = printer.doprint
    lines = [_signature(fname, inputs, outputs), '\t{']
    for n, syms in inputs:
        for i, s in enumerate(syms):
            lines.append('\t\tconst double %s = %s[%d];' % (printer._print(s), n, i))
    ops = 0
    nTemp = 0
    t_cse = 0.

    def assign(sym, expr):
        nonlocal ops, nTemp
        lines.append('\t\tconst double %s = %s;' % (printer._print(sym), P(expr)))
        ops += int(sp.count_ops(expr))
        nTemp += 1

    def local(name, expr):
        expr = sp.sympify(expr)
        if expr.is_Symbol or expr.is_Number:
            return expr
        s = sp.Symbol(name)
        assign(s, expr)
        return s

    # 1) 중간 변수 u 와 J = ∂u/∂x 의 상수 아닌 성분 (emit_staged 와 같다)
    Jv = [[None] * nDof for _ in range(m)]
    jlist = []
    for k in range(m):
        for b in range(nDof):
            e = J[k][b]
            if e.is_Number:
                Jv[k][b] = e
            else:
                s = sp.Symbol('J%d_%d' % (k, b))
                Jv[k][b] = s
                jlist.append((s, e))
    block = list(uex) + [e for _, e in jlist]
    t0 = time.perf_counter()
    repl, red = sp.cse(block, symbols=sp.numbered_symbols('a'), order='none')
    t_cse += time.perf_counter() - t0
    lines.append('\t\t// 중간 변수 u 와 J = du/dx')
    for s, e in repl:
        assign(s, e)
    for u, e in zip(us, red[:len(us)]):
        assign(u, e)
    for (s, _), e in zip(jlist, red[len(us):]):
        assign(s, e)

    lines.append('\t\t// 계수: -k, -kd' if damped else '\t\t// 계수: -k')
    p_mk = local('p_mk', -pk)
    p_mkd = local('p_mkd', -pkd) if damped else None
    gn = kind == 'gn'
    Pterms = [[[] for _ in range(nDof)] for _ in range(nDof)]   # Σ g gᵀ, 노드쌍 블록만
    Hterms = [[[] for _ in range(nDof)] for _ in range(nDof)]   # −Σ β' ∇²C, 노드쌍 블록만
    fterms = [[] for _ in range(nDof)]

    for c, (C, G, K) in enumerate(d['comps']):
        # 2) C(u), ∇ᵤC, ∇²ᵤC (GN 판은 ∇²ᵤC 를 쓰지 않는다)
        Gv = [None] * m
        Kv = [[None] * m for _ in range(m)]
        exprs, slots = [C], [('C', None)]
        for k in range(m):
            if G[k].is_Number:
                Gv[k] = G[k]
            else:
                exprs.append(G[k])
                slots.append(('G', k))
        for k in range(m if not gn else 0):
            for l in range(k, m):
                if K[k][l].is_Number:
                    Kv[k][l] = Kv[l][k] = K[k][l]
                else:
                    exprs.append(K[k][l])
                    slots.append(('K', (k, l)))
        t0 = time.perf_counter()
        repl, red = sp.cse(exprs, symbols=sp.numbered_symbols('b%d_' % c), order='none')
        t_cse += time.perf_counter() - t0
        lines.append('\t\t// 성분 %d: C(u), dC/du' % c if gn else '\t\t// 성분 %d: C(u), dC/du, d2C/du2' % c)
        for s, e in repl:
            assign(s, e)
        for (kd_, idx), e in zip(slots, red):
            if kd_ == 'C':
                Cc = sp.Symbol('Cv%d' % c)
                assign(Cc, e)
                lines.append('\t\tC[%d] = %s;' % (c, printer._print(Cc)))
            elif kd_ == 'G':
                s = sp.Symbol('G%d_%d' % (c, idx))
                assign(s, e)
                Gv[idx] = s
            else:
                k, l = idx
                s = sp.Symbol('K%d_%d_%d' % (c, k, l))
                assign(s, e)
                Kv[k][l] = Kv[l][k] = s
        # 3) g = Jᵀ ∇ᵤC, Ċ = g · v, s = k C + kd Ċ, f += −s g, Σ g gᵀ (노드쌍 블록만)
        lines.append('\t\t// g = J^T dC/du, Cdot = g.v, s = k C + kd Cdot' if damped else '\t\t// g = J^T dC/du, s = k C')
        gv = [local('g%d_%d' % (c, a), sp.Add(*[Jv[k][a] * Gv[k] for k in range(m) if Jv[k][a] != 0 and Gv[k] != 0]))
              for a in range(nDof)]
        if damped:
            Cd = local('Cd%d' % c, sp.Add(*[gv[a] * vs[a] for a in range(nDof) if gv[a] != 0]))
            sc = local('s%d' % c, pk * Cc + pkd * Cd)
        else:
            sc = local('s%d' % c, pk * Cc)
        msc = local('ms%d' % c, -sc)
        for a in range(nDof):
            fterms[a].append(msc * gv[a])
            for b in range(a + 1, nDof):
                if need(a, b) and gv[a] != 0 and gv[b] != 0:
                    Pterms[a][b].append(gv[a] * gv[b])
        if gn:
            continue
        # 4) J 의 ∇²C 몫: β'' = −β' (정확 −s, T −k max(C, 0)), K'' = β'' ∇²ᵤC, G'' = β'' ∇ᵤC
        if kind == 'tension':
            Cp = sp.Symbol('Cp%d' % c)
            lines.append('\t\tconst double %s = %s > 0.0 ? %s : 0.0;' % (printer._print(Cp), printer._print(Cc), printer._print(Cc)))
            ops += 1
            nTemp += 1
            bc = local('mbeta%d' % c, p_mk * Cp)
        elif kind == 'tension_rate':   # 엔진 손 코드의 처방: 감쇠의 Ċ∇²C 몫도 Ċ > 0 일 때만 넣는다
            Cp = sp.Symbol('Cp%d' % c)
            lines.append('\t\tconst double %s = %s > 0.0 ? %s : 0.0;' % (printer._print(Cp), printer._print(Cc), printer._print(Cc)))
            Cdp = sp.Symbol('Cdp%d' % c)
            lines.append('\t\tconst double %s = %s > 0.0 ? %s : 0.0;' % (printer._print(Cdp), printer._print(Cd), printer._print(Cd)))
            ops += 2
            nTemp += 2
            bc = local('mbeta%d' % c, p_mk * Cp + p_mkd * Cdp)
        else:
            bc = msc
        lines.append("\t\t// K'' = beta'' d2C/du2,  G'' = beta'' dC/du  (중간 변수 공간)")
        Gp = [local('Gp%d_%d' % (c, k), bc * Gv[k])
              if any(U2[k][a][b] != 0 for a in range(nDof) for b in range(nDof) if need(a, b)) else sp.S.Zero
              for k in range(m)]
        Kp = [[None] * m for _ in range(m)]
        for k in range(m):
            for l in range(k, m):
                Kp[k][l] = Kp[l][k] = local('Kp%d_%d_%d' % (c, k, l), bc * Kv[k][l])
        lines.append("\t\t// T = K'' J (노드쌍 블록에 쓰이는 열만)")
        T = [[0] * nDof for _ in range(m)]
        for k in range(m):
            for b in range(nDof):
                if not any(Jv[k][a] != 0 for a in range(nDof) if need(a, b)):
                    continue
                e = sp.Add(*[Kp[k][l] * Jv[l][b] for l in range(m) if Jv[l][b] != 0 and Kp[k][l] != 0])
                if e == 0:
                    continue
                if e.is_Symbol or e.is_Number:
                    T[k][b] = e
                else:
                    s = sp.Symbol('T%d_%d_%d' % (c, k, b))
                    assign(s, e)
                    T[k][b] = s
        for a in range(nDof):
            for b in range(a + 1, nDof):
                if not need(a, b):
                    continue
                Hterms[a][b] += [Jv[k][a] * T[k][b] for k in range(m) if Jv[k][a] != 0 and T[k][b] != 0]
                Hterms[a][b] += [Gp[k] * U2[k][a][b] for k in range(m) if U2[k][a][b] != 0 and Gp[k] != 0]

    # 5) f = Σ −s g,  블록: P = Σ g gᵀ,  J = −k P + Σ (−β' ∇²C),  V = −kd P
    lines.append('\t\t// f = sum -s g;  J = -k sum g g^T - sum beta\' d2C,  V = -kd sum g g^T (노드쌍 블록, 행 우선)' if damped
                 else '\t\t// f = sum -s g;  J = -k sum g g^T (노드쌍 블록, 행 우선)')
    for a in range(nDof):
        e = sp.Add(*fterms[a])
        lines.append('\t\tf[%d] = %s;' % (a, P(e)))
        ops += int(sp.count_ops(e))
    for p, (na, nb) in enumerate(pairs):
        for i in range(3):
            for j in range(3):
                a, b = 3 * na + i, 3 * nb + j
                Ps = local('P%d_%d' % (a, b), sp.Add(*Pterms[a][b]))
                eJ = p_mk * Ps + sp.Add(*Hterms[a][b])
                lines.append('\t\tJ[%d] = %s;' % (9 * p + 3 * i + j, P(eJ)))
                ops += int(sp.count_ops(eJ))
                if damped:
                    eV = p_mkd * Ps
                    lines.append('\t\tV[%d] = %s;' % (9 * p + 3 * i + j, P(eV)))
                    ops += int(sp.count_ops(eV))
    lines.append('\t}')
    return '\n'.join(lines), {'cse_seconds': round(t_cse, 2), 'temporaries': nTemp, 'ops': ops, 'stages': m}


def load_spec(path):
    s = importlib.util.spec_from_file_location('sxspec', path)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec')
    ap.add_argument('--out', default='generated')
    ap.add_argument('--results', default='results')
    a = ap.parse_args()

    spec = load_spec(a.spec)
    specText = open(a.spec, 'rb').read()
    printer = Printer()
    nDof = 3 * spec.N_NODES
    nHess = nDof * (nDof + 1) // 2
    nComp = spec.N_COMP

    funcs, stats = [], {'spec': os.path.basename(a.spec), 'spec_sha256': hashlib.sha256(specText).hexdigest()[:16],
                        'sympy': sp.__version__, 'variants': {}}
    for mode, fname in (('exact', 'Eval'), ('frozen', 'EvalFrozen')):
        xs, ev, C, g, H, t_diff = derive(spec, mode)
        inputs = [('x', list(xs))] + ([('ev', ev)] if ev else [])
        src, st = emit_function(fname, inputs, [('C', C), ('g', g), ('H', H)], printer)
        st['diff_seconds'] = round(t_diff, 2)
        st['body_sha256'] = hashlib.sha256(src.encode('utf-8')).hexdigest()
        stats['variants'][mode] = st
        funcs.append(src)
        print('%-7s diff %.1f s, cse %.1f s, temporaries %d, ops %d'
              % (mode, t_diff, st['cse_seconds'], st['temporaries'], st['ops']))
    d = derive_staged(spec, 'exact')
    if d is not None:
        src, st = emit_staged('EvalStaged', d, printer)
        st['diff_seconds'] = round(d['t_diff'], 2)
        st['body_sha256'] = hashlib.sha256(src.encode('utf-8')).hexdigest()
        stats['variants']['staged'] = st
        funcs.append(src)
        print('%-7s diff %.1f s, cse %.1f s, temporaries %d, ops %d, stages %d'
              % ('staged', d['t_diff'], st['cse_seconds'], st['temporaries'], st['ops'], st['stages']))
        src, st = emit_staged('EvalImplicit', d, printer, fused=True)
        st['diff_seconds'] = round(d['t_diff'], 2)
        st['body_sha256'] = hashlib.sha256(src.encode('utf-8')).hexdigest()
        stats['variants']['implicit'] = st
        funcs.append('\t// EvalImplicit: v[kDof] = 꼭짓점 속도, prm = {k, kd, h} → C[kComp], f[kDof] = f + d,\n'
                     '\t//   A[kHess] = −h ∂f/∂v − h² ∂f/∂x (위 삼각, 행 우선), b[kDof] = h(f + d) + h² ∂f/∂x v.\n'
                     '\t//   원문 식 (6) 의 요소 기여, 감쇠는 비대칭 항을 버린 원문 꼴. docs/GATE_ASSEMBLY.md 2 차.')
        funcs.append(src)
        print('%-8s cse %.1f s, temporaries %d, ops %d' % ('implicit', st['cse_seconds'], st['temporaries'], st['ops']))
        src, st = emit_staged('EvalImplicitGN', d, printer, fused=True, gn=True)
        st['body_sha256'] = hashlib.sha256(src.encode('utf-8')).hexdigest()
        stats['variants']['implicit_gn'] = st
        funcs.append('\t// EvalImplicitGN: EvalImplicit 와 입출력이 같고 ∇²C 를 0 으로 둔 Gauss–Newton 판.\n'
                     '\t//   f + d 는 같고, A = (h kd + h² k) Σ ∇C∇Cᵀ (양의 준정부호), b = h(f + d) − h² k Σ ∇C (∇C·v).\n'
                     '\t//   docs/GATE_GAUSS_NEWTON.md.')
        funcs.append(src)
        print('%-8s cse %.1f s, temporaries %d, ops %d' % ('impl_gn', st['cse_seconds'], st['temporaries'], st['ops']))
        if getattr(spec, 'IMPLICIT_TENSION', False):
            src, st = emit_staged('EvalImplicitTension', d, printer, fused=True, tension=True)
            st['body_sha256'] = hashlib.sha256(src.encode('utf-8')).hexdigest()
            stats['variants']['implicit_tension'] = st
            funcs.append('\t// EvalImplicitTension: EvalImplicit 와 입출력이 같은 인장 C > 0 판. f + d 는 같고,\n'
                         '\t//   A = (h kd + h² k) Σ ∇C∇Cᵀ + h² k Σ max(C, 0) ∇²C (∇²C 가 준정부호면 준정부호), ∂d/∂x 는 버린다.\n'
                         '\t//   b = h(f + d) + h² K_T v. docs/GATE_TENSION.md.')
            funcs.append(src)
            print('%-8s cse %.1f s, temporaries %d, ops %d' % ('impl_ten', st['cse_seconds'], st['temporaries'], st['ops']))

    os.makedirs(a.out, exist_ok=True)
    nodeVars = ''.join('L"%s", ' % n for n, _ in spec.NODE_VARS)
    nodeCols = ''.join('%d, ' % c for _, c in spec.NODE_VARS)
    elemVars = ''.join('L"%s", ' % n for n, _ in spec.ELEM_DOUBLE_VARS) or '0, '
    elemCols = ''.join('%d, ' % c for _, c in spec.ELEM_DOUBLE_VARS) or '0, '
    nEv = sum(c for _, c in spec.ELEM_DOUBLE_VARS)
    hdr = []
    hdr.append('// 자동 생성 — 손으로 고치지 않는다.')
    hdr.append('// 생성기 codegen/sxgen.py · 명세 codegen/specs/%s (sha256 %s) · SymPy %s'
               % (stats['spec'], stats['spec_sha256'], sp.__version__))
    hdr.append('// %s' % spec.DESCRIPTION)
    hdr.append('#pragma once')
    hdr.append('#include <cmath>')
    hdr.append('')
    hdr.append('namespace sxgen {')
    hdr.append('')
    hdr.append('struct %s {' % spec.NAME)
    hdr.append('\tenum { kNodes = %d, kDof = %d, kComp = %d, kHess = %d, kElemValues = %d };'
               % (spec.N_NODES, nDof, nComp, nHess, nEv))
    hdr.append('\t// SXMesh 변수 이름. 노드 쪽 double 변수 · 스텐실(정수, 행마다 꼭짓점 %d개) · 스텐실마다 double 변수' % spec.N_NODES)
    hdr.append('\tstatic const wchar_t* PositionVar() { return L"%s"; }' % spec.POSITION_VAR)
    hdr.append('\tstatic const wchar_t* VelocityVar() { return L"%s"; }' % spec.VELOCITY_VAR)
    hdr.append('\tstatic const wchar_t* StencilVar()  { return L"%s"; }' % spec.STENCIL_VAR)
    hdr.append('\tenum { kNodeVars = %d, kElemVars = %d };' % (len(spec.NODE_VARS), len(spec.ELEM_DOUBLE_VARS)))
    hdr.append('\tstatic const wchar_t* NodeVar(int i) { static const wchar_t* n[] = { %s}; return n[i]; }' % nodeVars)
    hdr.append('\tstatic int NodeVarCols(int i)      { static const int c[] = { %s}; return c[i]; }' % nodeCols)
    hdr.append('\tstatic const wchar_t* ElemVar(int i) { static const wchar_t* n[] = { %s}; return n[i]; }' % elemVars)
    hdr.append('\tstatic int ElemVarCols(int i)      { static const int c[] = { %s}; return c[i]; }' % elemCols)
    hdr.append('')
    hdr.append('\t// x[kDof] = 스텐실 꼭짓점 좌표(노드 순서대로 xyz), ev[kElemValues] = 스텐실 값(선언 순서)')
    hdr.append('\t// C[kComp], g[kComp*kDof], H[kComp*kHess] (∇²C 의 위 삼각, 행 우선)')
    hdr.extend(funcs)
    hdr.append('};')
    hdr.append('')
    hdr.append('} // namespace sxgen')
    outPath = os.path.join(a.out, spec.NAME + '.h')
    data = ('﻿' + '\n'.join(hdr) + '\n').encode('utf-8')
    open(outPath, 'wb').write(data)
    stats['header'] = outPath.replace('\\', '/')
    stats['header_bytes'] = len(data)
    stats['header_lines'] = data.count(b'\n')
    os.makedirs(a.results, exist_ok=True)
    print('wrote %s (%d lines, %d bytes)' % (outPath, stats['header_lines'], len(data)))

    # 엔진(cfmsDrape CPU)용 변형은 따로 된 헤더에 쓴다 — 위 하네스 헤더는 그대로 둔다. docs/GATE_LINK_L1.md.
    # 하네스 헤더를 다 쓴 뒤에 만든다(기존 변형의 출력에 어떤 식으로도 끼어들지 않도록).
    if d is not None:
        kinds = [('exact', 'EvalEngine', 'engine',
                  'J = −k Σ ∇C∇Cᵀ − Σ (kC + kdĊ) ∇²C. EvalImplicit 의 짝.'),
                 ('gn', 'EvalEngineGN', 'engine_gn',
                  'J = −k Σ ∇C∇Cᵀ (Gauss–Newton). EvalImplicitGN 의 짝.')]
        if getattr(spec, 'IMPLICIT_TENSION', False):
            kinds.append(('tension', 'EvalEngineTension', 'engine_tension',
                          'J = −k Σ ∇C∇Cᵀ − k Σ max(C, 0) ∇²C (인장 C > 0 판). EvalImplicitTension 의 짝.'))
        # 감쇠 없는 GN 은 기존 변형 뒤에 붙인다(앞 함수들의 본문 · 해시는 그대로). docs/GATE_LINK_L4.md 「2 차 (다)」.
        kinds.append(('gn_undamped', 'EvalEngineGNUndamped', 'engine_gn_undamped',
                      'J = −k Σ ∇C∇Cᵀ, 감쇠 없음 — (x, ev, prm = {k}) → C · f · J. EvalEngineGN(v = 0, kd = 0) 과 같다.'))
        # 엔진 처방 T 판도 기존 변형 뒤에 붙인다(앞 함수들의 본문 · 해시는 그대로). docs/GATE_LINK_L1.md L1-5.
        if getattr(spec, 'IMPLICIT_TENSION', False):
            kinds.append(('tension_rate', 'EvalEngineTensionRate', 'engine_tension_rate',
                          'J = −k Σ ∇C∇Cᵀ − Σ (k max(C, 0) + kd max(Ċ, 0)) ∇²C — 엔진 cfmsDrape 인장 · 전단 손 코드의 처방'
                          '(감쇠의 Ċ∇²C 몫을 Ċ > 0 일 때만 넣는다). C · f · V 는 EvalEngineTension 과 같다.'))
        efuncs = []
        for kind, fname, key, note in kinds:
            damped = kind != 'gn_undamped'
            src, st = emit_engine(fname, d, printer, kind if damped else 'gn', damped=damped)
            st['body_sha256'] = hashlib.sha256(src.encode('utf-8')).hexdigest()
            stats['variants'][key] = st
            efuncs.append('\t// %s: %s' % (fname, note))
            efuncs.append(src)
            print('%-8s cse %.1f s, temporaries %d, ops %d' % (key, st['cse_seconds'], st['temporaries'], st['ops']))
        nPairs = spec.N_NODES * (spec.N_NODES - 1) // 2
        pa = [a for a in range(spec.N_NODES) for b in range(a + 1, spec.N_NODES)]
        pb = [b for a in range(spec.N_NODES) for b in range(a + 1, spec.N_NODES)]
        eh = []
        eh.append('// 자동 생성 — 손으로 고치지 않는다.')
        eh.append('// 생성기 codegen/sxgen.py · 명세 codegen/specs/%s (sha256 %s) · SymPy %s'
                  % (stats['spec'], stats['spec_sha256'], sp.__version__))
        eh.append('// %s' % spec.DESCRIPTION)
        eh.append('// 엔진(cfmsDrape CPU)용 변형 — f + d 와 노드쌍 블록 J · V. docs/GATE_LINK_L1.md.')
        eh.append('#pragma once')
        eh.append('#include <cmath>')
        eh.append('')
        eh.append('namespace sxgen {')
        eh.append('')
        eh.append('struct %sEngine {' % spec.NAME)
        eh.append('\tenum { kNodes = %d, kDof = %d, kComp = %d, kPairs = %d, kElemValues = %d };'
                  % (spec.N_NODES, nDof, nComp, nPairs, nEv))
        eh.append('\t// 노드쌍 p = 지역 꼭짓점 (PairA(p), PairB(p)), a < b, 사전순')
        eh.append('\tstatic int PairA(int p) { static const int t[] = { %s}; return t[p]; }' % ''.join('%d, ' % x for x in pa))
        eh.append('\tstatic int PairB(int p) { static const int t[] = { %s}; return t[p]; }' % ''.join('%d, ' % x for x in pb))
        eh.append('')
        eh.append('\t// x[kDof] = 스텐실 꼭짓점 좌표(노드 순서대로 xyz), ev[kElemValues] = 스텐실 값, v[kDof] = 꼭짓점 속도,')
        eh.append('\t// prm = {k, kd} → C[kComp], f[kDof] = f + d, J[9 kPairs] · V[9 kPairs] = 노드쌍 p 의 3 × 3 블록(행 우선):')
        eh.append('\t//   J[9p + 3i + j] = ∂f_{3a+i}/∂x_{3b+j}, V 는 ∂/∂v (a = PairA(p), b = PairB(p)), V = −kd Σ ∇C∇Cᵀ.')
        eh.append('\t//   대각 블록은 내지 않는다 — −(그 행의 이웃 블록 합)이다(엔진 행 합 규약). 요소 행렬은 대칭이라 (b, a) 블록 = (a, b)ᵀ.')
        eh.append('\t//   엔진 A = M − h V − h² J, b = h (f + d) + h² J v 가 짝이 되는 EvalImplicit* 의 A · b 와 같다.')
        eh.append('\t// EvalEngineGNUndamped 는 감쇠 없는 GN: v · prm[1] = kd · V 가 없고 (x, [ev,] prm = {k}) → C · f · J.')
        eh.extend(efuncs)
        eh.append('};')
        eh.append('')
        eh.append('} // namespace sxgen')
        ePath = os.path.join(a.out, spec.NAME + '_Engine.h')
        edata = ('﻿' + '\n'.join(eh) + '\n').encode('utf-8')
        open(ePath, 'wb').write(edata)
        stats['engine_header'] = ePath.replace('\\', '/')
        stats['engine_header_bytes'] = len(edata)
        stats['engine_header_lines'] = edata.count(b'\n')
        print('wrote %s (%d lines, %d bytes)' % (ePath, stats['engine_header_lines'], len(edata)))
    rp = os.path.join(a.results, 'codegen_%s.json' % spec.NAME)
    json.dump(stats, open(rp, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
    print('wrote %s' % rp)


if __name__ == '__main__':
    sys.exit(main())
