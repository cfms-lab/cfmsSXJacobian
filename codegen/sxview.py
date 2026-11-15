"""조건 함수 명세를 텐서 꼴로 보여 주는 시제품 (주피터용). 생성기 sxgen.py 는 건드리지 않고 읽지도 않는다.

명세의 condition(ctx, X, E) 를 「보기용 ctx」로 한 번 부른다. ctx.stage 가 받은 중간 벡터 d_i 가 좌표의 등방 1 차식
(d_i = Σ_m c_im x_m, c_im 은 스칼라)이고, 조건이 그 길이 |d_i| 만의 함수이면 미분이 꼭짓점 쌍 3 × 3 블록의 닫힌 꼴로 나온다:

  ∂C/∂x_m      = Σ_i c_im (∂C/∂|d_i|) d̂_i
  ∂²C/∂x_m∂x_n = Σ_i c_im c_in [ (∂C/∂|d_i|)/|d_i| (I − d̂_i d̂_iᵀ) + (∂²C/∂|d_i|²) d̂_i d̂_iᵀ ]
                + Σ_{i≠j} c_im c_jn (∂²C/∂|d_i|∂|d_j|) d̂_i d̂_jᵀ

view(): 중간 벡터의 길이 · 내적만 쓰는 명세 — 인장(baraff_stretch.py) · Provot 전단(provot_shear.py) · Baraff–Witkin 전단(baraff_shear.py).
bend_view(): θ = atan2 를 한 번 쓰는 굽힘 명세 — baraff_bend.py · bridson_bend.py (아래 Hinge · BendView).

check() 는 화면의 텐서 식을 숫자로 평가해, 같은 condition() 을 성분으로 풀어 SymPy 로 미분한 값과 견준다.

사용(노트북):  import sxview;  v = sxview.view('../codegen/specs/baraff_stretch.py', names=['w_u', 'w_v']);  v;  v.check()
"""
import importlib.util
import random

import sympy as sp


class StageVec:
    """보기용 중간 벡터 d_i = Σ_m c_im x_m. 길이는 기호 하나(ell)로 둔다."""

    def __init__(self, index, coeffs, expr):
        self.index, self.coeffs, self.expr = index, coeffs, expr
        self.ell = sp.Symbol('ell%d' % index, positive=True)

    def dot(self, other):
        if not isinstance(other, StageVec):
            raise ValueError('내적은 ctx.stage 로 표시한 벡터끼리만 텐서 꼴로 보일 수 있다')
        if other.index == self.index:
            return self.ell ** 2
        key = tuple(sorted((self.index, other.index)))
        return self.ctx.dots.setdefault(key, sp.Symbol('dot%d_%d' % key, real=True))


class ViewCtx:
    def __init__(self, X):
        self.X = X
        self.stages = []
        self.dots = {}

    def stage(self, v):
        if not (isinstance(v, sp.MatrixBase) and v.shape == (3, 1)):
            raise ValueError('ctx.stage: 3 × 1 벡터만 텐서 꼴로 보일 수 있다')
        coeffs = []
        for xm in self.X:
            c = sp.diff(v[0], xm[0])
            for r in range(3):
                for s in range(3):
                    if sp.simplify(sp.diff(v[r], xm[s]) - (c if r == s else 0)) != 0:
                        raise ValueError('ctx.stage: 좌표의 등방 1 차식(Σ c_m x_m)이 아니다')
            coeffs.append(c)
        rest = v - sum((c * xm for c, xm in zip(coeffs, self.X)), sp.zeros(3, 1))
        if any(sp.simplify(e) != 0 for e in rest):
            raise ValueError('ctx.stage: 상수항이 있다')
        sv = StageVec(len(self.stages), coeffs, v)
        sv.ctx = self
        self.stages.append(sv)
        return sv

    def norm(self, v):
        if not isinstance(v, StageVec):
            raise ValueError('ctx.norm: ctx.stage 로 표시한 벡터의 길이만 텐서 꼴로 보일 수 있다')
        return v.ell

    def unit(self, v):
        raise ValueError('ctx.unit 은 시제품이 다루지 않는다')


class _FlatCtx:   # 성분으로 푸는 기준 (check 용)
    def stage(self, v):
        return v

    def norm(self, v):
        return sp.sqrt((v.T * v)[0])

    def unit(self, v):
        return v / self.norm(v)


def _load(path):
    s = importlib.util.spec_from_file_location('sxspec_view', path)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def _symbols(spec):
    X = [sp.Matrix(sp.symbols('x%d_0:3' % m, real=True)) for m in range(spec.N_NODES)]
    E, names = {}, {}
    for name, ncol in spec.ELEM_DOUBLE_VARS:
        if name == 'pUV' and ncol == 6:
            syms = sp.symbols('u0 v0 u1 v1 u2 v2', real=True)
        elif name == 'pDiag' and ncol == 2:
            syms = sp.symbols('b_p b_m', positive=True)
            names.update({syms[0]: 'b_{+}', syms[1]: 'b_{-}'})
        else:
            syms = sp.symbols('%s_0:%d' % (name, ncol), real=True)
        E[name] = list(syms)
    return X, E, names


class TensorView:
    def __init__(self, spec, path, names):
        self.spec, self.path = spec, path
        X, E, self.sym_names = _symbols(spec)
        self.X, self.E = X, E
        ctx = ViewCtx(X)
        self.C = [sp.sympify(c) for c in spec.condition(ctx, X, E)]
        self.stages = ctx.stages
        self._dots = dict(ctx.dots)
        self.names = list(names or ['d_{%d}' % (i + 1) for i in range(len(self.stages))])
        # 삼각형 명세: 쉬는 틀의 행렬식 Δ 와 넓이 a = |Δ|/2 로 줄여 쓴다
        self.abbrev = []
        if 'pUV' in E:
            u0, v0, u1, v1, u2, v2 = E['pUV']
            det = (u1 - u0) * (v2 - v0) - (u2 - u0) * (v1 - v0)
            self.D, self.a = sp.Symbol('Delta', real=True, nonzero=True), sp.Symbol('a', positive=True)
            self.det = det
            self.C = [c.subs(sp.Abs(det), 2 * self.a) for c in self.C]
            for sv in self.stages:
                sv.coeffs_show = [sp.simplify(c.subs(det, self.D)) for c in sv.coeffs]
            self.abbrev = [(r'\Delta', sp.latex(det)), ('a', r'\tfrac12\,|\Delta|')]
        else:
            for sv in self.stages:
                sv.coeffs_show = [sp.simplify(c) for c in sv.coeffs]
        self.ells = [sv.ell for sv in self.stages]
        # 불변량: 길이 |d_i| 와 내적 d_i · d_j. 조건은 이것들과 스텐실 값만의 함수다
        self.invs = [('norm', i, sv.ell) for i, sv in enumerate(self.stages)] + [('dot', k, q) for k, q in sorted(self._dots.items())]

    # ---- LaTeX 조각 ----
    def _split(self, i):   # 'w_u' → ('w', 'u'),  'd_{+}' → ('d', '{+}')
        base, _, sub = self.names[i].partition('_')
        return base, ('_' + sub if sub else '')

    def _vec(self, i):
        base, sub = self._split(i)
        return r'\mathbf{%s}%s' % (base, sub)

    def _hat(self, i):
        base, sub = self._split(i)
        return r'\hat{\mathbf{%s}}%s' % (base, sub)

    def _c(self, i, m):
        base, sub = self._split(i)
        return r'c^{%s}_{%s}' % (sub[1:] if sub else base, m)

    def _tex(self, e):
        nm = dict(self.sym_names)
        for i, sv in enumerate(self.stages):
            nm[sv.ell] = r'|%s|' % self._vec(i)
        for (i, j), q in self._dots.items():
            nm[q] = r'%s\cdot%s' % (self._vec(i), self._vec(j))
        if hasattr(self, 'D'):
            nm[self.D] = r'\Delta'
        return sp.latex(sp.simplify(e), symbol_names=nm)

    def _u_tex(self, kind, idx, m):
        """불변량의 꼭짓점 m 에 대한 미분(3-벡터)."""
        if kind == 'norm':
            return r'%s\,%s' % (self._c(idx, m), self._hat(idx))
        i, j = idx
        return r'\left(%s\,%s + %s\,%s\right)' % (self._c(i, m), self._vec(j), self._c(j, m), self._vec(i))

    def _b_tex(self, kind, idx):
        """불변량의 꼭짓점 쌍 (m, n) 에 대한 2 차 미분(3 × 3). 길이는 1/|d| 를 뺀 몫."""
        if kind == 'norm':
            return r'%s\,%s\left(\mathbf{I} - %s%s^{\mathsf T}\right)' % (self._c(idx, 'm'), self._c(idx, 'n'), self._hat(idx), self._hat(idx))
        i, j = idx
        return r'\left(%s\,%s + %s\,%s\right)\mathbf{I}' % (self._c(i, 'm'), self._c(j, 'n'), self._c(j, 'm'), self._c(i, 'n'))

    def _scaled(self, e, body):
        e = sp.simplify(e)
        if e == 1:
            return body
        if e == -1:
            return '-' + body
        return r'%s\,%s' % (self._tex(e), body)

    def markdown(self):
        spec, S = self.spec, self.stages
        nC = len(self.C)
        L = ['### %s — 텐서 꼴' % spec.NAME, '', spec.DESCRIPTION, '']
        L += ['**조건 함수** (꼭짓점 %d 개, 성분 %d 개)' % (spec.N_NODES, nC), '']
        for c, e in enumerate(self.C):
            L += ['$$C_{%d} = %s$$' % (c + 1, self._tex(e)), '']
        L += ['**중간 벡터** — 좌표의 1 차식이다. $m$ 은 꼭짓점 번호(0 부터).', '']
        for i, sv in enumerate(S):
            L += ['$$%s = \\sum_m %s\\,\\mathbf{x}_m,\\qquad %s$$' % (
                self._vec(i), self._c(i, 'm'),
                r',\quad '.join('%s = %s' % (self._c(i, m), self._tex(c)) for m, c in enumerate(sv.coeffs_show))), '']
        if self.abbrev:
            L += ['여기서 $' + r',\quad '.join('%s = %s' % ab for ab in self.abbrev) + '$ 이고, 계수는 쉬는 패턴만으로 정해진다'
                  '($\\sum_m c_m = 0$ — 병진 불변).', '']
        L += ['**1 차 미분** — 꼭짓점 $m$ 에 대한 3-벡터', '']
        for c, e in enumerate(self.C):
            terms = [self._scaled(sp.diff(e, q), self._u_tex(kind, idx, 'm')) for kind, idx, q in self.invs
                     if sp.simplify(sp.diff(e, q)) != 0]
            L += ['$$\\frac{\\partial C_{%d}}{\\partial \\mathbf{x}_m} = %s$$' % (c + 1, ' + '.join(terms) or '0'), '']
        L += ['**2 차 미분** — 꼭짓점 쌍 $(m, n)$ 의 3 × 3 블록', '']
        for c, e in enumerate(self.C):
            terms = []
            for kind, idx, q in self.invs:
                g = sp.simplify(sp.diff(e, q))
                if g != 0:
                    terms.append(self._scaled(g / q if kind == 'norm' else g, self._b_tex(kind, idx)))
            for k1, i1, q1 in self.invs:
                for k2, i2, q2 in self.invs:
                    h = sp.simplify(sp.diff(e, q1, q2))
                    if h != 0:
                        terms.append(self._scaled(h, r'\left(%s\right)\left(%s\right)^{\mathsf T}'
                                                  % (self._u_tex(k1, i1, 'm'), self._u_tex(k2, i2, 'n'))))
            L += ['$$\\frac{\\partial^2 C_{%d}}{\\partial \\mathbf{x}_m\\,\\partial \\mathbf{x}_n} = %s$$' % (c + 1, ' + '.join(terms) or '0'), '']
        L += ['**지배방정식에 들어가는 항** — 꼭짓점 블록으로. $\\nabla_m C \\equiv \\partial C/\\partial\\mathbf{x}_m$, '
              '$\\nabla^2_{mn} C \\equiv \\partial^2 C/\\partial\\mathbf{x}_m\\partial\\mathbf{x}_n$, 합은 조건 성분 $c$ 에 대해.', '',
              '$$\\dot C_c = \\sum_n \\nabla_n C_c\\cdot\\mathbf{v}_n$$', '',
              '$$\\mathbf{f}_m = -k\\sum_c C_c\\,\\nabla_m C_c,\\qquad \\mathbf{d}_m = -k_d\\sum_c \\dot C_c\\,\\nabla_m C_c$$', '',
              '$$\\mathbf{K}_{mn} = \\frac{\\partial \\mathbf{f}_m}{\\partial \\mathbf{x}_n} = -k\\sum_c\\Big[\\nabla_m C_c\\,\\nabla_n C_c^{\\mathsf T} '
              '+ \\beta_c\\,\\nabla^2_{mn} C_c\\Big],\\qquad \\mathbf{V}_{mn} = \\frac{\\partial \\mathbf{d}_m}{\\partial \\mathbf{v}_n} '
              '= -k_d\\sum_c \\nabla_m C_c\\,\\nabla_n C_c^{\\mathsf T}$$', '',
              '$$\\mathbf{A} = \\mathbf{M} - h\\,\\mathbf{V} - h^2\\,\\mathbf{K},\\qquad \\mathbf{b} = h\\,(\\mathbf{f} + \\mathbf{d}) + h^2\\,\\mathbf{K}\\,\\mathbf{v}$$', '',
              '| 강성 판 | $\\beta_c$ | 생성 함수 |', '|---|---|---|',
              '| 정확 | $C_c + (k_d/k)\\,\\dot C_c$ | `EvalImplicit` · `EvalEngine` |',
              '| Gauss–Newton | $0$ | `EvalImplicitGN` · `EvalEngineGN` |',
              '| 인장 C > 0 | $\\max(C_c, 0)$ | `EvalImplicitTension` · `EvalEngineTension` |',
              '| 엔진 처방 | $\\max(C_c, 0) + (k_d/k)\\max(\\dot C_c, 0)$ | `EvalEngineTensionRate` |', '']
        return '\n'.join(L)

    def _repr_markdown_(self):
        return self.markdown()

    def check(self, n=50, seed=1):
        """화면의 텐서 식을 숫자로 평가해, condition() 을 성분으로 풀어 미분한 값과 견준다. 돌려주는 것: 최대 상대 오차."""
        import numpy as np
        spec, X, E, S = self.spec, self.X, self.E, self.stages
        xs = [s for xm in X for s in xm]
        ev = [s for name, _ in spec.ELEM_DOUBLE_VARS for s in E[name]]
        flat = [sp.sympify(c) for c in spec.condition(_FlatCtx(), X, E)]
        f_ref = sp.lambdify(xs + ev, [[c] + [sp.diff(c, a) for a in xs] + [sp.diff(c, a, b) for a in xs for b in xs] for c in flat], 'numpy')
        sub = {self.a: sp.Abs(self.det) / 2} if hasattr(self, 'D') else {}
        qs = [q for _, _, q in self.invs]
        f_coef = sp.lambdify(ev, [[c for c in sv.coeffs] for sv in S], 'numpy')
        f_C = sp.lambdify(qs + ev, [c.subs(sub) for c in self.C], 'numpy')
        f_g = sp.lambdify(qs + ev, [[sp.diff(c.subs(sub), q) for q in qs] for c in self.C], 'numpy')
        f_h = sp.lambdify(qs + ev, [[[sp.diff(c.subs(sub), q1, q2) for q2 in qs] for q1 in qs] for c in self.C], 'numpy')
        rnd = random.Random(seed)
        N, I3 = spec.N_NODES, np.eye(3)
        worst = {'C': 0.0, 'grad': 0.0, 'hess': 0.0}
        for _ in range(n):
            xv = [rnd.uniform(-1, 1) for _ in xs]
            evv = []
            for name, ncol in spec.ELEM_DOUBLE_VARS:
                if name == 'pUV':
                    evv += [0.0, 0.0, rnd.uniform(0.6, 1.4), rnd.uniform(-0.3, 0.3), rnd.uniform(-0.3, 0.3), rnd.uniform(0.6, 1.4)]
                else:
                    evv += [rnd.uniform(1.0, 1.6) for _ in range(ncol)]
            ref = f_ref(*(xv + evv))
            coef = np.array(f_coef(*evv), dtype=float)                       # (stages, nodes)
            d = coef @ np.array(xv).reshape(N, 3)                            # (stages, 3)
            ell = np.linalg.norm(d, axis=1)
            hat = d / ell[:, None]
            qv, U, B = [], [], []        # 불변량의 값, ∂q/∂x_m (N × 3), ∂²q/∂x_m∂x_n (N × N × 3 × 3)
            for kind, idx, _q in self.invs:
                if kind == 'norm':
                    i = idx
                    qv.append(ell[i])
                    U.append(coef[i][:, None] * hat[i][None, :])
                    B.append(np.einsum('m,n,ab->mnab', coef[i], coef[i], (I3 - np.outer(hat[i], hat[i])) / ell[i]))
                else:
                    i, j = idx
                    qv.append(d[i] @ d[j])
                    U.append(coef[i][:, None] * d[j][None, :] + coef[j][:, None] * d[i][None, :])
                    B.append(np.einsum('mn,ab->mnab', np.outer(coef[i], coef[j]) + np.outer(coef[j], coef[i]), I3))
            args = qv + evv
            Cv, g, h = f_C(*args), np.array(f_g(*args), dtype=float), np.array(f_h(*args), dtype=float)
            for c in range(len(self.C)):
                r = np.array(ref[c], dtype=float)
                C_ref, g_ref, H_ref = r[0], r[1:1 + 3 * N], r[1 + 3 * N:].reshape(3 * N, 3 * N)
                grad = sum(g[c, a] * U[a] for a in range(len(qs))).reshape(3 * N)
                H4 = sum(g[c, a] * B[a] for a in range(len(qs)))
                for a in range(len(qs)):
                    for b in range(len(qs)):
                        H4 = H4 + h[c, a, b] * np.einsum('ma,nb->mnab', U[a], U[b])
                H = H4.transpose(0, 2, 1, 3).reshape(3 * N, 3 * N)
                worst['C'] = max(worst['C'], abs(Cv[c] - C_ref) / abs(C_ref) if C_ref != 0 else abs(Cv[c]))
                worst['grad'] = max(worst['grad'], np.abs(grad - g_ref).max() / np.abs(g_ref).max())
                worst['hess'] = max(worst['hess'], np.abs(H - H_ref).max() / np.abs(H_ref).max())
        return {'samples': n, 'max_rel_C': worst['C'], 'max_rel_grad': worst['grad'], 'max_rel_hess': worst['hess']}

def view(spec_path, names=None):
    return TensorView(_load(spec_path), spec_path, names)


# ---------------------------------------------------------------------------
# 원자 수준의 전개 — 2007 년 Mathematica 노트북(e_Star_자동화버전.nb)의 표기
#   원자(법선 n₁, 변 ẽ 등)를 꼭짓점의 함수로만 두고 미분을 풀지 않는다. 곱의 법칙으로 펼친 꼴을 ∂n₁/∂p_h 같은 편미분 분수로 보인다.
#   곱은 순서를 지킨다(비가환) — Dot 의 자리가 그대로 남는다. 원자가 의존하지 않는 꼭짓점의 미분은 저절로 0 이 된다.
# ---------------------------------------------------------------------------
def nodes(names):
    """꼭짓점 기호. nodes('h i j k')"""
    return sp.symbols(names, real=True)


def atom(tex, *deps):
    """꼭짓점 deps 에만 의존하는 원자. tex 는 표시할 LaTeX (예: 'n_1', r'\tilde{e}')."""
    f = sp.Function('atom_%d' % len(_ATOMS), commutative=False)
    _ATOMS[f] = tex
    return f(*deps)


_ATOMS = {}


def _leibniz_tex(e):
    from sympy.core.function import AppliedUndef
    names = {}

    def atom_tex(a):
        return _ATOMS.get(a.func, str(a.func))

    def rep_d(d):
        base = atom_tex(d.expr)
        vs = [v for v, n in d.variable_count for _ in range(int(n))]
        den = r'\,'.join(r'\partial p_{%s}' % v for v in vs)
        num = r'\partial %s' % base if len(vs) == 1 else r'\partial^{%d} %s' % (len(vs), base)
        sym = sp.Symbol('D%d' % len(names), commutative=False)
        names[sym] = r'\frac{%s}{%s}' % (num, den)
        return sym

    def rep_a(a):
        sym = sp.Symbol('A%d' % len(names), commutative=False)
        names[sym] = atom_tex(a)
        return sym

    e = e.replace(lambda x: isinstance(x, sp.Derivative), rep_d)
    e = e.replace(lambda x: isinstance(x, AppliedUndef), rep_a)
    return sp.latex(e, symbol_names=names, mul_symbol='dot')


class LeibnizView:
    def __init__(self, C, pts, name='C', upper=False):
        self.C, self.pts, self.name, self.upper = C, list(pts), name, upper

    def markdown(self):
        L = ['$$%s = %s$$' % (self.name, _leibniz_tex(self.C)), '', '**1 차 미분**', '']
        for m in self.pts:
            L += ['$$\\frac{\\partial %s}{\\partial p_{%s}} = %s$$' % (self.name, m, _leibniz_tex(sp.expand(sp.diff(self.C, m)))), '']
        L += ['**2 차 미분** — 꼭짓점 쌍마다(행 = 첫 미분의 꼭짓점)', '']
        for m in self.pts:
            L += ['(%s, ·)' % m, '']
            for n in self.pts:
                if self.upper and self.pts.index(n) < self.pts.index(m):
                    continue
                d = sp.expand(sp.diff(self.C, m, n))
                L += ['$$\\frac{\\partial^2 %s}{\\partial p_{%s}\\,\\partial p_{%s}} = %s$$' % (self.name, m, n, _leibniz_tex(d)), '']
        return '\n'.join(L)

    def _repr_markdown_(self):
        return self.markdown()


def leibniz(C, pts, name='C', upper=False):
    return LeibnizView(C, pts, name, upper)


# ---------------------------------------------------------------------------
# 힌지 원자의 미분을 텐서 꼴로 — 원자 수준 전개(leibniz)에 넣을 값과 표
#   힌지 (h, i, j, k): 삼각형 A = (h, i, j), B = (i, j, k), 공유 변 e = p_j − p_i.
#   N₁ = (p_i − p_h) × (p_j − p_h),  N₂ = (p_j − p_k) × (p_i − p_k),  ẽ = [e]× (ẽ a = e × a)
#   unit=True 면 n = N / |N| 을 원자로 쓴다.
# ---------------------------------------------------------------------------
def _skew(a):
    import numpy as np
    return np.array([[0., -a[2], a[1]], [a[2], 0., -a[0]], [-a[1], a[0], 0.]])


class Hinge:
    def __init__(self, unit=True):
        self.unit = unit
        self.pts = nodes('h i j k')
        h, i, j, k = self.pts
        nm = 'n' if unit else 'N'
        self.n1, self.n2 = atom('%s_1' % nm, h, i, j), atom('%s_2' % nm, i, j, k)
        self.et = atom(r'\tilde{e}', i, j)
        self.eh = atom(r'\tilde{\hat{e}}', i, j)   # 단위 변 ê = e/|e| 의 반대칭 행렬 [ê]×
        # 원자마다 (H, P, Q) 역할: N = (P − H) × (Q − H)
        self._roles = {self.n1.func: (h, i, j), self.n2.func: (k, j, i)}

    def table(self):
        n = 'n' if self.unit else 'N'
        L = ['**원자의 정의** — 힌지 $(h, i, j, k)$', '',
             r'$$\mathbf{N}_1 = (\mathbf{p}_i - \mathbf{p}_h)\times(\mathbf{p}_j - \mathbf{p}_h),\qquad '
             r'\mathbf{N}_2 = (\mathbf{p}_j - \mathbf{p}_k)\times(\mathbf{p}_i - \mathbf{p}_k),\qquad '
             r'\tilde{e} = [\mathbf{p}_j - \mathbf{p}_i]_\times$$', '']
        if self.unit:
            L += [r'$$n_1 = \mathbf{N}_1/|\mathbf{N}_1|,\qquad n_2 = \mathbf{N}_2/|\mathbf{N}_2|$$', '']
        L += [r'$[\mathbf{a}]_\times$ 는 $\mathbf{a}\times$ 의 3 × 3 반대칭 행렬이다.', '',
              '**원자의 1 차 미분** (3 × 3 블록)', '',
              r'$$\frac{\partial \mathbf{N}_1}{\partial p_h} = [\mathbf{p}_j - \mathbf{p}_i]_\times = \tilde{e},\qquad '
              r'\frac{\partial \mathbf{N}_1}{\partial p_i} = [\mathbf{p}_h - \mathbf{p}_j]_\times,\qquad '
              r'\frac{\partial \mathbf{N}_1}{\partial p_j} = [\mathbf{p}_i - \mathbf{p}_h]_\times$$', '',
              r'$$\frac{\partial \mathbf{N}_2}{\partial p_k} = [\mathbf{p}_i - \mathbf{p}_j]_\times = -\tilde{e},\qquad '
              r'\frac{\partial \mathbf{N}_2}{\partial p_j} = [\mathbf{p}_k - \mathbf{p}_i]_\times,\qquad '
              r'\frac{\partial \mathbf{N}_2}{\partial p_i} = [\mathbf{p}_j - \mathbf{p}_k]_\times$$', '']
        if self.unit:
            L += [r'$$\frac{\partial n}{\partial p_m} = \frac{1}{|\mathbf{N}|}\left(\mathbf{I} - n\,n^{\mathsf T}\right)\frac{\partial \mathbf{N}}{\partial p_m}$$', '']
        L += [r'$$\frac{\partial \tilde{e}}{\partial p_j}[\boldsymbol\delta] = [\boldsymbol\delta]_\times,\qquad '
              r'\frac{\partial \tilde{e}}{\partial p_i}[\boldsymbol\delta] = -[\boldsymbol\delta]_\times$$', '',
              r'**원자의 2 차 미분** — 방향 $\mathbf{a}$(첫 꼭짓점) · $\mathbf{b}$(둘째 꼭짓점)에 대한 값. '
              r'$\mathbf{N} = (P - H)\times(Q - H)$ 로 쓰면($\mathbf{N}_1$: $H{=}h, P{=}i, Q{=}j$; $\mathbf{N}_2$: $H{=}k, P{=}j, Q{=}i$)', '',
              r'$$D^2\mathbf{N}[\mathbf{a}, \mathbf{b}] = \sigma\;\mathbf{a}\times\mathbf{b},\qquad '
              r'\sigma_{PQ} = \sigma_{HP} = \sigma_{QH} = +1,\quad \sigma_{QP} = \sigma_{PH} = \sigma_{HQ} = -1,\quad \sigma_{PP} = \sigma_{QQ} = \sigma_{HH} = 0$$', '']
        if self.unit:
            L += [r'$$D^2 n[\mathbf{a}, \mathbf{b}] = \frac{\mathbf{P}\,D^2\mathbf{N}[\mathbf{a}, \mathbf{b}]}{|\mathbf{N}|} - '
                  r'\frac{(n\cdot \mathbf{N}_b)\,\mathbf{P}\mathbf{N}_a + (n\cdot \mathbf{N}_a)\,\mathbf{P}\mathbf{N}_b + (\mathbf{N}_a\cdot\mathbf{P}\mathbf{N}_b)\,n}{|\mathbf{N}|^2},'
                  r'\qquad \mathbf{P} = \mathbf{I} - n\,n^{\mathsf T},\quad \mathbf{N}_a = D\mathbf{N}[\mathbf{a}]$$', '']
        L += [r'$\tilde{e}$ 는 좌표의 1 차식이라 2 차 미분이 0 이다.', '',
              r'**단위 변** $\hat{\mathbf{e}} = \mathbf{e}/|\mathbf{e}|$, $\tilde{\hat{e}} = [\hat{\mathbf{e}}]_\times$ — 단위 법선과 같은 꼴이다'
              r'($\mathbf{e}_a = D\mathbf{e}[\mathbf{a}]$: 꼭짓점 $j$ 이면 $\mathbf{a}$, $i$ 이면 $-\mathbf{a}$).', '',
              r'$$D\hat{\mathbf{e}}[\mathbf{a}] = \frac{\mathbf{P}_e\,\mathbf{e}_a}{|\mathbf{e}|},\qquad '
              r'D^2\hat{\mathbf{e}}[\mathbf{a}, \mathbf{b}] = -\frac{(\hat{\mathbf{e}}\cdot\mathbf{e}_b)\,\mathbf{P}_e\mathbf{e}_a + '
              r'(\hat{\mathbf{e}}\cdot\mathbf{e}_a)\,\mathbf{P}_e\mathbf{e}_b + (\mathbf{e}_a\cdot\mathbf{P}_e\mathbf{e}_b)\,\hat{\mathbf{e}}}{|\mathbf{e}|^2},'
              r'\qquad \mathbf{P}_e = \mathbf{I} - \hat{\mathbf{e}}\hat{\mathbf{e}}^{\mathsf T}$$', '']
        return '\n'.join(L)

    # ---- 숫자 값 (위 표의 식 그대로) ----
    def _N(self, f, P):
        import numpy as np
        H, A, B = (P[x] for x in self._roles[f])
        return np.cross(A - H, B - H)

    def _DN(self, f, P, node, d):
        import numpy as np
        H, A, B = self._roles[f]
        if node == H:
            return _skew(P[B] - P[A]) @ d
        if node == A:
            return -_skew(P[B] - P[H]) @ d
        if node == B:
            return _skew(P[A] - P[H]) @ d
        return np.zeros(3)

    def _D2N(self, f, P, n1, a, n2, b):
        import numpy as np
        H, A, B = self._roles[f]
        sig = {(A, B): 1, (H, A): 1, (B, H): 1, (B, A): -1, (A, H): -1, (H, B): -1}.get((n1, n2), 0)
        return sig * np.cross(a, b)

    def value(self, obj, P, steps):
        """원자나 그 미분의 숫자 값. steps = [(꼭짓점, 방향)] — Derivative 의 변수 순서대로."""
        import numpy as np
        f = obj.func
        if f == self.et.func:
            i, j = self.pts[1], self.pts[2]
            if not steps:
                return _skew(P[j] - P[i])
            if len(steps) == 1:
                node, d = steps[0]
                return _skew(d) if node == j else (-_skew(d) if node == i else np.zeros((3, 3)))
            return np.zeros((3, 3))
        if f == self.eh.func:
            i, j = self.pts[1], self.pts[2]
            e = P[j] - P[i]
            Le = np.linalg.norm(e)
            eh = e / Le
            Pe = np.eye(3) - np.outer(eh, eh)
            de = [(d if node == j else (-d if node == i else np.zeros(3))) for node, d in steps]
            if not steps:
                return _skew(eh)
            if len(steps) == 1:
                return _skew(Pe @ de[0] / Le)
            ea, eb = de
            return _skew(-((eh @ eb) * (Pe @ ea) + (eh @ ea) * (Pe @ eb) + (ea @ (Pe @ eb)) * eh) / Le ** 2)
        N = self._N(f, P)
        L = np.linalg.norm(N)
        if not self.unit:
            if not steps:
                return N
            if len(steps) == 1:
                return self._DN(f, P, *steps[0])
            return self._D2N(f, P, steps[0][0], steps[0][1], steps[1][0], steps[1][1])
        n = N / L
        Pm = np.eye(3) - np.outer(n, n)
        if not steps:
            return n
        Na = self._DN(f, P, *steps[0])
        if len(steps) == 1:
            return Pm @ Na / L
        Nb = self._DN(f, P, *steps[1])
        Nab = self._D2N(f, P, steps[0][0], steps[0][1], steps[1][0], steps[1][1])
        return Pm @ Nab / L - ((n @ Nb) * (Pm @ Na) + (n @ Na) * (Pm @ Nb) + (Na @ (Pm @ Nb)) * n) / L ** 2

    def _nc_eval(self, expr, P, assign):
        """assign = [(꼭짓점, 방향)] 한두 개 — 미분 단계의 방향."""
        total = 0.0
        for term in sp.Add.make_args(sp.expand(expr)):
            c, nc = term.args_cnc()
            ders = [fac for fac in nc if isinstance(fac, sp.Derivative)]
            same = len(assign) == 2 and assign[0][0] == assign[1][0] and len(ders) == 2
            orders = [assign, assign[::-1]] if same else [assign]
            acc = 0.0
            for order in orders:
                left = list(order)
                r = None
                for fac in nc:
                    if isinstance(fac, sp.Derivative):
                        vs = [v for v, cnt in fac.variable_count for _ in range(int(cnt))]
                        steps = []
                        for v in vs:
                            k = next(t for t in left if t[0] == v)
                            left.remove(k)
                            steps.append(k)
                        M = self.value(fac.expr, P, steps)
                    else:
                        M = self.value(fac, P, [])
                    r = M if r is None else r @ M
                acc += r
            total += float(sp.Mul(*c)) * acc / len(orders)
        return total


    def tensor_eval(self, C, P):
        """비가환 곱 식 C 의 값 · 기울기(12) · Hessian(12 × 12) — 곱의 법칙 전개 + 표의 원자 미분."""
        import itertools
        import numpy as np
        pts, E3 = list(self.pts), np.eye(3)
        val = self._nc_eval(C, P, [])
        g = np.array([self._nc_eval(sp.diff(C, m), P, [(m, E3[a])]) for m in pts for a in range(3)])
        H = np.zeros((12, 12))
        for (mi, m), (ni, q) in itertools.product(enumerate(pts), enumerate(pts)):
            d2 = sp.diff(C, m, q)
            for a in range(3):
                for b in range(3):
                    H[3 * mi + a, 3 * ni + b] = self._nc_eval(d2, P, [(m, E3[a]), (q, E3[b])])
        return val, g, H

    def check(self, C, n=20, seed=3):
        """펼친 식(곱의 법칙 + 위 표의 원자 미분)을 숫자로 평가해, 같은 C 를 12 개 좌표 성분으로 풀어 SymPy 로 미분한 값과 견준다."""
        import itertools
        import numpy as np
        from sympy.core.function import AppliedUndef
        pts = list(self.pts)
        xs = {q: sp.Matrix(sp.symbols('x%s_0:3' % q, real=True)) for q in pts}
        flat = [c for q in pts for c in xs[q]]

        def sym_atom(a):
            f = a.func
            if f in (self.et.func, self.eh.func):
                e = xs[pts[2]] - xs[pts[1]]
                if f == self.eh.func:
                    e = e / sp.sqrt((e.T * e)[0])
                return sp.Matrix([[0, -e[2], e[1]], [e[2], 0, -e[0]], [-e[1], e[0], 0]])
            H, A, B = (xs[q] for q in self._roles[f])
            N = (A - H).cross(B - H)
            return N / sp.sqrt((N.T * N)[0]) if self.unit else N

        def sym_eval(term):
            c, nc = term.args_cnc()
            r = None
            for fac in nc:
                M = sym_atom(fac)
                r = M if r is None else (r.T * M if r.shape == (3, 1) and M.shape == (3, 3) else
                                         (r * M if M.shape == (3, 1) and r.shape[1] == 3 else r.T * M))
            return sp.Mul(*c) * (r[0] if r.shape == (1, 1) else r)

        Cs = sp.Add(*[sym_eval(t) for t in sp.Add.make_args(sp.expand(C))])
        f_ref = sp.lambdify(flat, [Cs, [sp.diff(Cs, a) for a in flat], [[sp.diff(Cs, a, b) for b in flat] for a in flat]], 'numpy')

        num_eval = self._nc_eval

        rnd = random.Random(seed)
        worst = {'C': 0.0, 'grad': 0.0, 'hess': 0.0}
        E3 = np.eye(3)
        for _ in range(n):
            P = {pts[0]: np.array([0.3, 1.0, 0.2]), pts[1]: np.array([0.0, 0.0, 0.0]),
                 pts[2]: np.array([1.0, 0.1, 0.0]), pts[3]: np.array([0.4, -0.9, 0.5])}
            for q in pts:
                P[q] = P[q] + np.array([rnd.uniform(-0.2, 0.2) for _ in range(3)])
            ref = f_ref(*[c for q in pts for c in P[q]])
            C_ref, g_ref, H_ref = float(ref[0]), np.array(ref[1], dtype=float), np.array(ref[2], dtype=float)
            worst['C'] = max(worst['C'], abs(num_eval(C, P, []) - C_ref) / abs(C_ref))
            g = np.array([num_eval(sp.diff(C, m), P, [(m, E3[a])]) for m in pts for a in range(3)])
            worst['grad'] = max(worst['grad'], np.abs(g - g_ref).max() / np.abs(g_ref).max())
            H = np.zeros((12, 12))
            for (mi, m), (ni, q) in itertools.product(enumerate(pts), enumerate(pts)):
                d2 = sp.diff(C, m, q)
                for a in range(3):
                    for b in range(3):
                        H[3 * mi + a, 3 * ni + b] = num_eval(d2, P, [(m, E3[a]), (q, E3[b])])
            worst['hess'] = max(worst['hess'], np.abs(H - H_ref).max() / np.abs(H_ref).max())
        return {'samples': n, 'max_rel_C': worst['C'], 'max_rel_grad': worst['grad'], 'max_rel_hess': worst['hess']}


# ---------------------------------------------------------------------------
# 굽힘 명세 (codegen/specs/baraff_bend.py):  sin θ = (n̂A × n̂B) · ê,  cos θ = n̂A · n̂B,  C = atan2(sin θ, cos θ) − θ0
#   원자로 쓰면  s = −n₁ · [ê]× · n₂,  c = n₁ · n₂  (n₁ = n̂A, n₂ = n̂B; (a × b) · ê = −a · (ê × b)).
#   s · c 는 원자 수준으로 펼치고, θ = atan2(s, c) 는 스칼라 연쇄 법칙으로 잇는다. 단위 벡터에서는 s² + c² = 1 이 항등식이다.
# ---------------------------------------------------------------------------
class BendView:
    """굽힘 명세: C = F(θ; 스텐실 값), θ = atan2(s, c). F 는 명세의 condition() 에서 atan2 를 θ 로 바꿔 얻는다.
    명세가 정규화하지 않은 법선 · 변으로 atan2 를 써도(Bridson 명세) atan2 는 양의 배율에 불변이라 같은 θ 다."""

    def __init__(self, spec_path):
        self.spec = _load(spec_path)
        self.hg = Hinge(unit=True)
        self.s = -self.hg.n1 * self.hg.eh * self.hg.n2
        self.c = self.hg.n1 * self.hg.n2
        self.X = [sp.Matrix(sp.symbols('x%d_0:3' % m, real=True)) for m in range(4)]
        self.E, self.ev, self.ev_tex = {}, [], {}
        for name, ncol in self.spec.ELEM_DOUBLE_VARS:
            tex = {'pTheta0': r'\theta_0', 'pBendW': r'\bar{w}'}.get(name)
            syms = [sp.Symbol('%s_%d' % (name, k), **({'positive': True} if name == 'pBendW' else {'real': True})) for k in range(ncol)]
            self.E[name] = syms
            self.ev += syms
            if tex and ncol == 1:
                self.ev_tex[syms[0]] = tex
        self.Cs = sp.sympify(self.spec.condition(_FlatCtx(), self.X, self.E)[0])
        ats = list(self.Cs.atoms(sp.atan2))
        if len(ats) != 1:
            raise ValueError('굽힘 보기: 조건에 atan2 가 정확히 하나 있어야 한다')
        self.theta = sp.Symbol('theta', real=True)
        self.F = self.Cs.xreplace({ats[0]: self.theta})
        if self.F.free_symbols - set(self.ev) - {self.theta}:
            raise ValueError('굽힘 보기: atan2 밖에 좌표가 남아 있다')
        self.F1, self.F2 = sp.simplify(sp.diff(self.F, self.theta)), sp.simplify(sp.diff(self.F, self.theta, 2))

    def _tex(self, e):
        nm = dict(self.ev_tex)
        nm[self.theta] = r'\theta'
        return sp.latex(e, symbol_names=nm)

    def markdown(self, second=True):
        hg = self.hg
        L = ['### %s — 텐서 꼴' % self.spec.NAME, '', self.spec.DESCRIPTION, '',
             r'$$C = F(\theta) = %s,\qquad F^\prime = %s,\qquad F^{\prime\prime} = %s$$' % (self._tex(self.F), self._tex(self.F1), self._tex(self.F2)), '',
             r'$\theta$ 는 두 삼각형의 법선이 이루는 각이다. 단위 원자로 쓰면($n_1 = \hat{\mathbf{n}}_A$, $n_2 = \hat{\mathbf{n}}_B$; '
             r'$(\mathbf{a}\times\mathbf{b})\cdot\hat{\mathbf{e}} = -\mathbf{a}\cdot(\hat{\mathbf{e}}\times\mathbf{b})$)', '',
             r'$$s \equiv \sin\theta = %s,\qquad c \equiv \cos\theta = %s,\qquad \theta = \operatorname{atan2}(s, c)$$'
             % (_leibniz_tex(self.s), _leibniz_tex(self.c)), '',
             '명세가 정규화하지 않은 법선 · 변으로 atan2 를 쓰더라도 atan2 는 양의 배율에 불변이라 같은 $\\theta$ 다(`check()` 가 확인한다).', '',
             r'**연쇄 법칙** — $\nabla_m \equiv \partial/\partial p_m$, $\nabla^2_{mn} \equiv \partial^2/\partial p_m\partial p_n$. '
             r'단위 벡터에서는 $s^2 + c^2 = 1$ 이 항등식이다.', '',
             r'$$\nabla_m \theta = c\,\nabla_m s - s\,\nabla_m c$$', '',
             r'$$\nabla^2_{mn} \theta = c\,\nabla^2_{mn} s - s\,\nabla^2_{mn} c - 2sc\left(\nabla_m s\,\nabla_n s^{\mathsf T} - \nabla_m c\,\nabla_n c^{\mathsf T}\right) '
             r'+ (s^2 - c^2)\left(\nabla_m s\,\nabla_n c^{\mathsf T} + \nabla_m c\,\nabla_n s^{\mathsf T}\right)$$', '',
             r'$$\nabla_m C = F^\prime\,\nabla_m\theta,\qquad \nabla^2_{mn} C = F^\prime\,\nabla^2_{mn}\theta + F^{\prime\prime}\,\nabla_m\theta\,\nabla_n\theta^{\mathsf T}$$', '',
             hg.table(), '']
        if second:
            L += ['#### $s$ 의 전개', '', LeibnizView(self.s, hg.pts, 's', upper=True).markdown(), '',
                  '#### $c$ 의 전개', '', LeibnizView(self.c, hg.pts, 'c', upper=True).markdown(), '',
                  '2 차 미분은 위 삼각(행 ≤ 열)만 적었다. 아래 삼각은 그 전치다.', '']
        return '\n'.join(L)

    def _repr_markdown_(self):
        return self.markdown()

    def check(self, n=3, seed=5):
        """s · c 의 전개(곱의 법칙 + 원자 미분) → θ 의 연쇄 법칙 → F(θ) 로 C · ∇C · ∇²C 를 만들어, 명세의 condition() 을 성분으로 풀어 미분한 값과 견준다."""
        import numpy as np
        hg = self.hg
        flat = [c for xm in self.X for c in xm]
        g_s = [sp.diff(self.Cs, a) for a in flat]
        H_s = [sp.diff(g_s[a], flat[b]) for a in range(12) for b in range(a, 12)]
        f_ref = sp.lambdify(flat + self.ev, [self.Cs] + g_s + H_s, 'numpy', cse=True)
        f_F = sp.lambdify([self.theta] + self.ev, [self.F, self.F1, self.F2], 'numpy')
        rnd = random.Random(seed)
        worst = {'C': 0.0, 'grad': 0.0, 'hess': 0.0, 'r2': 0.0}
        pts = list(hg.pts)
        for _ in range(n):
            base = [np.array([0.3, 1.0, 0.2]), np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.1, 0.0]), np.array([0.4, -0.9, 0.5])]
            P = {q: b + np.array([rnd.uniform(-0.2, 0.2) for _ in range(3)]) for q, b in zip(pts, base)}
            evv = [rnd.uniform(0.5, 2.0) if sym.is_positive else rnd.uniform(-0.3, 0.3) for sym in self.ev]
            ref = np.array(f_ref(*([c for q in pts for c in P[q]] + evv)), dtype=float)
            C_ref, g_ref = ref[0], ref[1:13]
            H_ref = np.zeros((12, 12))
            k = 13
            for a in range(12):
                for b in range(a, 12):
                    H_ref[a, b] = H_ref[b, a] = ref[k]
                    k += 1
            sv, sg, sH = hg.tensor_eval(self.s, P)
            cv, cg, cH = hg.tensor_eval(self.c, P)
            r2 = sv * sv + cv * cv
            ts, tc = cv / r2, -sv / r2
            tss, tcc, tsc = -2 * sv * cv / r2 ** 2, 2 * sv * cv / r2 ** 2, (sv * sv - cv * cv) / r2 ** 2
            tg = ts * sg + tc * cg
            tH = ts * sH + tc * cH + tss * np.outer(sg, sg) + tcc * np.outer(cg, cg) + tsc * (np.outer(sg, cg) + np.outer(cg, sg))
            Fv, F1, F2 = (float(x) for x in f_F(*([np.arctan2(sv, cv)] + evv)))
            g, H = F1 * tg, F1 * tH + F2 * np.outer(tg, tg)
            worst['C'] = max(worst['C'], abs(Fv - C_ref) / abs(C_ref))
            worst['grad'] = max(worst['grad'], np.abs(g - g_ref).max() / np.abs(g_ref).max())
            worst['hess'] = max(worst['hess'], np.abs(H - H_ref).max() / np.abs(H_ref).max())
            worst['r2'] = max(worst['r2'], abs(r2 - 1.0))
        return {'samples': n, 'max_rel_C': worst['C'], 'max_rel_grad': worst['grad'], 'max_rel_hess': worst['hess'],
                'max_abs_s2_c2_minus_1': worst['r2']}

def bend_view(spec_path):
    return BendView(spec_path)
