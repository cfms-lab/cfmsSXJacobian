"""엔진의 나머지 힘 — 손 코드 대 생성기의 유도 (docs/GATE_ENGINE_FORCES.md). 엔진 소스는 고치지 않는다. run_gate.bat 밖.

엔진 cfmsDrape 의 CPU DLL(기준 378d20c, 리포 밖에서 빌드한 것)을 디버그 고리로 찔러 힘 f(X) 와 강성 J = ∂f/∂x 를 재고,
같은 모형을 조건 함수 한 줄로 쓴 명세에서 생성기가 유도한 C · ∇C · ∇²C (codegen/sxgen.py 의 derive) 로 만든 값과 견준다.
  힘       f = −k Σ C ∇C
  J 의 판  GN = −k Σ ∇C∇Cᵀ,  정확판 = GN − k Σ C ∇²C,  C > 0 판 = GN − k Σ max(C, 0) ∇²C
엔진 J 가 엔진 힘의 미분인지는 엔진 힘의 중심차분으로 따로 본다(분모 = 실제 간격).
생성 C++ 을 부르지 않고 생성기의 유도 식을 Python 에서 수치 평가한다 — 명세마다 C++ 기준 도구를 짓지 않으려는 것이다.
그 식 자체는 명세마다 유한차분으로 한 번 확인한다.

사용:  python tests/engine_forces.py <cfmsDrape.dll>      (리포 루트에서, 결과는 results/gate_engine_forces.json)
엔진의 Python 포장 cfmsdrape 는 CFMS_DRAPE_DIR(없으면 기준 워크트리 ..\\cfmsDrape_dev_base378d20c)의 drape 에서 읽는다.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'codegen'))
import sxgen  # noqa: E402

ENGINE_COMMIT = '378d20c'
_BASE_DIR = os.path.join(ROOT, '..', 'cfmsDrape_dev_base' + ENGINE_COMMIT)
ENGINE_DIR = os.path.abspath(os.environ.get('CFMS_DRAPE_DIR') or _BASE_DIR)

SEED = 20261114
FD_STEP = 1e-4                # cm (엔진 힘의 중심차분, L0 과 같다)
LIM_SAME = 1e-10              # 「같다」 (L0 의 인장 문턱과 같다)
LIM_SPEC_FD = 1e-6            # 명세의 유도 식 대 유한차분


def load_drape():
    sys.path.insert(0, os.path.join(ENGINE_DIR, 'drape'))
    from cfmsdrape import Drape
    return Drape


def rel(a, b):
    s = float(np.abs(b).max())
    d = float(np.abs(np.asarray(a) - np.asarray(b)).max())
    return d / s if s > 0 else d


# ---- 생성기의 유도 식 ---------------------------------------------------------------
class SpecEval:
    """명세 하나의 C · ∇C · ∇²C 를 생성기(sxgen.derive, 평탄 정확판)로 유도해 수치 함수로 만든다."""

    def __init__(self, name):
        self.path = os.path.join(ROOT, 'codegen', 'specs', name + '.py')
        self.spec = sxgen.load_spec(self.path)
        xs, ev, C, g, H, _ = sxgen.derive(self.spec, 'exact')
        self.nd, self.nc, self.ne = len(xs), len(C), len(ev)
        self.nh = self.nd * (self.nd + 1) // 2
        self.fn = sp.lambdify(list(xs) + list(ev), [C, g, H], modules='numpy', cse=True)
        self.iu = np.triu_indices(self.nd)

    def __call__(self, x, ev=()):
        C, g, H = self.fn(*np.asarray(x, float).ravel(), *np.atleast_1d(ev))
        C = np.array(C, float)
        g = np.array(g, float).reshape(self.nc, self.nd)
        Hu = np.array(H, float).reshape(self.nc, self.nh)
        Hf = np.zeros((self.nc, self.nd, self.nd))
        for c in range(self.nc):
            Hf[c][self.iu] = Hu[c]
            Hf[c] = Hf[c] + np.triu(Hf[c], 1).T
        return C, g, Hf

    def fd_check(self, x, ev=(), step=1e-6):
        """유도한 ∇C · ∇²C 대 중심차분(분모 = 실제 간격). (∇C 상대 오차, ∇²C 상대 오차)."""
        x = np.asarray(x, float).ravel()
        C, g, H = self(x, ev)
        gfd, Hfd = np.zeros_like(g), np.zeros_like(H)
        for b in range(self.nd):
            xp, xm = x.copy(), x.copy()
            xp[b] += step
            xm[b] -= step
            Cp, gp, _ = self(xp, ev)
            Cm, gm, _ = self(xm, ev)
            hx = xp[b] - xm[b]
            gfd[:, b] = (Cp - Cm) / hx
            Hfd[:, :, b] = (gp - gm) / hx
        return rel(g, gfd), rel(H, Hfd)


def assemble(se, stencils, evs, X, k):
    """전역(조밀) 조립. k 는 스칼라이거나 요소 × 성분 배열. f, Jgg, JHC, JHCp, C 를 낸다."""
    n3 = X.size
    f = np.zeros(n3)
    Jgg, JHC, JHCp = np.zeros((n3, n3)), np.zeros((n3, n3)), np.zeros((n3, n3))
    Cs = []
    kk = np.broadcast_to(np.asarray(k, float).reshape(-1, 1) if np.ndim(k) == 1 else np.asarray(k, float),
                         (len(stencils), se.nc))
    for e, (s, ev) in enumerate(zip(stencils, evs)):
        idx = np.array([3 * v + c for v in s for c in range(3)])
        C, g, H = se(X[list(s)].ravel(), ev)
        ix = np.ix_(idx, idx)
        for c in range(se.nc):
            f[idx] += -kk[e, c] * C[c] * g[c]
            Jgg[ix] += -kk[e, c] * np.outer(g[c], g[c])
            JHC[ix] += -kk[e, c] * C[c] * H[c]
            JHCp[ix] += -kk[e, c] * max(C[c], 0.0) * H[c]
            Cs.append(C[c])
    return f, Jgg, JHC, JHCp, np.array(Cs)


# ---- 엔진 재기 ----------------------------------------------------------------------
def engine_measure(d, X, fd=True):
    """디버그 고리로 f, J(단위 벡터 곱), 엔진 힘의 중심차분 J_fd 를 잰다."""
    n3 = X.size
    out = {'f': d.debug_internal_force(X).ravel()}
    J = np.empty((n3, n3))
    e = np.zeros(n3)
    for k in range(n3):
        e[:] = 0.0
        e[k] = 1.0
        J[:, k] = d.debug_stiffness_product(X, e.reshape(-1, 3)).ravel()
    out['J'] = J
    if fd:
        Xf = X.ravel()
        Jfd = np.empty((n3, n3))
        for k in range(n3):
            xp, xm = Xf.copy(), Xf.copy()
            xp[k] += FD_STEP
            xm[k] -= FD_STEP
            fp = d.debug_internal_force(xp.reshape(-1, 3)).ravel()
            fm = d.debug_internal_force(xm.reshape(-1, 3)).ravel()
            Jfd[:, k] = (fp - fm) / (xp[k] - xm[k])
        out['Jfd'] = Jfd
    return out


def classify(J, variants):
    """엔진 J 가 생성기의 어느 판과 같은가. (가장 가까운 판 이름, 그 상대 오차, 판마다의 상대 오차)."""
    errs = {name: rel(J, Jv) for name, Jv in variants.items()}
    best = min(errs, key=errs.get)
    return best, errs[best], errs


def damping_force(se, stencils, evs, X, vel, kd):
    """d = −kd Σ Ċ ∇C,  Ċ = ∇C·v (전역, 조밀). kd 는 스칼라이거나 요소마다의 배열."""
    d = np.zeros(X.size)
    kk = np.broadcast_to(np.asarray(kd, float), (len(stencils),))
    for e, (s, ev) in enumerate(zip(stencils, evs)):
        idx = np.array([3 * v + c for v in s for c in range(3)])
        _, g, _ = se(X[list(s)].ravel(), ev)
        vl = vel[list(s)].ravel()
        for c in range(se.nc):
            d[idx] += -kk[e] * float(g[c] @ vl) * g[c]
    return d


# ---- 장면 ---------------------------------------------------------------------------
def grid_scene(m=4, s=2.0, x0=0.0):
    n = m + 1
    V = np.array([[x0 + s * i, s * j, 0.0] for j in range(n) for i in range(n)])
    F = []
    for j in range(m):
        for i in range(m):
            a = j * n + i
            F += [[a, a + 1, a + n + 1], [a, a + n + 1, a + n]]
    return V, np.array(F, np.int32)


def deformed(V, rng):
    x, y = V[:, 0], V[:, 1]
    X = np.stack([1.05 * x + 0.08 * y, 0.02 * x + 0.97 * y,
                  0.5 * np.sin(np.pi * x / 8.0) * np.sin(np.pi * y / 8.0)], axis=1)
    return X + rng.uniform(-0.15, 0.15, X.shape)


def new_drape(Drape, dll, V, F, uv, mat, fn=None, panels=None):
    d = Drape(dll)
    d.set_mesh(V, F)
    d.set_mesh_2d(uv)
    if panels:
        d.set_panels(*panels)
    d.set_material(**mat)
    d.set_gravity((0.0, 0.0, 0.0))
    d.set_timestep(0.005)
    if fn:
        fn(d)
    return d


def measure(Drape, dll, V, F, uv, X, mat, fn=None, vel=None, fd=True):
    d = new_drape(Drape, dll, V, F, uv, mat, fn)
    try:
        if vel is not None:
            d.set_velocities(vel)
        return engine_measure(d, X, fd=fd)
    finally:
        d.close()


def compare(name, f_eng, J_eng, Jfd, f, Jgg, JHC, JHCp):
    """한 경우의 기록: 힘, J 의 판, J 대 제 힘의 유한차분."""
    best, err, errs = classify(J_eng, {"GN": Jgg, "exact": Jgg + JHC, "tension": Jgg + JHCp})
    fr, jf = rel(f_eng, f), rel(J_eng, Jfd)
    r = {"f_rel": fr, "J_best": best, "J_best_rel": err, "J_rel": errs, "J_vs_own_fd": jf,
         "spec_exact_vs_engine_fd": rel(Jgg + JHC, Jfd), "f_same": bool(fr <= LIM_SAME), "J_same": bool(err <= LIM_SAME),
         "J_is_derivative": bool(jf <= LIM_SPEC_FD)}
    print("  %-36s 힘 %.1e %s | J %s %.1e %s (GN %.1e · 정확 %.1e · C>0 %.1e) | J 대 제 힘의 차분 %.1e %s"
          % (name, fr, "같다" if r["f_same"] else "다르다", best, err, "같다" if r["J_same"] else "어느 판과도 다르다",
             errs["GN"], errs["exact"], errs["tension"], jf, "미분이다" if r["J_is_derivative"] else "미분이 아니다"))
    return r


def save(out):
    with open(os.path.join(ROOT, "results", "gate_engine_forces.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
        fh.write("\n")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    dll = os.path.abspath(sys.argv[1])
    if not os.path.isfile(dll) or not os.path.isdir(os.path.join(ENGINE_DIR, "drape")):
        print("DLL · 엔진 리포 중 없는 것이 있음: %s, %s" % (dll, ENGINE_DIR))
        return 2
    import hashlib
    import itertools
    import warnings
    warnings.simplefilter("ignore")   # Piecewise 의 안 쓰는 갈래에서 나는 sqrt 경고
    Drape = load_drape()
    rng = np.random.default_rng(SEED)
    V, F = grid_scene()
    uv = V[:, :2].copy()
    X = deformed(V, rng)
    tri = [tuple(int(v) for v in t) for t in F]
    uvs = [uv[t].ravel() for t in F]
    E, SE, VIS = 1e5, 4e4, 1e5
    B = float(np.float32(1.41421356237))   # 엔진의 b_su = b_sv (float 리터럴)
    base = dict(mass=0.1, E=E, V=VIS, SE=0.0, BE_U=0.0)
    out = {"dll_sha256": hashlib.sha256(open(dll, "rb").read()).hexdigest(), "engine_commit": ENGINE_COMMIT,
           "nodes": int(len(V)), "triangles": int(len(F)), "lim_same": LIM_SAME, "lim_fd": LIM_SPEC_FD}
    d0 = new_drape(Drape, dll, V, F, uv, base)
    out["build_info"] = d0.build_info
    out["threads"] = int(d0.openmp_threads)
    d0.close()
    print("엔진의 나머지 힘 — 노드 %d, 삼각형 %d, 엔진 빌드 %s" % (len(V), len(F), out["build_info"]))

    # ---- 전제 H0 · H1 --------------------------------------------------------
    specs = {n: SpecEval(n) for n in ("baraff_stretch", "seam_spring", "point_tether", "stretch_nonlinear",
                                       "provot_shear_nonlinear", "shear_polynomial")}
    print("전제")
    m_lin = measure(Drape, dll, V, F, uv, X, base)
    f, Jgg, JHC, JHCp, _ = assemble(specs["baraff_stretch"], tri, uvs, X, E)
    out["H0"] = compare("H0 기본 인장 (L0 을 다시)", m_lin["f"], m_lin["J"], m_lin["Jfd"], f, Jgg, JHC, JHCp)
    h0 = out["H0"]["f_same"] and out["H0"]["J_same"] and out["H0"]["J_best"] == "tension"
    x3 = X[list(tri[3])].ravel()
    fd_in = {"seam_spring": (X[[0, 4]].ravel(), ()),
             "point_tether": (X[3], X[3] + 0.4),
             "stretch_nonlinear": (x3, np.concatenate([uvs[3], [0.3, 0.02, 2.0]])),
             "provot_shear_nonlinear": (x3, np.concatenate([uvs[3], [B, B], [0.02, 2.0]])),
             "shear_polynomial": (x3, np.concatenate([uvs[3], [4e4, 3e5, 2e6]]))}
    out["H1"] = {}
    h1 = True
    for n, (xx, ev) in fd_in.items():
        g_err, h_err = specs[n].fd_check(xx, ev)
        ok = g_err <= LIM_SPEC_FD and h_err <= LIM_SPEC_FD
        h1 = h1 and ok
        out["H1"][n] = {"grad_rel": g_err, "hess_rel": h_err, "ok": bool(ok)}
        print("  H1 %-24s ∇C %.1e · ∇²C %.1e  %s" % (n, g_err, h_err, "PASS" if ok else "FAIL"))
    out["premise"] = bool(h0 and h1)
    if not out["premise"]:
        print("전제가 깨져 나머지를 적지 않는다")
        save(out)
        return 1

    # ---- A 봉제 -----------------------------------------------------------------
    print("A 봉제 스프링 (봉제를 켠 것 − 끈 것)")
    n = 5
    pairs = np.array([[j * n, j * n + n - 1] for j in range(n)], np.int32)
    st = [tuple(int(v) for v in p) for p in pairs]
    no_ev = [()] * len(st)
    w = np.array([0.5, 1.0, 2.0, 0.0, 1.5])
    vel = rng.uniform(-5.0, 5.0, X.shape)

    def sew(d):
        d.set_sewing(pairs)

    def sew_w(d):
        d.set_sewing(pairs)
        d.set_sewing_weights(w)

    a = measure(Drape, dll, V, F, uv, X, base, sew)
    fs, Js, Jfds = a["f"] - m_lin["f"], a["J"] - m_lin["J"], a["Jfd"] - m_lin["Jfd"]
    f, Jgg, JHC, JHCp, _ = assemble(specs["seam_spring"], st, no_ev, X, 10.0 * E)
    out["A1"] = compare("A1 봉제 (k = 10 E)", fs, Js, Jfds, f, Jgg, JHC, JHCp)
    aw = measure(Drape, dll, V, F, uv, X, base, sew_w)
    f, Jgg, JHC, JHCp, _ = assemble(specs["seam_spring"], st, no_ev, X, 10.0 * E * w)
    out["A2"] = compare("A2 봉제 가중 (k = 10 E w)", aw["f"] - m_lin["f"], aw["J"] - m_lin["J"], aw["Jfd"] - m_lin["Jfd"],
                        f, Jgg, JHC, JHCp)
    av = measure(Drape, dll, V, F, uv, X, base, sew, vel=vel, fd=False)
    bv = measure(Drape, dll, V, F, uv, X, base, None, vel=vel, fd=False)
    d_eng = (av["f"] - bv["f"]) - fs
    d_spec = damping_force(specs["seam_spring"], st, no_ev, X, vel, 2.0 * 10.0 * VIS)
    out["A3"] = {"d_rel": rel(d_eng, d_spec), "same": bool(rel(d_eng, d_spec) <= LIM_SAME),
                 "velocity_changes_membrane_force": float(np.abs(bv["f"] - m_lin["f"]).max())}
    print("  %-36s 속도로 생긴 봉제 힘 %.1e %s" % ("A3 봉제 감쇠 (kd = 2 · 10 · V)", out["A3"]["d_rel"],
                                            "같다" if out["A3"]["same"] else "다르다"))

    # ---- B · C 비선형 막 법칙 ---------------------------------------------------
    print("B 비선형 막 법칙 — 인장")
    out["B"] = {}
    for cb, ec, kh in [(0.3, 0.0, 0.0), (1.0, 0.02, 2e5), (0.3, 0.02, 2e5)]:
        def law(d, cb=cb, ec=ec, kh=kh):
            d.set_membrane_nonlinear(cb=cb, eps_crit=ec, k_hard=kh)
        m = measure(Drape, dll, V, F, uv, X, base, law)
        evs = [np.concatenate([u, [cb, ec, kh / E]]) for u in uvs]
        f, Jgg, JHC, JHCp, Cs = assemble(specs["stretch_nonlinear"], tri, evs, X, E)
        key = "cb=%g,eps_c=%g,k_hard=%g" % (cb, ec, kh)
        out["B"][key] = compare("B " + key, m["f"], m["J"], m["Jfd"], f, Jgg, JHC, JHCp)
        out["B"][key].update(C_pos=int((Cs > 0).sum()), C_neg=int((Cs < 0).sum()))
    print("C 비선형 막 법칙 — Provot 전단 (전단을 켠 것 − 끈 것)")
    out["C"] = {}
    for cb, ec, kh in [(1.0, 0.02, 2e5), (0.3, 0.05, 2e5)]:
        def law(d, cb=cb, ec=ec, kh=kh):
            d.set_membrane_nonlinear(cb=cb, eps_crit=ec, k_hard=kh)
        a = measure(Drape, dll, V, F, uv, X, dict(base, SE=SE), law)
        b = measure(Drape, dll, V, F, uv, X, base, law)
        evs = [np.concatenate([u, [B, B], [ec, kh / E]]) for u in uvs]
        f, Jgg, JHC, JHCp, Cs = assemble(specs["provot_shear_nonlinear"], tri, evs, X, SE)
        key = "cb=%g,eps_c=%g,k_hard=%g" % (cb, ec, kh)
        out["C"][key] = compare("C " + key, a["f"] - b["f"], a["J"] - b["J"], a["Jfd"] - b["Jfd"], f, Jgg, JHC, JHCp)
        out["C"][key].update(C_pos=int((Cs > 0).sum()), C_neg=int((Cs < 0).sum()))

    # ---- D 전단 접선 다항식 -----------------------------------------------------
    print("D 전단 접선 다항식")
    coef = [4e4, 3e5, 2e6]
    off = dict(mass=0.1, E=0.0, V=0.0, SE=0.0, SV=0.0, BE_U=0.0, BV=0.0)

    def poly(d):
        d.set_shear_tangent_polynomial(coef, 1.2)

    m = measure(Drape, dll, V, F, uv, X, off, poly)
    evs = [np.concatenate([u, coef]) for u in uvs]
    f, Jgg, JHC, JHCp, _ = assemble(specs["shear_polynomial"], tri, evs, X, 1.0)
    out["D1"] = compare("D1 격자", m["f"], m["J"], m["Jfd"], f, Jgg, JHC, JHCp)
    out["D1"]["J_asymmetry"] = float(np.abs(m["J"] - m["J"].T).max() / np.abs(m["J"]).max())
    V1 = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.5, 1.7, 0.0]])
    X1 = np.array([[0.1, -0.05, 0.02], [2.25, 0.2, 0.3], [0.95, 1.6, -0.25]])
    out["D2"] = {}
    for perm in itertools.permutations(range(3)):
        F1 = np.array([list(perm)], np.int32)
        m = measure(Drape, dll, V1, F1, V1[:, :2].copy(), X1, off, poly)
        ev = np.concatenate([V1[list(perm), :2].ravel(), coef])
        f, Jgg, JHC, JHCp, _ = assemble(specs["shear_polynomial"], [perm], [ev], X1, 1.0)
        H = Jgg + JHC
        key = "(%d, %d, %d)" % perm
        r = compare("D2 삼각형 " + key, m["f"], m["J"], m["Jfd"], f, Jgg, JHC, JHCp)
        blocks = {}
        for p in range(3):
            for q in range(p + 1, 3):
                Jb, Hb = m["J"][3 * p:3 * p + 3, 3 * q:3 * q + 3], H[3 * p:3 * p + 3, 3 * q:3 * q + 3]
                blocks["%d-%d" % (p, q)] = {"vs_exact": rel(Jb, Hb), "vs_exact_transposed": rel(Jb, Hb.T),
                                            "local_order_reversed": bool(perm.index(p) > perm.index(q))}
        r["upper_blocks"] = blocks
        out["D2"][key] = r
    tr = [b for r in out["D2"].values() for b in r["upper_blocks"].values()]
    agree = sum(1 for b in tr if (b["vs_exact_transposed"] <= LIM_SAME) == b["local_order_reversed"]
                and (b["vs_exact"] <= LIM_SAME) == (not b["local_order_reversed"]))
    out["D2_block_rule"] = {"blocks": len(tr), "transposed_exactly_when_local_order_reversed": agree}
    print("  D2 노드쌍 블록 %d 개 가운데 %d 개가 「요소 안의 순서가 노드 번호와 거꾸로일 때에만 전치돼 있다」에 맞는다" % (len(tr), agree))

    # ---- E 테더 (한 단계) ---------------------------------------------------------
    print("E 꼭짓점 테더 (한 단계 적분, 내부 힘 · 감쇠 · 공기 저항 0)")
    Vt, Ft = grid_scene(2)
    uvt = Vt[:, :2].copy()
    Xt = Vt + rng.uniform(-0.2, 0.2, Vt.shape)
    V0 = rng.uniform(-3.0, 3.0, Vt.shape)
    idx = np.array([1, 4, 7], np.int32)
    tgt = Xt[idx] + rng.uniform(-0.5, 0.5, (3, 3))
    h, mass = 0.005, 0.1
    quiet = dict(mass=mass, E=0.0, V=0.0, SE=0.0, SV=0.0, BE_U=0.0, BV=0.0, air_damping=0.0, internal_damping=0.0)
    out["E"] = {}
    se = specs["point_tether"]
    for k, c in [(1e3, 0.0), (1e3, 2.0), (1e5, 50.0)]:
        d = new_drape(Drape, dll, Vt, Ft, uvt, quiet)
        try:
            d.set_timestep(h)
            d.set_air_drag(0.0)
            d.set_damping(0.0)
            d.step(1)                      # 설정 뒤 첫 단계는 준비 단계다 — 쉬는 모양에서 먼저 한 번 밟는다(L3 과 같다)
            d.set_positions(Xt)
            d.set_velocities(V0)
            d.set_vertex_tethers(idx, tgt, stiffness=k, damping=c)
            d.step(1)
            x1, v1 = d.positions().copy(), d.velocities().copy()
        finally:
            d.close()
        v_pred = np.zeros((len(idx), 3))
        for i, node in enumerate(idx):     # 식 (6): A = M − h ∂f/∂v − h² ∂f/∂x, b = h(f + d) + h² ∂f/∂x v
            C, g, _ = se(Xt[node], tgt[i])
            gg = sum(np.outer(g[cc], g[cc]) for cc in range(se.nc))
            fd_ = -k * (C @ g) - c * (gg @ V0[node])
            A = mass * np.eye(3) + h * c * gg + h * h * k * gg
            b = h * fd_ - h * h * k * (gg @ V0[node])
            v_pred[i] = V0[node] + np.linalg.solve(A, b)
        free = np.setdiff1d(np.arange(len(Vt)), idx)
        key = "k=%g,c=%g" % (k, c)
        vr = rel(v1[idx], v_pred)
        out["E"][key] = {"v_rel": vr, "x_abs": float(np.abs(x1[idx] - (Xt[idx] + h * v_pred)).max()),
                         "free_nodes_dv": float(np.abs(v1[free] - V0[free]).max()), "same": bool(vr <= 1e-9)}
        print("  %-36s 한 단계 뒤 속도 %.1e %s (테더 없는 노드의 속도 변화 %.1e)"
              % ("E " + key, vr, "같다" if vr <= 1e-9 else "다르다", out["E"][key]["free_nodes_dv"]))

    # ---- F 패널별 재질 (보고) ------------------------------------------------------
    print("F 패널별 재질이 막 힘에 먹는가 (보고)")
    Va, Fa = grid_scene(2)
    Vb, Fb = grid_scene(2, x0=10.0)
    Vp = np.vstack([Va, Vb])
    Fp = np.vstack([Fa, Fb + len(Va)]).astype(np.int32)
    uvp = Vp[:, :2].copy()
    Xp = Vp * np.array([1.06, 0.95, 1.0]) + rng.uniform(-0.1, 0.1, Vp.shape)
    trip = [tuple(int(v) for v in t) for t in Fp]
    uvsp = [uvp[t].ravel() for t in Fp]
    E1 = 4e5
    cand = {"uniform_E0": E, "uniform_E1": E1,
            "per_panel": np.array([[E, E]] * len(Fa) + [[E1, E1]] * len(Fb)),
            "per_panel_swapped": np.array([[E1, E1]] * len(Fa) + [[E, E]] * len(Fb))}
    fc = {n_: assemble(specs["baraff_stretch"], trip, uvsp, Xp, k_)[0] for n_, k_ in cand.items()}
    out["F"] = {}
    for label, which in (("패널 1 만 E 를 4 배", 1), ("패널 0 만 E 를 4 배", 0)):
        def panel(dd, which=which):
            dd.set_material_panel(which, E=E1)
        d = new_drape(Drape, dll, Vp, Fp, uvp, base, panel, panels=([len(Va), len(Vb)], [len(Fa), len(Fb)]))
        try:
            fe = d.debug_internal_force(Xp).ravel()
        finally:
            d.close()
        rr = {n_: rel(fe, f_) for n_, f_ in fc.items()}
        out["F"]["panel_%d_x4" % which] = rr
        print("  %-22s 엔진 힘 대 — 균일 E %.1e · 균일 4E %.1e · 패널별 %.1e · 패널별(뒤바꿈) %.1e"
              % (label, rr["uniform_E0"], rr["uniform_E1"], rr["per_panel"], rr["per_panel_swapped"]))

    save(out)
    print("\n엔진의 나머지 힘: 전제 PASS → results/gate_engine_forces.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
