"""패널별 재질이 막 힘에 먹게 한 수정 — 전 · 뒤 견주기 (docs/GATE_ENGINE_FORCES.md 「수정 2」). run_gate.bat 밖.

엔진 브랜치 sxjacobian-link 의 막 힘 커널(손 인장 · 손 Provot 전단 · 다항식 전단의 점성 · 생성 인장 · 생성 전단)이 재질을
노드 0 이 아니라 요소의 첫 노드에서 읽는지를 본다. 팔 넷을 저마다 다른 프로세스로 돌린다:
  before · after = 고치기 전(브랜치 3a35c64) · 고친 DLL, 기본 환경(생성 커널)
  before_kill · after_kill = 같은 DLL + CFMS_CPU_SXJ_BEND=off · CFMS_CPU_SXJ_MEMBRANE=off (손 커널)
기준 M0 – M7 은 판정문에 실행 전에 고정했다. 생성기 쪽 값은 tests/engine_forces.py 와 같은 길(sxgen.derive 를 Python 에서 수치 평가)이다.

사용:  python tests/engine_panel_material_fix.py <고치기 전 cfmsDrape.dll 의 사본> <고친 cfmsDrape.dll 의 사본>
엔진의 Python 포장 · 시험 스크립트는 CFMS_DRAPE_LINK_DIR(없으면 ..\\cfmsDrape_dev_sxjacobian-link)에서 읽는다(읽기만).
M7 은 엔진의 가상 환경(CFMS_DRAPE_MAIN_DIR, 없으면 ..\\cfmsDrape_dev 의 drape\\.venv)의 python 으로 엔진 시험을 돌린다.
"""
import sys
sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import engine_forces as ef  # noqa: E402
import engine_shear_poly_fix as fx  # noqa: E402

LINK_DIR, ENGINE_PY = fx.LINK_DIR, fx.ENGINE_PY
BRANCH_COMMIT = '3a35c64'
FIXED_FILES = {'drape/src/sde_core/SDECore.cpp', 'drape/src/sde_core/SDECore_SXJMembrane.cpp'}

E0, SE0, V0, SV0 = 1e5, 4e4, 1e3, 5e2          # 바탕 재질
KMUL, KDMUL = 4.0, 3.0                          # 재지정: 강성 4 배, 점성 3 배
B = float(np.float32(1.41421356237))            # 엔진의 b_su = b_sv (float 리터럴)
LAW = (0.3, 0.02, 2e5)                          # 비선형 막 법칙 cb · ε_c · k_hard
LIM_STEP = fx.LIM_STEP
ARMS = ('before', 'before_kill', 'after', 'after_kill')


def sha(a):
    return fx.sha(a)


def membrane_terms(se, stencils, evs, X, vel, k, kd):
    """f = −k C∇C − kd Ċ∇C,  J = −k ∇C∇Cᵀ − (k max(C, 0) + kd max(Ċ, 0)) ∇²C,  V = −kd ∇C∇Cᵀ.  k · kd 는 요소마다."""
    n3 = X.size
    f, J, Vm = np.zeros(n3), np.zeros((n3, n3)), np.zeros((n3, n3))
    kk = np.broadcast_to(np.asarray(k, float), (len(stencils),))
    kdd = np.broadcast_to(np.asarray(kd, float), (len(stencils),))
    for e, (s, ev) in enumerate(zip(stencils, evs)):
        idx = np.array([3 * v + c for v in s for c in range(3)])
        C, g, H = se(X[list(s)].ravel(), ev)
        vl = vel[list(s)].ravel()
        ix = np.ix_(idx, idx)
        for c in range(se.nc):
            cd = float(g[c] @ vl)
            f[idx] += -(kk[e] * C[c] + kdd[e] * cd) * g[c]
            J[ix] += -kk[e] * np.outer(g[c], g[c]) - (kk[e] * max(C[c], 0.0) + kdd[e] * max(cd, 0.0)) * H[c]
            Vm[ix] += -kdd[e] * np.outer(g[c], g[c])
    return f, J, Vm


def rate_terms_e(se, stencils, evs, X, vel, kd):
    """다항식 전단의 감쇠 몫: d = −kd Ċ∇C,  J_rate = −kd Ċ∇²C (max 없음).  kd 는 요소마다."""
    n3 = X.size
    d, Jr = np.zeros(n3), np.zeros((n3, n3))
    kdd = np.broadcast_to(np.asarray(kd, float), (len(stencils),))
    for e, (s, ev) in enumerate(zip(stencils, evs)):
        idx = np.array([3 * v + c for v in s for c in range(3)])
        _, g, H = se(X[list(s)].ravel(), ev)
        vl = vel[list(s)].ravel()
        for c in range(se.nc):
            cd = float(g[c] @ vl)
            d[idx] += -kdd[e] * cd * g[c]
            Jr[np.ix_(idx, idx)] += -kdd[e] * cd * H[c]
    return d, Jr


class Scene:
    """삼각형 메시 하나와 그 패널 나눔. 요소의 패널 = 첫 노드의 패널."""

    def __init__(self, V, F, node_counts, elem_counts, X, vel):
        self.V, self.F, self.X, self.vel = V, F, X, vel
        self.uv = V[:, :2].copy()
        self.panels = (list(node_counts), list(elem_counts))
        self.tri = [tuple(int(v) for v in t) for t in F]
        self.uvs = [self.uv[list(t)].ravel() for t in self.tri]
        edge = np.cumsum(node_counts)
        self.panel_of_first_node = np.array([int(np.searchsorted(edge, t[0], side='right')) for t in self.tri])
        self.panel_of_elem_range = np.repeat(np.arange(len(elem_counts)), elem_counts)

    def per_elem(self, base, over, which='first_node'):
        """재지정 {패널: 값} 을 요소마다의 배열로."""
        pe = self.panel_of_first_node if which == 'first_node' else self.panel_of_elem_range
        return np.array([over.get(int(p), base) for p in pe], float)


def spec_membrane(specs, sc, vel, kE, kdE, kS, kdS):
    f1, J1, V1 = membrane_terms(specs['stretch'], sc.tri, sc.uvs, sc.X, vel, kE, kdE)
    f2, J2, V2 = membrane_terms(specs['shear'], sc.tri, [np.concatenate([u, [B, B]]) for u in sc.uvs], sc.X, vel, kS, kdS)
    return f1 + f2, J1 + J2, V1 + V2


def hook(Drape, dll, sc, mat, fn=None, vel=None):
    d = ef.new_drape(Drape, dll, sc.V, sc.F, sc.uv, mat, fn, panels=sc.panels)
    try:
        if vel is not None:
            d.set_velocities(vel)
        return ef.engine_measure(d, sc.X, fd=False)
    finally:
        d.close()


def one_step(Drape, dll, sc, mat, fn=None):
    d = ef.new_drape(Drape, dll, sc.V, sc.F, sc.uv, mat, fn, panels=sc.panels)
    try:
        d.set_timestep(fx.STEP_H)
        d.set_air_drag(0.0)
        d.set_damping(0.0)
        d.set_dynamic_pcg(True, fx.PCG_REL, fx.PCG_CAP)
        d.set_dynamic_pcg_absolute_tolerance(fx.PCG_ABS)
        d.step(1)                          # 설정 뒤 첫 단계는 준비 단계다
        d.set_positions(sc.X)
        d.set_velocities(sc.vel)
        d.step(1)
        return d.velocities().copy(), d.get_pcg_solve_history()
    finally:
        d.close()


def step_case(Drape, dll, specs, sc, mat, fn, kE, kdE, kS, kdS):
    v1, hist = one_step(Drape, dll, sc, mat, fn)
    f, J, Vm = spec_membrane(specs, sc, sc.vel, kE, kdE, kS, kdS)
    h, m = fx.STEP_H, fx.STEP_MASS
    A = m * np.eye(sc.X.size) - h * Vm - h * h * J
    b = h * f + h * h * (J @ sc.vel.ravel())
    v_pred = sc.vel.ravel() + np.linalg.solve(A, b)
    eig = np.linalg.eigvalsh(0.5 * (A + A.T))
    reasons = [r['exit_reason'] for r in hist]
    return {'v_rel': ef.rel(v1.ravel(), v_pred), 'v_sha': sha(v1), 'pcg_exit': reasons,
            'pcg_iterations': [int(r['iterations']) for r in hist],
            'pcg_ok': bool(len(reasons) > 0 and all(r == 'relative' for r in reasons)),
            'A_min_eig': float(eig[0]), 'A_max_eig': float(eig[-1]),
            'max_hV': float(np.abs(h * Vm).max()), 'max_h2J': float(np.abs(h * h * J).max())}


def child(dll, out_path):
    import warnings
    warnings.simplefilter('ignore')
    ef.ENGINE_DIR = LINK_DIR
    Drape = ef.load_drape()
    rng = np.random.default_rng(ef.SEED)
    Va, Fa = ef.grid_scene(2)
    Vb, Fb = ef.grid_scene(2, x0=10.0)
    Vp = np.vstack([Va, Vb])
    Fp = np.vstack([Fa, Fb + len(Va)]).astype(np.int32)
    Xp = Vp * np.array([1.06, 0.95, 1.0]) + rng.uniform(-0.1, 0.1, Vp.shape)
    two = Scene(Vp, Fp, [len(Va), len(Vb)], [len(Fa), len(Fb)], Xp, rng.uniform(-5.0, 5.0, Vp.shape))
    Vg, Fg = ef.grid_scene()
    Xg = ef.deformed(Vg, rng)
    straddle = Scene(Vg, Fg, [13, 12], [16, 16], Xg, np.zeros_like(Xg))     # 한 격자를 노드 범위 둘로 — 요소가 걸친다
    V17, F17 = ef.grid_scene(16)
    X17 = ef.deformed(V17, rng)
    zero = np.zeros_like(Xp)
    specs = {'stretch': ef.SpecEval('baraff_stretch'), 'shear': ef.SpecEval('provot_shear'),
             'stretch_nl': ef.SpecEval('stretch_nonlinear'), 'shear_nl': ef.SpecEval('provot_shear_nonlinear'),
             'poly': ef.SpecEval('shear_polynomial'), 'q': fx.q_spec_eval()}
    base = dict(mass=0.1, E=E0, V=V0, SE=SE0, SV=SV0, BE_U=0.0, BV=0.0)

    out = {'dll_sha256': hashlib.sha256(open(dll, 'rb').read()).hexdigest()}
    d0 = ef.new_drape(Drape, dll, Vp, Fp, two.uv, base, panels=two.panels)
    out['build_info'] = d0.build_info
    out['threads'] = int(d0.openmp_threads)
    d0.close()

    def panel(which, **kw):
        def fn(d):
            d.set_material_panel(which, **kw)
        return fn

    def errs(m, cands):
        return {n: {'f_rel': ef.rel(m['f'], f), 'J_rel': ef.rel(m['J'], J)} for n, (f, J) in cands.items()}

    # ---- M1 · M2: 속도 0 에서 힘 · J ----
    out['static'] = {}
    for name, which, over in (('p1_E', 1, {'E': KMUL * E0}), ('p1_SE', 1, {'SE': KMUL * SE0}),
                              ('p1_E_SE', 1, {'E': KMUL * E0, 'SE': KMUL * SE0}),
                              ('p0_E_SE', 0, {'E': KMUL * E0, 'SE': KMUL * SE0})):
        m = hook(Drape, dll, two, base, panel(which, **over))
        e1, s1 = over.get('E', E0), over.get('SE', SE0)
        cands = {'per_panel': spec_membrane(specs, two, zero, two.per_elem(E0, {which: e1}), 0.0, two.per_elem(SE0, {which: s1}), 0.0)[:2],
                 'uniform_base': spec_membrane(specs, two, zero, E0, 0.0, SE0, 0.0)[:2],
                 'uniform_override': spec_membrane(specs, two, zero, e1, 0.0, s1, 0.0)[:2]}
        out['static'][name] = errs(m, cands)
        r = out['static'][name]
        print('  %-10s 힘 — 패널별 %.1e · 균일(바탕) %.1e · 균일(재지정) %.1e | J — %.1e · %.1e · %.1e'
              % (name, r['per_panel']['f_rel'], r['uniform_base']['f_rel'], r['uniform_override']['f_rel'],
                 r['per_panel']['J_rel'], r['uniform_base']['J_rel'], r['uniform_override']['J_rel']))

    # ---- M3: 속도를 준 상태 — 감쇠 계수도 패널별 ----
    over = {'E': KMUL * E0, 'SE': KMUL * SE0, 'V': KDMUL * V0, 'SV': KDMUL * SV0}
    kE, kS = two.per_elem(E0, {1: over['E']}), two.per_elem(SE0, {1: over['SE']})
    kdE, kdS = two.per_elem(V0, {1: over['V']}), two.per_elem(SV0, {1: over['SV']})
    m = hook(Drape, dll, two, base, panel(1, **over), vel=two.vel)
    f_el = spec_membrane(specs, two, zero, kE, 0.0, kS, 0.0)[0]
    f_pp = spec_membrane(specs, two, two.vel, kE, kdE, kS, kdS)[0]
    out['velocity'] = errs(m, {'per_panel': spec_membrane(specs, two, two.vel, kE, kdE, kS, kdS)[:2],
                               'uniform_base': spec_membrane(specs, two, two.vel, E0, V0, SE0, SV0)[:2]})
    out['velocity']['damping_share_of_f'] = float(np.abs(f_pp - f_el).max() / np.abs(f_pp).max())
    print('  %-10s 힘 — 패널별 %.1e · 균일(바탕) %.1e | J — %.1e · %.1e' % ('속도 · 점성', out['velocity']['per_panel']['f_rel'],
          out['velocity']['uniform_base']['f_rel'], out['velocity']['per_panel']['J_rel'], out['velocity']['uniform_base']['J_rel']))

    # ---- M4: 한 단계 (V 블록) — 균일 재질(전제 · M6) 과 패널별 ----
    quiet = dict(base, mass=fx.STEP_MASS, air_damping=0.0, internal_damping=0.0)
    out['step'] = {'uniform': step_case(Drape, dll, specs, two, quiet, None, E0, V0, SE0, SV0),
                   'per_panel': step_case(Drape, dll, specs, two, quiet, panel(1, **over), kE, kdE, kS, kdS)}
    for key, r in out['step'].items():
        print('  %-10s 한 단계 뒤 속도 %.1e (PCG %s, 반복 %s)' % ('step ' + key, r['v_rel'], r['pcg_exit'], r['pcg_iterations']))

    # ---- M5ⓐ: 비선형 막 법칙 + 패널별 E · SE ----
    cb, ec, kh = LAW

    def nl(d):
        d.set_membrane_nonlinear(cb=cb, eps_crit=ec, k_hard=kh)
        d.set_material_panel(1, E=KMUL * E0, SE=KMUL * SE0)

    def spec_nl(kEe, kSe):
        ev1 = [np.concatenate([u, [cb, ec, kh / e]]) for u, e in zip(two.uvs, np.broadcast_to(kEe, (len(two.tri),)))]
        ev2 = [np.concatenate([u, [B, B], [ec, kh / e]]) for u, e in zip(two.uvs, np.broadcast_to(kEe, (len(two.tri),)))]
        f1, g1, _, p1, _ = ef.assemble(specs['stretch_nl'], two.tri, ev1, two.X, np.broadcast_to(kEe, (len(two.tri),)).copy())
        f2, g2, _, p2, _ = ef.assemble(specs['shear_nl'], two.tri, ev2, two.X, np.broadcast_to(kSe, (len(two.tri),)).copy())
        return f1 + f2, g1 + p1 + g2 + p2

    m = hook(Drape, dll, two, base, nl)
    out['nonlinear'] = errs(m, {'per_panel': spec_nl(kE, kS), 'uniform_base': spec_nl(E0, SE0)})
    print('  %-10s 힘 — 패널별 %.1e · 균일(바탕) %.1e | J — %.1e · %.1e' % ('비선형', out['nonlinear']['per_panel']['f_rel'],
          out['nonlinear']['uniform_base']['f_rel'], out['nonlinear']['per_panel']['J_rel'], out['nonlinear']['uniform_base']['J_rel']))

    # ---- M5ⓑ: 전단 접선 다항식 + 패널별 SV + 속도 ----
    off_sv = dict(mass=0.1, E=0.0, V=0.0, SE=0.0, SV=SV0, BE_U=0.0, BV=0.0)

    def poly_sv(d):
        d.set_shear_tangent_polynomial(fx.COEF, fx.MAX_ANGLE)
        d.set_material_panel(1, SV=KDMUL * SV0)

    fpo, Jgg, JHC, _, _ = ef.assemble(specs['poly'], two.tri, [np.concatenate([u, fx.COEF]) for u in two.uvs], two.X, 1.0)
    m = hook(Drape, dll, two, off_sv, poly_sv, vel=two.vel)
    cands = {}
    for n, kd in (('per_panel', two.per_elem(SV0, {1: KDMUL * SV0})), ('uniform_base', SV0)):
        dmp, Jr = rate_terms_e(specs['q'], two.tri, two.uvs, two.X, two.vel, kd)
        cands[n] = (fpo + dmp, Jgg + JHC + Jr)
    out['polynomial_sv'] = errs(m, cands)
    print('  %-10s 힘 — 패널별 %.1e · 균일(바탕) %.1e | J — %.1e · %.1e' % ('다항식 SV', out['polynomial_sv']['per_panel']['f_rel'],
          out['polynomial_sv']['uniform_base']['f_rel'], out['polynomial_sv']['per_panel']['J_rel'], out['polynomial_sv']['uniform_base']['J_rel']))

    # ---- M6: 재질이 균일할 때 (비트 견주기용 해시) ----
    m0 = hook(Drape, dll, two, base)
    mv = hook(Drape, dll, two, base, vel=two.vel)
    out['uniform'] = {'f_sha': sha(m0['f']), 'J_sha': sha(m0['J']), 'f_vel_sha': sha(mv['f']), 'J_vel_sha': sha(mv['J']),
                      'free': {'grid 5x5': fx.free_run(Drape, dll, Vg, Fg, Xg), 'grid 17x17': fx.free_run(Drape, dll, V17, F17, X17)}}

    # ---- 보고: 요소가 두 재질에 걸칠 때 ----
    m = hook(Drape, dll, straddle, base, panel(1, E=KMUL * E0, SE=KMUL * SE0))
    zs = np.zeros_like(Xg)
    cands = {}
    for n, which in (('first_node', 'first_node'), ('element_range', 'elem_range')):
        cands[n] = spec_membrane(specs, straddle, zs, straddle.per_elem(E0, {1: KMUL * E0}, which), 0.0,
                                 straddle.per_elem(SE0, {1: KMUL * SE0}, which), 0.0)[:2]
    cands['uniform_base'] = spec_membrane(specs, straddle, zs, E0, 0.0, SE0, 0.0)[:2]
    out['straddle'] = errs(m, cands)
    out['straddle']['elements_where_rules_differ'] = int((straddle.panel_of_first_node != straddle.panel_of_elem_range).sum())
    with open(out_path, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    return 0


# ---- 부모 ---------------------------------------------------------------------------------------------
def run_arm(arm, dll, tmp):
    out = os.path.join(tmp, arm + '.json')
    print('팔 %s' % arm, flush=True)
    proc = subprocess.run([sys.executable, os.path.abspath(__file__), '--child', dll, out],
                          env=fx.clean_env(fx.KILL if arm.endswith('_kill') else None), cwd=ROOT)
    if proc.returncode != 0 or not os.path.isfile(out):
        raise RuntimeError('팔 %s 이(가) 끝나지 못했다 (종료 코드 %s)' % (arm, proc.returncode))
    with open(out, encoding='utf-8') as fh:
        return json.load(fh)


def engine_test(dll):
    """엔진 시험 py_pcg_parallel_test.py 다섯 장면 (M7). 장면마다 1 스레드 끝 위치 해시를 적는다."""
    proc = subprocess.run([ENGINE_PY, os.path.join('scripts', 'py_pcg_parallel_test.py')],
                          env=fx.clean_env({'CFMS_DRAPE_DLL': dll}), cwd=os.path.join(LINK_DIR, 'drape'),
                          capture_output=True, text=True, encoding='utf-8', errors='replace')
    hashes = dict(re.findall(r'^(\w+)\s+1T off ([0-9a-f]{16}) \|', proc.stdout, flags=re.M))
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    return {'returncode': proc.returncode, 'hashes': hashes, 'stdout_tail': lines[-1:], 'pass': bool(proc.returncode == 0)}


def git(*args):
    return subprocess.run(['git', '-C', LINK_DIR] + list(args), capture_output=True).stdout


def ok_per_panel(r):
    return bool(r['per_panel']['f_rel'] <= ef.LIM_SAME and r['per_panel']['J_rel'] <= ef.LIM_SAME)


def judge(arms, out):
    j = {}
    changed = sorted(x for x in git('diff', '--name-only', BRANCH_COMMIT).decode().split('\n') if x)
    head = git('rev-parse', '--short', 'HEAD').decode().strip()
    tree_dirty = bool(git('status', '--porcelain', '--untracked-files=no').strip())
    out['engine_branch'] = {'head': head, 'worktree_dirty': tree_dirty, 'changed_since_' + BRANCH_COMMIT: changed,
                            'diff_sha256': hashlib.sha256(git('diff', BRANCH_COMMIT, '--', *sorted(FIXED_FILES))).hexdigest()}
    bi, ai = arms['before']['build_info'], arms['after']['build_info']
    after_build = bool(ai['git_commit'] == head and ai['dirty'] is tree_dirty and tree_dirty == (head == BRANCH_COMMIT))
    j['M0'] = bool(bi['git_commit'] == BRANCH_COMMIT and bi['dirty'] is False and set(changed) - {'CLAUDE.md'} == FIXED_FILES
                   and after_build and arms['after']['dll_sha256'] != arms['before']['dll_sha256']
                   and arms['before_kill']['dll_sha256'] == arms['before']['dll_sha256']
                   and arms['after_kill']['dll_sha256'] == arms['after']['dll_sha256'])
    j['M1'] = {}
    for arm in ('before', 'before_kill'):
        s = arms[arm]['static']
        j['M1'][arm] = bool(all(s[n]['uniform_base']['f_rel'] <= ef.LIM_SAME and s[n]['uniform_base']['J_rel'] <= ef.LIM_SAME
                                for n in ('p1_E', 'p1_SE', 'p1_E_SE'))
                            and s['p0_E_SE']['uniform_override']['f_rel'] <= ef.LIM_SAME
                            and s['p0_E_SE']['uniform_override']['J_rel'] <= ef.LIM_SAME)
    j['M4_premise_a'] = bool(all(arms[a]['step']['uniform']['pcg_ok'] and arms[a]['step']['uniform']['v_rel'] <= LIM_STEP
                                 for a in ('before', 'before_kill')))
    for arm in ('after', 'after_kill'):
        r, b = arms[arm], arms['before' + arm[5:]]
        st = r['step']['per_panel']
        m6 = (r['uniform']['f_sha'] == b['uniform']['f_sha'] and r['uniform']['J_sha'] == b['uniform']['J_sha']
              and r['uniform']['f_vel_sha'] == b['uniform']['f_vel_sha'] and r['uniform']['J_vel_sha'] == b['uniform']['J_vel_sha']
              and all(r['uniform']['free'][k]['sha'] == b['uniform']['free'][k]['sha'] and r['uniform']['free'][k]['finite']
                      for k in r['uniform']['free'])
              and r['step']['uniform']['v_sha'] == b['step']['uniform']['v_sha'])
        j[arm] = {'M2': bool(all(ok_per_panel(v) for v in r['static'].values())),
                  'M3': ok_per_panel(r['velocity']),
                  'M4': bool(j['M4_premise_a'] and st['pcg_ok'] and st['v_rel'] <= LIM_STEP), 'M4_judged': bool(st['pcg_ok']),
                  'M5': bool(ok_per_panel(r['nonlinear']) and ok_per_panel(r['polynomial_sv'])),
                  'M6': bool(m6)}
    x = out['M7']
    j['M7'] = bool(x['after']['pass'] and len(x['after']['hashes']) == 5 and x['after']['hashes'] == x['before']['hashes'])
    j['pass'] = bool(j['M0'] and all(j['M1'].values()) and j['M7']
                     and all(j[arm][m] for arm in ('after', 'after_kill') for m in ('M2', 'M3', 'M4', 'M5', 'M6')))
    return j


def main():
    if len(sys.argv) == 4 and sys.argv[1] == '--child':
        return child(os.path.abspath(sys.argv[2]), sys.argv[3])
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    before, after = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    dlls = {'before': before, 'before_kill': before, 'after': after, 'after_kill': after}
    if not (os.path.isfile(before) and os.path.isfile(after) and os.path.isdir(os.path.join(LINK_DIR, 'drape'))):
        print('DLL · 엔진 브랜치 워크트리 중 없는 것이 있음: %s, %s, %s' % (before, after, LINK_DIR))
        return 2
    out = {'branch_commit': BRANCH_COMMIT, 'base': {'E': E0, 'SE': SE0, 'V': V0, 'SV': SV0}, 'k_mul': KMUL, 'kd_mul': KDMUL, 'law': list(LAW),
           'step': {'h': fx.STEP_H, 'mass': fx.STEP_MASS, 'pcg_rel': fx.PCG_REL, 'pcg_cap': fx.PCG_CAP, 'pcg_abs': fx.PCG_ABS, 'lim': LIM_STEP},
           'lim_same': ef.LIM_SAME, 'kill': fx.KILL}
    with tempfile.TemporaryDirectory() as tmp:
        arms = {arm: run_arm(arm, dlls[arm], tmp) for arm in ARMS}
    out['arms'] = arms
    print('엔진 시험 py_pcg_parallel_test.py (다섯 장면)', flush=True)
    out['M7'] = {'python': ENGINE_PY.replace('\\', '/'), 'before': engine_test(before), 'after': engine_test(after)}
    for arm in ('before', 'after'):
        print('  %-7s %s  %s' % (arm, 'PASS' if out['M7'][arm]['pass'] else 'FAIL', out['M7'][arm]['hashes']))
    j = out['judge'] = judge(arms, out)
    with open(os.path.join(ROOT, 'results', 'gate_engine_panel_material_fix.json'), 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
        fh.write('\n')
    print('\n판정')
    print('  M0 빌드 전제 %s | M1 결함이 있다 %s | M4 전제 ⓐ %s' % (j['M0'], j['M1'], j['M4_premise_a']))
    for arm in ('after', 'after_kill'):
        print('  %-10s M2 힘 · J %s | M3 속도 · 점성 %s | M4 한 단계 %s | M5 비선형 · 다항식 %s | M6 균일하면 그대로 %s'
              % (arm, j[arm]['M2'], j[arm]['M3'], j[arm]['M4'], j[arm]['M5'], j[arm]['M6']))
    print('  M7 엔진 시험 %s' % j['M7'])
    print('\n패널별 재질 수정: %s → results/gate_engine_panel_material_fix.json' % ('PASS' if j['pass'] else 'FAIL'))
    return 0 if j['pass'] else 1


if __name__ == '__main__':
    sys.exit(main())
