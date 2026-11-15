"""전단 접선 다항식의 노드쌍 블록 방향 수정 — 전 · 뒤 견주기 (docs/GATE_ENGINE_FORCES.md 「수정」). run_gate.bat 밖.

엔진 브랜치 sxjacobian-link 의 `_GetImplicitShear_TangentPolynomial_WithInvariant` 가 노드쌍 블록(J · V)을 슬롯의 방향
(행 = 작은 노드 번호)에 맞춰 넣는지를 본다. 팔 셋을 저마다 다른 프로세스로 돌린다:
  before = 고치기 전 DLL(브랜치 961765d),  after = 고친 DLL,  after_kill = 고친 DLL + CFMS_CPU_SXJ_BEND=off · CFMS_CPU_SXJ_MEMBRANE=off
기준 X0 – X6 은 판정문에 실행 전에 고정했다. 생성기 쪽 값은 tests/engine_forces.py 와 같은 길(sxgen.derive 를 Python 에서 수치 평가)이다.

사용:  python tests/engine_shear_poly_fix.py <고치기 전 cfmsDrape.dll 의 사본> <고친 cfmsDrape.dll 의 사본>
엔진의 Python 포장 · 시험 스크립트는 CFMS_DRAPE_LINK_DIR(없으면 ..\\cfmsDrape_dev_sxjacobian-link)에서 읽는다(읽기만).
X6 은 엔진의 가상 환경(CFMS_DRAPE_MAIN_DIR, 없으면 ..\\cfmsDrape_dev 의 drape\\.venv)의 python 으로 엔진 시험을 돌린다.
"""
import sys
sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import types  # noqa: E402

import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import engine_forces as ef  # noqa: E402
import sxgen  # noqa: E402  (engine_forces 가 codegen 을 경로에 넣는다)

LINK_DIR = os.path.abspath(os.environ.get('CFMS_DRAPE_LINK_DIR') or os.path.join(ROOT, '..', 'cfmsDrape_dev_sxjacobian-link'))
MAIN_DIR = os.path.abspath(os.environ.get('CFMS_DRAPE_MAIN_DIR') or os.path.join(ROOT, '..', 'cfmsDrape_dev'))
ENGINE_PY = os.path.join(MAIN_DIR, 'drape', '.venv', 'Scripts', 'python.exe')
BRANCH_COMMIT = '961765d'
FIXED_FILE = 'drape/src/sde_core/SDECore.cpp'

COEF = [4e4, 3e5, 2e6]
MAX_ANGLE = 1.2
SV = 1e3                      # X3b · X4 의 전단 점성
STEP_H, STEP_MASS = 1e-3, 1.0
PCG_REL, PCG_CAP, PCG_ABS = 1e-22, 500, 1e-40
LIM_STEP = 1e-9
KILL = {'CFMS_CPU_SXJ_BEND': 'off', 'CFMS_CPU_SXJ_MEMBRANE': 'off'}
ARMS = ('before', 'after', 'after_kill')


def sha(a):
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.float64).tobytes()).hexdigest()


# ---- 감쇠 몫의 조건 C_q = √a q (엔진: d = −α SV q̇ ∇q, V = −α SV ∇q∇qᵀ, J 에 −α SV q̇ ∇²q) ------------------
def _q_condition(ctx, X, E):
    x0, x1, x2 = X
    u0, v0, u1, v1, u2, v2 = E['pUV']
    du1, dv1, du2, dv2 = u1 - u0, v1 - v0, u2 - u0, v2 - v0
    det = du1 * dv2 - du2 * dv1
    dx1, dx2 = x1 - x0, x2 - x0
    wu = ctx.stage((dx1 * dv2 - dx2 * dv1) / det)
    wv = ctx.stage((-dx1 * du2 + dx2 * du1) / det)
    return [sp.sqrt(sp.Abs(det) / 2) * 2 * wu.dot(wv) / (wu.dot(wu) + wv.dot(wv))]


def q_spec_eval():
    spec = types.SimpleNamespace(NAME='ShearQ', N_NODES=3, N_COMP=1, ELEM_DOUBLE_VARS=[('pUV', 6)], condition=_q_condition)
    se = ef.SpecEval.__new__(ef.SpecEval)
    se.path, se.spec = None, spec
    xs, ev, C, g, H, _ = sxgen.derive(spec, 'exact')
    se.nd, se.nc, se.ne = len(xs), len(C), len(ev)
    se.nh = se.nd * (se.nd + 1) // 2
    se.fn = sp.lambdify(list(xs) + list(ev), [C, g, H], modules='numpy', cse=True)
    se.iu = np.triu_indices(se.nd)
    return se


def rate_terms(se, stencils, evs, X, vel, kd):
    """d = −kd Σ Ċ ∇C,  V = −kd Σ ∇C∇Cᵀ,  J_rate = −kd Σ Ċ ∇²C  (전역, 조밀)."""
    n3 = X.size
    d, Vm, Jr = np.zeros(n3), np.zeros((n3, n3)), np.zeros((n3, n3))
    for s, ev in zip(stencils, evs):
        idx = np.array([3 * v + c for v in s for c in range(3)])
        _, g, H = se(X[list(s)].ravel(), ev)
        vl = vel[list(s)].ravel()
        ix = np.ix_(idx, idx)
        for c in range(se.nc):
            cd = float(g[c] @ vl)
            d[idx] += -kd * cd * g[c]
            Vm[ix] += -kd * np.outer(g[c], g[c])
            Jr[ix] += -kd * cd * H[c]
    return d, Vm, Jr


def poly(d):
    d.set_shear_tangent_polynomial(COEF, MAX_ANGLE)


def one_step(Drape, dll, V, F, uv, X0, V0):
    """다항식 전단만 켜고(SV > 0) 속도를 준 상태에서 한 단계. (v1, x1, PCG 기록)."""
    quiet = dict(mass=STEP_MASS, E=0.0, V=0.0, SE=0.0, SV=SV, BE_U=0.0, BV=0.0, air_damping=0.0, internal_damping=0.0)
    d = ef.new_drape(Drape, dll, V, F, uv, quiet, poly)
    try:
        d.set_timestep(STEP_H)
        d.set_air_drag(0.0)
        d.set_damping(0.0)
        d.set_dynamic_pcg(True, PCG_REL, PCG_CAP)
        d.set_dynamic_pcg_absolute_tolerance(PCG_ABS)
        d.step(1)                          # 설정 뒤 첫 단계는 준비 단계다 — 쉬는 모양에서 먼저 한 번 밟는다
        d.set_positions(X0)
        d.set_velocities(V0)
        d.step(1)
        return d.velocities().copy(), d.positions().copy(), d.get_pcg_solve_history()
    finally:
        d.close()


def step_case(Drape, dll, sp_poly, sp_q, V, F, X0, V0):
    uv = V[:, :2].copy()
    tri = [tuple(int(v) for v in t) for t in F]
    uvs = [uv[list(t)].ravel() for t in tri]
    v1, x1, hist = one_step(Drape, dll, V, F, uv, X0, V0)
    f, Jgg, JHC, _, _ = ef.assemble(sp_poly, tri, [np.concatenate([u, COEF]) for u in uvs], X0, 1.0)
    dmp, Vm, Jr = rate_terms(sp_q, tri, uvs, X0, V0, SV)
    J = Jgg + JHC + Jr
    n3 = X0.size
    A = STEP_MASS * np.eye(n3) - STEP_H * Vm - STEP_H * STEP_H * J        # 식 (6)
    b = STEP_H * (f + dmp) + STEP_H * STEP_H * (J @ V0.ravel())
    v_pred = V0.ravel() + np.linalg.solve(A, b)
    eig = np.linalg.eigvalsh(0.5 * (A + A.T))
    reasons = [r['exit_reason'] for r in hist]
    return {'v_rel': ef.rel(v1.ravel(), v_pred), 'dv_rel_to_v': float(np.abs(v_pred - V0.ravel()).max() / np.abs(V0).max()),
            'x_abs': float(np.abs(x1.ravel() - (X0.ravel() + STEP_H * v_pred)).max()),
            'pcg_exit': reasons, 'pcg_iterations': [int(r['iterations']) for r in hist],
            'pcg_ok': bool(len(reasons) > 0 and all(r == 'relative' for r in reasons)),
            'A_min_eig': float(eig[0]), 'A_max_eig': float(eig[-1]), 'A_asym': float(np.abs(A - A.T).max())}


def free_run(Drape, dll, V, F, X):
    """다항식을 켜지 않은 기본 설정 — 흩뜨린 천을 40 단계."""
    d = Drape(dll)
    try:
        d.set_mesh(V, F)
        d.set_mesh_2d(V[:, :2].copy())
        d.set_timestep(0.005)
        d.step(1)
        d.set_positions(X)
        d.step(40)
        p = np.ascontiguousarray(d.positions(), dtype=np.float64)
        return {'sha': sha(p), 'finite': bool(np.isfinite(p).all()), 'moved': float(np.abs(p - X).max())}
    finally:
        d.close()


def child(dll, out_path):
    import warnings
    warnings.simplefilter('ignore')
    ef.ENGINE_DIR = LINK_DIR
    Drape = ef.load_drape()
    rng = np.random.default_rng(ef.SEED)
    V, F = ef.grid_scene()
    uv = V[:, :2].copy()
    X = ef.deformed(V, rng)                      # engine_forces 의 D1 과 같은 상태
    vel = rng.uniform(-5.0, 5.0, X.shape)
    V1 = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.5, 1.7, 0.0]])
    X1 = np.array([[0.1, -0.05, 0.02], [2.25, 0.2, 0.3], [0.95, 1.6, -0.25]])
    vel1 = rng.uniform(-5.0, 5.0, X1.shape)
    V17, F17 = ef.grid_scene(16)
    X17 = ef.deformed(V17, rng)
    vel17 = rng.uniform(-5.0, 5.0, X17.shape)
    tri = [tuple(int(v) for v in t) for t in F]
    uvs = [uv[list(t)].ravel() for t in tri]
    off = dict(mass=0.1, E=0.0, V=0.0, SE=0.0, SV=0.0, BE_U=0.0, BV=0.0)
    off_sv = dict(off, SV=SV)
    sp_poly, sp_q = ef.SpecEval('shear_polynomial'), q_spec_eval()

    out = {'dll_sha256': hashlib.sha256(open(dll, 'rb').read()).hexdigest()}
    d0 = ef.new_drape(Drape, dll, V, F, uv, off)
    out['build_info'] = d0.build_info
    out['threads'] = int(d0.openmp_threads)
    d0.close()
    x3 = X[list(tri[3])].ravel()
    out['H1'] = {'shear_polynomial': list(sp_poly.fd_check(x3, np.concatenate([uvs[3], COEF]))),
                 'shear_q': list(sp_q.fd_check(x3, uvs[3]))}

    # ---- D1 격자: 속도 0 (X2 · X3), 속도 · SV (X3b) ----
    evs = [np.concatenate([u, COEF]) for u in uvs]
    f, Jgg, JHC, JHCp, _ = ef.assemble(sp_poly, tri, evs, X, 1.0)
    m = ef.measure(Drape, dll, V, F, uv, X, off, poly)
    out['D1'] = ef.compare('D1 격자', m['f'], m['J'], m['Jfd'], f, Jgg, JHC, JHCp)
    out['D1'].update(f_sha=sha(m['f']), J_asymmetry=float(np.abs(m['J'] - m['J'].T).max() / np.abs(m['J']).max()))
    dmp, _, Jr = rate_terms(sp_q, tri, uvs, X, vel, SV)
    m = ef.measure(Drape, dll, V, F, uv, X, off_sv, poly, vel=vel, fd=False)
    out['D1v'] = {'f_rel': ef.rel(m['f'], f + dmp), 'J_rel': ef.rel(m['J'], Jgg + JHC + Jr),
                  'damping_share_of_f': float(np.abs(dmp).max() / np.abs(f + dmp).max()),
                  'rate_share_of_J': float(np.abs(Jr).max() / np.abs(Jgg + JHC + Jr).max())}
    print('  %-36s 힘 %.1e | J %.1e' % ('D1 격자, 속도 · SV', out['D1v']['f_rel'], out['D1v']['J_rel']))

    # ---- D2 삼각형 하나, 꼭짓점 순서 여섯 ----
    out['D2'], out['D2v'] = {}, {}
    for perm in itertools.permutations(range(3)):
        F1 = np.array([list(perm)], np.int32)
        uv1 = V1[:, :2].copy()
        ev = np.concatenate([V1[list(perm), :2].ravel(), COEF])
        f, Jgg, JHC, JHCp, _ = ef.assemble(sp_poly, [perm], [ev], X1, 1.0)
        H = Jgg + JHC
        key = '(%d, %d, %d)' % perm
        m = ef.measure(Drape, dll, V1, F1, uv1, X1, off, poly)
        r = ef.compare('D2 삼각형 ' + key, m['f'], m['J'], m['Jfd'], f, Jgg, JHC, JHCp)
        r['f_sha'] = sha(m['f'])
        blocks = {}
        for p in range(3):
            for q in range(p + 1, 3):
                Jb, Hb = m['J'][3 * p:3 * p + 3, 3 * q:3 * q + 3], H[3 * p:3 * p + 3, 3 * q:3 * q + 3]
                blocks['%d-%d' % (p, q)] = {'vs_exact': ef.rel(Jb, Hb), 'vs_exact_transposed': ef.rel(Jb, Hb.T),
                                            'local_order_reversed': bool(perm.index(p) > perm.index(q))}
        r['upper_blocks'] = blocks
        out['D2'][key] = r
        dmp, _, Jr = rate_terms(sp_q, [perm], [ev[:6]], X1, vel1, SV)
        m = ef.measure(Drape, dll, V1, F1, uv1, X1, off_sv, poly, vel=vel1, fd=False)
        out['D2v'][key] = {'f_rel': ef.rel(m['f'], f + dmp), 'J_rel': ef.rel(m['J'], H + Jr)}
        print('  %-36s 힘 %.1e | J %.1e' % ('D2 삼각형 ' + key + ', 속도 · SV', out['D2v'][key]['f_rel'], out['D2v'][key]['J_rel']))

    # ---- X4 한 단계 (V 블록) ----
    out['step'] = {}
    for perm in itertools.permutations(range(3)):
        key = 'tri (%d, %d, %d)' % perm
        out['step'][key] = step_case(Drape, dll, sp_poly, sp_q, V1, np.array([list(perm)], np.int32), X1, vel1)
    out['step']['grid 5x5'] = step_case(Drape, dll, sp_poly, sp_q, V, F, X, vel)
    out['step']['grid 17x17'] = step_case(Drape, dll, sp_poly, sp_q, V17, F17, X17, vel17)
    for key, r in out['step'].items():
        print('  %-36s 한 단계 뒤 속도 %.1e (PCG %s, 반복 %s)' % ('X4 ' + key, r['v_rel'], r['pcg_exit'], r['pcg_iterations']))

    # ---- X5 옵션을 끈 기본 설정 ----
    out['free'] = {'grid 5x5': free_run(Drape, dll, V, F, X), 'grid 17x17': free_run(Drape, dll, V17, F17, X17)}
    with open(out_path, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    return 0


# ---- 부모: 팔을 돌리고 판정한다 -------------------------------------------------------------------------
def clean_env(extra=None):
    env = {k: v for k, v in os.environ.items()
           if not (k.startswith('CFMS_') and k not in ('CFMS_DRAPE_LINK_DIR', 'CFMS_DRAPE_MAIN_DIR'))
           and k not in ('OMP_NUM_THREADS',)}
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env.update(extra or {})
    return env


def run_arm(arm, dll, tmp):
    out = os.path.join(tmp, arm + '.json')
    print('팔 %s' % arm, flush=True)
    proc = subprocess.run([sys.executable, os.path.abspath(__file__), '--child', dll, out],
                          env=clean_env(KILL if arm == 'after_kill' else None), cwd=ROOT)
    if proc.returncode != 0 or not os.path.isfile(out):
        raise RuntimeError('팔 %s 이(가) 끝나지 못했다 (종료 코드 %s)' % (arm, proc.returncode))
    with open(out, encoding='utf-8') as fh:
        return json.load(fh)


def engine_test(dll):
    """엔진 시험 py_pcg_parallel_test.py 의 장면 shear_poly (X6)."""
    proc = subprocess.run([ENGINE_PY, os.path.join('scripts', 'py_pcg_parallel_test.py'), '--scenes', 'shear_poly'],
                          env=clean_env({'CFMS_DRAPE_DLL': dll}), cwd=os.path.join(LINK_DIR, 'drape'),
                          capture_output=True, text=True, encoding='utf-8', errors='replace')
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    return {'returncode': proc.returncode, 'stdout_tail': lines[-4:], 'pass': bool(proc.returncode == 0)}


def git(*args):
    return subprocess.run(['git', '-C', LINK_DIR] + list(args), capture_output=True).stdout


def judge(arms, out):
    b, a = arms['before'], arms['after']
    j = {}
    # 뒤 DLL: 961765d 와 다른 빌드 대상 파일이 SDECore.cpp 하나(엔진 CLAUDE.md 는 문서라 센다 — 판정문 「커밋 뒤」),
    # 그리고 DLL 이 지금 워크트리 상태에서 빌드된 것(커밋 전 = 961765d · dirty, 커밋 뒤 = HEAD · 깨끗).
    changed = sorted(x for x in git('diff', '--name-only', BRANCH_COMMIT).decode().split('\n') if x)
    head = git('rev-parse', '--short', 'HEAD').decode().strip()
    tree_dirty = bool(git('status', '--porcelain', '--untracked-files=no').strip())
    out['engine_branch'] = {'head': head, 'worktree_dirty': tree_dirty, 'changed_since_' + BRANCH_COMMIT: changed,
                            'diff_sha256': hashlib.sha256(git('diff', BRANCH_COMMIT, '--', FIXED_FILE)).hexdigest()}
    ai = a['build_info']
    after_build = bool(ai['git_commit'] == head and ai['dirty'] is tree_dirty and tree_dirty == (head == BRANCH_COMMIT))
    j['X0'] = bool(b['build_info']['git_commit'] == BRANCH_COMMIT and b['build_info']['dirty'] is False
                   and set(changed) - {'CLAUDE.md'} == {FIXED_FILE} and after_build and a['dll_sha256'] != b['dll_sha256']
                   and arms['after_kill']['dll_sha256'] == a['dll_sha256'])
    rule = [bl for r in b['D2'].values() for bl in r['upper_blocks'].values()]
    agree = sum(1 for bl in rule if (bl['vs_exact_transposed'] <= ef.LIM_SAME) == bl['local_order_reversed']
                and (bl['vs_exact'] <= ef.LIM_SAME) == (not bl['local_order_reversed']))
    out['before_block_rule'] = {'blocks': len(rule), 'transposed_exactly_when_local_order_reversed': agree}
    j['X1'] = bool(not b['D1']['J_is_derivative'] and agree == len(rule) == 18)
    step_premise = bool(b['step']['tri (0, 1, 2)']['pcg_ok'] and b['step']['tri (0, 1, 2)']['v_rel'] <= LIM_STEP)
    out['X4_premise_a'] = step_premise
    for arm in ('after', 'after_kill'):
        r = arms[arm]
        x2 = (r['D1']['f_sha'] == b['D1']['f_sha'] and r['D1']['f_same']
              and all(r['D2'][k]['f_sha'] == b['D2'][k]['f_sha'] and r['D2'][k]['f_same'] for k in r['D2']))
        x3 = (r['D1']['J_best'] == 'exact' and r['D1']['J_same'] and r['D1']['J_is_derivative']
              and all(v['J_rel']['exact'] <= ef.LIM_SAME for v in r['D2'].values())
              and all(bl['vs_exact'] <= ef.LIM_SAME for v in r['D2'].values() for bl in v['upper_blocks'].values()))
        x3b = (r['D1v']['f_rel'] <= ef.LIM_SAME and r['D1v']['J_rel'] <= ef.LIM_SAME
               and all(v['f_rel'] <= ef.LIM_SAME and v['J_rel'] <= ef.LIM_SAME for v in r['D2v'].values()))
        unjudged = [k for k, v in r['step'].items() if not v['pcg_ok']]
        x4 = bool(step_premise and not unjudged and all(v['v_rel'] <= LIM_STEP for v in r['step'].values()))
        j[arm] = {'X2': bool(x2), 'X3': bool(x3), 'X3b': bool(x3b), 'X4': x4, 'X4_unjudged': unjudged}
    j['X5'] = bool(all(a['free'][k]['sha'] == b['free'][k]['sha'] and a['free'][k]['finite'] for k in a['free']))
    j['X6'] = bool(out['X6']['after']['pass'])
    j['pass'] = bool(j['X0'] and j['X1'] and j['X5'] and j['X6']
                     and all(j[arm][x] for arm in ('after', 'after_kill') for x in ('X2', 'X3', 'X3b', 'X4')))
    return j


def main():
    if len(sys.argv) == 4 and sys.argv[1] == '--child':
        return child(os.path.abspath(sys.argv[2]), sys.argv[3])
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    dlls = {'before': os.path.abspath(sys.argv[1]), 'after': os.path.abspath(sys.argv[2])}
    dlls['after_kill'] = dlls['after']
    if not all(os.path.isfile(p) for p in dlls.values()) or not os.path.isdir(os.path.join(LINK_DIR, 'drape')):
        print('DLL · 엔진 브랜치 워크트리 중 없는 것이 있음: %s, %s' % (dlls, LINK_DIR))
        return 2
    out = {'branch_commit': BRANCH_COMMIT, 'coef': COEF, 'max_angle': MAX_ANGLE, 'SV': SV,
           'step': {'h': STEP_H, 'mass': STEP_MASS, 'pcg_rel': PCG_REL, 'pcg_cap': PCG_CAP, 'pcg_abs': PCG_ABS, 'lim': LIM_STEP},
           'lim_same': ef.LIM_SAME, 'lim_fd': ef.LIM_SPEC_FD, 'kill': KILL}
    with tempfile.TemporaryDirectory() as tmp:
        arms = {arm: run_arm(arm, dlls[arm], tmp) for arm in ARMS}
    out['arms'] = arms
    print('엔진 시험 py_pcg_parallel_test.py --scenes shear_poly', flush=True)
    out['X6'] = {'python': ENGINE_PY.replace('\\', '/'), 'before': engine_test(dlls['before']), 'after': engine_test(dlls['after'])}
    for arm in ('before', 'after'):
        print('  %-7s %s  %s' % (arm, 'PASS' if out['X6'][arm]['pass'] else 'FAIL', out['X6'][arm]['stdout_tail'][-1:]))
    j = out['judge'] = judge(arms, out)
    with open(os.path.join(ROOT, 'results', 'gate_engine_shear_poly_fix.json'), 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
        fh.write('\n')
    print('\n판정')
    print('  X0 빌드 전제 %s | X1 결함이 브랜치에도 있다 %s (블록 %d / %d) | X4 전제 ⓐ %s'
          % (j['X0'], j['X1'], out['before_block_rule']['transposed_exactly_when_local_order_reversed'],
             out['before_block_rule']['blocks'], out['X4_premise_a']))
    for arm in ('after', 'after_kill'):
        print('  %-10s X2 힘 %s | X3 J %s | X3b J(속도) %s | X4 한 단계 %s%s'
              % (arm, j[arm]['X2'], j[arm]['X3'], j[arm]['X3b'], j[arm]['X4'],
                 ' (판정 못 한 회: %s)' % j[arm]['X4_unjudged'] if j[arm]['X4_unjudged'] else ''))
    print('  X5 옵션을 끄면 그대로 %s | X6 엔진 시험 %s' % (j['X5'], j['X6']))
    print('\n전단 접선 다항식 수정: %s → results/gate_engine_shear_poly_fix.json' % ('PASS' if j['pass'] else 'FAIL'))
    return 0 if j['pass'] else 1


if __name__ == '__main__':
    sys.exit(main())
