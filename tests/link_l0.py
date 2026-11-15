"""엔진 연결 L0 (docs/GATE_LINK_L0.md): cfmsDrape CPU DLL 의 디버그 고리로 엔진 규약을 잰다. 엔진 소스는 고치지 않는다.

사용:  python tests/link_l0.py <cfmsDrape.dll> <link_l0_ref.exe>
       (내부) python tests/link_l0.py --dump <cfmsDrape.dll> <out.npz> — L0-4 (ii) 의 OMP_NUM_THREADS = 1 프로세스
엔진의 Python 포장 cfmsdrape 는 CFMS_DRAPE_DIR(없으면 기준 워크트리 ..\\cfmsDrape_dev_base378d20c, 그것도 없으면 ..\\cfmsDrape_dev)의 drape 에서 읽어 들인다. 엔진 리포에 바이트코드를
쓰지 않는다(sys.dont_write_bytecode). 리포 루트에서 실행. 결과는 results/gate_link_l0.json. 통과 0, 실패 1, 준비 안 됨 2.
"""
import sys
sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402

import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENGINE_COMMIT = '378d20c'
# 기준 엔진(손 코드) = 378d20c. 엔진 main 은 2026-11-13 에 생성 커널을 합쳐 그 뒤로 나아갔다. 그래서 그 커밋에 고정한
# 워크트리 ..\cfmsDrape_dev_base378d20c 가 있으면 그것을 읽는다(없으면 종전대로 ..\cfmsDrape_dev). CFMS_DRAPE_DIR 이 먼저다.
_BASE_DIR = os.path.join(ROOT, '..', 'cfmsDrape_dev_base' + ENGINE_COMMIT)
ENGINE_DIR = os.path.abspath(os.environ.get('CFMS_DRAPE_DIR')
                             or (_BASE_DIR if os.path.isdir(os.path.join(_BASE_DIR, 'drape')) else os.path.join(ROOT, '..', 'cfmsDrape_dev')))

# 장면 (판정문 「장면」에 고정)
M, S = 4, 2.0                 # 4 × 4 칸, 간격 2 cm
SEED = 20261113
N_RAND = 5
FD_STEP = 1e-4                # cm
MATERIALS = [('stretch', dict(E=1e5, SE=0.0, BE_U=0.0)),
             ('shear', dict(E=0.0, SE=1e5, BE_U=0.0)),
             ('bend', dict(E=0.0, SE=0.0, BE_U=1.0)),
             ('all', dict(E=1e5, SE=1e5, BE_U=1.0))]
K_STRETCH = 1e5               # = E
K_BEND = 1.0                  # = BE_U
LIM_PRODUCT = 1e-12           # L0-1
LIM_STRETCH = 1e-10           # L0-2 · L0-3
MIN_BRANCH = 10               # L0-3 갈래


def load_drape():
    sys.path.insert(0, os.path.join(ENGINE_DIR, 'drape'))
    from cfmsdrape import Drape
    return Drape


def scene():
    n = M + 1
    V = np.array([[S * i, S * j, 0.0] for j in range(n) for i in range(n)])
    F = []
    for j in range(M):
        for i in range(M):
            a = j * n + i
            F += [[a, a + 1, a + n + 1], [a, a + n + 1, a + n]]
    F = np.array(F, np.int32)
    uv = V[:, :2].copy()
    rng = np.random.default_rng(SEED)
    x, y = V[:, 0], V[:, 1]
    X = np.stack([1.08 * x + 0.10 * y, 0.03 * x + 0.94 * y,
                  0.6 * np.sin(np.pi * x / 8.0) * np.sin(np.pi * y / 8.0)], axis=1)
    X = X + rng.uniform(-0.15, 0.15, X.shape)
    R = rng.standard_normal((N_RAND,) + X.shape)
    return V, F, uv, X, R


def hinges(F):
    """MakeGrid(tests/gate_common.h)와 같은 순서의 힌지 (h, i, j, k)."""
    third = {}
    for t in F:
        for e in range(3):
            third[(int(t[e]), int(t[(e + 1) % 3]))] = int(t[(e + 2) % 3])
    H = []
    for t in F:
        for e in range(3):
            p, q = int(t[e]), int(t[(e + 1) % 3])
            if p > q:
                continue
            k = third.get((q, p))
            if k is None:
                continue
            H.append((int(t[(e + 2) % 3]), p, q, k))
    return np.array(H, np.int64)


def engine_measure(Drape, dll, V, F, uv, X, R, mat, full):
    d = Drape(dll)
    try:
        d.set_mesh(V, F)
        d.set_mesh_2d(uv)
        d.set_material(mass=0.1, **mat)
        d.set_gravity((0.0, 0.0, 0.0))
        d.set_timestep(0.005)
        out = {'threads': int(d.openmp_threads), 'build_info': d.build_info}
        n3 = X.size
        out['f'] = d.debug_internal_force(X).ravel()
        J = np.empty((n3, n3))
        e = np.zeros(n3)
        for k in range(n3):
            e[:] = 0.0
            e[k] = 1.0
            J[:, k] = d.debug_stiffness_product(X, e.reshape(-1, 3)).ravel()
        out['J'] = J
        if full:
            out['f2'] = d.debug_internal_force(X).ravel()
            out['P'] = np.array([d.debug_stiffness_product(X, r).ravel() for r in R])
            out['P2'] = np.array([d.debug_stiffness_product(X, r).ravel() for r in R])
            Xf = X.ravel()
            Jfd = np.empty((n3, n3))
            for k in range(n3):
                xp, xm = Xf.copy(), Xf.copy()
                xp[k] += FD_STEP
                xm[k] -= FD_STEP
                fp = d.debug_internal_force(xp.reshape(-1, 3)).ravel()
                fm = d.debug_internal_force(xm.reshape(-1, 3)).ravel()
                Jfd[:, k] = (fp - fm) / (xp[k] - xm[k])   # 실제로 표현된 간격
            out['Jfd'] = Jfd
        return out
    finally:
        d.close()


def run_ref(ref_exe, kind, rows):
    with tempfile.TemporaryDirectory() as td:
        inp, outp = os.path.join(td, 'in.txt'), os.path.join(td, 'out.txt')
        with open(inp, 'w') as fh:
            for r in rows:
                fh.write(' '.join('%.17g' % v for v in r) + '\n')
        subprocess.run([ref_exe, kind, inp, outp], check=True, capture_output=True)
        return np.loadtxt(outp, ndmin=2)


def upper_to_full(h, n):
    Hf = np.zeros((n, n))
    i = 0
    for a in range(n):
        for b in range(a, n):
            Hf[a, b] = Hf[b, a] = h[i]
            i += 1
    return Hf


def reference(ref_exe, kind, stencils, evs, X, k):
    """이 리포 기준: f = −k Σ C∇C, Σ −k ∇C∇Cᵀ, Σ −k C ∇²C, Σ −k C⁺ ∇²C (전역, 조밀)."""
    nd = 3 * stencils.shape[1]
    nc = 2 if kind == 'stretch' else 1
    nh = nd * (nd + 1) // 2
    rows = [np.concatenate([X[s].ravel(), np.atleast_1d(ev)]) for s, ev in zip(stencils, evs)]
    out = run_ref(ref_exe, kind, rows)
    n3 = X.size
    f = np.zeros(n3)
    Jgg, JHC, JHCp = np.zeros((n3, n3)), np.zeros((n3, n3)), np.zeros((n3, n3))
    Cs = []
    for s, row in zip(stencils, out):
        idx = np.array([3 * v + c for v in s for c in range(3)])
        C = row[:nc]
        g = row[nc:nc + nc * nd].reshape(nc, nd)
        H = row[nc + nc * nd:].reshape(nc, nh)
        ix = np.ix_(idx, idx)
        for c in range(nc):
            Hf = upper_to_full(H[c], nd)
            f[idx] += -k * C[c] * g[c]
            Jgg[ix] += -k * np.outer(g[c], g[c])
            JHC[ix] += -k * C[c] * Hf
            JHCp[ix] += -k * max(C[c], 0.0) * Hf
            Cs.append(C[c])
    return f, Jgg, JHC, JHCp, np.array(Cs)


def bridson_w(X, uv, h):
    """엔진처럼 그 상태의 w = 3 L0² / (½|nA| + ½|nB|), L0 = 쉬는 변 길이(uv)."""
    xh, xi, xj, xk = X[h[0]], X[h[1]], X[h[2]], X[h[3]]
    nA = np.cross(xi - xh, xj - xh)
    nB = np.cross(xj - xk, xi - xk)
    L0 = np.linalg.norm(uv[h[1]] - uv[h[2]])
    return 3.0 * L0 * L0 / (0.5 * np.linalg.norm(nA) + 0.5 * np.linalg.norm(nB))


def rel(a, b):
    s = float(np.abs(b).max())
    d = float(np.abs(a - b).max())
    return d / s if s > 0 else d


def same_bits(a, b):
    a, b = np.ascontiguousarray(a), np.ascontiguousarray(b)
    return a.shape == b.shape and a.tobytes() == b.tobytes()


def dump(dll, path):
    Drape = load_drape()
    V, F, uv, X, R = scene()
    arrs = {}
    for name, mat in MATERIALS:
        m = engine_measure(Drape, dll, V, F, uv, X, R, mat, full=False)
        arrs['f_' + name], arrs['J_' + name] = m['f'], m['J']
        arrs['threads'] = np.array([m['threads']])
    np.savez(path, **arrs)
    return 0


def gate(fail, ok, label, fmt, value, limit):
    print('  %-5s %s  %s  (기준 %s)' % (label, 'PASS' if ok else 'FAIL', fmt % value, limit))
    return fail + (0 if ok else 1)


def main():
    if len(sys.argv) == 4 and sys.argv[1] == '--dump':
        return dump(sys.argv[2], sys.argv[3])
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    dll, ref_exe = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    if not os.path.isfile(dll) or not os.path.isfile(ref_exe) or not os.path.isdir(os.path.join(ENGINE_DIR, 'drape')):
        print('DLL · 기준 도구 · 엔진 리포 중 없는 것이 있음: %s, %s, %s' % (dll, ref_exe, ENGINE_DIR))
        return 2
    Drape = load_drape()
    V, F, uv, X, R = scene()
    H = hinges(F)
    print('엔진 연결 L0 — 노드 %d, 삼각형 %d, 힌지 %d, DLL %s' % (len(V), len(F), len(H), dll))

    meas = {}
    for name, mat in MATERIALS:
        meas[name] = engine_measure(Drape, dll, V, F, uv, X, R, mat, full=True)
    bi = meas['stretch']['build_info']
    threads = meas['stretch']['threads']

    # L0-4 (ii): OMP_NUM_THREADS = 1 인 다른 프로세스
    with tempfile.TemporaryDirectory() as td:
        npz = os.path.join(td, 'single.npz')
        env = dict(os.environ, OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
        subprocess.run([sys.executable, os.path.abspath(__file__), '--dump', dll, npz], check=True, env=env, cwd=ROOT)
        single = dict(np.load(npz))
    threads1 = int(single['threads'][0])

    # 이 리포 기준
    tri = F.astype(np.int64)
    fS, JggS, JHCS, JHCpS, CS = reference(ref_exe, 'stretch', tri, [uv[t].ravel() for t in tri], X, K_STRETCH)
    wB = np.array([bridson_w(X, uv, h) for h in H])
    fB, JggB, JHCB, _, CB = reference(ref_exe, 'bridson', H, wB, X, K_BEND)
    J_T = JggS + JHCpS
    J_E = JggS + JHCS

    st = subprocess.run(['git', '-C', ENGINE_DIR, 'status', '--porcelain'], capture_output=True, text=True, check=True).stdout
    dll_sha = hashlib.sha256(open(dll, 'rb').read()).hexdigest()

    fail = 0
    print(' L0-0 전제')
    pre_ok = str(bi.get('git_commit', '')).startswith(ENGINE_COMMIT) and bi.get('dirty') is False
    fail = gate(fail, pre_ok, 'DLL', '%s', 'git_commit %s · dirty %s' % (bi.get('git_commit'), bi.get('dirty')), '%s · dirty false' % ENGINE_COMMIT)
    fail = gate(fail, st.strip() == '', '리포', '%s', 'git status --porcelain %s' % ('비어 있음' if not st.strip() else '줄 %d 개' % len(st.splitlines())), '비어 있음')

    print(' L0-1 곱의 일관성 (되살린 J 로 곱한 값 대 엔진 곱, 무작위 벡터 %d 개)' % N_RAND)
    l01 = {}
    for name, _ in MATERIALS:
        m = meas[name]
        l01[name] = max(rel(m['J'] @ r.ravel(), p) for r, p in zip(R, m['P']))
        fail = gate(fail, l01[name] <= LIM_PRODUCT, name, 'rel max %.2e', l01[name], '%.0e' % LIM_PRODUCT)

    print(' L0-2 · L0-3 인장 (k = E = %g)' % K_STRETCH)
    l02 = rel(meas['stretch']['f'], fS)
    l03 = rel(meas['stretch']['J'], J_T)
    nPos, nNeg = int((CS > 0).sum()), int((CS <= 0).sum())
    fail = gate(fail, l02 <= LIM_STRETCH, 'f', 'rel max %.2e', l02, '%.0e' % LIM_STRETCH)
    fail = gate(fail, l03 <= LIM_STRETCH, 'J', 'J 대 J_T rel max %.2e', l03, '%.0e' % LIM_STRETCH)
    fail = gate(fail, nPos >= MIN_BRANCH and nNeg >= MIN_BRANCH, '갈래', '%s', 'C > 0 %d 개 · C ≤ 0 %d 개' % (nPos, nNeg), '각 %d 개 이상' % MIN_BRANCH)

    print(' L0-4 결정성 (기본 스레드 %d · 단일 스레드 프로세스 %d)' % (threads, threads1))
    l04 = {}
    for name, _ in MATERIALS:
        m = meas[name]
        same_proc = same_bits(m['f'], m['f2']) and same_bits(m['P'], m['P2'])
        across = same_bits(m['f'], single['f_' + name]) and same_bits(m['J'], single['J_' + name])
        l04[name] = {'same_process': same_proc, 'threads_1_vs_default': across}
        fail = gate(fail, same_proc and across, name, '%s', '같은 프로세스 %s · 스레드 1 대 기본 %s' % (
            '같음' if same_proc else '다름', '같음' if across else '다름'), '비트까지 같음')

    rep = {}
    for name, _ in MATERIALS:
        m = meas[name]
        J = m['J']
        rep[name] = {'asym': rel(J, J.T) if np.abs(J).max() > 0 else 0.0, 'J_vs_fd': rel(J, m['Jfd']),
                     'J_max': float(np.abs(J).max()), 'f_max': float(np.abs(m['f']).max())}
    rep['stretch']['J_vs_exact'] = rel(meas['stretch']['J'], J_E)
    rep['stretch']['exact_vs_fd'] = rel(J_E, meas['stretch']['Jfd'])
    rep['bend']['f_vs_bridson_state_w'] = rel(meas['bend']['f'], fB)
    rep['bend']['J_vs_gn_state_w'] = rel(meas['bend']['J'], JggB)
    rep['bend']['J_vs_exact_state_w'] = rel(meas['bend']['J'], JggB + JHCB)
    rep['bend']['exact_state_w_vs_fd'] = rel(JggB + JHCB, meas['bend']['Jfd'])
    print(' 보고')
    for name, _ in MATERIALS:
        r = rep[name]
        print('  %-7s 비대칭 %.2e · J 대 유한차분 %.2e (max|J| %.3g, max|f| %.3g)' % (name, r['asym'], r['J_vs_fd'], r['J_max'], r['f_max']))
    print('  인장   엔진 J 대 이 리포 정확판 %.2e · 정확판 대 유한차분 %.2e' % (rep['stretch']['J_vs_exact'], rep['stretch']['exact_vs_fd']))
    print('  굽힘   엔진 f 대 Bridson 힘(w = 그 상태) %.2e · 엔진 J 대 GN %.2e · 대 정확판(w 고정) %.2e · 정확판(w 고정) 대 유한차분 %.2e' % (
        rep['bend']['f_vs_bridson_state_w'], rep['bend']['J_vs_gn_state_w'], rep['bend']['J_vs_exact_state_w'], rep['bend']['exact_state_w_vs_fd']))
    print('  DLL SHA-256 %s' % dll_sha)

    res = {
        'engine_dir': ENGINE_DIR.replace('\\', '/'), 'engine_commit': ENGINE_COMMIT, 'dll_sha256': dll_sha, 'build_info': bi,
        'threads': {'default': threads, 'single_process': threads1},
        'scene': {'nodes': int(len(V)), 'triangles': int(len(F)), 'hinges': int(len(H)), 'cells': M, 'spacing_cm': S,
                  'seed': SEED, 'random_vectors': N_RAND, 'fd_step_cm': FD_STEP,
                  'materials': {n: m for n, m in MATERIALS}, 'k_stretch': K_STRETCH, 'k_bend': K_BEND},
        'L0_0': {'build_ok': pre_ok, 'engine_status_clean': st.strip() == ''},
        'L0_1': {'limit': LIM_PRODUCT, 'rel': l01},
        'L0_2': {'limit': LIM_STRETCH, 'f_rel': l02},
        'L0_3': {'limit': LIM_STRETCH, 'J_vs_J_T_rel': l03, 'components_C_pos': nPos, 'components_C_nonpos': nNeg, 'min_branch': MIN_BRANCH},
        'L0_4': l04,
        'report': rep,
        'gates_failed': fail,
    }
    with open(os.path.join(ROOT, 'results', 'gate_link_l0.json'), 'w', encoding='utf-8') as fh:
        json.dump(res, fh, indent=2, ensure_ascii=False)
    print('\n엔진 연결 L0 %s (실패 %d) → results/gate_link_l0.json' % ('PASS' if fail == 0 else 'FAIL', fail))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
