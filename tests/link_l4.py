"""엔진 연결 L4 (굽힘) (docs/GATE_LINK_L4.md): 엔진 브랜치 DLL 로 천을 끝까지 적분해 엔진 Oh 굽힘과 생성 Bridson GN 을 견준다.

사용:  python tests/link_l4.py <브랜치 cfmsDrape.dll> [--out <결과 json>]
       --out 을 주면 그 파일에 쓰고(2 차 (다) 는 results/gate_link_l4_2.json), 회마다 결정적 값을 1 차
       results/gate_link_l4.json 과 견준다 — 끔 회가 하나라도 다르면 전제 (e) 실패(판정문 「2 차 (다)」).
       (내부) python tests/link_l4.py --run <dll> <장면> <BE_U> <단계/초> <conv|fixed20> <out 접두> <samples 0|1>
              — 한 회를 다른 프로세스에서 돌린다. 굽힘 모드는 그 프로세스의 환경 변수 CFMS_CPU_SXJ_BEND 로 준다.
엔진의 Python 포장 cfmsdrape 는 tests/link_l0.py 와 같이 CFMS_DRAPE_DIR(없으면 ..\\cfmsDrape_dev)의 것을 읽는다. 브랜치 워크트리는
CFMS_DRAPE_LINK_DIR(없으면 ..\\cfmsDrape_dev_sxjacobian-link). editbin 은 PATH 에서 찾고, 없으면 VS 2022 Community 의 vcvars64.bat 를
불러 쓴다. 리포 루트에서 실행. 결과는 results/gate_link_l4.json. 통과 0, 실패 1, 준비 안 됨 2.
"""
import sys
sys.dont_write_bytecode = True

import concurrent.futures  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import statistics  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import link_l0 as L0  # noqa: E402
import link_l2 as L2  # noqa: E402

ROOT = L0.ROOT
OUT_JSON = os.path.join(ROOT, 'results', 'gate_link_l4.json')
FIRST_JSON = OUT_JSON                   # 1 차 결과(--out 으로 다른 곳에 쓸 때 견줄 상대)
FP_KEYS = ('usable', 'stop', 'steps', 'pcg_iters_sum', 'strain_max', 'vmax_cm_s')   # 회의 결정적 값

# 판정문 「장면」 · 「변형」 · 「측정」 · 「판정 기준」에 고정한 값
SCENES = ('P', 'R')
BENDS = (1.0, 1e2, 1e4)                 # BE_U, erg
SPS = (400, 200, 100, 50)               # 단계 / 초 (h = 1 / SPS)
PCGS = ('conv', 'fixed20')
MODES = ('off', 'gn_rest', 'gn_state')
GN_MODES = ('gn_rest', 'gn_state')
T_END = 2.0
SAMPLE_DT = 0.02
N_SAMPLES = 100
GRAVITY = 980.0
VMAX = 5000.0                           # cm/s
STRAIN_MAX = 0.25
RUN_TIMEOUT = 900.0                     # s, 한 회의 벽시계
TIME_RATIO = 1.05
REF_TIME = ('P', 1.0, 200, 'conv')      # (다) 기준 짝 = 엔진 기본 설정
DET = ('R', 1e2, 200, 'conv')           # L4-0 (d)
DET_MODES = ('off', 'gn_rest')
BASES = ['0x180000000', '0x180010000', '0x180020000', '0x180030000',
         '0x190000000', '0x1A0000000', '0x1E0000000', '0x280000000']
VCVARS = r'C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat'
PARALLEL = max(1, min(16, (os.cpu_count() or 2) // 2))


def scene(name):
    if name == 'P':
        n = 81
        xs = np.linspace(-15.0, 15.0, n)
        gx, gz = np.meshgrid(xs, xs)
        V = np.column_stack([gx.ravel(), np.full(gx.size, 20.0), gz.ravel()])
        pins = list(range(n))                          # z = −15 의 한 변 전부 (perf_bench)
    else:
        n = 21
        xs = np.linspace(0.0, 50.0, n)
        gx, gz = np.meshgrid(xs, xs)
        V = np.column_stack([gx.ravel(), np.zeros(gx.size), gz.ravel()])
        pins = [n * (n - 1), n * n - 1]                # z = 50 의 한 변의 두 끝
    F = []
    for j in range(n - 1):
        for i in range(n - 1):
            a = j * n + i
            F += [[a, a + 1, (j + 1) * n + i + 1], [a, (j + 1) * n + i + 1, (j + 1) * n + i]]
    uv = np.column_stack([gx.ravel(), gz.ravel()])
    return V, np.array(F, np.int32), uv, pins


def strain_fn(F, uv):
    Dm = np.stack([uv[F[:, 1]] - uv[F[:, 0]], uv[F[:, 2]] - uv[F[:, 0]]], axis=2)   # (T, 2, 2), 열 = 변
    Dinv = np.linalg.inv(Dm)

    def strain(X):
        Ds = np.stack([X[F[:, 1]] - X[F[:, 0]], X[F[:, 2]] - X[F[:, 0]]], axis=2)    # (T, 3, 2)
        W = np.einsum('tij,tjk->tik', Ds, Dinv)                                      # 열 = w_u, w_v
        return float(np.abs(np.linalg.norm(W, axis=1) - 1.0).max())
    return strain


def child_run(dll, sc, bend, sps, pcg, prefix, save_samples, pcg_cap=None, timeout=RUN_TIMEOUT):
    # pcg_cap · timeout 은 tests/link_l4r.py 의 기준 회만 준다(수렴 PCG 의 상한만 바꿈). L4 의 회는 기본값 그대로다.
    Drape = L0.load_drape()
    V, F, uv, pins = scene(sc)
    strain = strain_fn(F, uv)
    nsteps, every = int(round(T_END * sps)), int(round(SAMPLE_DT * sps))
    d = Drape(dll)
    info = {'base': hex(int(d._lib._handle)), 'build_info': d.build_info, 'threads': int(d.openmp_threads),
            'mode': os.environ.get('CFMS_CPU_SXJ_BEND', '') or 'off'}
    samples = np.full((N_SAMPLES, len(V), 3), np.nan) if save_samples else None
    stop, t_stop, step_time, steps_done = '', None, 0.0, 0
    iters, reasons, vmax, smax = [], {}, 0.0, 0.0
    try:
        d.set_mesh(V, F)
        d.set_mesh_2d(uv)
        d.set_material(BE_U=float(bend))
        d.set_gravity((0.0, -GRAVITY, 0.0))
        d.set_timestep(1.0 / sps)
        if pcg == 'fixed20':
            d.set_dynamic_pcg(False)
        elif pcg_cap is not None:
            d.set_dynamic_pcg(True, 1e-6, int(pcg_cap))
        for p in pins:
            d.pin_vertex(p)
        wall0 = time.perf_counter()
        for s in range(1, nsteps + 1):
            t0 = time.perf_counter()
            try:
                d.step(1)
            except Exception as e:   # StepDiverged 등
                stop = 'diverged: %s: %s' % (type(e).__name__, str(e)[:200])
                break
            step_time += time.perf_counter() - t0
            steps_done = s
            X = np.array(d.positions(), dtype=float)
            if not np.isfinite(X).all():
                stop = 'nonfinite'
                break
            hist = d.get_pcg_solve_history()
            failed, retried = False, False
            for r in hist:
                iters.append(int(r['iterations']))
                reasons[r['exit_reason']] = reasons.get(r['exit_reason'], 0) + 1
                failed |= r['exit_reason'] == 'numerical_failure'
                retried |= int(r['attempt']) > 0
            if failed:
                stop = 'pcg_numerical_failure'
                break
            if retried:
                stop = 'pcg_retry'
                break
            v = float(np.linalg.norm(np.array(d.velocities(), dtype=float), axis=1).max())
            vmax = max(vmax, v)
            if not (v <= VMAX):
                stop = 'speed'
                break
            if s % every == 0:
                e = strain(X)
                smax = max(smax, e)
                if save_samples:
                    samples[s // every - 1] = X
                if not (e <= STRAIN_MAX):
                    stop = 'strain'
                    break
            if time.perf_counter() - wall0 > timeout:
                stop = 'timeout'
                break
        if stop:
            t_stop = steps_done / sps
    finally:
        d.close()
    info.update({'usable': stop == '', 'stop': stop, 't_stop': t_stop, 'steps': steps_done, 'steps_planned': nsteps,
                 'step_seconds': step_time, 'ms_per_sim_second': 1e3 * step_time / T_END,
                 'pcg_iters_sum': int(sum(iters)), 'pcg_iters_mean': float(np.mean(iters)) if iters else 0.0,
                 'pcg_iters_max': int(max(iters)) if iters else 0, 'pcg_exit_reasons': reasons,
                 'vmax_cm_s': vmax, 'strain_max': smax})
    if save_samples:
        np.save(prefix + '.npy', samples)
    with open(prefix + '.json', 'w', encoding='utf-8') as fh:
        json.dump(info, fh)
    return 0


def spawn(dll, mode, cfg, prefix, save_samples):
    sc, bend, sps, pcg = cfg
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1')
    # 끔도 'off' 로 명시한다 — 브랜치의 기본(환경 변수 없음)이 gn_state 가 된 뒤(2026-11-13)로는 없음 ≠ 끔이다.
    env['CFMS_CPU_SXJ_BEND'] = mode
    # 인장 · 전단은 엔진 손 커널로 고정한다 — 브랜치의 기본이 on 이 된 뒤(2026-11-13)로는 없음 ≠ 끔이고, 이 판정은 굽힘만 본다.
    env['CFMS_CPU_SXJ_MEMBRANE'] = 'off'
    args = [sys.executable, os.path.abspath(__file__), '--run', dll, sc, repr(float(bend)), str(sps), pcg, prefix,
            '1' if save_samples else '0']
    try:
        r = subprocess.run(args, env=env, cwd=os.path.dirname(prefix), capture_output=True, text=True,
                           encoding='utf-8', errors='replace', timeout=RUN_TIMEOUT + 120)
    except subprocess.TimeoutExpired:
        return {'mode': mode, 'usable': False, 'stop': 'timeout (killed)', 'steps': None}
    if r.returncode != 0 or not os.path.isfile(prefix + '.json'):
        raise RuntimeError('회 실패 (%s, %s): %s' % (mode, cfg, r.stderr[-2000:]))
    with open(prefix + '.json', encoding='utf-8') as fh:
        info = json.load(fh)
    info['stderr'] = r.stderr.strip()[-500:]
    return info


def key(mode, cfg):
    sc, bend, sps, pcg = cfg
    return '%s|%s|%g|%d|%s' % (mode, sc, bend, sps, pcg)


def dist(pa, pb):
    A, B = np.load(pa), np.load(pb)
    return float(np.sqrt(((A - B) ** 2).sum(axis=2).mean(axis=1)).max())


def find_editbin():
    if shutil.which('editbin'):
        return None
    if os.path.isfile(VCVARS):
        return VCVARS
    raise RuntimeError('editbin 도 vcvars64.bat 도 없다')


def make_rebased(dll, td):
    vc = find_editbin()
    out = []
    for i, base in enumerate(BASES):
        dst_dir = os.path.join(td, 'base%d' % i)
        os.makedirs(dst_dir)
        dst = os.path.join(dst_dir, 'cfmsDrape.dll')
        shutil.copyfile(dll, dst)
        bat = os.path.join(dst_dir, 'rebase.bat')
        lines = ['@echo off']                       # cmd 는 \" 를 모르므로 .bat 로 부른다
        if vc:
            lines.append('call "%s" >nul' % vc)
        lines.append('editbin /NOLOGO /REBASE:BASE=%s /DYNAMICBASE:NO "%s"' % (base, dst))
        lines.append('exit /b %ERRORLEVEL%')
        with open(bat, 'w', encoding='ascii', newline='\r\n') as fh:
            fh.write('\n'.join(lines) + '\n')
        r = subprocess.run(['cmd', '/d', '/c', bat], capture_output=True, text=True, encoding='mbcs', errors='replace')
        if r.returncode != 0:
            raise RuntimeError('editbin 실패 (%s): %s %s' % (base, r.stdout, r.stderr))
        out.append(dst)
    return out


def gate(fail, ok, label, text, limit):
    print('  %-10s %s  %s  (기준 %s)' % (label, 'PASS' if ok else 'FAIL', text, limit))
    return fail + (0 if ok else 1)


def main():
    if len(sys.argv) == 9 and sys.argv[1] == '--run':
        a = sys.argv
        return child_run(a[2], a[3], float(a[4]), int(a[5]), a[6], a[7], a[8] == '1')
    out_json = OUT_JSON
    if len(sys.argv) == 4 and sys.argv[2] == '--out':
        out_json = os.path.abspath(sys.argv[3])
    elif len(sys.argv) != 2:
        print(__doc__)
        return 2
    dll = os.path.abspath(sys.argv[1])
    first = None
    if os.path.normcase(out_json) != os.path.normcase(OUT_JSON):
        if not os.path.isfile(FIRST_JSON):
            print('없는 것: %s (견줄 1 차 결과)' % FIRST_JSON)
            return 2
        with open(FIRST_JSON, encoding='utf-8') as fh:
            first = json.load(fh)
    missing = [p for p in (dll, L2.HEADER_GEN) if not os.path.isfile(p)]
    missing += [p for p in (os.path.join(L0.ENGINE_DIR, 'drape'), os.path.join(L2.LINK_DIR, 'drape')) if not os.path.isdir(p)]
    if missing:
        print('없는 것: %s' % ', '.join(missing))
        return 2

    configs = [(sc, b, s, p) for sc in SCENES for b in BENDS for s in SPS for p in PCGS]
    print('엔진 연결 L4 (굽힘) — 장면 %s, 회 %d (+ 결정성 %d), 동시 프로세스 %d\n  브랜치 DLL %s (워크트리 %s)'
          % ('·'.join(SCENES), len(configs) * len(MODES), len(DET_MODES), PARALLEL, dll, L2.LINK_DIR))
    t_all = time.perf_counter()
    td = tempfile.mkdtemp(prefix='cfmsSXJ_L4_')
    try:
        # 적분 회 — 오래 걸리는 P 를 먼저 넣는다
        jobs = [(m, c, os.path.join(td, 'run_%03d' % i)) for i, (m, c) in
                enumerate((m, c) for c in configs for m in MODES)]
        jobs += [(m, DET, os.path.join(td, 'det_%s' % m)) for m in DET_MODES]
        jobs.sort(key=lambda j: (j[1][0] != 'P', -j[1][2]))
        runs = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=PARALLEL) as ex:
            futs = {ex.submit(spawn, dll, m, c, pre, True): (m, c, pre) for m, c, pre in jobs}
            for n_done, fu in enumerate(concurrent.futures.as_completed(futs), 1):
                m, c, pre = futs[fu]
                info = fu.result()
                info['prefix'] = pre
                runs[('det:' if os.path.basename(pre).startswith('det_') else '') + key(m, c)] = info
                if n_done % 24 == 0 or n_done == len(jobs):
                    print('  적분 회 %d / %d (%.0f s)' % (n_done, len(jobs), time.perf_counter() - t_all))
        t_integ = time.perf_counter() - t_all

        # 시간 회 — 고정 주소 8 곳, 혼자, 주소마다 세 모드를 번갈아
        copies = make_rebased(dll, td)
        timing = {m: [] for m in MODES}
        for i, cp in enumerate(copies):
            for m in MODES:
                info = spawn(cp, m, REF_TIME, os.path.join(os.path.dirname(cp), 't_%s' % m), False)
                info['base_expected'] = BASES[i]
                timing[m].append(info)
        t_time = time.perf_counter() - t_all - t_integ

        # 전제
        st = L2.branch_state()
        sha_gen = L2.sha256_file(L2.HEADER_GEN)
        sha_eng = L2.sha256_file(L2.HEADER_ENGINE) if os.path.isfile(L2.HEADER_ENGINE) else ''
        bis = [r['build_info'] for r in runs.values() if 'build_info' in r]
        bi = bis[0] if bis else {}
        same_bi = all(b.get('git_commit') == bi.get('git_commit') and b.get('dirty') == bi.get('dirty') for b in bis)
        fail0 = 0
        print(' L4-0 전제')
        a_ok = (same_bi and st['head'] != '' and str(bi.get('git_commit', '')).startswith(st['head'][:7])
                and bi.get('dirty') is False and st['descends_from_base'] and not st['outside_allowed'] and not st['uncommitted'])
        fail0 = gate(fail0, a_ok, '(a)', 'DLL git_commit %s · 워크트리 HEAD %s · dirty %s · 미커밋 %s · 허용 밖 %s'
                     % (bi.get('git_commit'), st['head'], bi.get('dirty'), st['uncommitted'], st['outside_allowed'] or '없음'),
                     'HEAD 일치 · dirty false · 허용 목록 안')
        fail0 = gate(fail0, sha_gen == sha_eng, '(b)', '헤더 사본 %s' % ('같음' if sha_gen == sha_eng else '다름'), 'SHA-256 같음')
        c_bad = [(t['base_expected'], t.get('base')) for m in MODES for t in timing[m]
                 if t.get('base') is None or int(t['base'], 16) != int(t['base_expected'], 16)
                 or t.get('build_info', {}).get('git_commit') != bi.get('git_commit')
                 or t.get('build_info', {}).get('dirty') != bi.get('dirty')]
        fail0 = gate(fail0, not c_bad, '(c)', '사본 8 개 %s' % ('모두 제 주소 · 같은 빌드' if not c_bad else '어긋남 %s' % c_bad),
                     '정한 주소에 올라감 · git_commit · dirty 같음')
        det = {}
        for m in DET_MODES:
            ra, rb = runs[key(m, DET)], runs['det:' + key(m, DET)]
            same = (ra['usable'] and rb['usable']
                    and np.load(ra['prefix'] + '.npy').tobytes() == np.load(rb['prefix'] + '.npy').tobytes())
            det[m] = bool(same)
        fail0 = gate(fail0, all(det.values()), '(d)', '두 번 돌린 표본 위치 %s'
                     % ' · '.join('%s %s' % (m, '같음' if det[m] else '다름') for m in DET_MODES), '비트까지 같음')
        vs_first = None
        if first is not None:   # 2 차: 회마다 결정적 값을 1 차와 견준다
            fr = first['runs']

            def mode_of(k):
                return k.split('|')[0].split(':')[-1]
            same = {k: (k in fr and all(fr[k].get(f) == r.get(f) for f in FP_KEYS)) for k, r in runs.items()}
            vs_first = {'first_json': os.path.relpath(FIRST_JSON, ROOT), 'keys': list(FP_KEYS),
                        'same': {m: sum(1 for k, v in same.items() if v and mode_of(k) == m) for m in MODES},
                        'runs': {m: sum(1 for k in same if mode_of(k) == m) for m in MODES},
                        'different': sorted(k for k, v in same.items() if not v)}
            off_bad = [k for k, v in same.items() if not v and mode_of(k) == 'off']
            fail0 = gate(fail0, not off_bad, '(e)', '끔 회 결정적 값이 1 차와 %s' % (
                '모두 같음 (%d 회)' % vs_first['runs']['off'] if not off_bad else '다름 %d 회 %s' % (len(off_bad), off_bad[:4])),
                '끔 회 모두 같음')

        # 짝 · 거리
        def run(m, c):
            return runs[key(m, c)]
        pairs = []
        for c in configs:
            sc, b, s, p = c
            oh = run('off', c)
            half = (sc, b, 2 * s, p) if 2 * s in SPS else None
            for m in GN_MODES:
                g = run(m, c)
                rec = {'mode': m, 'scene': sc, 'BE_U': b, 'sps': s, 'pcg': p,
                       'off_usable': oh['usable'], 'gn_usable': g['usable'], 'violation': bool(oh['usable'] and not g['usable'])}
                if oh['usable'] and g['usable']:
                    rec['D_gn_off'] = dist(g['prefix'] + '.npy', oh['prefix'] + '.npy')
                if half is not None:
                    oh2 = run('off', half)
                    rec['off_half_usable'] = oh2['usable']
                    if oh['usable'] and oh2['usable']:
                        rec['D_off_half'] = dist(oh['prefix'] + '.npy', oh2['prefix'] + '.npy')
                    if 'D_gn_off' in rec and 'D_off_half' in rec:
                        rec['b_counted'] = True
                        rec['b_ok'] = bool(rec['D_gn_off'] <= rec['D_off_half'])
                pairs.append(rec)
        rest_state = {}
        for c in configs:
            ga, gb = run('gn_rest', c), run('gn_state', c)
            if ga['usable'] and gb['usable']:
                rest_state[key('rest_vs_state', c)] = dist(ga['prefix'] + '.npy', gb['prefix'] + '.npy')
    finally:
        shutil.rmtree(td, ignore_errors=True)

    # 판정
    print(' 모드별 판정 (가) 안정 · (나) 차이 · (다) 시간')
    def t_ok_all(ts):
        return all(t['usable'] for t in ts)
    tmin = {m: (min(t['ms_per_sim_second'] for t in timing[m]) if t_ok_all(timing[m]) else None) for m in MODES}
    tmed = {m: (statistics.median(t['ms_per_sim_second'] for t in timing[m]) if t_ok_all(timing[m]) else None) for m in MODES}
    t_usable = {m: all(t['usable'] for t in timing[m]) for m in MODES}
    verdict, fail = {}, 0
    for m in GN_MODES:
        pm = [r for r in pairs if r['mode'] == m]
        viol = [r for r in pm if r['violation']]
        counted = [r for r in pm if r.get('b_counted')]
        b_bad = [r for r in counted if not r['b_ok']]
        a_ok = not viol
        b_ok = bool(counted) and not b_bad
        if t_usable['off'] and t_usable[m]:
            ratio = tmin[m] / tmin['off']
            c_ok, c_text = ratio <= TIME_RATIO, 'min %.1f / %.1f = %.3f' % (tmin[m], tmin['off'], ratio)
        elif t_usable[m]:
            ratio, c_ok, c_text = None, True, '끔 이 시간 회에서 쓸 수 없고 %s 는 쓸 수 있음' % m
        else:
            ratio, c_ok, c_text = None, False, ('판정 불가(둘 다 쓸 수 없음)' if not t_usable['off'] else '%s 가 시간 회에서 쓸 수 없음' % m)
        print('  %s' % m)
        fail = gate(fail, a_ok, '(가)', '위반 %d / 짝 %d' % (len(viol), len(pm)), '위반 0')
        worst = max((r['D_gn_off'] / r['D_off_half'] for r in counted if r['D_off_half'] > 0), default=None)
        fail = gate(fail, b_ok, '(나)', '센 짝 %d, 넘음 %d · 못 넘음 %d, 최대 D(M,Oh)/D(Oh,Oh½) %s'
                    % (len(counted), len(counted) - len(b_bad), len(b_bad), '-' if worst is None else '%.3g' % worst),
                    '센 짝 모두 ≤, 하나 이상')
        fail = gate(fail, c_ok, '(다)', c_text, '≤ %.2f' % TIME_RATIO)
        verdict[m] = {'a_ok': a_ok, 'violations': len(viol), 'pairs': len(pm), 'b_ok': b_ok, 'b_counted': len(counted),
                      'b_failed': len(b_bad), 'b_worst_ratio': worst, 'c_ok': bool(c_ok), 'time_ratio': ratio,
                      'usable_mode': bool(a_ok and b_ok and c_ok)}
    if verdict['gn_rest']['usable_mode']:
        rec = 'gn_rest'
    elif verdict['gn_state']['usable_mode']:
        rec = 'gn_state'
    else:
        rec = 'off'
    rec_text = {'gn_rest': '엔진 굽힘을 gn_rest 로 권함', 'gn_state': '엔진 굽힘을 gn_state 로 권함',
                'off': '엔진 Oh 를 그대로 둠'}[rec]

    # 보고
    usable_counts = {}
    for m in MODES:
        for sc in SCENES:
            for p in PCGS:
                usable_counts['%s|%s|%s' % (m, sc, p)] = sum(1 for c in configs if c[0] == sc and c[3] == p and run(m, c)['usable'])
    print(' 보고 — 쓸 수 있는 회 (장면 · PCG 마다 BE_U 3 × h 4 = 12 회 중)')
    for sc in SCENES:
        for p in PCGS:
            print('  %s · %-7s ' % (sc, p) + ' · '.join('%s %d' % (m, usable_counts['%s|%s|%s' % (m, sc, p)]) for m in MODES))
    stops = {}
    for k, r in runs.items():
        if not k.startswith('det:') and not r['usable']:
            stops.setdefault(r['stop'].split(':')[0], []).append(k)
    print('  멈춘 까닭: ' + (' · '.join('%s %d' % (s, len(v)) for s, v in sorted(stops.items())) or '없음'))
    print('  시간 회(P · BE_U 1 · 1/200 · 수렴) ms / 모의 1 초 — 최소 · 중앙: ' + ' · '.join(
        '%s %s · %s' % (m, '-' if tmin[m] is None else '%.1f' % tmin[m], '-' if tmed[m] is None else '%.1f' % tmed[m])
        for m in MODES))
    if vs_first is not None:
        print('  1 차와 결정적 값이 같은 회: ' + ' · '.join('%s %d / %d' % (m, vs_first['same'][m], vs_first['runs'][m])
                                                     for m in MODES))
    print('  실행 시간: 적분 %.0f s · 시간 %.0f s' % (t_integ, t_time))

    status = 'PASS' if (fail0 == 0 and rec != 'off') else ('NOT_JUDGED' if fail0 else 'FAIL')
    res = {'gate': 'link_L4_bend', 'status': status, 'recommendation': rec,
           'L4_0': {'branch': st, 'build_info': bi, 'header_sha256': {'generated': sha_gen, 'engine_copy': sha_eng},
                    'rebased_bases_ok': not c_bad, 'determinism': det, 'fail': fail0},
           'verdict': verdict, 'pairs': pairs, 'usable_counts': usable_counts, 'D_rest_vs_state': rest_state,
           'timing': {'config': list(REF_TIME), 'bases': BASES, 'min_ms_per_sim_second': tmin, 'median_ms_per_sim_second': tmed,
                      'per_base': {m: [{k: t.get(k) for k in ('base', 'base_expected', 'usable', 'stop', 'ms_per_sim_second',
                                                              'pcg_iters_sum', 'steps')} for t in timing[m]] for m in MODES}},
           'runs': {k: {kk: vv for kk, vv in r.items() if kk not in ('prefix', 'build_info')} for k, r in sorted(runs.items())},
           'config': {'scenes': SCENES, 'BE_U': BENDS, 'sps': SPS, 'pcg': PCGS, 'T': T_END, 'sample_dt': SAMPLE_DT,
                      'vmax_cm_s': VMAX, 'strain_max': STRAIN_MAX, 'run_timeout_s': RUN_TIMEOUT, 'time_ratio': TIME_RATIO,
                      'parallel': PARALLEL},
           'seconds': {'integration': t_integ, 'timing': t_time}}
    if vs_first is not None:
        res['vs_first'] = vs_first
    with open(out_json, 'w', encoding='utf-8') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
    if fail0:
        print('\n엔진 연결 L4 (굽힘) 전제 실패 %d — 판정하지 않음 → %s' % (fail0, os.path.relpath(out_json, ROOT)))
        return 2
    print('\n엔진 연결 L4 (굽힘) %s — 권고: %s → %s' % (status, rec_text, os.path.relpath(out_json, ROOT)))
    return 0 if status == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
