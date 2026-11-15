"""판정 밖 진단 (docs/GATE_LINK_L4.md 「(가) 판정 뒤 진단」): L4 (가) 위반 넷(장면 P · BE_U 1e4)은 왜 생겼는가.

A. 엔진 기본 풀이(수렴 PCG, 상대 1e-6 · 상한 300). 같은 회를 상한 2000 으로도 돌려, 상한만 넘은 것인지 계가 풀리지 않는 것인지
   본다(1/200 · 1/100 · 1/50). 20 회 고정(1/200)은 풀이가 끝난 잔차 비(전처리 잔차 제곱의 비 δ_끝 / δ0)를 본다.
B. 물리 속력. 작은 h(1/1600, 수렴 PCG)에서 노드 최대 속력을 본다. 이 장면은 한 변에 매달린 30 cm 판이 수평에서 내려
   흔들리는 것이라, 굽지 않는 판이면 끝 속력이 √(3 g L) 다(경첩에 매달린 막대 — 무게 중심이 L/2 내려가고 관성 모멘트 m L² / 3).
C. 엔진 AMG 전처리(MGCG, 적응 끔 — 환경 변수 CFMS_MGCG_ADAPTIVE=0)로 P · 1e4 의 세 h · 두 PCG 와 엔진 기본 설정
   (P · BE_U 1 · 1/200 · 수렴)을 돌린다. 적응을 켜면(엔진 기본) 반복 12 회가 넘는 풀이에서 Jacobi 로 물러서 효과가 없다.

사용:  python tests/diag_link_l4_a.py <브랜치 cfmsDrape.dll>     (리포 루트에서). 결과 results/diag_link_l4_a.json.
회는 tests/link_l4.py 와 같은 장면 · 측정 · 단일 스레드 프로세스다. 판정에 쓰지 않는다.
"""
import sys
sys.dont_write_bytecode = True

import concurrent.futures  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import link_l0 as L0  # noqa: E402
import link_l4 as L4  # noqa: E402

OUT = os.path.join(L4.ROOT, 'results', 'diag_link_l4_a.json')
PLATE_L = 30.0                              # 장면 P 의 판 길이(cm)
MODES = ('off', 'gn_state')
# (절, 장면, BE_U, 단계/초, PCG, 전처리) — PCG: conv = 엔진 기본(상한 300), conv2000 = 상한만 2000, fixed20 = 옛 경로
JOBS = ([('A', 'P', 1e4, 200, 'conv', 'jacobi')]
        + [('A', 'P', 1e4, s, 'conv2000', 'jacobi') for s in (200, 100, 50)]
        + [('A', 'P', 1e4, 200, 'fixed20', 'jacobi')]
        + [('B', 'P', 1e4, 1600, 'conv', 'jacobi')]
        + [('C', 'P', 1e4, s, p, 'mgcg') for s in (200, 100, 50) for p in ('conv', 'fixed20')]
        + [('C', 'P', 1.0, 200, 'conv', 'mgcg')])


def child(dll, sc, bend, sps, pcg, prec, out):
    Drape = L0.load_drape()
    V, F, uv, pins = L4.scene(sc)
    strain = L4.strain_fn(F, uv)
    nsteps = int(round(L4.T_END * sps))
    d = Drape(dll)
    info = {'mode': os.environ.get('CFMS_CPU_SXJ_BEND', '') or 'off', 'build_info': d.build_info}
    stop, steps, iters, first, ratios, vmax, smax, wall = '', 0, [], [], [], 0.0, 0.0, 0.0
    try:
        d.set_mesh(V, F)
        d.set_mesh_2d(uv)
        d.set_material(BE_U=float(bend))
        d.set_gravity((0.0, -L4.GRAVITY, 0.0))
        d.set_timestep(1.0 / sps)
        if pcg == 'fixed20':
            d.set_dynamic_pcg(False)
        elif pcg == 'conv2000':
            d.set_dynamic_pcg(True, 1e-6, 2000)
        if prec == 'mgcg':
            d.set_mgcg(True, 0)
        for p in pins:
            d.pin_vertex(p)
        for s in range(1, nsteps + 1):
            t0 = time.perf_counter()
            try:
                d.step(1)
            except Exception as e:   # StepDiverged 등
                stop = 'diverged: %s' % str(e)[:120]
                break
            wall += time.perf_counter() - t0
            steps = s
            X = np.array(d.positions(), dtype=float)
            if not np.isfinite(X).all():
                stop = 'nonfinite'
                break
            hist = d.get_pcg_solve_history()
            n = sum(int(r['iterations']) for r in hist)
            iters.append(n)
            if len(first) < 6:
                first.append(n)
            last = hist[-1] if hist else None
            if last is not None and last['delta0'] > 0:
                ratios.append(last['delta_final'] / last['delta0'])
            if any(r['exit_reason'] == 'numerical_failure' for r in hist):
                stop = 'pcg_numerical_failure'
                break
            if any(int(r['attempt']) > 0 for r in hist):
                stop = 'pcg_retry'
                break
            v = float(np.linalg.norm(np.array(d.velocities(), dtype=float), axis=1).max())
            vmax = max(vmax, v)
            if not (v <= L4.VMAX):
                stop = 'speed'
                break
            if s % int(round(L4.SAMPLE_DT * sps)) == 0:
                e = strain(X)
                smax = max(smax, e)
                if not (e <= L4.STRAIN_MAX):
                    stop = 'strain'
                    break
    finally:
        d.close()
    info.update({'usable': stop == '', 'stop': stop, 't_stop': (steps / sps) if stop else None, 'steps': steps,
                 'pcg_iters_mean': float(np.mean(iters)) if iters else 0.0, 'pcg_iters_max': int(max(iters)) if iters else 0,
                 'pcg_iters_first': first, 'exit_ratio_median': float(np.median(ratios)) if ratios else None,
                 'vmax_cm_s': vmax, 'strain_max': smax, 'ms_per_step': 1e3 * wall / max(steps, 1)})
    with open(out, 'w', encoding='utf-8') as fh:
        json.dump(info, fh)


def spawn(dll, mode, job, td):
    sec, sc, bend, sps, pcg, prec = job
    out = os.path.join(td, '%s_%s_%s_%g_%d_%s_%s.json' % (sec, mode, sc, bend, sps, pcg, prec))
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1')
    env.pop('CFMS_MGCG_ADAPTIVE', None)
    # 끔도 'off' 로 명시한다 — 브랜치의 기본(환경 변수 없음)이 gn_state 가 된 뒤(2026-11-13)로는 없음 ≠ 끔이다.
    env['CFMS_CPU_SXJ_BEND'] = mode
    # 인장 · 전단은 엔진 손 커널로 고정한다 — 브랜치의 기본이 on 이 된 뒤(2026-11-13)로는 없음 ≠ 끔이고, 이 판정은 굽힘만 본다.
    env['CFMS_CPU_SXJ_MEMBRANE'] = 'off'
    if prec == 'mgcg':
        env['CFMS_MGCG_ADAPTIVE'] = '0'
    r = subprocess.run([sys.executable, os.path.abspath(__file__), '--run', dll, sc, repr(float(bend)), str(sps), pcg,
                        prec, out], env=env, cwd=td, capture_output=True, text=True, encoding='utf-8', errors='replace')
    if r.returncode != 0 or not os.path.isfile(out):
        raise RuntimeError('회 실패 (%s, %s): %s' % (mode, job, r.stderr[-2000:]))
    with open(out, encoding='utf-8') as fh:
        return json.load(fh)


def main():
    dll = os.path.abspath(sys.argv[1])
    td = tempfile.mkdtemp(prefix='cfmsSXJ_L4a_')
    t0 = time.perf_counter()
    rows = []
    try:
        todo = [(m, j) for j in JOBS for m in MODES]
        with concurrent.futures.ThreadPoolExecutor(max_workers=L4.PARALLEL) as ex:
            futs = {ex.submit(spawn, dll, m, j, td): (m, j) for m, j in todo}
            got = {}
            for fu in concurrent.futures.as_completed(futs):
                got[futs[fu]] = fu.result()
        for m, j in todo:
            sec, sc, bend, sps, pcg, prec = j
            info = got[(m, j)]
            rows.append({'section': sec, 'mode': m, 'scene': sc, 'BE_U': bend, 'sps': sps, 'pcg': pcg, 'prec': prec,
                         **{k: info[k] for k in ('usable', 'stop', 't_stop', 'steps', 'pcg_iters_mean', 'pcg_iters_max',
                                                 'pcg_iters_first', 'exit_ratio_median', 'vmax_cm_s', 'strain_max',
                                                 'ms_per_step')},
                         'git_commit': info['build_info'].get('git_commit'), 'dirty': info['build_info'].get('dirty')})
    finally:
        shutil.rmtree(td, ignore_errors=True)
    res = {'rigid_plate_tip_speed_cm_s': math.sqrt(3.0 * L4.GRAVITY * PLATE_L), 'runs': rows,
           'seconds': time.perf_counter() - t0}
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
        fh.write('\n')

    print('판정 밖 진단 — L4 (가) (장면 P). 굽지 않는 판의 끝 속력 √(3 g L) = %.0f cm/s' % res['rigid_plate_tip_speed_cm_s'])
    for r in rows:
        print('  %s %-8s BE_U %-5g 1/%-4d %-8s %-6s %-16s PCG 평균 %6.1f 최대 %4d 처음 %-26s 끝 잔차 비 %-8s 최대 속력 %4.0f cm/s'
              ' 변형률 %.3f  %5.0f ms/단계'
              % (r['section'], r['mode'], r['BE_U'], r['sps'], r['pcg'], r['prec'],
                 'U' if r['usable'] else ('%s t=%.2f' % (r['stop'][:9], r['t_stop'])),
                 r['pcg_iters_mean'], r['pcg_iters_max'], ','.join(str(x) for x in r['pcg_iters_first']),
                 '-' if r['exit_ratio_median'] is None else '%.1e' % r['exit_ratio_median'], r['vmax_cm_s'],
                 r['strain_max'], r['ms_per_step']))
    print('결과 %s (%.0f s)' % (os.path.relpath(OUT, L4.ROOT), res['seconds']))
    return 0


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--run':
        a = sys.argv[2:]
        child(a[0], a[1], float(a[2]), int(a[3]), a[4], a[5], a[6])
        sys.exit(0)
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main())
