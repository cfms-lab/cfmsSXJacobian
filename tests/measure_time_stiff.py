"""논문용 측정 2 (docs/GATE_JACOBIAN_DAMPING.md 「논문용 측정 2」, 판정 아님): 깨끗한 설정의 낙하에서 두 Jacobian 의 계산 시간 — 굽힘 강성 1e3 – 1e5 erg.

사용:  python tests/measure_time_stiff.py <브랜치 cfmsDrape.dll>     (리포 루트에서, run_gate 밖 — 약 25 분, **혼자 돌린다**)
       (내부) python tests/measure_time_stiff.py --child <인자 json>
회는 tests/gate_jacobian_damping.py 의 깨끗한 설정(make)을 그대로 쓴다. 시간은 이 리포의 규칙대로 고정 로드 주소 8 곳(tests/link_l4.py 의
BASES, editbin 으로 옮긴 DLL 사본)에서 한 프로세스씩 돌려 팔마다 최솟값으로 견준다. 잰 것은 모의 1 초 동안 step 호출의 벽시계 합이다.
전제: 사본이 제 주소에 올라가고, 시간 회의 결정적 값이 results/measure_jacobian_drop.json 의 같은 회와 끝자리까지 같다.
결과 results/measure_time_stiff.json. 전제가 서면 0, 깨지면 1, 준비 안 됨 2.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import statistics  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gate_jacobian_damping as G  # noqa: E402
import link_l4 as L4  # noqa: E402

OUT = os.path.join(G.ROOT, 'results', 'measure_time_stiff.json')
REC = os.path.join(G.ROOT, 'results', 'measure_jacobian_drop.json')
N = 81
CASES = [(1e3, 200), (1e4, 200), (1e5, 200), (1e4, 50), (1e5, 50)]     # (BE_U, 단계 / 초)
ARMS = G.ARMS                                                          # off = 블록 대칭, gn_state = GN
COMMIT = '441b685'
ITER_COST = 1.05                                                       # 예측에 쓴 반복 한 번의 비용 비(L4)
FP = ('usable', 'v_tip_peak', 'pcg_iters_mean', 'pcg_iters_max')       # 기록과 견주는 결정적 값


def child(a):
    n, sps = a['n'], a['sps']
    d, V, s, row, skipped = G.make(a['dll'], n, a['bend'], 1.0 / sps, True, G.GRAV)
    tip = row == n - 1
    res = {'base': hex(int(d._lib._handle)), 'build_info': d.build_info, 'threads': int(d.openmp_threads), 'skipped': skipped}
    stop, steps, vtip, iters, step_time = '', 0, [], [], 0.0
    try:
        for k in range(1, int(round(G.B2_T * sps)) + 1):
            t0 = time.perf_counter()
            try:
                d.step(1)
            except Exception as e:   # StepDiverged 등
                stop = 'diverged: %s' % str(e)[:100]
                break
            step_time += time.perf_counter() - t0
            X = np.asarray(d.positions(), float)
            v = np.asarray(d.velocities(), float)
            if not (np.isfinite(X).all() and np.isfinite(v).all()):
                stop = 'nonfinite'
                break
            it, r, f = G.pcg_flags(d)
            iters.append(it)
            if f:
                stop = 'pcg_numerical_failure'
                break
            if r:
                stop = 'pcg_retry'
                break
            steps = k
            vtip.append(float(np.linalg.norm(v, axis=1)[tip].mean()))
    finally:
        d.close()
    res.update({'usable': stop == '', 'stop': stop, 'steps': steps, 'step_seconds': step_time,
                'v_tip_peak': max(vtip) if vtip else 0.0, 'pcg_iters_sum': int(sum(iters)),
                'pcg_iters_mean': float(np.mean(iters)) if iters else 0.0, 'pcg_iters_max': int(max(iters)) if iters else 0})
    with open(a['out'], 'w', encoding='utf-8') as fh:
        json.dump(res, fh)
    return 0


def spawn(td, dll, arm, bend, sps, tag):
    out = os.path.join(td, tag + '.json')
    arg = os.path.join(td, tag + '.arg.json')
    with open(arg, 'w', encoding='utf-8') as fh:
        json.dump({'dll': dll, 'n': N, 'bend': bend, 'sps': sps, 'out': out}, fh)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1')
    env.pop('CFMS_MGCG_ADAPTIVE', None)
    env['CFMS_CPU_SXJ_BEND'] = arm           # 끔도 'off' 로 명시한다(브랜치의 기본은 gn_state)
    env['CFMS_CPU_SXJ_MEMBRANE'] = 'off'     # 인장 · 전단은 엔진 손 커널(논문용 측정과 같다)
    r = subprocess.run([sys.executable, os.path.abspath(__file__), '--child', arg], env=env, cwd=td, capture_output=True,
                       text=True, encoding='utf-8', errors='replace')
    if r.returncode != 0 or not os.path.isfile(out):
        raise RuntimeError('회 실패 (%s, %g, %d): %s' % (arm, bend, sps, r.stderr[-1500:]))
    with open(out, encoding='utf-8') as fh:
        return json.load(fh)


def main():
    dll = os.path.abspath(sys.argv[1])
    if not (os.path.isfile(dll) and os.path.isfile(REC)):
        print('없는 것:', dll if not os.path.isfile(dll) else REC)
        return 2
    with open(REC, encoding='utf-8') as fh:
        rec = {(r['n'], r['BE_U'], r['sps']): r for r in json.load(fh)['runs']}
    name = {'off': 'oh', 'gn_state': 'gn'}
    print('논문용 측정 2 — 뻣뻣한 굽힘에서의 계산 시간 (n %d, 깨끗한 설정, 주소 %d 곳 × 회 %d × 팔 %d)' % (N, len(L4.BASES), len(CASES), len(ARMS)))
    print(' 예측(실행 전): 시간 비 ≈ 반복 수 평균의 비 × %.2f' % ITER_COST)
    pred = {}
    for b, s in CASES:
        r = rec[(N, b, s)]
        it = r['gn']['pcg_iters_mean'] / r['oh']['pcg_iters_mean']
        pred[(b, s)] = it * ITER_COST
        print('  BE_U %-6g 1/%-3d 반복 GN %.2f · 블록 대칭 %.2f (비 %.3f) → 시간 비 예측 %.3f' % (
            b, s, r['gn']['pcg_iters_mean'], r['oh']['pcg_iters_mean'], it, pred[(b, s)]))
    t0 = time.perf_counter()
    td = tempfile.mkdtemp(prefix='cfmsSXJ_tstiff_')
    runs = {(arm, b, s): [] for b, s in CASES for arm in ARMS}
    try:
        copies = L4.make_rebased(dll, td)
        for i, cp in enumerate(copies):
            for b, s in CASES:
                for arm in ARMS:
                    info = spawn(os.path.dirname(cp), cp, arm, b, s, 't_%s_%g_%d' % (arm, b, s))
                    info['base_expected'] = L4.BASES[i]
                    runs[(arm, b, s)].append(info)
            print('  주소 %d / %d (%s) 끝 (%.0f s)' % (i + 1, len(copies), L4.BASES[i], time.perf_counter() - t0))
    finally:
        shutil.rmtree(td, ignore_errors=True)

    # 전제
    allr = [t for ts in runs.values() for t in ts]
    bi = allr[0]['build_info']
    a_ok = all(str(t['build_info'].get('git_commit', '')).startswith(COMMIT) and t['build_info'].get('dirty') is False for t in allr)
    b_bad = [(t['base_expected'], t.get('base')) for t in allr if int(t['base'], 16) != int(t['base_expected'], 16)]
    c_bad = []
    for (arm, b, s), ts in runs.items():
        want = rec[(N, b, s)][name[arm]]
        for t in ts:
            if any(t[k] != want[k] for k in FP):
                c_bad.append('%s|%g|%d@%s' % (arm, b, s, t['base_expected']))
    d_ok = all(t['threads'] == 1 for t in allr)
    premise = {'a_build': bool(a_ok), 'b_bases': not b_bad, 'c_deterministic': not c_bad, 'd_threads': bool(d_ok)}
    print(' 전제')
    print('  (a) %s  DLL git_commit %s · dirty %s' % ('PASS' if a_ok else 'FAIL', bi.get('git_commit'), bi.get('dirty')))
    print('  (b) %s  사본 %s' % ('PASS' if not b_bad else 'FAIL', '모두 제 주소' if not b_bad else '어긋남 %s' % b_bad))
    print('  (c) %s  시간 회 %d 개의 결정적 값이 기록과 %s' % ('PASS' if not c_bad else 'FAIL', len(allr),
                                                    '모두 같음' if not c_bad else '다름 %d 회 %s' % (len(c_bad), c_bad[:4])))
    print('  (d) %s  스레드 1' % ('PASS' if d_ok else 'FAIL'))
    valid = all(premise.values())

    rows = []
    print(' 모의 1 초의 step 시간 (ms) — 주소 8 곳의 최소 · 중앙, GN / 블록 대칭')
    for b, s in CASES:
        row = {'n': N, 'BE_U': b, 'sps': s, 'predicted_ratio': pred[(b, s)]}
        for arm in ARMS:
            ts = runs[(arm, b, s)]
            ms = [1e3 * t['step_seconds'] for t in ts]
            row[name[arm]] = {'ms_min': min(ms), 'ms_median': statistics.median(ms), 'ms_per_base': ms,
                              'pcg_iters_sum': ts[0]['pcg_iters_sum'], 'pcg_iters_mean': ts[0]['pcg_iters_mean'],
                              'us_per_iter_min': 1e3 * min(ms) / ts[0]['pcg_iters_sum'], 'v_tip_peak': ts[0]['v_tip_peak'],
                              'usable': all(t['usable'] for t in ts)}
        per = [g / o for g, o in zip(row['gn']['ms_per_base'], row['oh']['ms_per_base'])]
        row.update({'ratio_min': row['gn']['ms_min'] / row['oh']['ms_min'], 'ratio_median': row['gn']['ms_median'] / row['oh']['ms_median'],
                    'ratio_per_base_range': [min(per), max(per)], 'iters_ratio': row['gn']['pcg_iters_mean'] / row['oh']['pcg_iters_mean'],
                    'iter_cost_ratio': row['gn']['us_per_iter_min'] / row['oh']['us_per_iter_min']})
        rows.append(row)
        print('  BE_U %-6g 1/%-3d  블록 대칭 %8.1f · %8.1f   GN %8.1f · %8.1f   비(최소) %.3f (중앙 %.3f, 주소별 %.3f – %.3f)   '
              '예측 %.3f (잰 값 / 예측 %.3f)   반복 비 %.3f · 반복 한 번의 시간 비 %.3f' % (
                  b, s, row['oh']['ms_min'], row['oh']['ms_median'], row['gn']['ms_min'], row['gn']['ms_median'], row['ratio_min'],
                  row['ratio_median'], per and min(per), per and max(per), pred[(b, s)], row['ratio_min'] / pred[(b, s)],
                  row['iters_ratio'], row['iter_cost_ratio']))
    res = {'dll': dll, 'engine': {'git_commit': bi.get('git_commit'), 'dirty': bi.get('dirty')}, 'n': N, 'bases': L4.BASES,
           'config': 'clean (measure_jacobian_drop 과 같다), 수렴 PCG %g · 상한 %d, 단일 스레드, 모의 1 초의 step 벽시계 합' % (G.PCG_TOL, G.PCG_CAP),
           'premise': premise, 'valid': valid, 'premise_bad': {'bases': b_bad, 'deterministic': c_bad},
           'iter_cost_assumed': ITER_COST, 'runs': rows, 'seconds': time.perf_counter() - t0}
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
        fh.write('\n')
    print('전제 %s → %s (%.0f s)' % ('모두 섬' if valid else '깨짐 — 값을 쓰지 않는다', os.path.relpath(OUT, G.ROOT), res['seconds']))
    return 0 if valid else 1


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--child':
        with open(sys.argv[2], encoding='utf-8') as fh:
            sys.exit(child(json.load(fh)))
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main())
