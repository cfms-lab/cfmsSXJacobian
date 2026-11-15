"""논문용 측정 3 (docs/GATE_JACOBIAN_DAMPING.md 「논문용 측정 3」, 판정 아님): 정확한 굽힘 행렬과 그 준정부호 투영으로 같은 낙하를 돌린다.

사용:  python tests/measure_jacobian_drop_exact.py <실험 cfmsDrape.dll>     (리포 루트에서, run_gate 밖 — 약 30 분)
실험 DLL = 엔진 브랜치에 exact_state · psd_state 모드를 넣은 빌드(%TEMP%\\cfmsDrape_build_sxjexact). 먼저 tests/engine_bend_exact.py 가 통과해야 한다.
회는 tests/gate_jacobian_damping.py 의 낙하 회(깨끗한 설정)다. 팔 다섯: off(블록 대칭) · gn_state · gn_rest · exact_state · psd_state.
전제: off · gn_state 회가 results/measure_jacobian_drop.json 의 같은 회와 결정적 값이 끝자리까지 같다.
결과 results/measure_jacobian_drop_exact.json. 전제가 서면 0, 깨지면 1, 준비 안 됨 2.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gate_jacobian_damping as G  # noqa: E402

OUT = os.path.join(G.ROOT, 'results', 'measure_jacobian_drop_exact.json')
REC = os.path.join(G.ROOT, 'results', 'measure_jacobian_drop.json')
GATE = os.path.join(G.ROOT, 'results', 'gate_engine_bend_exact.json')
N, BEND = 81, 1e4
SPS_ALL = (50, 100, 200, 400, 800, 1600, 3200)
SPS_PSD = (50, 100, 200, 400, 800, 1600)             # psd_state 는 힌지마다 고윳값 분해라 1/3200 을 뺀다
ARMS = ('off', 'gn_state', 'gn_rest', 'exact_state', 'psd_state')
EXTRA = [(arm, b, 200) for arm in ('psd_state', 'gn_rest') for b in (1e3, 1e5)]
FP = ('usable', 'v_tip_peak', 'pcg_iters_mean', 'pcg_iters_max')
KEEP = ('usable', 'stop', 'steps', 'v_tip_peak', 't_peak', 'v_node_max', 'chord_min', 'chord_at_peak', 'pcg_iters_mean', 'pcg_iters_max',
        'skipped', 'v_tip_series')


def halvings(v):
    """h 를 반으로 줄일 때마다의 증가분(쓸 수 있는 이웃 회 사이만)."""
    return [(b - a) if (a is not None and b is not None) else None for a, b in zip(v[:-1], v[1:])]


def main():
    dll = os.path.abspath(sys.argv[1])
    if not all(os.path.isfile(p) for p in (dll, REC, GATE)):
        print('없는 것:', [p for p in (dll, REC, GATE) if not os.path.isfile(p)])
        return 2
    with open(GATE, encoding='utf-8') as fh:
        gate = json.load(fh)
    if gate.get('gates_failed') != 0 or os.path.normcase(gate['dll']['new']) != os.path.normcase(dll):
        print('커널 확인(tests/engine_bend_exact.py)이 이 DLL 로 통과하지 않았다')
        return 2
    with open(REC, encoding='utf-8') as fh:
        rec = {(r['n'], r['BE_U'], r['sps']): r for r in json.load(fh)['runs']}
    jobs = {}
    for arm in ARMS:
        for s in (SPS_PSD if arm == 'psd_state' else SPS_ALL):
            jobs[(arm, BEND, s)] = (arm, 'drop', dict(dll=dll, n=N, bend=BEND, sps=s, clean=True))
    for arm, b, s in EXTRA:
        jobs[(arm, b, s)] = (arm, 'drop', dict(dll=dll, n=N, bend=b, sps=s, clean=True))
    print('논문용 측정 3 — 정확한 행렬 · 준정부호 투영의 낙하 (n %d, 깨끗한 설정, 회 %d, 동시 프로세스 %d)' % (N, len(jobs), G.PARALLEL))
    td = tempfile.mkdtemp(prefix='cfmsSXJ_jdropx_')
    t0 = time.perf_counter()
    try:
        got = G.run_all(td, jobs)
    finally:
        shutil.rmtree(td, ignore_errors=True)
    info = got[('off', BEND, 200)]['build_info']

    # 전제: 기존 두 팔이 기록과 같다
    name = {'off': 'oh', 'gn_state': 'gn'}
    bad = []
    for arm in name:
        for s in SPS_ALL:
            want, r = rec[(N, BEND, s)][name[arm]], got[(arm, BEND, s)]
            if any(r[k] != want[k] for k in FP):
                bad.append('%s|1/%d' % (arm, s))
    valid = not bad
    print(' 전제 %s  off · gn_state %d 회의 결정적 값이 기록과 %s (엔진 %s, dirty %s)' % (
        'PASS' if valid else 'FAIL', 2 * len(SPS_ALL), '모두 같음' if valid else '다름 %s' % bad, info.get('git_commit'), info.get('dirty')))

    runs = {'%s|%g|%d' % k: {kk: r.get(kk) for kk in KEEP} for k, r in sorted(got.items())}
    peak = {arm: [(got[(arm, BEND, s)]['v_tip_peak'] if (arm, BEND, s) in got and got[(arm, BEND, s)]['usable'] else None) for s in SPS_ALL]
            for arm in ARMS}
    inc = {arm: halvings(peak[arm]) for arm in ARMS}
    print(' 끝 줄 평균 속력의 1 초 안 최댓값 (cm/s), BE_U %g — 괄호는 멈춘 까닭' % BEND)
    print('  %-12s' % '1/h' + ''.join('%12d' % s for s in SPS_ALL))
    for arm in ARMS:
        cells = []
        for s in SPS_ALL:
            r = got.get((arm, BEND, s))
            cells.append('%12s' % ('-' if r is None else ('%.1f' % r['v_tip_peak'] if r['usable'] else '(%s@%d)' % (r['stop'][:8], r['steps']))))
        print('  %-12s' % arm + ''.join(cells))
    print(' h 를 반으로 줄일 때의 증가분 (cm/s)')
    for arm in ARMS:
        print('  %-12s' % arm + ''.join('%12s' % ('-' if v is None else '%.1f' % v) for v in inc[arm]))
    print(' gn_state 의 증가분에 대한 비')
    ratio = {}
    for arm in ('gn_rest', 'exact_state', 'psd_state'):
        ratio[arm] = [(a / g if (a is not None and g) else None) for a, g in zip(inc[arm], inc['gn_state'])]
        print('  %-12s' % arm + ''.join('%12s' % ('-' if v is None else '%.3f' % v) for v in ratio[arm]))
    print(' 같은 h 에서 gn_state 에 대한 최대 속력의 비 · PCG 반복 평균')
    rel = {}
    for arm in ('gn_rest', 'exact_state', 'psd_state'):
        rel[arm] = [(a / g if (a is not None and g is not None) else None) for a, g in zip(peak[arm], peak['gn_state'])]
        its = [got[(arm, BEND, s)]['pcg_iters_mean'] if (arm, BEND, s) in got else None for s in SPS_ALL]
        print('  %-12s' % arm + ''.join('%12s' % ('-' if v is None else '%.4f' % v) for v in rel[arm]))
        print('  %-12s' % '  (PCG)' + ''.join('%12s' % ('-' if v is None else '%.0f' % v) for v in its))
    print(' 굽힘 강성 1e3 · 1e5, h 1/200')
    for arm, b, s in EXTRA:
        r = got[(arm, b, s)]
        g = rec[(N, b, s)]['gn']
        print('  %-12s BE_U %-6g %s (gn_state 기록 %.1f, PCG %.0f 대 %.0f)' % (
            arm, b, '%.1f' % r['v_tip_peak'] if r['usable'] else '(%s@%d)' % (r['stop'][:8], r['steps']), g['v_tip_peak'],
            r['pcg_iters_mean'], g['pcg_iters_mean']))

    res = {'dll': dll, 'engine': {'git_commit': info.get('git_commit'), 'dirty': info.get('dirty')}, 'n': N, 'BE_U': BEND,
           'config': 'clean (measure_jacobian_drop 과 같다)', 'sps': list(SPS_ALL), 'sps_psd': list(SPS_PSD), 'arms': list(ARMS),
           'premise': {'old_arms_match_record': valid, 'different': bad}, 'valid': valid,
           'peak': peak, 'increment_per_halving': inc, 'increment_over_gn_state': ratio, 'peak_over_gn_state': rel,
           'runs': runs, 'seconds': time.perf_counter() - t0}
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
        fh.write('\n')
    print('전제 %s → %s (%.0f s)' % ('섬' if valid else '깨짐 — 값을 쓰지 않는다', os.path.relpath(OUT, G.ROOT), res['seconds']))
    return 0 if valid else 1


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main())
