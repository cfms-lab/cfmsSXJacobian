"""논문용 측정 (docs/GATE_JACOBIAN_DAMPING.md 「논문용 측정」, 판정 아님): 깨끗한 설정에서 시간 간격 · 굽힘 강성 · 메시에 따른 두 Jacobian 의 낙하 속력.

사용:  python tests/measure_jacobian_drop.py <브랜치 cfmsDrape.dll>     (리포 루트에서, run_gate 밖 — 약 5 분)
회는 tests/gate_jacobian_damping.py 의 낙하 회(깨끗한 설정: 공기 저항 · 보조 감쇠 · 굽힘 감쇠 끔)를 부른다. 두 팔(off = Oh 꼴, gn_state = GN).
끝 줄 평균 속력의 1 초 안 최댓값과, Oh 팔에 대해 식의 종단 속도 a / (hλ)(평평한 상태의 행렬, k_d = 0)를 함께 적는다.
결과 results/measure_jacobian_drop.json.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gate_jacobian_damping as G  # noqa: E402

OUT = os.path.join(G.ROOT, 'results', 'measure_jacobian_drop.json')
SPS_ALL = (50, 100, 200, 400, 800, 1600, 3200)
SPS_SHORT = (50, 100, 200, 400, 800)
CASES = ([(81, 1e4, s) for s in SPS_ALL] + [(81, b, s) for b in (1e3, 1e5) for s in SPS_SHORT]
         + [(n, 1e4, s) for n in (41, 21) for s in SPS_SHORT])


def main():
    dll = os.path.abspath(sys.argv[1])
    if not os.path.isfile(dll):
        print('DLL 없음:', dll)
        return 2
    td = tempfile.mkdtemp(prefix='cfmsSXJ_jdrop_')
    t0 = time.perf_counter()
    try:
        keys = sorted({(n, b) for n, b, _ in CASES})
        mats = G.run_all(td, {k: ('off', 'matrix', dict(dll=dll, n=k[0], bend=k[1])) for k in keys})
        info = mats[keys[0]]['build_info']
        jobs = {(arm, n, b, s): (arm, 'drop', dict(dll=dll, n=n, bend=b, sps=s, clean=True)) for n, b, s in CASES for arm in G.ARMS}
        got = G.run_all(td, jobs)
    finally:
        shutil.rmtree(td, ignore_errors=True)
    rows = []
    print('논문용 측정 — 깨끗한 설정의 낙하 (엔진 %s, dirty %s). 끝 줄 평균 속력의 1 초 안 최댓값, cm/s' % (info.get('git_commit'), info.get('dirty')))
    for n, b, s in CASES:
        K = mats[(n, b)]['K']
        lam, _, mr = G.eig(K, n)
        vss = float(np.linalg.solve(0.5 * (K + K.T), np.full(n - 1, mr * G.GRAV))[-1]) * s
        ro, rg = got[('off', n, b, s)], got[('gn_state', n, b, s)]
        row = {'n': n, 'BE_U': b, 'sps': s, 'h_lambda1': float(lam[0]) / s, 'v_ss': vss}
        for name, r in (('oh', ro), ('gn', rg)):
            row[name] = {k: r[k] for k in ('usable', 'stop', 'v_tip_peak', 't_peak', 'v_node_max', 'chord_min', 'chord_at_peak',
                                           'pcg_iters_mean', 'pcg_iters_max', 'skipped')}
        rows.append(row)
        print('  n %-2d BE_U %-6g 1/%-4d hλ1 %8.4g  식 %8.4g · Oh %8.4g%s (비 %.3f) · GN %8.4g%s  (PCG 평균 Oh %.0f GN %.0f)' % (
            n, b, s, row['h_lambda1'], vss, ro['v_tip_peak'], '' if ro['usable'] else '[' + ro['stop'][:10] + ']',
            ro['v_tip_peak'] / vss, rg['v_tip_peak'], '' if rg['usable'] else '[' + rg['stop'][:10] + ']',
            ro['pcg_iters_mean'], rg['pcg_iters_mean']))
    res = {'dll': dll, 'engine': {'git_commit': info.get('git_commit'), 'dirty': info.get('dirty')}, 'gravity': G.GRAV, 'plate_cm': G.PLATE,
           'config': 'clean (set_damping 0, air drag 0, bend damping 0, BV = air_damping = internal_damping = 0)',
           'rigid_plate_tip_speed': float(np.sqrt(3.0 * G.GRAV * G.PLATE)), 'runs': rows, 'seconds': time.perf_counter() - t0}
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
        fh.write('\n')
    print('결과 %s (%.0f s)' % (os.path.relpath(OUT, G.ROOT), res['seconds']))
    return 0


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main())
