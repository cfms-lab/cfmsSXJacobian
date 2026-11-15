"""한 모드 식 시험 2 차 (docs/GATE_JACOBIAN_DAMPING.md 「2 차」): 식에 엔진 기본의 법선 공기 저항 k_d 를 넣고, 1 차에서 돌리지 않은 회로 다시 본다.

사용:  python tests/gate_jacobian_damping_2.py <브랜치 cfmsDrape.dll>     (리포 루트에서, run_gate 밖 — 약 3 분)
법선 모드의 단계 = (v + h a) / (1 + h² λ + h k_d), 종단 속도 V_ss′ = [(K̃_rows + (k_d / h) M_rows)⁻¹ M_rows g]_끝 / h.
2-B0 한 모드의 감쇠(단계 수 규칙을 고침) · 2-D 공기 저항 항만(무중력, 원래 설정) · 2-B2 새 회의 종단 속도. 회는 1 차 하네스의 것을 부른다.
끝에 1 차 원래 설정의 종단 영역 아홉 회를 고친 식으로 다시 견준다(사후 — 판정에 쓰지 않는다).
결과 results/gate_jacobian_damping_2.json. 통과 0, 미달 1, 준비 안 됨 2.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gate_jacobian_damping as G  # noqa: E402

OUT = os.path.join(G.ROOT, 'results', 'gate_jacobian_damping_2.json')
FIRST = os.path.join(G.ROOT, 'results', 'gate_jacobian_damping.json')
KD = 5.0                        # 엔진 기본 법선 공기 저항(1/s) — CFMS_CPU_AIR_DRAG_DEFAULT
BASE = (81, 1e4)
D_SPS, D_RATE, D_STEPS_GN = 200, 0.05, 40
TERMINAL = 25.0                 # hλ1 + k_d ≥ 25 / s
# (n, BE_U, 단계/초, 공기 저항 — None 이면 엔진 기본 5)
B2_CASES = ([(81, 3e4, s, None) for s in (400, 200, 100, 50)] + [(41, 1e5, s, None) for s in (400, 200, 100, 50)]
            + [(21, 1e6, s, None) for s in (200, 100, 50)] + [(81, 1e4, s, 20.0) for s in (200, 100, 50)] + [(81, 1e4, 100, 0.0)])
MATRIX_CASES = sorted({(n, b) for n, b, _, _ in B2_CASES} | {BASE, (81, 1e5), (41, 1e4)})


def v_ss(K, n, sps, kd):
    mr = n * G.node_mass(n)
    Ks = 0.5 * (K + K.T)
    return float(np.linalg.solve(Ks + kd * sps * mr * np.eye(n - 1), np.full(n - 1, mr * G.GRAV))[-1]) * sps


def main():
    dll = os.path.abspath(sys.argv[1])
    if not os.path.isfile(dll):
        print('DLL 없음:', dll)
        return 2
    td = tempfile.mkdtemp(prefix='cfmsSXJ_jdamp2_')
    t0 = time.perf_counter()
    fail, out = [], {'dll': dll, 'k_d_default': KD, 'pcg': [G.PCG_TOL, G.PCG_CAP]}
    try:
        mats = G.run_all(td, {(n, b): ('off', 'matrix', dict(dll=dll, n=n, bend=b)) for n, b in MATRIX_CASES})
        info = mats[BASE]['build_info']
        out['engine'] = {'git_commit': info.get('git_commit'), 'dirty': info.get('dirty')}
        n0, b0 = BASE
        lam, U, mrow = G.eig(mats[BASE]['K'], n0)
        srow = np.linspace(0.0, G.PLATE, n0)[1:]

        # ---- 예측 (돌리기 전에 찍는다)
        pred, jobs = {'lambda': [float(x) for x in lam[:6]], 'B0': {}, 'D': {}, 'B2': {}}, {}
        for sps in G.B0_SPS:
            for k in G.B0_MODES:
                fac = 1.0 / (1.0 + lam[k - 1] / sps ** 2)
                N = int(max(1, min(200, math.floor(math.log(0.01) / math.log(fac)))))
                pred['B0']['%d|%d' % (sps, k)] = {'factor': fac, 'steps': N}
                u = U[:, k - 1]
                w0 = os.path.join(td, 'w0_%d_%d.npy' % (sps, k))
                np.save(w0, u / np.abs(u).max())
                jobs[('B0', 'off', sps, k)] = ('off', 'decay', dict(dll=dll, n=n0, bend=b0, sps=sps, steps=N, w0=w0))
                if k == 1:
                    jobs[('B0', 'gn_state', sps, k)] = ('gn_state', 'decay', dict(dll=dll, n=n0, bend=b0, sps=sps, steps=N, w0=w0))
        w0r = D_RATE * srow
        bk = U.T @ (mrow * w0r)
        dec = 1.0 / (1.0 + lam / D_SPS ** 2 + KD / D_SPS)
        ser = np.array([float((bk * U[-1, :] * dec ** k).sum()) for k in range(0, 2001)]) / float((bk * U[-1, :]).sum())
        n_oh = int(np.argmax(ser <= 0.2))
        pred['D'] = {'gn_steps': D_STEPS_GN, 'gn_ratio': (1.0 + KD / D_SPS) ** (-D_STEPS_GN), 'oh_steps': n_oh, 'oh_ratio': float(ser[n_oh])}
        w0f = os.path.join(td, 'w0_rigid.npy')
        np.save(w0f, w0r)
        jobs[('D', 'gn_state')] = ('gn_state', 'decay', dict(dll=dll, n=n0, bend=b0, sps=D_SPS, steps=D_STEPS_GN, w0=w0f, clean=False))
        jobs[('D', 'off')] = ('off', 'decay', dict(dll=dll, n=n0, bend=b0, sps=D_SPS, steps=n_oh, w0=w0f, clean=False))
        for n, b, sps, air in B2_CASES:
            kd = KD if air is None else air
            l1 = float(G.eig(mats[(n, b)]['K'], n)[0][0])
            pred['B2']['%d|%g|%d|%g' % (n, b, sps, kd)] = {'v_ss': v_ss(mats[(n, b)]['K'], n, sps, kd), 'h_lambda1': l1 / sps, 'k_d': kd,
                                                           'terminal': l1 / sps + kd >= TERMINAL}
            jobs[('B2', n, b, sps, kd)] = ('off', 'drop', dict(dll=dll, n=n, bend=b, sps=sps, clean=False, air_drag=air))
        out['predictions'] = pred
        print('2 차 예측 (엔진 %s, dirty %s; Oh 팔, k_d = 기본 %g 1/s)' % (info.get('git_commit'), info.get('dirty'), KD))
        for key, p in pred['B0'].items():
            print('  2-B0 %s: 단계당 인자 %.5f, N = %d' % (key, p['factor'], p['steps']))
        print('  2-D: GN %d 단계 비 %.4f · Oh %d 단계 비 %.4f' % (D_STEPS_GN, pred['D']['gn_ratio'], n_oh, pred['D']['oh_ratio']))
        for key, p in pred['B2'].items():
            print('  2-B2 %s: V_ss′ %.4g cm/s, hλ1 + k_d = %.4g 1/s%s' % (key, p['v_ss'], p['h_lambda1'] + p['k_d'], ' (종단)' if p['terminal'] else ''))
        with open(OUT + '.pred', 'w', encoding='utf-8') as fh:
            json.dump(pred, fh, indent=1)

        got = G.run_all(td, jobs)

        print('2-B0. 한 모드의 감쇠 (무중력, 깨끗한 설정)')
        B0, ok = {}, True
        for sps in G.B0_SPS:
            for k in G.B0_MODES:
                p = pred['B0']['%d|%d' % (sps, k)]
                u, N = U[:, k - 1], p['steps']
                r = got[('B0', 'off', sps, k)]
                c = (np.array(r['W']) * (mrow * u)[None, :]).sum(axis=1)
                f = float((c[N] / c[0]) ** (1.0 / N)) if c[N] / c[0] > 0 else float('nan')
                ratio = math.log(f) / math.log(p['factor']) if f == f else float('nan')
                good = 0.95 <= ratio <= 1.05 and not (r['pcg_retry'] or r['pcg_failure'])
                row = {'predicted': p['factor'], 'steps': N, 'oh': f, 'log_ratio': ratio, 'pcg_retry': r['pcg_retry'], 'skipped': r['skipped']}
                text = ''
                if k == 1:
                    rg = got[('B0', 'gn_state', sps, k)]
                    cg = (np.array(rg['W']) * (mrow * u)[None, :]).sum(axis=1)
                    fg = float((cg[N] / cg[0]) ** (1.0 / N))
                    gnr = abs(math.log(fg)) / abs(math.log(p['factor']))
                    good = good and gnr <= 0.02
                    row.update({'gn': fg, 'gn_exponent_over_oh_predicted': gnr})
                    text = ' · GN %.6f (지수 비 %.2e)' % (fg, gnr)
                ok = ok and good
                B0['%d|%d' % (sps, k)] = row
                print('  1/%-3d 모드 %d (N = %d): 예측 %.5f · Oh %.5f (ln 비 %.4f)%s%s' % (sps, k, N, p['factor'], f, ratio, text, '' if good else '  ← 벗어남'))
        out['B0'] = B0
        G.gate(fail, ok, '2-B0', 'Oh 의 ln 비 ∈ [0.95, 1.05] 여섯 모두, GN 대조(모드 1) 지수 ≤ 2 %')

        print('2-D. 공기 저항 항만 (무중력, 원래 설정, 강체 회전, 1/%d)' % D_SPS)
        Wg, Wo = np.array(got[('D', 'gn_state')]['W']), np.array(got[('D', 'off')]['W'])
        rg, ro = float(Wg[D_STEPS_GN, -1] / Wg[0, -1]), float(Wo[n_oh, -1] / Wo[0, -1])
        relg, relo = rg / pred['D']['gn_ratio'], ro / pred['D']['oh_ratio']
        print('  GN %d 단계: 예측 %.4f · 측정 %.4f (측정 / 예측 %.4f)' % (D_STEPS_GN, pred['D']['gn_ratio'], rg, relg))
        print('  Oh %d 단계: 예측 %.4f · 측정 %.4f (측정 / 예측 %.4f)' % (n_oh, pred['D']['oh_ratio'], ro, relo))
        out['D'] = {'gn_measured': rg, 'gn_rel': relg, 'oh_measured': ro, 'oh_rel': relo,
                    'gn_per_step_rate': float(-math.log(rg) / D_STEPS_GN * D_SPS)}
        G.gate(fail, 0.95 <= relg <= 1.05 and 0.95 <= relo <= 1.05, '2-D', '두 팔의 측정 비 / 예측 비 ∈ [0.95, 1.05]')

        print('2-B2. 중력 아래의 종단 속도 — 새 회 (Oh 팔, 끝 줄 평균 속력의 1 초 안 최댓값, cm/s)')
        B2, ok, nterm = {}, True, 0
        for n, b, sps, air in B2_CASES:
            kd = KD if air is None else air
            p = pred['B2']['%d|%g|%d|%g' % (n, b, sps, kd)]
            r = got[('B2', n, b, sps, kd)]
            ratio = r['v_tip_peak'] / p['v_ss'] if r['usable'] else float('nan')
            judged = p['terminal'] and r['usable']
            good = (0.9 <= ratio <= 1.1) if judged else True
            nterm += int(judged)
            ok = ok and good and (r['usable'] or not p['terminal'])
            B2['%d|%g|%d|%g' % (n, b, sps, kd)] = {'v_ss': p['v_ss'], 'h_lambda1': p['h_lambda1'], 'k_d': kd, 'terminal': p['terminal'],
                                                   'judged': judged, 'ratio': ratio,
                                                   'oh': {k: r[k] for k in ('usable', 'stop', 'v_tip_peak', 't_peak', 'v_node_max',
                                                                            'pcg_iters_mean', 'pcg_iters_max')}}
            print('  n %-2d BE_U %-6g 1/%-3d k_d %-4g hλ1+k_d %8.4g %s 예측 %8.4g · 측정 %8.4g%s (비 %.3f%s)' % (
                n, b, sps, kd, p['h_lambda1'] + kd, '종단' if p['terminal'] else '    ', p['v_ss'], r['v_tip_peak'],
                '' if r['usable'] else '[' + r['stop'][:14] + ']', ratio, '' if good else ' ← 벗어남'))
        out['B2'] = B2
        G.gate(fail, ok and nterm > 0, '2-B2', '종단 영역의 새 회 %d 개 모두 0.9 ≤ V_peak / V_ss′ ≤ 1.1' % nterm)

        print('사후(판정 밖): 1 차 원래 설정의 종단 영역 회를 고친 식으로 다시 견줌')
        post = {}
        with open(FIRST, encoding='utf-8') as fh:
            first = json.load(fh)
        for key, e in first['B2'].items():
            n, b, sps, cfg = key.split('|')
            n, b, sps = int(n), float(b), int(sps)
            if not (e['judged'] and (n, b) in mats):
                continue
            v2 = v_ss(mats[(n, b)]['K'], n, sps, KD)
            post[key] = {'v_ss_first': e['v_ss'], 'v_ss_with_drag': v2, 'measured': e['oh']['v_tip_peak'],
                         'ratio_first': e['ratio'], 'ratio_with_drag': e['oh']['v_tip_peak'] / v2}
            print('  %-22s 측정 %8.4g · 1 차 예측 %8.4g (비 %.3f) · 고친 예측 %8.4g (비 %.3f)' % (
                key, e['oh']['v_tip_peak'], e['v_ss'], e['ratio'], v2, e['oh']['v_tip_peak'] / v2))
        out['post_first_round'] = post
    finally:
        shutil.rmtree(td, ignore_errors=True)
    out['fail'] = fail
    out['verdict'] = 'PASS' if not fail else 'FAIL'
    out['seconds'] = time.perf_counter() - t0
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
        fh.write('\n')
    if os.path.isfile(OUT + '.pred'):
        os.remove(OUT + '.pred')
    print('2 차 종합: %s%s  (%.0f s, 결과 %s)' % (out['verdict'], '' if not fail else ' — ' + ' '.join(fail), out['seconds'],
                                             os.path.relpath(OUT, G.ROOT)))
    return 0 if not fail else 1


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main())
