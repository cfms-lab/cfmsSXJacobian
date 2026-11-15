"""엔진 연결 L4R (굽힘) — 엔진 Oh 대신 작은 h 해를 기준으로 다시 본 안정 · 차이. 판정문 docs/GATE_LINK_L4R.md.

사용:  python tests/link_l4r.py <브랜치 cfmsDrape.dll>     (리포 루트에서)
       python tests/link_l4r.py --ref <dll> <장면> <BE_U> <단계/초> <접두>   — 기준 회 하나(이 스크립트가 부른다)

144 회 + 결정성 2 회는 tests/link_l4.py 와 같은 계산이다(전제 (c) 가 results/gate_link_l4_2.json 과 견준다). 기준 회는 같은 장면 ·
측정에 수렴 PCG 의 상한만 2000 으로 올리고(상대 1e-6 그대로) 벽시계 상한을 3600 s 로 둔다. 굽힘 모드는 L4 와 같이 프로세스마다
환경 변수 CFMS_CPU_SXJ_BEND 로 준다. (다) 시간은 다시 재지 않고 results/gate_link_l4_2.json 의 것을 쓴다.
통과 0, 실패 1, 준비 안 됨 2. 결과 results/gate_link_l4r.json.
"""
import sys
sys.dont_write_bytecode = True

import concurrent.futures  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import link_l0 as L0  # noqa: E402
import link_l2 as L2  # noqa: E402
import link_l4 as L4  # noqa: E402

ROOT = L0.ROOT
OUT_JSON = os.path.join(ROOT, 'results', 'gate_link_l4r.json')
L4_2_JSON = os.path.join(ROOT, 'results', 'gate_link_l4_2.json')

# 판정문 「기준 해」 · 「판정 기준」에 고정한 값
REF_SPS = 6400
CHECK_SPS = 3200
REF_PCG_CAP = 2000
REF_TIMEOUT = 3600.0
REF_JOBS = (('gn_state', REF_SPS), ('gn_state', CHECK_SPS), ('off', REF_SPS), ('off', CHECK_SPS), ('gn_rest', REF_SPS))
SIDE_CM = {'P': 30.0, 'R': 50.0}
ACC_FRAC = 0.25                          # 정확히 쓸 수 있음: D(회, 기준) ≤ 0.25 L
EREF_FRAC = 0.02                         # 전제 (d): e_ref ≤ 0.02 L
OFF_APPROACH = 0.75                      # 전제 (d): D(끔 1/6400, 기준) ≤ 0.75 × D(끔 1/3200, 기준)
# L4R-3 (다) = L4 2 차 판정값 그대로(results/gate_link_l4_2.json, 커밋 b3fbd4a). run_gate.bat 20 단계가 그 파일을 다시 쓰면 시간이
# 다시 재여 통과 · 미달이 퍼짐 안에서 바뀔 수 있으므로, 판정은 이 고정값으로 하고 지금 파일의 값은 보고만 한다.
L4_2_TIME = {'gn_rest': (False, 1.0508761525618886), 'gn_state': (True, 1.0477111100581422)}


def ref_spawn(dll, mode, sc, bend, sps, prefix):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1')
    # 끔도 'off' 로 명시한다 — 브랜치의 기본(환경 변수 없음)이 gn_state 가 된 뒤(2026-11-13)로는 없음 ≠ 끔이다.
    env['CFMS_CPU_SXJ_BEND'] = mode
    # 인장 · 전단은 엔진 손 커널로 고정한다 — 브랜치의 기본이 on 이 된 뒤(2026-11-13)로는 없음 ≠ 끔이고, 이 판정은 굽힘만 본다.
    env['CFMS_CPU_SXJ_MEMBRANE'] = 'off'
    args = [sys.executable, os.path.abspath(__file__), '--ref', dll, sc, repr(float(bend)), str(sps), prefix]
    try:
        r = subprocess.run(args, env=env, cwd=os.path.dirname(prefix), capture_output=True, text=True,
                           encoding='utf-8', errors='replace', timeout=REF_TIMEOUT + 120)
    except subprocess.TimeoutExpired:
        return {'mode': mode, 'usable': False, 'stop': 'timeout (killed)', 'steps': None}
    if r.returncode != 0 or not os.path.isfile(prefix + '.json'):
        raise RuntimeError('기준 회 실패 (%s, %s, %g, %d): %s' % (mode, sc, bend, sps, r.stderr[-2000:]))
    with open(prefix + '.json', encoding='utf-8') as fh:
        info = json.load(fh)
    info['stderr'] = r.stderr.strip()[-500:]
    return info


def fmt(x, f='%.2f'):
    return '-' if x is None else f % x


def main():
    if len(sys.argv) == 7 and sys.argv[1] == '--ref':
        a = sys.argv
        return L4.child_run(a[2], a[3], float(a[4]), int(a[5]), 'conv', a[6], True,
                            pcg_cap=REF_PCG_CAP, timeout=REF_TIMEOUT)
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    dll = os.path.abspath(sys.argv[1])
    missing = [p for p in (dll, L2.HEADER_GEN, L4_2_JSON) if not os.path.isfile(p)]
    missing += [p for p in (os.path.join(L0.ENGINE_DIR, 'drape'), os.path.join(L2.LINK_DIR, 'drape')) if not os.path.isdir(p)]
    if missing:
        print('없는 것: %s' % ', '.join(missing))
        return 2
    with open(L4_2_JSON, encoding='utf-8') as fh:
        l42 = json.load(fh)

    configs = [(sc, b, s, p) for sc in L4.SCENES for b in L4.BENDS for s in L4.SPS for p in L4.PCGS]
    refsets = [(sc, b) for sc in L4.SCENES for b in L4.BENDS]
    print('엔진 연결 L4R (굽힘) — 회 %d (+ 결정성 %d) + 기준 회 %d, 동시 프로세스 %d\n  브랜치 DLL %s (워크트리 %s)'
          % (len(configs) * len(L4.MODES), len(L4.DET_MODES), len(refsets) * len(REF_JOBS), L4.PARALLEL, dll, L2.LINK_DIR))
    t_all = time.perf_counter()
    td = tempfile.mkdtemp(prefix='cfmsSXJ_L4R_')
    try:
        jobs = [('run', m, c, os.path.join(td, 'run_%03d' % i)) for i, (m, c) in
                enumerate((m, c) for c in configs for m in L4.MODES)]
        jobs += [('det', m, L4.DET, os.path.join(td, 'det_%s' % m)) for m in L4.DET_MODES]
        jobs += [('ref', m, (sc, b, s), os.path.join(td, 'ref_%s_%s_%g_%d' % (m, sc, b, s)))
                 for sc, b in refsets for m, s in REF_JOBS]
        # 오래 걸리는 것(P 의 기준 회)을 먼저 넣는다
        jobs.sort(key=lambda j: (j[2][0] != 'P', j[0] != 'ref', -j[2][2]))
        runs, refs = {}, {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=L4.PARALLEL) as ex:
            futs = {}
            for kind, m, c, pre in jobs:
                if kind == 'ref':
                    fu = ex.submit(ref_spawn, dll, m, c[0], c[1], c[2], pre)
                else:
                    fu = ex.submit(L4.spawn, dll, m, c, pre, True)
                futs[fu] = (kind, m, c, pre)
            for n_done, fu in enumerate(concurrent.futures.as_completed(futs), 1):
                kind, m, c, pre = futs[fu]
                info = fu.result()
                info['prefix'] = pre
                if kind == 'ref':
                    refs[(m, c)] = info
                else:
                    runs[('det:' if kind == 'det' else '') + L4.key(m, c)] = info
                if n_done % 24 == 0 or n_done == len(jobs):
                    print('  회 %d / %d (%.0f s)' % (n_done, len(jobs), time.perf_counter() - t_all))
        t_runs = time.perf_counter() - t_all

        def dref(pre_a, pre_b):
            return L4.dist(pre_a + '.npy', pre_b + '.npy')

        # 전제
        st = L2.branch_state()
        sha_gen = L2.sha256_file(L2.HEADER_GEN)
        sha_eng = L2.sha256_file(L2.HEADER_ENGINE) if os.path.isfile(L2.HEADER_ENGINE) else ''
        bis = [r['build_info'] for r in list(runs.values()) + list(refs.values()) if 'build_info' in r]
        bi = bis[0] if bis else {}
        same_bi = all(b.get('git_commit') == bi.get('git_commit') and b.get('dirty') == bi.get('dirty') for b in bis)
        l42_commit = l42.get('L4_0', {}).get('build_info', {}).get('git_commit')
        fail0 = 0
        print(' L4R-0 전제')
        a_ok = (same_bi and st['head'] != '' and str(bi.get('git_commit', '')).startswith(st['head'][:7])
                and bi.get('dirty') is False and st['descends_from_base'] and not st['outside_allowed'] and not st['uncommitted']
                and bi.get('git_commit') == l42_commit)
        fail0 = L4.gate(fail0, a_ok, '(a)', 'DLL git_commit %s · 워크트리 HEAD %s · dirty %s · 미커밋 %s · 허용 밖 %s · L4 2 차 %s'
                        % (bi.get('git_commit'), st['head'], bi.get('dirty'), st['uncommitted'], st['outside_allowed'] or '없음',
                           l42_commit), 'HEAD 일치 · dirty false · 허용 목록 안 · L4 2 차와 같은 커밋')
        fail0 = L4.gate(fail0, sha_gen == sha_eng, '(b)', '헤더 사본 %s' % ('같음' if sha_gen == sha_eng else '다름'),
                        'SHA-256 같음')
        fr = l42['runs']
        diff = sorted(k for k, r in runs.items() if not (k in fr and all(fr[k].get(f) == r.get(f) for f in L4.FP_KEYS)))
        fail0 = L4.gate(fail0, not diff, '(c)', '%d 회 중 L4 2 차와 결정적 값이 다른 회 %d %s'
                        % (len(runs), len(diff), diff[:4] if diff else ''), '모두 같음')
        refinfo, d_bad = {}, []
        for sc, b in refsets:
            L = SIDE_CM[sc]
            rs = {(m, s): refs[(m, (sc, b, s))] for m, s in REF_JOBS}
            ok_all = all(r['usable'] for r in rs.values())
            rec = {'usable': {'%s|%d' % k: r['usable'] for k, r in rs.items()},
                   'stop': {'%s|%d' % k: r['stop'] for k, r in rs.items() if not r['usable']}}
            if ok_all:
                ref_pre = rs[('gn_state', REF_SPS)]['prefix']
                rec['e_ref'] = dref(rs[('gn_state', CHECK_SPS)]['prefix'], ref_pre)
                rec['D_off_3200'] = dref(rs[('off', CHECK_SPS)]['prefix'], ref_pre)
                rec['D_off_6400'] = dref(rs[('off', REF_SPS)]['prefix'], ref_pre)
                rec['D_rest_6400'] = dref(rs[('gn_rest', REF_SPS)]['prefix'], ref_pre)
                rec['ok'] = bool(rec['e_ref'] <= EREF_FRAC * L and rec['D_off_6400'] <= OFF_APPROACH * rec['D_off_3200'])
                rec['ref_prefix'] = ref_pre
            else:
                rec['ok'] = False
            refinfo['%s|%g' % (sc, b)] = rec
            if not rec['ok']:
                d_bad.append('%s · %g' % (sc, b))
        fail0 = L4.gate(fail0, not d_bad, '(d)', '기준 해 여섯 %s' % ('모두 섬' if not d_bad else '안 섬 %s' % d_bad),
                        '기준 회 쓸 수 있음 · e_ref ≤ %.2f L · 끔 1/%d ≤ %.2f × 끔 1/%d' % (EREF_FRAC, REF_SPS, OFF_APPROACH, CHECK_SPS))
        for sc, b in refsets:
            r = refinfo['%s|%g' % (sc, b)]
            print('      %s · BE_U %-5g e_ref %s cm · 끔 1/%d %s → 1/%d %s cm · gn_rest 1/%d %s cm%s'
                  % (sc, b, fmt(r.get('e_ref'), '%.3f'), CHECK_SPS, fmt(r.get('D_off_3200')), REF_SPS, fmt(r.get('D_off_6400')),
                     REF_SPS, fmt(r.get('D_rest_6400')), '' if 'e_ref' in r else '  멈춤 %s' % r['stop']))

        # 회마다 기준 거리 · 정확히 쓸 수 있음
        dist_ref, acc = {}, {}
        if not d_bad:
            for k, r in runs.items():
                if k.startswith('det:'):
                    continue
                m, sc, b, s, p = k.split('|')
                ref_pre = refinfo['%s|%s' % (sc, b)]['ref_prefix']
                dist_ref[k] = dref(r['prefix'], ref_pre) if r['usable'] else None
                acc[k] = bool(r['usable'] and dist_ref[k] <= ACC_FRAC * SIDE_CM[sc])
    finally:
        shutil.rmtree(td, ignore_errors=True)

    res = {'gate': 'link_L4R_bend', 'status': 'NOT_JUDGED', 'recommendation': None,
           'L4R_0': {'branch': st, 'build_info': bi, 'l4_2_commit': l42_commit,
                     'header_sha256': {'generated': sha_gen, 'engine_copy': sha_eng},
                     'same_as_l4_2_different': diff, 'reference': {k: {kk: vv for kk, vv in v.items() if kk != 'ref_prefix'}
                                                                    for k, v in refinfo.items()}, 'fail': fail0},
           'config': {'ref_sps': REF_SPS, 'check_sps': CHECK_SPS, 'ref_pcg_cap': REF_PCG_CAP, 'ref_timeout_s': REF_TIMEOUT,
                      'side_cm': SIDE_CM, 'acc_frac': ACC_FRAC, 'eref_frac': EREF_FRAC, 'off_approach': OFF_APPROACH,
                      'parallel': L4.PARALLEL},
           'refs': {'%s|%s|%g|%d' % (m, c[0], c[1], c[2]): {kk: vv for kk, vv in r.items() if kk not in ('prefix', 'build_info')}
                    for (m, c), r in sorted(refs.items())},
           'seconds': {'runs': t_runs, 'all': time.perf_counter() - t_all}}
    if fail0:
        with open(OUT_JSON, 'w', encoding='utf-8', newline='\n') as fh:
            json.dump(res, fh, ensure_ascii=False, indent=1)
        print('\n엔진 연결 L4R (굽힘) 전제 실패 %d — 판정하지 않음 → %s' % (fail0, os.path.relpath(OUT_JSON, ROOT)))
        return 2

    # 판정
    print(' 모드별 판정 (가R) 안정 · (나R) 차이 · (다) 시간(L4 2 차)')
    pairs, verdict, fail = [], {}, 0
    for c in configs:
        sc, b, s, p = c
        ko = L4.key('off', c)
        eref = refinfo['%s|%g' % (sc, b)]['e_ref']
        for m in L4.GN_MODES:
            km = L4.key(m, c)
            rec = {'mode': m, 'scene': sc, 'BE_U': b, 'sps': s, 'pcg': p,
                   'off_usable': runs[ko]['usable'], 'gn_usable': runs[km]['usable'],
                   'off_D_ref': dist_ref[ko], 'gn_D_ref': dist_ref[km], 'off_acc': acc[ko], 'gn_acc': acc[km],
                   'violation': bool(acc[ko] and not acc[km]), 'e_ref': eref}
            if acc[ko] and acc[km]:
                rec['b_counted'] = True
                rec['b_ok'] = bool(dist_ref[km] <= dist_ref[ko] + 2.0 * eref)
            pairs.append(rec)
    for m in L4.GN_MODES:
        pm = [r for r in pairs if r['mode'] == m]
        viol = [r for r in pm if r['violation']]
        counted = [r for r in pm if r.get('b_counted')]
        b_bad = [r for r in counted if not r['b_ok']]
        a_ok = not viol
        b_ok = bool(counted) and not b_bad
        c_ok, ratio = L4_2_TIME[m]
        v2 = l42['verdict'][m]
        if v2.get('time_ratio') != ratio:
            print('  (보고) 지금 %s 의 %s (다) 는 %s — 판정은 L4 2 차 고정값 %.3f'
                  % (os.path.relpath(L4_2_JSON, ROOT), m, v2.get('time_ratio'), ratio))
        worst = max((r['gn_D_ref'] - r['off_D_ref'] - 2.0 * r['e_ref'] for r in counted), default=None)
        print('  %s' % m)
        fail = L4.gate(fail, a_ok, '(가R)', '위반 %d / 짝 %d%s' % (
            len(viol), len(pm), '' if not viol else ' — %s' % ', '.join('%s·%g·1/%d·%s' % (r['scene'], r['BE_U'], r['sps'], r['pcg'])
                                                                       for r in viol[:6])), '위반 0')
        fail = L4.gate(fail, b_ok, '(나R)', '센 짝 %d, 넘음 %d · 못 넘음 %d, 최대 D(M)−D(끔)−2e_ref %s cm'
                       % (len(counted), len(counted) - len(b_bad), len(b_bad), fmt(worst, '%.3f')), '센 짝 모두 ≤, 하나 이상')
        fail = L4.gate(fail, c_ok, '(다)', 'L4 2 차 그대로 %s' % ('-' if ratio is None else '%.3f' % ratio), '≤ %.2f' % L4.TIME_RATIO)
        verdict[m] = {'a_ok': a_ok, 'violations': len(viol), 'pairs': len(pm), 'b_ok': b_ok, 'b_counted': len(counted),
                      'b_failed': len(b_bad), 'b_worst_excess_cm': worst, 'c_ok': c_ok, 'time_ratio': ratio,
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
    acc_counts, use_counts = {}, {}
    for m in L4.MODES:
        for sc in L4.SCENES:
            for p in L4.PCGS:
                ks = [L4.key(m, c) for c in configs if c[0] == sc and c[3] == p]
                acc_counts['%s|%s|%s' % (m, sc, p)] = sum(1 for k in ks if acc[k])
                use_counts['%s|%s|%s' % (m, sc, p)] = sum(1 for k in ks if runs[k]['usable'])
    print(' 보고 — 정확히 쓸 수 있는 회 / 쓸 수 있는 회 (장면 · PCG 마다 BE_U 3 × h 4 = 12 회 중)')
    for sc in L4.SCENES:
        for p in L4.PCGS:
            print('  %s · %-7s ' % (sc, p) + ' · '.join('%s %d / %d' % (m, acc_counts['%s|%s|%s' % (m, sc, p)],
                                                                       use_counts['%s|%s|%s' % (m, sc, p)]) for m in L4.MODES))
    print(' 보고 — 기준 거리 D(회, 기준) cm (쓸 수 있는 회만, * = 정확히 쓸 수 있음)')
    for sc in L4.SCENES:
        for b in L4.BENDS:
            for p in L4.PCGS:
                cells = []
                for s in L4.SPS:
                    cells.append('1/%d ' % s + ' '.join(
                        '%s%s' % (fmt(dist_ref[L4.key(m, (sc, b, s, p))], '%.1f'), '*' if acc[L4.key(m, (sc, b, s, p))] else '')
                        for m in L4.MODES))
                print('  %s · %-5g · %-7s (끔 rest state)  %s' % (sc, b, p, ' | '.join(cells)))
    print('  실행 시간: 회 %.0f s' % t_runs)

    status = 'PASS' if rec != 'off' else 'FAIL'
    res.update({'status': status, 'recommendation': rec, 'verdict': verdict, 'pairs': pairs,
                'accurate_counts': acc_counts, 'usable_counts': use_counts,
                'runs': {k: {**{kk: vv for kk, vv in r.items() if kk not in ('prefix', 'build_info')},
                             'D_ref': dist_ref.get(k), 'accurate': acc.get(k)} for k, r in sorted(runs.items())}})
    with open(OUT_JSON, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
    print('\n엔진 연결 L4R (굽힘) %s — 권고: %s → %s' % (status, rec_text, os.path.relpath(OUT_JSON, ROOT)))
    return 0 if status == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
