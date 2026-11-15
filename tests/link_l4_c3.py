"""엔진 연결 L4 3 차 (다) — 시간 묶음을 다섯 번 재서 비의 중앙값으로 판정한다. 판정문 docs/GATE_LINK_L4.md 「3 차 (다)」.

사용:  python tests/link_l4_c3.py <브랜치 cfmsDrape.dll>     (리포 루트에서, 다른 작업이 없을 때)

한 묶음 = tests/link_l4.py 의 시간 회 그대로(고정 주소 8 곳의 DLL 사본, 기준 짝, 단일 스레드, 혼자). 주소마다 세 모드의 순서만
한 칸씩 돌린다. 적분 144 회는 돌지 않는다. 통과 0, 실패 1, 준비 안 됨 2. 결과 results/gate_link_l4_c3.json.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import statistics  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import link_l0 as L0  # noqa: E402
import link_l2 as L2  # noqa: E402
import link_l4 as L4  # noqa: E402

ROOT = L0.ROOT
OUT_JSON = os.path.join(ROOT, 'results', 'gate_link_l4_c3.json')
BATCHES = 5
PCG_SUM = {'off': 39042, 'gn_rest': 38882, 'gn_state': 38919}     # 2 차 시간 회의 PCG 반복 합(같은 계산인지 본다)
KNOWN = {'gn_rest': [1.0508761525618886, 1.0534648341233714], 'gn_state': [1.0477111100581422, 1.060857402919238]}


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    dll = os.path.abspath(sys.argv[1])
    if not os.path.isfile(dll) or not os.path.isdir(os.path.join(L2.LINK_DIR, 'drape')):
        print('없는 것: DLL 또는 브랜치 워크트리')
        return 2
    print('엔진 연결 L4 3 차 (다) — 묶음 %d × 주소 %d × 모드 %d, 혼자\n  브랜치 DLL %s' % (BATCHES, len(L4.BASES), len(L4.MODES), dll))
    t0 = time.perf_counter()
    td = tempfile.mkdtemp(prefix='cfmsSXJ_L4c3_')
    batches = []
    try:
        copies = L4.make_rebased(dll, td)
        for b in range(BATCHES):
            rec = {m: [] for m in L4.MODES}
            for i, cp in enumerate(copies):
                order = [L4.MODES[(i + k) % len(L4.MODES)] for k in range(len(L4.MODES))]
                for m in order:
                    info = L4.spawn(cp, m, L4.REF_TIME, os.path.join(os.path.dirname(cp), 't%d_%s' % (b, m)), False)
                    info['base_expected'] = L4.BASES[i]
                    rec[m].append(info)
            batches.append(rec)
            mins = {m: min(t['ms_per_sim_second'] for t in rec[m]) for m in L4.MODES}
            print('  묶음 %d (%.0f s): 최소 ' % (b + 1, time.perf_counter() - t0)
                  + ' · '.join('%s %.1f' % (m, mins[m]) for m in L4.MODES)
                  + ' — 비 ' + ' · '.join('%s %.4f' % (m, mins[m] / mins['off']) for m in L4.GN_MODES), flush=True)
    finally:
        shutil.rmtree(td, ignore_errors=True)

    st = L2.branch_state()
    allt = [t for rec in batches for m in L4.MODES for t in rec[m]]
    bi = allt[0].get('build_info', {})
    fail0 = 0
    print(' 전제')
    a_ok = (all(t.get('build_info', {}).get('git_commit') == bi.get('git_commit') and t.get('build_info', {}).get('dirty') is False
                for t in allt) and st['head'] != '' and str(bi.get('git_commit', '')).startswith(st['head'][:7])
            and st['descends_from_base'] and not st['outside_allowed'] and not st['uncommitted'])
    fail0 = L4.gate(fail0, a_ok, '(a)', 'DLL git_commit %s · 워크트리 HEAD %s · 미커밋 %s' % (bi.get('git_commit'), st['head'], st['uncommitted']),
                    'HEAD 일치 · dirty false · 허용 목록 안')
    c_bad = [t['base_expected'] for t in allt if t.get('base') is None or int(t['base'], 16) != int(t['base_expected'], 16)]
    fail0 = L4.gate(fail0, not c_bad, '(c)', '사본 %s' % ('모두 제 주소' if not c_bad else '어긋남 %s' % c_bad[:4]), '정한 주소에 올라감')
    u_bad = sum(1 for t in allt if not t['usable'])
    fail0 = L4.gate(fail0, u_bad == 0, '쓸 수 있음', '시간 회 %d 중 쓸 수 없는 회 %d' % (len(allt), u_bad), '0')
    p_bad = sum(1 for rec in batches for m in L4.MODES for t in rec[m] if t.get('pcg_iters_sum') != PCG_SUM[m])
    fail0 = L4.gate(fail0, p_bad == 0, 'PCG 합', '2 차와 다른 회 %d' % p_bad, '모두 같음 %s' % PCG_SUM)

    mins = [{m: min(t['ms_per_sim_second'] for t in rec[m]) for m in L4.MODES} for rec in batches]
    ratios = {m: [mn[m] / mn['off'] for mn in mins] for m in L4.GN_MODES}
    pooled = {m: min(mn[m] for mn in mins) for m in L4.MODES}
    verdict, fail = {}, 0
    print(' 판정 — 다섯 묶음의 비의 중앙값')
    for m in L4.GN_MODES:
        med = statistics.median(ratios[m])
        ok = med <= L4.TIME_RATIO
        fail = L4.gate(fail, ok, m, '중앙 %.4f (묶음 %s, 최소 %.4f · 최대 %.4f) · 합친 최소의 비 %.4f · 이미 잰 값 %s'
                       % (med, ' '.join('%.4f' % r for r in ratios[m]), min(ratios[m]), max(ratios[m]),
                          pooled[m] / pooled['off'], ' '.join('%.3f' % r for r in KNOWN[m])), '≤ %.2f' % L4.TIME_RATIO)
        verdict[m] = {'c_ok': bool(ok), 'median_ratio': med, 'ratios': ratios[m], 'pooled_min_ratio': pooled[m] / pooled['off']}
    res = {'gate': 'link_L4_bend_time_3', 'status': 'NOT_JUDGED' if fail0 else ('PASS' if fail == 0 else 'FAIL'),
           'prereq_fail': fail0, 'build_info': bi, 'branch': st, 'batches': BATCHES, 'bases': L4.BASES,
           'config': list(L4.REF_TIME), 'time_ratio': L4.TIME_RATIO, 'verdict': verdict, 'batch_min_ms_per_sim_second': mins,
           'pooled_min_ms_per_sim_second': pooled, 'known_before': KNOWN,
           'per_batch': [{m: [{k: t.get(k) for k in ('base', 'usable', 'ms_per_sim_second', 'pcg_iters_sum')} for t in rec[m]]
                          for m in L4.MODES} for rec in batches],
           'seconds': time.perf_counter() - t0}
    with open(OUT_JSON, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
    if fail0:
        print('\n3 차 (다) 전제 실패 %d — 판정하지 않음 → %s' % (fail0, os.path.relpath(OUT_JSON, ROOT)))
        return 2
    print('\n3 차 (다) — ' + ' · '.join('%s %s' % (m, 'PASS' if verdict[m]['c_ok'] else 'FAIL') for m in L4.GN_MODES)
          + ' → %s (%.0f s)' % (os.path.relpath(OUT_JSON, ROOT), res['seconds']))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
