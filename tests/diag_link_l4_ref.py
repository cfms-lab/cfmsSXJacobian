"""판정 밖 진단 (docs/GATE_LINK_L4.md 「판정 뒤 진단」): L4 가 왜 걸렸는가.

A. 작은 h 기준. L4 (나) 는 D(GN, Oh) 를 Oh 의 시간 간격 오차와 견줬다. 그러나 Oh 의 J 도 참 Jacobian 이 아니라서, D(GN, Oh) 가
   크다는 것만으로는 어느 쪽이 참 궤적에서 먼지 모른다. 힘이 같으면(끔 · gn_state) 선형화가 달라도 h → 0 에서 같은 궤적으로
   간다 — 그 궤적을 h = 1/1600 로 근사해 기준으로 둔다(장면 R · P, BE_U 1e2 · 1e4, 수렴 PCG).
B. 뻣뻣한 굽힘(R · BE_U 1e4)에서는 1/1600 도 수렴하지 않았다 — 1/3200 · 1/6400 으로 수렴을 보고, gn_state 1/6400 을 기준으로
   네 h · 두 PCG 의 세 모드 거리를 잰다.
C. 시간. 엔진 기본 설정(P · BE_U 1 · 1/200 · 수렴, 단일 스레드)을 엔진 위상 프로파일러(CFMS_PROF)로 한 번씩 돌려 굽힘 · PCG
   고리 시간을 나눈다(고정 로드 주소가 아닌 한 번 잰 값 — 나눠 보기만 한다).

사용:  python tests/diag_link_l4_ref.py <브랜치 cfmsDrape.dll>     (리포 루트에서). 결과 results/diag_link_l4_ref.json.
회는 tests/link_l4.py 와 같은 장면 · 측정 · 단일 스레드 프로세스다.
"""
import sys
sys.dont_write_bytecode = True

import concurrent.futures  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import link_l4 as L4  # noqa: E402

REF_SPS = 1600
SCENES = ('R', 'P')
BENDS = (1e2, 1e4)
B_SPS = (1600, 3200, 6400)
B_MODES = ('off', 'gn_state')
OUT = os.path.join(L4.ROOT, 'results', 'diag_link_l4_ref.json')


def fmt(x, f='%.3g'):
    return '-' if x is None else f % x


def profile(dll, mode, td):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1', CFMS_PROF='400')
    # 끔도 'off' 로 명시한다 — 브랜치의 기본(환경 변수 없음)이 gn_state 가 된 뒤(2026-11-13)로는 없음 ≠ 끔이다.
    env['CFMS_CPU_SXJ_BEND'] = mode
    # 인장 · 전단은 엔진 손 커널로 고정한다 — 브랜치의 기본이 on 이 된 뒤(2026-11-13)로는 없음 ≠ 끔이고, 이 판정은 굽힘만 본다.
    env['CFMS_CPU_SXJ_MEMBRANE'] = 'off'
    sc, b, s, p = L4.REF_TIME
    r = subprocess.run([sys.executable, os.path.join(L4.ROOT, 'tests', 'link_l4.py'), '--run', dll, sc, repr(b), str(s), p,
                        os.path.join(td, 'prof_' + mode), '0'], env=env, cwd=td, capture_output=True, text=True,
                       encoding='utf-8', errors='replace')
    out = {}
    for line in r.stderr.splitlines():
        m = re.match(r'\s*(\S.*?)\s+([0-9.]+) ms\s+[0-9.]+%$', line)
        if m:
            out[m.group(1).strip()] = float(m.group(2))
        m = re.search(r'last (\d+) steps, ([0-9.]+) ms/step', line)
        if m:
            out['ms_per_step'] = float(m.group(2))
        m = re.search(r'PCG: ([0-9.]+) iters/solve', line)
        if m:
            out['pcg_iters_per_solve'] = float(m.group(1))
    return out


def main():
    dll = os.path.abspath(sys.argv[1])
    td = tempfile.mkdtemp(prefix='cfmsSXJ_L4ref_')
    t0 = time.perf_counter()
    res = {'ref_sps': REF_SPS}
    try:
        jobs = [(m, (sc, b, s, 'conv')) for sc in SCENES for b in BENDS for s in (REF_SPS,) + L4.SPS for m in L4.MODES]
        jobs += [(m, ('R', 1e4, s, 'conv')) for s in B_SPS[1:] for m in B_MODES]
        jobs += [(m, ('R', 1e4, s, 'fixed20')) for s in L4.SPS for m in L4.MODES]
        jobs.sort(key=lambda j: (-j[1][2] * (10 if j[1][0] == 'P' else 1)))
        runs = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=L4.PARALLEL) as ex:
            futs = {}
            for i, (m, c) in enumerate(jobs):
                pre = os.path.join(td, 'r%03d' % i)
                futs[ex.submit(L4.spawn, dll, m, c, pre, True)] = (m, c, pre)
            for fu in concurrent.futures.as_completed(futs):
                m, c, pre = futs[fu]
                info = fu.result()
                info['prefix'] = pre
                runs[(m, c)] = info

        def D(a, b):
            ra, rb = runs[a], runs[b]
            return L4.dist(ra['prefix'] + '.npy', rb['prefix'] + '.npy') if ra['usable'] and rb['usable'] else None

        # A
        rows = []
        for sc in SCENES:
            for b in BENDS:
                ref = {m: ('%s' % m, (sc, b, REF_SPS, 'conv')) for m in L4.MODES}
                head = {'scene': sc, 'BE_U': b, 'ref_usable': {m: runs[ref[m]]['usable'] for m in L4.MODES},
                        'D_ref_state_vs_off': D(ref['gn_state'], ref['off']), 'D_ref_rest_vs_off': D(ref['gn_rest'], ref['off'])}
                rows.append(head)
                for s in L4.SPS:
                    for m in L4.MODES:
                        k = (m, (sc, b, s, 'conv'))
                        rows.append({'scene': sc, 'BE_U': b, 'sps': s, 'mode': m, 'usable': runs[k]['usable'],
                                     'stop': runs[k]['stop'], 'pcg_iters_mean': runs[k]['pcg_iters_mean'],
                                     'D_to_ref_off': D(k, ref['off']), 'D_to_ref_rest': D(k, ref['gn_rest'])})
        res['A'] = rows
        # B
        conv = {'D_off_vs_gn_state': {s: D(('off', ('R', 1e4, s, 'conv')), ('gn_state', ('R', 1e4, s, 'conv'))) for s in B_SPS}}
        for m in B_MODES:
            conv['D_%s_halving' % m] = {'%d-%d' % (a, 2 * a): D((m, ('R', 1e4, a, 'conv')), (m, ('R', 1e4, 2 * a, 'conv')))
                                        for a in B_SPS[:-1]}
        ref6400 = ('gn_state', ('R', 1e4, B_SPS[-1], 'conv'))
        table = []
        for p in L4.PCGS:
            for s in L4.SPS:
                table.append({'pcg': p, 'sps': s, **{m: D((m, ('R', 1e4, s, p)), ref6400) for m in L4.MODES}})
        res['B'] = {'convergence': conv, 'D_to_gn_state_6400': table}
        # C (혼자 돌린다)
        res['C'] = {'config': list(L4.REF_TIME), 'profile_ms_400_steps': {m: profile(dll, m, td) for m in ('off', 'gn_state')}}
    finally:
        shutil.rmtree(td, ignore_errors=True)
    res['seconds'] = time.perf_counter() - t0

    print('판정 밖 진단 A — 기준 = 같은 장면 · BE_U 의 h = 1/%d (수렴 PCG), D 는 cm' % REF_SPS)
    for r in res['A']:
        if 'sps' not in r:
            print(' %s · BE_U %g  D(gn_state 기준, 끔 기준) %s  D(gn_rest 기준, 끔 기준) %s'
                  % (r['scene'], r['BE_U'], fmt(r['D_ref_state_vs_off']), fmt(r['D_ref_rest_vs_off'])))
        else:
            print('   1/%-3d %-8s %s  D(·, 끔 기준) %s  D(·, gn_rest 기준) %s  PCG 평균 %.1f'
                  % (r['sps'], r['mode'], 'U' if r['usable'] else r['stop'][:14].ljust(14),
                     fmt(r['D_to_ref_off']), fmt(r['D_to_ref_rest']), r['pcg_iters_mean']))
    c = res['B']['convergence']
    print('B — R · BE_U 1e4 수렴: D(끔, gn_state) ' + ' · '.join('1/%d %s' % (s, fmt(v)) for s, v in c['D_off_vs_gn_state'].items()))
    for m in B_MODES:
        print('    %s 반으로: ' % m + ' · '.join('%s %s' % (k, fmt(v)) for k, v in c['D_%s_halving' % m].items()))
    print('    기준 gn_state 1/%d 까지의 거리' % B_SPS[-1])
    for r in res['B']['D_to_gn_state_6400']:
        print('    %-7s 1/%-3d ' % (r['pcg'], r['sps']) + '  '.join('%s %s' % (m, fmt(r[m], '%.2f')) for m in L4.MODES))
    print('C — 엔진 위상 프로파일러 (400 단계, ms)')
    for m, pr in res['C']['profile_ms_400_steps'].items():
        print('    %-8s 단계 %s ms · 굽힘 %s · PCG 고리 %s · 인장 %s · 전단 %s · PCG %s 회/풀이'
              % (m, fmt(pr.get('ms_per_step'), '%.2f'), fmt(pr.get('bending'), '%.0f'), fmt(pr.get('solve:pcg_loop'), '%.0f'),
                 fmt(pr.get('tension'), '%.0f'), fmt(pr.get('shear'), '%.0f'), fmt(pr.get('pcg_iters_per_solve'), '%.1f')))
    with open(OUT, 'w', encoding='utf-8') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
    print('→ %s (%.0f s)' % (os.path.relpath(OUT, L4.ROOT), res['seconds']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
