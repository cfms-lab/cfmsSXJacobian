"""측정(판정 밖): 옷깃의 비선형 평형 — 정확한 굽힘 힘과 길이 고정 굽힘 힘으로 따로 평형시킨 차 (docs/MEASURE_FROZEN_FORCE.md 「M4c」).

사용:  python tests/measure_frozen_garment_nl.py [<옷 폴더 이름> ...]
자료 · 옷깃 물성은 M4b(tests/measure_frozen_garment_clo.py)와 같다 — data/garments/<옷>/ 의 garment.npz · garment.json · materials.json.
tests/measure_frozen_garment_nl.cpp 의 실행 파일(%TEMP%\\cfmsSXJacobian_build\\measure_frozen_garment_nl.exe)이 푼다:
  A. 옷깃만 다시 드레이프(봉제 노드 고정 + 중력, 접촉 없음): 정확한 힘의 평형 → 거기서 길이 고정 힘의 평형, 그 차와 선형 응답 둘
  B. CLO 모양에서의 죽은 하중 문제(보조)
전제: 새 프로그램이 CLO 모양에서 낸 준정부호 선형 응답이 M4b 의 기록(results/measure_frozen_garment_clo.json)과 상대 1e-6 안에서 같다.
결과: results/measure_frozen_garment_nl.json. 전제가 서면 0, 깨지면 1, 준비 안 됨 2. 엔진 · 테스트베드는 쓰지 않는다.
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
import measure_frozen_garment_clo as M  # noqa: E402

ROOT = M.ROOT
TMP = os.environ.get('TEMP') or tempfile.gettempdir()
WORK = os.path.join(TMP, 'cfmsSXJ_frozen_garment_nl')
EXE = os.path.join(TMP, 'cfmsSXJacobian_build', 'measure_frozen_garment_nl.exe')
REC = os.path.join(ROOT, 'results', 'measure_frozen_garment_clo.json')
OUT = os.path.join(ROOT, 'results', 'measure_frozen_garment_nl.json')
DENSITY = 0.02          # g/cm² — 주름진 천 · 드레이프 시험과 같은 면밀도
GRAVITY = 980.0         # cm/s², −y
LIM_CHECK = 1e-6


def run(gid):
    gdir, meta, mat, data = M.load_garment(gid)
    i, name, n2, n3, el, pinned = M.collar_of(meta, mat, data)
    collar = dict(mat['default'])
    collar.update(mat.get('collar') or {})
    E, SE, bes = float(collar['E']), float(collar['SE']), M.as_list(collar['BE'])
    inp = os.path.join(WORK, gid + '.txt')
    res = os.path.join(WORK, gid + '.json')
    M.write_input(inp, E, SE, n2, n3, el, pinned)
    t0 = time.perf_counter()
    r = subprocess.run([EXE, inp, res, '%.17g' % DENSITY, '%.17g' % GRAVITY] + ['%.17g' % b for b in bes], capture_output=True)
    log = r.stdout.decode('utf-8', 'replace')
    if r.returncode != 0 or not os.path.isfile(res):
        raise RuntimeError('옷깃 풀이 실패 (%s): %s %s' % (gid, log[-600:], r.stderr.decode('utf-8', 'replace')[-600:]))
    out = json.load(open(res, encoding='utf-8'))
    return {'garment': gid, 'pattern': name, 'pattern_index': i, 'nodes': int(len(n2)), 'triangles': int(len(el)), 'pinned': int(pinned.sum()),
            'collar_material': {'E': E, 'SE': SE, 'BE': bes}, 'seconds': time.perf_counter() - t0, 'result': out, 'log': log}


def length(cm):
    return ('%.3g mm' % (cm * 10.0)) if cm >= 0.01 else ('%.3g µm' % (cm * 1.0e4))


def main():
    names = sys.argv[1:] or sorted(d for d in os.listdir(M.GARMENTS) if os.path.isfile(os.path.join(M.GARMENTS, d, 'garment.npz')))
    if not (os.path.isfile(EXE) and os.path.isfile(REC)):
        print('없는 것:', EXE if not os.path.isfile(EXE) else REC)
        return 2
    if os.path.isdir(WORK):
        shutil.rmtree(WORK)
    os.makedirs(WORK)
    rec = {g['garment']: {c['BE']: c for c in g['results']['collar']['cases']} for g in json.load(open(REC, encoding='utf-8'))['garments']}
    t0 = time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(names)) as ex:
        garments = list(ex.map(run, names))
    bad = []
    for g in garments:
        r = g['result']
        print('%s — 옷깃 패턴 「%s」, 노드 %d(봉제 %d) · 삼각형 %d · 반띠폭 %d (%.0f s)' % (
            g['garment'], g['pattern'], g['nodes'], g['pinned'], g['triangles'], r['half_bandwidth'], g['seconds']))
        for c in r['cases']:
            a, b = c['A'], c['B']
            want = rec[g['garment']][c['BE']]['u_max_cm']
            rel = abs(b['u_psd_at_x0_max_cm'] - want) / want
            c['check_vs_m4b'] = {'u_max_cm_m4b': want, 'relative_difference': rel}
            if not rel <= LIM_CHECK:
                bad.append('%s|%g' % (g['garment'], c['BE']))
            nl = a['u_nl_max_cm']
            print('  굽힘 %-6g A: x_E %s(반복 %d, CLO 에서 %.2f cm, 각 rms %.1f°, 참 접선 %s) | u_NL %s%s · Newton 첫 걸음 %s · 준정부호 선형 %s → NL/T %.3f · NL/L %.3f · cos %.4f | 보정 %.1e cm'
                  ' | M4b(CLO 모양의 선형) %s · 확인 %.1e | B: x0 참 접선 %s, 길이 고정 %s%s' % (
                      c['BE'], '수렴' if a['exact']['converged'] else '미수렴', a['exact']['iterations'], a['moved_from_clo_max_cm'], a['theta_rms_deg'],
                      '정부호' if a['true_tangent_pd_at_xE'] else '부정치', length(nl), '' if a['frozen']['converged'] else ' [미수렴]',
                      length(a['u_true_tangent_max_cm']), length(a['u_psd_max_cm']),
                      nl / a['u_true_tangent_max_cm'] if a['u_true_tangent_max_cm'] else float('nan'),
                      nl / a['u_psd_max_cm'] if a['u_psd_max_cm'] else float('nan'), a['cos_nl_true'], a['corrected']['u_max_cm'],
                      length(want), rel, '정부호' if b['true_tangent_pd_at_x0'] else '부정치',
                      '수렴' if b['frozen']['converged'] else '미수렴', (' u %s' % length(b['frozen']['u_max_cm'])) if b['frozen']['converged'] else ''))
    valid = not bad
    print('전제(CLO 모양의 준정부호 선형 응답 = M4b 기록, 상대 %.0e): %s' % (LIM_CHECK, 'PASS' if valid else 'FAIL %s' % bad))
    for g in garments:
        g.pop('log')
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump({'data': 'data/garments', 'density_g_cm2': DENSITY, 'gravity_cm_s2': GRAVITY, 'premise': {'matches_m4b': valid, 'different': bad},
                   'valid': valid, 'garments': garments, 'seconds': time.perf_counter() - t0}, fh, indent=1, ensure_ascii=False)
        fh.write('\n')
    print('→ %s (%.0f s)' % (os.path.relpath(OUT, ROOT), time.perf_counter() - t0))
    return 0 if valid else 1


if __name__ == '__main__':
    sys.exit(main())
