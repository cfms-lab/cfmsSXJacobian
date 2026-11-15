"""data/garments/<옷>/source/*.zprj (CLO 프로젝트의 사본)를 풀어 garment.npz · garment.json 으로 저장하고, 옷깃 물성 파일 materials.json 을 만든다.

사용:  python tests/export_clo_garment.py [<옷 폴더 이름> ...] [--testbed <테스트베드 drape 폴더>]
  - CLO 파일은 테스트베드의 ZprjReader 로 읽는다(읽기만; 그 venv 의 python 으로 자식 프로세스를 띄운다). 한 번 풀어 두면 측정 · 그림은
    테스트베드 없이 이 폴더의 자료만 쓴다.
  - garment.npz: 패턴마다 n2_<i>(2 차원 메시, cm) · n3_<i>(CLO 가 시뮬레이션한 3 차원 위치, cm) · el_<i>(삼각형), stitches(패턴A, 노드A, 패턴B, 노드B).
  - garment.json: 출처(파일 이름 · SHA-256 · CLO 판) · 패턴 목록(이름 · 크기 · 배치점 · 원단) · 원단 물성(CLO 계수 그대로) · 옷깃 패턴.
  - materials.json: **사용자가 고치는 파일.** 옷깃 패턴의 인장 · 전단 · 굽힘 강성을 나머지 패턴과 따로 준다. 있으면 덮어쓰지 않는다.
자료는 CLO-SET CONNECT 의 무료 라이선스 항목이라 리포에 넣지 않는다(.gitignore). docs/MEASURE_FROZEN_FORCE.md 「M4b」.
"""
import sys
sys.dont_write_bytecode = True

import glob  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import subprocess  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GARMENTS = os.path.join(ROOT, 'data', 'garments')
TESTBED = os.path.normpath(os.path.join(ROOT, '..', 'cfmsMiindo_drape_improve_dev', 'drape'))
# 옷깃 패턴 이름(CLO 패턴 이름은 번호다). 훑어서 고른 값 — docs/MEASURE_FROZEN_FORCE.md 「M4b」
COLLAR = {'kofoti_024_classic_dress_shirt': '10', 'kofoti_370_button_down_shirt': '8', 'kofoti_244_wide_fit_open_collar_shirt': '6'}
DEFAULT_MATERIAL = {'E': 1.0e6, 'SE': 5.0e4, 'BE': 10.0}        # 인장 · 전단(dyn/cm) · 굽힘(dyn·cm) — 셔츠 b83 실행과 같은 값
COLLAR_BE_SWEEP = [10.0, 100.0, 1000.0, 10000.0]


def child(spec_path):
    import numpy as np
    spec = json.load(open(spec_path, encoding='utf-8'))
    sys.path.insert(0, spec['testbed'])
    from ZprjReader import ZprjReader

    for gdir in spec['garments']:
        hits = glob.glob(os.path.join(gdir, 'source', '*.zprj'))
        assert len(hits) == 1, (gdir, hits)
        r = ZprjReader(hits[0], verbose=False)
        arrays, patterns = {}, []
        for i, p in enumerate(r.patterns):
            n2, n3, el = np.asarray(p.nodes2d, float), np.asarray(p.nodes3d, float), np.asarray(p.elements, np.int32)
            arrays['n2_%d' % i], arrays['n3_%d' % i], arrays['el_%d' % i] = n2, n3, el
            clo = p.info['clo']
            ext = n2.max(axis=0) - n2.min(axis=0)
            patterns.append({'index': i, 'name': str(p.name), 'nodes': int(len(n2)), 'triangles': int(len(el)),
                             'size_cm': [float(ext.max()), float(ext.min())], 'height_mean_cm': float(n3[:, 1].mean()),
                             'arrange_point': clo.get('arrange_point'), 'layer': clo.get('layer'), 'fabric_index': clo.get('fabric_index')})
        arrays['stitches'] = np.asarray(r.stitches, np.int32).reshape(-1, 4)
        np.savez_compressed(os.path.join(gdir, 'garment.npz'), **arrays)
        fabrics = [{'name': f.get('name'), 'physical': {k: (float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else v)
                                                        for k, v in (f.get('physical') or {}).items()}} for f in r.fabrics]
        gid = os.path.basename(gdir)
        meta = {'id': gid, 'source_file': os.path.basename(hits[0]), 'source_sha256': hashlib.sha256(open(hits[0], 'rb').read()).hexdigest(),
                'clo_version': r.clo_version, 'units': 'cm', 'patterns': patterns, 'skipped_patterns': len(r.skipped_patterns),
                'stitch_pairs': int(len(arrays['stitches'])), 'fabrics': fabrics, 'collar_pattern': spec['collar'].get(gid)}
        src = os.path.join(gdir, 'source', 'source.json')
        if os.path.isfile(src):
            meta['source'] = json.load(open(src, encoding='utf-8'))
        json.dump(meta, open(os.path.join(gdir, 'garment.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    return 0


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == '--child':
        return child(sys.argv[2])
    args = sys.argv[1:]
    testbed = TESTBED
    if '--testbed' in args:
        i = args.index('--testbed')
        testbed = os.path.abspath(args[i + 1])
        del args[i:i + 2]
    names = args or sorted(d for d in os.listdir(GARMENTS) if os.path.isdir(os.path.join(GARMENTS, d, 'source')))
    gdirs = [os.path.join(GARMENTS, n) for n in names]
    spec = {'testbed': testbed, 'garments': gdirs, 'collar': COLLAR}
    sp = os.path.join(GARMENTS, '_export_spec.json')
    json.dump(spec, open(sp, 'w', encoding='utf-8'), ensure_ascii=False)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONUTF8='1')
    python = os.path.join(testbed, '.venv', 'Scripts', 'python.exe')
    r = subprocess.run([python, os.path.abspath(__file__), '--child', sp], env=env, capture_output=True, timeout=3600)
    os.remove(sp)
    if r.returncode != 0:
        print('풀기 실패:\n' + (r.stdout.decode('utf-8', 'replace') + r.stderr.decode('utf-8', 'replace'))[-1500:])
        return 1
    for gdir in gdirs:
        meta = json.load(open(os.path.join(gdir, 'garment.json'), encoding='utf-8'))
        collar = meta.get('collar_pattern')
        mpath = os.path.join(gdir, 'materials.json')
        made = False
        if not os.path.isfile(mpath):
            mat = {
                '_읽는 법': '옷깃 패턴의 물성을 나머지 패턴과 따로 준다. 단위: 길이 cm, 인장 E · 전단 SE 는 dyn/cm, 굽힘 BE 는 dyn·cm. '
                         'collar 의 값은 숫자 하나 또는 목록(목록이면 값마다 한 번씩 푼다). 이 파일은 손으로 고친다 — export_clo_garment.py 가 덮어쓰지 않는다.',
                'collar_pattern': collar,
                'default': dict(DEFAULT_MATERIAL),
                'collar': {'E': DEFAULT_MATERIAL['E'], 'SE': DEFAULT_MATERIAL['SE'], 'BE': list(COLLAR_BE_SWEEP)},
            }
            json.dump(mat, open(mpath, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
            made = True
        cp = [p for p in meta['patterns'] if p['name'] == collar]
        print('%s — CLO %s, 패턴 %d · 봉제 쌍 %d · 원단 %d | 옷깃 패턴 「%s」%s | materials.json %s' % (
            meta['id'], meta['clo_version'], len(meta['patterns']), meta['stitch_pairs'], len(meta['fabrics']), collar,
            (' %.1f × %.1f cm, 노드 %d' % (cp[0]['size_cm'][0], cp[0]['size_cm'][1], cp[0]['nodes'])) if cp else ' (없음!)',
            '새로 만듦' if made else '그대로 둠'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
