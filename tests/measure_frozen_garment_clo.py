"""측정(판정 밖): CLO 가 드레이프시킨 옷깃에서 길이 고정 근사의 가짜 힘이 만드는 변위 (docs/MEASURE_FROZEN_FORCE.md 「M4b」).

사용:  python tests/measure_frozen_garment_clo.py [<옷 폴더 이름> ...] [--collar-be 10,100,1000] [--collar-e 1e6] [--collar-se 5e4]
자료는 이 리포의 data/garments/<옷>/ 에 있다(tests/export_clo_garment.py 가 CLO 프로젝트에서 풀어 둔 garment.npz · garment.json).
**옷깃의 물성은 그 폴더의 materials.json 에서 읽는다** — collar 의 E(인장) · SE(전단) · BE(굽힘)를 나머지 패턴과 따로 준다. BE 가 목록이면
값마다 한 번씩 푼다. 명령줄의 --collar-* 는 그 파일의 값을 이번 실행에서만 덮는다.
tests/measure_frozen_garment.cpp 의 실행 파일(%TEMP%\\cfmsSXJacobian_build\\measure_frozen_garment.exe)이 옷깃에서 가짜 힘과 선형 응답을 푼다.
결과: results/measure_frozen_garment_clo.json(숫자) · data/garments/<옷>/collar_response.npz(노드별 변위 · 변 길이 변화 — 그림용, 리포에 넣지 않는다).
엔진 · 테스트베드는 쓰지 않는다.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402

import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GARMENTS = os.path.join(ROOT, 'data', 'garments')
TMP = os.environ.get('TEMP') or tempfile.gettempdir()
WORK = os.path.join(TMP, 'cfmsSXJ_frozen_garment_clo')
EXE = os.path.join(TMP, 'cfmsSXJacobian_build', 'measure_frozen_garment.exe')


def as_list(v):
    return [float(x) for x in v] if isinstance(v, (list, tuple)) else [float(v)]


def load_garment(gid):
    gdir = os.path.join(GARMENTS, gid)
    meta = json.load(open(os.path.join(gdir, 'garment.json'), encoding='utf-8'))
    mat = json.load(open(os.path.join(gdir, 'materials.json'), encoding='utf-8'))
    data = np.load(os.path.join(gdir, 'garment.npz'))
    return gdir, meta, mat, data


def collar_of(meta, mat, data):
    name = str(mat.get('collar_pattern') or meta.get('collar_pattern'))
    idx = [p['index'] for p in meta['patterns'] if p['name'] == name]
    assert len(idx) == 1, ('옷깃 패턴을 찾을 수 없다', name)
    i = idx[0]
    n2, n3, el = data['n2_%d' % i], data['n3_%d' % i], data['el_%d' % i]
    pinned = np.zeros(len(n2), int)
    for pa, na, pb, nb in data['stitches']:
        if pa == i:
            pinned[na] = 1
        if pb == i:
            pinned[nb] = 1
    return i, name, n2, n3, el, pinned


def write_input(path, E, SE, n2, n3, el, pinned):
    with open(path, 'w', encoding='ascii') as fh:
        fh.write('%.17g %.17g\n%d %d\n' % (E, SE, len(n2), len(el)))
        for k in range(len(n2)):
            fh.write('%.17g %.17g %.17g %.17g %.17g %d\n' % (n3[k, 0], n3[k, 1], n3[k, 2], n2[k, 0], n2[k, 1], pinned[k]))
        for a, b, c in el:
            fh.write('%d %d %d\n' % (a, b, c))


def solve(tag, E, SE, bes, n2, n3, el, pinned, dump):
    inp = os.path.join(WORK, tag + '.txt')
    res = os.path.join(WORK, tag + '.json')
    write_input(inp, E, SE, n2, n3, el, pinned)
    cmd = [EXE, inp, res] + (['--dump', os.path.join(WORK, tag)] if dump else []) + ['%.17g' % b for b in bes]
    r = subprocess.run(cmd, capture_output=True)
    sys.stdout.write(r.stdout.decode('utf-8', 'replace'))
    if r.returncode != 0 or not os.path.isfile(res):
        raise RuntimeError('옷깃 풀이 실패: ' + r.stderr.decode('utf-8', 'replace')[-600:])
    out = json.load(open(res, encoding='utf-8'))
    fields = [np.loadtxt(os.path.join(WORK, '%s_%d.txt' % (tag, k))).reshape(len(n2), 6) for k in range(len(bes))] if dump else []
    return out, fields


def main():
    args = sys.argv[1:]
    opt = {}
    for flag in ('--collar-be', '--collar-e', '--collar-se'):
        if flag in args:
            i = args.index(flag)
            opt[flag] = args[i + 1]
            del args[i:i + 2]
    names = args or sorted(d for d in os.listdir(GARMENTS) if os.path.isfile(os.path.join(GARMENTS, d, 'garment.npz')))
    if not os.path.isfile(EXE):
        print('먼저 tests/measure_frozen_garment.cpp 를 빌드한다: %s 없음' % EXE)
        return 1
    if os.path.isdir(WORK):
        shutil.rmtree(WORK)
    os.makedirs(WORK)
    garments = []
    for gid in names:
        gdir, meta, mat, data = load_garment(gid)
        i, name, n2, n3, el, pinned = collar_of(meta, mat, data)
        collar = dict(mat['default'])
        collar.update(mat.get('collar') or {})
        E = float(opt.get('--collar-e', collar['E']))
        SE = float(opt.get('--collar-se', collar['SE']))
        bes = [float(x) for x in opt['--collar-be'].split(',')] if '--collar-be' in opt else as_list(collar['BE'])
        ext = n2.max(axis=0) - n2.min(axis=0)
        print('%s — CLO %s, 옷깃 패턴 「%s」 %.1f × %.1f cm, 노드 %d(봉제 %d) · 삼각형 %d' % (
            gid, meta['clo_version'], name, ext.max(), ext.min(), len(n2), int(pinned.sum()), len(el)))
        print(' ① 옷깃 물성(materials.json%s): 인장 %g · 전단 %g · 굽힘 %s' % (
            ' + 명령줄' if opt else '', E, SE, ', '.join('%g' % b for b in bes)))
        res, fields = solve(gid + '_collar', E, SE, bes, n2, n3, el, pinned, True)
        np.savez_compressed(os.path.join(gdir, 'collar_response.npz'), pattern_index=i, E=E, SE=SE, BE=np.asarray(bes),
                            u=np.asarray([f[:, 0:3] for f in fields]), strain_mean_abs=np.asarray([f[:, 3] for f in fields]),
                            strain_mean_signed=np.asarray([f[:, 4] for f in fields]), strain_max_abs=np.asarray([f[:, 5] for f in fields]),
                            pinned=pinned)
        g = {'garment': gid, 'source_file': meta['source_file'], 'clo_version': meta['clo_version'], 'pattern': name, 'pattern_index': i,
             'size_cm': [float(ext.max()), float(ext.min())], 'nodes': int(len(n2)), 'triangles': int(len(el)), 'pinned': int(pinned.sum()),
             'collar_material': {'E': E, 'SE': SE, 'BE': bes}, 'results': {'collar': res}}
        # 참고: CLO 파일의 원단 물성을 단위만 환산 (g/s² = dyn/cm, g·mm²/s² = 0.01 dyn·cm). 뜻이 같다는 근거는 없다
        fi = [p['fabric_index'] for p in meta['patterns'] if p['index'] == i][0]
        phys = meta['fabrics'][fi]['physical'] if isinstance(fi, int) and 0 <= fi < len(meta['fabrics']) else {}
        if phys:
            e_clo = 0.5 * (float(phys.get('fSuK', 0.)) + float(phys.get('fSvK', 0.)))
            se_clo = float(phys.get('fLeftShearK', 0.))
            be_clo = 0.5 * (float(phys.get('fBuK', 0.)) + float(phys.get('fBvK', 0.))) * 0.01
            print(' ② 참고 — CLO 물성을 단위만 환산: 인장 %g · 전단 %g · 굽힘 %g' % (e_clo, se_clo, be_clo))
            res2, _ = solve(gid + '_clo', e_clo, se_clo, [be_clo], n2, n3, el, pinned, False)
            g['clo_converted'] = {'E': e_clo, 'SE': se_clo, 'BE': be_clo}
            g['results']['clo'] = res2
        garments.append(g)
    with open(os.path.join(ROOT, 'results', 'measure_frozen_garment_clo.json'), 'w', encoding='utf-8') as fh:
        json.dump({'data': 'data/garments', 'garments': garments}, fh, indent=1, ensure_ascii=False)
    print('→ results/measure_frozen_garment_clo.json · data/garments/<옷>/collar_response.npz')
    return 0


if __name__ == '__main__':
    sys.exit(main())
