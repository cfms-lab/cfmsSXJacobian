"""측정(판정 밖): 실제 옷의 옷깃에서 길이 고정 근사의 가짜 힘이 만드는 변위 (docs/MEASURE_FROZEN_FORCE.md 「M4」).

사용:  python tests/measure_frozen_garment.py [--testbed <테스트베드 drape 폴더>] [--archive <옷 zip 이름>] [--dll <엔진 DLL>] [--steps 2000]
  1. 테스트베드 로더(SCADReader · main.build_sim)로 옷을 엔진에 올려 드레이프시킨다(그 venv 의 python 으로 자식 프로세스를 띄운다).
     옷깃 패널만 굽힘 강성을 패턴 파일의 값(10)으로 둔다. 테스트베드 · 엔진 리포는 읽기만 하고, DLL 은 사본(%TEMP%\\cfmsSXJ_frozen_garment)을 쓴다.
  2. 드레이프된 옷깃 패널(위치 · 패턴 좌표 · 봉제 노드 · 삼각형)을 글 파일로 뽑는다.
  3. tests/measure_frozen_garment.cpp 를 빌드한 실행 파일(%TEMP%\\cfmsSXJacobian_build\\measure_frozen_garment.exe)이 가짜 힘과 선형 응답을 푼다.
결과는 results/measure_frozen_garment.json.
"""
import sys
sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP = os.environ.get('TEMP') or tempfile.gettempdir()
WORK = os.path.join(TMP, 'cfmsSXJ_frozen_garment')
TESTBED = os.path.normpath(os.path.join(ROOT, '..', 'cfmsMiindo_drape_improve_dev', 'drape'))
ARCHIVE = '봉제연b83(A)_NBS_neckfix.zip'
DLL = os.path.join(TMP, 'cfmsDrape_build_sxjlink', 'Release', 'cfmsDrape.dll')
EXE = os.path.join(TMP, 'cfmsSXJacobian_build', 'measure_frozen_garment.exe')
COLLAR_NAME = 'Collar_UP'
COLLAR_BE = 10.0                      # 패턴 파일(b38_collar.2dd)의 재질 줄에서 읽은 옷깃 굽힘 강성
BE_CASES = (10.0, 1000.0, 10000.0)    # 파일의 값 · 심지 든 옷깃 수준 둘
BLOCK = 100
SNAPSHOTS = (0.5, 0.75, 1.0)             # 드레이프가 접촉에서 떨며 멈추지 않으므로 세 시각의 모양에서 따로 푼다 (끝 걸음 수의 배)


def child(spec_path):
    import numpy as np
    spec = json.load(open(spec_path, encoding='utf-8'))
    os.environ['CFMS_DRAPE_DLL'] = spec['dll']
    sys.path.insert(0, spec['testbed'])
    from main import build_sim, scad_pattern_names
    from SCADReader import SCADReader

    reader = SCADReader(os.path.join(spec['testbed'], 'MeshData', spec['archive']))
    sim, info, body = build_sim(reader, 'triangle')
    names = list(scad_pattern_names(reader))
    counts = [int(c) for c in info['panel_node_counts']]
    ecounts = [int(c) for c in info['panel_elem_counts']]
    pi = names.index(COLLAR_NAME)
    n0, e0 = sum(counts[:pi]), sum(ecounts[:pi])
    nn, ne = counts[pi], ecounts[pi]
    sim.set_material_panel(pi, BE_U=COLLAR_BE, BE_V=COLLAR_BE)
    out = {'build_info': sim.build_info, 'threads': sim.openmp_threads, 'panel_names': names, 'panel_node_counts': counts,
           'collar_panel': pi, 'global_material': [float(v) for v in sim._material], 'collar_BE_U': COLLAR_BE,
           'body': body is not None, 'verts': int(np.asarray(info['verts']).shape[0])}
    uv = np.asarray(info['uv'], float)
    F = np.asarray(info['faces'], int)[e0:e0 + ne] - n0
    assert F.min() >= 0 and F.max() < nn, '옷깃 삼각형이 옷깃 노드 범위를 벗어난다'
    sew = np.asarray(info['sewing'], int).reshape(-1, 2)
    pinned = np.zeros(nn, int)
    for a, b in sew:
        for g in (a, b):
            if n0 <= g < n0 + nn:
                pinned[g - n0] = 1
    out['collar_nodes'] = nn
    out['collar_pinned'] = int(pinned.sum())
    E, SE = float(sim._material[1]), float(sim._material[6])
    marks = sorted(set(max(1, int(round(f * spec['steps'] / BLOCK))) for f in SNAPSHOTS))
    t0 = time.perf_counter()
    speed, exports = [], []
    for blk in range(1, spec['steps'] // BLOCK + 1):
        sim.step(BLOCK)
        v = sim.velocities()
        speed.append([float(np.abs(v).max()), float(np.abs(v[n0:n0 + nn]).max())])
        if blk in marks:
            P = sim.positions()
            assert np.isfinite(P).all()
            path = spec['export'] % (blk * BLOCK)
            with open(path, 'w', encoding='ascii') as fh:
                fh.write('%.17g %.17g\n%d %d\n' % (E, SE, nn, ne))
                for i in range(nn):
                    p, q = P[n0 + i], uv[n0 + i]
                    fh.write('%.17g %.17g %.17g %.17g %.17g %d\n' % (p[0], p[1], p[2], q[0], q[1], pinned[i]))
                for a, b, c in F:
                    fh.write('%d %d %d\n' % (a, b, c))
            exports.append({'step': blk * BLOCK, 'path': path, 'y_min': float(P[:, 1].min()), 'y_max': float(P[:, 1].max()),
                            'collar_y_mean': float(P[n0:n0 + nn, 1].mean())})
    out['seconds'] = time.perf_counter() - t0
    out['max_speed_all_collar_per_block'] = speed
    out['exports'] = exports
    sim.close()
    json.dump(out, open(spec['out'], 'w', encoding='utf-8'), ensure_ascii=False)
    return 0


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == '--child':
        return child(sys.argv[2])
    args = sys.argv[1:]
    opt = {}
    for flag in ('--testbed', '--archive', '--dll', '--steps'):
        if flag in args:
            i = args.index(flag)
            opt[flag] = args[i + 1]
            del args[i:i + 2]
    testbed = os.path.abspath(opt.get('--testbed', TESTBED))
    archive = opt.get('--archive', ARCHIVE)
    steps = int(opt.get('--steps', 2000))
    if not os.path.isfile(EXE):
        print('먼저 tests/measure_frozen_garment.cpp 를 빌드한다 (docs/MEASURE_FROZEN_FORCE.md 「재현」): %s 없음' % EXE)
        return 1
    if os.path.isdir(WORK):
        shutil.rmtree(WORK)
    os.makedirs(WORK)
    dll = os.path.join(WORK, 'cfmsDrape.dll')
    shutil.copyfile(os.path.abspath(opt.get('--dll', DLL)), dll)
    spec = {'testbed': testbed, 'archive': archive, 'dll': dll, 'steps': steps, 'export': os.path.join(WORK, 'collar_%d.txt'),
            'out': os.path.join(WORK, 'drape.json')}
    sp = os.path.join(WORK, 'spec.json')
    json.dump(spec, open(sp, 'w', encoding='utf-8'), ensure_ascii=False)
    env = {k: v for k, v in os.environ.items()
           if k not in ('CFMS_CPU_SXJ_BEND', 'CFMS_CPU_SXJ_MEMBRANE', 'CFMS_DRAPE_DLL', 'CFMS_PCG_PARALLEL', 'OMP_NUM_THREADS')}
    env.update({'OMP_NUM_THREADS': '1', 'PYTHONDONTWRITEBYTECODE': '1', 'PYTHONUTF8': '1'})
    python = os.path.join(testbed, '.venv', 'Scripts', 'python.exe')
    r = subprocess.run([python, os.path.abspath(__file__), '--child', sp], env=env, cwd=WORK, capture_output=True, timeout=3600)
    if r.returncode != 0 or not os.path.isfile(spec['out']):
        print('드레이프 실패:\n' + (r.stdout.decode('utf-8', 'replace') + r.stderr.decode('utf-8', 'replace'))[-1500:])
        return 1
    drape = json.load(open(spec['out'], encoding='utf-8'))
    sp_all = drape['max_speed_all_collar_per_block']
    print('옷 %s — 엔진 %s(dirty %s), 단일 스레드, 노드 %d, %d 걸음 %.0f s. 옷깃 = 패널 %d(%s), 노드 %d 중 봉제 %d' % (
        archive, drape['build_info'].get('git_commit'), drape['build_info'].get('dirty'), drape['verts'], steps, drape['seconds'],
        drape['collar_panel'], COLLAR_NAME, drape['collar_nodes'], drape['collar_pinned']))
    print('  최대 속력 cm/s (옷 전체 · 옷깃): 처음 묶음 %.2f · %.2f → 마지막 묶음 %.3f · %.3f' % (sp_all[0][0], sp_all[0][1], sp_all[-1][0], sp_all[-1][1]))
    snaps = []
    for ex in drape['exports']:
        res = os.path.join(WORK, 'collar_%d.json' % ex['step'])
        r = subprocess.run([EXE, ex['path'], res] + ['%g' % b for b in BE_CASES], capture_output=True)
        print(' 걸음 %d (옷 전체 y %.1f – %.1f cm, 옷깃 평균 높이 %.1f cm)' % (ex['step'], ex['y_min'], ex['y_max'], ex['collar_y_mean']))
        sys.stdout.write(r.stdout.decode('utf-8', 'replace'))
        if r.returncode != 0 or not os.path.isfile(res):
            print('옷깃 풀이 실패:\n' + r.stderr.decode('utf-8', 'replace')[-800:])
            return 1
        snaps.append({'step': ex['step'], 'y_min': ex['y_min'], 'y_max': ex['y_max'], 'collar_y_mean': ex['collar_y_mean'],
                      'collar': json.load(open(res, encoding='utf-8'))})
    for e in drape['exports']:
        e.pop('path', None)
    out = {'archive': archive, 'testbed': testbed.replace('\\', '/'), 'dll_sha256': hashlib.sha256(open(dll, 'rb').read()).hexdigest(),
           'steps': steps, 'timestep': 0.005, 'collar_name': COLLAR_NAME, 'be_cases': list(BE_CASES), 'drape': drape, 'snapshots': snaps}
    with open(os.path.join(ROOT, 'results', 'measure_frozen_garment.json'), 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print('→ results/measure_frozen_garment.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
