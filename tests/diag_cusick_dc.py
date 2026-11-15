"""판정 밖 진단: 엔진 합성 시험에서 새 기본값(N)과 킬스위치(K)의 Cusick 드레이프 계수가 다른 까닭 (docs/GATE_ENGINE_SUITE.md 「판정 뒤 진단 6」).

사용:  python tests/diag_cusick_dc.py [<엔진 DLL>] [--engine <엔진 리포>] [--jobs 8]
엔진 리포의 시험 `drape/scripts/py_cusick_test.py` 의 장면을 그대로 쓴다(시편 반경 15 cm · 받침 9 cm · 노드 1009 · 자가 충돌 켬 · 단일 스레드 ·
동적 10 초 뒤 정적 모드로 안착). 메시와 드레이프 계수는 그 스크립트의 함수를 불러 쓴다(엔진 코드를 이 리포에 싣지 않는다). 엔진 리포는 읽기만
하고, DLL 은 사본(%TEMP%\\cfmsSXJ_diag_cusick)을 쓴다. 굽힘 강성은 합성 시험에서 값이 갈린 BE_U 50 하나다.

묻는 것 (돌리기 전에 적었다 — 판정이 아니라 진단이라 통과 기준은 없다)
  Q1  어느 스위치가 값을 바꾸나 — N(둘 다 켬) · K(둘 다 끔) · 굽힘만 끔 · 인장 · 전단만 끔.
  Q2  두 끝 상태가 같은 힘의 평형인가 — 끝 위치에서 내부 힘 + 중력의 잔차. 두 커널의 힘이 그 위치에서 같은가.
  Q3  끝 상태를 바꿔 넣으면 머무는가 — N 의 끝 위치에서 K 를 정적 모드로 이어 돌리고, 그 반대도. 머물면 평형이 둘이고 어느 쪽에 닿느냐는 길의 문제다.
  Q4  반올림 수준의 섭동에 값이 얼마나 흔들리나 — 처음 위치에 1e-10 · 1e-6 cm 의 난수를 더한다(팔마다 씨 셋).
      한 팔 안의 흔들림이 두 팔의 차이만큼 크면, 차이는 커널이 체계적으로 만든 것이 아니다.
  Q5  「안착」 판정 뒤에도 값이 기어가나 — 안착 뒤 정적 2000 걸음을 더 밟는다.
  Q6  동적 단계의 시간 간격을 줄이면 두 팔이 모이나 — 1/2 · 1/4 (동적 10 초는 그대로).
      처음 실행에서 1/4 에 두 팔이 모이는 것을 본 뒤 더했다: 1/8, 그리고 1/4 에서의 섭동(1e-10 cm · 씨 둘) — 모인 값이 굳은가.
결과는 results/diag_cusick_dc.json.
"""
import sys
sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ThreadPoolExecutor  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP = os.environ.get('TEMP') or tempfile.gettempdir()
WORK = os.path.join(TMP, 'cfmsSXJ_diag_cusick')
ENGINE_DEFAULT = os.path.normpath(os.path.join(ROOT, '..', 'cfmsDrape_dev'))
DLL_DEFAULT = os.path.join(TMP, 'cfmsDrape_build_sxjlink', 'Release', 'cfmsDrape.dll')
ARMS = {
    'N': {},
    'K': {'CFMS_CPU_SXJ_BEND': 'off', 'CFMS_CPU_SXJ_MEMBRANE': 'off'},
    'bend_off': {'CFMS_CPU_SXJ_BEND': 'off'},
    'membrane_off': {'CFMS_CPU_SXJ_MEMBRANE': 'off'},
}
BE = 50.0
R_S, R_D, RHO_A = 15.0, 9.0, 0.02
DT0, DYN_SECONDS = 0.005, 10.0
BLOCK, MAX_BLOCKS, MIN_BLOCKS, SETTLE = 50, 80, 6, 2e-3      # 엔진 시험의 안착 규칙
EXTRA_BLOCKS = 40                                             # Q5: 안착 뒤 정적 2000 걸음


# ----------------------------------------------------------------------------------------------------------------- 자식
def child(spec_path):
    import numpy as np
    spec = json.load(open(spec_path, encoding='utf-8'))
    sys.path.insert(0, os.path.join(spec['engine'], 'drape', 'scripts'))
    sys.path.insert(0, os.path.join(spec['engine'], 'drape'))
    import py_cusick_test as cus                 # 엔진 시험의 메시 · 지표 함수 (OMP_NUM_THREADS 기본 1 도 여기서 건다)
    from cfmsdrape import Drape

    V, F, ring_idx = cus.polar_disk(R_S)
    nv = len(V)
    edge_idx = ring_idx[-1]
    m = RHO_A * math.pi * R_S ** 2 / nv
    pins = [0]
    for ring in ring_idx:
        for i in ring:
            if math.hypot(V[i, 0], V[i, 2]) <= R_D + 1e-6:
                pins.append(i)
    free = np.ones(nv, bool)
    free[list(set(pins))] = False
    V3 = V.copy()
    if spec.get('amp', 0.0) > 0.0:
        rng = np.random.default_rng(spec['seed'])
        V3[free] += spec['amp'] * rng.uniform(-1.0, 1.0, (int(free.sum()), 3))

    d = Drape(spec['dll'])
    out = {'id': spec['id'], 'arm': spec['arm'], 'build_info': d.build_info, 'threads': d.openmp_threads, 'nodes': nv,
           'free_nodes': int(free.sum()), 'mass_per_node': m}
    d.set_mesh(V3, F)
    d.set_mesh_2d(V[:, [0, 2]].copy())
    for p in set(pins):
        d.pin_vertex(p)
    d.set_material(mass=m, E=1e6, BE_U=BE, SE=5e4)
    d.set_gravity((0.0, -cus.G_ACC, 0.0))
    d.set_self_collision(True)
    d.set_timestep(spec['dt'])

    def metrics(P):
        dc, nodes, area = cus.drape_metrics(P, F, edge_idx, R_S, R_D)
        return {'dc': dc, 'nodes': int(nodes), 'area': area, 'y_min': float(P[:, 1].min()), 'y_edge_mean': float(P[edge_idx, 1].mean())}

    if spec['mode'] == 'force':                  # Q2: 주어진 위치에서 이 팔의 커널이 내는 내부 힘
        # 걸음을 밟지 않고 묻는다 — 평평한 상태에서 한 걸음 밟은 뒤 먼 위치의 힘을 물으면 평형 위치에서도 잔차가 무게의 수천 배가 된다
        # (엔진에 남은 속도의 감쇠 힘이 섞이는 것으로 보인다 — 확인하지 않았다). 처음 실행에서 그렇게 물었다가 두 커널의 힘이 18 % 다르다는
        # 틀린 값을 얻었다 — 하네스 잘못.
        for src, dst in spec['force_jobs']:
            np.save(dst, d.debug_internal_force(np.load(src)))
        out['force_jobs'] = len(spec['force_jobs'])
    else:
        t0 = time.perf_counter()
        if spec.get('start'):                    # Q3: 다른 팔의 끝 위치에서 정적 모드로 잇는다
            P0 = np.load(spec['start'])
            d.set_positions(P0)
            out['start'] = metrics(P0)
        else:
            d.step(int(round(DYN_SECONDS / spec['dt'])))
            P = d.positions()
            out['after_dynamic'] = dict(metrics(P), v_max=float(np.abs(d.velocities()).max()))
        d.set_static_mode(True)
        prev = d.positions().copy()
        blocks, settled, cur = [], None, prev
        for k in range(MAX_BLOCKS):
            d.step(BLOCK)
            cur = d.positions()
            if not np.isfinite(cur).all():
                out['diverged'] = True
                break
            delta = float(np.abs(cur - prev).max())
            blocks.append(delta)
            if k + 1 >= MIN_BLOCKS and delta < SETTLE:
                settled = k + 1
                break
            prev = cur.copy()
        out['static_delta'] = blocks
        out['settled_blocks'] = settled
        out['at_settle'] = metrics(cur)
        extra = []
        prev = cur.copy()
        for j in range(EXTRA_BLOCKS):
            d.step(BLOCK)
            cur = d.positions()
            delta = float(np.abs(cur - prev).max())
            prev = cur.copy()
            if (j + 1) % 5 == 0:
                extra.append(dict(metrics(cur), block=j + 1, delta=delta))
        out['extra'] = extra
        out['final'] = metrics(cur)
        f = d.debug_internal_force(cur)
        r = f.copy()
        r[:, 1] -= m * cus.G_ACC
        rn = np.linalg.norm(r[free], axis=1) / (m * cus.G_ACC)      # 노드 무게 단위
        out['residual_over_weight'] = {'max': float(rn.max()), 'rms': float(np.sqrt((rn ** 2).mean())), 'median': float(np.median(rn)),
                                       'nodes_over_0.05': int((rn > 0.05).sum()), 'nodes_over_1': int((rn > 1.0).sum())}
        out['gravity_pe_over_weight_cm'] = float(cur[free, 1].sum())   # Σ y (자유 노드) — 작을수록 더 처졌다
        np.save(spec['npy'], cur)
        out['seconds'] = time.perf_counter() - t0
    d.close()
    json.dump(out, open(spec['out'], 'w', encoding='utf-8'), ensure_ascii=False)
    return 0


# ----------------------------------------------------------------------------------------------------------------- 부모
def launch(spec, python):
    p = os.path.join(WORK, spec['id'] + '.spec.json')
    spec['out'] = os.path.join(WORK, spec['id'] + '.out.json')
    spec['npy'] = os.path.join(WORK, spec['id'] + '.npy')
    json.dump(spec, open(p, 'w', encoding='utf-8'))
    env = {k: v for k, v in os.environ.items()
           if k not in ('CFMS_CPU_SXJ_BEND', 'CFMS_CPU_SXJ_MEMBRANE', 'CFMS_DRAPE_DLL', 'CFMS_PCG_PARALLEL', 'OMP_NUM_THREADS')}
    env.update(ARMS[spec['arm']])
    env.update({'OMP_NUM_THREADS': '1', 'PYTHONDONTWRITEBYTECODE': '1', 'PYTHONUTF8': '1'})
    r = subprocess.run([python, os.path.abspath(__file__), '--child', p], env=env, capture_output=True, timeout=3600)
    if r.returncode != 0 or not os.path.isfile(spec['out']):
        return {'id': spec['id'], 'arm': spec['arm'], 'error': (r.stderr.decode('utf-8', 'replace') + r.stdout.decode('utf-8', 'replace'))[-600:]}
    return json.load(open(spec['out'], encoding='utf-8'))


def row(r):
    if 'error' in r:
        return '%-22s 오류 %s' % (r['id'], r['error'][-200:].replace('\n', ' | '))
    a, s, f = r.get('after_dynamic') or r.get('start'), r['at_settle'], r['final']
    res = r['residual_over_weight']
    return '%-22s %-12s 처음 %5.1f %%  안착 %5.1f %% (%s 묶음)  +2000 걸음 %5.1f %%  노드 %d  잔차/무게 중앙 %.1e · 최대 %.1e · >0.05 인 노드 %d' % (
        r['id'], r['arm'], 100 * a['dc'], 100 * s['dc'], r['settled_blocks'], 100 * f['dc'], f['nodes'], res['median'], res['max'],
        res['nodes_over_0.05'])


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == '--child':
        return child(sys.argv[2])
    import numpy as np
    args = sys.argv[1:]
    opt = {}
    for flag in ('--engine', '--jobs'):
        if flag in args:
            i = args.index(flag)
            opt[flag] = args[i + 1]
            del args[i:i + 2]
    engine = os.path.abspath(opt.get('--engine', ENGINE_DEFAULT))
    jobs = int(opt.get('--jobs', 8))
    dll_src = os.path.abspath(args[0]) if args else DLL_DEFAULT
    python = os.path.join(engine, 'drape', '.venv', 'Scripts', 'python.exe')
    if os.path.isdir(WORK):
        shutil.rmtree(WORK)
    os.makedirs(WORK)
    dll = os.path.join(WORK, 'cfmsDrape.dll')
    shutil.copyfile(dll_src, dll)
    sha = hashlib.sha256(open(dll, 'rb').read()).hexdigest()
    base = {'engine': engine, 'dll': dll, 'mode': 'run', 'dt': DT0, 'amp': 0.0, 'seed': 0}
    t0 = time.perf_counter()

    wave1 = [dict(base, id='base_' + a, arm=a) for a in ARMS]                                             # Q1 · Q5
    for a in ('N', 'K'):
        for amp in (1e-10, 1e-6):
            for seed in (1, 2, 3):
                wave1.append(dict(base, id='perturb_%s_%g_s%d' % (a, amp, seed), arm=a, amp=amp, seed=seed))     # Q4
        for div in (2, 4, 8):
            wave1.append(dict(base, id='dt_%s_div%d' % (a, div), arm=a, dt=DT0 / div))                          # Q6
        for seed in (1, 2):
            wave1.append(dict(base, id='perturb_dt4_%s_1e-10_s%d' % (a, seed), arm=a, dt=DT0 / 4, amp=1e-10, seed=seed))
    print('Cusick 드레이프 계수 진단 — BE_U %g, DLL %s (sha256 %s…), 프로세스 %d 개 동시' % (BE, dll_src, sha[:12], jobs))
    with ThreadPoolExecutor(jobs) as ex:
        res1 = list(ex.map(lambda s: launch(s, python), wave1))
    by = {r['id']: r for r in res1}
    for r in res1:
        print(' ' + row(r))

    wave2 = []
    if all('error' not in by['base_' + a] for a in ('N', 'K')):
        end = {a: os.path.join(WORK, 'base_%s.npy' % a) for a in ('N', 'K')}
        for a in ('N', 'K'):
            for src in ('N', 'K'):
                wave2.append(dict(base, id='cross_%s_from_%s' % (a, src), arm=a, start=end[src]))                # Q3
            wave2.append(dict(base, id='force_' + a, arm=a, mode='force',
                              force_jobs=[[end[s], os.path.join(WORK, 'force_%s_at_%s.npy' % (a, s))] for s in ('N', 'K')]))   # Q2
    with ThreadPoolExecutor(jobs) as ex:
        res2 = list(ex.map(lambda s: launch(s, python), wave2))
    force = {}
    for r in res2:
        if r.get('force_jobs') or 'error' in r and r['id'].startswith('force_'):
            if 'error' in r:
                print(' %-22s 오류 %s' % (r['id'], r['error'][-200:].replace('\n', ' | ')))
            continue
        print(' ' + row(r))
    if wave2 and all(os.path.isfile(os.path.join(WORK, 'force_%s_at_%s.npy' % (a, s))) for a in 'NK' for s in 'NK'):
        for s in ('N', 'K'):
            fn = np.load(os.path.join(WORK, 'force_N_at_%s.npy' % s))
            fk = np.load(os.path.join(WORK, 'force_K_at_%s.npy' % s))
            force['at_%s_end' % s] = {'max_abs_force': float(np.abs(fk).max()), 'max_abs_diff_N_minus_K': float(np.abs(fn - fk).max()),
                                      'relative': float(np.abs(fn - fk).max() / np.abs(fk).max())}
            print(' 두 커널의 내부 힘, %s 의 끝 위치에서: 차이 최대 %.2e (힘 최대 %.3e 의 %.1e 배)' % (
                s, force['at_%s_end' % s]['max_abs_diff_N_minus_K'], force['at_%s_end' % s]['max_abs_force'], force['at_%s_end' % s]['relative']))
    dist = {}
    for a, b in (('base_N', 'base_K'), ('base_N', 'base_membrane_off'), ('base_K', 'base_bend_off'), ('dt_N_div4', 'dt_K_div4'),
                 ('dt_N_div8', 'dt_K_div8'), ('dt_N_div4', 'dt_N_div8'), ('dt_K_div4', 'dt_K_div8'), ('cross_K_from_N', 'base_N'),
                 ('cross_N_from_K', 'base_K'), ('cross_N_from_N', 'base_N'), ('cross_K_from_K', 'base_K')):
        fa, fb = os.path.join(WORK, a + '.npy'), os.path.join(WORK, b + '.npy')
        if os.path.isfile(fa) and os.path.isfile(fb):
            pa, pb = np.load(fa), np.load(fb)
            dist['%s|%s' % (a, b)] = {'max': float(np.abs(pa - pb).max()), 'rms': float(np.sqrt(((pa - pb) ** 2).sum(axis=1).mean()))}
            print(' 끝 위치의 거리 %-34s 최대 %.3e cm · rms %.3e cm' % (a + ' · ' + b, dist['%s|%s' % (a, b)]['max'], dist['%s|%s' % (a, b)]['rms']))

    out = {'be_u': BE, 'dll_source': dll_src, 'dll_sha256': sha, 'engine': engine, 'dt0': DT0, 'dyn_seconds': DYN_SECONDS,
           'settle_rule': {'block': BLOCK, 'max_blocks': MAX_BLOCKS, 'min_blocks': MIN_BLOCKS, 'delta_cm': SETTLE},
           'extra_static_steps': EXTRA_BLOCKS * BLOCK, 'arms': ARMS, 'runs': res1 + res2, 'force': force, 'end_position_distance_cm': dist,
           'seconds': time.perf_counter() - t0}
    with open(os.path.join(ROOT, 'results', 'diag_cusick_dc.json'), 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print('→ results/diag_cusick_dc.json (%.0f s)' % out['seconds'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
