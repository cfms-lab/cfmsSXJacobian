"""시간 판정 — 시험 exe 를 고정 로드 주소 8 곳에 다시 링크해 돌리고, 변형마다 최솟값으로 판정한다.

왜: 같은 기계어라도 실행 파일이 놓이는 로드 주소(ASLR)에 따라 함수 하나가 1.8 – 3.3 배 느려진다.
    한 exe 안에서 5 회 중 최소를 해도, 같은 exe 를 다시 돌려도 없어지지 않는다(그 exe 의 주소가 그대로라서).
    고정 주소로 링크하면 실행마다 값이 같고, 주소를 바꾸면 느려지는 함수가 바뀐다. docs/GATE_SPEED.md
    「측정 방법 변경」.

판정(기준은 각 판정문에 실행 전 고정한 그대로, 바뀐 것은 시간을 모으는 방법뿐):
  S2   (docs/GATE_SPEED.md)    min ns(EvalStaged) ≤ 2.0 × min ns(EvalFrozen)
  R3′  (docs/GATE_ASSEMBLY.md) min ns(EvalImplicit + b · f) − min ns(EvalStaged)
                               ≤ [min ns(EvalStaged + 기존 조립) − min ns(EvalStaged)] / 3
  E4   (docs/GATE_ENGINE.md)   min ns(모으기 + EvalStaged + 출력 저장) ≤ min ns(엔진 CheckBendingChange_Macri_SOA)
  B8   (docs/GATE_BRIDSON.md)  min ns(모으기 + EvalImplicit + 출력 저장) ≤ min ns(엔진 Bridson: Macri θ + GetImpBending_Oh_SOA,
                               쌍 슬롯 환원 없음)
  N2   (docs/GATE_GAUSS_NEWTON.md) min ns(모으기 + BridsonBend::EvalImplicitGN + 출력 저장) ≤ min ns(B8 의 엔진)
  (다) (docs/GATE_INTEGRATION.md) 힘마다 기준 짝의 min 시간(GN) ≤ min 시간(정확판), 시간 = 조립 + 풀이 / 모의 1 초.
                               (가) · (나) 는 test_integration.exe 가 results/gate_integration.json 에 낸 것을 읽어 권고를 낸다
  (다′) (docs/GATE_TENSION.md) 기준 짝(전단 E · Bridson G · 뻣뻣한 굽힘 · 1/60)의 min 시간(인장 T) ≤ 1.05 × min 시간(인장 E).
                               (가′) · (나′) · (라′) 는 results/gate_tension_integration.json 에서 읽어 권고를 낸다
  (다″) (docs/GATE_SHEAR_GN.md) 인장 T 에서 기준 짝 둘(Bridson G · 1/60 · 부드러운 · 뻣뻣한 굽힘)마다
                               min 시간(전단 G) ≤ 1.05 × min 시간(전단 E). (가″) · (나″) · (라″) 는
                               results/gate_shear_gn_integration.json 에서 읽어 권고를 낸다
  min =고정 로드 주소 8 곳에서 잰 값(각각 한 exe 안 5 회 중 최소)의 최솟값. 중앙값도 보고한다.

사용:  python tests/timing_bases.py <빌드 폴더> [core|engine|bridson|gn|integ|tension|shear_gn]   (리포 루트에서, vcvars64 환경 — run_gate.bat 이 부른다)
       core(기본) = S2 · R3′, 빌드 폴더에 test_speed.obj · test_assembly.obj · SXMesh.obj 가 있어야 한다.
       engine     = E4, test_engine.obj · engine_macri.obj · SDEmysse.obj · SXMesh.obj.
       bridson    = B8, test_bridson.obj(/DSX_WITH_ENGINE) · engine_macri.obj · SDEmysse.obj · SXMesh.obj.
       gn         = N2, test_gauss_newton.obj(/DSX_WITH_ENGINE) · engine_macri.obj · SDEmysse.obj · SXMesh.obj.
       integ      = (다), test_integration.obj · SXMesh.obj(엔진과 함께 빌드했으면 engine_macri.obj · SDEmysse.obj 도),
                    exe 를 timing 인자로 돌린다. 먼저 인자 없이 한 번 돌려 results/gate_integration.json 이 있어야 한다.
       tension    = (다′), integ 와 같은 obj, exe 를 timing_tension 인자로 돌린다. 먼저 tension 인자로 한 번 돌려
                    results/gate_tension_integration.json 이 있어야 한다. 숫자는 results/gate_timing_tension.json.
       shear_gn   = (다″), integ 와 같은 obj, exe 를 timing_shear_gn 인자로 돌린다. 먼저 shear_gn 인자로 한 번 돌려
                    results/gate_shear_gn_integration.json 이 있어야 한다. 숫자는 results/gate_timing_shear_gn.json.
숫자는 results/gate_timing.json(core) · results/gate_timing_engine.json(engine) · results/gate_timing_bridson.json(bridson)
· results/gate_timing_gn.json(gn) · results/gate_timing_integration.json(integ).
각 exe 가 쓰는 results/gate_*.json 의 ns_per_hinge(gn 은 ns_per_stencil, integ 는 ms_per_sim_second)는 마지막 실행 한 번의 값이다.
"""
import json
import os
import statistics
import subprocess
import sys

# 로드 주소(64 KiB 정렬). 0x140000000 = x64 exe 기본값. 실행 전에 고정한 목록이다.
BASES = ['0x140000000', '0x140010000', '0x140020000', '0x140030000',
         '0x150000000', '0x160000000', '0x1A0000000', '0x240000000']
TESTS = [('test_speed', 'results/gate_speed.json'), ('test_assembly', 'results/gate_assembly.json')]
ENGINE_TEST = ('test_engine', 'results/gate_engine.json', ['engine_macri.obj', 'SDEmysse.obj'])
BRIDSON_TEST = ('test_bridson', 'results/gate_bridson.json', ['engine_macri.obj', 'SDEmysse.obj'])
GN_TEST = ('test_gauss_newton', 'results/gate_gauss_newton.json', ['engine_macri.obj', 'SDEmysse.obj'])
INTEG_TEST = ('test_integration', 'results/gate_integration_time.json', ['engine_macri.obj', 'SDEmysse.obj'])
TENSION_TEST = ('test_integration', 'results/gate_tension_time.json', ['engine_macri.obj', 'SDEmysse.obj'])
SHEAR_GN_TEST = ('test_integration', 'results/gate_shear_gn_time.json', ['engine_macri.obj', 'SDEmysse.obj'])


def run_bases(out, stem, jpath, extra=(), key='ns_per_hinge', args=()):
    per = []
    objs = [os.path.join(out, o) for o in [stem + '.obj'] + list(extra) + ['SXMesh.obj']]
    for i, base in enumerate(BASES):
        exe = os.path.join(out, '%s_base%d.exe' % (stem, i))
        subprocess.run(['link', '/nologo'] + objs + ['/DYNAMICBASE:NO', '/FIXED', '/BASE:' + base, '/OUT:' + exe],
                       check=True, stdout=subprocess.DEVNULL)
        r = subprocess.run([exe] + list(args), stdout=subprocess.DEVNULL)
        if r.returncode != 0 and i == 0:   # 정확도 판정은 exe 의 몫이다. 시간은 그래도 잰다
            print('  주의: %s 의 정확도 판정이 실패했다(종료 코드 %d) — 시간만 잰다' % (stem, r.returncode))
        ns = json.load(open(jpath, encoding='utf-8'))[key]
        per.append(ns)
        print('  %-13s base %s  %s' % (stem, base, '  '.join('%s %.0f' % (k, v) for k, v in ns.items())))
    keys = list(per[0])
    return {'per_base': per,
            'min': {k: min(p[k] for p in per) for k in keys},
            'median': {k: statistics.median(p[k] for p in per) for k in keys}}


def engine(out):
    stem, jpath, extra = ENGINE_TEST
    print('엔진 손 코드 비교 E4 시간 — 고정 로드 주소 %d 곳, 변형마다 최솟값' % len(BASES))
    agg = run_bases(out, stem, jpath, extra)
    res = {'bases': BASES, 'raw': agg, 'gates': {}}
    for stat in ('min', 'median'):
        t = agg[stat]
        res[stat] = {'E4_gen_store_over_engine': round(t['gen_staged_store'] / t['engine_macri'], 4),
                     'gen_staged_over_engine': round(t['gen_staged'] / t['engine_macri'], 4),
                     'gen_implicit_over_engine': round(t['gen_implicit'] / t['engine_macri'], 4),
                     'compact_store_over_engine': round(t['compact_staged_store'] / t['engine_macri'], 4),
                     'small_gen_store_over_engine': round(t['small_gen_staged_store'] / t['small_engine_macri'], 4)}
    m = agg['min']
    e4 = m['gen_staged_store'] / m['engine_macri']
    ok = e4 <= 1.0
    res['gates']['E4'] = {'value': round(e4, 4), 'limit': 1.0, 'pass': ok}
    print('판정 (최솟값)')
    print('  %-4s %s  생성(EvalStaged + 저장) / 엔진 Macri = %.3f  (기준 ≤ 1.000)' % ('E4', 'PASS' if ok else 'FAIL', e4))
    print('  보고 — 힌지당 ns: 엔진 %.0f, 생성 + 저장 %.0f, 저장 없음 %.0f, EvalFrozen %.0f, 압축 + 저장 %.0f, EvalImplicit %.0f'
          % (m['engine_macri'], m['gen_staged_store'], m['gen_staged'], m['gen_frozen'], m['compact_staged_store'], m['gen_implicit']))
    print('  보고 — 작은 격자: 엔진 %.0f, 생성 + 저장 %.0f ns → 비 %.3f'
          % (m['small_engine_macri'], m['small_gen_staged_store'], m['small_gen_staged_store'] / m['small_engine_macri']))
    md = res['median']
    print('보고 (중앙값): E4 %.3f' % md['E4_gen_store_over_engine'])
    pb = [round(p['gen_staged_store'] / p['engine_macri'], 3) for p in agg['per_base']]
    res['per_base_ratios'] = {'E4': pb}
    print('보고 (로드 주소마다 따로 낸 비): E4 %s' % pb)
    res['gates_failed'] = 0 if ok else 1
    json.dump(res, open('results/gate_timing_engine.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print('\n시간 판정 E4: %s → results/gate_timing_engine.json' % ('GO' if ok else 'NO-GO'))
    return 0 if ok else 1


def bridson(out):
    stem, jpath, extra = BRIDSON_TEST
    print('Bridson 굽힘 B8 시간 — 고정 로드 주소 %d 곳, 변형마다 최솟값' % len(BASES))
    agg = run_bases(out, stem, jpath, extra)
    res = {'bases': BASES, 'raw': agg, 'gates': {}}
    for stat in ('min', 'median'):
        t = agg[stat]
        e = t['engine_oh']
        res[stat] = {'B8_gen_implicit_store_over_engine': round(t['gen_implicit_store'] / e, 4),
                     'gen_staged_store_over_engine': round(t['gen_staged_store'] / e, 4),
                     'gen_implicit_scatter_over_engine': round(t['gen_implicit_scatter'] / e, 4),
                     'engine_serial_over_engine': round(t['engine_oh_serial'] / e, 4),
                     'gen_implicit_store_over_engine_serial': round(t['gen_implicit_store'] / t['engine_oh_serial'], 4),
                     'small_gen_over_engine': round(t['small_gen_implicit_store'] / t['small_engine_oh'], 4)}
    m = agg['min']
    b8 = m['gen_implicit_store'] / m['engine_oh']
    ok = b8 <= 1.0
    res['gates']['B8'] = {'value': round(b8, 4), 'limit': 1.0, 'pass': ok}
    print('판정 (최솟값)')
    print('  %-4s %s  생성(EvalImplicit + 저장) / 엔진 Bridson(환원 없음) = %.3f  (기준 ≤ 1.000)' % ('B8', 'PASS' if ok else 'FAIL', b8))
    print('  보고 — 힌지당 ns: 엔진 %.0f (기본값 환원 포함 %.0f), 생성 EvalImplicit + 저장 %.0f, EvalStaged + 저장 %.0f, EvalImplicit + 흩뿌리기 %.0f'
          % (m['engine_oh'], m['engine_oh_serial'], m['gen_implicit_store'], m['gen_staged_store'], m['gen_implicit_scatter']))
    print('  보고 — 엔진 기본값 대비 %.3f, EvalStaged + 저장 / 엔진 %.3f; 작은 격자: 엔진 %.0f, 생성 %.0f ns → 비 %.3f'
          % (m['gen_implicit_store'] / m['engine_oh_serial'], m['gen_staged_store'] / m['engine_oh'],
             m['small_engine_oh'], m['small_gen_implicit_store'], m['small_gen_implicit_store'] / m['small_engine_oh']))
    md = res['median']
    print('보고 (중앙값): B8 %.3f' % md['B8_gen_implicit_store_over_engine'])
    pb = [round(p['gen_implicit_store'] / p['engine_oh'], 3) for p in agg['per_base']]
    res['per_base_ratios'] = {'B8': pb}
    print('보고 (로드 주소마다 따로 낸 비): B8 %s' % pb)
    res['gates_failed'] = 0 if ok else 1
    json.dump(res, open('results/gate_timing_bridson.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print('\n시간 판정 B8: %s → results/gate_timing_bridson.json' % ('GO' if ok else 'NO-GO'))
    return 0 if ok else 1


def gn(out):
    stem, jpath, extra = GN_TEST
    print('Gauss–Newton N2 시간 — 고정 로드 주소 %d 곳, 변형마다 최솟값' % len(BASES))
    agg = run_bases(out, stem, jpath, extra, key='ns_per_stencil')
    res = {'bases': BASES, 'raw': agg, 'gates': {}}
    specs = [('bridson', 'Bridson'), ('baraff', 'Baraff 굽힘'), ('compact', '압축 명세'), ('stretch', '인장'), ('shear', '전단')]
    for stat in ('min', 'median'):
        t = agg[stat]
        e = t['engine_oh']
        r = {'N2_bridson_gn_store_over_engine': round(t['bridson_gn_store'] / e, 4),
             'bridson_gn_store_over_engine_serial': round(t['bridson_gn_store'] / t['engine_oh_serial'], 4),
             'bridson_gn_scatter_over_engine': round(t['bridson_gn_scatter'] / e, 4),
             'bridson_implicit_store_over_engine': round(t['bridson_implicit_store'] / e, 4),
             'small_bridson_gn_over_engine': round(t['small_bridson_gn_store'] / t['small_engine_oh'], 4)}
        for k, _ in specs:
            r['%s_gn_over_implicit' % k] = round(t['%s_gn_store' % k] / t['%s_implicit_store' % k], 4)
        res[stat] = r
    m = agg['min']
    n2 = m['bridson_gn_store'] / m['engine_oh']
    ok = n2 <= 1.0
    res['gates']['N2'] = {'value': round(n2, 4), 'limit': 1.0, 'pass': ok}
    print('판정 (최솟값)')
    print('  %-4s %s  생성 GN(EvalImplicitGN + 저장) / 엔진 Bridson(환원 없음) = %.3f  (기준 ≤ 1.000)' % ('N2', 'PASS' if ok else 'FAIL', n2))
    print('  보고 — 힌지당 ns: 엔진 %.0f (기본값 환원 포함 %.0f), 생성 GN + 저장 %.0f, 정확판 + 저장 %.0f, GN + 흩뿌리기 %.0f'
          % (m['engine_oh'], m['engine_oh_serial'], m['bridson_gn_store'], m['bridson_implicit_store'], m['bridson_gn_scatter']))
    print('  보고 — 엔진 기본값 대비 %.3f; 작은 격자: 엔진 %.0f, 생성 GN %.0f ns → 비 %.3f'
          % (res['min']['bridson_gn_store_over_engine_serial'], m['small_engine_oh'], m['small_bridson_gn_store'],
             res['min']['small_bridson_gn_over_engine']))
    print('  보고 — GN / 정확판(모으기 + 저장): ' + ', '.join('%s %.0f / %.0f = %.3f' % (lab, m['%s_gn_store' % k], m['%s_implicit_store' % k],
                                                                                      res['min']['%s_gn_over_implicit' % k]) for k, lab in specs))
    print('보고 (중앙값): N2 %.3f' % res['median']['N2_bridson_gn_store_over_engine'])
    pb = [round(p['bridson_gn_store'] / p['engine_oh'], 3) for p in agg['per_base']]
    res['per_base_ratios'] = {'N2': pb}
    print('보고 (로드 주소마다 따로 낸 비): N2 %s' % pb)
    res['gates_failed'] = 0 if ok else 1
    json.dump(res, open('results/gate_timing_gn.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print('\n시간 판정 N2: %s → results/gate_timing_gn.json' % ('GO' if ok else 'NO-GO'))
    return 0 if ok else 1


def integ(out):
    stem, jpath, extra = INTEG_TEST
    full = json.load(open('results/gate_integration.json', encoding='utf-8'))
    if not full['engine']:
        extra = []
    print('시간 적분 비교 (다) 시간 — 고정 로드 주소 %d 곳, 기준 짝마다 최솟값 (ms / 모의 1 초)' % len(BASES))
    agg = run_bases(out, stem, jpath, extra, key='ms_per_sim_second', args=['timing'])
    tj = json.load(open(jpath, encoding='utf-8'))   # 마지막 실행 — 안정 · 반복 수는 실행마다 같다
    stable = tj['stable']
    asm_last = tj['asm_ms_per_sim_second']
    iters = tj['pcg_iters']
    m = agg['min']
    forces = [('baraff', 'Baraff 굽힘', 'baraff', 'EE_baraffE', 'EE_baraffG'),
              ('compact', '압축 명세', 'baraff', 'EE_compactE', 'EE_compactG'),
              ('stretch', '인장', 'stretch', 'EE_bridsonG', 'GE_bridsonG'),
              ('shear', '전단', 'shear', 'EE_bridsonG', 'EG_bridsonG'),
              ('bridson', 'Bridson 굽힘 (확인)', 'bridson', 'EE_bridsonE', 'EE_bridsonG')]
    res = {'bases': BASES, 'raw': agg, 'stable': stable, 'pcg_iters': iters, 'asm_ms_per_sim_second_last': asm_last,
           'I0_pass': full['I0']['pass'], 'forces': {}}
    print('권고 (가) 안정 · (나) 차이 = results/gate_integration.json, (다) 시간 = 위 최솟값')
    for key, label, src, kE, kG in forces:
        rule = full['rules'][src]
        tE, tG = m[kE], m[kG]
        sE, sG = bool(stable[kE]), bool(stable[kG])
        c_ok = (tG <= tE) if sE else sG
        c_na = not sE and not sG   # 둘 다 불안정 — 판정문이 정하지 않은 경우, 넘지 못한 것으로 센다
        a_ok, b_ok = rule['A_ok'], rule['B_ok']
        rec = 'GN' if (a_ok and b_ok and c_ok) else 'exact'
        res['forces'][key] = {'A_ok': a_ok, 'B_ok': b_ok, 'C_ok': c_ok, 'C_not_measurable': c_na, 'time_E': round(tE, 2),
                              'time_G': round(tG, 2), 'G_over_E': round(tG / tE, 4), 'iters_E': iters[kE], 'iters_G': iters[kG],
                              'stable_E': sE, 'stable_G': sG, 'recommend': rec}
        if key == 'bridson':
            verdict = 'GN 기본 유지' if rec == 'GN' else ('확인 못 함' if c_na else '결정을 다시 볼 근거')
        else:
            verdict = 'GN 을 기본으로 권함' if rec == 'GN' else '정확판 유지'
        c_txt = '판정 불가(기준 짝 둘 다 불안정)' if c_na else ('넘음' if c_ok else '못 넘음')
        print('  %-18s (가) %s · (나) %s · (다) %s  시간 G / E = %.1f / %.1f = %.3f, PCG 반복 %d / %d  →  %s'
              % (label, '넘음' if a_ok else '못 넘음', '넘음' if b_ok else '못 넘음', c_txt,
                 tG, tE, tG / tE, iters[kG], iters[kE], verdict))
    if 'EE_engine' in m:
        print('  보고 — 엔진 Bridson %.1f ms/s (생성 Bridson G %.1f, E %.1f), PCG 반복 %d'
              % (m['EE_engine'], m['EE_bridsonG'], m['EE_bridsonE'], iters['EE_engine']))
    # 보고(판정 밖, 첫 실행 뒤 더함): 뻣뻣한 굽힘의 같은 짝 — 기준 짝이 양쪽 다 불안정이라 시간을 견줄 곳이 없어서
    res['stiff_report'] = {}
    print('보고 — 뻣뻣한 굽힘(k = 1e-3, h = 1/60)의 같은 짝 (판정 밖)')
    for key, label, src, kE, kG in forces:
        sE_, sG_ = kE + '_stiff', kG + '_stiff'
        if sE_ not in m:
            continue
        tE, tG = m[sE_], m[sG_]
        a_tE, a_tG = asm_last[sE_], asm_last[sG_]
        res['stiff_report'][key] = {'time_E': round(tE, 2), 'time_G': round(tG, 2), 'G_over_E': round(tG / tE, 4),
                                    'stable_E': bool(stable[sE_]), 'stable_G': bool(stable[sG_]),
                                    'iters_E': iters[sE_], 'iters_G': iters[sG_],
                                    'asm_last_E': round(a_tE, 2), 'asm_last_G': round(a_tG, 2)}
        print('  %-18s 안정 E %s · G %s  시간 G / E = %.1f / %.1f = %.3f, PCG 반복 %d / %d, 조립(마지막 실행) %.1f / %.1f'
              % (label, '예' if stable[sE_] else '아니오', '예' if stable[sG_] else '아니오', tG, tE, tG / tE,
                 iters[sG_], iters[sE_], a_tG, a_tE))
    if 'EE_engine_stiff' in m:
        res['stiff_report']['engine'] = {'time': round(m['EE_engine_stiff'], 2), 'iters': iters['EE_engine_stiff'],
                                         'over_bridsonG': round(m['EE_engine_stiff'] / m['EE_bridsonG_stiff'], 4)}
        print('  엔진 Bridson %.1f ms/s (생성 Bridson G 대비 %.3f), PCG 반복 %d'
              % (m['EE_engine_stiff'], m['EE_engine_stiff'] / m['EE_bridsonG_stiff'], iters['EE_engine_stiff']))
    pb = {k: [round(p[kG] / p[kE], 3) for p in agg['per_base']] for k, _, _, kE, kG in forces}
    pb.update({k + '_stiff': [round(p[kG + '_stiff'] / p[kE + '_stiff'], 3) for p in agg['per_base']]
               for k, _, _, kE, kG in forces if kE + '_stiff' in m})
    res['per_base_ratios'] = pb
    print('보고 (로드 주소마다 따로 낸 G / E): ' + ', '.join('%s %s' % (k, v) for k, v in pb.items()))
    json.dump(res, open('results/gate_timing_integration.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    ok = full['I0']['pass']
    print('\n시간 적분 비교 (다) 와 권고 → results/gate_timing_integration.json%s' % ('' if ok else ' (I0 실패 — 권고 무효)'))
    return 0 if ok else 1


def tension(out):
    """인장 C > 0 판 게이트 (다′) 와 권고(docs/GATE_TENSION.md). (가′) · (나′) · (라′) 는 test_integration.exe tension 의 것."""
    stem, jpath, extra = TENSION_TEST
    full = json.load(open('results/gate_tension_integration.json', encoding='utf-8'))
    if not full['engine']:
        extra = []
    print('인장 C > 0 판 (다′) 시간 — 고정 로드 주소 %d 곳, 최솟값 (ms / 모의 1 초)' % len(BASES))
    agg = run_bases(out, stem, jpath, extra, key='ms_per_sim_second', args=['timing_tension'])
    tj = json.load(open(jpath, encoding='utf-8'))   # 쓸 수 있음 · 반복 수는 실행마다 같다
    usable, iters, asm_last = tj['usable'], tj['pcg_iters'], tj['asm_ms_per_sim_second']
    m = agg['min']
    kE, kT = 'EE_bridsonG_stiff', 'TE_bridsonG_stiff'
    uE, uT = bool(usable[kE]), bool(usable[kT])
    tE, tT = m[kE], m[kT]
    c_ok = (tT <= 1.05 * tE) if (uE and uT) else (uT and not uE)
    rules = full['rules']
    a_ok, b_ok, d_ok = rules['A']['ok'], rules['B']['ok'], rules['D']['ok']
    rec = 'tension' if (a_ok and b_ok and c_ok and d_ok) else 'exact'
    res = {'bases': BASES, 'raw': agg, 'usable': usable, 'pcg_iters': iters, 'asm_ms_per_sim_second_last': asm_last,
           'C': {'time_E': round(tE, 2), 'time_T': round(tT, 2), 'T_over_E': round(tT / tE, 4), 'limit': 1.05,
                 'usable_E': uE, 'usable_T': uT, 'ok': c_ok},
           'A_ok': a_ok, 'B_ok': b_ok, 'D_ok': d_ok, 'recommend': rec}
    print('판정 (최솟값)')
    print('  (다′) %s  T / E = %.1f / %.1f = %.3f  (기준 ≤ 1.050), 쓸 수 있음 E %s · T %s, PCG 반복 %d / %d'
          % ('넘음' if c_ok else '못 넘음', tT, tE, tT / tE, '예' if uE else '아니오', '예' if uT else '아니오', iters[kT], iters[kE]))
    print('  권고 (가′) %s · (나′) %s · (다′) %s · (라′) %s  →  %s'
          % ('넘음' if a_ok else '못 넘음', '넘음' if b_ok else '못 넘음', '넘음' if c_ok else '못 넘음', '넘음' if d_ok else '못 넘음',
             '인장 기본을 EvalImplicitTension 으로 권함' if rec == 'tension' else '정확판 유지'))
    sE, sT = 'EE_bridsonG', 'TE_bridsonG'   # 보고: 부드러운 굽힘의 같은 짝(둘 다 쓸 수 있을 때만 뜻이 있다)
    res['soft_report'] = {'time_E': round(m[sE], 2), 'time_T': round(m[sT], 2), 'usable_E': bool(usable[sE]),
                          'usable_T': bool(usable[sT]), 'iters_E': iters[sE], 'iters_T': iters[sT]}
    print('  보고 — 부드러운 굽힘 같은 짝: 쓸 수 있음 E %s · T %s, 시간 T %.1f · E %.1f, PCG 반복 %d / %d'
          % ('예' if usable[sE] else '아니오', '예' if usable[sT] else '아니오', m[sT], m[sE], iters[sT], iters[sE]))
    pb = [round(p[kT] / p[kE], 3) for p in agg['per_base']]
    res['per_base_ratios'] = {'C': pb}
    print('보고 (로드 주소마다 따로 낸 T / E): %s' % pb)
    json.dump(res, open('results/gate_timing_tension.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print('\n인장 C > 0 판 (다′) 와 권고 → results/gate_timing_tension.json')
    return 0


def shear_gn(out):
    """전단 GN 재비교 (다″) 와 권고(docs/GATE_SHEAR_GN.md). (가″) · (나″) · (라″) 는 test_integration.exe shear_gn 의 것."""
    stem, jpath, extra = SHEAR_GN_TEST
    full = json.load(open('results/gate_shear_gn_integration.json', encoding='utf-8'))
    if not full['engine']:
        extra = []
    print('전단 GN 재비교 (다″) 시간 — 인장 T, 고정 로드 주소 %d 곳, 최솟값 (ms / 모의 1 초)' % len(BASES))
    agg = run_bases(out, stem, jpath, extra, key='ms_per_sim_second', args=['timing_shear_gn'])
    tj = json.load(open(jpath, encoding='utf-8'))   # 쓸 수 있음 · 반복 수는 실행마다 같다
    usable, iters, asm_last = tj['usable'], tj['pcg_iters'], tj['asm_ms_per_sim_second']
    m = agg['min']
    refs = {'soft': ('TE_bridsonG', 'TG_bridsonG'), 'stiff': ('TE_bridsonG_stiff', 'TG_bridsonG_stiff')}
    res = {'bases': BASES, 'raw': agg, 'usable': usable, 'pcg_iters': iters, 'asm_ms_per_sim_second_last': asm_last, 'C': {}}
    c_all = True
    print('판정 (최솟값)')
    for name, (kE, kG) in refs.items():
        uE, uG = bool(usable[kE]), bool(usable[kG])
        tE, tG = m[kE], m[kG]
        ok = (tG <= 1.05 * tE) if (uE and uG) else (uG and not uE)
        c_all = c_all and ok
        res['C'][name] = {'time_E': round(tE, 2), 'time_G': round(tG, 2), 'G_over_E': round(tG / tE, 4), 'limit': 1.05,
                          'usable_E': uE, 'usable_G': uG, 'iters_E': iters[kE], 'iters_G': iters[kG], 'ok': ok}
        print('  (다″) %-5s %s  G / E = %.1f / %.1f = %.3f  (기준 ≤ 1.050), 쓸 수 있음 E %s · G %s, PCG 반복 %d / %d'
              % (name, '넘음' if ok else '못 넘음', tG, tE, tG / tE, '예' if uE else '아니오', '예' if uG else '아니오', iters[kG], iters[kE]))
    rules = full['rules']
    a_ok, b_ok, d_ok = rules['A']['ok'], rules['B']['ok'], rules['D']['ok']
    rec = 'gn' if (a_ok and b_ok and c_all and d_ok) else 'exact'
    res.update({'C_ok': c_all, 'A_ok': a_ok, 'B_ok': b_ok, 'D_ok': d_ok, 'recommend': rec})
    print('  권고 (가″) %s · (나″) %s · (다″) %s · (라″) %s  →  %s'
          % ('넘음' if a_ok else '못 넘음', '넘음' if b_ok else '못 넘음', '넘음' if c_all else '못 넘음', '넘음' if d_ok else '못 넘음',
             '전단 기본을 EvalImplicitGN 으로 권함' if rec == 'gn' else '정확판 유지'))
    pb = {name: [round(p[kG] / p[kE], 3) for p in agg['per_base']] for name, (kE, kG) in refs.items()}
    res['per_base_ratios'] = pb
    print('보고 (로드 주소마다 따로 낸 G / E): 부드러움 %s' % pb['soft'])
    print('                                  뻣뻣함   %s' % pb['stiff'])
    json.dump(res, open('results/gate_timing_shear_gn.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print('\n전단 GN 재비교 (다″) 와 권고 → results/gate_timing_shear_gn.json')
    return 0


def main():
    if len(sys.argv) not in (2, 3) or (len(sys.argv) == 3 and sys.argv[2] not in ('core', 'engine', 'bridson', 'gn', 'integ', 'tension', 'shear_gn')):
        print(__doc__)
        return 2
    out = sys.argv[1]
    if len(sys.argv) == 3 and sys.argv[2] == 'shear_gn':
        return shear_gn(out)
    if len(sys.argv) == 3 and sys.argv[2] == 'tension':
        return tension(out)
    if len(sys.argv) == 3 and sys.argv[2] == 'integ':
        return integ(out)
    if len(sys.argv) == 3 and sys.argv[2] == 'gn':
        return gn(out)
    if len(sys.argv) == 3 and sys.argv[2] == 'engine':
        return engine(out)
    if len(sys.argv) == 3 and sys.argv[2] == 'bridson':
        return bridson(out)
    print('시간 판정 — 고정 로드 주소 %d 곳, 변형마다 최솟값' % len(BASES))
    agg = {stem: run_bases(out, stem, jpath) for stem, jpath in TESTS}

    res = {'bases': BASES, 'raw': agg, 'gates': {}}
    fail = 0

    def gate(name, value, limit, label):
        nonlocal fail
        ok = value <= limit
        fail += 0 if ok else 1
        print('  %-4s %s  %s = %.3f  (기준 ≤ %.3f)' % (name, 'PASS' if ok else 'FAIL', label, value, limit))
        res['gates'][name] = {'value': round(value, 4), 'limit': round(limit, 4), 'pass': ok}

    for stat in ('min', 'median'):
        sp = agg['test_speed'][stat]
        a = agg['test_assembly'][stat]
        s2 = sp['eval_staged'] / sp['eval_frozen']
        inc_old = a['staged_old_assembly'] - a['staged_eval']
        inc_imp = a['eval_implicit'] - a['staged_eval']
        inc_new = a['staged_new_assembly'] - a['staged_eval']
        res[stat] = {'S2_staged_over_frozen': round(s2, 4),
                     'R3p_implicit_over_old': round(inc_imp / inc_old, 4),
                     'R3_first_attempt_over_old': round(inc_new / inc_old, 4),
                     'assembly_old_ns': round(inc_old, 1), 'implicit_increment_ns': round(inc_imp, 1),
                     'first_attempt_increment_ns': round(inc_new, 1)}
        if stat == 'min':
            print('판정 (최솟값)')
            gate('S2', s2, 2.0, 'EvalStaged / EvalFrozen')
            gate('R3′', inc_imp / inc_old, 1. / 3., 'EvalImplicit 조립 몫 / 기존 조립')
            print('  보고 — 1 차 AssembleImplicit / 기존 조립 = %.3f' % (inc_new / inc_old))
            print('  보고 — 힌지당 ns: 평가 %.0f, + 기존 조립 %.0f (조립 %.0f), EvalImplicit %.0f (조립 몫 %.0f), 1 차 %.0f (조립 %.0f)'
                  % (a['staged_eval'], a['staged_old_assembly'], inc_old, a['eval_implicit'], inc_imp,
                     a['staged_new_assembly'], inc_new))
        else:
            print('보고 (중앙값): S2 %.3f, R3′ %.3f, 1 차 %.3f' % (s2, inc_imp / inc_old, inc_new / inc_old))
    # 보고: 로드 주소마다 따로 낸 비(같은 실행 안의 차) — 함수마다 느려짐이 달라 흩어진다
    pb = {'S2': [], 'R3p': [], 'R3_first': []}
    for sp, a in zip(agg['test_speed']['per_base'], agg['test_assembly']['per_base']):
        io = a['staged_old_assembly'] - a['staged_eval']
        pb['S2'].append(round(sp['eval_staged'] / sp['eval_frozen'], 3))
        pb['R3p'].append(round((a['eval_implicit'] - a['staged_eval']) / io, 3))
        pb['R3_first'].append(round((a['staged_new_assembly'] - a['staged_eval']) / io, 3))
    res['per_base_ratios'] = pb
    print('보고 (로드 주소마다 따로 낸 비): S2 %s' % pb['S2'])
    print('                                R3′ %s' % pb['R3p'])
    print('                                1 차 %s' % pb['R3_first'])
    res['gates_failed'] = fail
    json.dump(res, open('results/gate_timing.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print('\n시간 판정 S2 · R3′: %s (실패 %d) → results/gate_timing.json' % ('GO' if fail == 0 else 'NO-GO', fail))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
