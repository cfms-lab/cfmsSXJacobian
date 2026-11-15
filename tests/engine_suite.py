"""엔진 합성 시험 (docs/GATE_ENGINE_SUITE.md): 엔진 리포의 합성 시험을 세 팔로 돌려 새 기본값이 통과를 잃는지 본다.

사용:  python tests/engine_suite.py <새 엔진 DLL> <종전 엔진 DLL(378d20c)> [--run 2] [--resume] [--only a.py,b.py] [--out <json>]
--run:    판정의 차수(기본 1). 차수마다 새 엔진 커밋 · 시험 트리 폴더 · 결과 파일 · 측정 방법의 차이가 아래 RUNS 에 고정돼 있다
          (1 = 엔진 92d8838, 2 = 엔진 브랜치 961765d — docs/GATE_ENGINE_SUITE.md 「2 차」).
--resume: 시험 하나가 끝날 때마다 %TEMP%\\cfmsSXJ_engine_suite\\progress.json 에 기록을 쓴다. 기계가 내려가 끊겼으면 같은 명령을
          --resume 으로 다시 부르면 끝난 시험은 건너뛰고 이어 돌린다(같은 DLL SHA-256 · 같은 커밋일 때만; 아니면 처음부터).
          같은 시험을 돌리다 기계가 두 번 내려가면 그 시험은 건너뛰고 「판정 불가(기계)」로 적는다.
팔: N = 새 DLL · 환경 변수 없음, K = 새 DLL + 킬스위치(CFMS_CPU_SXJ_BEND=off · CFMS_CPU_SXJ_MEMBRANE=off), B = 종전 DLL.
시험 트리는 엔진 리포 밖에 둔다 — `git archive 92d8838 drape` 를 %TEMP%\\cfmsSXJ_engine_suite\\<팔> 에 풀고, 래퍼의 기본 후보 자리
drape/build/0.2.2/Release/cfmsDrape.dll 에 그 팔의 DLL 을 놓는다. 엔진 리포(main · 브랜치 · 기준 워크트리)는 읽기만 한다.
python 은 엔진의 가상 환경(cfmsDrape_dev\\drape\\.venv)을 쓴다. 엔진 main 은 CFMS_DRAPE_MAIN_DIR(없으면 ..\\cfmsDrape_dev).
--only 를 주면 그 시험만 돌리고, --out 이 없으면 결과 json 을 쓰지 않는다(하네스 확인용).
리포 루트에서 실행. 결과는 results/gate_engine_suite.json(2 차는 gate_engine_suite_2.json). 통과 0, 실패 1, 준비 안 됨 2.
"""
import sys
sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import io  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import tarfile  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENGINE_MAIN = os.path.abspath(os.environ.get('CFMS_DRAPE_MAIN_DIR') or os.path.join(ROOT, '..', 'cfmsDrape_dev'))
LINK_DIR = os.path.abspath(os.environ.get('CFMS_DRAPE_LINK_DIR') or os.path.join(ROOT, '..', 'cfmsDrape_dev_sxjacobian-link'))
BASE_DIR = os.path.abspath(os.path.join(ROOT, '..', 'cfmsDrape_dev_base378d20c'))
NEW, BASE = '92d8838', '378d20c'
PY = os.path.join(ENGINE_MAIN, 'drape', '.venv', 'Scripts', 'python.exe')
SUITE = os.path.join(os.environ.get('TEMP') or tempfile.gettempdir(), 'cfmsSXJ_engine_suite')   # 문서대로 %TEMP%
ARMS = ('N', 'K', 'B')
ARM_LABEL = {'N': '새 기본값', 'K': '킬스위치', 'B': '종전 DLL'}
TIMEOUT = 900.0                                   # 엔진 10 차 검증 기록과 같다
ARGS = {'py_velocity_readback_test.py': ['--backend', 'cpu']}   # 10 차 검증 기록의 부르는 꼴
KILL = {'CFMS_CPU_SXJ_BEND': 'off', 'CFMS_CPU_SXJ_MEMBRANE': 'off'}
DROP = ('OMP_NUM_THREADS', 'CFMS_OMP_THREADS', 'CFMS_DRAPE_DLL')
SAME_PATHS = ('drape/scripts', 'drape/cfmsdrape', 'drape/pyproject.toml')
OUT_JSON = os.path.join(ROOT, 'results', 'gate_engine_suite.json')
PROGRESS = os.path.join(SUITE, 'progress.json')   # 시험마다 쓰는 진행 기록 — --resume 이 읽는다
MAX_INTERRUPTS = 2                                # 같은 시험을 돌리다 기계가 이만큼 내려가면 그 시험은 건너뛴다(판정 불가 — 기계)
# 판정의 차수마다 고정한 것(docs/GATE_ENGINE_SUITE.md). main() 이 --run 으로 골라 위 전역(NEW · SUITE · PROGRESS · OUT_JSON)을 바꾼다.
#   new           새 엔진 커밋(N · K 팔의 DLL, 시험 트리의 git archive)
#   suffix        시험 트리 폴더 · 결과 파일 이름에 붙는 꼬리(앞 차수의 로그를 지우지 않게)
#   changed       SAME_PATHS 안에서 378d20c 와 달라도 되는 파일 — S0 (b) 는 「다른 파일이 정확히 이것뿐」
#   known_timeout 앞 차수에서 세 팔 모두 하나씩 돌려도 제한 시간을 넘은 시험 — 본 실행에서 세 팔 모두 제한 시간이면 다시 돌리지 않는다
#   extra_dll     래퍼의 기본 자리 말고도 팔의 DLL 을 놓는 자리(시험 스크립트가 박아 둔 곳)
RUNS = {
    '1': {'new': '92d8838', 'suffix': '', 'changed': (), 'known_timeout': (), 'extra_dll': ()},
    '2': {'new': '961765d', 'suffix': '_2', 'changed': ('drape/scripts/py_membrane_nonlinear_test.py',),
          'known_timeout': ('py_central_load_spin_test.py', 'py_open_shell_contact_test.py'),
          'extra_dll': ('build/Release',)},
}
CHANGED, KNOWN_TIMEOUT, EXTRA_DLL = (), (), ()


def git(*args, cwd):
    return subprocess.run(['git', '-C', cwd] + list(args), capture_output=True, text=True, encoding='utf-8', errors='replace')


def arm_env(arm):
    env = {k: v for k, v in os.environ.items() if k not in DROP and not k.startswith('CFMS_CPU_SXJ_')}
    env.update(PYTHONUTF8='1', PYTHONDONTWRITEBYTECODE='1')
    if arm == 'K':
        env.update(KILL)
    return env


def arm_dir(arm):
    return os.path.join(SUITE, arm, 'drape')


def prepare(dll_new, dll_base):
    if os.path.isdir(SUITE):
        shutil.rmtree(SUITE)
    r = subprocess.run(['git', '-C', ENGINE_MAIN, 'archive', '--format=tar', NEW, 'drape'], capture_output=True)
    if r.returncode != 0:
        raise RuntimeError('git archive 실패: %s' % r.stderr.decode('utf-8', 'replace')[-500:])
    for arm in ARMS:
        os.makedirs(os.path.join(SUITE, arm))
        with tarfile.open(fileobj=io.BytesIO(r.stdout)) as tf:
            tf.extractall(os.path.join(SUITE, arm), filter='data')
        for sub in ('build/0.2.2/Release',) + tuple(EXTRA_DLL):
            rel = os.path.join(arm_dir(arm), *sub.split('/'))
            os.makedirs(rel)
            shutil.copyfile(dll_base if arm == 'B' else dll_new, os.path.join(rel, 'cfmsDrape.dll'))
        os.makedirs(os.path.join(SUITE, 'logs', arm))


def probe(arm):
    """사본 안에서 래퍼가 기본으로 여는 DLL 과 래퍼의 위치."""
    code = ("import json, cfmsdrape; d = cfmsdrape.Drape(); "
            "print(json.dumps({'wrapper': cfmsdrape.__file__, 'dll': str(cfmsdrape.find_dll()), 'build_info': d.build_info})); d.close()")
    r = subprocess.run([PY, '-c', code], cwd=arm_dir(arm), env=arm_env(arm), capture_output=True, text=True,
                       encoding='utf-8', errors='replace', timeout=120)
    if r.returncode != 0:
        return {'error': r.stderr[-600:]}
    return json.loads(r.stdout.strip().splitlines()[-1])


def command(name):
    if name.startswith('test_'):
        return [PY, '-m', 'pytest', 'scripts/' + name, '-q']
    return [PY, 'scripts/' + name] + ARGS.get(name, [])


def kill_tree(pid):
    subprocess.run(['taskkill', '/F', '/T', '/PID', str(pid)], capture_output=True)


def finish(arm, name, proc, t0, fo, fe, timed_out):
    fo.close()
    fe.close()
    base = os.path.join(SUITE, 'logs', arm, name)
    out = open(base + '.out', 'rb').read()
    err = open(base + '.err', 'rb').read()
    rec = {'exit': None if timed_out else proc.returncode, 'timeout': timed_out, 'seconds': round(time.perf_counter() - t0, 1),
           'stdout_sha256': hashlib.sha256(out).hexdigest(), 'stdout_bytes': len(out)}
    if timed_out or proc.returncode != 0:
        rec['stdout_tail'] = out.decode('utf-8', 'replace')[-500:]
        rec['stderr_tail'] = err.decode('utf-8', 'replace')[-500:]
    return rec


def run_test(name, arms, suffix=''):
    """arms 의 팔을 동시에 띄운다(하나만 주면 혼자 돈다). 팔마다 기록을 돌려준다."""
    procs = {}
    for arm in arms:
        base = os.path.join(SUITE, 'logs', arm, name + suffix)
        fo, fe = open(base + '.out', 'wb'), open(base + '.err', 'wb')
        p = subprocess.Popen(command(name), cwd=arm_dir(arm), env=arm_env(arm), stdout=fo, stderr=fe, stdin=subprocess.DEVNULL,
                             creationflags=subprocess.CREATE_NEW_PROCESS_GROUP)
        procs[arm] = (p, time.perf_counter(), fo, fe)
    recs = {}
    while len(recs) < len(arms):
        for arm, (p, t0, fo, fe) in procs.items():
            if arm in recs:
                continue
            if p.poll() is not None:
                recs[arm] = finish(arm, name + suffix, p, t0, fo, fe, False)
            elif time.perf_counter() - t0 > TIMEOUT:
                kill_tree(p.pid)
                try:
                    p.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    pass
                recs[arm] = finish(arm, name + suffix, p, t0, fo, fe, True)
        time.sleep(0.2)
    return recs


def status(rec):
    return 'timeout' if rec['timeout'] else rec['exit']


def show(st):
    return '시간' if st == 'timeout' else ('통과' if st == 0 else '실패(%s)' % st)


def sha256_file(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def load_progress(key):
    """이어 돌리기: 같은 DLL(SHA-256) · 같은 커밋 · 같은 시험 고르기의 진행 기록이 있으면 돌려준다."""
    if not os.path.isfile(PROGRESS):
        return None
    try:
        with open(PROGRESS, encoding='utf-8') as fh:
            pr = json.load(fh)
    except ValueError:
        return None
    return pr if pr.get('key') == key else None


def save_progress(pr):
    """시험 하나가 끝날 때마다 디스크에 쓴다(기계가 내려가도 끝난 시험은 남는다)."""
    tmp = PROGRESS + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(pr, fh, ensure_ascii=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, PROGRESS)


def clean_states():
    return {n: git('status', '--porcelain', cwd=d).stdout.strip() == '' for n, d in
            (('engine_main', ENGINE_MAIN), ('link_worktree', LINK_DIR), ('base_worktree', BASE_DIR))}


def gate(fail, ok, label, text, limit):
    print('  %-8s %s  %s  (기준 %s)' % (label, 'PASS' if ok else 'FAIL', text, limit), flush=True)
    return fail + (0 if ok else 1)


def main():
    global NEW, SUITE, PROGRESS, OUT_JSON, CHANGED, KNOWN_TIMEOUT, EXTRA_DLL
    argv = sys.argv[1:]
    only, out_json, run = None, None, '1'
    resume = '--resume' in argv
    if resume:
        argv.remove('--resume')
    for flag in ('--only', '--out', '--run'):
        if flag in argv:
            i = argv.index(flag)
            val = argv[i + 1]
            del argv[i:i + 2]
            if flag == '--only':
                only = [s for s in val.split(',') if s]
            elif flag == '--out':
                out_json = os.path.abspath(val)
            else:
                run = val
    if len(argv) != 2 or run not in RUNS:
        print(__doc__)
        return 2
    cfg = RUNS[run]
    NEW, CHANGED, KNOWN_TIMEOUT, EXTRA_DLL = cfg['new'], cfg['changed'], cfg['known_timeout'], cfg['extra_dll']
    SUITE = SUITE + cfg['suffix']
    PROGRESS = os.path.join(SUITE, 'progress.json')
    OUT_JSON = os.path.join(ROOT, 'results', 'gate_engine_suite%s.json' % cfg['suffix'])
    dll_new, dll_base = (os.path.abspath(a) for a in argv)
    missing = [p for p in (dll_new, dll_base, PY) if not os.path.isfile(p)]
    missing += [p for p in (ENGINE_MAIN, LINK_DIR, BASE_DIR) if not os.path.isdir(os.path.join(p, 'drape'))]
    if missing:
        print('없는 것: %s' % ', '.join(missing))
        return 2
    if only is None and out_json is None:
        out_json = OUT_JSON

    t_all = time.perf_counter()
    dq = git('diff', '--name-only', BASE, NEW, '--', *SAME_PATHS, cwd=ENGINE_MAIN)
    changed_found = sorted(s for s in dq.stdout.splitlines() if s.strip()) if dq.returncode == 0 else ['(git diff 실패)']
    same = changed_found == sorted(CHANGED)
    key = {'new_commit': NEW, 'base_commit': BASE, 'dll_new_sha256': sha256_file(dll_new), 'dll_base_sha256': sha256_file(dll_base),
           'only': only, 'run': run}
    pr = load_progress(key) if resume else None
    resumed = pr is not None
    if not resumed:
        prepare(dll_new, dll_base)
        pr = {'key': key, 'segments': [], 'in_progress': None, 'interrupted': {}, 'results': {}, 'seconds': 0.0}
    elif pr.get('in_progress'):
        # 앞 구간이 이 시험을 돌리다 끊겼다(기계가 내려갔다)
        ip = pr['in_progress']
        pr['interrupted'][ip] = pr['interrupted'].get(ip, 0) + 1
        pr['in_progress'] = None
    pr['segments'].append({'started': time.strftime('%Y-%m-%dT%H:%M:%S'), 'resumed': resumed, 'clean_at_start': clean_states()})
    save_progress(pr)
    clean_before = pr['segments'][0]['clean_at_start']
    probes = {arm: probe(arm) for arm in ARMS}
    scripts = os.path.join(arm_dir('N'), 'scripts')
    tests = sorted(f for f in os.listdir(scripts) if f.endswith('.py') and
                   ((f.startswith('py_') and f.endswith('_test.py')) or f.startswith('test_')))
    if only is not None:
        unknown = [n for n in only if n not in tests]
        if unknown:
            print('모르는 시험: %s' % ', '.join(unknown))
            return 2
        tests = [n for n in tests if n in only]
    print('엔진 합성 시험 %s 차 — 시험 %d 개, 팔 %s, 제한 시간 %d s\n  시험 트리 %s (엔진 %s 의 git archive)\n  새 DLL %s\n  종전 DLL %s'
          % (run, len(tests), ' · '.join('%s = %s' % (a, ARM_LABEL[a]) for a in ARMS), TIMEOUT, SUITE, NEW, dll_new, dll_base), flush=True)
    if CHANGED or KNOWN_TIMEOUT or EXTRA_DLL:
        print('  378d20c 와 달라도 되는 파일 %s\n  세 팔 모두 제한 시간이면 다시 돌리지 않는 시험 %s\n  DLL 을 더 놓는 자리 %s'
              % (list(CHANGED) or '없음', list(KNOWN_TIMEOUT) or '없음', list(EXTRA_DLL) or '없음'), flush=True)

    results = pr['results']
    if resumed:
        print('  이어 돌리기 — 앞 구간에서 끝난 시험 %d 개, 구간 %d 번째, 돌리다 끊긴 시험: %s' % (
            sum(1 for n in tests if 'first' in results.get(n, {})), len(pr['segments']), pr['interrupted'] or '없음'), flush=True)

    def guarded(job, fn):
        """job 을 돌리는 동안 기계가 내려가면 다음 --resume 이 그것을 안다. 같은 job 에서 MAX_INTERRUPTS 번 내려갔으면 돌리지 않는다."""
        if pr['interrupted'].get(job, 0) >= MAX_INTERRUPTS:
            return None
        pr['in_progress'] = job
        save_progress(pr)
        out = fn()
        pr['in_progress'] = None
        return out

    for i, name in enumerate(tests, 1):
        if 'first' in results.get(name, {}) or results.get(name, {}).get('machine'):
            continue
        recs = guarded(name, lambda: run_test(name, ARMS))
        if recs is None:
            results[name] = {'machine': True}
            save_progress(pr)
            print('  [%2d/%d] %-44s 건너뜀 — 이 시험을 돌리다 기계가 %d 번 내려갔다' % (i, len(tests), name, MAX_INTERRUPTS), flush=True)
            continue
        results[name] = {'first': recs}
        save_progress(pr)
        print('  [%2d/%d] %-44s %s' % (i, len(tests), name, '  '.join(
            '%s %-9s %6.1fs' % (a, show(status(recs[a])), recs[a]['seconds']) for a in ARMS)), flush=True)

    # 제한 시간 규칙: 어느 팔이든 제한 시간에 걸린 시험은 세 팔을 하나씩 다시 돌려 그 결과로 판정한다.
    for name in tests:
        r = results[name]
        if r.get('machine') or 'rerun' in r or not any(r['first'][a]['timeout'] for a in ARMS):
            continue
        if name in KNOWN_TIMEOUT and all(r['first'][a]['timeout'] for a in ARMS):
            # 앞 차수에서 세 팔 모두 혼자 돌려도 제한 시간을 넘은 시험 — 이번에도 세 팔 모두 걸렸으면 다시 돌리지 않는다(판정 불가 — 시간)
            if not r.get('rerun_skipped'):
                r['rerun_skipped'] = True
                save_progress(pr)
            print('  다시(하나씩) %-36s 안 돌림 — 앞 차수에서 하나씩 돌려도 세 팔 모두 제한 시간, 이번 본 실행도 세 팔 모두 제한 시간' % name, flush=True)
            continue
        rer = {}
        for a in ARMS:
            one = guarded('%s#rerun' % name, lambda: run_test(name, (a,), suffix='.rerun'))
            if one is None:
                break
            rer.update(one)
        if len(rer) < len(ARMS):
            r['machine'] = True
            save_progress(pr)
            print('  다시(하나씩) %-36s 건너뜀 — 다시 돌리다 기계가 %d 번 내려갔다' % (name, MAX_INTERRUPTS), flush=True)
            continue
        r['rerun'] = rer
        save_progress(pr)
        print('  다시(하나씩) %-36s %s' % (name, '  '.join(
            '%s %-9s %6.1fs' % (a, show(status(rer[a])), rer[a]['seconds']) for a in ARMS)), flush=True)
    machine = [n for n in tests if results[n].get('machine')]
    tests_run = [n for n in tests if not results[n].get('machine')]
    for name in tests_run:
        final = results[name].get('rerun', results[name]['first'])
        results[name]['final'] = {a: status(final[a]) for a in ARMS}
        results[name]['undecidable_time'] = all(final[a]['timeout'] for a in ARMS)
        results[name]['stdout_K_eq_B'] = final['K']['stdout_sha256'] == final['B']['stdout_sha256']
        results[name]['stdout_N_eq_K'] = final['N']['stdout_sha256'] == final['K']['stdout_sha256']
    clean_after = clean_states()
    pr['seconds'] += time.perf_counter() - t_all
    save_progress(pr)

    fail = 0
    print(' S0 전제')
    a_ok = True
    for arm in ARMS:
        pb = probes[arm]
        bi = pb.get('build_info', {})
        want = BASE if arm == 'B' else NEW
        inside = os.path.normcase(os.path.abspath(pb.get('wrapper', ''))).startswith(os.path.normcase(arm_dir(arm))) and \
            os.path.normcase(os.path.abspath(pb.get('dll', ''))).startswith(os.path.normcase(arm_dir(arm)))
        ok = str(bi.get('git_commit', '')).startswith(want) and bi.get('dirty') is False and inside
        a_ok &= ok
        fail = gate(fail, ok, '(a) ' + arm, 'DLL git_commit %s · dirty %s · 래퍼 · DLL 이 사본 안 %s%s' % (
            bi.get('git_commit'), bi.get('dirty'), inside, (' · ' + pb['error'][-200:]) if 'error' in pb else ''),
            '%s · dirty false · 사본 안' % want)
    fail = gate(fail, same, '(b)', '시험 스크립트 · 래퍼에서 다른 파일 %s' % (changed_found or '없음'),
                '%s 와 %s 에서 다른 파일이 %s' % (BASE, NEW, ('정확히 ' + ', '.join(CHANGED)) if CHANGED else '없음'))
    starts_clean = all(all(s['clean_at_start'].values()) for s in pr['segments'])
    c_ok = starts_clean and all(clean_after.values())
    fail = gate(fail, c_ok, '(c)', '엔진 워크트리 셋 git status — 구간 %d 개의 시작 %s · 끝 %s' % (
        len(pr['segments']), '비어 있음' if starts_clean else [s['clean_at_start'] for s in pr['segments']],
        '비어 있음' if all(clean_after.values()) else clean_after), '비어 있음')

    judged = [n for n in tests_run if not results[n]['undecidable_time']]
    s1_bad = [n for n in judged if results[n]['final']['K'] != results[n]['final']['B']]
    s2_base = [n for n in judged if results[n]['final']['K'] == 0]
    s2_bad = [n for n in s2_base if results[n]['final']['N'] != 0]
    print(' S1 킬스위치 = 종전 (K 의 종료 상태 = B 의 종료 상태)')
    fail = gate(fail, not s1_bad, 'S1', '다른 시험 %d / %d%s' % (len(s1_bad), len(judged), (': ' + ', '.join(s1_bad)) if s1_bad else ''), '0')
    print(' S2 새 기본값이 통과를 잃지 않는다 (K 가 통과한 시험에서 N 도 통과)')
    fail = gate(fail, not s2_bad, 'S2', 'K 통과 %d 개 중 N 이 못 넘은 것 %d%s' % (
        len(s2_base), len(s2_bad), (': ' + ', '.join(s2_bad)) if s2_bad else ''), '0')

    all_fail = [n for n in judged if all(results[n]['final'][a] != 0 for a in ARMS)]
    n_only = [n for n in judged if results[n]['final']['N'] == 0 and results[n]['final']['K'] != 0]
    undec = [n for n in tests_run if results[n]['undecidable_time']]
    print(' 보고')
    print('  세 팔 모두 통과 %d · 세 팔 모두 실패 %d · N 만 통과 %d · 판정 불가(시간) %d · 판정 불가(기계) %d' % (
        sum(1 for n in judged if all(results[n]['final'][a] == 0 for a in ARMS)), len(all_fail), len(n_only), len(undec), len(machine)))
    for n in machine:
        print('   판정 불가(기계)  %s' % n)
    if len(pr['segments']) > 1 or pr['interrupted']:
        print('  구간 %d 개(기계가 내려가 이어 돌림) — 돌리다 끊긴 시험: %s' % (len(pr['segments']), pr['interrupted'] or '없음'))
    for n in all_fail:
        f = results[n].get('rerun', results[n]['first'])['K']
        tail = (f.get('stderr_tail') or f.get('stdout_tail') or '').strip().splitlines()
        print('   모두 실패  %-42s K %s  %s' % (n, show(results[n]['final']['K']), tail[-1][:110] if tail else ''))
    for n in n_only:
        print('   N 만 통과  %s (K %s · B %s)' % (n, show(results[n]['final']['K']), show(results[n]['final']['B'])))
    for n in undec:
        print('   판정 불가(시간)  %s' % n)
    print('  표준 출력: K = B 인 시험 %d / %d · N = K 인 시험 %d / %d' % (
        sum(1 for n in tests_run if results[n]['stdout_K_eq_B']), len(tests_run),
        sum(1 for n in tests_run if results[n]['stdout_N_eq_K']), len(tests_run)))
    secs = pr['seconds']
    print('  실행 시간 %.1f 분 (구간 %d 개의 합)' % (secs / 60, len(pr['segments'])))

    verdict = 'PASS' if fail == 0 else 'FAIL'
    res = {'gate': 'engine_suite', 'run': run, 'verdict': verdict, 'gates_failed': fail,
           'run_config': {'changed_allowed': list(CHANGED), 'changed_found': changed_found, 'known_timeout': list(KNOWN_TIMEOUT),
                          'extra_dll': list(EXTRA_DLL), 'suite_dir': SUITE.replace('\\', '/')},
           'engine_main': ENGINE_MAIN.replace('\\', '/'), 'new_commit': NEW, 'base_commit': BASE, 'timeout_s': TIMEOUT,
           'python': PY.replace('\\', '/'), 'arms': {'N': 'new DLL, no switch', 'K': 'new DLL + ' + ' '.join('%s=%s' % kv for kv in KILL.items()),
                                                     'B': 'base DLL'},
           'S0': {'probes': probes, 'a_ok': a_ok, 'scripts_same': same, 'clean_before': clean_before, 'clean_after': clean_after},
           'S1': {'judged': len(judged), 'differ': s1_bad}, 'S2': {'k_passed': len(s2_base), 'n_lost': s2_bad},
           'report': {'all_fail': all_fail, 'n_only_pass': n_only, 'undecidable_time': undec, 'undecidable_machine': machine,
                      'stdout_K_eq_B': sum(1 for n in tests_run if results[n]['stdout_K_eq_B']),
                      'stdout_N_eq_K': sum(1 for n in tests_run if results[n]['stdout_N_eq_K']), 'minutes': round(secs / 60, 1)},
           'segments': pr['segments'], 'interrupted': pr['interrupted'], 'dll_sha256': {'new': key['dll_new_sha256'], 'base': key['dll_base_sha256']},
           'only': only, 'tests': results}
    if out_json:
        with open(out_json, 'w', encoding='utf-8') as fh:
            json.dump(res, fh, indent=2, ensure_ascii=False)
    print('\n엔진 합성 시험 %s (실패 %d)%s' % (verdict, fail, (' → ' + os.path.relpath(out_json, ROOT)) if out_json else ' (결과 json 안 씀 — 하네스 확인용)'))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
