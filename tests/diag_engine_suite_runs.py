"""판정 밖 진단: 엔진 합성 시험 1 차 · 2 차에서 팔마다 종료 상태 · 표준 출력이 같은가 (docs/GATE_ENGINE_SUITE.md 「2 차」의 「보고만(더함)」).

사용:  python tests/diag_engine_suite_runs.py
두 실행의 결과 json(results/gate_engine_suite.json · gate_engine_suite_2.json)과 로그(%TEMP%\\cfmsSXJ_engine_suite\\logs ·
%TEMP%\\cfmsSXJ_engine_suite_2\\logs)를 읽는다. 엔진 리포는 읽지 않는다.
K · B 는 고친 경로(생성 굽힘의 병렬 gather)를 밟지 않으므로 두 차수에서 같은 값을 내야 한다. N 은 2 스레드 이상에서 달라질 수 있다.
견주기 전에 맞추는 것(값이 아닌 줄): 시험 트리 폴더 이름, 빌드 정보 줄(커밋 · DLL 해시), 걸린 시간(`1.23s` · `(x0.95)`),
고친 시험(py_membrane_nonlinear_test)의 기준 2 이름표. 어느 차수에서든 제한 시간에 걸린 팔은 견주지 않는다(출력이 잘린다).
결과는 results/diag_engine_suite_runs.json.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402
import tempfile  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP = os.environ.get('TEMP') or tempfile.gettempdir()
RUNS = {'1': ('gate_engine_suite.json', 'cfmsSXJ_engine_suite'), '2': ('gate_engine_suite_2.json', 'cfmsSXJ_engine_suite_2')}
ARMS = ('N', 'K', 'B')
TIME = re.compile(r'\d+\.\d+ ?s\b')
RATIO = re.compile(r'\(x\d+\.\d+\)')


def norm(raw):
    out = []
    for line in raw.decode('utf-8', 'replace').replace('\r', '').split('\n'):
        if 'git_commit' in line:
            continue
        line = line.replace('cfmsSXJ_engine_suite_2', 'cfmsSXJ_engine_suite')
        line = line.replace(' bitwise = OFF', ' = OFF to 1e-12')
        line = RATIO.sub('(xR)', TIME.sub('Ts', line))
        out.append(line)
    return out


def main():
    res, logs = {}, {}
    for r, (js, d) in RUNS.items():
        with open(os.path.join(ROOT, 'results', js), encoding='utf-8') as fh:
            res[r] = json.load(fh)['tests']
        logs[r] = os.path.join(TMP, d, 'logs')
    names = sorted(set(res['1']) & set(res['2']))
    have_logs = all(os.path.isdir(p) for p in logs.values())
    out = {'tests': len(names), 'logs_found': have_logs, 'arms': {}}
    for arm in ARMS:
        a = {'compared': 0, 'skipped_timeout': [], 'exit_differs': [], 'stdout_same_raw': 0, 'stdout_same_normalized': 0,
             'stdout_differs': {}}
        for n in names:
            f1, f2 = res['1'][n].get('first', {}).get(arm), res['2'][n].get('first', {}).get(arm)
            if not f1 or not f2:
                continue
            if f1['timeout'] or f2['timeout']:
                a['skipped_timeout'].append(n)
                continue
            a['compared'] += 1
            if f1['exit'] != f2['exit']:
                a['exit_differs'].append({'test': n, 'run1': f1['exit'], 'run2': f2['exit']})
            if f1['stdout_sha256'] == f2['stdout_sha256']:
                a['stdout_same_raw'] += 1
                a['stdout_same_normalized'] += 1
                continue
            if not have_logs:
                continue
            l1 = norm(open(os.path.join(logs['1'], arm, n + '.out'), 'rb').read())
            l2 = norm(open(os.path.join(logs['2'], arm, n + '.out'), 'rb').read())
            if l1 == l2:
                a['stdout_same_normalized'] += 1
                continue
            pairs = [(x, y) for x, y in zip(l1, l2) if x != y]
            a['stdout_differs'][n] = {'lines_run1': len(l1), 'lines_run2': len(l2), 'differing_lines': len(pairs) + abs(len(l1) - len(l2)),
                                      'first': [{'run1': x[:200], 'run2': y[:200]} for x, y in pairs[:3]]}
        out['arms'][arm] = a
        print('%s: 견준 시험 %d (제한 시간이라 뺀 것 %d) · 종료 상태 다른 것 %d · 표준 출력 그대로 같음 %d · 맞춘 뒤 같음 %d · 다른 것 %d' % (
            arm, a['compared'], len(a['skipped_timeout']), len(a['exit_differs']), a['stdout_same_raw'], a['stdout_same_normalized'],
            len(a['stdout_differs'])))
        for e in a['exit_differs']:
            print('   종료 상태  %-44s 1 차 %s · 2 차 %s' % (e['test'], e['run1'], e['run2']))
        for n, d in a['stdout_differs'].items():
            print('   출력 다름  %-44s 다른 줄 %d' % (n, d['differing_lines']))
            for p in d['first'][:2]:
                print('      1 차: %s\n      2 차: %s' % (p['run1'][:150], p['run2'][:150]))
    with open(os.path.join(ROOT, 'results', 'diag_engine_suite_runs.json'), 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)
    print('→ results/diag_engine_suite_runs.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
