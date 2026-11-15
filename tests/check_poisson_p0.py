"""Poisson 힘 게이트 Q0 (docs/GATE_POISSON.md): 생성기에 log 출력(_print_log)을 더해도 기존 명세의 출력이 그대로인지 본다.

- 명세 여섯의 results/codegen_<이름>.json 에서 모든 변형의 body_sha256 과 header_bytes · engine_header_bytes 가
  기준 커밋(4798bc3)의 것과 같다.
- SulPoisson 에는 인장 판 변형(implicit_tension · engine_tension · engine_tension_rate)이 없다.
- 생성 헤더 generated/SulPoisson.h · SulPoisson_Engine.h 에 std::log 가 있고, std:: 없는 log( 가 없다.
기준 json 은 git show 로 읽는다(리포 루트에서 실행). 통과 0, 실패 1.
"""
import json
import re
import subprocess
import sys

BASE = '4798bc3'
SPECS = ['BaraffBend', 'BaraffBendCompact', 'BaraffStretch', 'BaraffShear', 'BridsonBend', 'ProvotShear']
NEW = 'SulPoisson'
NEW_VARIANTS = ['exact', 'frozen', 'staged', 'implicit', 'implicit_gn', 'engine', 'engine_gn', 'engine_gn_undamped']
BARE_LOG = re.compile(r'(?<![:\w])log\(')


def main():
    fail = 0
    for name in SPECS:
        path = 'results/codegen_%s.json' % name
        old = json.loads(subprocess.run(['git', 'show', '%s:%s' % (BASE, path)], capture_output=True, check=True).stdout)
        new = json.load(open(path, encoding='utf-8'))
        bad = [v for v in old['variants'] if v not in new['variants']
               or old['variants'][v]['body_sha256'] != new['variants'][v]['body_sha256']]
        extra = [v for v in new['variants'] if v not in old['variants']]
        size_ok = all(old[k] == new.get(k) for k in ('header_bytes', 'engine_header_bytes'))
        ok = not bad and not extra and size_ok
        fail += 0 if ok else 1
        print('  %-18s %s  변형 %d 개 %s · 헤더 크기 %s' % (
            name, 'PASS' if ok else 'FAIL', len(old['variants']),
            '같음' if not bad and not extra else '다름 ' + ','.join(bad + extra), '같음' if size_ok else '다름'))

    new = json.load(open('results/codegen_%s.json' % NEW, encoding='utf-8'))
    var_ok = list(new['variants']) == NEW_VARIANTS
    n_std, n_bare = 0, 0
    for hdr in ('generated/%s.h' % NEW, 'generated/%s_Engine.h' % NEW):
        text = open(hdr, encoding='utf-8').read()
        n_std += text.count('std::log(')
        n_bare += len(BARE_LOG.findall(text))
    log_ok = n_std > 0 and n_bare == 0
    ok = var_ok and log_ok
    fail += 0 if ok else 1
    print('  %-18s %s  변형 %s · std::log %d 곳 · std:: 없는 log( %d 곳' % (
        NEW, 'PASS' if ok else 'FAIL', '판정문대로' if var_ok else '다름 ' + ','.join(new['variants']), n_std, n_bare))
    print('Q0 %s (기준 %s)' % ('PASS' if fail == 0 else 'FAIL', BASE))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
