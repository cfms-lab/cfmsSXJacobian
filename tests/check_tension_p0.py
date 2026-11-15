"""인장 C > 0 판 게이트 P0 (docs/GATE_TENSION.md): 생성기에 EvalImplicitTension 을 더해도 기존 변형이 그대로인지 본다.

명세 다섯의 기존 변형(exact · frozen · staged · implicit · implicit_gn) body_sha256 이 기준 커밋(ada7e61)의
results/codegen_*.json 과 같고, 인장이 아닌 넷은 header_bytes 도 같아야 한다. 인장 명세에만 implicit_tension 이 있어야 한다.
기준값은 git show 로 읽는다(리포 루트에서 실행). 통과 0, 실패 1.
"""
import json
import subprocess
import sys

BASE = 'ada7e61'
SPECS = ['BaraffBend', 'BaraffBendCompact', 'BaraffStretch', 'BaraffShear', 'BridsonBend']
OLD_VARIANTS = ['exact', 'frozen', 'staged', 'implicit', 'implicit_gn']


def main():
    fail = 0
    for name in SPECS:
        path = 'results/codegen_%s.json' % name
        old = json.loads(subprocess.run(['git', 'show', '%s:%s' % (BASE, path)], capture_output=True, check=True).stdout)
        new = json.load(open(path, encoding='utf-8'))
        bad = [v for v in OLD_VARIANTS if old['variants'][v]['body_sha256'] != new['variants'][v]['body_sha256']]
        has_t = 'implicit_tension' in new['variants']
        want_t = name == 'BaraffStretch'
        hdr_ok = want_t or old['header_bytes'] == new['header_bytes']
        ok = not bad and has_t == want_t and hdr_ok
        fail += 0 if ok else 1
        print('  %-18s %s  기존 변형 %s%s · implicit_tension %s · header %s' % (
            name, 'PASS' if ok else 'FAIL', '같음' if not bad else '다름 ' + ','.join(bad),
            '', '있음' if has_t else '없음', '같음' if old['header_bytes'] == new['header_bytes'] else
            '%d → %d 바이트' % (old['header_bytes'], new['header_bytes'])))
    print('P0 %s (기준 %s)' % ('PASS' if fail == 0 else 'FAIL', BASE))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
