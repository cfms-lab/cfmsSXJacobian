"""엔진 연결 L1-0 (docs/GATE_LINK_L1.md): 생성기에 엔진 변형(EvalEngine*)을 더해도 기존 출력이 그대로인지 본다.

- 하네스 헤더 generated/<이름>.h 다섯의 SHA-256 이 판정문에 고정한 값과 같다.
- results/codegen_<이름>.json 의 기존 변형 body_sha256 이 기준 커밋(55c9745)의 것과 같다.
- 엔진 헤더 generated/<이름>_Engine.h 가 있고, json 에 engine · engine_gn 이 있으며, engine_tension 은 인장에만 있다.
- (2026-11-13 덧붙임, docs/GATE_LINK_L4.md 「2 차 (다)」) 감쇠 없는 GN engine_gn_undamped 가 다섯 모두 있고, 기존 엔진 변형
  (engine · engine_gn · engine_tension)의 body_sha256 이 그 변형을 더하기 전 커밋(df547d8)의 것과 같다.
- (2026-11-13 덧붙임, docs/GATE_LINK_L1.md 「L1-5」) 엔진 처방 T 판 engine_tension_rate 가 인장에만 있고(다섯 가운데), 그것을 더하기 전
  커밋(816a73d)의 엔진 변형(engine_gn_undamped 포함)의 body_sha256 이 그대로다.
기준 json 은 git show 로 읽는다(리포 루트에서 실행). 통과 0, 실패 1.
"""
import hashlib
import json
import os
import subprocess
import sys

BASE = '55c9745'
BASE_ENGINE = 'df547d8'   # engine_gn_undamped 를 더하기 전
BASE_RATE = '816a73d'     # engine_tension_rate 를 더하기 전
HEADER_SHA = {
    'BaraffBend': '38b552fdeaad573dee7a475373a68391c5777f4b0e8a470c0fee3068f16d5649',
    'BaraffBendCompact': '90fb1230301a5ea67d8d96624a44d35fdd68cdcc19b6ca2d7d48071abf60fad0',
    'BaraffStretch': 'b28f4fe6b8ccb9bca602cc3b08c86b2ba45e5c2281a990fae3561313d737bd0c',
    'BaraffShear': '5d086b5e84c914d0f2f786337dea728f014d1fc40ac98538b62cc5c51f556a60',
    'BridsonBend': '2d9a7bd5a47b1b5c285ba8e7c9b258652fda14427db84efdab9270619d63794d',
}
NEW_VARIANTS = ['engine', 'engine_gn', 'engine_tension', 'engine_gn_undamped', 'engine_tension_rate']


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def main():
    fail = 0
    for name, want in HEADER_SHA.items():
        hdr = 'generated/%s.h' % name
        eng = 'generated/%s_Engine.h' % name
        path = 'results/codegen_%s.json' % name
        old = json.loads(subprocess.run(['git', 'show', '%s:%s' % (BASE, path)], capture_output=True, check=True).stdout)
        old_eng = json.loads(subprocess.run(['git', 'show', '%s:%s' % (BASE_ENGINE, path)], capture_output=True,
                                            check=True).stdout)
        new = json.load(open(path, encoding='utf-8'))
        hdr_ok = os.path.isfile(hdr) and sha(hdr) == want
        old_vars = list(old['variants'])
        bad = [v for v in old_vars if v not in new['variants'] or old['variants'][v]['body_sha256'] != new['variants'][v]['body_sha256']]
        old_rate = json.loads(subprocess.run(['git', 'show', '%s:%s' % (BASE_RATE, path)], capture_output=True,
                                             check=True).stdout)
        want_eng = (['engine', 'engine_gn'] + (['engine_tension'] if name == 'BaraffStretch' else []) + ['engine_gn_undamped']
                    + (['engine_tension_rate'] if name == 'BaraffStretch' else []))
        eng_old = [v for v in old_eng['variants'] if v.startswith('engine')]
        bad += [v for v in eng_old if v not in new['variants']
                or old_eng['variants'][v]['body_sha256'] != new['variants'][v]['body_sha256']]
        bad += [v for v in old_rate['variants'] if v.startswith('engine') and v not in bad and (
            v not in new['variants'] or old_rate['variants'][v]['body_sha256'] != new['variants'][v]['body_sha256'])]
        has_eng = [v for v in NEW_VARIANTS if v in new['variants']]
        eng_ok = os.path.isfile(eng) and new.get('engine_header_bytes') == os.path.getsize(eng) and has_eng == want_eng
        extra = [v for v in new['variants'] if v not in old_vars and v not in NEW_VARIANTS]
        ok = hdr_ok and not bad and eng_ok and not extra
        fail += 0 if ok else 1
        print('  %-18s %s  헤더 %s · 기존 변형 %d 개 %s · 엔진 변형 %s%s' % (
            name, 'PASS' if ok else 'FAIL', '같음' if hdr_ok else '다름', len(old_vars) + len(eng_old),
            '같음' if not bad else '다름 ' + ','.join(bad), ','.join(has_eng) or '없음',
            '' if eng_ok else ' (엔진 헤더 · 변형이 판정문과 다름)'))
    print('L1-0 %s (기준 %s · 엔진 변형 %s · %s)' % ('PASS' if fail == 0 else 'FAIL', BASE, BASE_ENGINE, BASE_RATE))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
