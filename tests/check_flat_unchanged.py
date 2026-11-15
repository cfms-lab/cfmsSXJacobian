"""속도 게이트 S3 — `ctx.stage` 표시가 평탄판을 바꾸지 않았는가(docs/GATE_SPEED.md).

생성기가 results/codegen_BaraffBend.json 에 적은 함수 본문 SHA-256 을, 착수 게이트(커밋 64a2f95) 때
생성된 generated/BaraffBend.h 의 Eval · EvalFrozen 본문에서 잰 값과 견준다.
사용: python tests/check_flat_unchanged.py   (실패하면 종료 코드 1)
"""
import json
import sys

GATE0 = {   # 2026-11-12, 착수 게이트 판 generated/BaraffBend.h (파일 sha256 2fa13f73...)
    'exact': '06fadd36f3efca38c82786346ee7631a02032ca4bcbbc542050420258e970093',    # Eval
    'frozen': 'd46888e7fedd6d5b9219ca0df6c2057f3964ee60a47868a75da5fa9116e64c49',   # EvalFrozen
}


def main(path='results/codegen_BaraffBend.json'):
    st = json.load(open(path, encoding='utf-8'))['variants']
    fail = 0
    print('S3 평탄판 불변 (함수 본문 SHA-256, 착수 게이트 판과 비교)')
    for mode, want in GATE0.items():
        got = st.get(mode, {}).get('body_sha256')
        ok = got == want
        fail += not ok
        print('  %-6s %s  %s' % (mode, 'PASS' if ok else 'FAIL', (got or 'none')[:16]))
    return 1 if fail else 0


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
