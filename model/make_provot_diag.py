"""Provot 전단 게이트 표본: 삼각형마다 쉬는 대각선 길이 b₊ b₋ 를 결정적으로 만든다(시드 고정). docs/GATE_PROVOT_SHEAR.md.

results/triangles.txt(model/make_triangles.py)의 줄과 짝이다. 짝수 줄(0 부터 셈)과 처음 4 개(고정 표본)는 √2 √2 —
직교 · 단위 틀의 쉬는 대각선 길이다. 홀수 줄은 √2 × [0.9, 1.1] 의 무작위 값 둘이다(b 가 실제 변수로 들어가는지 본다).

한 줄 = b₊ b₋, repr 표기(왕복 정밀도).
사용: python model/make_provot_diag.py results/provot_diag.txt
"""
import math
import random
import sys


def main(path, n=1000, seed=20261114):
    rnd = random.Random(seed)
    r2 = math.sqrt(2.0)
    out = []
    for s in range(n):
        bp, bm = r2 * rnd.uniform(0.9, 1.1), r2 * rnd.uniform(0.9, 1.1)   # 줄마다 뽑아 두고(순서 고정) 홀수 줄만 쓴다
        if s < 4 or s % 2 == 0:
            bp, bm = r2, r2
        out.append('%r %r' % (bp, bm))
    with open(path, 'w', newline='\n') as f:
        f.write('\n'.join(out) + '\n')
    print('wrote %d diagonal rest lengths to %s' % (n, path))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'results/provot_diag.txt')
