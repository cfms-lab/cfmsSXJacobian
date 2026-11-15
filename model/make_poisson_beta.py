"""Poisson 힘 게이트 표본: 삼각형마다 Poisson 계수 β 를 결정적으로 만든다(시드 고정). docs/GATE_POISSON.md.

results/triangles.txt(model/make_triangles.py)의 줄과 짝이다. 처음 4 개(고정 표본)는 1 · 0 · −1 · 0.6,
그 뒤 짝수 줄(0 부터 셈)은 1 · 0 · −1 을 돌려 쓰고 홀수 줄은 [−1, 1] 의 무작위 값이다(β 가 실제 변수로 들어가는지 본다).
원문(Sul 외 2010)의 β = 1 · 0 · −1 은 넓이 보존 · 가로 변형 없음 · 가로로 같이 늘어남이다.

한 줄 = β, repr 표기(왕복 정밀도).
사용: python model/make_poisson_beta.py results/poisson_beta.txt
"""
import random
import sys


def main(path, n=1000, seed=20261115):
    rnd = random.Random(seed)
    fixed = [1.0, 0.0, -1.0, 0.6]
    cycle = [1.0, 0.0, -1.0]
    out = []
    for s in range(n):
        beta = rnd.uniform(-1.0, 1.0)            # 줄마다 뽑아 두고(순서 고정) 홀수 줄만 쓴다
        if s < len(fixed):
            beta = fixed[s]
        elif s % 2 == 0:
            beta = cycle[(s // 2) % 3]
        out.append('%r' % beta)
    with open(path, 'w', newline='\n') as f:
        f.write('\n'.join(out) + '\n')
    print('wrote %d Poisson coefficients to %s' % (n, path))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'results/poisson_beta.txt')
