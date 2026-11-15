"""게이트 표본: 무작위 힌지 1000개를 결정적으로 만든다(시드 고정).

실제 의복 스케일을 넣는다(규칙 §1 「실제 치수」): 변 길이 0.5 – 3 cm, 좌표는 원점에서 최대 1 m 떨어진 곳
(인체 위 의복). 접힘각 θ 는 −170° … 170° (atan2 의 ±π 가지 자름을 유한차분이 넘지 않게).
처음 4개는 평평(θ = 0) · ±90° · 150° 로 고정한다.

한 줄 = x_h x_i x_j x_k (12) · θ0 (1) · v_h v_i v_j v_k (12), repr 표기(왕복 정밀도).
사용: python model/make_hinges.py results/hinges.txt
"""
import math
import random
import sys


def rot_from_quat(q):
    w, x, y, z = q
    return [[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
            [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
            [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]]


def main(path, n=1000, seed=20261112):
    rnd = random.Random(seed)
    fixed = [0.0, 90.0, -90.0, 150.0]
    out = []
    for s in range(n):
        L = rnd.uniform(0.005, 0.03)
        u1, w1 = rnd.uniform(-0.2, 1.2) * L, rnd.uniform(0.3, 1.2) * L
        u2, w2 = rnd.uniform(-0.2, 1.2) * L, rnd.uniform(0.3, 1.2) * L
        phi = math.radians(fixed[s] if s < len(fixed) else rnd.uniform(-170.0, 170.0))
        pts = [(u1, -w1, 0.0), (0.0, 0.0, 0.0), (L, 0.0, 0.0), (u2, w2 * math.cos(phi), w2 * math.sin(phi))]
        q = [rnd.gauss(0, 1) for _ in range(4)]
        qn = math.sqrt(sum(c * c for c in q))
        R = rot_from_quat([c / qn for c in q])
        t = [rnd.uniform(-1.0, 1.0) for _ in range(3)]
        xs = []
        for p in pts:
            for r in range(3):
                xs.append(R[r][0] * p[0] + R[r][1] * p[1] + R[r][2] * p[2] + t[r])
        theta0 = rnd.uniform(-0.3, 0.3)
        vs = [rnd.uniform(-0.2, 0.2) for _ in range(12)]
        out.append(' '.join(repr(v) for v in xs + [theta0] + vs))
    with open(path, 'w', newline='\n') as f:
        f.write('\n'.join(out) + '\n')
    print('wrote %d hinges to %s' % (n, path))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'results/hinges.txt')
