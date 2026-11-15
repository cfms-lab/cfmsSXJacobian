"""인장 · 전단 게이트 표본: 무작위 삼각형 1000 개를 결정적으로 만든다(시드 고정). docs/GATE_STRETCH_SHEAR.md.

쉬는 uv 삼각형: 첫 변 0.5 – 3 cm, 둘째 변은 첫 변의 0.5 – 1.5 배, 사잇각 30 – 150°. uv 평면에서 무작위로 돌리고
±1 m 안으로 옮긴다(패턴 좌표). 절반은 v 를 뒤집어 감김 방향을 바꾼다.
변형: x = R F (uv − uv0) + t,  F = [[s_u, γ], [0, s_v], [0, 0]],  s_u · s_v ∈ [0.8, 1.3],  γ ∈ [−0.4, 0.4],
R 무작위 회전, t ∈ ±1 m. 그러면 w = R F. 처음 4 개는 쉬는 모양 · u 방향 1.2 배 · 순수 전단 γ = 0.3 · 압축(0.85, 0.9).

한 줄 = x0 x1 x2 (9) · u0 v0 u1 v1 u2 v2 (6) · v0 v1 v2 (9), repr 표기(왕복 정밀도).
사용: python model/make_triangles.py results/triangles.txt
"""
import math
import random
import sys


def rot_from_quat(q):
    w, x, y, z = q
    return [[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
            [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
            [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]]


def main(path, n=1000, seed=20261113):
    rnd = random.Random(seed)
    fixed = [(1.0, 1.0, 0.0), (1.2, 1.0, 0.0), (1.0, 1.0, 0.3), (0.85, 0.9, 0.0)]   # (s_u, s_v, γ)
    out = []
    for s in range(n):
        L = rnd.uniform(0.005, 0.03)
        r = rnd.uniform(0.5, 1.5) * L
        phi = math.radians(rnd.uniform(30.0, 150.0))
        tri = [(0.0, 0.0), (L, 0.0), (r * math.cos(phi), r * math.sin(phi))]
        psi = rnd.uniform(0.0, 2.0 * math.pi)
        c, sn = math.cos(psi), math.sin(psi)
        off = (rnd.uniform(-1.0, 1.0), rnd.uniform(-1.0, 1.0))
        flip = -1.0 if rnd.random() < 0.5 else 1.0
        uv = [(c * p[0] - sn * p[1] + off[0], flip * (sn * p[0] + c * p[1]) + off[1]) for p in tri]

        su, sv, gam = fixed[s] if s < len(fixed) else (rnd.uniform(0.8, 1.3), rnd.uniform(0.8, 1.3), rnd.uniform(-0.4, 0.4))
        q = [rnd.gauss(0, 1) for _ in range(4)]
        qn = math.sqrt(sum(c * c for c in q))
        R = rot_from_quat([c / qn for c in q])
        t = [rnd.uniform(-1.0, 1.0) for _ in range(3)]
        xs = []
        for p in uv:
            du, dv = p[0] - uv[0][0], p[1] - uv[0][1]
            loc = (su * du + gam * dv, sv * dv, 0.0)          # F (uv − uv0)
            for k in range(3):
                xs.append(R[k][0] * loc[0] + R[k][1] * loc[1] + R[k][2] * loc[2] + t[k])
        uvs = [c for p in uv for c in p]
        vs = [rnd.uniform(-0.2, 0.2) for _ in range(9)]
        out.append(' '.join(repr(v) for v in xs + uvs + vs))
    with open(path, 'w', newline='\n') as f:
        f.write('\n'.join(out) + '\n')
    print('wrote %d triangles to %s' % (n, path))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'results/triangles.txt')
