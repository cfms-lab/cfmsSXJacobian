"""논문용 측정 4 의 판정 뒤 진단 (docs/GATE_JACOBIAN_DAMPING.md 「논문용 측정 4 의 결과」): 주름의 방향.

측정 4 의 접힌 천은 주름 선이 고정 변과 수직(y = 20 + A sin(2π·2·x/L))이었고, 블록 대칭 Gauss–Newton 이 고정 변을 축으로 한
강체 회전을 막지 않았다(R ≈ 1e-13). 식으로는 그 방향이라 대칭화한 항이 사라진다. 주름 선을 고정 변과 나란히 둔 판
(y = 20 + A sin(2π·2·(z + L/2)/L))에서 같은 것을 잰다. 적분 없음.

사용:  python tests/diag_symgn_folded.py <실험 cfmsDrape.dll>     (리포 루트에서, run_gate 밖 — 약 1 분)
결과 results/diag_symgn_folded.json.
"""
import sys
sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import tempfile  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import measure_jacobian_symgn as M  # noqa: E402

OUT = os.path.join(M.G.ROOT, 'results', 'diag_symgn_folded.json')


def main():
    dll = os.path.abspath(sys.argv[1])
    td = tempfile.mkdtemp(prefix='cfmsSXJ_symgn_fold_')
    res = {'dll': dll, 'n': M.FOLD_N, 'BE_U': M.BEND, 'amps_cm': list(M.FOLD_AMP), 'orient': {}}
    try:
        for orient in ('x', 'z'):
            res['orient'][orient] = {arm: M.spawn_folded(td, arm, dict(dll=dll, n=M.FOLD_N, bend=M.BEND, amps=list(M.FOLD_AMP), orient=orient))
                                     for arm in M.ARMS}
    finally:
        shutil.rmtree(td, ignore_errors=True)
    for orient, label in (('x', '주름 선이 고정 변과 수직 (측정 4)'), ('z', '주름 선이 고정 변과 나란히')):
        print(label)
        for k, amp in enumerate(M.FOLD_AMP):
            c0 = res['orient'][orient]['off']['cases'][k]
            print('  A %.2f cm (접힌 각 rms %.1f° · 최대 %.1f°)' % (amp, c0['theta_rms_deg'], c0['theta_max_deg']))
            for arm in M.ARMS:
                c = res['orient'][orient][arm]['cases'][k]
                print('     %-12s R %.6g   고윳값 비 %.3e' % (arm, c['rayleigh_rot'], c['eig_min'] / c['eig_max']))
    res['engine'] = {kk: res['orient']['x']['off']['build_info'].get(kk) for kk in ('git_commit', 'dirty')}
    for orient in res['orient'].values():
        for arm in orient.values():
            arm.pop('build_info', None)
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
        fh.write('\n')
    print('→ %s' % os.path.relpath(OUT, M.G.ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
