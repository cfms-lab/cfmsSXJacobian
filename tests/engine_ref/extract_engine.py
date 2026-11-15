"""엔진 손 코드 비교 게이트 E1 — 엔진 Macri 굽힘 코드를 바이트 그대로 뽑는다(docs/GATE_ENGINE.md).

엔진 리포(cfmsDrape_dev)는 읽기만 한다: `git show <커밋>:<경로>` 로 고정 커밋의 원본을 읽고, 필요한 함수 정의 ·
선언 · 블록을 **이름과 표지 줄로 찾아** 그 바이트를 그대로 generated/engine_ref/ 에 쓴다(git 에 넣지 않는다 —
엔진 코드 사본을 이 리포에 싣지 않기 위해서). 조각마다 SHA-256 을 tests/engine_ref/manifest.json 의 고정값과
견준다.

함수 정의 = 0 열에서 시작하는 정의 줄부터 그 뒤 처음 나오는 「}」 한 줄까지(엔진 코드의 함수 끝은 모두 0 열 「}」).
조각 사이에는 이쪽이 단 구분 주석 한 줄만 들어간다. 해시는 조각 바이트(구분 주석 제외)로 잰다.

사용:  python tests/engine_ref/extract_engine.py          뽑고 고정값과 견준다(E1)
       python tests/engine_ref/extract_engine.py --pin    뽑고 고정값을 새로 쓴다(처음 한 번, 또는 커밋을 바꿀 때 —
                                                          그때는 판정문에 이유를 적는다)
엔진 리포 위치 = 환경 변수 CFMS_DRAPE_DIR, 없으면 이 리포 옆의 cfmsDrape_dev.
종료 코드: 0 통과, 1 E1 실패, 3 엔진 리포 없음(게이트를 건너뛴다).
"""
import hashlib
import json
import os
import re
import subprocess
import sys

COMMIT = '378d20c'
SRC = 'drape/src/sde_core/'
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, 'generated', 'engine_ref')
MANIFEST = os.path.join(ROOT, 'tests', 'engine_ref', 'manifest.json')

# SDECore.cpp 에서 정의를 뽑을 멤버 함수와 정의 개수(오버로드). 파일 순서대로 쓴다.
MEMBER_FUNCS = [
    ('CheckBendingChange_Macri_SOA', 1),
    ('_Input_b012', 1),
    ('_Tilde_Operation', 1),
    ('_Dot31_M33_by_V31', 3),
    ('_TripleProduct31_V31_M33_W31', 1),
    ('_Dot33_M33_by_N33', 1),
    ('_Input_CosD2', 1),
    ('_Input_SinD2', 1),
    ('_Outer33_V31_W31', 1),
    ('_Add33_M33_N33', 2),
    ('_AXPY33_M33_N33', 2),
    ('_Dot33_T333_M33_V3', 1),
    ('_Factor33', 2),
    ('_Dot33_T333_V3', 1),
    ('__StartCheckTEMPmem', 1),
    ('__EndCheckTEMPmem', 1),
]


def engine_dir():
    d = os.environ.get('CFMS_DRAPE_DIR')
    return d if d else os.path.join(os.path.dirname(ROOT), 'cfmsDrape_dev')


def show(repo, path):
    return subprocess.run(['git', '-C', repo, 'show', '%s:%s' % (COMMIT, SRC + path)],
                          check=True, stdout=subprocess.PIPE).stdout


def lines_of(blob):
    return blob.splitlines(keepends=True)


def text(line):
    return line.decode('latin-1').rstrip('\r\n')


def span(lines, start_pat, end_pat, start_from=0):
    """start_pat 와 맞는 첫 줄부터 end_pat 와 맞는 첫 줄까지(포함). 정규식은 latin-1 로 읽은 줄에 건다."""
    a = next(i for i in range(start_from, len(lines)) if re.search(start_pat, text(lines[i])))
    b = next(i for i in range(a, len(lines)) if re.search(end_pat, text(lines[i])))
    return a, b


def func_defs(lines, qual_name, count):
    """0 열 정의 줄들(qual_name 다음 '(')과 각자의 끝 「}」."""
    pat = re.compile(r'^\S.*\b' + re.escape(qual_name) + r'\s*\(')
    out = []
    for i, ln in enumerate(lines):
        if pat.search(text(ln)) and not text(ln).rstrip().endswith(';'):
            j = next(k for k in range(i + 1, len(lines)) if text(lines[k]).rstrip() == '}')
            out.append((i, j))
    if len(out) != count:
        raise SystemExit('정의 개수가 다르다: %s %d (기대 %d)' % (qual_name, len(out), count))
    return out


def piece(name, path, lines, a, b):
    data = b''.join(lines[a:b + 1])
    return {'name': name, 'source': path, 'first_line': a + 1, 'last_line': b + 1,
            'sha256': hashlib.sha256(data).hexdigest(), 'data': data}


def collect(repo):
    core = lines_of(show(repo, 'SDECore.cpp'))
    hdr = lines_of(show(repo, 'SDECore.h'))
    iv = lines_of(show(repo, 'SDEImplicitVariables.h'))
    files = {}

    # 1) 전역: BEST_SSE_* 포인터 선언 · 상수 · my_sse_crossproduct
    g = []
    a, b = span(core, r'^void \(\*BEST_SSE_AX\)\(', r'^dFloat\* \(\*BEST_SSE_MALLOC\)\(')
    g.append(piece('BEST_SSE pointers', 'SDECore.cpp', core, a, b))
    a, b = span(core, r'static dFloat MINUS = -1\., PLUS = 1\.', r'static int _K = 3;')
    g.append(piece('constants', 'SDECore.cpp', core, a, b))
    (a, b), = func_defs(core, 'my_sse_crossproduct', 1)
    g.append(piece('my_sse_crossproduct', 'SDECore.cpp', core, a, b))
    files['SDECore_globals.inc'] = g

    # 2) 실행 때 포인터 대입(_DO_NOT_USE_SSE 경로) — 하네스가 함수 본문 안에 #include 한다
    a, b = span(core, r'^\s*BEST_SSE_AX\s*=\s*my_SSE_AX;', r'^\s*BEST_SSE_MALLOC\s*=\s*my_SSE_MALLOC;')
    files['SDECore_assign.inc'] = [piece('BEST_SSE assignment (__PrepareImplicit)', 'SDECore.cpp', core, a, b)]

    # 3) 선언(기본 인수가 여기 있다) — 하네스가 SDECore 껍데기 클래스 안에 #include 한다
    d = []
    a, b = span(hdr, r'void\s+_Input_b012\(void\);', r'SortBendingTheta0\(void\);')
    d.append(piece('declarations: bending', 'SDECore.h', hdr, a, b))
    a, b = span(hdr, r'void\s+_Tilde_Operation\(', r'void _Input_F_and_J\(int I0, int I\);')
    d.append(piece('declarations: SoA helpers', 'SDECore.h', hdr, a, b))
    a, b = span(hdr, r'__StartCheckTEMPmem\(', r'__EndCheckTEMPmem\(')
    d.append(piece('declarations: TEMP check', 'SDECore.h', hdr, a, b))
    files['SDECore_decl.inc'] = d

    # 4) sbv 클래스 선언과 NUM_TEMP
    v = []
    a, b = span(iv, r'#define NUM_TEMP\b', r'#define NUM_TEMP\b')
    v.append(piece('NUM_TEMP', 'SDEImplicitVariables.h', iv, a, b))
    a, b = span(iv, r'class _SerialBendingVariables \{', r'^\s*\};\s*$')
    p = piece('class _SerialBendingVariables', 'SDEImplicitVariables.h', iv, a, b)
    if b'_Na_outer9' not in p['data'] or b'BendVars' in p['data']:
        raise SystemExit('_SerialBendingVariables 범위가 예상과 다르다')
    v.append(p)
    files['SDEImplicitVariables_sbv.inc'] = v

    # 5) 멤버 함수 정의
    defs = []
    for name, count in MEMBER_FUNCS:
        for k, (a, b) in enumerate(func_defs(core, 'SDECore::' + name, count)):
            defs.append(piece('SDECore::%s%s' % (name, '' if count == 1 else ' #%d' % (k + 1)), 'SDECore.cpp', core, a, b))
    defs.sort(key=lambda q: q['first_line'])
    files['SDECore_macri.inc'] = defs

    # 6) Bridson(Oh) 경로 — docs/GATE_BRIDSON.md: 전역 OuterProduct, s_BendSkip 선언, GetImpBending_Oh_SOA
    oh = []
    a, b = span(core, r'^void\s+OuterProduct\(dFloat\*a, dFloat\*b, dFloat \*outer\)', r'^}\s*$')   # 0 열 주석에도 이름이 있다
    oh.append(piece('OuterProduct', 'SDECore.cpp', core, a, b))
    a, b = span(core, r'^static std::vector<char> s_BendSkip;', r'^static std::vector<char> s_BendSkip;')
    oh.append(piece('s_BendSkip', 'SDECore.cpp', core, a, b))
    (a, b), = func_defs(core, 'SDECore::GetImpBending_Oh_SOA', 1)
    p = piece('SDECore::GetImpBending_Oh_SOA', 'SDECore.cpp', core, a, b)
    if b'pSerialJacobianJ[ pCGSolver->piBendingJ[i*6+4]]' not in p['data']:
        raise SystemExit('GetImpBending_Oh_SOA 범위가 예상과 다르다(직렬 환원이 없다)')
    oh.append(p)
    files['SDECore_oh.inc'] = oh

    # 7) SDEMatrix3f 의 두 정의(Oh 경로가 부르는 오버로드만) — 하네스의 SDEMatrix3f 껍데기가 선언한다
    mat = lines_of(show(repo, 'SDEMatrix3f.cpp'))
    m3 = []
    a, b = span(mat, r'^void\s+SDEMatrix3f::AXPY\( double a, double x\[3\]\[3\], double b\)', r'^}\s*$')
    m3.append(piece('SDEMatrix3f::AXPY(double, double[3][3], double)', 'SDEMatrix3f.cpp', mat, a, b))
    a, b = span(mat, r'^void\s+SDEMatrix3f::Set\( double\* d,double alpha\)', r'^}\s*$')
    m3.append(piece('SDEMatrix3f::Set(double*, double)', 'SDEMatrix3f.cpp', mat, a, b))
    files['SDEMatrix3f_oh.inc'] = m3

    # 8) my_SSE_* 파일 전체
    for fn in ('SDEmysse.h', 'SDEmysse.cpp'):
        blob = show(repo, fn)
        files[fn] = [{'name': fn + ' (whole file)', 'source': fn, 'first_line': 1, 'last_line': len(lines_of(blob)),
                      'sha256': hashlib.sha256(blob).hexdigest(), 'data': blob}]
    return files


def write(files):
    os.makedirs(OUT, exist_ok=True)
    for fn, pieces in files.items():
        with open(os.path.join(OUT, fn), 'wb') as f:
            if fn.startswith('SDEmysse'):
                f.write(pieces[0]['data'])      # 파일 전체 — 구분 주석도 넣지 않는다
                continue
            for p in pieces:
                f.write(('// ==== %s — %s:%d-%d (cfmsDrape_dev %s, 바이트 그대로) ====\n'
                         % (p['name'], p['source'], p['first_line'], p['last_line'], COMMIT)).encode('utf-8'))
                f.write(p['data'])
                if not p['data'].endswith(b'\n'):
                    f.write(b'\n')


def main():
    pin = '--pin' in sys.argv[1:]
    repo = engine_dir()
    if not os.path.isdir(os.path.join(repo, '.git')):
        print('엔진 리포 없음(%s) — 엔진 손 코드 비교 게이트를 건너뛴다' % repo)
        return 3
    files = collect(repo)
    write(files)
    table = {fn: [{k: p[k] for k in ('name', 'source', 'first_line', 'last_line', 'sha256')} for p in ps]
             for fn, ps in files.items()}
    n = sum(len(ps) for ps in table.values())
    if pin:
        with open(MANIFEST, 'w', encoding='utf-8', newline='\n') as f:
            json.dump({'engine_commit': COMMIT, 'files': table}, f, ensure_ascii=False, indent=1)
            f.write('\n')
        print('E1 고정값을 썼다: 조각 %d 개 → %s' % (n, os.path.relpath(MANIFEST, ROOT)))
        return 0
    if not os.path.exists(MANIFEST):
        print('E1 FAIL: 고정값 파일이 없다 — 처음이면 --pin')
        return 1
    ref = json.load(open(MANIFEST, encoding='utf-8'))
    bad = []
    if ref.get('engine_commit') != COMMIT:
        bad.append('커밋 %s ≠ 고정값 %s' % (COMMIT, ref.get('engine_commit')))
    for fn in sorted(set(ref['files']) | set(table)):
        a = [(p['name'], p['sha256']) for p in ref['files'].get(fn, [])]
        b = [(p['name'], p['sha256']) for p in table.get(fn, [])]
        if a != b:
            bad.append(fn)
    if bad:
        print('E1 FAIL: 고정값과 다르다 — %s' % ', '.join(bad))
        return 1
    print('E1 PASS: 엔진 %s 조각 %d 개(파일 %d 개)가 고정값 SHA-256 과 같다 → generated/engine_ref/'
          % (COMMIT, n, len(table)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
