//---------------------------------------------------------------------------
// link_l0_ref.cpp — 엔진 연결 L0(docs/GATE_LINK_L0.md)의 이 리포 기준: 스텐실마다 평탄판 Eval 의 C · ∇C · ∇²C.
// 사용: link_l0_ref.exe stretch|bridson <입력> <출력>
//   입력 한 줄 = x(9) uv(6) (stretch, BaraffStretch) 또는 x(12) w̄(1) (bridson, BridsonBend).
//   출력 한 줄 = C(kComp) g(kComp·kDof) H(kComp·kHess, ∇²C 의 위 삼각 행 우선), %.17g.
// 조립(전역 힘 · J)은 tests/link_l0.py 가 한다.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffStretch.h"
#include "../generated/BridsonBend.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

template <class S>
static int Run(const char* inPath, const char* outPath)
{
	enum { N = S::kDof, NE = S::kElemValues, NC = S::kComp, NH = S::kHess };
	std::ifstream in(inPath);
	std::FILE* out = std::fopen(outPath, "w");
	if (!in || !out) { std::printf("파일을 열 수 없음\n"); return 2; }
	std::string line;
	int n = 0;
	while (std::getline(in, line)) {
		if (line.empty()) continue;
		std::istringstream is(line);
		double x[N], ev[NE];
		for (int i = 0; i < N; ++i) is >> x[i];
		for (int i = 0; i < NE; ++i) is >> ev[i];
		if (!is) { std::printf("%d 번째 줄을 읽을 수 없음\n", n + 1); std::fclose(out); return 2; }
		double C[NC], g[NC * N], H[NC * NH];
		S::Eval(x, ev, C, g, H);
		for (int i = 0; i < NC; ++i) std::fprintf(out, "%.17g ", C[i]);
		for (int i = 0; i < NC * N; ++i) std::fprintf(out, "%.17g ", g[i]);
		for (int i = 0; i < NC * NH; ++i) std::fprintf(out, "%.17g%s", H[i], i + 1 < NC * NH ? " " : "\n");
		++n;
	}
	std::fclose(out);
	std::printf("link_l0_ref: %d 개\n", n);
	return 0;
}

int main(int argc, char** argv)
{
	if (argc != 4) { std::printf("사용: link_l0_ref.exe stretch|bridson <입력> <출력>\n"); return 2; }
	if (!std::strcmp(argv[1], "stretch")) return Run<sxgen::BaraffStretch>(argv[2], argv[3]);
	if (!std::strcmp(argv[1], "bridson")) return Run<sxgen::BridsonBend>(argv[2], argv[3]);
	std::printf("명세는 stretch 또는 bridson\n");
	return 2;
}
