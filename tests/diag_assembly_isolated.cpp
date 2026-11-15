//---------------------------------------------------------------------------
// diag_assembly_isolated.cpp — 조립 게이트의 진단(판정 밖, docs/GATE_ASSEMBLY.md 「결과」).
// 평가를 빼고 조립만 잰다. 격자 힌지 4096 개의 C · g · H · v 를 미리 계산해 두고 조립만 30 회 반복, 7 회 중 최소.
// 빌드(리포 루트, vcvars64 뒤):
//   cl /nologo /EHsc /O2 /MD /DNDEBUG /utf-8 /std:c++14 /bigobj tests\diag_assembly_isolated.cpp third_party\sxmesh\SXMesh.cpp
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffBend.h"
#include "gate_common.h"
#include <chrono>
#include <cstdio>
using namespace gate;
typedef sxgen::BaraffBend Spec;
static const int N = 12, NH = 78;
struct Pre { double C, g[N], H[NH], v[N]; int nd[4]; };
int main() {
	std::mt19937_64 rng(7);
	const Grid grid = MakeGrid(200, rng);
	const int M = 4096;
	std::vector<Pre> P(M);
	for (int s = 0; s < M; ++s) {
		double x[N]; const double t0 = 0.;
		for (int a = 0; a < 4; ++a) { P[s].nd[a] = grid.hinges[4*s+a]; for (int c = 0; c < 3; ++c) { x[3*a+c] = grid.X[3*P[s].nd[a]+c]; P[s].v[3*a+c] = grid.V[3*P[s].nd[a]+c]; } }
		Spec::EvalStaged(x, &t0, &P[s].C, P[s].g, P[s].H);
	}
	std::vector<double> Fn(3*grid.nNode), Bn(3*grid.nNode);
	typedef std::chrono::steady_clock Clock;
	double best[3] = {1e30,1e30,1e30}, cs = 0.;
	const int R = 30;
	for (int rep = 0; rep < 7; ++rep) {
		for (int mode = 0; mode < 3; ++mode) {
			const Clock::time_point t0 = Clock::now();
			for (int r = 0; r < R; ++r) for (int s = 0; s < M; ++s) {
				const Pre& p = P[s];
				if (mode == 0) { cs += p.H[0] + p.g[0]; }
				else if (mode == 1) {
					sx::ConditionForce<N> o; sx::AssembleConditionForce<N>(1, &p.C, p.g, p.H, p.v, 1., 0.1, false, o);
					for (int a = 0; a < 4; ++a) for (int c = 0; c < 3; ++c) Fn[3*p.nd[a]+c] += o.f[3*a+c] + o.d[3*a+c];
					cs += o.K[13] + o.Dx[7] + o.Dv[5];
				} else {
					sx::ImplicitBlock<N> o; sx::AssembleImplicit<N>(1, &p.C, p.g, p.H, p.v, 1., 0.1, 0.01, o);
					for (int a = 0; a < 4; ++a) for (int c = 0; c < 3; ++c) { Fn[3*p.nd[a]+c] += o.f[3*a+c]; Bn[3*p.nd[a]+c] += o.b[3*a+c]; }
					cs += o.A[13];
				}
			}
			best[mode] = std::min(best[mode], std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / (double(R) * M));
		}
	}
	std::printf("isolated ns/hinge: read-only %.1f, old assembly %.1f, new assembly %.1f -> (new-base)/(old-base) = %.3f  cs %.3e\n",
		best[0], best[1], best[2], (best[2]-best[0])/(best[1]-best[0]), cs);
}
