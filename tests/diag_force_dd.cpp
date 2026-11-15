// diag_force_dd.cpp — 진단(판정 밖, docs/GATE_ASSEMBLY.md 2 차 결과). R1′ f+d 차이가 누구의 반올림인가.
// 같은 C · g · v 에서 double-double 로 f = −(kC + kd g·v) g 를 만들어 기준(조밀 조립) · 1 차 AssembleImplicit ·
// 2 차 EvalImplicit 를 각각 견준다. 빌드: cl ... tests\diag_force_dd.cpp third_party\sxmesh\SXMesh.cpp (리포 루트에서 실행)
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffBend.h"
#include "gate_common.h"
#include <cmath>
using namespace gate;
typedef sxgen::BaraffBend Spec;
static const int N = Spec::kDof, NH = Spec::kHess;
static const double kK = 1.0, kKd = 0.1, kH = 0.01;

struct DD { double hi, lo; };
static DD TwoSum(double a, double b) { const double s = a + b, bb = s - a; return { s, (a - (s - bb)) + (b - bb) }; }
static DD TwoProd(double a, double b) { const double p = a * b; return { p, std::fma(a, b, -p) }; }
static DD Add(DD a, DD b) { DD s = TwoSum(a.hi, b.hi); s.lo += a.lo + b.lo; return TwoSum(s.hi, s.lo); }
static DD Mul(DD a, double b) { DD p = TwoProd(a.hi, b); p.lo += a.lo * b; return TwoSum(p.hi, p.lo); }

int main()
{
	const std::vector<Hinge> hs = ReadHinges("results/hinges.txt");
	const double prm[3] = { kK, kKd, kH };
	Stat eRef, eFirst, eImp, cond;
	double worst = 0.; std::size_t ws = 0;
	for (std::size_t s = 0; s < hs.size(); ++s) {
		const Hinge& hg = hs[s];
		double C, g[N], H[NH];
		Spec::EvalStaged(hg.x, &hg.theta0, &C, g, H);
		sx::ConditionForce<N> dense;
		sx::AssembleConditionForce<N>(1, &C, g, H, hg.v, kK, kKd, false, dense);
		sx::ImplicitBlock<N> blk;
		sx::AssembleImplicit<N>(1, &C, g, H, hg.v, kK, kKd, kH, blk);
		double Ci, fi[N], Ai[NH], bi[N];
		Spec::EvalImplicit(hg.x, &hg.theta0, hg.v, prm, &Ci, fi, Ai, bi);

		DD cd = { 0., 0. };
		for (int a = 0; a < N; ++a) cd = Add(cd, TwoProd(g[a], hg.v[a]));
		const DD sdd = Add(TwoProd(kK, C), Mul(cd, kKd));
		double fdd[N], fref[N], mx = 0.;
		for (int a = 0; a < N; ++a) {
			const DD f = Mul(sdd, -g[a]);
			fdd[a] = f.hi + f.lo;
			fref[a] = dense.f[a] + dense.d[a];
			mx = std::max(mx, std::fabs(fdd[a]));
		}
		double mEl = 0.;   // 탄성력 크기 / 전체 힘 크기 = 상쇄 정도
		for (int a = 0; a < N; ++a) mEl = std::max(mEl, std::fabs(kK * C * g[a]));
		cond.Add(mEl / mx);
		const double r = MaxAbsDiff(fref, fdd, N) / mx;
		eRef.Add(r);
		eFirst.Add(MaxAbsDiff(blk.f, fdd, N) / mx);
		const double ri = MaxAbsDiff(fi, fdd, N) / mx;
		eImp.Add(ri);
		if (ri > worst) { worst = ri; ws = s; }
	}
	std::printf("f+d 대 double-double (같은 C · g · v), 1000 힌지, 상대오차 max / 중앙\n");
	std::printf("  기준(조밀 조립 f + d)    %.2e / %.2e\n", eRef.Max(), eRef.Median());
	std::printf("  1 차 AssembleImplicit    %.2e / %.2e\n", eFirst.Max(), eFirst.Median());
	std::printf("  2 차 EvalImplicit        %.2e / %.2e\n", eImp.Max(), eImp.Median());
	std::printf("  상쇄 max|kC g| / max|f+d|: max %.2e, 중앙 %.2e;  EvalImplicit 최악 힌지 %zu\n", cond.Max(), cond.Median(), ws);
	return 0;
}
