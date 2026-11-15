//---------------------------------------------------------------------------
// diag_poisson_b.cpp — Poisson 힘 게이트 1 차의 Q6 b 미달(1.26e-12, 기준 1e-12)을 들여다본다. 판정 밖 진단(docs/GATE_POISSON.md).
// 성분이 하나라 Gauss–Newton 의 b 는 스칼라 하나에 ∇C 를 곱한 것이다:
//   b_GN = −h s ∇C,   s = k C + (kd + h k) Ċ
// k C 와 (kd + h k) Ċ 가 서로 지우는 표본에서는 |b| 가 작아져 상대 오차가 커진다. 그것이 까닭인지를 표본마다
//   상대 오차 = max|Δb| / max|b_기준|,   지움 = |s| / S,   지움 없는 척도 = max|Δb| / (h S max|∇C|),
//   S = k|C| + (kd + h k) Σ_a |∇_a C v_a|     (Ċ = Σ_a ∇_a C v_a 안에서 지우는 몫도 뺀다)
// 로 본다. 정확판(Q4)도 같은 식으로 본다 — b_E = −h s ∇C − h² (k C + kd Ċ) ∇²C v, 척도에 h² (k|C| + kd Σ|∇C v|) max_a Σ_b |∇²C_ab v_b| 를 더한다.
// 실행 디렉터리 = 리포 루트. results/diag_poisson_b.json 을 쓴다.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/SulPoisson.h"
#include "gate_common.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;
typedef sxgen::SulPoisson S;

enum { N = S::kDof, NC = S::kComp, NH = S::kHess };

static const double kK = 1.0, kKd = 0.1, kH = 0.01;
static const double kPrm[3] = { kK, kKd, kH };

struct Row { double rel, cancel, robust; int idx; };

static void Summary(const char* name, const std::vector<Row>& rows, std::FILE* js, const char* tail)
{
	Row worst = rows[0];
	double robustMax = 0., prodMax = 0.;
	int over13 = 0, over12 = 0;
	Stat rel, robust;
	for (std::size_t i = 0; i < rows.size(); ++i) {
		const Row& r = rows[i];
		if (r.rel > worst.rel) worst = r;
		robustMax = std::max(robustMax, r.robust);
		prodMax = std::max(prodMax, r.rel * r.cancel);
		if (r.rel > 1e-13) ++over13;
		if (r.rel > 1e-12) ++over12;
		rel.Add(r.rel); robust.Add(r.robust);
	}
	std::printf("%s\n", name);
	std::printf("  상대 오차 최대 %.3e (표본 %d, 그 표본의 지움 %.3e · 지움 없는 척도 %.3e), 중앙 %.3e\n",
		worst.rel, worst.idx, worst.cancel, worst.robust, rel.Median());
	std::printf("  상대 오차 > 1e-13 %d 개 · > 1e-12 %d 개\n", over13, over12);
	std::printf("  지움 없는 척도 최대 %.3e · 중앙 %.3e,  상대 오차 × 지움 의 최대 %.3e\n", robustMax, robust.Median(), prodMax);
	std::fprintf(js, "  \"%s\": {\"rel_max\": %.3e, \"rel_median\": %.3e, \"worst_sample\": %d, \"worst_cancel\": %.3e, \"worst_robust\": %.3e, "
		"\"over_1e-13\": %d, \"over_1e-12\": %d, \"robust_max\": %.3e, \"robust_median\": %.3e, \"rel_times_cancel_max\": %.3e}%s\n",
		name, worst.rel, rel.Median(), worst.idx, worst.cancel, worst.robust, over13, over12, robustMax, robust.Median(), prodMax, tail);
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Tri> tris = ReadTriangles("results/triangles.txt");
	const std::vector<std::vector<double> > beta = ReadRows("results/poisson_beta.txt");
	if (tris.empty() || beta.size() < tris.size()) { std::printf("표본 없음\n"); return 2; }
	std::vector<Row> gn, ex;
	for (std::size_t s = 0; s < tris.size(); ++s) {
		double ev[7];
		std::memcpy(ev, tris[s].uv, 6 * sizeof(double));
		ev[6] = beta[s][0];
		const double* x = tris[s].x;
		const double* v = tris[s].v;
		double C[NC], g[NC * N], H[NC * NH], H0[NC * NH];
		S::Eval(x, ev, C, g, H);
		for (int i = 0; i < NC * NH; ++i) H0[i] = 0.;
		double Cd = 0., CdAbs = 0.;   // Ċ = ∇C·v 와 그 항들의 크기 합(Ċ 안의 지움을 뺀 척도)
		for (int a = 0; a < N; ++a) { Cd += g[a] * v[a]; CdAbs += std::fabs(g[a] * v[a]); }
		double Hv[N];
		for (int a = 0; a < N; ++a) Hv[a] = 0.;
		for (int a = 0; a < N; ++a)
			for (int b = 0; b < N; ++b) Hv[a] += std::fabs(H[a <= b ? sx::UpperIndex<N>(a, b) : sx::UpperIndex<N>(b, a)] * v[b]);
		const double sAbs = kK * std::fabs(C[0]) + (kKd + kH * kK) * CdAbs;
		const double sVal = kK * C[0] + (kKd + kH * kK) * Cd;
		const double gmax = MaxAbs(g, N);

		sx::ImplicitBlock<N> refG, refE;
		sx::AssembleImplicit<N>(NC, C, g, H0, v, kK, kKd, kH, refG);
		sx::AssembleImplicit<N>(NC, C, g, H, v, kK, kKd, kH, refE);
		double Co[NC], fo[N], Ao[NH], bo[N];
		S::EvalImplicitGN(x, ev, v, kPrm, Co, fo, Ao, bo);
		Row r;
		r.idx = static_cast<int>(s);
		r.rel = Rel(MaxAbsDiff(bo, refG.b, N), MaxAbs(refG.b, N));
		r.cancel = sAbs > 0. ? std::fabs(sVal) / sAbs : 1.;
		r.robust = MaxAbsDiff(bo, refG.b, N) / (kH * sAbs * gmax);
		gn.push_back(r);
		S::EvalImplicit(x, ev, v, kPrm, Co, fo, Ao, bo);
		const double scaleE = kH * sAbs * gmax + kH * kH * (kK * std::fabs(C[0]) + kKd * CdAbs) * MaxAbs(Hv, N);
		r.rel = Rel(MaxAbsDiff(bo, refE.b, N), MaxAbs(refE.b, N));
		r.robust = MaxAbsDiff(bo, refE.b, N) / scaleE;
		ex.push_back(r);
	}
	std::FILE* js = std::fopen("results/diag_poisson_b.json", "w");
	if (!js) return 2;
	std::fprintf(js, "{\n  \"samples\": %zu, \"k\": %g, \"kd\": %g, \"h\": %g,\n", tris.size(), kK, kKd, kH);
	std::printf("Poisson 힘 게이트 판정 뒤 진단 — b 의 상대 오차와 지움 (%zu 표본)\n", tris.size());
	Summary("gn_Q6", gn, js, ",");
	Summary("exact_Q4", ex, js, "");
	std::fprintf(js, "}\n");
	std::fclose(js);
	return 0;
}
