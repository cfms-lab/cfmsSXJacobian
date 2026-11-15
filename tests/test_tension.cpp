//---------------------------------------------------------------------------
// test_tension.cpp — 인장 C > 0 판 게이트의 정의 쪽(docs/GATE_TENSION.md).
// P1:  EvalImplicitTension 이 기준 조립(Eval 의 C · ∇C · ∇²C → AssembleImplicit(∇²C = 0) + h² k Σ max(C, 0) ∇²C)과 같은지,
//      표본에 C > 0 · C ≤ 0 성분이 모두 있는지.
// P1b: 두 성분이 모두 늘어난 표본(v = 0)에서 EvalImplicit 와, 모두 눌린 표본에서 EvalImplicitGN 과 같은지.
// P2:  A_T 가 모든 표본에서 양의 준정부호인지.
// P0(기존 변형 불변)은 tests/check_tension_p0.py, 적분 쪽은 test_integration.exe tension.
// 실행 디렉터리 = 리포 루트. results/triangles.txt 를 읽고 results/gate_tension.json 을 쓴다.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffStretch.h"
#include "gate_common.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;
typedef sxgen::BaraffStretch S;

static const double kPrmK = 1.0, kPrmKd = 0.1, kPrmH = 0.01;   // Gauss–Newton 게이트와 같다
static const double kPrm[3] = { kPrmK, kPrmKd, kPrmH };
static const double kLimit = 1e-12;     // P1 · P1b
static const double kPsd = -1e-12;      // P2: 최소 고윳값 / 최대 고윳값 크기
static const int kMinBranch = 100;      // P1 갈래

// 대칭 행렬(n × n 전체, 행 우선)의 고윳값 — 순환 Jacobi(Gauss–Newton 게이트와 같은 방법). M 은 망가진다
static void SymEig(double* M, int n, double* lam)
{
	for (int sweep = 0; sweep < 64; ++sweep) {
		double off = 0., all = 0.;
		for (int p = 0; p < n; ++p)
			for (int q = 0; q < n; ++q) { all += M[p * n + q] * M[p * n + q]; if (q != p) off += M[p * n + q] * M[p * n + q]; }
		if (off <= 1e-32 * all) break;
		for (int p = 0; p < n; ++p)
			for (int q = p + 1; q < n; ++q) {
				const double apq = M[p * n + q];
				if (apq == 0.) continue;
				const double th = (M[q * n + q] - M[p * n + p]) / (2. * apq);
				const double t = (th >= 0. ? 1. : -1.) / (std::fabs(th) + std::sqrt(th * th + 1.));
				const double c = 1. / std::sqrt(t * t + 1.), s = t * c;
				for (int k = 0; k < n; ++k) {
					const double mp = M[k * n + p], mq = M[k * n + q];
					M[k * n + p] = c * mp - s * mq; M[k * n + q] = s * mp + c * mq;
				}
				for (int k = 0; k < n; ++k) {
					const double mp = M[p * n + k], mq = M[q * n + k];
					M[p * n + k] = c * mp - s * mq; M[q * n + k] = s * mp + c * mq;
				}
			}
	}
	for (int p = 0; p < n; ++p) lam[p] = M[p * n + p];
}

// 위 삼각 저장 → 최소 고윳값 / 최대 고윳값 크기
static double MinEigRatioUpper(const double* U)
{
	enum { N = S::kDof };
	double M[N * N], lam[N];
	for (int a = 0; a < N; ++a)
		for (int b = 0; b < N; ++b) M[a * N + b] = U[a <= b ? sx::UpperIndex<N>(a, b) : sx::UpperIndex<N>(b, a)];
	SymEig(M, N, lam);
	double lo = lam[0], big = 0.;
	for (int p = 0; p < N; ++p) { lo = std::min(lo, lam[p]); big = std::max(big, std::fabs(lam[p])); }
	return big > 0. ? lo / big : 0.;
}

// 기준 조립: AssembleImplicit(∇²C = 0) 에 h² k Σ C⁺ ∇²C 를 A 에 더하고 b 에서 h² k Σ C⁺ ∇²C v 를 뺀다
static void Reference(const double* C, const double* g, const double* H, const double* v, sx::ImplicitBlock<S::kDof>& ref)
{
	enum { N = S::kDof, NH = S::kHess, NC = S::kComp };
	double H0[NC * NH];
	for (int i = 0; i < NC * NH; ++i) H0[i] = 0.;
	sx::AssembleImplicit<N>(NC, C, g, H0, v, kPrmK, kPrmKd, kPrmH, ref);
	const double hhk = kPrmH * kPrmH * kPrmK;
	for (int c = 0; c < NC; ++c) {
		const double cp = C[c] > 0. ? C[c] : 0.;
		const double* Hc = H + c * NH;
		double Hv[N];
		for (int a = 0; a < N; ++a) Hv[a] = 0.;
		int i = 0;
		for (int a = 0; a < N; ++a) {
			Hv[a] += Hc[i] * v[a];
			ref.A[i] += hhk * cp * Hc[i];
			++i;
			for (int j = a + 1; j < N; ++j, ++i) {
				Hv[a] += Hc[i] * v[j];
				Hv[j] += Hc[i] * v[a];
				ref.A[i] += hhk * cp * Hc[i];
			}
		}
		for (int a = 0; a < N; ++a) ref.b[a] -= hhk * cp * Hv[a];
	}
}

struct Cmp { Stat c, A, b, f; int n; Cmp() : n(0) {} };

static void Compare(const double* C0, const double* Cs, double cscale, const double* A0, const double* A1, const double* b0,
	const double* b1, const double* f0, const double* f1, Cmp& m)
{
	enum { N = S::kDof, NH = S::kHess, NC = S::kComp };
	m.c.Add(MaxAbsDiff(Cs, C0, NC) / cscale);
	m.A.Add(Rel(MaxAbsDiff(A1, A0, NH), MaxAbs(A0, NH)));
	m.b.Add(Rel(MaxAbsDiff(b1, b0, N), MaxAbs(b0, N)));
	m.f.Add(Rel(MaxAbsDiff(f1, f0, N), MaxAbs(f0, N)));
	++m.n;
}

static void GateCmp(Report& rpt, const char* title, const Cmp& m, bool withC)
{
	std::printf(" %s (%d 개)\n", title, m.n);
	if (withC) rpt.Gate("C", m.n > 0 && m.c.Max() <= kLimit, "|ΔC|/√a %.2e", m.c.Max(), kLimit);
	rpt.Gate("A", m.n > 0 && m.A.Max() <= kLimit, "rel max %.2e", m.A.Max(), kLimit);
	rpt.Gate("b", m.n > 0 && m.b.Max() <= kLimit, "rel max %.2e", m.b.Max(), kLimit);
	rpt.Gate("f+d", m.n > 0 && m.f.Max() <= kLimit, "rel max %.2e", m.f.Max(), kLimit);
}

static void JsonCmp(std::FILE* js, const char* name, const Cmp& m, const char* tail)
{
	std::fprintf(js, "    \"%s\": {\"samples\": %d, \"C_abs_over_sqrt_a\": %.3e, \"A_rel\": %.3e, \"b_rel\": %.3e, \"f_rel\": %.3e, "
		"\"A_rel_median\": %.3e, \"b_rel_median\": %.3e, \"f_rel_median\": %.3e}%s\n",
		name, m.n, m.c.Max(), m.A.Max(), m.b.Max(), m.f.Max(), m.A.Median(), m.b.Median(), m.f.Median(), tail);
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	enum { N = S::kDof, NH = S::kHess, NC = S::kComp };
	const std::vector<Tri> ts = ReadTriangles("results/triangles.txt");
	if (ts.empty()) { std::printf("results/triangles.txt 없음\n"); return 2; }
	Report rpt("results/gate_tension.json");
	std::FILE* js = rpt.json;
	std::printf("인장 C > 0 판 게이트 (정의) — 삼각형 %zu 개, k = %g, kd = %g, h = %g\n", ts.size(), kPrmK, kPrmKd, kPrmH);

	Cmp p1, pPos, pNeg;
	int nPosComp = 0, nNegComp = 0, nAllPos = 0, nAllNeg = 0, nIndefE = 0, nIndefT = 0;
	double tMin = 1e30, gnMin = 1e30;
	Stat diffTE, eigE;
	const double zero[N] = { 0. };
	for (std::size_t s = 0; s < ts.size(); ++s) {
		const Tri& t = ts[s];
		const double sa = std::sqrt(UvArea(t.uv));
		double C[NC], g[NC * N], H[NC * NH];
		S::Eval(t.x, t.uv, C, g, H);
		sx::ImplicitBlock<N> ref;
		Reference(C, g, H, t.v, ref);
		double Ct[NC], ft[N], At[NH], bt[N];
		S::EvalImplicitTension(t.x, t.uv, t.v, kPrm, Ct, ft, At, bt);
		Compare(C, Ct, sa, ref.A, At, ref.b, bt, ref.f, ft, p1);
		bool allPos = true, allNeg = true;
		for (int c = 0; c < NC; ++c) {
			if (C[c] > 0.) { ++nPosComp; allNeg = false; } else { ++nNegComp; allPos = false; }
		}
		if (allPos) {   // P1b (i): v = 0 이면 s = kC 라 정확판과 같다
			++nAllPos;
			double C0[NC], f0[N], A0[NH], b0[N], C1[NC], f1[N], A1[NH], b1[N];
			S::EvalImplicit(t.x, t.uv, zero, kPrm, C0, f0, A0, b0);
			S::EvalImplicitTension(t.x, t.uv, zero, kPrm, C1, f1, A1, b1);
			Compare(C0, C1, sa, A0, A1, b0, b1, f0, f1, pPos);
		}
		if (allNeg) {   // P1b (ii): 눌린 성분만이면 GN 판과 같다
			++nAllNeg;
			double C0[NC], f0[N], A0[NH], b0[N];
			S::EvalImplicitGN(t.x, t.uv, t.v, kPrm, C0, f0, A0, b0);
			Compare(C0, Ct, sa, A0, At, b0, bt, f0, ft, pNeg);
		}
		// P2 와 보고
		const double rT = MinEigRatioUpper(At);
		tMin = std::min(tMin, rT);
		if (rT < kPsd) ++nIndefT;
		double Ce[NC], fe[N], Ae[NH], be[N], Cg[NC], fg[N], Ag[NH], bg[N];
		S::EvalImplicit(t.x, t.uv, t.v, kPrm, Ce, fe, Ae, be);
		S::EvalImplicitGN(t.x, t.uv, t.v, kPrm, Cg, fg, Ag, bg);
		const double rE = MinEigRatioUpper(Ae);
		if (rE < kPsd) { ++nIndefE; eigE.Add(-rE); }
		gnMin = std::min(gnMin, MinEigRatioUpper(Ag));
		diffTE.Add(Rel(MaxAbsDiff(At, Ae, NH), MaxAbs(Ae, NH)));
	}

	GateCmp(rpt, "P1 EvalImplicitTension 대 기준 조립(∇²C = 0 + h² k Σ C⁺ ∇²C)", p1, true);
	std::printf(" P1 갈래 — 인장 성분 C > 0 %d 개 · C ≤ 0 %d 개\n", nPosComp, nNegComp);
	rpt.Gate("C>0", nPosComp >= kMinBranch, "%.0f 개", nPosComp, kMinBranch);
	rpt.Gate("C≤0", nNegComp >= kMinBranch, "%.0f 개", nNegComp, kMinBranch);
	GateCmp(rpt, "P1b (i) 두 성분 모두 C > 0 · v = 0 — EvalImplicit 와", pPos, false);
	GateCmp(rpt, "P1b (ii) 두 성분 모두 C ≤ 0 — EvalImplicitGN 과", pNeg, false);
	std::printf(" P2 준정부호 (%zu 개)\n", ts.size());
	rpt.Gate("A_T", tMin >= kPsd, "최소 고윳값 / 최대 크기의 최솟값 %.2e", tMin, kPsd);
	std::printf("  보고 — 정확판 A 부정치 %d 개(|최소|/최대 중앙 %.2e), A_T 부정치 %d 개, GN A 최솟값 %.2e, |A_T − A_E|/|A_E| 중앙 %.2e · max %.2e\n",
		nIndefE, eigE.Median(), nIndefT, gnMin, diffTE.Median(), diffTE.Max());

	std::fprintf(js, "{\n  \"triangles\": %zu,\n  \"k\": %g, \"kd\": %g, \"h\": %g,\n", ts.size(), kPrmK, kPrmKd, kPrmH);
	std::fprintf(js, "  \"P1\": {\n");
	JsonCmp(js, "tension_vs_reference", p1, ",");
	std::fprintf(js, "    \"components_C_pos\": %d, \"components_C_nonpos\": %d, \"min_branch\": %d\n  },\n", nPosComp, nNegComp, kMinBranch);
	std::fprintf(js, "  \"P1b\": {\n");
	JsonCmp(js, "all_pos_v0_vs_implicit", pPos, ",");
	JsonCmp(js, "all_nonpos_vs_gn", pNeg, "");
	std::fprintf(js, "  },\n");
	std::fprintf(js, "  \"P2\": {\"min_eig_ratio\": %.3e, \"limit\": %g, \"indefinite\": %d},\n", tMin, kPsd, nIndefT);
	std::fprintf(js, "  \"report\": {\"indefinite_A_exact\": %d, \"indefinite_A_exact_eig_ratio_median\": %.3e, \"gn_min_eig_ratio\": %.3e, "
		"\"A_T_vs_A_E_rel_median\": %.3e, \"A_T_vs_A_E_rel_max\": %.3e},\n", nIndefE, eigE.Median(), gnMin, diffTE.Median(), diffTE.Max());
	std::fprintf(js, "  \"gates_failed\": %d\n}\n", rpt.fail);
	std::printf("\n인장 C > 0 판 정의 %s (실패 %d) → results/gate_tension.json\n", rpt.fail == 0 ? "PASS" : "FAIL", rpt.fail);
	return rpt.fail == 0 ? 0 : 1;
}
