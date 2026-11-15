//---------------------------------------------------------------------------
// test_poisson.cpp — Poisson 힘 게이트 Q1 – Q7 (docs/GATE_POISSON.md) 와 보고용 측정.
// Provot 전단 게이트(test_provot_shear.cpp)의 V1 – V5 와 Gauss–Newton 게이트(test_gauss_newton.cpp)의 N1 을 명세 SulPoisson 에 옮긴 것이다.
// 실행 디렉터리 = 리포 루트. results/triangles.txt · results/poisson_beta.txt · results/mma_ref_poisson.txt 를 읽고
// results/gate_poisson.json 을 쓴다. Q0(기존 출력 불변)은 check_flat_unchanged.py · check_tension_p0.py · check_link_l1_p0.py · check_poisson_p0.py.
// 2 차(판정문 「2 차 기준」): Q4 · Q6 의 b 는 상대 오차가 아니라 지움을 뺀 척도로 판정한다 — 성분이 하나라 b 가 스칼라 × ∇C 이고,
// b 가 0 에 가까운 표본에서는 상대 오차가 지움의 크기를 잰다(tests/diag_poisson_b.cpp). 1 차의 상대 오차는 보고로 남긴다.
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

enum { N = S::kDof, NC = S::kComp, NH = S::kHess, NE = S::kElemValues };
static_assert(S::kNodes == 3 && N == 9 && NC == 1 && NE == 7 && S::kElemVars == 2, "삼각형 · 성분 1 · uv 6 칸 + β 1 칸");

static const double kK = 1.0, kKd = 0.1, kH = 0.01;   // 판정문에 고정한 계수(조립 · Provot 전단 게이트와 같음)
static const double kPrm[3] = { kK, kKd, kH };
static const int kMmaCols = 1 + 9 + 45 + 9 + 45;     // C · g · H · 원문 닫힌 꼴 g · H
static const double kLimit = 1e-12;                  // Q3 · Q4 · Q6
static const double kPsd = -1e-12;                   // Q7: 최소 고윳값 / 최대 고윳값 크기
static const int kMinBranch = 100;                   // Q6 갈래

struct Sample { double x[9]; double ev[7]; double v[9]; };   // ev = uv(6) · β

static double CScale(const double* ev) { return std::sqrt(UvArea(ev)); }   // C 오차의 척도 √a

static double Len3(const double* a, const double* b)
{
	const double dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2];
	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// C · ∇C · ∇²C 오차: C 는 √a 로 나눈 절대, ∇C · ∇²C 는 삼각형마다 max|Δ|/max|기준|
struct Err {
	Stat c, g, h;
	void Add(const double* C, const double* g_, const double* H, const double* Cr, const double* gr, const double* Hr, double scale)
	{
		c.Add(MaxAbsDiff(C, Cr, NC) / scale);
		g.Add(Rel(MaxAbsDiff(g_, gr, NC * N), MaxAbs(gr, NC * N)));
		h.Add(Rel(MaxAbsDiff(H, Hr, NC * NH), MaxAbs(Hr, NC * NH)));
	}
};

// 중심차분: gFD[b] = ∂C/∂x_b, HFD[(a ≤ b)] = ∂g_a/∂x_b. 분모 = 실제 간격
static void FiniteDiffN(EvalFn f, const double* x, const double* ev, double step, double* gFD, double* HFD)
{
	double xp[N], xm[N], Cp[NC], Cm[NC], gp[NC * N], gm[NC * N], Hp[NC * NH], Hm[NC * NH];
	for (int b = 0; b < N; ++b) {
		std::memcpy(xp, x, sizeof xp); std::memcpy(xm, x, sizeof xm);
		xp[b] += step; xm[b] -= step;
		const double h2 = xp[b] - xm[b];
		f(xp, ev, Cp, gp, Hp);
		f(xm, ev, Cm, gm, Hm);
		for (int c = 0; c < NC; ++c) {
			gFD[c * N + b] = (Cp[c] - Cm[c]) / h2;
			for (int a = 0; a <= b; ++a) HFD[c * NH + sx::UpperIndex<N>(a, b)] = (gp[c * N + a] - gm[c * N + a]) / h2;
		}
	}
}

static void JsonErr(std::FILE* js, const char* key, const Err& e, const char* tail)
{
	std::fprintf(js, "\"%s\": {\"C_abs_over_sqrt_a\": %.3e, \"grad_rel\": %.3e, \"hess_rel\": %.3e, \"grad_rel_median\": %.3e, \"hess_rel_median\": %.3e}%s",
		key, e.c.Max(), e.g.Max(), e.h.Max(), e.g.Median(), e.h.Median(), tail);
}

// 대칭 행렬(n × n 전체, 행 우선)의 고윳값 — 순환 Jacobi(인장 판 · Provot 전단 게이트와 같은 방법). M 은 망가진다
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

// 위 삼각 저장 → 최소 · 최대 고윳값을 최대 고윳값 크기로 나눈 값
static void EigRangeUpper(const double* U, double* lo, double* hi)
{
	double M[N * N], lam[N];
	for (int a = 0; a < N; ++a)
		for (int b = 0; b < N; ++b) M[a * N + b] = U[a <= b ? sx::UpperIndex<N>(a, b) : sx::UpperIndex<N>(b, a)];
	SymEig(M, N, lam);
	double l = lam[0], h = lam[0], big = 0.;
	for (int p = 0; p < N; ++p) { l = std::min(l, lam[p]); h = std::max(h, lam[p]); big = std::max(big, std::fabs(lam[p])); }
	*lo = big > 0. ? l / big : 0.;
	*hi = big > 0. ? h / big : 0.;
}

static double MinEigRatioUpper(const double* U)
{
	double lo, hi;
	EigRangeUpper(U, &lo, &hi);
	return lo;
}

struct Cmp { Stat c, A, b, f; int n; Cmp() : n(0) {} };

static void Compare(const double* C0, const double* Cs, double cscale, const double* A0, const double* A1, const double* b0,
	const double* b1, const double* f0, const double* f1, Cmp& m)
{
	m.c.Add(MaxAbsDiff(Cs, C0, NC) / cscale);
	m.A.Add(Rel(MaxAbsDiff(A1, A0, NH), MaxAbs(A0, NH)));
	m.b.Add(Rel(MaxAbsDiff(b1, b0, N), MaxAbs(b0, N)));
	m.f.Add(Rel(MaxAbsDiff(f1, f0, N), MaxAbs(f0, N)));
	++m.n;
}

// 보고용: w_u, w_v 의 길이(명세와 같은 식을 손으로)
static void FrameLengths(const Sample& t, double* lenU, double* lenV)
{
	const double* uv = t.ev;
	const double du1 = uv[2] - uv[0], dv1 = uv[3] - uv[1], du2 = uv[4] - uv[0], dv2 = uv[5] - uv[1];
	const double det = du1 * dv2 - du2 * dv1;
	double su = 0., sv = 0.;
	for (int k = 0; k < 3; ++k) {
		const double dx1 = t.x[3 + k] - t.x[k], dx2 = t.x[6 + k] - t.x[k];
		const double wu = (dx1 * dv2 - dx2 * dv1) / det, wv = (-dx1 * du2 + dx2 * du1) / det;
		su += wu * wu; sv += wv * wv;
	}
	*lenU = std::sqrt(su); *lenV = std::sqrt(sv);
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Tri> tris = ReadTriangles("results/triangles.txt");
	const std::vector<std::vector<double> > beta = ReadRows("results/poisson_beta.txt");
	if (tris.empty() || beta.size() < tris.size()) {
		std::printf("results/triangles.txt · results/poisson_beta.txt 없음 — model/make_triangles.py · model/make_poisson_beta.py 를 먼저 돌린다\n");
		return 2;
	}
	std::vector<Sample> ts(tris.size());
	for (std::size_t s = 0; s < tris.size(); ++s) {
		if (beta[s].size() != 1) { std::printf("results/poisson_beta.txt %zu 번째 줄이 한 칸이 아니다\n", s + 1); return 2; }
		std::memcpy(ts[s].x, tris[s].x, sizeof ts[s].x);
		std::memcpy(ts[s].ev, tris[s].uv, 6 * sizeof(double));
		ts[s].ev[6] = beta[s][0];
		std::memcpy(ts[s].v, tris[s].v, sizeof ts[s].v);
	}
	const std::vector<std::vector<double> > ref = ReadRows("results/mma_ref_poisson.txt");
	std::printf("Poisson 힘 게이트 — 표본 %zu 삼각형, Mathematica 기준 %zu 줄 (성분 %d, 자유도 %d, 요소 값 %d)\n",
		ts.size(), ref.size(), static_cast<int>(NC), static_cast<int>(N), static_cast<int>(NE));
	Report rpt("results/gate_poisson.json");
	std::FILE* js = rpt.json;
	std::fprintf(js, "{\n  \"attempt\": 2, \"samples\": %zu, \"k\": %g, \"kd\": %g, \"h\": %g,\n", ts.size(), kK, kKd, kH);

	// ---- Q1 Mathematica · Q1b 원문 닫힌 꼴 -----------------------------------
	{
		Err eFlat, eStaged;
		Stat pG, pH;
		int n = 0;
		for (std::size_t s = 0; s < ref.size() && s < ts.size(); ++s) {
			if (static_cast<int>(ref[s].size()) != kMmaCols) continue;
			const double* Cr = &ref[s][0];
			const double* gr = Cr + NC;
			const double* Hr = gr + NC * N;
			const double* gP = Hr + NC * NH;
			const double* HP = gP + NC * N;
			const double scale = CScale(ts[s].ev);
			double C[NC], g[NC * N], H[NC * NH];
			S::Eval(ts[s].x, ts[s].ev, C, g, H);
			eFlat.Add(C, g, H, Cr, gr, Hr, scale);
			pG.Add(Rel(MaxAbsDiff(g, gP, NC * N), MaxAbs(gP, NC * N)));
			pH.Add(Rel(MaxAbsDiff(H, HP, NC * NH), MaxAbs(HP, NC * NH)));
			S::EvalStaged(ts[s].x, ts[s].ev, C, g, H);
			eStaged.Add(C, g, H, Cr, gr, Hr, scale);
			++n;
		}
		const bool ran = n > 0;
		std::printf("Q1 Mathematica 대조 (%d 삼각형)\n", n);
		const Err* es[2] = { &eFlat, &eStaged };
		const char* tag[2] = { "Eval", "EvalStaged" };
		for (int i = 0; i < 2; ++i) {
			std::printf(" %s\n", tag[i]);
			rpt.Gate("C/√a", ran && es[i]->c.Max() <= 1e-12, "|ΔC|/√a %.2e", es[i]->c.Max(), 1e-12);
			rpt.Gate("∇C", ran && es[i]->g.Max() <= 1e-10, "rel %.2e", es[i]->g.Max(), 1e-10);
			rpt.Gate("∇²C", ran && es[i]->h.Max() <= 1e-10, "rel %.2e", es[i]->h.Max(), 1e-10);
		}
		std::printf("Q1b 원문 식 (7) · (8) 의 닫힌 꼴 대 생성 Eval (%d 삼각형)\n", n);
		rpt.Gate("∇C", ran && pG.Max() <= 1e-10, "rel %.2e", pG.Max(), 1e-10);
		rpt.Gate("∇²C", ran && pH.Max() <= 1e-10, "rel %.2e", pH.Max(), 1e-10);
		std::fprintf(js, "  \"Q1_mathematica\": {\"triangles\": %d, ", n);
		JsonErr(js, "eval", eFlat, ", ");
		JsonErr(js, "staged", eStaged, "},\n");
		std::fprintf(js, "  \"Q1b_paper_closed_form\": {\"triangles\": %d, \"grad_rel\": %.3e, \"hess_rel\": %.3e, \"grad_rel_median\": %.3e, \"hess_rel_median\": %.3e},\n",
			n, pG.Max(), pH.Max(), pG.Median(), pH.Median());
	}

	// ---- Q2 유한차분 · Q3 단계판 대 평탄판 · 보고 ---------------------------
	Err fdFlat, fdStaged, st;
	Stat sumG, torqueG;
	int frozenDiff = 0;
	double uMin = 1e30, uMax = 0., vMin = 1e30, vMax = 0.;
	for (std::size_t s = 0; s < ts.size(); ++s) {
		const Sample& t = ts[s];
		const double scale = CScale(t.ev);
		const double step = 1e-6 * Len3(t.x, t.x + 3);
		double C[NC], g[NC * N], H[NC * NH], Cs[NC], gs[NC * N], Hs[NC * NH], gFD[NC * N], HFD[NC * NH];
		S::Eval(t.x, t.ev, C, g, H);
		S::EvalStaged(t.x, t.ev, Cs, gs, Hs);
		FiniteDiffN(&S::Eval, t.x, t.ev, step, gFD, HFD);
		fdFlat.Add(C, g, H, C, gFD, HFD, scale);
		FiniteDiffN(&S::EvalStaged, t.x, t.ev, step, gFD, HFD);
		fdStaged.Add(Cs, gs, Hs, Cs, gFD, HFD, scale);
		st.Add(Cs, gs, Hs, C, g, H, scale);

		// 보고: Σᵢ ∇ᵢC = 0(병진), Σᵢ (xᵢ − x̄) × ∇ᵢC = 0(회전)
		double xc[3];
		for (int k = 0; k < 3; ++k) xc[k] = (t.x[k] + t.x[3 + k] + t.x[6 + k]) / 3.;
		double rmax = 0.;
		for (int a = 0; a < 3; ++a) for (int k = 0; k < 3; ++k) rmax = std::max(rmax, std::fabs(t.x[3 * a + k] - xc[k]));
		const double gm = MaxAbs(g, N);
		double sg[3] = { 0., 0., 0. }, tq[3] = { 0., 0., 0. };
		for (int a = 0; a < 3; ++a) {
			const double r[3] = { t.x[3 * a] - xc[0], t.x[3 * a + 1] - xc[1], t.x[3 * a + 2] - xc[2] };
			const double* q = g + 3 * a;
			for (int k = 0; k < 3; ++k) sg[k] += q[k];
			tq[0] += r[1] * q[2] - r[2] * q[1];
			tq[1] += r[2] * q[0] - r[0] * q[2];
			tq[2] += r[0] * q[1] - r[1] * q[0];
		}
		sumG.Add(Rel(MaxAbs(sg, 3), gm));
		torqueG.Add(Rel(MaxAbs(tq, 3), gm * rmax));
		double Cf[NC], gf[NC * N], Hf[NC * NH];
		S::EvalFrozen(t.x, t.ev, Cf, gf, Hf);
		if (std::memcmp(Cf, C, sizeof Cf) != 0 || std::memcmp(gf, g, sizeof gf) != 0 || std::memcmp(Hf, H, sizeof Hf) != 0) ++frozenDiff;
		double lu, lv;
		FrameLengths(t, &lu, &lv);
		uMin = std::min(uMin, lu); uMax = std::max(uMax, lu);
		vMin = std::min(vMin, lv); vMax = std::max(vMax, lv);
	}
	std::printf("Q2 유한차분 (%zu 삼각형)\n", ts.size());
	rpt.Gate("∇C", fdFlat.g.Max() <= 1e-6, "Eval rel max %.2e", fdFlat.g.Max(), 1e-6);
	rpt.Gate("∇²C", fdFlat.h.Max() <= 1e-6, "Eval rel max %.2e", fdFlat.h.Max(), 1e-6);
	rpt.Gate("∇C", fdStaged.g.Max() <= 1e-6, "EvalStaged rel max %.2e", fdStaged.g.Max(), 1e-6);
	rpt.Gate("∇²C", fdStaged.h.Max() <= 1e-6, "EvalStaged rel max %.2e", fdStaged.h.Max(), 1e-6);
	std::printf("Q3 단계판 대 평탄판 (%zu 삼각형)\n", ts.size());
	rpt.Gate("C/√a", st.c.Max() <= kLimit, "|ΔC|/√a %.2e", st.c.Max(), kLimit);
	rpt.Gate("∇C", st.g.Max() <= kLimit, "rel %.2e", st.g.Max(), kLimit);
	rpt.Gate("∇²C", st.h.Max() <= kLimit, "rel %.2e", st.h.Max(), kLimit);
	std::fprintf(js, "  \"Q2_finite_difference\": {\"eval\": {\"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e}, \"staged\": {\"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e}},\n",
		fdFlat.g.Max(), fdFlat.h.Max(), fdStaged.g.Max(), fdStaged.h.Max());
	std::fprintf(js, "  ");
	JsonErr(js, "Q3_staged_vs_flat", st, ",\n");

	// ---- Q4 EvalImplicit 대 AssembleImplicit · Q6 · Q7 · 보고 ---------------
	Stat eA, eB, eF, eBrel, gnBs;   // eB · gnBs = 지움을 뺀 척도의 b(2 차 판정), eBrel · gn.b = 1 차의 상대 오차(보고)
	Cmp gn;
	int nPos = 0, nNeg = 0, nBetaPos = 0, nBetaNeg = 0, nBetaZero = 0;
	int nIndefH = 0, nIndefEPos = 0, nIndefENeg = 0, nIndefTPos = 0, nIndefGN = 0;
	double gnMin = 1e30;
	Stat diffGE, eigE;
	for (std::size_t s = 0; s < ts.size(); ++s) {
		const Sample& t = ts[s];
		const double sa = CScale(t.ev);
		double C[NC], g[NC * N], H[NC * NH];
		S::Eval(t.x, t.ev, C, g, H);
		// 2 차의 b 척도: b 를 이루는 항들의 크기 합. S = k|C| + (kd + h k) Σ_a|∇_a C v_a|
		double CdAbs = 0., HvAbs = 0.;
		for (int a = 0; a < N; ++a) CdAbs += std::fabs(g[a] * t.v[a]);
		for (int a = 0; a < N; ++a) {
			double r = 0.;
			for (int b = 0; b < N; ++b) r += std::fabs(H[a <= b ? sx::UpperIndex<N>(a, b) : sx::UpperIndex<N>(b, a)] * t.v[b]);
			HvAbs = std::max(HvAbs, r);
		}
		const double bScaleGN = kH * (kK * std::fabs(C[0]) + (kKd + kH * kK) * CdAbs) * MaxAbs(g, N);
		const double bScaleE = bScaleGN + kH * kH * (kK * std::fabs(C[0]) + kKd * CdAbs) * HvAbs;
		// Q4
		sx::ImplicitBlock<N> blk;
		sx::AssembleImplicit<N>(NC, C, g, H, t.v, kK, kKd, kH, blk);
		double Ce[NC], fe[N], Ae[NH], be[N];
		S::EvalImplicit(t.x, t.ev, t.v, kPrm, Ce, fe, Ae, be);
		eA.Add(Rel(MaxAbsDiff(Ae, blk.A, NH), MaxAbs(blk.A, NH)));
		eB.Add(Rel(MaxAbsDiff(be, blk.b, N), bScaleE));
		eBrel.Add(Rel(MaxAbsDiff(be, blk.b, N), MaxAbs(blk.b, N)));
		eF.Add(Rel(MaxAbsDiff(fe, blk.f, N), MaxAbs(blk.f, N)));
		// Q6: 기준 = AssembleImplicit(C, ∇C, ∇²C = 0)
		double H0[NC * NH];
		for (int i = 0; i < NC * NH; ++i) H0[i] = 0.;
		sx::ImplicitBlock<N> refG;
		sx::AssembleImplicit<N>(NC, C, g, H0, t.v, kK, kKd, kH, refG);
		double Cg[NC], fg[N], Ag[NH], bg[N];
		S::EvalImplicitGN(t.x, t.ev, t.v, kPrm, Cg, fg, Ag, bg);
		Compare(C, Cg, sa, refG.A, Ag, refG.b, bg, refG.f, fg, gn);
		gnBs.Add(Rel(MaxAbsDiff(bg, refG.b, N), bScaleGN));
		const bool pos = C[0] > 0.;
		if (pos) ++nPos; else if (C[0] < 0.) ++nNeg;
		const double bt = t.ev[6];
		if (bt > 0.) ++nBetaPos; else if (bt < 0.) ++nBetaNeg; else ++nBetaZero;
		// Q7 과 보고
		const double rG = MinEigRatioUpper(Ag);
		gnMin = std::min(gnMin, rG);
		if (rG < kPsd) ++nIndefGN;
		double lo, hi;
		EigRangeUpper(H, &lo, &hi);
		if (lo < kPsd && hi > -kPsd) ++nIndefH;               // ∇²C 에 음 · 양 고윳값이 다 있다
		const double rE = MinEigRatioUpper(Ae);
		if (rE < kPsd) { if (pos) ++nIndefEPos; else ++nIndefENeg; eigE.Add(-rE); }
		if (pos) {   // 인장 C > 0 판의 처방을 여기에 썼다면: A_GN + h² k C ∇²C
			double At[NH];
			for (int i = 0; i < NH; ++i) At[i] = refG.A[i] + kH * kH * kK * C[0] * H[i];
			if (MinEigRatioUpper(At) < kPsd) ++nIndefTPos;
		}
		diffGE.Add(Rel(MaxAbsDiff(Ag, Ae, NH), MaxAbs(Ae, NH)));
	}
	std::printf("Q4 EvalImplicit 대 AssembleImplicit (%zu 삼각형, k = %g, kd = %g, h = %g)\n", ts.size(), kK, kKd, kH);
	rpt.Gate("A", eA.Max() <= kLimit, "rel max %.2e", eA.Max(), kLimit);
	rpt.Gate("b", eB.Max() <= kLimit, "지움 뺀 척도 max %.2e", eB.Max(), kLimit);
	rpt.Gate("f+d", eF.Max() <= kLimit, "rel max %.2e", eF.Max(), kLimit);
	std::printf("  보고 — b 의 상대 오차(1 차의 척도) max %.2e\n", eBrel.Max());
	std::fprintf(js, "  \"Q4_implicit\": {\"A_rel\": %.3e, \"b_scaled\": %.3e, \"f_rel\": %.3e, \"A_rel_median\": %.3e, \"b_scaled_median\": %.3e, \"f_rel_median\": %.3e, "
		"\"report_b_rel_first_attempt\": %.3e},\n",
		eA.Max(), eB.Max(), eF.Max(), eA.Median(), eB.Median(), eF.Median(), eBrel.Max());

	// ---- Q5 SXMesh 삼각형 루프 (요소 값 변수 둘: uv · β) ---------------------
	{
		std::mt19937_64 rng(7);
		const int G = 200;
		const Grid grid = MakeGrid(G, rng);
		const int nNode = grid.nNode, nTri = grid.nTri;
		const double cyc[4] = { 1.0, 0.0, -1.0, 0.6 };
		std::vector<double> EV(static_cast<std::size_t>(NE) * nTri);   // 쉬는 uv = 평평한 격자(간격 1 cm), β = 1 · 0 · −1 · 0.6 을 돌려 쓴다
		for (int s = 0; s < nTri; ++s) {
			for (int a = 0; a < 3; ++a) {
				const int nd = grid.tris[3 * s + a];
				EV[NE * s + 2 * a] = 0.01 * (nd % G);
				EV[NE * s + 2 * a + 1] = 0.01 * (nd / G);
			}
			EV[NE * s + 6] = cyc[s % 4];
		}
		SXMesh mesh;
		sx::CreateVariables<S>(mesh, nNode, nTri);
		SXNodeNamedContainer& px = mesh.GetNodeV(S::PositionVar());
		SXNodeNamedContainer& pv = mesh.GetNodeV(S::VelocityVar());
		SXElementNamedContainer& pt = mesh.GetElementV(S::StencilVar());
		SXNodeNamedContainer& puv = mesh.GetNodeV(S::ElemVar(0));
		SXNodeNamedContainer& pbt = mesh.GetNodeV(S::ElemVar(1));
		for (int n = 0; n < nNode; ++n) { px.Add(n, &grid.X[3 * n], 3); pv.Add(n, &grid.V[3 * n], 3); }
		for (int s = 0; s < nTri; ++s) {
			pt.Add(s, &grid.tris[3 * s], 3);
			for (int c = 0; c < 6; ++c) puv[s][c] = EV[NE * s + c];
			pbt[s][0] = EV[NE * s + 6];
		}
		sx::StencilView<S> view(mesh);
		int mismatch = 0;
		for (int s = 0; s < nTri; ++s) {
			double xa[N], va[N], eva[NE], xb[N], vb[N];
			if (!view.Gather(s, xa, va, eva)) { ++mismatch; continue; }
			for (int a = 0; a < 3; ++a) for (int c = 0; c < 3; ++c) { xb[3 * a + c] = grid.X[3 * grid.tris[3 * s + a] + c]; vb[3 * a + c] = grid.V[3 * grid.tris[3 * s + a] + c]; }
			double Ca[NC], ga[NC * N], Ha[NC * NH], Cb[NC], gb[NC * N], Hb[NC * NH];
			S::Eval(xa, eva, Ca, ga, Ha);
			S::Eval(xb, &EV[NE * s], Cb, gb, Hb);
			if (std::memcmp(Ca, Cb, sizeof Ca) != 0 || std::memcmp(ga, gb, sizeof ga) != 0 || std::memcmp(Ha, Hb, sizeof Ha) != 0 ||
				std::memcmp(va, vb, sizeof va) != 0 || std::memcmp(eva, &EV[NE * s], sizeof eva) != 0) ++mismatch;
		}
		std::printf("Q5 SXMesh 삼각형 루프 (노드 %d, 삼각형 %d, 요소 값 변수 %d 개)\n", nNode, nTri, static_cast<int>(S::kElemVars));
		rpt.Gate("루프", mismatch == 0, "불일치 %.0f", static_cast<double>(mismatch), 0.);
		std::fprintf(js, "  \"Q5_stencil_loop\": {\"nodes\": %d, \"triangles\": %d, \"mismatch\": %d},\n", nNode, nTri, mismatch);
	}

	// ---- Q6 · Q7 판정 --------------------------------------------------------
	std::printf("Q6 EvalImplicitGN 대 AssembleImplicit(C, ∇C, ∇²C = 0) (%d 개)\n", gn.n);
	rpt.Gate("C", gn.n > 0 && gn.c.Max() <= kLimit, "|ΔC|/√a %.2e", gn.c.Max(), kLimit);
	rpt.Gate("A", gn.n > 0 && gn.A.Max() <= kLimit, "rel max %.2e", gn.A.Max(), kLimit);
	rpt.Gate("b", gn.n > 0 && gnBs.Max() <= kLimit, "지움 뺀 척도 max %.2e", gnBs.Max(), kLimit);
	rpt.Gate("f+d", gn.n > 0 && gn.f.Max() <= kLimit, "rel max %.2e", gn.f.Max(), kLimit);
	std::printf("  보고 — b 의 상대 오차(1 차의 척도) max %.2e\n", gn.b.Max());
	std::printf(" Q6 갈래 — C > 0 %d 개 · C < 0 %d 개, β > 0 %d 개 · β < 0 %d 개 · β = 0 %d 개\n", nPos, nNeg, nBetaPos, nBetaNeg, nBetaZero);
	rpt.Gate("C>0", nPos >= kMinBranch, "%.0f 개", nPos, kMinBranch);
	rpt.Gate("C<0", nNeg >= kMinBranch, "%.0f 개", nNeg, kMinBranch);
	rpt.Gate("β>0", nBetaPos >= kMinBranch, "%.0f 개", nBetaPos, kMinBranch);
	rpt.Gate("β<0", nBetaNeg >= kMinBranch, "%.0f 개", nBetaNeg, kMinBranch);
	rpt.Gate("β=0", nBetaZero >= kMinBranch, "%.0f 개", nBetaZero, kMinBranch);
	std::printf("Q7 GN 준정부호 (%zu 개)\n", ts.size());
	rpt.Gate("A_GN", gnMin >= kPsd, "최소 고윳값 / 최대 크기의 최솟값 %.2e", gnMin, kPsd);

	// 보고
	double Crest[NC], grest[NC * N], Hrest[NC * NH];
	S::Eval(ts[0].x, ts[0].ev, Crest, grest, Hrest);
	const double restC = MaxAbs(Crest, NC) / CScale(ts[0].ev);
	std::printf("  보고 — 병진 |Σ∇C|/max|∇C| max %.1e, 회전 |Σ r×∇C| max %.1e, 쉬는 모양 |C|/√a %.1e, EvalFrozen ≠ Eval %d 개, |w_u| %.3f – %.3f, |w_v| %.3f – %.3f\n",
		sumG.Max(), torqueG.Max(), restC, frozenDiff, uMin, uMax, vMin, vMax);
	std::printf("  보고 — ∇²C 부정치 %d 개, 정확판 A 부정치 C > 0 %d / %d · C ≤ 0 %d / %d (|최소|/최대 중앙 %.2e), C > 0 판 처방이면 부정치 %d / %d, GN A 부정치 %d 개, |A_GN − A_E|/|A_E| 중앙 %.2e · max %.2e\n",
		nIndefH, nIndefEPos, nPos, nIndefENeg, static_cast<int>(ts.size()) - nPos, eigE.Median(), nIndefTPos, nPos, nIndefGN,
		diffGE.Median(), diffGE.Max());

	std::fprintf(js, "  \"Q6\": {\"samples\": %d, \"C_abs_over_sqrt_a\": %.3e, \"A_rel\": %.3e, \"b_scaled\": %.3e, \"f_rel\": %.3e, "
		"\"A_rel_median\": %.3e, \"b_scaled_median\": %.3e, \"f_rel_median\": %.3e, \"report_b_rel_first_attempt\": %.3e,\n"
		"    \"C_pos\": %d, \"C_neg\": %d, \"beta_pos\": %d, \"beta_neg\": %d, \"beta_zero\": %d, \"min_branch\": %d},\n",
		gn.n, gn.c.Max(), gn.A.Max(), gnBs.Max(), gn.f.Max(), gn.A.Median(), gnBs.Median(), gn.f.Median(), gn.b.Max(),
		nPos, nNeg, nBetaPos, nBetaNeg, nBetaZero, kMinBranch);
	std::fprintf(js, "  \"Q7\": {\"min_eig_ratio\": %.3e, \"limit\": %g, \"indefinite\": %d},\n", gnMin, kPsd, nIndefGN);
	std::fprintf(js, "  \"report\": {\"sum_grad_rel_max\": %.3e, \"torque_grad_rel_max\": %.3e, \"rest_C_over_sqrt_a\": %.3e, \"frozen_differs\": %d, "
		"\"wu_len_min\": %.6f, \"wu_len_max\": %.6f, \"wv_len_min\": %.6f, \"wv_len_max\": %.6f,\n"
		"    \"indefinite_hessian_C\": %d, \"indefinite_A_exact_C_pos\": %d, \"indefinite_A_exact_C_nonpos\": %d, "
		"\"indefinite_A_exact_eig_ratio_median\": %.3e, \"indefinite_A_tension_rule_C_pos\": %d, "
		"\"A_GN_vs_A_E_rel_median\": %.3e, \"A_GN_vs_A_E_rel_max\": %.3e},\n",
		sumG.Max(), torqueG.Max(), restC, frozenDiff, uMin, uMax, vMin, vMax, nIndefH, nIndefEPos, nIndefENeg, eigE.Median(), nIndefTPos,
		diffGE.Median(), diffGE.Max());
	std::fprintf(js, "  \"gates_failed\": %d\n}\n", rpt.fail);
	std::printf("\nPoisson 힘 게이트 Q1 – Q7: %s (실패 %d) → results/gate_poisson.json\n", rpt.fail == 0 ? "GO" : "NO-GO", rpt.fail);
	return rpt.fail == 0 ? 0 : 1;
}
