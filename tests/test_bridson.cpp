//---------------------------------------------------------------------------
// test_bridson.cpp — Bridson 굽힘 게이트(docs/GATE_BRIDSON.md).
// 명세 쪽 B1 – B4 를 판정한다. 엔진이 있으면(/DSX_WITH_ENGINE, engine_macri.obj · SDEmysse.obj 와 링크) B5 의 θ 를
// 판정하고 B6 · B7 을 분류하며, B8 시간은 이 실행 한 번의 값으로 보고만 한다.
// B5 의 옮김(SHA-256)은 tests/engine_ref/extract_engine.py, B8 판정은 tests/timing_bases.py bridson(고정 로드 주소 8 곳).
// 실행 디렉터리 = 리포 루트. results/hinges.txt · mma_ref_bridson.txt · mma_ref.txt 를 읽고 results/gate_bridson.json 을 쓴다.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BridsonBend.h"
#include "../generated/BaraffBendCompact.h"
#include "gate_common.h"
#ifdef SX_WITH_ENGINE
#include "engine_ref/engine_macri.h"
#endif

#include <chrono>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;

typedef sxgen::BridsonBend Spec;
typedef sxgen::BaraffBendCompact Theta;   // θ · ∇θ (θ0 = 0) — 보고용: 힘의 척도 · 오차 나눔 · 강성 비
static const int N = Spec::kDof;          // 12
static const int NH = Spec::kHess;        // 78
static const int NM = 1 + N + NH;         // Mathematica 한 줄 91
static const int NI = 1 + N + NH + N;     // EvalImplicit 출력 103 개(C · f + d · A · b)
static const double kPrmK = 1.0, kPrmKd = 0.1, kPrmH = 0.01;   // 판정문에 고정한 계수(조립 게이트와 같음)
static const double kPrm[3] = { kPrmK, kPrmKd, kPrmH };

// 위 삼각 78 → 12 × 12 전체(대칭)
static void Full(const double* Hu, double* Hf)
{
	for (int a = 0; a < N; ++a)
		for (int b = 0; b < N; ++b) Hf[a * N + b] = Hu[a <= b ? sx::UpperIndex<N>(a, b) : sx::UpperIndex<N>(b, a)];
}

static void Sub3(const double* a, const double* b, double* r) { r[0] = a[0] - b[0]; r[1] = a[1] - b[1]; r[2] = a[2] - b[2]; }
static double Dot3(const double* a, const double* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
static double Norm3(const double* a) { return std::sqrt(Dot3(a, a)); }
static void Cross3(const double* a, const double* b, double* r)
{
	r[0] = a[1] * b[2] - a[2] * b[1];
	r[1] = a[2] * b[0] - a[0] * b[2];
	r[2] = a[0] * b[1] - a[1] * b[0];
}

// 명세의 법선 nA = (xi − xh) × (xj − xh), nB = (xj − xk) × (xi − xk), 변 e = xj − xi
static void Normals(const double* x, double* nA, double* nB, double* e)
{
	const double *xh = x, *xi = x + 3, *xj = x + 6, *xk = x + 9;
	double a[3], b[3];
	Sub3(xi, xh, a); Sub3(xj, xh, b); Cross3(a, b, nA);
	Sub3(xj, xk, a); Sub3(xi, xk, b); Cross3(a, b, nB);
	Sub3(xj, xi, e);
}

// 표본 힌지 가중치 w̄ = 3|e|² / (½|nA| + ½|nB|) (판정문 「판정 기준」, Mathematica 기준과 같은 식)
static double Weight(const double* x)
{
	double nA[3], nB[3], e[3];
	Normals(x, nA, nB, e);
	return 3. * Dot3(e, e) / (0.5 * Norm3(nA) + 0.5 * Norm3(nB));
}

#ifdef SX_WITH_ENGINE
// ---- 읽어서 세운 엔진 Jacobian 예측(판정문 「읽어서 세운 예측」) ---------------------------------
// 엔진의 u_m = α_m n̂A + β_m n̂B (명세 순서 h, i, j, k). α · β 는 Oh 함수의 식을 명세 꼭짓점으로 옮겨 쓴 것
static void BridsonU(const double* x, double* alpha, double* beta, double* nhA, double* nhB)
{
	const double *xh = x, *xi = x + 3, *xj = x + 6, *xk = x + 9;
	double nA[3], nB[3], e[3], t[3];
	Normals(x, nA, nB, e);
	const double la = Norm3(nA), lb = Norm3(nB), le = Norm3(e);
	for (int c = 0; c < 3; ++c) { nhA[c] = nA[c] / la; nhB[c] = nB[c] / lb; }
	Sub3(xi, xh, t); const double dot0332 = Dot3(e, t);      // qA[0]·qA[2] = (xJ − xI)·(xI − xH)
	Sub3(xk, xi, t); const double dot1332 = -Dot3(e, t);     // −qA[0]·qB[2]
	Sub3(xh, xj, t); const double dot2032 = Dot3(e, t);      // qA[0]·qA[1]
	Sub3(xj, xk, t); const double dot2132 = -Dot3(e, t);     // −qA[0]·qB[1]
	alpha[0] = le / la; alpha[1] = dot2032 / (le * la); alpha[2] = dot0332 / (le * la); alpha[3] = 0.;
	beta[0] = 0.; beta[1] = dot2132 / (le * lb); beta[2] = dot1332 / (le * lb); beta[3] = le / lb;
}

// 같은 법선끼리의 항만 남긴 −k w (Σ α_mα_n n̂An̂Aᵀ + β_mβ_n n̂Bn̂Bᵀ). dropIJ 면 (i, j) 블록을 지우고 대각 = −행 합
static void SameNormal(const double* alpha, const double* beta, const double* nhA, const double* nhB, double kw, bool dropIJ, double* M)
{
	for (int m = 0; m < 4; ++m) for (int n = 0; n < 4; ++n) for (int s = 0; s < 3; ++s) for (int t = 0; t < 3; ++t)
		M[(3 * m + s) * N + 3 * n + t] = -kw * (alpha[m] * alpha[n] * nhA[s] * nhA[t] + beta[m] * beta[n] * nhB[s] * nhB[t]);
	if (!dropIJ) return;
	for (int s = 0; s < 3; ++s) for (int t = 0; t < 3; ++t) { M[(3 + s) * N + 6 + t] = 0.; M[(6 + s) * N + 3 + t] = 0.; }
	for (int m = 0; m < 4; ++m) for (int s = 0; s < 3; ++s) for (int t = 0; t < 3; ++t) {
		double sum = 0.;
		for (int n = 0; n < 4; ++n) if (n != m) sum += M[(3 * m + s) * N + 3 * n + t];
		M[(3 * m + s) * N + 3 * m + t] = -sum;
	}
}

static double Quad(const double* M, const double* u)
{
	double q = 0.;
	for (int a = 0; a < N; ++a) for (int b = 0; b < N; ++b) q += u[a] * M[a * N + b] * u[b];
	return q;
}

// ---- 시간 ----------------------------------------------------------------------
typedef std::chrono::steady_clock Clock;

static double TimeEngine(engine_ref::Bridson* eb, double& checksum)
{
	const Clock::time_point t0 = Clock::now();
	engine_ref::RunBridson(eb);
	const double ns = std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / engine_ref::SizeBridson(eb);
	checksum += engine_ref::ChecksumBridson(eb);
	return ns;
}

// 모으기 + EvalImplicit + 출력 103 개를 힌지별 배열에 저장 (B8 의 생성 쪽)
static double TimeImplicitStore(const sx::StencilView<Spec>& view, std::vector<double>& out, double& checksum)
{
	const int nHinge = view.size();
	const Clock::time_point t0 = Clock::now();
	for (int s = 0; s < nHinge; ++s) {
		double x[N], v[N], ev[1];
		view.Gather(s, x, v, ev);
		double* o = &out[static_cast<std::size_t>(s) * NI];
		Spec::EvalImplicit(x, ev, v, kPrm, o, o + 1, o + 1 + N, o + 1 + N + NH);
	}
	const double ns = std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / nHinge;
	for (int s = 0; s < nHinge; s += 97) checksum += out[static_cast<std::size_t>(s) * NI + 20];
	return ns;
}

// 보고: 모으기 + EvalStaged + C · ∇C · ∇²C(91 개) 저장
static double TimeStagedStore(const sx::StencilView<Spec>& view, std::vector<double>& out, double& checksum)
{
	const int nHinge = view.size();
	const Clock::time_point t0 = Clock::now();
	for (int s = 0; s < nHinge; ++s) {
		double x[N], v[N], ev[1];
		view.Gather(s, x, v, ev);
		double* o = &out[static_cast<std::size_t>(s) * NI];
		Spec::EvalStaged(x, ev, o, o + 1, o + 1 + N);
	}
	const double ns = std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / nHinge;
	for (int s = 0; s < nHinge; s += 97) checksum += out[static_cast<std::size_t>(s) * NI + 20];
	return ns;
}

// 보고: 모으기 + EvalImplicit + f · b 를 노드 배열로 흩뿌리기(엔진 비교 게이트의 EvalImplicit 보고와 같은 꼴)
static double TimeImplicitScatter(const sx::StencilView<Spec>& view, std::vector<double>& Fn, std::vector<double>& Bn, double& checksum)
{
	std::fill(Fn.begin(), Fn.end(), 0.);
	std::fill(Bn.begin(), Bn.end(), 0.);
	const int nHinge = view.size();
	const Clock::time_point t0 = Clock::now();
	for (int s = 0; s < nHinge; ++s) {
		double x[N], v[N], ev[1], C, f[N], A[NH], b[N];
		view.Gather(s, x, v, ev);
		Spec::EvalImplicit(x, ev, v, kPrm, &C, f, A, b);
		const int* nd = view.Nodes(s);
		for (int a = 0; a < 4; ++a) for (int c = 0; c < 3; ++c) { Fn[3 * nd[a] + c] += f[3 * a + c]; Bn[3 * nd[a] + c] += b[3 * a + c]; }
		checksum += A[13];
	}
	return std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / nHinge;
}

struct Timing {
	double engine, engineSerial, genImplicit, genStaged, genScatter;
	int nHinge, nPair;
	double engineBytes;
};

static Timing TimeGrid(int G, bool full)
{
	std::mt19937_64 rng(7);
	const Grid grid = MakeGrid(G, rng);
	std::vector<double> W(grid.nHinge), L0(grid.nHinge);
	for (int s = 0; s < grid.nHinge; ++s) {   // 격자 모양을 쉬는 모양으로 본다: w̄ · Length0 를 여기서 정한다
		double x[N];
		for (int a = 0; a < 4; ++a) for (int c = 0; c < 3; ++c) x[3 * a + c] = grid.X[3 * grid.hinges[4 * s + a] + c];
		W[s] = Weight(x);
		L0[s] = EdgeLength(x);
	}
	SXMesh mesh;
	FillMesh<Spec>(mesh, grid, W);
	const sx::StencilView<Spec> view(mesh);

	std::vector<double> X(grid.X), V(grid.V);   // 엔진은 노드 배열을 가리키기만 한다
	engine_ref::Bridson* eb = engine_ref::CreateBridson(grid.nHinge, &grid.hinges[0], grid.nNode, &X[0], &V[0], &L0[0]);

	std::vector<double> out(static_cast<std::size_t>(grid.nHinge) * NI, 0.);
	std::vector<double> Fn(3 * grid.nNode), Bn(3 * grid.nNode);
	Timing t;
	t.engine = t.engineSerial = t.genImplicit = t.genStaged = t.genScatter = 1e30;
	t.nHinge = grid.nHinge;
	t.nPair = engine_ref::PairsBridson(eb);
	t.engineBytes = engine_ref::BytesBridson(eb);
	double checksum = 0.;
	for (int rep5 = 0; rep5 < 5; ++rep5) {   // 변형을 번갈아 5 회, 변형마다 최소
		engine_ref::SetSerialBridson(eb, false);
		t.engine      = std::min(t.engine,      TimeEngine(eb, checksum));
		t.genImplicit = std::min(t.genImplicit, TimeImplicitStore(view, out, checksum));
		if (!full) continue;
		engine_ref::SetSerialBridson(eb, true);
		t.engineSerial = std::min(t.engineSerial, TimeEngine(eb, checksum));
		t.genStaged    = std::min(t.genStaged,    TimeStagedStore(view, out, checksum));
		t.genScatter   = std::min(t.genScatter,   TimeImplicitScatter(view, Fn, Bn, checksum));
	}
	std::printf("  (격자 G = %d, checksum %.3e)\n", G, checksum);
	engine_ref::DestroyBridson(eb);
	return t;
}
#endif   // SX_WITH_ENGINE

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Hinge> hs = ReadHinges("results/hinges.txt");
	const std::vector<std::vector<double> > ref = ReadRows("results/mma_ref_bridson.txt");
	if (hs.empty() || ref.empty()) { std::printf("results/hinges.txt · mma_ref_bridson.txt 없음\n"); return 2; }
	const int nS = static_cast<int>(hs.size());
	Report rpt("results/gate_bridson.json");
	std::FILE* js = rpt.json;
	std::printf("Bridson 굽힘 게이트 — 표본 %d 힌지, Mathematica 기준 %zu 개\n", nS, ref.size());
	std::fprintf(js, "{\n  \"engine_commit\": \"378d20c\",\n  \"samples\": %d,\n", nS);

	// ---- B1 Mathematica ----------------------------------------------------------
	{
		Stat c[2], gr[2], hr[2];
		int n = 0;
		for (int s = 0; s < nS && s < static_cast<int>(ref.size()); ++s) {
			const std::vector<double>& r = ref[s];
			if (static_cast<int>(r.size()) != NM) continue;
			const double w = Weight(hs[s].x), sw = std::sqrt(w);
			for (int k = 0; k < 2; ++k) {
				double C, g[N], H[NH];
				(k == 0 ? &Spec::Eval : &Spec::EvalStaged)(hs[s].x, &w, &C, g, H);
				c[k].Add(std::fabs(C - r[0]) / sw);
				gr[k].Add(Rel(MaxAbsDiff(g, &r[1], N), MaxAbs(&r[1], N)));
				hr[k].Add(Rel(MaxAbsDiff(H, &r[1 + N], NH), MaxAbs(&r[1 + N], NH)));
			}
			++n;
		}
		const bool ran = n > 0;
		std::printf("B1 Mathematica 대조 (%d 힌지)\n", n);
		static const char* tag[2] = { "Eval", "EvalStaged" };
		for (int k = 0; k < 2; ++k) {
			std::printf(" %s\n", tag[k]);
			rpt.Gate("C/√w", ran && c[k].Max() <= 1e-12, "|ΔC|/√w̄ %.2e", c[k].Max(), 1e-12);
			rpt.Gate("∇C", ran && gr[k].Max() <= 1e-10, "rel %.2e", gr[k].Max(), 1e-10);
			rpt.Gate("∇²C", ran && hr[k].Max() <= 1e-10, "rel %.2e", hr[k].Max(), 1e-10);
		}
		std::fprintf(js, "  \"B1_mathematica\": {\"hinges\": %d", n);
		for (int k = 0; k < 2; ++k)
			std::fprintf(js, ", \"%s\": {\"C_abs_over_sqrt_w\": %.3e, \"grad_rel\": %.3e, \"hess_rel\": %.3e, \"grad_rel_median\": %.3e, \"hess_rel_median\": %.3e}",
				k == 0 ? "eval" : "staged", c[k].Max(), gr[k].Max(), hr[k].Max(), gr[k].Median(), hr[k].Median());
		std::fprintf(js, "},\n");
	}

	// ---- B2 유한차분 · B3 단계판 대 평탄판 · B4 조립 ------------------------------------
	{
		Stat fdG[2], fdH[2], stC, stG, stH, eA, eB, eF, sumG;
		for (int s = 0; s < nS; ++s) {
			const Hinge& h = hs[s];
			const double w = Weight(h.x), sw = std::sqrt(w);
			const double step = 1e-6 * EdgeLength(h.x);
			double C, g[N], H[NH], Cs, gs[N], Hs[NH], gFD[N], HFD[NH];
			Spec::Eval(h.x, &w, &C, g, H);
			Spec::EvalStaged(h.x, &w, &Cs, gs, Hs);
			FiniteDiff<N>(&Spec::Eval, h.x, &w, step, gFD, HFD);
			fdG[0].Add(Rel(MaxAbsDiff(g, gFD, N), MaxAbs(g, N)));
			fdH[0].Add(Rel(MaxAbsDiff(H, HFD, NH), MaxAbs(H, NH)));
			FiniteDiff<N>(&Spec::EvalStaged, h.x, &w, step, gFD, HFD);
			fdG[1].Add(Rel(MaxAbsDiff(gs, gFD, N), MaxAbs(gs, N)));
			fdH[1].Add(Rel(MaxAbsDiff(Hs, HFD, NH), MaxAbs(Hs, NH)));
			stC.Add(std::fabs(Cs - C) / sw);
			stG.Add(Rel(MaxAbsDiff(gs, g, N), MaxAbs(g, N)));
			stH.Add(Rel(MaxAbsDiff(Hs, H, NH), MaxAbs(H, NH)));
			double sg[3] = { 0., 0., 0. };
			for (int a = 0; a < 4; ++a) for (int k = 0; k < 3; ++k) sg[k] += g[3 * a + k];
			sumG.Add(Rel(MaxAbs(sg, 3), MaxAbs(g, N)));

			sx::ImplicitBlock<N> blk;
			sx::AssembleImplicit<N>(1, &C, g, H, h.v, kPrmK, kPrmKd, kPrmH, blk);
			double Ci, f[N], A[NH], b[N];
			Spec::EvalImplicit(h.x, &w, h.v, kPrm, &Ci, f, A, b);
			eA.Add(Rel(MaxAbsDiff(A, blk.A, NH), MaxAbs(blk.A, NH)));
			eB.Add(Rel(MaxAbsDiff(b, blk.b, N), MaxAbs(blk.b, N)));
			eF.Add(Rel(MaxAbsDiff(f, blk.f, N), MaxAbs(blk.f, N)));
		}
		std::printf("B2 유한차분 (%d 힌지)\n", nS);
		rpt.Gate("∇C", fdG[0].Max() <= 1e-6, "Eval rel max %.2e", fdG[0].Max(), 1e-6);
		rpt.Gate("∇²C", fdH[0].Max() <= 1e-6, "Eval rel max %.2e", fdH[0].Max(), 1e-6);
		rpt.Gate("∇C", fdG[1].Max() <= 1e-6, "EvalStaged rel max %.2e", fdG[1].Max(), 1e-6);
		rpt.Gate("∇²C", fdH[1].Max() <= 1e-6, "EvalStaged rel max %.2e", fdH[1].Max(), 1e-6);
		std::printf("B3 단계판 대 평탄판 (%d 힌지)\n", nS);
		rpt.Gate("C/√w", stC.Max() <= 1e-12, "|ΔC|/√w̄ %.2e", stC.Max(), 1e-12);
		rpt.Gate("∇C", stG.Max() <= 1e-12, "rel %.2e", stG.Max(), 1e-12);
		rpt.Gate("∇²C", stH.Max() <= 1e-12, "rel %.2e", stH.Max(), 1e-12);
		std::printf("B4 EvalImplicit 대 AssembleImplicit (%d 힌지, k = %g, kd = %g, h = %g)\n", nS, kPrmK, kPrmKd, kPrmH);
		rpt.Gate("A", eA.Max() <= 1e-12, "rel max %.2e", eA.Max(), 1e-12);
		rpt.Gate("b", eB.Max() <= 1e-12, "rel max %.2e", eB.Max(), 1e-12);
		rpt.Gate("f+d", eF.Max() <= 1e-12, "rel max %.2e", eF.Max(), 1e-12);
		std::printf("  보고 — 병진 |Σ∇C|/max|∇C| max %.1e; B4 중앙값 A %.1e · b %.1e · f+d %.1e\n",
			sumG.Max(), eA.Median(), eB.Median(), eF.Median());
		std::fprintf(js, "  \"B2_finite_difference\": {\"eval\": {\"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e}, \"staged\": {\"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e}},\n",
			fdG[0].Max(), fdH[0].Max(), fdG[1].Max(), fdH[1].Max());
		std::fprintf(js, "  \"B3_staged_vs_flat\": {\"C_abs_over_sqrt_w\": %.3e, \"grad_rel\": %.3e, \"hess_rel\": %.3e},\n", stC.Max(), stG.Max(), stH.Max());
		std::fprintf(js, "  \"B4_implicit\": {\"A_rel\": %.3e, \"b_rel\": %.3e, \"f_rel\": %.3e, \"A_rel_median\": %.3e, \"b_rel_median\": %.3e, \"f_rel_median\": %.3e},\n",
			eA.Max(), eB.Max(), eF.Max(), eA.Median(), eB.Median(), eF.Median());
		std::fprintf(js, "  \"report_sum_grad_rel_max\": %.3e,\n", sumG.Max());
	}
	const int specFail = rpt.fail;

#ifdef SX_WITH_ENGINE
	// ---- 엔진: 표본마다 자기 노드 4 개, 그리고 유한차분용 사본 24 개(좌표 하나를 ± 간격) -----------
	const std::vector<std::vector<double> > refTheta = ReadRows("results/mma_ref.txt");   // 착수 게이트 기준(θ − θ0)
	const int kCopies = 1 + 2 * N;
	const int nE = nS * kCopies;
	std::vector<double> X(static_cast<std::size_t>(12) * nE), V(static_cast<std::size_t>(12) * nE), L0(nE), H2(static_cast<std::size_t>(nS) * N);
	std::vector<int> vid(static_cast<std::size_t>(4) * nE);
	for (int s = 0; s < nS; ++s) {
		const Hinge& h = hs[s];
		const double step = 1e-6 * EdgeLength(h.x);
		for (int c = 0; c < kCopies; ++c) {
			const int idx = s * kCopies + c;
			double* xc = &X[static_cast<std::size_t>(12) * idx];
			for (int i = 0; i < 12; ++i) { xc[i] = h.x[i]; V[static_cast<std::size_t>(12) * idx + i] = h.v[i]; }
			if (c > 0) {
				const int b = (c - 1) / 2;
				xc[b] = (c - 1) % 2 == 0 ? h.x[b] + step : h.x[b] - step;
			}
			L0[idx] = EdgeLength(h.x);   // 쉬는 변 길이 = 표본의 |e| (사본도 같다)
			for (int m = 0; m < 4; ++m) vid[static_cast<std::size_t>(4) * idx + m] = 4 * idx + m;
		}
		for (int b = 0; b < N; ++b) H2[static_cast<std::size_t>(s) * N + b] = (h.x[b] + step) - (h.x[b] - step);
	}
	engine_ref::Bridson* eb = engine_ref::CreateBridson(nE, &vid[0], 4 * nE, &X[0], &V[0], &L0[0]);
	engine_ref::RunBridson(eb);

	// ---- B5 θ ---------------------------------------------------------------------
	Stat b5;
	int nSkip = 0;
	for (int i = 0; i < nE; ++i) if (engine_ref::SkippedBridson(eb, i)) ++nSkip;
	for (int s = 0; s < nS && s < static_cast<int>(refTheta.size()); ++s)
		if (static_cast<int>(refTheta[s].size()) == NM) b5.Add(std::fabs(engine_ref::ThetaBridson(eb, s * kCopies) + (refTheta[s][0] + hs[s].theta0)));
	std::printf("B5 엔진 θ — Mathematica %zu 힌지 (옮김 SHA-256 은 extract_engine.py)\n", b5.v.size());
	rpt.Gate("B5", !b5.v.empty() && b5.Max() <= 1e-12, "max |θ_엔진 + θ| %.2e rad", b5.Max(), 1e-12);
	std::printf("  보고 — 퇴화 가드로 건너뛴 힌지 %d 개 (사본 포함 %d 개 중)\n", nSkip, nE);
	std::fprintf(js, "  \"B5_theta\": {\"hinges\": %zu, \"abs_max\": %.3e, \"skipped\": %d},\n", b5.v.size(), b5.Max(), nSkip);

	// ---- B6 · B7 --------------------------------------------------------------------
	Stat fRel, fRelPlain, fRelPlainRest, fMom;  // B6 (Rest = 평평한 표본 0 을 뺀 것)
	double fPlain0 = 0.;
	Stat jK, jGN, jKflat, jP, rho, rhoFlat;     // B7
	Stat bud[4];                                // 오차 나눔: K→Ka, Ka→Kb, Kb→Kc, Kc→Kd(= 예측)
	Stat uSign, fdEng, fdK, fdKw, wTerm;        // u = ±∇θ, 엔진 힘의 유한차분 Jacobian
	int uPlus = 0, nFlat = 0;
	double rho0 = 0.;
	for (int s = 0; s < nS; ++s) {
		const Hinge& h = hs[s];
		const int idx = s * kCopies;
		const double w = engine_ref::WeightBridson(eb, idx);
		double C, g[N], Hu[NH], Hf[N * N], K[N * N], KGN[N * N], f[N];
		Spec::Eval(h.x, &w, &C, g, Hu);
		Full(Hu, Hf);
		for (int a = 0; a < N; ++a) {
			f[a] = -kPrmK * C * g[a];
			for (int b = 0; b < N; ++b) { KGN[a * N + b] = -kPrmK * g[a] * g[b]; K[a * N + b] = KGN[a * N + b] - kPrmK * C * Hf[a * N + b]; }
		}
		const double zero = 0.;
		double th, gt[N], Ht[NH];
		Theta::Eval(h.x, &zero, &th, gt, Ht);
		const double fScale = kPrmK * w * MaxAbs(gt, N);   // θ = π 에서의 힘 크기 척도(평평한 표본은 f = 0)
		const double kScale = MaxAbs(K, N * N);

		// B6
		double fe[N];
		engine_ref::ForceBridson(eb, idx, fe);
		fRel.Add(MaxAbsDiff(fe, f, N) / fScale);
		const double plain = Rel(MaxAbsDiff(fe, f, N), MaxAbs(f, N));
		fRelPlain.Add(plain);
		if (s == 0) fPlain0 = plain; else fRelPlainRest.Add(plain);
		double sf[3] = { 0., 0., 0. };
		for (int a = 0; a < 4; ++a) for (int k = 0; k < 3; ++k) sf[k] += fe[3 * a + k];
		fMom.Add(MaxAbs(sf, 3) / fScale);

		// B7
		double Je[N * N];
		engine_ref::JacobianBridson(eb, idx, Je);
		const double rK = MaxAbsDiff(Je, K, N * N) / kScale;
		jK.Add(rK);
		jGN.Add(Rel(MaxAbsDiff(Je, KGN, N * N), MaxAbs(KGN, N * N)));
		double alpha[4], beta[4], nhA[3], nhB[3], u[N];
		BridsonU(h.x, alpha, beta, nhA, nhB);
		for (int m = 0; m < 4; ++m) for (int k = 0; k < 3; ++k) u[3 * m + k] = alpha[m] * nhA[k] + beta[m] * nhB[k];
		double mgt[N];
		for (int a = 0; a < N; ++a) mgt[a] = -gt[a];
		const double dPlus = MaxAbsDiff(u, gt, N), dMinus = MaxAbsDiff(u, mgt, N);
		if (dPlus <= dMinus) ++uPlus;
		uSign.Add(std::min(dPlus, dMinus) / MaxAbs(gt, N));
		double Ka[N * N], Kb[N * N], Kc[N * N], Kd[N * N];
		for (int a = 0; a < N; ++a) for (int b = 0; b < N; ++b) {
			Kb[a * N + b] = -kPrmK * w * gt[a] * gt[b];
			Ka[a * N + b] = 0.5 * std::cos(0.5 * th) * Kb[a * N + b];
		}
		SameNormal(alpha, beta, nhA, nhB, kPrmK * w, false, Kc);
		SameNormal(alpha, beta, nhA, nhB, kPrmK * w, true, Kd);
		bud[0].Add(MaxAbsDiff(K, Ka, N * N) / kScale);
		bud[1].Add(MaxAbsDiff(Ka, Kb, N * N) / kScale);
		bud[2].Add(MaxAbsDiff(Kb, Kc, N * N) / kScale);
		bud[3].Add(MaxAbsDiff(Kc, Kd, N * N) / kScale);
		jP.Add(MaxAbsDiff(Je, Kd, N * N) / kScale);
		const double r = Quad(Je, gt) / Quad(K, gt);
		rho.Add(r);
		if (std::fabs(th) < 10. * 3.14159265358979323846 / 180.) { jKflat.Add(rK); rhoFlat.Add(r); ++nFlat; }
		if (s == 0) rho0 = r;

		// 보고: 엔진 힘의 유한차분 Jacobian(지금 넓이로 다시 잰 w 의 미분 포함) — ∂f/∂x = K + f ∇wᵀ / w 여야 한다
		double JFD[N * N], Kw[N * N], gw[N];
		for (int b = 0; b < N; ++b) {
			double fp[N], fm[N];
			engine_ref::ForceBridson(eb, idx + 1 + 2 * b, fp);
			engine_ref::ForceBridson(eb, idx + 2 + 2 * b, fm);
			const double h2 = H2[static_cast<std::size_t>(s) * N + b];
			for (int a = 0; a < N; ++a) JFD[a * N + b] = (fp[a] - fm[a]) / h2;
			gw[b] = (engine_ref::WeightBridson(eb, idx + 1 + 2 * b) - engine_ref::WeightBridson(eb, idx + 2 + 2 * b)) / h2;
		}
		double wt[N * N];
		for (int a = 0; a < N; ++a) for (int b = 0; b < N; ++b) { wt[a * N + b] = f[a] * gw[b] / w; Kw[a * N + b] = K[a * N + b] + wt[a * N + b]; }
		fdEng.Add(MaxAbsDiff(JFD, Je, N * N) / kScale);
		fdK.Add(MaxAbsDiff(JFD, K, N * N) / kScale);
		fdKw.Add(MaxAbsDiff(JFD, Kw, N * N) / kScale);
		wTerm.Add(MaxAbs(wt, N * N) / kScale);
	}
	engine_ref::DestroyBridson(eb);

	const bool fExact = fRel.Max() <= 1e-10;
	std::printf("B6 엔진 힘 — %s (분류, 기준 1e-10; 척도 = k w̄ max|∇θ|)\n", fExact ? "정확" : "다름");
	std::printf("  max|f_엔진 − (−kC∇C)| / 척도: 중앙 %.2e · max %.2e; 운동량 |Σf|/척도 max %.1e\n",
		fRel.Median(), fRel.Max(), fMom.Max());
	std::printf("  보고 — max|f| 대비: 중앙 %.2e, 평평한 표본 0 %.2e, 그 밖 max %.2e\n",
		fRelPlain.Median(), fPlain0, fRelPlainRest.Max());
	std::fprintf(js, "  \"B6_class\": \"%s\",\n", fExact ? "exact" : "differs");
	std::fprintf(js, "  \"B6_force\": {\"rel_scale_median\": %.3e, \"rel_scale_max\": %.3e, \"rel_plain_median\": %.3e, \"rel_plain_sample0\": %.3e, "
		"\"rel_plain_max_others\": %.3e, \"momentum_max\": %.3e},\n",
		fRel.Median(), fRel.Max(), fRelPlain.Median(), fPlain0, fRelPlainRest.Max(), fMom.Max());

	const bool jExact = jK.Max() <= 1e-10, jIsGN = jGN.Max() <= 1e-10;
	const char* jCls = jExact ? "정확" : (jIsGN ? "Gauss–Newton" : "어느 쪽도 아님");
	std::printf("B7 엔진 Jacobian — %s (분류, 기준 1e-10)\n", jCls);
	std::printf("  K 대비 rel: 중앙 %.2e · max %.2e; |θ| < 10° (%d 개) 중앙 %.2e · max %.2e; −k∇C∇Cᵀ 대비 중앙 %.2e\n",
		jK.Median(), jK.Max(), nFlat, jKflat.Median(), jKflat.Max(), jGN.Median());
	std::printf("  읽은 예측 대비 max %.2e (≤ 1e-10 이면 읽은 대로); 엔진 u = %s∇θ (%d / %d 표본), |u ∓ ∇θ|/max|∇θ| max %.1e\n",
		jP.Max(), uPlus * 2 >= nS ? "+" : "−", uPlus * 2 >= nS ? uPlus : nS - uPlus, nS, uSign.Max());
	std::printf("  오차 나눔(|K| 대비 중앙값 · max): ∇²θ 항 %.2e · %.2e, 계수 ½cos(θ/2)→1 %.2e · %.2e, 교차 항 %.2e · %.2e, (i, j) 옮김 %.2e · %.2e\n",
		bud[0].Median(), bud[0].Max(), bud[1].Median(), bud[1].Max(), bud[2].Median(), bud[2].Max(), bud[3].Median(), bud[3].Max());
	std::printf("  ∇θ 방향 강성 비 ∇θᵀJ∇θ / ∇θᵀK∇θ: 중앙 %.3f, |θ| < 10° 중앙 %.3f (범위 %.3f – %.3f), 평평한 표본 0 %.3f\n",
		rho.Median(), rhoFlat.Median(),
		rhoFlat.v.empty() ? 0. : *std::min_element(rhoFlat.v.begin(), rhoFlat.v.end()),
		rhoFlat.v.empty() ? 0. : *std::max_element(rhoFlat.v.begin(), rhoFlat.v.end()), rho0);
	std::printf("  엔진 힘의 유한차분 Jacobian(|K| 대비 max): 대 엔진 J %.2e, 대 K %.2e, 대 K + f∇wᵀ/w %.2e; w 미분 항 크기 중앙 %.2e · max %.2e\n",
		fdEng.Max(), fdK.Max(), fdKw.Max(), wTerm.Median(), wTerm.Max());
	std::fprintf(js, "  \"B7_class\": \"%s\",\n", jExact ? "exact" : (jIsGN ? "gauss_newton" : "neither"));
	std::fprintf(js, "  \"B7_jacobian\": {\"vs_K_median\": %.3e, \"vs_K_max\": %.3e, \"flat_hinges\": %d, \"flat_vs_K_median\": %.3e, \"flat_vs_K_max\": %.3e, "
		"\"vs_GN_median\": %.3e, \"vs_GN_max\": %.3e, \"vs_prediction_max\": %.3e, \"u_sign_plus\": %d, \"u_vs_grad_theta_max\": %.3e,\n",
		jK.Median(), jK.Max(), nFlat, jKflat.Median(), jKflat.Max(), jGN.Median(), jGN.Max(), jP.Max(), uPlus, uSign.Max());
	std::fprintf(js, "    \"budget_median\": {\"hess_theta_term\": %.3e, \"factor_half_cos\": %.3e, \"cross_terms\": %.3e, \"ij_block\": %.3e},\n",
		bud[0].Median(), bud[1].Median(), bud[2].Median(), bud[3].Median());
	std::fprintf(js, "    \"budget_max\": {\"hess_theta_term\": %.3e, \"factor_half_cos\": %.3e, \"cross_terms\": %.3e, \"ij_block\": %.3e},\n",
		bud[0].Max(), bud[1].Max(), bud[2].Max(), bud[3].Max());
	std::fprintf(js, "    \"stiffness_ratio_along_grad_theta\": {\"median\": %.4f, \"flat_median\": %.4f, \"flat_min\": %.4f, \"flat_max\": %.4f, \"sample0\": %.4f},\n",
		rho.Median(), rhoFlat.Median(),
		rhoFlat.v.empty() ? 0. : *std::min_element(rhoFlat.v.begin(), rhoFlat.v.end()),
		rhoFlat.v.empty() ? 0. : *std::max_element(rhoFlat.v.begin(), rhoFlat.v.end()), rho0);
	std::fprintf(js, "    \"engine_force_fd_jacobian\": {\"vs_engine_J_max\": %.3e, \"vs_K_max\": %.3e, \"vs_K_plus_w_term_max\": %.3e, \"w_term_median\": %.3e, \"w_term_max\": %.3e}},\n",
		fdEng.Max(), fdK.Max(), fdKw.Max(), wTerm.Median(), wTerm.Max());

	// ---- B8 시간 — 보고만, 이 실행 한 번 ------------------------------------------------
	std::printf("B8 시간 — 보고만, 이 실행 한 번 (모으기 포함, 단일 스레드, 5 회 중 최소; 판정은 timing_bases.py bridson)\n");
	const Timing t = TimeGrid(200, true);
	const Timing ts = TimeGrid(60, false);
	std::printf("  힌지 %d: 엔진 Oh(환원 없음) %.0f ns, 생성 EvalImplicit + 저장 %.0f ns → 비 %.3f\n",
		t.nHinge, t.engine, t.genImplicit, t.genImplicit / t.engine);
	std::printf("  보고 — 엔진 기본값(쌍 슬롯 %d 개로 환원) %.0f, 생성 EvalStaged + 저장 %.0f, EvalImplicit + 흩뿌리기 %.0f ns; 엔진 메모리 %.0f MB\n",
		t.nPair, t.engineSerial, t.genStaged, t.genScatter, t.engineBytes / 1048576.);
	std::printf("  보고 — 작은 격자 힌지 %d: 엔진 %.0f, 생성 %.0f ns → 비 %.3f\n", ts.nHinge, ts.engine, ts.genImplicit, ts.genImplicit / ts.engine);
	std::fprintf(js, "  \"B8_time\": {\"hinges\": %d, \"pairs\": %d, \"gen_over_engine\": %.3f, \"engine_bytes\": %.0f, \"small_hinges\": %d, \"small_gen_over_engine\": %.3f},\n",
		t.nHinge, t.nPair, t.genImplicit / t.engine, t.engineBytes, ts.nHinge, ts.genImplicit / ts.engine);
	std::fprintf(js, "  \"ns_per_hinge\": {\"engine_oh\": %.1f, \"engine_oh_serial\": %.1f, \"gen_implicit_store\": %.1f, \"gen_staged_store\": %.1f, "
		"\"gen_implicit_scatter\": %.1f, \"small_engine_oh\": %.1f, \"small_gen_implicit_store\": %.1f},\n",
		t.engine, t.engineSerial, t.genImplicit, t.genStaged, t.genScatter, ts.engine, ts.genImplicit);
#else
	std::printf("엔진 없이 빌드 — B5 – B8 건너뜀\n");
#endif

	std::fprintf(js, "  \"spec_gates_failed\": %d,\n  \"gates_failed\": %d\n}\n", specFail, rpt.fail);
	std::printf("\nBridson 명세 B1 – B4: %s (실패 %d)", specFail == 0 ? "통과" : "실패", specFail);
#ifdef SX_WITH_ENGINE
	std::printf(", 엔진 θ B5: %s", rpt.fail == specFail ? "통과" : "실패");
#endif
	std::printf(" — 옮김은 extract_engine.py, B8 은 timing_bases.py bridson\n");
	return rpt.fail == 0 ? 0 : 1;
}
