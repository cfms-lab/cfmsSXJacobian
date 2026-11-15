//---------------------------------------------------------------------------
// test_engine.cpp — 엔진 손 코드 비교 게이트 E2 · E3(docs/GATE_ENGINE.md). E1 은 tests/engine_ref/extract_engine.py,
// E4 시간 판정은 tests/timing_bases.py(고정 로드 주소 8 곳의 최솟값). 이 exe 는 E2 를 판정하고 E3 를 분류하며,
// 시간은 이 실행 한 번의 값으로 보고만 한다.
// 실행 디렉터리 = 리포 루트. results/gate_engine.json 을 쓴다.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffBend.h"
#include "../generated/BaraffBendCompact.h"
#include "gate_common.h"
#include "engine_ref/engine_macri.h"

#include <chrono>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;

typedef sxgen::BaraffBend Spec;
typedef sxgen::BaraffBendCompact Compact;
static const int N = Spec::kDof;      // 12
static const int NH = Spec::kHess;    // 78
static const int NO = 1 + N + NH;     // 저장하는 출력 91 개
static const double kPrm[3] = { 1., 0.1, 0.01 };   // EvalImplicit 보고용 {k, kd, h} — 조립 게이트와 같다

// 위 삼각 78 → 12 × 12 전체(대칭)
static void Full(const double* Hu, double* Hf)
{
	for (int a = 0; a < N; ++a)
		for (int b = 0; b < N; ++b) Hf[a * N + b] = Hu[a <= b ? sx::UpperIndex<N>(a, b) : sx::UpperIndex<N>(b, a)];
}

static void Negate(double* p, int n) { for (int i = 0; i < n; ++i) p[i] = -p[i]; }

// ---- 시간 ----------------------------------------------------------------------
typedef std::chrono::steady_clock Clock;

static double TimeEngine(engine_ref::Macri* em, double& checksum)
{
	const Clock::time_point t0 = Clock::now();
	engine_ref::Run(em);
	const double ns = std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / engine_ref::Size(em);
	checksum += engine_ref::Checksum(em);
	return ns;
}

// 모으기 + 평가. store = C · ∇C · ∇²C(91 개)를 힌지별 출력 배열에 저장(E4 의 생성 쪽), 아니면 checksum 만
template <void (*F)(const double*, const double*, double*, double*, double*)>
static double TimeGen(const sx::StencilView<Spec>& view, bool store, std::vector<double>& out, double& checksum)
{
	const int nHinge = view.size();
	const Clock::time_point t0 = Clock::now();
	for (int s = 0; s < nHinge; ++s) {
		double x[N], v[N], ev[1];
		view.Gather(s, x, v, ev);
		if (store) {
			double* o = &out[static_cast<std::size_t>(s) * NO];
			F(x, ev, o, o + 1, o + 1 + N);
		} else {
			double C, g[N], H[NH];
			F(x, ev, &C, g, H);
			checksum += H[0];
		}
	}
	const double ns = std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / nHinge;
	if (store) for (int s = 0; s < nHinge; s += 97) checksum += out[static_cast<std::size_t>(s) * NO + 13];
	return ns;
}

template <void (*F)(const double*, const double*, const double*, const double*, double*, double*, double*, double*)>
static double TimeImplicit(const sx::StencilView<Spec>& view, std::vector<double>& Fn, std::vector<double>& Bn, double& checksum)
{
	std::fill(Fn.begin(), Fn.end(), 0.);
	std::fill(Bn.begin(), Bn.end(), 0.);
	const int nHinge = view.size();
	const Clock::time_point t0 = Clock::now();
	for (int s = 0; s < nHinge; ++s) {
		double x[N], v[N], ev[1], C, f[N], A[NH], b[N];
		view.Gather(s, x, v, ev);
		F(x, ev, v, kPrm, &C, f, A, b);
		const int* nd = view.Nodes(s);
		for (int a = 0; a < 4; ++a) for (int c = 0; c < 3; ++c) { Fn[3 * nd[a] + c] += f[3 * a + c]; Bn[3 * nd[a] + c] += b[3 * a + c]; }
		checksum += A[13];
	}
	return std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / nHinge;
}

struct Timing {
	double engine, genStore, gen, frozen, compactStore, implicit;
	int nHinge;
	double engineBytes;
};

static Timing TimeGrid(int G, bool full)
{
	std::mt19937_64 rng(7);
	const Grid grid = MakeGrid(G, rng);
	std::uniform_real_distribution<double> noise(-1., 1.);
	std::vector<double> T0(grid.nHinge);
	for (int s = 0; s < grid.nHinge; ++s) T0[s] = 0.01 * noise(rng);
	SXMesh mesh;
	FillMesh<Spec>(mesh, grid, T0);
	const sx::StencilView<Spec> view(mesh);

	std::vector<double> X(grid.X), V(grid.V);   // 엔진은 노드 배열을 가리키기만 한다
	engine_ref::Macri* em = engine_ref::Create(grid.nHinge, &grid.hinges[0], grid.nNode, &X[0], &V[0]);

	std::vector<double> out(static_cast<std::size_t>(grid.nHinge) * NO, 0.);
	std::vector<double> Fn(3 * grid.nNode), Bn(3 * grid.nNode);
	Timing t;
	t.engine = t.genStore = t.gen = t.frozen = t.compactStore = t.implicit = 1e30;
	t.nHinge = grid.nHinge;
	t.engineBytes = engine_ref::Bytes(em);
	double checksum = 0.;
	for (int rep5 = 0; rep5 < 5; ++rep5) {   // 변형을 번갈아 5 회, 변형마다 최소
		t.engine   = std::min(t.engine,   TimeEngine(em, checksum));
		t.genStore = std::min(t.genStore, TimeGen<&Spec::EvalStaged>(view, true, out, checksum));
		if (!full) continue;
		t.gen          = std::min(t.gen,          TimeGen<&Spec::EvalStaged>(view, false, out, checksum));
		t.frozen       = std::min(t.frozen,       TimeGen<&Spec::EvalFrozen>(view, false, out, checksum));
		t.compactStore = std::min(t.compactStore, TimeGen<&Compact::EvalStaged>(view, true, out, checksum));
		t.implicit     = std::min(t.implicit,     TimeImplicit<&Spec::EvalImplicit>(view, Fn, Bn, checksum));
	}
	std::printf("  (격자 G = %d, checksum %.3e)\n", G, checksum);
	engine_ref::Destroy(em);
	return t;
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Hinge> hs = ReadHinges("results/hinges.txt");
	const std::vector<std::vector<double> > ref = ReadRows("results/mma_ref.txt");
	if (hs.empty() || ref.empty()) { std::printf("results/hinges.txt · mma_ref.txt 없음\n"); return 2; }
	const int nS = static_cast<int>(hs.size());
	Report rpt("results/gate_engine.json");
	std::FILE* js = rpt.json;
	std::printf("엔진 손 코드 비교 게이트 — 표본 %d 힌지, Mathematica 기준 %zu 개\n", nS, ref.size());
	std::fprintf(js, "{\n  \"engine_commit\": \"378d20c\",\n  \"samples\": %d,\n", nS);

	// ---- 엔진 손 코드를 표본에 한 번 돌린다: 힌지마다 자기 노드 4 개 --------------------
	std::vector<double> X(12 * nS), V(12 * nS);
	std::vector<int> vid(4 * nS);
	for (int s = 0; s < nS; ++s) {
		for (int i = 0; i < 12; ++i) { X[12 * s + i] = hs[s].x[i]; V[12 * s + i] = hs[s].v[i]; }
		for (int m = 0; m < 4; ++m) vid[4 * s + m] = 4 * s + m;
	}
	engine_ref::Macri* em = engine_ref::Create(nS, &vid[0], 4 * nS, &X[0], &V[0]);
	engine_ref::Run(em);

	// ---- E2 · E3 ------------------------------------------------------------------
	Stat e2, e2all;                       // |θ_엔진 − (−θ)|
	Stat gM, hM, gF, hF;                  // Mathematica · EvalFrozen 대비(200)
	Stat gMs[2], hMs[2];                  // 가지별(0 = |cos| ≥ |sin|, 1 = |sin| > |cos|)
	Stat gE, hE;                          // Eval 대비(1000, 보고)
	Stat asym, ratioG[2];                 // 엔진 ∇²C 비대칭, max|∇C_엔진| / max|∇C| (가지별)
	int nRef = 0, nBranch[2] = { 0, 0 };
	for (int s = 0; s < nS; ++s) {
		const Hinge& h = hs[s];
		const double th = engine_ref::Theta(em, s);
		double ge[N], He[N * N];
		engine_ref::Grad(em, s, ge);
		engine_ref::Hess(em, s, He);
		Negate(ge, N); Negate(He, N * N);     // θ_엔진 = −θ (판정문 「꼭짓점 · 부호 대응」)
		const int br = engine_ref::SinBigger(em, s) ? 1 : 0;

		double C, g[N], Hu[NH], Hf[N * N];
		Spec::Eval(h.x, &h.theta0, &C, g, Hu);
		Full(Hu, Hf);
		e2all.Add(std::fabs(th + (C + h.theta0)));
		gE.Add(Rel(MaxAbsDiff(ge, g, N), MaxAbs(g, N)));
		hE.Add(Rel(MaxAbsDiff(He, Hf, N * N), MaxAbs(Hf, N * N)));
		ratioG[br].Add(MaxAbs(ge, N) / MaxAbs(g, N));
		double a = 0.;
		for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) a = std::max(a, std::fabs(He[i * N + j] - He[j * N + i]));
		asym.Add(Rel(a, MaxAbs(He, N * N)));

		if (s < static_cast<int>(ref.size()) && ref[s].size() == static_cast<std::size_t>(NO)) {
			const std::vector<double>& r = ref[s];
			double Hr[N * N];
			Full(&r[1 + N], Hr);
			e2.Add(std::fabs(th + (r[0] + h.theta0)));
			const double rg = Rel(MaxAbsDiff(ge, &r[1], N), MaxAbs(&r[1], N));
			const double rh = Rel(MaxAbsDiff(He, Hr, N * N), MaxAbs(Hr, N * N));
			gM.Add(rg); hM.Add(rh);
			gMs[br].Add(rg); hMs[br].Add(rh);
			++nBranch[br];
			double Cf, gf[N], Hfu[NH], Hff[N * N];
			Spec::EvalFrozen(h.x, &h.theta0, &Cf, gf, Hfu);
			Full(Hfu, Hff);
			gF.Add(Rel(MaxAbsDiff(ge, gf, N), MaxAbs(gf, N)));
			hF.Add(Rel(MaxAbsDiff(He, Hff, N * N), MaxAbs(Hff, N * N)));
			++nRef;
		}
	}
	engine_ref::Destroy(em);

	std::printf("E2 θ — Mathematica %d 힌지\n", nRef);
	rpt.Gate("E2", nRef > 0 && e2.Max() <= 1e-12, "max |θ_엔진 + θ| %.2e rad", e2.Max(), 1e-12);
	std::printf("  보고 — 1000 힌지 Eval 대비 max %.2e rad\n", e2all.Max());

	const bool exact = gM.Max() <= 1e-10 && hM.Max() <= 1e-10;
	const bool frozen = gF.Max() <= 1e-10 && hF.Max() <= 1e-10;
	const char* cls = exact ? "정확" : (frozen ? "길이 고정 근사" : "어느 쪽도 아님");
	std::printf("E3 분류 — %s (게이트 아님, 기준 1e-10)\n", cls);
	std::printf("  Mathematica 대비: ∇C rel 중앙 %.2e · max %.2e, ∇²C rel 중앙 %.2e · max %.2e\n",
		gM.Median(), gM.Max(), hM.Median(), hM.Max());
	std::printf("  EvalFrozen 대비:  ∇C rel 중앙 %.2e · max %.2e, ∇²C rel 중앙 %.2e · max %.2e\n",
		gF.Median(), gF.Max(), hF.Median(), hF.Max());
	static const char* brName[2] = { "|cos| ≥ |sin|", "|sin| > |cos|" };
	for (int b = 0; b < 2; ++b)
		std::printf("  가지 %s (%d 개): ∇C rel 중앙 %.2e · max %.2e, ∇²C rel 중앙 %.2e · max %.2e\n",
			brName[b], nBranch[b], gMs[b].Median(), gMs[b].Max(), hMs[b].Median(), hMs[b].Max());
	std::printf("  보고 — 1000 힌지 Eval 대비: ∇C rel max %.2e, ∇²C rel max %.2e; 엔진 ∇²C 비대칭 rel max %.2e\n",
		gE.Max(), hE.Max(), asym.Max());
	std::printf("  보고 — 크기 비 max|∇C_엔진| / max|∇C| 중앙값: %s %.3e (%zu 개), %s %.3e (%zu 개)\n",
		brName[0], ratioG[0].Median(), ratioG[0].v.size(), brName[1], ratioG[1].Median(), ratioG[1].v.size());
	std::fprintf(js, "  \"E2_theta\": {\"hinges\": %d, \"abs_max\": %.3e, \"abs_max_1000_vs_eval\": %.3e},\n", nRef, e2.Max(), e2all.Max());
	std::fprintf(js, "  \"E3_class\": \"%s\",\n", exact ? "exact" : (frozen ? "length_frozen" : "neither"));
	std::fprintf(js, "  \"E3_vs_mathematica\": {\"grad_rel_median\": %.3e, \"grad_rel_max\": %.3e, \"hess_rel_median\": %.3e, \"hess_rel_max\": %.3e},\n",
		gM.Median(), gM.Max(), hM.Median(), hM.Max());
	std::fprintf(js, "  \"E3_vs_frozen\": {\"grad_rel_median\": %.3e, \"grad_rel_max\": %.3e, \"hess_rel_median\": %.3e, \"hess_rel_max\": %.3e},\n",
		gF.Median(), gF.Max(), hF.Median(), hF.Max());
	std::fprintf(js, "  \"E3_by_branch\": {\"cos_ge_sin\": {\"hinges\": %d, \"grad_rel_median\": %.3e, \"grad_rel_max\": %.3e, \"hess_rel_median\": %.3e, \"hess_rel_max\": %.3e}, "
		"\"sin_gt_cos\": {\"hinges\": %d, \"grad_rel_median\": %.3e, \"grad_rel_max\": %.3e, \"hess_rel_median\": %.3e, \"hess_rel_max\": %.3e}},\n",
		nBranch[0], gMs[0].Median(), gMs[0].Max(), hMs[0].Median(), hMs[0].Max(),
		nBranch[1], gMs[1].Median(), gMs[1].Max(), hMs[1].Median(), hMs[1].Max());
	std::fprintf(js, "  \"report_1000_vs_eval\": {\"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e, \"engine_hess_asym_rel_max\": %.3e, "
		"\"grad_scale_median_cos_ge_sin\": %.3e, \"grad_scale_median_sin_gt_cos\": %.3e},\n",
		gE.Max(), hE.Max(), asym.Max(), ratioG[0].Median(), ratioG[1].Median());

	// ---- E4 시간 — 보고만, 이 실행 한 번 ------------------------------------------------
	std::printf("E4 시간 — 보고만, 이 실행 한 번 (모으기 포함, 단일 스레드, 5 회 중 최소; 판정은 timing_bases.py)\n");
	const Timing t = TimeGrid(200, true);
	const Timing ts = TimeGrid(60, false);
	std::printf("  힌지 %d: 엔진 Macri %.0f ns, 생성 EvalStaged + 저장 %.0f ns → 비 %.3f\n",
		t.nHinge, t.engine, t.genStore, t.genStore / t.engine);
	std::printf("  보고 — 저장 없음 %.0f, EvalFrozen %.0f, 압축 + 저장 %.0f, EvalImplicit(조립까지) %.0f ns; 엔진 sbv %.0f MB (%.0f B/힌지)\n",
		t.gen, t.frozen, t.compactStore, t.implicit, t.engineBytes / 1048576., t.engineBytes / t.nHinge);
	std::printf("  보고 — 작은 격자 힌지 %d: 엔진 %.0f, 생성 + 저장 %.0f ns → 비 %.3f\n",
		ts.nHinge, ts.engine, ts.genStore, ts.genStore / ts.engine);
	std::fprintf(js, "  \"E4_time\": {\"hinges\": %d, \"gen_store_over_engine\": %.3f, \"engine_sbv_bytes\": %.0f, \"small_hinges\": %d, \"small_gen_store_over_engine\": %.3f},\n",
		t.nHinge, t.genStore / t.engine, t.engineBytes, ts.nHinge, ts.genStore / ts.engine);
	std::fprintf(js, "  \"ns_per_hinge\": {\"engine_macri\": %.1f, \"gen_staged_store\": %.1f, \"gen_staged\": %.1f, \"gen_frozen\": %.1f, "
		"\"compact_staged_store\": %.1f, \"gen_implicit\": %.1f, \"small_engine_macri\": %.1f, \"small_gen_staged_store\": %.1f},\n",
		t.engine, t.genStore, t.gen, t.frozen, t.compactStore, t.implicit, ts.engine, ts.genStore);

	std::fprintf(js, "  \"gates_failed\": %d\n}\n", rpt.fail);
	std::printf("\n엔진 비교 E2: %s (실패 %d) — E1 은 extract_engine.py, E4 는 timing_bases.py\n", rpt.fail == 0 ? "통과" : "실패", rpt.fail);
	return rpt.fail == 0 ? 0 : 1;
}
