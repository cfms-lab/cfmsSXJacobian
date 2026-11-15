//---------------------------------------------------------------------------
// test_speed.cpp — 속도 게이트 S1 · S2 (docs/GATE_SPEED.md). S3(평탄판 불변)은 tests/check_flat_unchanged.py.
// S2 시간 판정은 tests/timing_bases.py(고정 로드 주소 8 곳의 최솟값). 이 exe 는 S1 을 판정하고 시간은 보고만 한다.
// 실행 디렉터리 = 리포 루트. results/gate_speed.json 을 쓴다.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffBend.h"
#include "../generated/BaraffBendCompact.h"
#include "gate_common.h"

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
static_assert(Compact::kDof == N && Compact::kHess == NH && Compact::kElemValues == Spec::kElemValues, "같은 힌지 배치");

// 변형 f 를 기준 ref 와 견준다: C 절대, ∇C · ∇²C 상대(힌지마다 max|Δ|/max|기준|)
static void CompareTo(EvalFn f, EvalFn ref, const std::vector<Hinge>& hs, Stat& sc, Stat& sg, Stat& sh)
{
	for (std::size_t s = 0; s < hs.size(); ++s) {
		double C, g[N], H[NH], Cr, gr[N], Hr[NH];
		f(hs[s].x, &hs[s].theta0, &C, g, H);
		ref(hs[s].x, &hs[s].theta0, &Cr, gr, Hr);
		sc.Add(std::fabs(C - Cr));
		sg.Add(Rel(MaxAbsDiff(g, gr, N), MaxAbs(gr, N)));
		sh.Add(Rel(MaxAbsDiff(H, Hr, NH), MaxAbs(Hr, NH)));
	}
}

// Mathematica 기준(results/mma_ref.txt 한 줄 = C − θ0, g, H 위 삼각)
static int CompareMma(EvalFn f, const std::vector<Hinge>& hs, const std::vector<std::vector<double> >& ref,
	double& eC, double& eG, double& eH)
{
	eC = eG = eH = 0.;
	int n = 0;
	for (std::size_t s = 0; s < ref.size() && s < hs.size(); ++s) {
		if (ref[s].size() != 1 + N + NH) continue;
		double C, g[N], H[NH];
		f(hs[s].x, &hs[s].theta0, &C, g, H);
		eC = std::max(eC, std::fabs(C - ref[s][0]));
		eG = std::max(eG, Rel(MaxAbsDiff(g, &ref[s][1], N), MaxAbs(&ref[s][1], N)));
		eH = std::max(eH, Rel(MaxAbsDiff(H, &ref[s][1 + N], NH), MaxAbs(&ref[s][1 + N], NH)));
		++n;
	}
	return n;
}

static void CompareFD(EvalFn f, const std::vector<Hinge>& hs, Stat& sg, Stat& sh)
{
	for (std::size_t s = 0; s < hs.size(); ++s) {
		const Hinge& h = hs[s];
		double C, g[N], H[NH], gFD[N], HFD[NH];
		f(h.x, &h.theta0, &C, g, H);
		FiniteDiff<N>(f, h.x, &h.theta0, 1e-6 * EdgeLength(h.x), gFD, HFD);
		sg.Add(Rel(MaxAbsDiff(g, gFD, N), MaxAbs(g, N)));
		sh.Add(Rel(MaxAbsDiff(H, HFD, NH), MaxAbs(H, NH)));
	}
}

// 격자 힌지 전체를 한 번 돈다. F 는 컴파일 시점 상수라 인라인될 수 있다(착수 게이트의 직접 호출과 같은 조건).
// assemble = 힘 · Jacobian 조립과 힘 흩뿌리기까지
template <void (*F)(const double*, const double*, double*, double*, double*)>
static double TimeLoop(const sx::StencilView<Spec>& view, bool assemble, std::vector<double>& Fext, double& checksum)
{
	typedef std::chrono::steady_clock Clock;
	std::fill(Fext.begin(), Fext.end(), 0.);
	const int nHinge = view.size();
	const Clock::time_point t0 = Clock::now();
	for (int s = 0; s < nHinge; ++s) {
		double x[N], v[N], ev[1], C, g[N], H[NH];
		view.Gather(s, x, v, ev);
		F(x, ev, &C, g, H);
		if (!assemble) { checksum += H[0]; continue; }
		sx::ConditionForce<N> out;
		sx::AssembleConditionForce<N>(1, &C, g, H, v, 1., 0.1, false, out);
		const int* nd = view.Nodes(s);
		for (int a = 0; a < 4; ++a) for (int c = 0; c < 3; ++c) Fext[3 * nd[a] + c] += out.f[3 * a + c] + out.d[3 * a + c];
		checksum += out.K[13] + out.Dx[7];
	}
	return std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / nHinge;
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Hinge> hs = ReadHinges("results/hinges.txt");
	if (hs.empty()) { std::printf("results/hinges.txt 없음 — model/make_hinges.py 를 먼저 돌린다\n"); return 2; }
	const std::vector<std::vector<double> > mma = ReadRows("results/mma_ref.txt");
	Report rep("results/gate_speed.json");
	std::FILE* js = rep.json;
	std::printf("속도 게이트 — 표본 %zu 힌지\n", hs.size());
	std::fprintf(js, "{\n  \"samples\": %zu,\n", hs.size());

	// ---- S1 단계판이 정확하다 -------------------------------------------------
	{
		double eC, eG, eH;
		const int n = CompareMma(&Spec::EvalStaged, hs, mma, eC, eG, eH);
		std::printf("S1a Mathematica 대조, EvalStaged (%d 힌지)\n", n);
		rep.Gate("C", n > 0 && eC <= 1e-12, "|ΔC| %.2e", eC, 1e-12);
		rep.Gate("∇C", n > 0 && eG <= 1e-10, "rel %.2e", eG, 1e-10);
		rep.Gate("∇²C", n > 0 && eH <= 1e-10, "rel %.2e", eH, 1e-10);
		std::fprintf(js, "  \"S1a_mathematica\": {\"hinges\": %d, \"C_abs\": %.3e, \"grad_rel\": %.3e, \"hess_rel\": %.3e},\n", n, eC, eG, eH);

		Stat fg, fh;
		CompareFD(&Spec::EvalStaged, hs, fg, fh);
		std::printf("S1b 유한차분, EvalStaged (%zu 힌지)\n", hs.size());
		rep.Gate("∇C", fg.Max() <= 1e-6, "rel max %.2e", fg.Max(), 1e-6);
		rep.Gate("∇²C", fh.Max() <= 1e-6, "rel max %.2e", fh.Max(), 1e-6);
		std::fprintf(js, "  \"S1b_finite_difference\": {\"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e},\n", fg.Max(), fh.Max());

		Stat sc, sg, sh;
		CompareTo(&Spec::EvalStaged, &Spec::Eval, hs, sc, sg, sh);
		std::printf("S1c 평탄판 Eval 과 비교, EvalStaged (%zu 힌지)\n", hs.size());
		rep.Gate("C", sc.Max() <= 1e-12, "|ΔC| max %.2e", sc.Max(), 1e-12);
		rep.Gate("∇C", sg.Max() <= 1e-12, "rel max %.2e", sg.Max(), 1e-12);
		rep.Gate("∇²C", sh.Max() <= 1e-12, "rel max %.2e", sh.Max(), 1e-12);
		std::fprintf(js, "  \"S1c_vs_flat\": {\"C_abs_max\": %.3e, \"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e, \"hess_rel_median\": %.3e},\n",
			sc.Max(), sg.Max(), sh.Max(), sh.Median());
	}

	// ---- 보고: 압축 명세 -------------------------------------------------------
	{
		double eC, eG, eH;
		const int n = CompareMma(&Compact::EvalStaged, hs, mma, eC, eG, eH);
		Stat sc, sg, sh, fg, fh;
		CompareTo(&Compact::EvalStaged, &Spec::Eval, hs, sc, sg, sh);
		CompareFD(&Compact::EvalStaged, hs, fg, fh);
		std::printf("  보고 — 압축 명세 EvalStaged: Mathematica(%d) |ΔC| %.1e · ∇C %.1e · ∇²C %.1e, 평탄 Eval 과 ∇²C rel max %.1e, 유한차분 ∇²C %.1e\n",
			n, eC, eG, eH, sh.Max(), fh.Max());
		std::fprintf(js, "  \"compact\": {\"mma_C_abs\": %.3e, \"mma_grad_rel\": %.3e, \"mma_hess_rel\": %.3e, \"vs_flat_C_abs_max\": %.3e, \"vs_flat_grad_rel_max\": %.3e, \"vs_flat_hess_rel_max\": %.3e, \"fd_grad_rel_max\": %.3e, \"fd_hess_rel_max\": %.3e},\n",
			eC, eG, eH, sc.Max(), sg.Max(), sh.Max(), fg.Max(), fh.Max());
	}

	// ---- S2 시간 ---------------------------------------------------------------
	{
		std::mt19937_64 rng(7);
		const Grid grid = MakeGrid(200, rng);
		std::uniform_real_distribution<double> noise(-1., 1.);
		std::vector<double> T0(grid.nHinge);
		for (int s = 0; s < grid.nHinge; ++s) T0[s] = 0.01 * noise(rng);
		SXMesh mesh;
		FillMesh<Spec>(mesh, grid, T0);
		const sx::StencilView<Spec> view(mesh);

		enum { kFlat, kFrozen, kStaged, kCompactStaged, kCompactFlat, kStagedAsm, kFrozenAsm, kFlatAsm, kVariants };
		static const char* names[kVariants] = { "eval_flat", "eval_frozen", "eval_staged", "compact_staged", "compact_flat",
			"staged_assemble", "frozen_assemble", "flat_assemble" };
		double best[kVariants];
		for (int i = 0; i < kVariants; ++i) best[i] = 1e30;
		std::vector<double> Fext(3 * grid.nNode);
		double checksum = 0.;
		for (int rep5 = 0; rep5 < 5; ++rep5) {   // 변형을 번갈아 5 회, 변형마다 최소
			best[kFlat]          = std::min(best[kFlat],          TimeLoop<&Spec::Eval>(view, false, Fext, checksum));
			best[kFrozen]        = std::min(best[kFrozen],        TimeLoop<&Spec::EvalFrozen>(view, false, Fext, checksum));
			best[kStaged]        = std::min(best[kStaged],        TimeLoop<&Spec::EvalStaged>(view, false, Fext, checksum));
			best[kCompactStaged] = std::min(best[kCompactStaged], TimeLoop<&Compact::EvalStaged>(view, false, Fext, checksum));
			best[kCompactFlat]   = std::min(best[kCompactFlat],   TimeLoop<&Compact::Eval>(view, false, Fext, checksum));
			best[kStagedAsm]     = std::min(best[kStagedAsm],     TimeLoop<&Spec::EvalStaged>(view, true, Fext, checksum));
			best[kFrozenAsm]     = std::min(best[kFrozenAsm],     TimeLoop<&Spec::EvalFrozen>(view, true, Fext, checksum));
			best[kFlatAsm]       = std::min(best[kFlatAsm],       TimeLoop<&Spec::Eval>(view, true, Fext, checksum));
		}
		const double ratio = best[kStaged] / best[kFrozen];
		// S2 판정은 tests/timing_bases.py 가 고정 로드 주소 8 곳의 최솟값으로 한다. 여기 값은 이 실행(로드 주소 하나)뿐이다.
		std::printf("S2 시간 — 보고만, 이 실행 한 번 (힌지 %d, 모으기 포함, 단일 스레드, 5 회 중 최소)\n", grid.nHinge);
		std::printf("  보고 — EvalStaged / EvalFrozen = %.2f (판정은 timing_bases.py)\n", ratio);
		std::printf("  보고 — 힌지당 ns: 평탄 정확 %.0f, 길이 고정 %.0f, 단계 %.0f (평탄 대비 %.1f 배 빠름), 압축 단계 %.0f, 압축 평탄 %.0f\n",
			best[kFlat], best[kFrozen], best[kStaged], best[kFlat] / best[kStaged], best[kCompactStaged], best[kCompactFlat]);
		std::printf("  보고 — 조립 포함 ns: 단계 %.0f, 길이 고정 %.0f, 평탄 정확 %.0f  (checksum %.3e)\n",
			best[kStagedAsm], best[kFrozenAsm], best[kFlatAsm], checksum);
		std::fprintf(js, "  \"S2_time\": {\"hinges\": %d, \"staged_over_frozen\": %.3f, \"flat_over_staged\": %.3f},\n",
			grid.nHinge, ratio, best[kFlat] / best[kStaged]);
		std::fprintf(js, "  \"ns_per_hinge\": {");
		for (int i = 0; i < kVariants; ++i) std::fprintf(js, "%s\"%s\": %.1f", i ? ", " : "", names[i], best[i]);
		std::fprintf(js, "},\n");
	}

	std::fprintf(js, "  \"gates_failed\": %d\n}\n", rep.fail);
	std::printf("\n속도 게이트 S1: %s (실패 %d) — S2 는 timing_bases.py\n", rep.fail == 0 ? "GO" : "NO-GO", rep.fail);
	return rep.fail == 0 ? 0 : 1;
}
