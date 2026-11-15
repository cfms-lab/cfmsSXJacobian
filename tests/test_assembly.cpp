//---------------------------------------------------------------------------
// test_assembly.cpp — 조립 게이트 (docs/GATE_ASSEMBLY.md).
//   2 차(판정): R1′ EvalImplicit 정확도 · R2 기존 K · R3′ 시간(판정은 tests/timing_bases.py, 이 exe 는 보고만)
//   1 차(보고만): R1 AssembleImplicit 정확도 · R3 시간 비
// 실행 디렉터리 = 리포 루트. results/gate_assembly.json 을 쓴다.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffBend.h"
#include "../generated/BaraffBendCompact.h"
#include "gate_common.h"

#include <chrono>
#include <type_traits>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;

typedef sxgen::BaraffBend Spec;
typedef sxgen::BaraffBendCompact Compact;
static const int N = Spec::kDof;      // 12
static const int NH = Spec::kHess;    // 78

static const double kK = 1.0, kKd = 0.1, kH = 0.01;   // 판정문에 고정한 계수
static const double kPrm[3] = { kK, kKd, kH };        // EvalImplicit 의 prm = {k, kd, h}

// 기존 조립 결과에서 식 (6) 의 요소 기여를 만든다(기준): A = −h Dv − h²(K + Dx), b = h(f + d) + h²(K + Dx) v
static void ReferenceFromDense(const sx::ConditionForce<N>& o, const double* v, double h, double* Aup, double* b, double* f)
{
	for (int a = 0; a < N; ++a) {
		double Jv = 0.;
		for (int j = 0; j < N; ++j) Jv += (o.K[a * N + j] + o.Dx[a * N + j]) * v[j];
		f[a] = o.f[a] + o.d[a];
		b[a] = h * f[a] + h * h * Jv;
		for (int j = a; j < N; ++j) Aup[sx::UpperIndex<N>(a, j)] = -h * o.Dv[a * N + j] - h * h * (o.K[a * N + j] + o.Dx[a * N + j]);
	}
}

// 탄성력 f = −k C ∇C
static void SpringForce(const double* x, const double* ev, double* f)
{
	double C, g[N], H[NH];
	Spec::EvalStaged(x, ev, &C, g, H);
	for (int a = 0; a < N; ++a) f[a] = -kK * C * g[a];
}

// 힌지 하나의 평가 뒤 단계. 0 = 평가만, 1 = 기존 조립 + f+d 흩뿌리기, 2 = 1 차 새 조립 + b · f 흩뿌리기
typedef std::integral_constant<int, 0> EvalOnly;
typedef std::integral_constant<int, 1> OldAssembly;
typedef std::integral_constant<int, 2> NewAssembly;

static inline void Step(EvalOnly, const int*, double, const double*, const double* H, const double*,
	std::vector<double>&, std::vector<double>&, double& checksum)
{
	checksum += H[0];
}
static inline void Step(OldAssembly, const int* nd, double C, const double* g, const double* H, const double* v,
	std::vector<double>& Fn, std::vector<double>&, double& checksum)
{
	sx::ConditionForce<N> out;
	sx::AssembleConditionForce<N>(1, &C, g, H, v, kK, kKd, false, out);
	for (int a = 0; a < 4; ++a) for (int c = 0; c < 3; ++c) Fn[3 * nd[a] + c] += out.f[3 * a + c] + out.d[3 * a + c];
	checksum += out.K[13] + out.Dx[7] + out.Dv[5];
}
static inline void Step(NewAssembly, const int* nd, double C, const double* g, const double* H, const double* v,
	std::vector<double>& Fn, std::vector<double>& Bn, double& checksum)
{
	sx::ImplicitBlock<N> out;
	sx::AssembleImplicit<N>(1, &C, g, H, v, kK, kKd, kH, out);
	for (int a = 0; a < 4; ++a) for (int c = 0; c < 3; ++c) { Fn[3 * nd[a] + c] += out.f[3 * a + c]; Bn[3 * nd[a] + c] += out.b[3 * a + c]; }
	checksum += out.A[13];
}

// 격자 힌지 전체를 한 번 돈다. F · ASM 은 컴파일 시점 상수라 인라인될 수 있다
template <void (*F)(const double*, const double*, double*, double*, double*), int ASM>
static double TimeLoop(const sx::StencilView<Spec>& view, std::vector<double>& Fn, std::vector<double>& Bn, double& checksum)
{
	typedef std::chrono::steady_clock Clock;
	std::fill(Fn.begin(), Fn.end(), 0.);
	std::fill(Bn.begin(), Bn.end(), 0.);
	const int nHinge = view.size();
	const Clock::time_point t0 = Clock::now();
	for (int s = 0; s < nHinge; ++s) {
		double x[N], v[N], ev[1], C, g[N], H[NH];
		view.Gather(s, x, v, ev);
		F(x, ev, &C, g, H);
		Step(std::integral_constant<int, ASM>(), view.Nodes(s), C, g, H, v, Fn, Bn, checksum);
	}
	return std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / nHinge;
}

// 2 차: 생성된 EvalImplicit(평가 + 조립) + b · f 흩뿌리기. 1 차 NewAssembly 와 같은 일을 한다
typedef void (*ImplicitFn)(const double*, const double*, const double*, const double*, double*, double*, double*, double*);
template <ImplicitFn F>
static double TimeLoopImplicit(const sx::StencilView<Spec>& view, std::vector<double>& Fn, std::vector<double>& Bn, double& checksum)
{
	typedef std::chrono::steady_clock Clock;
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

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Hinge> hs = ReadHinges("results/hinges.txt");
	if (hs.empty()) { std::printf("results/hinges.txt 없음 — model/make_hinges.py 를 먼저 돌린다\n"); return 2; }
	Report rep("results/gate_assembly.json");
	std::FILE* js = rep.json;
	std::printf("조립 게이트 — 표본 %zu 힌지, k = %g, kd = %g, h = %g\n", hs.size(), kK, kKd, kH);
	std::fprintf(js, "{\n  \"samples\": %zu, \"k\": %g, \"kd\": %g, \"h\": %g,\n", hs.size(), kK, kKd, kH);

	// ---- R1′ (2 차, 판정) EvalImplicit = 기존 조립에서 만든 식 (6) 기여 / R1 (1 차, 보고) -----
	{
		Stat eA, eB, eF, eC, nA, nB, nF, cA, cB, cF, asym;
		for (std::size_t s = 0; s < hs.size(); ++s) {
			const Hinge& hg = hs[s];
			double C, g[N], H[NH];
			Spec::EvalStaged(hg.x, &hg.theta0, &C, g, H);
			sx::ConditionForce<N> dense;
			sx::AssembleConditionForce<N>(1, &C, g, H, hg.v, kK, kKd, false, dense);
			double Aref[NH], bref[N], fref[N];
			ReferenceFromDense(dense, hg.v, kH, Aref, bref, fref);
			const double mA = MaxAbs(Aref, NH), mB = MaxAbs(bref, N), mF = MaxAbs(fref, N);

			double Ci, fi[N], Ai[NH], bi[N];
			Spec::EvalImplicit(hg.x, &hg.theta0, hg.v, kPrm, &Ci, fi, Ai, bi);
			eA.Add(Rel(MaxAbsDiff(Ai, Aref, NH), mA));
			eB.Add(Rel(MaxAbsDiff(bi, bref, N), mB));
			eF.Add(Rel(MaxAbsDiff(fi, fref, N), mF));
			eC.Add(std::fabs(Ci - C));

			sx::ImplicitBlock<N> blk;   // 1 차
			sx::AssembleImplicit<N>(1, &C, g, H, hg.v, kK, kKd, kH, blk);
			nA.Add(Rel(MaxAbsDiff(blk.A, Aref, NH), mA));
			nB.Add(Rel(MaxAbsDiff(blk.b, bref, N), mB));
			nF.Add(Rel(MaxAbsDiff(blk.f, fref, N), mF));

			double Cc, fc[N], Ac[NH], bc[N];   // 압축 명세(보고)
			Compact::EvalImplicit(hg.x, &hg.theta0, hg.v, kPrm, &Cc, fc, Ac, bc);
			cA.Add(Rel(MaxAbsDiff(Ac, Aref, NH), mA));
			cB.Add(Rel(MaxAbsDiff(bc, bref, N), mB));
			cF.Add(Rel(MaxAbsDiff(fc, fref, N), mF));

			// 기준 행렬이 정말 대칭인가(원문 꼴의 감쇠라면 그래야 한다)
			double m = 0., ma = 0.;
			for (int a = 0; a < N; ++a) for (int j = 0; j < N; ++j) {
				const double Aaj = -kH * dense.Dv[a * N + j] - kH * kH * (dense.K[a * N + j] + dense.Dx[a * N + j]);
				const double Aja = -kH * dense.Dv[j * N + a] - kH * kH * (dense.K[j * N + a] + dense.Dx[j * N + a]);
				m = std::max(m, std::fabs(Aaj - Aja)); ma = std::max(ma, std::fabs(Aaj));
			}
			asym.Add(Rel(m, ma));
		}
		std::printf("R1′ EvalImplicit 대 기존 조립 (%zu 힌지)\n", hs.size());
		rep.Gate("A", eA.Max() <= 1e-13, "rel max %.2e", eA.Max(), 1e-13);
		rep.Gate("b", eB.Max() <= 1e-13, "rel max %.2e", eB.Max(), 1e-13);
		rep.Gate("f+d", eF.Max() <= 1e-13, "rel max %.2e", eF.Max(), 1e-13);
		std::printf("  보고 — C 차 max %.2e (EvalStaged 와)\n", eC.Max());
		std::printf("  보고 — 1 차 AssembleImplicit: A %.2e, b %.2e, f+d %.2e\n", nA.Max(), nB.Max(), nF.Max());
		std::printf("  보고 — 압축 명세 EvalImplicit: A %.2e, b %.2e, f+d %.2e\n", cA.Max(), cB.Max(), cF.Max());
		std::printf("  보고 — 기준 행렬의 비대칭 max|A − Aᵀ|/max|A|: %.2e\n", asym.Max());
		std::fprintf(js, "  \"R1p_implicit_vs_dense\": {\"A_rel_max\": %.3e, \"b_rel_max\": %.3e, \"f_rel_max\": %.3e, \"C_abs_max\": %.3e},\n",
			eA.Max(), eB.Max(), eF.Max(), eC.Max());
		std::fprintf(js, "  \"R1_first_attempt_vs_dense\": {\"A_rel_max\": %.3e, \"b_rel_max\": %.3e, \"f_rel_max\": %.3e, \"ref_asym_rel_max\": %.3e},\n",
			nA.Max(), nB.Max(), nF.Max(), asym.Max());
		std::fprintf(js, "  \"compact_implicit_vs_dense\": {\"A_rel_max\": %.3e, \"b_rel_max\": %.3e, \"f_rel_max\": %.3e},\n",
			cA.Max(), cB.Max(), cF.Max());
	}

	// ---- R2 기존 K 대 탄성력의 유한차분 ----------------------------------------
	{
		Stat eK;
		for (std::size_t s = 0; s < hs.size(); ++s) {
			const Hinge& hg = hs[s];
			double C, g[N], H[NH];
			Spec::EvalStaged(hg.x, &hg.theta0, &C, g, H);
			sx::ConditionForce<N> dense;
			sx::AssembleConditionForce<N>(1, &C, g, H, hg.v, kK, kKd, false, dense);
			const double step = 1e-6 * EdgeLength(hg.x);
			double KFD[N * N];
			for (int j = 0; j < N; ++j) {
				double xp[N], xm[N], fp[N], fm[N];
				std::memcpy(xp, hg.x, sizeof xp); std::memcpy(xm, hg.x, sizeof xm);
				xp[j] += step; xm[j] -= step;
				SpringForce(xp, &hg.theta0, fp);
				SpringForce(xm, &hg.theta0, fm);
				const double h2 = xp[j] - xm[j];
				for (int a = 0; a < N; ++a) KFD[a * N + j] = (fp[a] - fm[a]) / h2;
			}
			eK.Add(Rel(MaxAbsDiff(dense.K, KFD, N * N), MaxAbs(dense.K, N * N)));
		}
		std::printf("R2 기존 K 대 탄성력 유한차분 (%zu 힌지)\n", hs.size());
		rep.Gate("K", eK.Max() <= 1e-6, "rel max %.2e", eK.Max(), 1e-6);
		std::fprintf(js, "  \"R2_K_vs_fd\": {\"rel_max\": %.3e, \"rel_median\": %.3e},\n", eK.Max(), eK.Median());
	}

	// ---- R3′ (2 차, 판정) 시간 / R3 (1 차, 보고) ---------------------------------
	{
		std::mt19937_64 rng(7);
		const Grid grid = MakeGrid(200, rng);
		std::uniform_real_distribution<double> noise(-1., 1.);
		std::vector<double> T0(grid.nHinge);
		for (int s = 0; s < grid.nHinge; ++s) T0[s] = 0.01 * noise(rng);
		SXMesh mesh;
		FillMesh<Spec>(mesh, grid, T0);
		const sx::StencilView<Spec> view(mesh);

		enum { kE, kEOld, kENew, kImp, kCompactNew, kCompactImp, kFlatOld, kV };
		static const char* names[kV] = { "staged_eval", "staged_old_assembly", "staged_new_assembly", "eval_implicit",
			"compact_staged_new_assembly", "compact_eval_implicit", "flat_old_assembly" };
		double best[kV];
		for (int i = 0; i < kV; ++i) best[i] = 1e30;
		std::vector<double> Fn(3 * grid.nNode), Bn(3 * grid.nNode);
		double checksum = 0.;
		for (int r = 0; r < 5; ++r) {   // 변형을 번갈아 5 회, 변형마다 최소
			best[kE]          = std::min(best[kE],          TimeLoop<&Spec::EvalStaged, 0>(view, Fn, Bn, checksum));
			best[kEOld]       = std::min(best[kEOld],       TimeLoop<&Spec::EvalStaged, 1>(view, Fn, Bn, checksum));
			best[kENew]       = std::min(best[kENew],       TimeLoop<&Spec::EvalStaged, 2>(view, Fn, Bn, checksum));
			best[kImp]        = std::min(best[kImp],        TimeLoopImplicit<&Spec::EvalImplicit>(view, Fn, Bn, checksum));
			best[kCompactNew] = std::min(best[kCompactNew], TimeLoop<&Compact::EvalStaged, 2>(view, Fn, Bn, checksum));
			best[kCompactImp] = std::min(best[kCompactImp], TimeLoopImplicit<&Compact::EvalImplicit>(view, Fn, Bn, checksum));
			best[kFlatOld]    = std::min(best[kFlatOld],    TimeLoop<&Spec::Eval, 1>(view, Fn, Bn, checksum));
		}
		const double incOld = best[kEOld] - best[kE];
		const double incNew = best[kENew] - best[kE];
		const double incImp = best[kImp] - best[kE];
		// R3′ 판정은 tests/timing_bases.py 가 고정 로드 주소 8 곳의 최솟값으로 한다. 여기 값은 이 실행(로드 주소 하나)뿐이다.
		std::printf("R3′ 시간 — 보고만, 이 실행 한 번 (힌지 %d, 모으기 포함, 단일 스레드, 5 회 중 최소)\n", grid.nHinge);
		std::printf("  보고 — EvalImplicit 조립 몫 / 기존 조립 = %.3f (판정은 timing_bases.py)\n", incImp / incOld);
		std::printf("  보고 — 힌지당 ns: 평가 %.0f, + 기존 조립 %.0f (조립 %.0f), EvalImplicit %.0f (조립 몫 %.0f)\n",
			best[kE], best[kEOld], incOld, best[kImp], incImp);
		std::printf("  보고 — 1 차 AssembleImplicit: %.0f (조립 %.0f), 비 %.3f\n", best[kENew], incNew, incNew / incOld);
		std::printf("  보고 — 압축 명세: EvalStaged + 1 차 조립 %.0f, EvalImplicit %.0f ns\n", best[kCompactNew], best[kCompactImp]);
		std::printf("  보고 — 처음부터 끝까지: 착수 게이트 꼴(평탄 정확 + 기존 조립) %.0f → 압축 명세 EvalImplicit %.0f ns (%.1f 배)  (checksum %.3e)\n",
			best[kFlatOld], best[kCompactImp], best[kFlatOld] / best[kCompactImp], checksum);
		std::fprintf(js, "  \"R3p_time\": {\"hinges\": %d, \"assembly_old_ns\": %.1f, \"implicit_increment_ns\": %.1f, \"implicit_over_old\": %.3f, \"end_to_end_speedup\": %.2f},\n",
			grid.nHinge, incOld, incImp, incImp / incOld, best[kFlatOld] / best[kCompactImp]);
		std::fprintf(js, "  \"R3_first_attempt\": {\"assembly_new_ns\": %.1f, \"new_over_old\": %.3f},\n", incNew, incNew / incOld);
		std::fprintf(js, "  \"ns_per_hinge\": {");
		for (int i = 0; i < kV; ++i) std::fprintf(js, "%s\"%s\": %.1f", i ? ", " : "", names[i], best[i]);
		std::fprintf(js, "},\n");
	}

	std::fprintf(js, "  \"gates_failed\": %d\n}\n", rep.fail);
	std::printf("\n조립 게이트 2 차 R1′ · R2: %s (실패 %d) — R3′ 는 timing_bases.py\n", rep.fail == 0 ? "GO" : "NO-GO", rep.fail);
	return rep.fail == 0 ? 0 : 1;
}
