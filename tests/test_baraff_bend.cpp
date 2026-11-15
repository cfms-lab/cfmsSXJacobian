//---------------------------------------------------------------------------
// test_baraff_bend.cpp — 착수 게이트 G1 – G4 (docs/MODEL_PROBLEM.md) 와 보고용 측정.
// 실행 디렉터리 = 리포 루트(results/hinges.txt, results/mma_ref.txt 를 읽고 results/gate_baraff_bend.json 을 쓴다).
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffBend.h"
#include "gate_common.h"

#include <chrono>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;

typedef sxgen::BaraffBend Spec;
static const int N = Spec::kDof;      // 12
static const int NH = Spec::kHess;    // 78

// 감쇠력 d(x) = −kd ∇C (∇C·v)
static void DampingForce(const double* x, const double* ev, const double* v, double kd, double* d)
{
	double C, g[12], H[NH];
	Spec::Eval(x, ev, &C, g, H);
	double Cdot = 0.;
	for (int a = 0; a < N; ++a) Cdot += g[a] * v[a];
	for (int a = 0; a < N; ++a) d[a] = -kd * Cdot * g[a];
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Hinge> hs = ReadHinges("results/hinges.txt");
	if (hs.empty()) { std::printf("results/hinges.txt 없음 — model/make_hinges.py 를 먼저 돌린다\n"); return 2; }
	std::printf("표본 %zu 힌지\n", hs.size());
	Report rpt("results/gate_baraff_bend.json");
	std::FILE* g_json = rpt.json;
	std::fprintf(g_json, "{\n  \"samples\": %zu,\n", hs.size());

	// ---- G1 Mathematica 대조 ------------------------------------------------
	{
		const std::vector<std::vector<double> > ref = ReadRows("results/mma_ref.txt");
		double eC = 0., eG = 0., eH = 0.;
		int n = 0;
		for (std::size_t s = 0; s < ref.size() && s < hs.size(); ++s) {
			if (ref[s].size() != 1 + N + NH) continue;
			double C, g[12], H[NH];
			Spec::Eval(hs[s].x, &hs[s].theta0, &C, g, H);
			eC = std::max(eC, std::fabs(C - ref[s][0]));
			eG = std::max(eG, Rel(MaxAbsDiff(g, &ref[s][1], N), MaxAbs(&ref[s][1], N)));
			eH = std::max(eH, Rel(MaxAbsDiff(H, &ref[s][1 + N], NH), MaxAbs(&ref[s][1 + N], NH)));
			++n;
		}
		std::printf("G1 Mathematica 대조 (%d 힌지)\n", n);
		const bool ran = n > 0;
		rpt.Gate("C", ran && eC <= 1e-12, "|ΔC| %.2e", eC, 1e-12);
		rpt.Gate("∇C", ran && eG <= 1e-10, "rel %.2e", eG, 1e-10);
		rpt.Gate("∇²C", ran && eH <= 1e-10, "rel %.2e", eH, 1e-10);
		std::fprintf(g_json, "  \"G1_mathematica\": {\"hinges\": %d, \"C_abs\": %.3e, \"grad_rel\": %.3e, \"hess_rel\": %.3e},\n", n, eC, eG, eH);
	}

	// ---- G2 유한차분 · 보고: 길이 고정 근사 -----------------------------------
	{
		Stat fdG, fdH, frG, frH, frFdG, frC;
		for (std::size_t s = 0; s < hs.size(); ++s) {
			const Hinge& h = hs[s];
			const double step = 1e-6 * EdgeLength(h.x);
			double C, g[12], H[NH], gFD[12], HFD[NH];
			Spec::Eval(h.x, &h.theta0, &C, g, H);
			FiniteDiff<N>(&Spec::Eval, h.x, &h.theta0, step, gFD, HFD);
			fdG.Add(Rel(MaxAbsDiff(g, gFD, N), MaxAbs(g, N)));
			fdH.Add(Rel(MaxAbsDiff(H, HFD, NH), MaxAbs(H, NH)));

			double Cf, gf[12], Hf[NH];
			Spec::EvalFrozen(h.x, &h.theta0, &Cf, gf, Hf);
			frC.Add(std::fabs(Cf - C));
			frG.Add(Rel(MaxAbsDiff(gf, g, N), MaxAbs(g, N)));
			frH.Add(Rel(MaxAbsDiff(Hf, H, NH), MaxAbs(H, NH)));
			frFdG.Add(Rel(MaxAbsDiff(gf, gFD, N), MaxAbs(gFD, N)));
		}
		std::printf("G2 유한차분 (%zu 힌지)\n", hs.size());
		rpt.Gate("∇C", fdG.Max() <= 1e-6, "rel max %.2e", fdG.Max(), 1e-6);
		rpt.Gate("∇²C", fdH.Max() <= 1e-6, "rel max %.2e", fdH.Max(), 1e-6);
		std::printf("  보고 — 길이 고정 근사(원문 §4.3) 대 정확: C 차 max %.1e, ∇C rel 중앙 %.2e · max %.2e, ∇²C rel 중앙 %.2e · max %.2e\n",
			frC.Max(), frG.Median(), frG.Max(), frH.Median(), frH.Max());
		std::fprintf(g_json, "  \"G2_finite_difference\": {\"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e, \"grad_rel_median\": %.3e, \"hess_rel_median\": %.3e},\n",
			fdG.Max(), fdH.Max(), fdG.Median(), fdH.Median());
		std::fprintf(g_json, "  \"frozen_vs_exact\": {\"C_abs_max\": %.3e, \"grad_rel_median\": %.3e, \"grad_rel_max\": %.3e, \"hess_rel_median\": %.3e, \"hess_rel_max\": %.3e, \"grad_vs_fd_rel_max\": %.3e},\n",
			frC.Max(), frG.Median(), frG.Max(), frH.Median(), frH.Max(), frFdG.Max());
	}

	// ---- G3 감쇠 ∂d/∂x 유한차분 · 보고: 버린 비대칭 항 ------------------------
	{
		Stat fdDx, fdDv, dropped;
		const double kd = 1.;
		for (std::size_t s = 0; s < hs.size(); ++s) {
			const Hinge& h = hs[s];
			const double step = 1e-6 * EdgeLength(h.x);
			double C, g[12], H[NH];
			Spec::Eval(h.x, &h.theta0, &C, g, H);
			sx::ConditionForce<N> full, baraff;
			sx::AssembleConditionForce<N>(1, &C, g, H, h.v, 1., kd, true, full);
			sx::AssembleConditionForce<N>(1, &C, g, H, h.v, 1., kd, false, baraff);
			double DxFD[N * N];
			for (int b = 0; b < N; ++b) {
				double xp[12], xm[12], dp[12], dm[12];
				std::memcpy(xp, h.x, sizeof xp); std::memcpy(xm, h.x, sizeof xm);
				xp[b] += step; xm[b] -= step;
				DampingForce(xp, &h.theta0, h.v, kd, dp);
				DampingForce(xm, &h.theta0, h.v, kd, dm);
				const double h2 = xp[b] - xm[b];
				for (int a = 0; a < N; ++a) DxFD[a * N + b] = (dp[a] - dm[a]) / h2;
			}
			fdDx.Add(Rel(MaxAbsDiff(full.Dx, DxFD, N * N), MaxAbs(full.Dx, N * N)));
			// ∂d/∂v 는 v 에 선형: −kd g gᵀ 와 정확히 같아야 한다
			double DvRef[N * N];
			for (int a = 0; a < N; ++a) for (int b = 0; b < N; ++b) DvRef[a * N + b] = -kd * g[a] * g[b];
			fdDv.Add(Rel(MaxAbsDiff(full.Dv, DvRef, N * N), MaxAbs(DvRef, N * N)));
			dropped.Add(Rel(MaxAbsDiff(full.Dx, baraff.Dx, N * N), MaxAbs(full.Dx, N * N)));
		}
		std::printf("G3 감쇠 Jacobian (%zu 힌지)\n", hs.size());
		rpt.Gate("∂d/∂x", fdDx.Max() <= 1e-6, "rel max %.2e", fdDx.Max(), 1e-6);
		rpt.Gate("∂d/∂v", fdDv.Max() <= 1e-15, "rel max %.2e", fdDv.Max(), 1e-15);
		std::printf("  보고 — 원 논문이 버린 비대칭 항 / ∂d/∂x: 중앙 %.2e, max %.2e\n", dropped.Median(), dropped.Max());
		std::fprintf(g_json, "  \"G3_damping\": {\"dx_rel_max\": %.3e, \"dv_rel_max\": %.3e, \"dropped_asym_rel_median\": %.3e, \"dropped_asym_rel_max\": %.3e},\n",
			fdDx.Max(), fdDv.Max(), dropped.Median(), dropped.Max());
	}

	// ---- G4 SXMesh 스텐실 루프 · 보고: 시간 ---------------------------------
	{
		std::mt19937_64 rng(7);
		const Grid grid = MakeGrid(200, rng);                  // 200 × 200 노드, 간격 1 cm
		const int nNode = grid.nNode, nHinge = grid.nHinge;
		const std::vector<double>& X = grid.X;
		const std::vector<double>& V = grid.V;
		const std::vector<int>& hinges = grid.hinges;
		std::uniform_real_distribution<double> noise(-1., 1.);
		std::vector<double> T0(nHinge);
		for (int s = 0; s < nHinge; ++s) T0[s] = 0.01 * noise(rng);
		SXMesh mesh;
		FillMesh<Spec>(mesh, grid, T0);   // 명세의 선언대로 SXMesh 변수를 만들고 채운다

		sx::StencilView<Spec> view(mesh);
		// 비트 비교: SXMesh 를 거친 값 대 직접 모은 값
		int mismatch = 0;
		for (int s = 0; s < nHinge; ++s) {
			double xa[12], va[12], eva[1], xb[12], vb[12];
			if (!view.Gather(s, xa, va, eva)) { ++mismatch; continue; }
			for (int a = 0; a < 4; ++a) for (int c = 0; c < 3; ++c) { xb[3 * a + c] = X[3 * hinges[4 * s + a] + c]; vb[3 * a + c] = V[3 * hinges[4 * s + a] + c]; }
			double Ca, ga[12], Ha[NH], Cb, gb[12], Hb[NH];
			Spec::Eval(xa, eva, &Ca, ga, Ha);
			Spec::Eval(xb, &T0[s], &Cb, gb, Hb);
			if (std::memcmp(&Ca, &Cb, sizeof Ca) != 0 || std::memcmp(ga, gb, sizeof ga) != 0 || std::memcmp(Ha, Hb, sizeof Ha) != 0 ||
				std::memcmp(va, vb, sizeof va) != 0) ++mismatch;
		}
		std::printf("G4 SXMesh 스텐실 루프 (노드 %d, 힌지 %d)\n", nNode, nHinge);
		rpt.Gate("루프", mismatch == 0, "불일치 %.0f", static_cast<double>(mismatch), 0.);

		// 시간(단일 스레드, 3회 중 최소)
		typedef std::chrono::steady_clock Clock;
		double best[4] = { 1e30, 1e30, 1e30, 1e30 };
		double checksum = 0.;
		std::vector<double> F(3 * nNode);
		for (int rep = 0; rep < 3; ++rep) {
			for (int mode = 0; mode < 4; ++mode) {
				std::fill(F.begin(), F.end(), 0.);
				const Clock::time_point t0 = Clock::now();
				for (int s = 0; s < nHinge; ++s) {
					double x[12], v[12], ev[1], C, g[12], H[NH];
					view.Gather(s, x, v, ev);
					if (mode == 0)      { Spec::Eval(x, ev, &C, g, H); checksum += H[0]; }
					else if (mode == 1) { Spec::EvalFrozen(x, ev, &C, g, H); checksum += H[0]; }
					else {
						if (mode == 2) Spec::Eval(x, ev, &C, g, H); else Spec::EvalFrozen(x, ev, &C, g, H);
						sx::ConditionForce<N> out;
						sx::AssembleConditionForce<N>(1, &C, g, H, v, 1., 0.1, false, out);
						const int* nd = view.Nodes(s);
						for (int a = 0; a < 4; ++a) for (int c = 0; c < 3; ++c) F[3 * nd[a] + c] += out.f[3 * a + c] + out.d[3 * a + c];
						checksum += out.K[13] + out.Dx[7];
					}
				}
				const double ns = std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / nHinge;
				best[mode] = std::min(best[mode], ns);
			}
		}
		std::printf("  보고 — 힌지당 ns (단일 스레드, 모으기 포함): Eval %.0f, EvalFrozen %.0f, Eval+조립 %.0f, EvalFrozen+조립 %.0f  (checksum %.3e)\n",
			best[0], best[1], best[2], best[3], checksum);
		std::fprintf(g_json, "  \"G4_stencil_loop\": {\"nodes\": %d, \"hinges\": %d, \"mismatch\": %d},\n", nNode, nHinge, mismatch);
		std::fprintf(g_json, "  \"ns_per_hinge\": {\"eval\": %.1f, \"eval_frozen\": %.1f, \"eval_assemble\": %.1f, \"eval_frozen_assemble\": %.1f},\n",
			best[0], best[1], best[2], best[3]);
	}

	std::fprintf(g_json, "  \"gates_failed\": %d\n}\n", rpt.fail);
	std::printf("\n게이트 G1 – G4: %s (실패 %d)\n", rpt.fail == 0 ? "GO" : "NO-GO", rpt.fail);
	return rpt.fail == 0 ? 0 : 1;
}
