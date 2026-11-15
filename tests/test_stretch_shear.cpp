//---------------------------------------------------------------------------
// test_stretch_shear.cpp — 인장 · 전단 게이트 T1 – T5 (docs/GATE_STRETCH_SHEAR.md) 와 보고용 측정.
// 실행 디렉터리 = 리포 루트(results/triangles.txt, results/mma_ref_tri.txt 를 읽고 results/gate_stretch_shear.json 을 쓴다).
// 시간은 이 exe 한 번의 값이라 보고만 한다(로드 주소 함정, docs/GATE_SPEED.md 「측정 방법 변경」).
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffStretch.h"
#include "../generated/BaraffShear.h"
#include "gate_common.h"

#include <chrono>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;

static const double kK = 1.0, kKd = 0.1, kH = 0.01;   // 판정문에 고정한 계수(조립 게이트와 같음)
static const double kPrm[3] = { kK, kKd, kH };
static const int kMmaCols = 165;                      // 인장 2 + 18 + 90, 전단 1 + 9 + 45

// C 오차의 척도 √a (C = √a × 차원 없는 양, docs/GATE_STRETCH_SHEAR.md 「측정 방법 변경」)
static double CScale(const double* uv) { return std::sqrt(UvArea(uv)); }

static double Len3(const double* a, const double* b)
{
	const double dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2];
	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// C · ∇C · ∇²C 오차: C 는 √a 로 나눈 절대, ∇C · ∇²C 는 삼각형마다 모든 성분에 걸친 max|Δ|/max|기준|
struct Err {
	Stat c, g, h;
	template <class Spec>
	void Add(const double* C, const double* g_, const double* H, const double* Cr, const double* gr, const double* Hr, double scale)
	{
		const int NC = Spec::kComp, N = Spec::kDof, NH = Spec::kHess;
		c.Add(MaxAbsDiff(C, Cr, NC) / scale);
		g.Add(Rel(MaxAbsDiff(g_, gr, NC * N), MaxAbs(gr, NC * N)));
		h.Add(Rel(MaxAbsDiff(H, Hr, NC * NH), MaxAbs(Hr, NC * NH)));
	}
};

// 성분 여럿 중심차분: gFD[c·N + b] = ∂C_c/∂x_b, HFD[c·NH + (a ≤ b)] = ∂g_{c,a}/∂x_b. 분모 = 실제 간격
template <class Spec>
static void FiniteDiffN(EvalFn f, const double* x, const double* ev, double step, double* gFD, double* HFD)
{
	enum { N = Spec::kDof, NC = Spec::kComp, NH = Spec::kHess };
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

// 명세 하나의 T1 – T5 와 보고. off = Mathematica 한 줄에서 이 명세 값이 시작하는 자리
template <class Spec>
static void RunSpec(const char* name, const std::vector<Tri>& ts, const std::vector<std::vector<double> >& ref, int off,
                    Report& rpt)
{
	enum { N = Spec::kDof, NC = Spec::kComp, NH = Spec::kHess };
	static_assert(Spec::kNodes == 3 && Spec::kElemValues == 6, "삼각형 · uv 6 칸");
	std::FILE* js = rpt.json;
	std::printf("\n== %s (성분 %d, 자유도 %d) ==\n", name, NC, N);
	std::fprintf(js, "  \"%s\": {\n", name);

	// ---- T1 Mathematica ------------------------------------------------------
	{
		Err eFlat, eStaged;
		int n = 0;
		for (std::size_t s = 0; s < ref.size() && s < ts.size(); ++s) {
			if (static_cast<int>(ref[s].size()) != kMmaCols) continue;
			const double* Cr = &ref[s][off];
			const double* gr = Cr + NC;
			const double* Hr = gr + NC * N;
			const double scale = CScale(ts[s].uv);
			double C[NC], g[NC * N], H[NC * NH];
			Spec::Eval(ts[s].x, ts[s].uv, C, g, H);
			eFlat.Add<Spec>(C, g, H, Cr, gr, Hr, scale);
			Spec::EvalStaged(ts[s].x, ts[s].uv, C, g, H);
			eStaged.Add<Spec>(C, g, H, Cr, gr, Hr, scale);
			++n;
		}
		const bool ran = n > 0;
		std::printf("T1 Mathematica 대조 (%d 삼각형)\n", n);
		const Err* es[2] = { &eFlat, &eStaged };
		const char* tag[2] = { "Eval", "EvalStaged" };
		for (int i = 0; i < 2; ++i) {
			std::printf(" %s\n", tag[i]);
			rpt.Gate("C/√a", ran && es[i]->c.Max() <= 1e-12, "|ΔC|/√a %.2e", es[i]->c.Max(), 1e-12);
			rpt.Gate("∇C", ran && es[i]->g.Max() <= 1e-10, "rel %.2e", es[i]->g.Max(), 1e-10);
			rpt.Gate("∇²C", ran && es[i]->h.Max() <= 1e-10, "rel %.2e", es[i]->h.Max(), 1e-10);
		}
		std::fprintf(js, "    \"T1_mathematica\": {\"triangles\": %d, ", n);
		JsonErr(js, "eval", eFlat, ", ");
		JsonErr(js, "staged", eStaged, "},\n");
	}

	// ---- T2 유한차분 · T3 단계판 대 평탄판 · 보고 ---------------------------
	Err fdFlat, fdStaged, st;
	Stat sumG, torqueG;
	int frozenDiff = 0;
	for (std::size_t s = 0; s < ts.size(); ++s) {
		const Tri& t = ts[s];
		const double scale = CScale(t.uv);
		const double step = 1e-6 * Len3(t.x, t.x + 3);
		double C[NC], g[NC * N], H[NC * NH], Cs[NC], gs[NC * N], Hs[NC * NH], gFD[NC * N], HFD[NC * NH];
		Spec::Eval(t.x, t.uv, C, g, H);
		Spec::EvalStaged(t.x, t.uv, Cs, gs, Hs);
		FiniteDiffN<Spec>(&Spec::Eval, t.x, t.uv, step, gFD, HFD);
		fdFlat.Add<Spec>(C, g, H, C, gFD, HFD, scale);
		FiniteDiffN<Spec>(&Spec::EvalStaged, t.x, t.uv, step, gFD, HFD);
		fdStaged.Add<Spec>(Cs, gs, Hs, Cs, gFD, HFD, scale);
		st.Add<Spec>(Cs, gs, Hs, C, g, H, scale);

		// 보고: 성분마다 Σᵢ ∇ᵢC = 0(병진), Σᵢ (xᵢ − x̄) × ∇ᵢC = 0(회전)
		double xc[3];
		for (int k = 0; k < 3; ++k) xc[k] = (t.x[k] + t.x[3 + k] + t.x[6 + k]) / 3.;
		double rmax = 0.;
		for (int a = 0; a < 3; ++a) for (int k = 0; k < 3; ++k) rmax = std::max(rmax, std::fabs(t.x[3 * a + k] - xc[k]));
		for (int c = 0; c < NC; ++c) {
			const double* gc = g + c * N;
			const double gm = MaxAbs(gc, N);
			double sg[3] = { 0., 0., 0. }, tq[3] = { 0., 0., 0. };
			for (int a = 0; a < 3; ++a) {
				const double r[3] = { t.x[3 * a] - xc[0], t.x[3 * a + 1] - xc[1], t.x[3 * a + 2] - xc[2] };
				const double* q = gc + 3 * a;
				for (int k = 0; k < 3; ++k) sg[k] += q[k];
				tq[0] += r[1] * q[2] - r[2] * q[1];
				tq[1] += r[2] * q[0] - r[0] * q[2];
				tq[2] += r[0] * q[1] - r[1] * q[0];
			}
			sumG.Add(Rel(MaxAbs(sg, 3), gm));
			torqueG.Add(Rel(MaxAbs(tq, 3), gm * rmax));
		}
		// 보고: EvalFrozen = Eval (ctx.unit 을 쓰지 않으므로 같아야 한다)
		double Cf[NC], gf[NC * N], Hf[NC * NH];
		Spec::EvalFrozen(t.x, t.uv, Cf, gf, Hf);
		if (std::memcmp(Cf, C, sizeof Cf) != 0 || std::memcmp(gf, g, sizeof gf) != 0 || std::memcmp(Hf, H, sizeof Hf) != 0) ++frozenDiff;
	}
	std::printf("T2 유한차분 (%zu 삼각형)\n", ts.size());
	rpt.Gate("∇C", fdFlat.g.Max() <= 1e-6, "Eval rel max %.2e", fdFlat.g.Max(), 1e-6);
	rpt.Gate("∇²C", fdFlat.h.Max() <= 1e-6, "Eval rel max %.2e", fdFlat.h.Max(), 1e-6);
	rpt.Gate("∇C", fdStaged.g.Max() <= 1e-6, "EvalStaged rel max %.2e", fdStaged.g.Max(), 1e-6);
	rpt.Gate("∇²C", fdStaged.h.Max() <= 1e-6, "EvalStaged rel max %.2e", fdStaged.h.Max(), 1e-6);
	std::printf("T3 단계판 대 평탄판 (%zu 삼각형)\n", ts.size());
	rpt.Gate("C/√a", st.c.Max() <= 1e-12, "|ΔC|/√a %.2e", st.c.Max(), 1e-12);
	rpt.Gate("∇C", st.g.Max() <= 1e-12, "rel %.2e", st.g.Max(), 1e-12);
	rpt.Gate("∇²C", st.h.Max() <= 1e-12, "rel %.2e", st.h.Max(), 1e-12);
	std::fprintf(js, "    \"T2_finite_difference\": {\"eval\": {\"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e}, \"staged\": {\"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e}},\n",
		fdFlat.g.Max(), fdFlat.h.Max(), fdStaged.g.Max(), fdStaged.h.Max());
	std::fprintf(js, "    ");
	JsonErr(js, "T3_staged_vs_flat", st, ",\n");

	// ---- T4 EvalImplicit 대 AssembleImplicit ---------------------------------
	{
		Stat eA, eB, eF;
		for (std::size_t s = 0; s < ts.size(); ++s) {
			const Tri& t = ts[s];
			double C[NC], g[NC * N], H[NC * NH];
			Spec::Eval(t.x, t.uv, C, g, H);
			sx::ImplicitBlock<N> blk;
			sx::AssembleImplicit<N>(NC, C, g, H, t.v, kK, kKd, kH, blk);
			double Ci[NC], f[N], A[NH], b[N];
			Spec::EvalImplicit(t.x, t.uv, t.v, kPrm, Ci, f, A, b);
			eA.Add(Rel(MaxAbsDiff(A, blk.A, NH), MaxAbs(blk.A, NH)));
			eB.Add(Rel(MaxAbsDiff(b, blk.b, N), MaxAbs(blk.b, N)));
			eF.Add(Rel(MaxAbsDiff(f, blk.f, N), MaxAbs(blk.f, N)));
		}
		std::printf("T4 EvalImplicit 대 AssembleImplicit (%zu 삼각형, k = %g, kd = %g, h = %g)\n", ts.size(), kK, kKd, kH);
		rpt.Gate("A", eA.Max() <= 1e-12, "rel max %.2e", eA.Max(), 1e-12);
		rpt.Gate("b", eB.Max() <= 1e-12, "rel max %.2e", eB.Max(), 1e-12);
		rpt.Gate("f+d", eF.Max() <= 1e-12, "rel max %.2e", eF.Max(), 1e-12);
		std::fprintf(js, "    \"T4_implicit\": {\"A_rel\": %.3e, \"b_rel\": %.3e, \"f_rel\": %.3e, \"A_rel_median\": %.3e, \"b_rel_median\": %.3e, \"f_rel_median\": %.3e},\n",
			eA.Max(), eB.Max(), eF.Max(), eA.Median(), eB.Median(), eF.Median());
	}

	// 보고: 쉬는 모양(표본 0, F = I 를 강체 이동)의 C/√a
	double Crest[NC], grest[NC * N], Hrest[NC * NH];
	Spec::Eval(ts[0].x, ts[0].uv, Crest, grest, Hrest);
	const double restC = MaxAbs(Crest, NC) / CScale(ts[0].uv);
	std::printf("  보고 — 병진 |Σ∇C|/max|∇C| max %.1e, 회전 |Σ r×∇C| max %.1e, 쉬는 모양 |C|/√a %.1e, EvalFrozen ≠ Eval %d 개\n",
		sumG.Max(), torqueG.Max(), restC, frozenDiff);
	std::fprintf(js, "    \"report\": {\"sum_grad_rel_max\": %.3e, \"torque_grad_rel_max\": %.3e, \"rest_C_over_sqrt_a\": %.3e, \"frozen_differs\": %d},\n",
		sumG.Max(), torqueG.Max(), restC, frozenDiff);

	// ---- T5 SXMesh 삼각형 루프 · 보고: 시간 ---------------------------------
	{
		std::mt19937_64 rng(7);
		const int G = 200;
		const Grid grid = MakeGrid(G, rng);
		const int nNode = grid.nNode, nTri = grid.nTri;
		std::vector<double> UV(6 * static_cast<std::size_t>(nTri));   // 쉬는 uv = 평평한 격자(간격 1 cm)
		for (int s = 0; s < nTri; ++s)
			for (int a = 0; a < 3; ++a) {
				const int nd = grid.tris[3 * s + a];
				UV[6 * s + 2 * a] = 0.01 * (nd % G);
				UV[6 * s + 2 * a + 1] = 0.01 * (nd / G);
			}
		SXMesh mesh;
		sx::CreateVariables<Spec>(mesh, nNode, nTri);
		SXNodeNamedContainer& px = mesh.GetNodeV(Spec::PositionVar());
		SXNodeNamedContainer& pv = mesh.GetNodeV(Spec::VelocityVar());
		SXElementNamedContainer& pt = mesh.GetElementV(Spec::StencilVar());
		SXNodeNamedContainer& puv = mesh.GetNodeV(Spec::ElemVar(0));
		for (int n = 0; n < nNode; ++n) { px.Add(n, &grid.X[3 * n], 3); pv.Add(n, &grid.V[3 * n], 3); }
		for (int s = 0; s < nTri; ++s) { pt.Add(s, &grid.tris[3 * s], 3); for (int c = 0; c < 6; ++c) puv[s][c] = UV[6 * s + c]; }

		sx::StencilView<Spec> view(mesh);
		int mismatch = 0;
		for (int s = 0; s < nTri; ++s) {
			double xa[N], va[N], eva[6], xb[N], vb[N];
			if (!view.Gather(s, xa, va, eva)) { ++mismatch; continue; }
			for (int a = 0; a < 3; ++a) for (int c = 0; c < 3; ++c) { xb[3 * a + c] = grid.X[3 * grid.tris[3 * s + a] + c]; vb[3 * a + c] = grid.V[3 * grid.tris[3 * s + a] + c]; }
			double Ca[NC], ga[NC * N], Ha[NC * NH], Cb[NC], gb[NC * N], Hb[NC * NH];
			Spec::Eval(xa, eva, Ca, ga, Ha);
			Spec::Eval(xb, &UV[6 * s], Cb, gb, Hb);
			if (std::memcmp(Ca, Cb, sizeof Ca) != 0 || std::memcmp(ga, gb, sizeof ga) != 0 || std::memcmp(Ha, Hb, sizeof Ha) != 0 ||
				std::memcmp(va, vb, sizeof va) != 0) ++mismatch;
		}
		std::printf("T5 SXMesh 삼각형 루프 (노드 %d, 삼각형 %d)\n", nNode, nTri);
		rpt.Gate("루프", mismatch == 0, "불일치 %.0f", static_cast<double>(mismatch), 0.);

		typedef std::chrono::steady_clock Clock;
		double best[3] = { 1e30, 1e30, 1e30 };
		double checksum = 0.;
		std::vector<double> F(3 * static_cast<std::size_t>(nNode));
		for (int rep = 0; rep < 5; ++rep) {
			for (int mode = 0; mode < 3; ++mode) {
				std::fill(F.begin(), F.end(), 0.);
				const Clock::time_point t0 = Clock::now();
				for (int s = 0; s < nTri; ++s) {
					double x[N], v[N], ev[6], C[NC], g[NC * N], H[NC * NH];
					view.Gather(s, x, v, ev);
					if (mode == 0)      { Spec::Eval(x, ev, C, g, H); checksum += H[0]; }
					else if (mode == 1) { Spec::EvalStaged(x, ev, C, g, H); checksum += H[0]; }
					else {
						double f[N], A[NH], b[N];
						Spec::EvalImplicit(x, ev, v, kPrm, C, f, A, b);
						const int* nd = view.Nodes(s);
						for (int a = 0; a < 3; ++a) for (int c = 0; c < 3; ++c) F[3 * nd[a] + c] += f[3 * a + c];
						checksum += A[0] + b[0];
					}
				}
				const double ns = std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / nTri;
				best[mode] = std::min(best[mode], ns);
			}
		}
		std::printf("  보고 — 삼각형당 ns (단일 스레드, 모으기 포함, 이 exe 한 번): Eval %.0f, EvalStaged %.0f, EvalImplicit+힘 흩뿌리기 %.0f  (checksum %.3e)\n",
			best[0], best[1], best[2], checksum);
		std::fprintf(js, "    \"T5_stencil_loop\": {\"nodes\": %d, \"triangles\": %d, \"mismatch\": %d},\n", nNode, nTri, mismatch);
		std::fprintf(js, "    \"ns_per_triangle\": {\"eval\": %.1f, \"staged\": %.1f, \"implicit\": %.1f}\n  },\n",
			best[0], best[1], best[2]);
	}
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Tri> ts = ReadTriangles("results/triangles.txt");
	if (ts.empty()) { std::printf("results/triangles.txt 없음 — model/make_triangles.py 를 먼저 돌린다\n"); return 2; }
	const std::vector<std::vector<double> > ref = ReadRows("results/mma_ref_tri.txt");
	std::printf("인장 · 전단 게이트 — 표본 %zu 삼각형, Mathematica 기준 %zu 줄\n", ts.size(), ref.size());
	Report rpt("results/gate_stretch_shear.json");
	std::fprintf(rpt.json, "{\n  \"samples\": %zu, \"k\": %g, \"kd\": %g, \"h\": %g,\n", ts.size(), kK, kKd, kH);

	RunSpec<sxgen::BaraffStretch>("BaraffStretch", ts, ref, 0, rpt);
	RunSpec<sxgen::BaraffShear>("BaraffShear", ts, ref, 2 + 18 + 90, rpt);

	std::fprintf(rpt.json, "  \"gates_failed\": %d\n}\n", rpt.fail);
	std::printf("\n게이트 T1 – T5: %s (실패 %d)\n", rpt.fail == 0 ? "GO" : "NO-GO", rpt.fail);
	return rpt.fail == 0 ? 0 : 1;
}
