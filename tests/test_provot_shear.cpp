//---------------------------------------------------------------------------
// test_provot_shear.cpp — Provot 전단 게이트 V1 – V7 (docs/GATE_PROVOT_SHEAR.md) 와 보고용 측정.
// 인장 · 전단 게이트(test_stretch_shear.cpp)의 T1 – T5 와 인장 판 게이트(test_tension.cpp)의 P1 – P2 를 명세 ProvotShear 에 옮긴 것이다.
// 실행 디렉터리 = 리포 루트. results/triangles.txt · results/provot_diag.txt · results/mma_ref_provot.txt 를 읽고
// results/gate_provot_shear.json 을 쓴다. V0(기존 변형 불변)은 check_flat_unchanged.py · check_tension_p0.py · check_link_l1_p0.py.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/ProvotShear.h"
#include "gate_common.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;
typedef sxgen::ProvotShear S;

enum { N = S::kDof, NC = S::kComp, NH = S::kHess, NE = S::kElemValues };
static_assert(S::kNodes == 3 && N == 9 && NC == 2 && NE == 8 && S::kElemVars == 2, "삼각형 · 성분 2 · uv 6 칸 + 대각선 2 칸");

static const double kK = 1.0, kKd = 0.1, kH = 0.01;   // 판정문에 고정한 계수(조립 · 인장 판 게이트와 같음)
static const double kPrm[3] = { kK, kKd, kH };
static const int kMmaCols = 2 + 18 + 90 + 9;         // C · g · H · 엔진 꼴 힘
static const double kLimit = 1e-12;                  // V3 · V4 · V6 · V6b
static const double kPsd = -1e-12;                   // V7: 최소 고윳값 / 최대 고윳값 크기
static const int kMinBranch = 100;                   // V6 갈래

struct Sample { double x[9]; double ev[8]; double v[9]; };   // ev = uv(6) · b₊ b₋

static double CScale(const double* ev) { return std::sqrt(UvArea(ev)); }   // C 오차의 척도 √a

static double Len3(const double* a, const double* b)
{
	const double dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2];
	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// C · ∇C · ∇²C 오차: C 는 √a 로 나눈 절대, ∇C · ∇²C 는 삼각형마다 모든 성분에 걸친 max|Δ|/max|기준|
struct Err {
	Stat c, g, h;
	void Add(const double* C, const double* g_, const double* H, const double* Cr, const double* gr, const double* Hr, double scale)
	{
		c.Add(MaxAbsDiff(C, Cr, NC) / scale);
		g.Add(Rel(MaxAbsDiff(g_, gr, NC * N), MaxAbs(gr, NC * N)));
		h.Add(Rel(MaxAbsDiff(H, Hr, NC * NH), MaxAbs(Hr, NC * NH)));
	}
};

// 성분 여럿 중심차분: gFD[c·N + b] = ∂C_c/∂x_b, HFD[c·NH + (a ≤ b)] = ∂g_{c,a}/∂x_b. 분모 = 실제 간격
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

// 대칭 행렬(n × n 전체, 행 우선)의 고윳값 — 순환 Jacobi(인장 판 게이트와 같은 방법). M 은 망가진다
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
	double M[N * N], lam[N];
	for (int a = 0; a < N; ++a)
		for (int b = 0; b < N; ++b) M[a * N + b] = U[a <= b ? sx::UpperIndex<N>(a, b) : sx::UpperIndex<N>(b, a)];
	SymEig(M, N, lam);
	double lo = lam[0], big = 0.;
	for (int p = 0; p < N; ++p) { lo = std::min(lo, lam[p]); big = std::max(big, std::fabs(lam[p])); }
	return big > 0. ? lo / big : 0.;
}

// T 판의 기준 조립: AssembleImplicit(∇²C = 0) 에 h² k Σ C⁺ ∇²C 를 A 에 더하고 b 에서 h² k Σ C⁺ ∇²C v 를 뺀다
static void Reference(const double* C, const double* g, const double* H, const double* v, sx::ImplicitBlock<N>& ref)
{
	double H0[NC * NH];
	for (int i = 0; i < NC * NH; ++i) H0[i] = 0.;
	sx::AssembleImplicit<N>(NC, C, g, H0, v, kK, kKd, kH, ref);
	const double hhk = kH * kH * kK;
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

// 보고용: 두 대각선 w_u ± w_v 의 길이(명세와 같은 식을 손으로)
static void DiagLengths(const Sample& t, double* lenP, double* lenM)
{
	const double* uv = t.ev;
	const double du1 = uv[2] - uv[0], dv1 = uv[3] - uv[1], du2 = uv[4] - uv[0], dv2 = uv[5] - uv[1];
	const double det = du1 * dv2 - du2 * dv1;
	double dp[3], dm[3];
	for (int k = 0; k < 3; ++k) {
		const double dx1 = t.x[3 + k] - t.x[k], dx2 = t.x[6 + k] - t.x[k];
		const double wu = (dx1 * dv2 - dx2 * dv1) / det, wv = (-dx1 * du2 + dx2 * du1) / det;
		dp[k] = wu + wv; dm[k] = wu - wv;
	}
	*lenP = std::sqrt(dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]);
	*lenM = std::sqrt(dm[0] * dm[0] + dm[1] * dm[1] + dm[2] * dm[2]);
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Tri> tris = ReadTriangles("results/triangles.txt");
	const std::vector<std::vector<double> > diag = ReadRows("results/provot_diag.txt");
	if (tris.empty() || diag.size() < tris.size()) {
		std::printf("results/triangles.txt · results/provot_diag.txt 없음 — model/make_triangles.py · model/make_provot_diag.py 를 먼저 돌린다\n");
		return 2;
	}
	std::vector<Sample> ts(tris.size());
	for (std::size_t s = 0; s < tris.size(); ++s) {
		if (diag[s].size() != 2) { std::printf("results/provot_diag.txt %zu 번째 줄이 두 칸이 아니다\n", s + 1); return 2; }
		std::memcpy(ts[s].x, tris[s].x, sizeof ts[s].x);
		std::memcpy(ts[s].ev, tris[s].uv, 6 * sizeof(double));
		ts[s].ev[6] = diag[s][0]; ts[s].ev[7] = diag[s][1];
		std::memcpy(ts[s].v, tris[s].v, sizeof ts[s].v);
	}
	const std::vector<std::vector<double> > ref = ReadRows("results/mma_ref_provot.txt");
	std::printf("Provot 전단 게이트 — 표본 %zu 삼각형, Mathematica 기준 %zu 줄 (성분 %d, 자유도 %d, 요소 값 %d)\n",
		ts.size(), ref.size(), static_cast<int>(NC), static_cast<int>(N), static_cast<int>(NE));
	Report rpt("results/gate_provot_shear.json");
	std::FILE* js = rpt.json;
	std::fprintf(js, "{\n  \"samples\": %zu, \"k\": %g, \"kd\": %g, \"h\": %g,\n", ts.size(), kK, kKd, kH);

	// ---- V1 Mathematica · V1b 엔진 꼴 힘 ------------------------------------
	{
		Err eFlat, eStaged;
		Stat eForce;
		int n = 0;
		for (std::size_t s = 0; s < ref.size() && s < ts.size(); ++s) {
			if (static_cast<int>(ref[s].size()) != kMmaCols) continue;
			const double* Cr = &ref[s][0];
			const double* gr = Cr + NC;
			const double* Hr = gr + NC * N;
			const double* fe = Hr + NC * NH;
			const double scale = CScale(ts[s].ev);
			double C[NC], g[NC * N], H[NC * NH];
			S::Eval(ts[s].x, ts[s].ev, C, g, H);
			eFlat.Add(C, g, H, Cr, gr, Hr, scale);
			double f[N];
			for (int b = 0; b < N; ++b) { f[b] = 0.; for (int c = 0; c < NC; ++c) f[b] -= C[c] * g[c * N + b]; }
			eForce.Add(MaxAbsDiff(f, fe, N) / (scale * MaxAbs(g, NC * N)));
			S::EvalStaged(ts[s].x, ts[s].ev, C, g, H);
			eStaged.Add(C, g, H, Cr, gr, Hr, scale);
			++n;
		}
		const bool ran = n > 0;
		std::printf("V1 Mathematica 대조 (%d 삼각형)\n", n);
		const Err* es[2] = { &eFlat, &eStaged };
		const char* tag[2] = { "Eval", "EvalStaged" };
		for (int i = 0; i < 2; ++i) {
			std::printf(" %s\n", tag[i]);
			rpt.Gate("C/√a", ran && es[i]->c.Max() <= 1e-12, "|ΔC|/√a %.2e", es[i]->c.Max(), 1e-12);
			rpt.Gate("∇C", ran && es[i]->g.Max() <= 1e-10, "rel %.2e", es[i]->g.Max(), 1e-10);
			rpt.Gate("∇²C", ran && es[i]->h.Max() <= 1e-10, "rel %.2e", es[i]->h.Max(), 1e-10);
		}
		std::printf("V1b 엔진 꼴 힘 −(1/a) Σ C_e ∇C_e 대 생성 −Σ C ∇C (%d 삼각형)\n", n);
		rpt.Gate("힘", ran && eForce.Max() <= 1e-10, "|Δf| / (√a · max|∇C|) %.2e", eForce.Max(), 1e-10);
		std::fprintf(js, "  \"V1_mathematica\": {\"triangles\": %d, ", n);
		JsonErr(js, "eval", eFlat, ", ");
		JsonErr(js, "staged", eStaged, "},\n");
		std::fprintf(js, "  \"V1b_engine_form_force\": {\"triangles\": %d, \"abs_over_sqrt_a_grad\": %.3e, \"median\": %.3e},\n", n, eForce.Max(), eForce.Median());
	}

	// ---- V2 유한차분 · V3 단계판 대 평탄판 · 보고 ---------------------------
	Err fdFlat, fdStaged, st;
	Stat sumG, torqueG;
	int frozenDiff = 0;
	double dMin = 1e30, dMax = 0.;
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
		double Cf[NC], gf[NC * N], Hf[NC * NH];
		S::EvalFrozen(t.x, t.ev, Cf, gf, Hf);
		if (std::memcmp(Cf, C, sizeof Cf) != 0 || std::memcmp(gf, g, sizeof gf) != 0 || std::memcmp(Hf, H, sizeof Hf) != 0) ++frozenDiff;
		double lp, lm;
		DiagLengths(t, &lp, &lm);
		dMin = std::min(dMin, std::min(lp, lm)); dMax = std::max(dMax, std::max(lp, lm));
	}
	std::printf("V2 유한차분 (%zu 삼각형)\n", ts.size());
	rpt.Gate("∇C", fdFlat.g.Max() <= 1e-6, "Eval rel max %.2e", fdFlat.g.Max(), 1e-6);
	rpt.Gate("∇²C", fdFlat.h.Max() <= 1e-6, "Eval rel max %.2e", fdFlat.h.Max(), 1e-6);
	rpt.Gate("∇C", fdStaged.g.Max() <= 1e-6, "EvalStaged rel max %.2e", fdStaged.g.Max(), 1e-6);
	rpt.Gate("∇²C", fdStaged.h.Max() <= 1e-6, "EvalStaged rel max %.2e", fdStaged.h.Max(), 1e-6);
	std::printf("V3 단계판 대 평탄판 (%zu 삼각형)\n", ts.size());
	rpt.Gate("C/√a", st.c.Max() <= kLimit, "|ΔC|/√a %.2e", st.c.Max(), kLimit);
	rpt.Gate("∇C", st.g.Max() <= kLimit, "rel %.2e", st.g.Max(), kLimit);
	rpt.Gate("∇²C", st.h.Max() <= kLimit, "rel %.2e", st.h.Max(), kLimit);
	std::fprintf(js, "  \"V2_finite_difference\": {\"eval\": {\"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e}, \"staged\": {\"grad_rel_max\": %.3e, \"hess_rel_max\": %.3e}},\n",
		fdFlat.g.Max(), fdFlat.h.Max(), fdStaged.g.Max(), fdStaged.h.Max());
	std::fprintf(js, "  ");
	JsonErr(js, "V3_staged_vs_flat", st, ",\n");

	// ---- V4 EvalImplicit 대 AssembleImplicit · V6 · V6b · V7 ----------------
	Stat eA, eB, eF;
	Cmp p1, pPos, pNeg;
	int nPosComp = 0, nNegComp = 0, nIndefE = 0, nIndefT = 0;
	double tMin = 1e30, gnMin = 1e30;
	Stat diffTE, diffGE, eigE;
	const double zero[N] = { 0. };
	for (std::size_t s = 0; s < ts.size(); ++s) {
		const Sample& t = ts[s];
		const double sa = CScale(t.ev);
		double C[NC], g[NC * N], H[NC * NH];
		S::Eval(t.x, t.ev, C, g, H);
		// V4
		sx::ImplicitBlock<N> blk;
		sx::AssembleImplicit<N>(NC, C, g, H, t.v, kK, kKd, kH, blk);
		double Ce[NC], fe[N], Ae[NH], be[N];
		S::EvalImplicit(t.x, t.ev, t.v, kPrm, Ce, fe, Ae, be);
		eA.Add(Rel(MaxAbsDiff(Ae, blk.A, NH), MaxAbs(blk.A, NH)));
		eB.Add(Rel(MaxAbsDiff(be, blk.b, N), MaxAbs(blk.b, N)));
		eF.Add(Rel(MaxAbsDiff(fe, blk.f, N), MaxAbs(blk.f, N)));
		// V6
		sx::ImplicitBlock<N> refT;
		Reference(C, g, H, t.v, refT);
		double Ct[NC], ft[N], At[NH], bt[N];
		S::EvalImplicitTension(t.x, t.ev, t.v, kPrm, Ct, ft, At, bt);
		Compare(C, Ct, sa, refT.A, At, refT.b, bt, refT.f, ft, p1);
		bool allPos = true, allNeg = true;
		for (int c = 0; c < NC; ++c) {
			if (C[c] > 0.) { ++nPosComp; allNeg = false; } else { ++nNegComp; allPos = false; }
		}
		double Cg[NC], fg[N], Ag[NH], bg[N];
		S::EvalImplicitGN(t.x, t.ev, t.v, kPrm, Cg, fg, Ag, bg);
		if (allPos) {   // V6b (i): v = 0 이면 s = kC 라 정확판과 같다
			double C0[NC], f0[N], A0[NH], b0[N], C1[NC], f1[N], A1[NH], b1[N];
			S::EvalImplicit(t.x, t.ev, zero, kPrm, C0, f0, A0, b0);
			S::EvalImplicitTension(t.x, t.ev, zero, kPrm, C1, f1, A1, b1);
			Compare(C0, C1, sa, A0, A1, b0, b1, f0, f1, pPos);
		}
		if (allNeg) Compare(Cg, Ct, sa, Ag, At, bg, bt, fg, ft, pNeg);   // V6b (ii): 눌린 성분만이면 GN 판과 같다
		// V7 과 보고
		const double rT = MinEigRatioUpper(At);
		tMin = std::min(tMin, rT);
		if (rT < kPsd) ++nIndefT;
		const double rE = MinEigRatioUpper(Ae);
		if (rE < kPsd) { ++nIndefE; eigE.Add(-rE); }
		gnMin = std::min(gnMin, MinEigRatioUpper(Ag));
		diffTE.Add(Rel(MaxAbsDiff(At, Ae, NH), MaxAbs(Ae, NH)));
		diffGE.Add(Rel(MaxAbsDiff(Ag, Ae, NH), MaxAbs(Ae, NH)));
	}
	std::printf("V4 EvalImplicit 대 AssembleImplicit (%zu 삼각형, k = %g, kd = %g, h = %g)\n", ts.size(), kK, kKd, kH);
	rpt.Gate("A", eA.Max() <= kLimit, "rel max %.2e", eA.Max(), kLimit);
	rpt.Gate("b", eB.Max() <= kLimit, "rel max %.2e", eB.Max(), kLimit);
	rpt.Gate("f+d", eF.Max() <= kLimit, "rel max %.2e", eF.Max(), kLimit);
	std::fprintf(js, "  \"V4_implicit\": {\"A_rel\": %.3e, \"b_rel\": %.3e, \"f_rel\": %.3e, \"A_rel_median\": %.3e, \"b_rel_median\": %.3e, \"f_rel_median\": %.3e},\n",
		eA.Max(), eB.Max(), eF.Max(), eA.Median(), eB.Median(), eF.Median());

	// ---- V5 SXMesh 삼각형 루프 (요소 값 변수 둘) ----------------------------
	{
		std::mt19937_64 rng(7);
		const int G = 200;
		const Grid grid = MakeGrid(G, rng);
		const int nNode = grid.nNode, nTri = grid.nTri;
		const double r2 = std::sqrt(2.0);
		std::vector<double> EV(static_cast<std::size_t>(NE) * nTri);   // 쉬는 uv = 평평한 격자(간격 1 cm), 대각선 √2 √2
		for (int s = 0; s < nTri; ++s) {
			for (int a = 0; a < 3; ++a) {
				const int nd = grid.tris[3 * s + a];
				EV[NE * s + 2 * a] = 0.01 * (nd % G);
				EV[NE * s + 2 * a + 1] = 0.01 * (nd / G);
			}
			EV[NE * s + 6] = r2; EV[NE * s + 7] = r2;
		}
		SXMesh mesh;
		sx::CreateVariables<S>(mesh, nNode, nTri);
		SXNodeNamedContainer& px = mesh.GetNodeV(S::PositionVar());
		SXNodeNamedContainer& pv = mesh.GetNodeV(S::VelocityVar());
		SXElementNamedContainer& pt = mesh.GetElementV(S::StencilVar());
		SXNodeNamedContainer& puv = mesh.GetNodeV(S::ElemVar(0));
		SXNodeNamedContainer& pdg = mesh.GetNodeV(S::ElemVar(1));
		for (int n = 0; n < nNode; ++n) { px.Add(n, &grid.X[3 * n], 3); pv.Add(n, &grid.V[3 * n], 3); }
		for (int s = 0; s < nTri; ++s) {
			pt.Add(s, &grid.tris[3 * s], 3);
			for (int c = 0; c < 6; ++c) puv[s][c] = EV[NE * s + c];
			for (int c = 0; c < 2; ++c) pdg[s][c] = EV[NE * s + 6 + c];
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
		std::printf("V5 SXMesh 삼각형 루프 (노드 %d, 삼각형 %d, 요소 값 변수 %d 개)\n", nNode, nTri, static_cast<int>(S::kElemVars));
		rpt.Gate("루프", mismatch == 0, "불일치 %.0f", static_cast<double>(mismatch), 0.);
		std::fprintf(js, "  \"V5_stencil_loop\": {\"nodes\": %d, \"triangles\": %d, \"mismatch\": %d},\n", nNode, nTri, mismatch);
	}

	// ---- V6 · V6b · V7 판정 --------------------------------------------------
	GateCmp(rpt, "V6 EvalImplicitTension 대 기준 조립(∇²C = 0 + h² k Σ C⁺ ∇²C)", p1, true);
	std::printf(" V6 갈래 — 성분 C > 0 %d 개 · C ≤ 0 %d 개\n", nPosComp, nNegComp);
	rpt.Gate("C>0", nPosComp >= kMinBranch, "%.0f 개", nPosComp, kMinBranch);
	rpt.Gate("C≤0", nNegComp >= kMinBranch, "%.0f 개", nNegComp, kMinBranch);
	GateCmp(rpt, "V6b (i) 두 성분 모두 C > 0 · v = 0 — EvalImplicit 와", pPos, false);
	GateCmp(rpt, "V6b (ii) 두 성분 모두 C ≤ 0 — EvalImplicitGN 과", pNeg, false);
	std::printf(" V7 준정부호 (%zu 개)\n", ts.size());
	rpt.Gate("A_T", tMin >= kPsd, "최소 고윳값 / 최대 크기의 최솟값 %.2e", tMin, kPsd);

	// 보고
	double Crest[NC], grest[NC * N], Hrest[NC * NH];
	S::Eval(ts[0].x, ts[0].ev, Crest, grest, Hrest);
	const double restC = MaxAbs(Crest, NC) / CScale(ts[0].ev);
	std::printf("  보고 — 병진 |Σ∇C|/max|∇C| max %.1e, 회전 |Σ r×∇C| max %.1e, 쉬는 모양 |C|/√a %.1e, EvalFrozen ≠ Eval %d 개, |w_u ± w_v| %.3f – %.3f\n",
		sumG.Max(), torqueG.Max(), restC, frozenDiff, dMin, dMax);
	std::printf("  보고 — 정확판 A 부정치 %d 개(|최소|/최대 중앙 %.2e), A_T 부정치 %d 개, GN A 최솟값 %.2e, |A_T − A_E|/|A_E| 중앙 %.2e · max %.2e, |A_GN − A_E|/|A_E| 중앙 %.2e · max %.2e\n",
		nIndefE, eigE.Median(), nIndefT, gnMin, diffTE.Median(), diffTE.Max(), diffGE.Median(), diffGE.Max());

	std::fprintf(js, "  \"V6\": {\n");
	JsonCmp(js, "tension_vs_reference", p1, ",");
	std::fprintf(js, "    \"components_C_pos\": %d, \"components_C_nonpos\": %d, \"min_branch\": %d\n  },\n", nPosComp, nNegComp, kMinBranch);
	std::fprintf(js, "  \"V6b\": {\n");
	JsonCmp(js, "all_pos_v0_vs_implicit", pPos, ",");
	JsonCmp(js, "all_nonpos_vs_gn", pNeg, "");
	std::fprintf(js, "  },\n");
	std::fprintf(js, "  \"V7\": {\"min_eig_ratio\": %.3e, \"limit\": %g, \"indefinite\": %d},\n", tMin, kPsd, nIndefT);
	std::fprintf(js, "  \"report\": {\"sum_grad_rel_max\": %.3e, \"torque_grad_rel_max\": %.3e, \"rest_C_over_sqrt_a\": %.3e, \"frozen_differs\": %d, "
		"\"diag_len_min\": %.6f, \"diag_len_max\": %.6f, \"indefinite_A_exact\": %d, \"indefinite_A_exact_eig_ratio_median\": %.3e, "
		"\"gn_min_eig_ratio\": %.3e, \"A_T_vs_A_E_rel_median\": %.3e, \"A_T_vs_A_E_rel_max\": %.3e, \"A_GN_vs_A_E_rel_median\": %.3e, \"A_GN_vs_A_E_rel_max\": %.3e},\n",
		sumG.Max(), torqueG.Max(), restC, frozenDiff, dMin, dMax, nIndefE, eigE.Median(), gnMin, diffTE.Median(), diffTE.Max(),
		diffGE.Median(), diffGE.Max());
	std::fprintf(js, "  \"gates_failed\": %d\n}\n", rpt.fail);
	std::printf("\nProvot 전단 게이트 V1 – V7: %s (실패 %d) → results/gate_provot_shear.json\n", rpt.fail == 0 ? "GO" : "NO-GO", rpt.fail);
	return rpt.fail == 0 ? 0 : 1;
}
