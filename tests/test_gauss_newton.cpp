//---------------------------------------------------------------------------
// test_gauss_newton.cpp — Gauss–Newton 게이트(docs/GATE_GAUSS_NEWTON.md).
// N1: 다섯 명세의 EvalImplicitGN 이 AssembleImplicit(C, ∇C, ∇²C = 0) 과 같은지 판정한다(C · ∇C 는 Eval).
// 보고: 버린 몫의 크기, 정부호, 격자에서의 시간. 엔진이 있으면(/DSX_WITH_ENGINE, engine_macri.obj · SDEmysse.obj 와 링크)
// 엔진 Bridson 시간도 잰다. N2 판정은 tests/timing_bases.py gn(고정 로드 주소 8 곳), 이 exe 의 시간은 한 번의 값이다.
// N0(기존 변형 불변)은 판정문 「재현」대로 생성기 쪽에서 본다.
// 실행 디렉터리 = 리포 루트. results/hinges.txt · triangles.txt 를 읽고 results/gate_gauss_newton.json 을 쓴다.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffBend.h"
#include "../generated/BaraffBendCompact.h"
#include "../generated/BaraffStretch.h"
#include "../generated/BaraffShear.h"
#include "../generated/BridsonBend.h"
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

static const double kPrmK = 1.0, kPrmKd = 0.1, kPrmH = 0.01;   // 판정문에 고정한 계수(조립 게이트와 같음)
static const double kPrm[3] = { kPrmK, kPrmKd, kPrmH };
static const double kLimit = 1e-12;                            // N1

// Bridson 표본 가중치 w̄ = 3|e|² / (½|nA| + ½|nB|) (Bridson 게이트와 같은 식)
static double BendWeight(const double* x)
{
	const double *xh = x, *xi = x + 3, *xj = x + 6, *xk = x + 9;
	double a[3], b[3], c[3], d[3], e[3];
	for (int t = 0; t < 3; ++t) { a[t] = xi[t] - xh[t]; b[t] = xj[t] - xh[t]; c[t] = xj[t] - xk[t]; d[t] = xi[t] - xk[t]; e[t] = xj[t] - xi[t]; }
	const double nA[3] = { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
	const double nB[3] = { c[1] * d[2] - c[2] * d[1], c[2] * d[0] - c[0] * d[2], c[0] * d[1] - c[1] * d[0] };
	const double la = std::sqrt(nA[0] * nA[0] + nA[1] * nA[1] + nA[2] * nA[2]);
	const double lb = std::sqrt(nB[0] * nB[0] + nB[1] * nB[1] + nB[2] * nB[2]);
	return 3. * (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]) / (0.5 * la + 0.5 * lb);
}

// 대칭 행렬(n × n 전체, 행 우선)의 고윳값 — 순환 Jacobi. M 은 망가진다
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

// 최소 고윳값 / 최대 고윳값 크기
static double MinEigRatio(double* M, int n)
{
	double lam[12];
	SymEig(M, n, lam);
	double lo = lam[0], big = 0.;
	for (int p = 0; p < n; ++p) { lo = std::min(lo, lam[p]); big = std::max(big, std::fabs(lam[p])); }
	return big > 0. ? lo / big : 0.;
}

template <int N>
static void FullFromUpper(const double* U, double* F)
{
	for (int a = 0; a < N; ++a)
		for (int b = 0; b < N; ++b) F[a * N + b] = U[a <= b ? sx::UpperIndex<N>(a, b) : sx::UpperIndex<N>(b, a)];
}

// ---- N1 과 보고(버린 몫 · 정부호) ------------------------------------------------
struct Acc {
	Stat c, A, b, f;            // N1
	Stat dropK, dropA;          // 보고: |C∇²C| / |∇C∇Cᵀ + C∇²C|, |A − A_GN| / |A|
	Stat dropKNear;             // 보고(판정문 고정 뒤 더함): 쉬는 모양 가까이(max|C|/척도 < 0.1)의 탄성 몫
	Stat eigK, eigA;            // 보고(판정문 고정 뒤 더함): 음일 때 |최소 고윳값| / 최대 고윳값 크기
	int n, indefK, indefA;      // 보고: 음의 고윳값을 갖는 표본 수
	double gnMin;               // 보고: GN A 의 (최소 고윳값 / 최대) 최솟값
	Acc() : n(0), indefK(0), indefA(0), gnMin(1e30) {}
};

static const double kIndef = -1e-12;   // 판정문 「보고만」: 최대 고윳값 크기의 1e-12 아래면 음
static const double kNear = 0.1;       // 쉬는 모양 가까이: max|C| / 척도 < 0.1 (굽힘 0.1 rad, 인장 · 전단 변형률 0.1)

template <class Spec>
static void Check(const double* x, const double* ev, const double* v, double cscale, Acc& acc)
{
	enum { N = Spec::kDof, NH = Spec::kHess, NC = Spec::kComp };
	double C[NC], g[NC * N], H[NC * NH], H0[NC * NH];
	Spec::Eval(x, ev, C, g, H);
	for (int i = 0; i < NC * NH; ++i) H0[i] = 0.;
	sx::ImplicitBlock<N> ref;
	sx::AssembleImplicit<N>(NC, C, g, H0, v, kPrmK, kPrmKd, kPrmH, ref);
	double Cg[NC], f[N], A[NH], b[N];
	Spec::EvalImplicitGN(x, ev, v, kPrm, Cg, f, A, b);
	acc.c.Add(MaxAbsDiff(Cg, C, NC) / cscale);
	acc.A.Add(Rel(MaxAbsDiff(A, ref.A, NH), MaxAbs(ref.A, NH)));
	acc.b.Add(Rel(MaxAbsDiff(b, ref.b, N), MaxAbs(ref.b, N)));
	acc.f.Add(Rel(MaxAbsDiff(f, ref.f, N), MaxAbs(ref.f, N)));

	// 보고: 탄성 몫 −K/k = Σ (∇C∇Cᵀ + C∇²C) 와 버린 C∇²C
	double Kf[N * N], CH[N * N];
	for (int i = 0; i < N * N; ++i) { Kf[i] = 0.; CH[i] = 0.; }
	for (int c = 0; c < NC; ++c)
		for (int a = 0; a < N; ++a)
			for (int bb = 0; bb < N; ++bb) {
				const double h = H[c * NH + (a <= bb ? sx::UpperIndex<N>(a, bb) : sx::UpperIndex<N>(bb, a))];
				CH[a * N + bb] += C[c] * h;
				Kf[a * N + bb] += g[c * N + a] * g[c * N + bb] + C[c] * h;
			}
	acc.dropK.Add(Rel(MaxAbs(CH, N * N), MaxAbs(Kf, N * N)));
	if (MaxAbs(C, NC) / cscale < kNear) acc.dropKNear.Add(Rel(MaxAbs(CH, N * N), MaxAbs(Kf, N * N)));
	double Ce[NC], fe[N], Ae[NH], be[N];
	Spec::EvalImplicit(x, ev, v, kPrm, Ce, fe, Ae, be);
	acc.dropA.Add(Rel(MaxAbsDiff(Ae, A, NH), MaxAbs(Ae, NH)));

	const double rK = MinEigRatio(Kf, N);
	if (rK < kIndef) { ++acc.indefK; acc.eigK.Add(-rK); }
	double M[N * N];
	FullFromUpper<N>(Ae, M);
	const double rA = MinEigRatio(M, N);
	if (rA < kIndef) { ++acc.indefA; acc.eigA.Add(-rA); }
	FullFromUpper<N>(A, M);
	acc.gnMin = std::min(acc.gnMin, MinEigRatio(M, N));
	++acc.n;
}

static void PrintN1(Report& rpt, const char* name, const char* scale, const Acc& a)
{
	std::printf(" %s (%d 개, C 척도 %s)\n", name, a.n, scale);
	rpt.Gate("C", a.n > 0 && a.c.Max() <= kLimit, "|ΔC|/척도 %.2e", a.c.Max(), kLimit);
	rpt.Gate("A", a.n > 0 && a.A.Max() <= kLimit, "rel max %.2e", a.A.Max(), kLimit);
	rpt.Gate("b", a.n > 0 && a.b.Max() <= kLimit, "rel max %.2e", a.b.Max(), kLimit);
	rpt.Gate("f+d", a.n > 0 && a.f.Max() <= kLimit, "rel max %.2e", a.f.Max(), kLimit);
	std::printf("  보고 — 버린 몫: 탄성 |C∇²C|/|−K/k| 중앙 %.2e · max %.2e, 조립 |A − A_GN|/|A| 중앙 %.2e · max %.2e\n",
		a.dropK.Median(), a.dropK.Max(), a.dropA.Median(), a.dropA.Max());
	std::printf("  보고 — 쉬는 모양 가까이(|C|/척도 < %.1f, %zu 개): 탄성 |C∇²C|/|−K/k| 중앙 %.2e · max %.2e\n",
		kNear, a.dropKNear.v.size(), a.dropKNear.Median(), a.dropKNear.Max());
	std::printf("  보고 — 음의 고윳값: −K/k %d 개(|최소|/최대 중앙 %.2e), 정확판 A %d 개(%.2e) (%d 개 중); GN A 최소/최대 고윳값 최솟값 %.1e\n",
		a.indefK, a.eigK.Median(), a.indefA, a.eigA.Median(), a.n, a.gnMin);
}

static void JsonN1(std::FILE* js, const char* name, const Acc& a, const char* tail)
{
	std::fprintf(js, "    \"%s\": {\"samples\": %d, \"C_abs_over_scale\": %.3e, \"A_rel\": %.3e, \"b_rel\": %.3e, \"f_rel\": %.3e, "
		"\"A_rel_median\": %.3e, \"b_rel_median\": %.3e, \"f_rel_median\": %.3e, "
		"\"drop_K_median\": %.3e, \"drop_K_max\": %.3e, \"drop_A_median\": %.3e, \"drop_A_max\": %.3e, "
		"\"near_rest_samples\": %zu, \"drop_K_near_median\": %.3e, \"drop_K_near_max\": %.3e, "
		"\"indefinite_K\": %d, \"indefinite_K_eig_ratio_median\": %.3e, \"indefinite_A_exact\": %d, \"indefinite_A_eig_ratio_median\": %.3e, "
		"\"gn_min_eig_ratio\": %.3e}%s\n",
		name, a.n, a.c.Max(), a.A.Max(), a.b.Max(), a.f.Max(), a.A.Median(), a.b.Median(), a.f.Median(),
		a.dropK.Median(), a.dropK.Max(), a.dropA.Median(), a.dropA.Max(),
		a.dropKNear.v.size(), a.dropKNear.Median(), a.dropKNear.Max(),
		a.indefK, a.eigK.Median(), a.indefA, a.eigA.Median(), a.gnMin, tail);
}

// ---- 시간 --------------------------------------------------------------------
typedef std::chrono::steady_clock Clock;

struct UseExact {
	template <class S> static void Run(const double* x, const double* ev, const double* v, double* o)
	{ S::EvalImplicit(x, ev, v, kPrm, o, o + S::kComp, o + S::kComp + S::kDof, o + S::kComp + S::kDof + S::kHess); }
};
struct UseGN {
	template <class S> static void Run(const double* x, const double* ev, const double* v, double* o)
	{ S::EvalImplicitGN(x, ev, v, kPrm, o, o + S::kComp, o + S::kComp + S::kDof, o + S::kComp + S::kDof + S::kHess); }
};

// 모으기 + 변형 + 출력(C · f + d · A · b)을 스텐실별 배열에 저장 (N2 의 생성 쪽 = Bridson · UseGN)
template <class Spec, class Var>
static double TimeStore(const sx::StencilView<Spec>& view, std::vector<double>& out, double& checksum)
{
	enum { N = Spec::kDof, NE = Spec::kElemValues, NO = Spec::kComp + 2 * Spec::kDof + Spec::kHess };
	const int n = view.size();
	const Clock::time_point t0 = Clock::now();
	for (int s = 0; s < n; ++s) {
		double x[N], v[N], ev[NE];
		view.Gather(s, x, v, ev);
		Var::template Run<Spec>(x, ev, v, &out[static_cast<std::size_t>(s) * NO]);
	}
	const double ns = std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / n;
	for (int s = 0; s < n; s += 97) checksum += out[static_cast<std::size_t>(s) * NO + 20];
	return ns;
}

// 보고: 모으기 + 변형 + f · b 를 노드 배열로 흩뿌리기
template <class Spec, class Var>
static double TimeScatter(const sx::StencilView<Spec>& view, std::vector<double>& Fn, std::vector<double>& Bn, double& checksum)
{
	enum { N = Spec::kDof, NE = Spec::kElemValues, NC = Spec::kComp, NO = Spec::kComp + 2 * Spec::kDof + Spec::kHess };
	std::fill(Fn.begin(), Fn.end(), 0.);
	std::fill(Bn.begin(), Bn.end(), 0.);
	const int n = view.size();
	const Clock::time_point t0 = Clock::now();
	for (int s = 0; s < n; ++s) {
		double x[N], v[N], ev[NE], o[NO];
		view.Gather(s, x, v, ev);
		Var::template Run<Spec>(x, ev, v, o);
		const int* nd = view.Nodes(s);
		for (int a = 0; a < Spec::kNodes; ++a)
			for (int c = 0; c < 3; ++c) { Fn[3 * nd[a] + c] += o[NC + 3 * a + c]; Bn[3 * nd[a] + c] += o[NC + N + Spec::kHess + 3 * a + c]; }
		checksum += o[NC + N + 13];
	}
	return std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / n;
}

// 삼각형 명세의 SXMesh(쉬는 uv = 평평한 격자, 간격 1 cm — 인장 · 전단 게이트 T5 와 같다)
template <class Spec>
static void FillTriMesh(SXMesh& mesh, const Grid& grid, int G)
{
	sx::CreateVariables<Spec>(mesh, grid.nNode, grid.nTri);
	SXNodeNamedContainer& px = mesh.GetNodeV(Spec::PositionVar());
	SXNodeNamedContainer& pv = mesh.GetNodeV(Spec::VelocityVar());
	SXElementNamedContainer& pt = mesh.GetElementV(Spec::StencilVar());
	SXNodeNamedContainer& puv = mesh.GetNodeV(Spec::ElemVar(0));
	for (int n = 0; n < grid.nNode; ++n) { px.Add(n, &grid.X[3 * n], 3); pv.Add(n, &grid.V[3 * n], 3); }
	for (int s = 0; s < grid.nTri; ++s) {
		pt.Add(s, &grid.tris[3 * s], 3);
		for (int a = 0; a < 3; ++a) {
			const int nd = grid.tris[3 * s + a];
			puv[s][2 * a] = 0.01 * (nd % G);
			puv[s][2 * a + 1] = 0.01 * (nd / G);
		}
	}
}

struct Timing {
	double bridsonGN, bridsonExact, bridsonGNScatter, baraffGN, baraffExact, compactGN, compactExact,
		stretchGN, stretchExact, shearGN, shearExact, engine, engineSerial;
	int nHinge, nTri;
};

#ifdef SX_WITH_ENGINE
static double TimeEngine(engine_ref::Bridson* eb, double& checksum)
{
	const Clock::time_point t0 = Clock::now();
	engine_ref::RunBridson(eb);
	const double ns = std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / engine_ref::SizeBridson(eb);
	checksum += engine_ref::ChecksumBridson(eb);
	return ns;
}
#endif

// full = false 이면 Bridson 의 GN · 엔진만(작은 격자 보고)
static Timing TimeGrid(int G, bool full)
{
	std::mt19937_64 rng(7);
	const Grid grid = MakeGrid(G, rng);
	std::vector<double> W(grid.nHinge), L0(grid.nHinge), Z(grid.nHinge, 0.);
	for (int s = 0; s < grid.nHinge; ++s) {   // 격자 모양을 쉬는 모양으로 본다(Bridson 게이트 B8 과 같다)
		double x[12];
		for (int a = 0; a < 4; ++a) for (int c = 0; c < 3; ++c) x[3 * a + c] = grid.X[3 * grid.hinges[4 * s + a] + c];
		W[s] = BendWeight(x);
		L0[s] = EdgeLength(x);
	}
	SXMesh mBridson, mBaraff, mCompact, mStretch, mShear;
	FillMesh<sxgen::BridsonBend>(mBridson, grid, W);
	const sx::StencilView<sxgen::BridsonBend> vBridson(mBridson);
	std::vector<double> out(static_cast<std::size_t>(std::max(grid.nHinge, grid.nTri)) * 103, 0.);
	std::vector<double> Fn(3 * grid.nNode), Bn(3 * grid.nNode);
#ifdef SX_WITH_ENGINE
	std::vector<double> X(grid.X), V(grid.V);
	engine_ref::Bridson* eb = engine_ref::CreateBridson(grid.nHinge, &grid.hinges[0], grid.nNode, &X[0], &V[0], &L0[0]);
#endif
	Timing t;
	double* all[] = { &t.bridsonGN, &t.bridsonExact, &t.bridsonGNScatter, &t.baraffGN, &t.baraffExact, &t.compactGN, &t.compactExact,
		&t.stretchGN, &t.stretchExact, &t.shearGN, &t.shearExact, &t.engine, &t.engineSerial };
	for (std::size_t i = 0; i < sizeof all / sizeof all[0]; ++i) *all[i] = 1e30;
	t.nHinge = grid.nHinge;
	t.nTri = grid.nTri;
	double checksum = 0.;
	for (int rep = 0; rep < 5; ++rep) {   // 변형을 번갈아 5 회, 변형마다 최소
#ifdef SX_WITH_ENGINE
		engine_ref::SetSerialBridson(eb, false);
		t.engine = std::min(t.engine, TimeEngine(eb, checksum));
#endif
		t.bridsonGN = std::min(t.bridsonGN, TimeStore<sxgen::BridsonBend, UseGN>(vBridson, out, checksum));
		if (!full) continue;
#ifdef SX_WITH_ENGINE
		engine_ref::SetSerialBridson(eb, true);
		t.engineSerial = std::min(t.engineSerial, TimeEngine(eb, checksum));
#endif
		t.bridsonExact = std::min(t.bridsonExact, TimeStore<sxgen::BridsonBend, UseExact>(vBridson, out, checksum));
		t.bridsonGNScatter = std::min(t.bridsonGNScatter, TimeScatter<sxgen::BridsonBend, UseGN>(vBridson, Fn, Bn, checksum));
	}
#ifdef SX_WITH_ENGINE
	engine_ref::DestroyBridson(eb);
#endif
	if (full) {   // 다른 명세: 엔진이 없으니 Bridson 시간과 섞지 않고 따로 번갈아 잰다
		FillMesh<sxgen::BaraffBend>(mBaraff, grid, Z);
		FillMesh<sxgen::BaraffBendCompact>(mCompact, grid, Z);
		FillTriMesh<sxgen::BaraffStretch>(mStretch, grid, G);
		FillTriMesh<sxgen::BaraffShear>(mShear, grid, G);
		const sx::StencilView<sxgen::BaraffBend> vBaraff(mBaraff);
		const sx::StencilView<sxgen::BaraffBendCompact> vCompact(mCompact);
		const sx::StencilView<sxgen::BaraffStretch> vStretch(mStretch);
		const sx::StencilView<sxgen::BaraffShear> vShear(mShear);
		for (int rep = 0; rep < 5; ++rep) {
			t.baraffGN     = std::min(t.baraffGN,     TimeStore<sxgen::BaraffBend, UseGN>(vBaraff, out, checksum));
			t.baraffExact  = std::min(t.baraffExact,  TimeStore<sxgen::BaraffBend, UseExact>(vBaraff, out, checksum));
			t.compactGN    = std::min(t.compactGN,    TimeStore<sxgen::BaraffBendCompact, UseGN>(vCompact, out, checksum));
			t.compactExact = std::min(t.compactExact, TimeStore<sxgen::BaraffBendCompact, UseExact>(vCompact, out, checksum));
			t.stretchGN    = std::min(t.stretchGN,    TimeStore<sxgen::BaraffStretch, UseGN>(vStretch, out, checksum));
			t.stretchExact = std::min(t.stretchExact, TimeStore<sxgen::BaraffStretch, UseExact>(vStretch, out, checksum));
			t.shearGN      = std::min(t.shearGN,      TimeStore<sxgen::BaraffShear, UseGN>(vShear, out, checksum));
			t.shearExact   = std::min(t.shearExact,   TimeStore<sxgen::BaraffShear, UseExact>(vShear, out, checksum));
		}
	}
	std::printf("  (격자 G = %d, checksum %.3e)\n", G, checksum);
	return t;
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Hinge> hs = ReadHinges("results/hinges.txt");
	const std::vector<Tri> ts = ReadTriangles("results/triangles.txt");
	if (hs.empty() || ts.empty()) { std::printf("results/hinges.txt · triangles.txt 없음\n"); return 2; }
	Report rpt("results/gate_gauss_newton.json");
	std::FILE* js = rpt.json;
	std::printf("Gauss–Newton 게이트 — 힌지 %zu 개, 삼각형 %zu 개, k = %g, kd = %g, h = %g\n", hs.size(), ts.size(), kPrmK, kPrmKd, kPrmH);
	std::fprintf(js, "{\n  \"hinges\": %zu,\n  \"triangles\": %zu,\n  \"k\": %g, \"kd\": %g, \"h\": %g,\n", hs.size(), ts.size(), kPrmK, kPrmKd, kPrmH);

	// ---- N1 ---------------------------------------------------------------------
	Acc aBaraff, aCompact, aStretch, aShear, aBridson;
	for (std::size_t s = 0; s < hs.size(); ++s) {
		const Hinge& h = hs[s];
		Check<sxgen::BaraffBend>(h.x, &h.theta0, h.v, 1., aBaraff);
		Check<sxgen::BaraffBendCompact>(h.x, &h.theta0, h.v, 1., aCompact);
		const double w = BendWeight(h.x);
		Check<sxgen::BridsonBend>(h.x, &w, h.v, std::sqrt(w), aBridson);
	}
	for (std::size_t s = 0; s < ts.size(); ++s) {
		const Tri& t = ts[s];
		const double sa = std::sqrt(UvArea(t.uv));
		Check<sxgen::BaraffStretch>(t.x, t.uv, t.v, sa, aStretch);
		Check<sxgen::BaraffShear>(t.x, t.uv, t.v, sa, aShear);
	}
	std::printf("N1 EvalImplicitGN 대 AssembleImplicit(C, ∇C, ∇²C = 0)\n");
	PrintN1(rpt, "BaraffBend", "1 rad", aBaraff);
	PrintN1(rpt, "BaraffBendCompact", "1 rad", aCompact);
	PrintN1(rpt, "BaraffStretch", "√a", aStretch);
	PrintN1(rpt, "BaraffShear", "√a", aShear);
	PrintN1(rpt, "BridsonBend", "√w̄", aBridson);
	std::fprintf(js, "  \"N1\": {\n");
	JsonN1(js, "BaraffBend", aBaraff, ",");
	JsonN1(js, "BaraffBendCompact", aCompact, ",");
	JsonN1(js, "BaraffStretch", aStretch, ",");
	JsonN1(js, "BaraffShear", aShear, ",");
	JsonN1(js, "BridsonBend", aBridson, "");
	std::fprintf(js, "  },\n");

	// ---- 시간 (N2 는 timing_bases.py gn) --------------------------------------------
	const Timing t = TimeGrid(200, true);
	const Timing ts60 = TimeGrid(60, false);
	std::printf("시간 (이 exe 한 번 — N2 판정은 tests/timing_bases.py gn, 힌지 %d · 삼각형 %d, ns / 스텐실)\n", t.nHinge, t.nTri);
	std::printf("  Bridson   GN + 저장 %.0f, 정확판 + 저장 %.0f, GN + 흩뿌리기 %.0f\n", t.bridsonGN, t.bridsonExact, t.bridsonGNScatter);
#ifdef SX_WITH_ENGINE
	std::printf("  엔진 Bridson %.0f (기본값 환원 포함 %.0f) → GN / 엔진 %.3f; 작은 격자 GN %.0f / 엔진 %.0f\n",
		t.engine, t.engineSerial, t.bridsonGN / t.engine, ts60.bridsonGN, ts60.engine);
#endif
	std::printf("  Baraff 굽힘 GN %.0f / 정확판 %.0f, 압축 명세 GN %.0f / %.0f, 인장 GN %.0f / %.0f, 전단 GN %.0f / %.0f\n",
		t.baraffGN, t.baraffExact, t.compactGN, t.compactExact, t.stretchGN, t.stretchExact, t.shearGN, t.shearExact);
	std::fprintf(js, "  \"grid\": {\"hinges\": %d, \"triangles\": %d},\n", t.nHinge, t.nTri);
	std::fprintf(js, "  \"ns_per_stencil\": {\"bridson_gn_store\": %.1f, \"bridson_implicit_store\": %.1f, \"bridson_gn_scatter\": %.1f, "
		"\"baraff_gn_store\": %.1f, \"baraff_implicit_store\": %.1f, \"compact_gn_store\": %.1f, \"compact_implicit_store\": %.1f, "
		"\"stretch_gn_store\": %.1f, \"stretch_implicit_store\": %.1f, \"shear_gn_store\": %.1f, \"shear_implicit_store\": %.1f, "
		"\"small_bridson_gn_store\": %.1f",
		t.bridsonGN, t.bridsonExact, t.bridsonGNScatter, t.baraffGN, t.baraffExact, t.compactGN, t.compactExact,
		t.stretchGN, t.stretchExact, t.shearGN, t.shearExact, ts60.bridsonGN);
#ifdef SX_WITH_ENGINE
	std::fprintf(js, ", \"engine_oh\": %.1f, \"engine_oh_serial\": %.1f, \"small_engine_oh\": %.1f", t.engine, t.engineSerial, ts60.engine);
#endif
	std::fprintf(js, "},\n");
	std::fprintf(js, "  \"gates_failed\": %d\n}\n", rpt.fail);
	std::printf("\nN1 %s (실패 %d)\n", rpt.fail == 0 ? "PASS — 생성 GO 의 조건" : "FAIL", rpt.fail);
	return rpt.fail == 0 ? 0 : 1;
}
