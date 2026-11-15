//---------------------------------------------------------------------------
// test_poisson_integration.cpp — Poisson 힘 단축 인장 적분(docs/GATE_POISSON_INTEGRATION.md).
// 늘려 놓은 작은 천을 Baraff–Witkin 식 (6)(단계마다 한 번, 거르개 PCG)으로 적분해 가로 변형이 따로 구한 평형으로 가는지 본다.
// 힘 = 인장(EvalImplicitTension) + 전단(EvalImplicitGN) + Poisson(정확판 EvalImplicit 또는 EvalImplicitGN). 굽힘 · 중력 없음.
//   U0a 정적 균형 · U0b 전역 조립 대 유한차분(전제) · 장면 U(양 끝의 x 만 고정, 판정) · 장면 C(양 끝 고정, 보고)
//   → results/gate_poisson_integration.json
// 조립 · 거르개 PCG 의 짜임은 tests/test_integration.cpp 와 같다(고정을 노드가 아니라 자유도마다 둔다).
// 2 차(판정문 「2 차 기준」): 장면 U 의 「정상」은 노드 평균 속도(y · z)를 뺀 속력으로 본다 — 풀린 강체 병진이 남기 때문이다.
// 실행 디렉터리 = 리포 루트. results/mma_ref_poisson_uniaxial.txt(model/sul_poisson_uniaxial_ref.wl)를 읽는다.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffStretch.h"
#include "../generated/BaraffShear.h"
#include "../generated/SulPoisson.h"
#include "gate_common.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;

// ---- 시나리오 (판정문에 고정) ------------------------------------------------------
static const int    kNU = 21, kNV = 11;            // 노드 수 (u × v)
static const double kEll = 0.005;                  // 노드 간격 5 mm → 10 cm × 5 cm
static const double kJitter = 0.2;                 // 안쪽 노드의 쉬는 uv 흩뜨림 (간격의 ±20 %)
static const unsigned long long kSeed = 20261116ULL;
static const double kRho = 0.2;                    // 면밀도 (kg/m²)
static const double kStretchK = 1000., kStretchKd = 0.1;
static const double kShearK = 100., kShearKd = 0.01;
static const double kPoissonKdRatio = 1e-4;        // kd_P = 1e-4 k_P (인장의 kd / k 와 같다)
enum { kNR = 2, kNB = 4, kNS = 3, kNH = 2, kNCase = kNR * kNB * kNS };
static const double kR[kNR] = { 0.3, 1.0 };                 // k_P / k_T
static const double kBeta[kNB] = { -0.6, 0.0, 0.6, 1.0 };
static const double kStrain[kNS] = { 0.02, 0.2, 1.0 };      // ε (u 방향 늘림)
static const int    kSps[kNH] = { 120, 30 };                // h = 1/120 · 1/30
static const double kT = 2.0;
static const double kTail = 0.5;                   // 「정상」을 보는 마지막 구간 (s)
static const double kPcgEps = 1e-6;                // PCG δ ≤ ε² δ₀
static const int    kCap = 2000;                   // PCG 반복 상한
static const double kVBlow = 50.;                  // 발산 속력 (m/s)
static const double kSteadySpeed = 1e-8;           // 정상: 최대 노드 속력 (m/s)
static const double kReach = 1e-6;                 // 도달: |w_u| · |w_v| · ŵ_u·ŵ_v 의 오차
static const double kU0aLimit = 1e-10, kU0bLimit = 1e-6;

enum { kSceneU = 0, kSceneC = 1 };   // U = 양 끝의 x 만 고정, C = 양 끝을 모두 고정
enum { kModeE = 0, kModeG = 1 };     // Poisson 의 강성 판
static const char* kModeName[2] = { "E", "G" };

// ---- 천: 위상 · 쉬는 uv · 블록 행렬 틀 -----------------------------------------------
struct Sheet {
	int nNode, nTri, nPair;
	std::vector<double> U, W, mass;     // 노드의 쉬는 좌표 (u, v) · 질량
	std::vector<int> tris;              // 삼각형(반시계) 3 개
	std::vector<double> uv;             // 삼각형마다 쉬는 uv 6 칸
	std::vector<int> pairI, pairJ;      // 노드쌍 블록 (i < j)
	std::vector<int> triSlot;           // 삼각형의 지역 쌍 (a < b) → 노드쌍 블록 번호
};

inline int LocalPair3(int a, int b) { return a * 3 - a * (a + 1) / 2 + (b - a - 1); }   // (0,1)(0,2)(1,2)

static Sheet MakeSheet()
{
	Sheet c;
	c.nNode = kNU * kNV;
	c.U.resize(c.nNode); c.W.resize(c.nNode);
	std::mt19937_64 rng(kSeed);
	std::uniform_real_distribution<double> jit(-1., 1.);
	for (int j = 0; j < kNV; ++j)
		for (int i = 0; i < kNU; ++i) {
			const int n = j * kNU + i;
			double u = kEll * i, v = kEll * j;
			if (i > 0 && i + 1 < kNU && j > 0 && j + 1 < kNV) {
				u += kJitter * kEll * jit(rng);
				v += kJitter * kEll * jit(rng);
			}
			c.U[n] = u; c.W[n] = v;
		}
	for (int j = 0; j + 1 < kNV; ++j)
		for (int i = 0; i + 1 < kNU; ++i) {
			const int a = j * kNU + i, b = a + 1, cc = a + kNU + 1, d = a + kNU;
			const int t1[3] = { a, b, cc }, t2[3] = { a, cc, d };
			c.tris.insert(c.tris.end(), t1, t1 + 3);
			c.tris.insert(c.tris.end(), t2, t2 + 3);
		}
	c.nTri = static_cast<int>(c.tris.size() / 3);
	c.mass.assign(c.nNode, 0.);
	c.uv.resize(6 * static_cast<std::size_t>(c.nTri));
	for (int t = 0; t < c.nTri; ++t) {
		for (int a = 0; a < 3; ++a) { const int n = c.tris[3 * t + a]; c.uv[6 * t + 2 * a] = c.U[n]; c.uv[6 * t + 2 * a + 1] = c.W[n]; }
		const double area = UvArea(&c.uv[6 * t]);
		for (int a = 0; a < 3; ++a) c.mass[c.tris[3 * t + a]] += kRho * area / 3.;
	}
	std::unordered_map<long long, int> slot;
	c.triSlot.resize(3 * static_cast<std::size_t>(c.nTri));
	for (int t = 0; t < c.nTri; ++t)
		for (int a = 0; a < 3; ++a) for (int b = a + 1; b < 3; ++b) {
			int p = c.tris[3 * t + a], q = c.tris[3 * t + b];
			if (p > q) { const int tmp = p; p = q; q = tmp; }
			const long long key = static_cast<long long>(p) * c.nNode + q;
			std::unordered_map<long long, int>::const_iterator it = slot.find(key);
			int s;
			if (it != slot.end()) s = it->second;
			else { s = static_cast<int>(c.pairI.size()); slot[key] = s; c.pairI.push_back(p); c.pairJ.push_back(q); }
			c.triSlot[3 * t + LocalPair3(a, b)] = s;
		}
	c.nPair = static_cast<int>(c.pairI.size());
	return c;
}

// 자유도마다의 고정: 양 끝(u = 0, u = L) 노드의 x (장면 C 는 y · z 도)
static std::vector<char> FixedMask(const Sheet& c, int scene)
{
	std::vector<char> f(3 * static_cast<std::size_t>(c.nNode), 0);
	for (int j = 0; j < kNV; ++j)
		for (int i = 0; i < kNU; i += kNU - 1) {
			const int n = j * kNU + i;
			f[3 * n] = 1;
			if (scene == kSceneC) { f[3 * n + 1] = 1; f[3 * n + 2] = 1; }
		}
	return f;
}

// 전체 행렬: 노드 대각 블록 D[9n] + 노드쌍 블록 O[9p](행 = pairI, 열 = pairJ, 아래쪽은 전치), 3 × 3 행 우선
struct BlockMat {
	std::vector<double> D, O;
	void Resize(const Sheet& c) { D.assign(9 * static_cast<std::size_t>(c.nNode), 0.); O.assign(9 * static_cast<std::size_t>(c.nPair), 0.); }
	void Zero() { std::fill(D.begin(), D.end(), 0.); std::fill(O.begin(), O.end(), 0.); }
};

// 삼각형 요소 A(위 삼각 행 우선, 9 자유도)를 전체 행렬에 더한다
inline void ScatterA(const double* A, const int* nd, const int* slot, BlockMat& M)
{
	int idx = 0;
	for (int p = 0; p < 9; ++p) {
		const int a = p / 3, s = p % 3;
		for (int q = p; q < 9; ++q, ++idx) {
			const int b = q / 3, t = q % 3;
			const double val = A[idx];
			if (a == b) {
				double* Db = &M.D[9 * static_cast<std::size_t>(nd[a])];
				Db[3 * s + t] += val;
				if (s != t) Db[3 * t + s] += val;
			}
			else {
				double* Ob = &M.O[9 * static_cast<std::size_t>(slot[LocalPair3(a, b)])];
				if (nd[a] < nd[b]) Ob[3 * s + t] += val; else Ob[3 * t + s] += val;
			}
		}
	}
}

static void MatVec(const Sheet& c, const BlockMat& M, const double* x, double* y)
{
	for (int n = 0; n < c.nNode; ++n) {
		const double* B = &M.D[9 * static_cast<std::size_t>(n)];
		const double* xn = x + 3 * n;
		for (int s = 0; s < 3; ++s) y[3 * n + s] = B[3 * s] * xn[0] + B[3 * s + 1] * xn[1] + B[3 * s + 2] * xn[2];
	}
	for (int p = 0; p < c.nPair; ++p) {
		const double* B = &M.O[9 * static_cast<std::size_t>(p)];
		const int i = c.pairI[p], j = c.pairJ[p];
		const double *xi = x + 3 * i, *xj = x + 3 * j;
		double* yi = y + 3 * i;
		double* yj = y + 3 * j;
		for (int s = 0; s < 3; ++s) {
			yi[s] += B[3 * s] * xj[0] + B[3 * s + 1] * xj[1] + B[3 * s + 2] * xj[2];
			yj[s] += B[s] * xi[0] + B[3 + s] * xi[1] + B[6 + s] * xi[2];
		}
	}
}

// ---- PCG (Baraff–Witkin 거르개, Jacobi 전처리) -----------------------------------------
struct PcgOut { int iters; bool converged, breakdown, precond; };
struct PcgWork { std::vector<double> r, c, q, s, inv; };

static PcgOut Pcg(const Sheet& cl, const BlockMat& M, const std::vector<double>& b, std::vector<double>& x, PcgWork& w,
	const std::vector<char>& fixed)
{
	const int n3 = 3 * cl.nNode;
	PcgOut o = { 0, false, false, false };
	w.r.resize(n3); w.c.resize(n3); w.q.resize(n3); w.s.resize(n3); w.inv.resize(n3);
	for (int n = 0; n < cl.nNode; ++n)
		for (int t = 0; t < 3; ++t) {
			const int d = 3 * n + t;
			if (fixed[d]) { w.inv[d] = 0.; continue; }
			const double diag = M.D[9 * static_cast<std::size_t>(n) + 4 * t];
			if (diag > 0.) w.inv[d] = 1. / diag;
			else { w.inv[d] = 1. / cl.mass[n]; o.precond = true; }
		}
	double d0 = 0., dn = 0.;
	for (int d = 0; d < n3; ++d) {
		x[d] = 0.;
		w.r[d] = w.inv[d] != 0. ? b[d] : 0.;   // 거르개: 고정 자유도는 inv = 0
		w.c[d] = w.inv[d] * w.r[d];
		d0 += w.r[d] * w.c[d];
		dn += w.r[d] * w.c[d];
	}
	const double stop = kPcgEps * kPcgEps * d0;
	while (dn > stop && o.iters < kCap) {
		MatVec(cl, M, &w.c[0], &w.q[0]);
		double pAp = 0.;
		for (int d = 0; d < n3; ++d) { if (w.inv[d] == 0.) w.q[d] = 0.; pAp += w.c[d] * w.q[d]; }
		if (!(pAp > 0.) || !std::isfinite(pAp)) { o.breakdown = true; break; }
		const double alpha = dn / pAp;
		double dnew = 0.;
		for (int d = 0; d < n3; ++d) {
			x[d] += alpha * w.c[d];
			w.r[d] -= alpha * w.q[d];
			w.s[d] = w.inv[d] * w.r[d];
			dnew += w.r[d] * w.s[d];
		}
		++o.iters;
		if (!std::isfinite(dnew)) { o.breakdown = true; break; }
		const double beta = dnew / dn;
		dn = dnew;
		for (int d = 0; d < n3; ++d) w.c[d] = w.s[d] + beta * w.c[d];
	}
	o.converged = !o.breakdown && dn <= stop;
	return o;
}

// ---- 조립 ---------------------------------------------------------------------------
typedef void (*ImplicitFn)(const double*, const double*, const double*, const double*, double*, double*, double*, double*);

struct Terms {
	ImplicitFn stretch, shear, poisson;   // 0 이면 그 힘을 뺀다
	double kP, kdT, kdS, kdP, h;
	bool mass;
};

struct Sim {
	const Sheet& c;
	SXMesh mesh;
	double* X;
	double* V;
	SXNodeNamedContainer* beta;
	BlockMat M;
	std::vector<double> F, b;
	explicit Sim(const Sheet& sh) : c(sh), X(0), V(0), beta(0)
	{
		// 명세가 선언한 변수를 모두 먼저 만들고(같은 이름은 다시 잡힌다) 채운다
		sx::CreateVariables<sxgen::BaraffStretch>(mesh, c.nNode, c.nTri);
		sx::CreateVariables<sxgen::BaraffShear>(mesh, c.nNode, c.nTri);
		sx::CreateVariables<sxgen::SulPoisson>(mesh, c.nNode, c.nTri);
		SXNodeNamedContainer& px = mesh.GetNodeV(sxgen::SulPoisson::PositionVar());
		SXElementNamedContainer& pt = mesh.GetElementV(sxgen::SulPoisson::StencilVar());
		SXNodeNamedContainer& puv = mesh.GetNodeV(sxgen::SulPoisson::ElemVar(0));
		beta = &mesh.GetNodeV(sxgen::SulPoisson::ElemVar(1));
		for (int n = 0; n < c.nNode; ++n) { const double p[3] = { c.U[n], c.W[n], 0. }; px.Add(n, p, 3); }
		for (int t = 0; t < c.nTri; ++t) {
			pt.Add(t, &c.tris[3 * t], 3);
			for (int k = 0; k < 6; ++k) puv[t][k] = c.uv[6 * t + k];
			(*beta)[t][0] = 0.;
		}
		X = px.Data();
		V = mesh.GetNodeV(sxgen::SulPoisson::VelocityVar()).Data();
		M.Resize(c);
		F.assign(3 * static_cast<std::size_t>(c.nNode), 0.);
		b.assign(3 * static_cast<std::size_t>(c.nNode), 0.);
	}
	void SetBeta(double bt) { for (int t = 0; t < c.nTri; ++t) (*beta)[t][0] = bt; }
	void SetState(double su, double sv)   // 균일 변형 x = s_u u, y = s_v v, 속도 0
	{
		for (int n = 0; n < c.nNode; ++n) { X[3 * n] = su * c.U[n]; X[3 * n + 1] = sv * c.W[n]; X[3 * n + 2] = 0.; }
		std::fill(V, V + 3 * c.nNode, 0.);
	}
};

struct Views {
	sx::StencilView<sxgen::BaraffStretch> st;
	sx::StencilView<sxgen::BaraffShear> sh;
	sx::StencilView<sxgen::SulPoisson> po;
	explicit Views(const SXMesh& m) : st(m), sh(m), po(m) {}
};

template <class Spec>
static void AddElements(const sx::StencilView<Spec>& view, const int* slots, const double* prm, ImplicitFn fn, Sim& sim)
{
	enum { N = Spec::kDof, NE = Spec::kElemValues, NC = Spec::kComp, NH = Spec::kHess };
	static_assert(Spec::kNodes == 3, "삼각형 스텐실");
	const int n = view.size();
	double* F = &sim.F[0];
	double* b = &sim.b[0];
	for (int s = 0; s < n; ++s) {
		double x[N], v[N], ev[NE], C[NC], f[N], A[NH], bb[N];
		view.Gather(s, x, v, ev);
		fn(x, ev, v, prm, C, f, A, bb);
		const int* nd = view.Nodes(s);
		for (int a = 0; a < 3; ++a)
			for (int t = 0; t < 3; ++t) { F[3 * nd[a] + t] += f[3 * a + t]; b[3 * nd[a] + t] += bb[3 * a + t]; }
		ScatterA(A, nd, slots + 3 * s, sim.M);
	}
}

// F = Σ(f + d), b = Σ b_e, M = Σ A_e (+ 질량)
static void Assemble(Sim& sim, const Views& vw, const Terms& tm)
{
	const Sheet& c = sim.c;
	sim.M.Zero();
	std::fill(sim.F.begin(), sim.F.end(), 0.);
	std::fill(sim.b.begin(), sim.b.end(), 0.);
	const double pSt[3] = { kStretchK, tm.kdT, tm.h }, pSh[3] = { kShearK, tm.kdS, tm.h }, pP[3] = { tm.kP, tm.kdP, tm.h };
	if (tm.stretch) AddElements(vw.st, &c.triSlot[0], pSt, tm.stretch, sim);
	if (tm.shear) AddElements(vw.sh, &c.triSlot[0], pSh, tm.shear, sim);
	if (tm.poisson) AddElements(vw.po, &c.triSlot[0], pP, tm.poisson, sim);
	if (tm.mass)
		for (int n = 0; n < c.nNode; ++n) {
			double* Db = &sim.M.D[9 * static_cast<std::size_t>(n)];
			Db[0] += c.mass[n]; Db[4] += c.mass[n]; Db[8] += c.mass[n];
		}
}

// 적분 회의 항: 인장 T 판 · 전단 GN · Poisson 은 mode
static Terms RunTerms(int mode, double r, double h)
{
	Terms t;
	t.stretch = &sxgen::BaraffStretch::EvalImplicitTension;
	t.shear = &sxgen::BaraffShear::EvalImplicitGN;
	t.poisson = mode == kModeG ? &sxgen::SulPoisson::EvalImplicitGN : &sxgen::SulPoisson::EvalImplicit;
	t.kP = r * kStretchK; t.kdT = kStretchKd; t.kdS = kShearKd; t.kdP = kPoissonKdRatio * t.kP;
	t.h = h; t.mass = true;
	return t;
}

// ---- 끝 상태 재기 ---------------------------------------------------------------------
struct Measure { double uMin, uMax, vMin, vMax, vMean, cosMax; int flipped; };

static Measure MeasureState(const Sheet& c, const double* X)
{
	Measure m = { 1e30, 0., 1e30, 0., 0., 0., 0 };
	double aSum = 0.;
	for (int t = 0; t < c.nTri; ++t) {
		const double* uv = &c.uv[6 * t];
		const int* nd = &c.tris[3 * t];
		const double du1 = uv[2] - uv[0], dv1 = uv[3] - uv[1], du2 = uv[4] - uv[0], dv2 = uv[5] - uv[1];
		const double det = du1 * dv2 - du2 * dv1;
		double wu[3], wv[3], dx1[3], dx2[3];
		for (int k = 0; k < 3; ++k) {
			dx1[k] = X[3 * nd[1] + k] - X[3 * nd[0] + k];
			dx2[k] = X[3 * nd[2] + k] - X[3 * nd[0] + k];
			wu[k] = (dx1[k] * dv2 - dx2[k] * dv1) / det;
			wv[k] = (-dx1[k] * du2 + dx2[k] * du1) / det;
		}
		const double lu = std::sqrt(wu[0] * wu[0] + wu[1] * wu[1] + wu[2] * wu[2]);
		const double lv = std::sqrt(wv[0] * wv[0] + wv[1] * wv[1] + wv[2] * wv[2]);
		const double cs = std::fabs(wu[0] * wv[0] + wu[1] * wv[1] + wu[2] * wv[2]) / (lu * lv);
		m.uMin = std::min(m.uMin, lu); m.uMax = std::max(m.uMax, lu);
		m.vMin = std::min(m.vMin, lv); m.vMax = std::max(m.vMax, lv);
		m.cosMax = std::max(m.cosMax, cs);
		const double nz = dx1[0] * dx2[1] - dx1[1] * dx2[0];   // 평면 안의 방향: 쉬는 삼각형과 같은 부호여야 한다
		if (!(nz * det > 0.)) ++m.flipped;
		const double area = 0.5 * std::fabs(det);
		m.vMean += area * lv; aSum += area;
	}
	m.vMean /= aSum;
	return m;
}

// ---- 한 회의 적분 --------------------------------------------------------------------
struct RunResult {
	int scene, ir, ib, is, ih, mode;
	bool stable, steady, diverged, reached;
	int steps, nSteps, breakdownSteps, capSteps, precondSteps, maxIters, stepsToSteady, flipped;
	long long iters;
	double endSpeed, endSpeedAbs, drift, errU, errV, cosMax, svMean, sMid;   // endSpeed = 강체 병진을 뺀 것, drift = 뺀 평균 속도의 크기
};

static RunResult Run(Sim& sim, const Views& vw, int scene, int ir, int ib, int is, int ih, int mode, double svRef)
{
	const Sheet& c = sim.c;
	RunResult r = RunResult();
	r.scene = scene; r.ir = ir; r.ib = ib; r.is = is; r.ih = ih; r.mode = mode;
	const double h = 1. / kSps[ih];
	const int nSteps = static_cast<int>(kT * kSps[ih] + 0.5);
	const int tailFrom = nSteps - static_cast<int>(kTail * kSps[ih] + 0.5);
	r.nSteps = nSteps;
	const std::vector<char> fixed = FixedMask(c, scene);
	const double su = 1. + kStrain[is];
	sim.SetBeta(kBeta[ib]);
	sim.SetState(su, 1.);   // 늘려 놓고 놓는다
	const Terms tm = RunTerms(mode, kR[ir], h);
	const int n3 = 3 * c.nNode;
	std::vector<double> dv(n3, 0.);
	PcgWork work;
	int lastAbove = -1;
	for (int step = 0; step < nSteps; ++step) {
		Assemble(sim, vw, tm);
		const PcgOut po = Pcg(c, sim.M, sim.b, dv, work, fixed);
		r.iters += po.iters;
		r.maxIters = std::max(r.maxIters, po.iters);
		if (po.breakdown) ++r.breakdownSteps;
		else if (!po.converged) ++r.capSteps;
		if (po.precond) ++r.precondSteps;
		double vmax = 0.;
		bool finite = true;
		for (int n = 0; n < c.nNode; ++n) {
			double v2 = 0.;
			for (int t = 0; t < 3; ++t) {
				const int d = 3 * n + t;
				if (!fixed[d]) { sim.V[d] += dv[d]; sim.X[d] += h * sim.V[d]; }
				v2 += sim.V[d] * sim.V[d];
				if (!std::isfinite(sim.X[d])) finite = false;
			}
			vmax = std::max(vmax, std::sqrt(v2));
		}
		r.steps = step + 1;
		if (!finite || !(vmax <= kVBlow)) { r.diverged = true; break; }
		// 정상은 강체 병진을 뺀 속력으로 본다(2 차 기준): 장면 U 는 y · z 병진이 풀려 있어 PCG 가 남긴 운동량이 지워지지 않는다.
		// x 는 양 끝이 고정이라 그대로 두고, 장면 C 는 풀린 강체 운동이 없어 빼지 않는다
		double mean[3] = { 0., 0., 0. };
		if (scene == kSceneU) {
			for (int n = 0; n < c.nNode; ++n) { mean[1] += sim.V[3 * n + 1]; mean[2] += sim.V[3 * n + 2]; }
			mean[1] /= c.nNode; mean[2] /= c.nNode;
		}
		double vrel = 0.;
		for (int n = 0; n < c.nNode; ++n) {
			const double ax = sim.V[3 * n], ay = sim.V[3 * n + 1] - mean[1], az = sim.V[3 * n + 2] - mean[2];
			vrel = std::max(vrel, std::sqrt(ax * ax + ay * ay + az * az));
		}
		if (vrel > kSteadySpeed) lastAbove = step;
		if (step >= tailFrom) { r.endSpeed = std::max(r.endSpeed, vrel); r.endSpeedAbs = std::max(r.endSpeedAbs, vmax); }
		r.drift = std::sqrt(mean[1] * mean[1] + mean[2] * mean[2]);
	}
	r.stable = !r.diverged && r.breakdownSteps == 0 && r.capSteps == 0 && r.precondSteps == 0;
	r.steady = !r.diverged && lastAbove < tailFrom;
	r.stepsToSteady = lastAbove + 1;
	if (!r.diverged) {
		const Measure m = MeasureState(c, sim.X);
		r.errU = std::max(std::fabs(m.uMin - su), std::fabs(m.uMax - su));
		r.errV = std::max(std::fabs(m.vMin - svRef), std::fabs(m.vMax - svRef));
		r.cosMax = m.cosMax;
		r.svMean = m.vMean;
		r.flipped = m.flipped;
		const int mid = (kNU - 1) / 2;   // 가운데 단면의 폭 비
		r.sMid = (sim.X[3 * ((kNV - 1) * kNU + mid) + 1] - sim.X[3 * mid + 1]) / (kEll * (kNV - 1));
		r.reached = r.stable && r.steady && r.errU <= kReach && r.errV <= kReach && r.cosMax <= kReach && r.flipped == 0;
	}
	return r;
}

// ---- 전제 ---------------------------------------------------------------------------
// U0a: 균일 변형 (s_u, s_v) 에서 자유 자유도의 힘 최댓값 / 인장만 조립한 절점 힘의 최댓값
static double StaticBalance(Sim& sim, const Views& vw, const std::vector<char>& fixed, double r, double bt, double su, double sv)
{
	const int n3 = 3 * sim.c.nNode;
	sim.SetBeta(bt);
	sim.SetState(su, sv);
	Terms tm = RunTerms(kModeG, r, 1.);
	tm.mass = false;
	Assemble(sim, vw, tm);
	double res = 0.;
	for (int d = 0; d < n3; ++d) if (!fixed[d]) res = std::max(res, std::fabs(sim.F[d]));
	tm.shear = 0; tm.poisson = 0;
	Assemble(sim, vw, tm);
	return Rel(res, MaxAbs(&sim.F[0], n3));
}

// U0b: 세 힘을 정확판으로, kd = 0 · 질량 없이 · h = 1 로 조립하면 A = −∂F/∂x 다. 열을 힘의 중심차분과 견준다
static double CheckAssembly(Sim& sim, const Views& vw, int* nColOut)
{
	const Sheet& c = sim.c;
	const int n3 = 3 * c.nNode;
	sim.SetBeta(0.6);
	sim.SetState(1.2, 0.95);
	std::mt19937_64 rng(7);
	std::uniform_real_distribution<double> u(-1., 1.);
	for (int n = 0; n < c.nNode; ++n) {   // 평면 밖으로도 흩뜨린다
		sim.X[3 * n] += 0.1 * kEll * u(rng);
		sim.X[3 * n + 1] += 0.1 * kEll * u(rng);
		sim.X[3 * n + 2] += 0.05 * kEll * u(rng);
	}
	Terms tm;
	tm.stretch = &sxgen::BaraffStretch::EvalImplicit;
	tm.shear = &sxgen::BaraffShear::EvalImplicit;
	tm.poisson = &sxgen::SulPoisson::EvalImplicit;
	tm.kP = kStretchK; tm.kdT = tm.kdS = tm.kdP = 0.; tm.h = 1.; tm.mass = false;
	Assemble(sim, vw, tm);
	const BlockMat A = sim.M;
	double amax = 0.;
	for (std::size_t i = 0; i < A.D.size(); ++i) amax = std::max(amax, std::fabs(A.D[i]));
	for (std::size_t i = 0; i < A.O.size(); ++i) amax = std::max(amax, std::fabs(A.O[i]));
	std::vector<double> e(n3, 0.), col(n3), Fp(n3), Fm(n3);
	const double dx = 1e-6 * kEll;
	double worst = 0.;
	int nCol = 0;
	for (int j = 0; j < n3; j += 7) {
		std::fill(e.begin(), e.end(), 0.);
		e[j] = 1.;
		MatVec(c, A, &e[0], &col[0]);   // 대칭이라 열 = 행
		const double x0 = sim.X[j];
		sim.X[j] = x0 + dx; Assemble(sim, vw, tm); Fp = sim.F;
		const double xp = sim.X[j];
		sim.X[j] = x0 - dx; Assemble(sim, vw, tm); Fm = sim.F;
		const double hx = xp - sim.X[j];
		sim.X[j] = x0;
		for (int d = 0; d < n3; ++d) worst = std::max(worst, std::fabs(col[d] + (Fp[d] - Fm[d]) / hx));
		++nCol;
	}
	*nColOut = nCol;
	return Rel(worst, amax);
}

// ---- 출력 ---------------------------------------------------------------------------
static int CaseIndex(int ir, int ib, int is) { return (ir * kNB + ib) * kNS + is; }
static int RunIndex(int ir, int ib, int is, int ih, int mode) { return (CaseIndex(ir, ib, is) * kNH + ih) * 2 + mode; }

static char CodeU(const RunResult& r) { return r.reached ? 'O' : r.diverged ? 'D' : !r.stable ? 'F' : !r.steady ? 'n' : 's'; }
static char CodeC(const RunResult& r) { return r.diverged ? 'D' : !r.stable ? 'F' : !r.steady ? 'n' : 'O'; }

static void JsonRun(std::FILE* js, const RunResult& r, char code, const char* tail)
{
	std::fprintf(js, "    {\"scene\": \"%s\", \"r\": %g, \"beta\": %g, \"eps\": %g, \"sps\": %d, \"mode\": \"%s\", \"code\": \"%c\", "
		"\"stable\": %s, \"steady\": %s, \"diverged\": %s, \"reached\": %s, \"steps\": %d, \"breakdown_steps\": %d, \"cap_steps\": %d, "
		"\"precond_steps\": %d, \"iters\": %lld, \"max_iters\": %d, \"steps_to_steady\": %d, \"end_speed\": %.3e, "
		"\"end_speed_abs\": %.3e, \"drift\": %.3e, "
		"\"err_wu\": %.3e, \"err_wv\": %.3e, \"cos_max\": %.3e, \"flipped\": %d, \"sv_mean\": %.12f, \"s_mid\": %.12f}%s\n",
		r.scene == kSceneU ? "U" : "C", kR[r.ir], kBeta[r.ib], kStrain[r.is], kSps[r.ih], kModeName[r.mode], code,
		r.stable ? "true" : "false", r.steady ? "true" : "false", r.diverged ? "true" : "false", r.reached ? "true" : "false",
		r.steps, r.breakdownSteps, r.capSteps, r.precondSteps, r.iters, r.maxIters, r.stepsToSteady, r.endSpeed,
		r.endSpeedAbs, r.drift, r.errU, r.errV, r.cosMax, r.flipped, r.svMean, r.sMid, tail);
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<std::vector<double> > ref = ReadRows("results/mma_ref_poisson_uniaxial.txt");
	if (static_cast<int>(ref.size()) != kNCase) {
		std::printf("results/mma_ref_poisson_uniaxial.txt 가 %d 줄이 아니다 — model/sul_poisson_uniaxial_ref.wl 을 먼저 돌린다\n", static_cast<int>(kNCase));
		return 2;
	}
	std::vector<double> svRef(kNCase), nuLin(kNCase);
	std::vector<char> valid(kNCase, 0);
	for (int ir = 0; ir < kNR; ++ir) for (int ib = 0; ib < kNB; ++ib) for (int is = 0; is < kNS; ++is) {
		const int k = CaseIndex(ir, ib, is);
		const std::vector<double>& row = ref[k];
		if (row.size() != 7 || std::fabs(row[0] - kR[ir]) > 1e-12 || std::fabs(row[1] - kBeta[ib]) > 1e-12 || std::fabs(row[2] - kStrain[is]) > 1e-12) {
			std::printf("기준 %d 번째 줄이 하네스의 경우(r %g, β %g, ε %g)와 다르다\n", k + 1, kR[ir], kBeta[ib], kStrain[is]);
			return 2;
		}
		svRef[k] = row[3];
		nuLin[k] = row[6];
		valid[k] = row[4] > 0. && row[5] == 1.;   // 극소이고 근이 하나
	}

	const Sheet sheet = MakeSheet();
	Sim sim(sheet);
	const Views vw(sim.mesh);
	const std::vector<char> fixedU = FixedMask(sheet, kSceneU);
	double aMin = 1e30, aMax = 0.;
	for (int t = 0; t < sheet.nTri; ++t) { const double a = UvArea(&sheet.uv[6 * t]); aMin = std::min(aMin, a); aMax = std::max(aMax, a); }
	std::printf("Poisson 힘 단축 인장 적분 — 노드 %d, 삼각형 %d (넓이 %.2f – %.2f mm²), 경우 %d, 회 %d × 2 장면\n",
		sheet.nNode, sheet.nTri, 1e6 * aMin, 1e6 * aMax, static_cast<int>(kNCase), static_cast<int>(kNCase * kNH * 2));
	Report rpt("results/gate_poisson_integration.json");
	std::FILE* js = rpt.json;
	std::fprintf(js, "{\n  \"attempt\": 2, \"nodes\": %d, \"triangles\": %d, \"T\": %g, \"k_T\": %g, \"k_shear\": %g, \"rho\": %g,\n",
		sheet.nNode, sheet.nTri, kT, kStretchK, kShearK, kRho);

	// ---- 전제 U0a · U0b -------------------------------------------------------
	double u0a = 0., u0aStart = 1e30;
	int nValid = 0;
	std::fprintf(js, "  \"U0a_static_balance\": [\n");
	for (int ir = 0; ir < kNR; ++ir) for (int ib = 0; ib < kNB; ++ib) for (int is = 0; is < kNS; ++is) {
		const int k = CaseIndex(ir, ib, is);
		const double su = 1. + kStrain[is];
		const double atRef = StaticBalance(sim, vw, fixedU, kR[ir], kBeta[ib], su, svRef[k]);
		const double atStart = StaticBalance(sim, vw, fixedU, kR[ir], kBeta[ib], su, 1.);
		u0a = std::max(u0a, atRef);
		if (kBeta[ib] != 0.) u0aStart = std::min(u0aStart, atStart);
		if (valid[k]) ++nValid;
		std::fprintf(js, "    {\"r\": %g, \"beta\": %g, \"eps\": %g, \"sv_ref\": %.15f, \"valid\": %s, \"at_ref\": %.3e, \"at_start\": %.3e}%s\n",
			kR[ir], kBeta[ib], kStrain[is], svRef[k], valid[k] ? "true" : "false", atRef, atStart, k + 1 < kNCase ? "," : "");
	}
	std::fprintf(js, "  ],\n");
	int nCol = 0;
	const double u0b = CheckAssembly(sim, vw, &nCol);
	std::printf("전제\n");
	rpt.Gate("U0a", u0a <= kU0aLimit, "기준 평형에서 자유 자유도의 힘 / 인장 절점 힘, 24 경우 최대 %.2e", u0a, kU0aLimit);
	std::printf("       (보고: 처음 상태 s_v = 1 에서는 β ≠ 0 인 경우의 최소 %.2e, 기준이 유효한 경우 %d / %d)\n", u0aStart, nValid, static_cast<int>(kNCase));
	rpt.Gate("U0b", u0b <= kU0bLimit, "max|A − A_FD| / max|A| %.2e", u0b, kU0bLimit);
	std::printf("       (열 %d 개)\n", nCol);
	const bool premise = rpt.fail == 0;
	std::fprintf(js, "  \"U0a_max\": %.3e, \"U0a_at_start_min_beta_nonzero\": %.3e, \"U0b\": %.3e, \"U0b_columns\": %d, \"valid_cases\": %d,\n",
		u0a, u0aStart, u0b, nCol, nValid);

	// ---- 적분 -----------------------------------------------------------------
	std::vector<RunResult> runU(kNCase * kNH * 2), runC(kNCase * kNH * 2);
	for (int ir = 0; ir < kNR; ++ir) for (int ib = 0; ib < kNB; ++ib) for (int is = 0; is < kNS; ++is)
		for (int ih = 0; ih < kNH; ++ih) for (int mode = 0; mode < 2; ++mode) {
			const int k = CaseIndex(ir, ib, is);
			runU[RunIndex(ir, ib, is, ih, mode)] = Run(sim, vw, kSceneU, ir, ib, is, ih, mode, svRef[k]);
			runC[RunIndex(ir, ib, is, ih, mode)] = Run(sim, vw, kSceneC, ir, ib, is, ih, mode, svRef[k]);
		}

	// ---- 장면 U: 부호표와 판정 --------------------------------------------------
	std::printf("\n장면 U (양 끝의 x 만 고정) — O 도달 · s 안정 · 정상이지만 기준과 다름 · n 정상이 아님 · F PCG 깨짐 · 상한 · 음의 대각 · D 발산\n");
	std::printf("    r     β     ε   | 1/120 E G | 1/30 E G |  s_v 기준     ν 기준   ν 선형   | G · 1/120 의 s_v 오차\n");
	int violations = 0, reachE = 0, reachG = 0, confirmed = 0, judged = 0;
	long long itE = 0, itG = 0;
	int bothPairs = 0;
	Stat stepsE, stepsG;
	double driftMax = 0., relMax = 0., absMax = 0.;   // 도달한 회의 미끄러지는 속력 · 평균을 뺀 속력 · 빼지 않은 속력 (마지막 0.5 s)
	for (std::size_t i = 0; i < runU.size(); ++i)
		if (runU[i].reached) {
			driftMax = std::max(driftMax, runU[i].drift);
			relMax = std::max(relMax, runU[i].endSpeed);
			absMax = std::max(absMax, runU[i].endSpeedAbs);
		}
	for (int ir = 0; ir < kNR; ++ir) for (int ib = 0; ib < kNB; ++ib) for (int is = 0; is < kNS; ++is) {
		const int k = CaseIndex(ir, ib, is);
		bool any = false;
		char code[kNH][2];
		for (int ih = 0; ih < kNH; ++ih) {
			const RunResult& e = runU[RunIndex(ir, ib, is, ih, kModeE)];
			const RunResult& g = runU[RunIndex(ir, ib, is, ih, kModeG)];
			code[ih][0] = CodeU(e); code[ih][1] = CodeU(g);
			if (!valid[k]) continue;
			if (e.reached && !g.reached) ++violations;
			if (e.reached) { ++reachE; any = true; }
			if (g.reached) { ++reachG; any = true; }
			if (e.reached && g.reached) {
				++bothPairs; itE += e.iters; itG += g.iters;
				stepsE.Add(e.stepsToSteady); stepsG.Add(g.stepsToSteady);
			}
		}
		if (valid[k]) { ++judged; if (any) ++confirmed; }
		const RunResult& g0 = runU[RunIndex(ir, ib, is, 0, kModeG)];
		std::printf("  %4.1f  %+4.1f  %4.2f  |    %c %c    |    %c %c   |  %.9f  %+.4f  %+.4f  | %.1e%s\n",
			kR[ir], kBeta[ib], kStrain[is], code[0][0], code[0][1], code[1][0], code[1][1],
			svRef[k], (1. - svRef[k]) / kStrain[is], nuLin[k], g0.diverged ? NAN : g0.errV, valid[k] ? "" : "  (판정에서 뺌)");
	}
	std::printf("\n판정 (장면 U, 판정 대상 경우 %d · 짝 %d)\n", judged, judged * kNH);
	if (premise) {
		rpt.Gate("(가)", violations == 0, "E 는 도달 · G 는 못 도달한 짝 %.0f", violations, 0.);
		rpt.Gate("(나)", reachG >= reachE, "도달한 회 G %.0f", reachG, reachE);
		std::printf("       (E 도달 %d / %d, G 도달 %d / %d)\n", reachE, judged * kNH, reachG, judged * kNH);
		rpt.Gate("(다)", confirmed == judged && judged > 0, "네 회 가운데 하나 이상 도달한 경우 %.0f", confirmed, judged);
	}
	else std::printf("  전제가 깨져 판정하지 않는다\n");
	const bool gnRule = premise && violations == 0 && reachG >= reachE;
	const char* recommend = !premise ? "판정하지 않음" : gnRule ? "EvalImplicitGN" : "EvalImplicit";
	std::printf("  권고(Poisson 힘의 기본 강성): %s\n", recommend);
	if (bothPairs > 0)
		std::printf("  보고 — 두 판이 모두 도달한 짝 %d 개: PCG 반복 합 E %lld · G %lld (G / E %.3f), 정상에 이른 단계 수 중앙 E %.0f · G %.0f\n",
			bothPairs, itE, itG, static_cast<double>(itG) / static_cast<double>(itE), stepsE.Median(), stepsG.Median());
	std::printf("  보고 — 도달한 회의 마지막 0.5 s: 강체 병진을 뺀 속력 최대 %.1e m/s, 빼지 않은 속력 최대 %.1e m/s, 끝의 미끄러지는 속력 최대 %.1e m/s\n",
		relMax, absMax, driftMax);
	std::fprintf(js, "  \"reached_rel_speed_max\": %.3e, \"reached_abs_speed_max\": %.3e, \"reached_drift_max\": %.3e,\n", relMax, absMax, driftMax);

	// ---- 장면 C: 보고 ------------------------------------------------------------
	std::printf("\n장면 C (양 끝을 모두 고정, 보고만) — O 안정 · 정상 · n 정상이 아님 · F PCG 깨짐 · 상한 · 음의 대각 · D 발산\n");
	std::printf("    r     β     ε   | 1/120 E G | 1/30 E G |  가운데 폭 비 (G · 1/120)   균일 변형의 s_v   수축의 비\n");
	int okE = 0, okG = 0;
	for (int ir = 0; ir < kNR; ++ir) for (int ib = 0; ib < kNB; ++ib) for (int is = 0; is < kNS; ++is) {
		const int k = CaseIndex(ir, ib, is);
		char code[kNH][2];
		for (int ih = 0; ih < kNH; ++ih) {
			const RunResult& e = runC[RunIndex(ir, ib, is, ih, kModeE)];
			const RunResult& g = runC[RunIndex(ir, ib, is, ih, kModeG)];
			code[ih][0] = CodeC(e); code[ih][1] = CodeC(g);
			if (code[ih][0] == 'O') ++okE;
			if (code[ih][1] == 'O') ++okG;
		}
		const RunResult& g0 = runC[RunIndex(ir, ib, is, 0, kModeG)];
		const double mid = g0.diverged ? NAN : g0.sMid;
		const double frac = svRef[k] != 1. ? (1. - mid) / (1. - svRef[k]) : NAN;
		std::printf("  %4.1f  %+4.1f  %4.2f  |    %c %c    |    %c %c   |  %.6f                  %.6f          %.3f\n",
			kR[ir], kBeta[ib], kStrain[is], code[0][0], code[0][1], code[1][0], code[1][1], mid, svRef[k], frac);
	}
	std::printf("  보고 — 안정 · 정상인 회: E %d / %d, G %d / %d\n", okE, static_cast<int>(kNCase * kNH), okG, static_cast<int>(kNCase * kNH));

	// ---- json -------------------------------------------------------------------
	std::fprintf(js, "  \"judged_cases\": %d, \"judged_pairs\": %d, \"violations\": %d, \"reached_E\": %d, \"reached_G\": %d, \"confirmed_cases\": %d,\n",
		judged, judged * kNH, violations, reachE, reachG, confirmed);
	std::fprintf(js, "  \"recommend\": \"%s\", \"both_reached_pairs\": %d, \"iters_E\": %lld, \"iters_G\": %lld,\n", recommend, bothPairs, itE, itG);
	std::fprintf(js, "  \"scene_C_stable_steady_E\": %d, \"scene_C_stable_steady_G\": %d,\n", okE, okG);
	std::fprintf(js, "  \"runs\": [\n");
	for (std::size_t i = 0; i < runU.size(); ++i) JsonRun(js, runU[i], CodeU(runU[i]), ",");
	for (std::size_t i = 0; i < runC.size(); ++i) JsonRun(js, runC[i], CodeC(runC[i]), i + 1 < runC.size() ? "," : "");
	std::fprintf(js, "  ],\n  \"gates_failed\": %d\n}\n", rpt.fail);
	const bool pass = premise && confirmed == judged && judged > 0;
	std::printf("\nPoisson 힘 단축 인장 적분: 전제 %s, (다) 예측 확인 %s, 권고 %s → results/gate_poisson_integration.json\n",
		premise ? "PASS" : "FAIL", !premise ? "—" : (confirmed == judged ? "PASS" : "FAIL"), recommend);
	return pass ? 0 : 1;
}
