//---------------------------------------------------------------------------
// test_integration.cpp — 시간 적분 비교(docs/GATE_INTEGRATION.md).
// 두 모서리에 매단 작은 천을 Baraff–Witkin 식 (6)(단계마다 한 번, 거르개 PCG)으로 적분한다.
// 요소 기여는 생성 코드 EvalImplicit(정확판) · EvalImplicitGN(GN)이, 엔진이 있으면(/DSX_WITH_ENGINE, engine_macri.obj ·
// SDEmysse.obj 와 링크) 엔진 Bridson 굽힘(힘 + Jacobian)도 낸다.
//   인자 없음  I0(전역 조립 대 유한차분) · I1 안정 · I2 차이 · 보고 → results/gate_integration.json
//   timing     기준 짝만 돌려 시간을 낸다 → results/gate_integration_time.json (I3 판정은 tests/timing_bases.py integ)
//   tension    인장 C > 0 판 게이트(docs/GATE_TENSION.md)의 적분 쪽: 인장 E · G · T 180 회, (가′) · (나′) · (라′)
//              → results/gate_tension_integration.json
//   timing_tension  그 게이트 (다′) 의 기준 짝 → results/gate_tension_time.json (판정은 tests/timing_bases.py tension)
//   shear_gn   전단 GN 재비교(docs/GATE_SHEAR_GN.md)의 적분 쪽: 인장 T · 전단 E · G 60 회, (가″) · (나″) · (라″)
//              → results/gate_shear_gn_integration.json
//   timing_shear_gn  그 재비교 (다″) 의 기준 짝 둘 → results/gate_shear_gn_time.json (판정은 tests/timing_bases.py shear_gn)
// 실행 디렉터리 = 리포 루트.
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
#include <cstdlib>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;

// ---- 시나리오 (판정문에 고정) ------------------------------------------------------
static const int    kG = 21;                       // 한 변의 노드 수
static const double kL = 0.5;                      // 한 변 길이 (m)
static const double kEll = kL / (kG - 1);          // 노드 간격 2.5 cm
static const double kRho = 0.2;                    // 면밀도 (kg/m²)
static const double kGrav = 9.8;                   // −z
static const double kStretchK = 1000., kStretchKd = 0.1;
static const double kShearK = 100., kShearKd = 0.01;
static const double kBendK[2] = { 1e-5, 1e-3 };   // 부드러움 · 뻣뻣함 (J), 굽힘 kd = 0
static const double kT = 2.0;
static const int    kNH = 3;
static const int    kStepsPerSec[kNH] = { 120, 60, 30 };   // h = 1/120 · 1/60 · 1/30
static const int    kSamplesPerSec = 30;                    // 표본 시각 k/30 s
static const double kEps = 1e-6;                   // PCG δ ≤ ε² δ₀
static const int    kCap = 2000;                   // PCG 반복 상한
static const double kVBlow = 50.;                  // 발산 속력 (m/s)
static const double kI0Limit = 1e-6;

enum { kE = 0, kGN = 1, kTen = 2 };   // kTen = 인장 C > 0 판(인장에만, docs/GATE_TENSION.md)
enum { kBaraffE = 0, kBaraffG, kBridsonE, kBridsonG, kEngine, kCompactE, kCompactG, kNBend };
static const char* kBendName[kNBend] = { "baraffE", "baraffG", "bridsonE", "bridsonG", "engine", "compactE", "compactG" };
static const char* kModeName[3] = { "E", "G", "T" };
static const double kUsableStrain = 0.25;          // 인장 C > 0 판 게이트의 「쓸 수 있음」: 안정 그리고 최대 변형률 ≤ 0.25

struct Config { int stretch, shear, bend, stiff, ih, sps; };   // sps > 0 이면 h = 1/sps (판정 밖 진단), 아니면 kStepsPerSec[ih]

static double StepH(int ih) { return 1. / kStepsPerSec[ih]; }
static int Sps(const Config& cf) { return cf.sps > 0 ? cf.sps : kStepsPerSec[cf.ih]; }

// ---- 천: 위상 · 쉬는 모양 · 블록 행렬 틀 ---------------------------------------------
struct Cloth {
	int nNode, nTri, nHinge, nPair;
	std::vector<double> X0, mass;
	std::vector<char> pinned;
	std::vector<int> tris, hinges;      // 삼각형(반시계) 3 개 · 힌지 (h, i, j, k) 4 개 — gate::MakeGrid 와 같은 짜임
	std::vector<double> uv, wbar, len0;
	std::vector<int> pairI, pairJ;      // 노드쌍 블록 (i < j)
	std::vector<int> triSlot, hingeSlot;   // 요소의 지역 쌍 (a < b) → 노드쌍 블록 번호
	double massFreeNorm;                // ‖M g‖ / g 의 자유 자유도 몫 = sqrt(Σ_자유 m²)
};

// a < b 인 지역 쌍의 번호: NN = 3 → (0,1)(0,2)(1,2), NN = 4 → (0,1)(0,2)(0,3)(1,2)(1,3)(2,3)
template <int NN>
inline int LocalPair(int a, int b) { return a * NN - a * (a + 1) / 2 + (b - a - 1); }

static Cloth MakeCloth()
{
	Cloth c;
	const int G = kG;
	c.nNode = G * G;
	c.X0.assign(3 * c.nNode, 0.);
	for (int j = 0; j < G; ++j)
		for (int i = 0; i < G; ++i) { const int n = j * G + i; c.X0[3 * n] = kEll * i; c.X0[3 * n + 1] = kEll * j; }
	for (int j = 0; j + 1 < G; ++j)
		for (int i = 0; i + 1 < G; ++i) {
			const int a = j * G + i, b = a + 1, cc = a + G + 1, d = a + G;
			const int t1[3] = { a, b, cc }, t2[3] = { a, cc, d };
			c.tris.insert(c.tris.end(), t1, t1 + 3);
			c.tris.insert(c.tris.end(), t2, t2 + 3);
		}
	c.nTri = static_cast<int>(c.tris.size() / 3);
	std::unordered_map<long long, int> third;
	for (int t = 0; t < c.nTri; ++t)
		for (int e = 0; e < 3; ++e) third[static_cast<long long>(c.tris[3 * t + e]) * c.nNode + c.tris[3 * t + (e + 1) % 3]] = c.tris[3 * t + (e + 2) % 3];
	for (int t = 0; t < c.nTri; ++t)
		for (int e = 0; e < 3; ++e) {
			const int p = c.tris[3 * t + e], q = c.tris[3 * t + (e + 1) % 3];
			if (p > q) continue;
			std::unordered_map<long long, int>::const_iterator it = third.find(static_cast<long long>(q) * c.nNode + p);
			if (it == third.end()) continue;
			const int hh[4] = { c.tris[3 * t + (e + 2) % 3], p, q, it->second };
			c.hinges.insert(c.hinges.end(), hh, hh + 4);
		}
	c.nHinge = static_cast<int>(c.hinges.size() / 4);

	c.mass.assign(c.nNode, 0.);
	c.uv.resize(6 * static_cast<std::size_t>(c.nTri));
	for (int t = 0; t < c.nTri; ++t) {
		for (int a = 0; a < 3; ++a) { const int n = c.tris[3 * t + a]; c.uv[6 * t + 2 * a] = c.X0[3 * n]; c.uv[6 * t + 2 * a + 1] = c.X0[3 * n + 1]; }
		const double area = UvArea(&c.uv[6 * t]);
		for (int a = 0; a < 3; ++a) c.mass[c.tris[3 * t + a]] += kRho * area / 3.;
	}
	c.wbar.resize(c.nHinge);
	c.len0.resize(c.nHinge);
	for (int s = 0; s < c.nHinge; ++s) {   // Bridson w̄ = 3|e|²/(A1 + A2) (쉬는 모양), Length0 = 쉬는 변 길이
		double x[12];
		for (int a = 0; a < 4; ++a) for (int t = 0; t < 3; ++t) x[3 * a + t] = c.X0[3 * c.hinges[4 * s + a] + t];
		const double *xh = x, *xi = x + 3, *xj = x + 6, *xk = x + 9;
		double u[3], w[3], p[3], q[3], e[3];
		for (int t = 0; t < 3; ++t) { u[t] = xi[t] - xh[t]; w[t] = xj[t] - xh[t]; p[t] = xj[t] - xk[t]; q[t] = xi[t] - xk[t]; e[t] = xj[t] - xi[t]; }
		const double nA[3] = { u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0] };
		const double nB[3] = { p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0] };
		const double la = std::sqrt(nA[0] * nA[0] + nA[1] * nA[1] + nA[2] * nA[2]);
		const double lb = std::sqrt(nB[0] * nB[0] + nB[1] * nB[1] + nB[2] * nB[2]);
		const double e2 = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
		c.wbar[s] = 3. * e2 / (0.5 * la + 0.5 * lb);
		c.len0[s] = std::sqrt(e2);
	}
	c.pinned.assign(c.nNode, 0);
	c.pinned[(G - 1) * G] = 1;           // (0, L, 0)
	c.pinned[(G - 1) * G + G - 1] = 1;   // (L, L, 0)
	double m2 = 0.;
	for (int n = 0; n < c.nNode; ++n) if (!c.pinned[n]) m2 += c.mass[n] * c.mass[n];
	c.massFreeNorm = std::sqrt(m2);

	std::unordered_map<long long, int> slot;
	struct Slot {
		static int Get(std::unordered_map<long long, int>& m, Cloth& cl, int p, int q)
		{
			if (p > q) { const int t = p; p = q; q = t; }
			const long long key = static_cast<long long>(p) * cl.nNode + q;
			std::unordered_map<long long, int>::const_iterator it = m.find(key);
			if (it != m.end()) return it->second;
			const int s = static_cast<int>(cl.pairI.size());
			m[key] = s;
			cl.pairI.push_back(p);
			cl.pairJ.push_back(q);
			return s;
		}
	};
	c.triSlot.resize(3 * static_cast<std::size_t>(c.nTri));
	for (int t = 0; t < c.nTri; ++t)
		for (int a = 0; a < 3; ++a) for (int b = a + 1; b < 3; ++b)
			c.triSlot[3 * t + LocalPair<3>(a, b)] = Slot::Get(slot, c, c.tris[3 * t + a], c.tris[3 * t + b]);
	c.hingeSlot.resize(6 * static_cast<std::size_t>(c.nHinge));
	for (int s = 0; s < c.nHinge; ++s)
		for (int a = 0; a < 4; ++a) for (int b = a + 1; b < 4; ++b)
			c.hingeSlot[6 * s + LocalPair<4>(a, b)] = Slot::Get(slot, c, c.hinges[4 * s + a], c.hinges[4 * s + b]);
	c.nPair = static_cast<int>(c.pairI.size());
	return c;
}

// 전체 행렬: 노드 대각 블록 D[9n] + 노드쌍 블록 O[9p](행 = pairI, 열 = pairJ, 아래쪽은 전치), 3 × 3 행 우선
struct BlockMat {
	std::vector<double> D, O;
	void Resize(const Cloth& c) { D.assign(9 * static_cast<std::size_t>(c.nNode), 0.); O.assign(9 * static_cast<std::size_t>(c.nPair), 0.); }
	void Zero() { std::fill(D.begin(), D.end(), 0.); std::fill(O.begin(), O.end(), 0.); }
};

// 요소 A(위 삼각 행 우선, 3·NN 자유도)를 전체 행렬에 더한다
template <int NN>
inline void ScatterA(const double* A, const int* nd, const int* slot, BlockMat& M)
{
	enum { N = 3 * NN };
	int idx = 0;
	for (int p = 0; p < N; ++p) {
		const int a = p / 3, s = p % 3;
		for (int q = p; q < N; ++q, ++idx) {
			const int b = q / 3, t = q % 3;
			const double val = A[idx];
			if (a == b) {
				double* Db = &M.D[9 * static_cast<std::size_t>(nd[a])];
				Db[3 * s + t] += val;
				if (s != t) Db[3 * t + s] += val;
			}
			else {
				double* Ob = &M.O[9 * static_cast<std::size_t>(slot[LocalPair<NN>(a, b)])];
				if (nd[a] < nd[b]) Ob[3 * s + t] += val; else Ob[3 * t + s] += val;
			}
		}
	}
}

static void MatVec(const Cloth& c, const BlockMat& M, const double* x, double* y)
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

static PcgOut Pcg(const Cloth& cl, const BlockMat& M, const std::vector<double>& b, std::vector<double>& x, PcgWork& w)
{
	const int n3 = 3 * cl.nNode;
	PcgOut o = { 0, false, false, false };
	w.r.resize(n3); w.c.resize(n3); w.q.resize(n3); w.s.resize(n3); w.inv.resize(n3);
	for (int n = 0; n < cl.nNode; ++n)
		for (int t = 0; t < 3; ++t) {
			const int d = 3 * n + t;
			if (cl.pinned[n]) { w.inv[d] = 0.; continue; }
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
	const double stop = kEps * kEps * d0;
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
struct Terms {
	int stretch, shear, bend;   // 판(kE · kGN, 인장은 kTen 도), 없으면 −1. bend 는 kBaraffE …, 없으면 −1
	double kb, h;
	bool mass;                  // 질량 · 중력
};

struct Sim {
	const Cloth& c;
	SXMesh mesh;
	double* X;
	double* V;
	BlockMat M;
	std::vector<double> F, b;
	double maxStrainC;          // 마지막 조립의 max|C_인장| (√a 로 나누기 전)
	double tEngine;             // 엔진 함수 시간 (s, 누적)
	double tExtract;            // 하네스가 엔진 출력으로 12 × 12 를 만들고 흩뿌린 시간 (s, 누적)
	int skipped;
#ifdef SX_WITH_ENGINE
	engine_ref::Bridson* eb;
#endif
	explicit Sim(const Cloth& cl) : c(cl), X(0), V(0), maxStrainC(0.), tEngine(0.), tExtract(0.), skipped(0)
	{
		// 명세가 선언한 변수를 모두 먼저 만들고(같은 이름은 다시 잡힌다) 채운다
		sx::CreateVariables<sxgen::BaraffStretch>(mesh, c.nNode, c.nTri);
		sx::CreateVariables<sxgen::BaraffShear>(mesh, c.nNode, c.nTri);
		sx::CreateVariables<sxgen::BaraffBend>(mesh, c.nNode, c.nHinge);
		sx::CreateVariables<sxgen::BaraffBendCompact>(mesh, c.nNode, c.nHinge);
		sx::CreateVariables<sxgen::BridsonBend>(mesh, c.nNode, c.nHinge);
		SXNodeNamedContainer& px = mesh.GetNodeV(sxgen::BaraffStretch::PositionVar());
		SXElementNamedContainer& pt = mesh.GetElementV(sxgen::BaraffStretch::StencilVar());
		SXNodeNamedContainer& puv = mesh.GetNodeV(sxgen::BaraffStretch::ElemVar(0));
		SXElementNamedContainer& ph = mesh.GetElementV(sxgen::BaraffBend::StencilVar());
		SXNodeNamedContainer& pth = mesh.GetNodeV(sxgen::BaraffBend::ElemVar(0));
		SXNodeNamedContainer& pw = mesh.GetNodeV(sxgen::BridsonBend::ElemVar(0));
		for (int n = 0; n < c.nNode; ++n) px.Add(n, &c.X0[3 * n], 3);
		for (int t = 0; t < c.nTri; ++t) { pt.Add(t, &c.tris[3 * t], 3); for (int k = 0; k < 6; ++k) puv[t][k] = c.uv[6 * t + k]; }
		for (int s = 0; s < c.nHinge; ++s) { ph.Add(s, &c.hinges[4 * s], 4); pth[s][0] = 0.; pw[s][0] = c.wbar[s]; }
		X = px.Data();
		V = mesh.GetNodeV(sxgen::BaraffStretch::VelocityVar()).Data();
		M.Resize(c);
		F.assign(3 * c.nNode, 0.);
		b.assign(3 * c.nNode, 0.);
#ifdef SX_WITH_ENGINE
		eb = engine_ref::CreateBridson(c.nHinge, &c.hinges[0], c.nNode, X, V, &c.len0[0]);
		engine_ref::SetSerialBridson(eb, false);   // 힌지 소유 슬롯까지 (B8 · N2 와 같다)
#endif
	}
	~Sim()
	{
#ifdef SX_WITH_ENGINE
		engine_ref::DestroyBridson(eb);
#endif
	}
	void Reset()
	{
		std::copy(c.X0.begin(), c.X0.end(), X);
		std::fill(V, V + 3 * c.nNode, 0.);
	}
};

typedef std::chrono::steady_clock Clock;
static double Sec(const Clock::time_point& a, const Clock::time_point& b) { return std::chrono::duration<double>(b - a).count(); }

// 인장 C > 0 판은 인장 명세에만 있다(명세의 IMPLICIT_TENSION). 다른 명세로 부르면 하네스 잘못이다
template <class Spec> struct TensionEval {
	static void Run(const double*, const double*, const double*, const double*, double*, double*, double*, double*) { std::abort(); }
};
template <> struct TensionEval<sxgen::BaraffStretch> {
	static void Run(const double* x, const double* ev, const double* v, const double* prm, double* C, double* f, double* A, double* b)
	{ sxgen::BaraffStretch::EvalImplicitTension(x, ev, v, prm, C, f, A, b); }
};

template <class Spec>
static void AddElements(const sx::StencilView<Spec>& view, const int* slots, const double* prm, int mode, Sim& sim, bool strain)
{
	enum { NN = Spec::kNodes, N = Spec::kDof, NE = Spec::kElemValues, NC = Spec::kComp, NH = Spec::kHess, NP = NN * (NN - 1) / 2 };
	const int n = view.size();
	double* F = &sim.F[0];
	double* b = &sim.b[0];
	double cmax = 0.;
	for (int s = 0; s < n; ++s) {
		double x[N], v[N], ev[NE], C[NC], f[N], A[NH], bb[N];
		view.Gather(s, x, v, ev);
		if (mode == kGN) Spec::EvalImplicitGN(x, ev, v, prm, C, f, A, bb);
		else if (mode == kTen) TensionEval<Spec>::Run(x, ev, v, prm, C, f, A, bb);
		else Spec::EvalImplicit(x, ev, v, prm, C, f, A, bb);
		const int* nd = view.Nodes(s);
		for (int a = 0; a < NN; ++a)
			for (int t = 0; t < 3; ++t) { F[3 * nd[a] + t] += f[3 * a + t]; b[3 * nd[a] + t] += bb[3 * a + t]; }
		ScatterA<NN>(A, nd, slots + NP * s, sim.M);
		if (strain) for (int k = 0; k < NC; ++k) cmax = std::max(cmax, std::fabs(C[k]));
	}
	if (strain) sim.maxStrainC = cmax;
}

#ifdef SX_WITH_ENGINE
// 엔진 Bridson: 힌지별 f · J(∂f/∂x, 12 × 12) → A = −h² J, b = h f + h² J v. 엔진 굽힘 감쇠는 꺼져 있다
static void AddEngine(Sim& sim, double h)
{
	const Cloth& c = sim.c;
	const Clock::time_point t0 = Clock::now();
	engine_ref::RunBridson(sim.eb);
	const Clock::time_point t1 = Clock::now();
	double* F = &sim.F[0];
	double* b = &sim.b[0];
	for (int s = 0; s < c.nHinge; ++s) {
		if (engine_ref::SkippedBridson(sim.eb, s)) { ++sim.skipped; continue; }
		double f[12], J[144], A[78], bb[12], v[12];
		engine_ref::ForceBridson(sim.eb, s, f);
		engine_ref::JacobianBridson(sim.eb, s, J);
		const int* nd = &c.hinges[4 * s];
		for (int a = 0; a < 4; ++a) for (int t = 0; t < 3; ++t) v[3 * a + t] = sim.V[3 * nd[a] + t];
		int idx = 0;
		for (int p = 0; p < 12; ++p) {
			double Jv = 0.;
			for (int q = 0; q < 12; ++q) Jv += J[12 * p + q] * v[q];
			bb[p] = h * f[p] + h * h * Jv;
			for (int q = p; q < 12; ++q) A[idx++] = -h * h * J[12 * p + q];
		}
		for (int a = 0; a < 4; ++a)
			for (int t = 0; t < 3; ++t) { F[3 * nd[a] + t] += f[3 * a + t]; b[3 * nd[a] + t] += bb[3 * a + t]; }
		ScatterA<4>(A, nd, &c.hingeSlot[6 * static_cast<std::size_t>(s)], sim.M);
	}
	const Clock::time_point t2 = Clock::now();
	sim.tEngine += Sec(t0, t1);
	sim.tExtract += Sec(t1, t2);
}
#endif

struct Views {
	sx::StencilView<sxgen::BaraffStretch> st;
	sx::StencilView<sxgen::BaraffShear> sh;
	sx::StencilView<sxgen::BaraffBend> ba;
	sx::StencilView<sxgen::BaraffBendCompact> co;
	sx::StencilView<sxgen::BridsonBend> br;
	explicit Views(const SXMesh& m) : st(m), sh(m), ba(m), co(m), br(m) {}
};

// F = Σ(f + d) (+ M g), b = Σ b_e (+ h M g), M = Σ A_e (+ M)
static void Assemble(Sim& sim, const Views& vw, const Terms& tm)
{
	const Cloth& c = sim.c;
	sim.M.Zero();
	std::fill(sim.F.begin(), sim.F.end(), 0.);
	std::fill(sim.b.begin(), sim.b.end(), 0.);
	const double pSt[3] = { kStretchK, kStretchKd, tm.h }, pSh[3] = { kShearK, kShearKd, tm.h }, pB[3] = { tm.kb, 0., tm.h };
	if (tm.stretch >= 0) AddElements(vw.st, &c.triSlot[0], pSt, tm.stretch, sim, true);
	if (tm.shear >= 0) AddElements(vw.sh, &c.triSlot[0], pSh, tm.shear == kGN ? kGN : kE, sim, false);
	switch (tm.bend) {
	case kBaraffE: case kBaraffG: AddElements(vw.ba, &c.hingeSlot[0], pB, tm.bend == kBaraffG ? kGN : kE, sim, false); break;
	case kCompactE: case kCompactG: AddElements(vw.co, &c.hingeSlot[0], pB, tm.bend == kCompactG ? kGN : kE, sim, false); break;
	case kBridsonE: case kBridsonG: AddElements(vw.br, &c.hingeSlot[0], pB, tm.bend == kBridsonG ? kGN : kE, sim, false); break;
#ifdef SX_WITH_ENGINE
	case kEngine: engine_ref::SetStiffnessBridson(sim.eb, tm.kb); AddEngine(sim, tm.h); break;
#endif
	default: break;
	}
	if (tm.mass)
		for (int n = 0; n < c.nNode; ++n) {
			double* Db = &sim.M.D[9 * static_cast<std::size_t>(n)];
			Db[0] += c.mass[n]; Db[4] += c.mass[n]; Db[8] += c.mass[n];
			sim.F[3 * n + 2] -= c.mass[n] * kGrav;
			sim.b[3 * n + 2] -= tm.h * c.mass[n] * kGrav;
		}
}

// ---- 한 회의 적분 --------------------------------------------------------------------
struct RunResult {
	Config cfg;
	bool ran, stable, diverged;
	int steps, nSteps, breakdownSteps, capSteps, precondSteps, maxIters, skipped;
	long long iters;
	double tAsm, tPcg, tEngine, tExtract;   // s
	double maxSpeed, finalKE, maxStrain;
	Stat beRes;
	int nSamples;
	std::vector<double> samples;            // 표본 시각마다 X (3 × 노드)
	RunResult() : ran(false), stable(false), diverged(false), steps(0), nSteps(0), breakdownSteps(0), capSteps(0), precondSteps(0),
		maxIters(0), skipped(0), iters(0), tAsm(0.), tPcg(0.), tEngine(0.), tExtract(0.), maxSpeed(0.), finalKE(0.), maxStrain(0.), nSamples(0) {}
};

static Terms TermsOf(const Config& cf)
{
	Terms t;
	t.stretch = cf.stretch; t.shear = cf.shear; t.bend = cf.bend;
	t.kb = kBendK[cf.stiff]; t.h = 1. / Sps(cf); t.mass = true;
	return t;
}

// keepAt > 0 이면 그 단계 수를 마친 상태(또는 그 전 마지막 유한 상태)를 keep 에 남긴다
static RunResult Run(Sim& sim, const Views& vw, const Config& cf, bool keepSamples, int keepAt, std::vector<double>* keep)
{
	const Cloth& c = sim.c;
	RunResult r;
	r.cfg = cf;
	r.ran = true;
	const Terms tm = TermsOf(cf);
	const double h = tm.h;
	const int nSteps = static_cast<int>(kT * Sps(cf) + 0.5);
	const int perSample = std::max(1, Sps(cf) / kSamplesPerSec);
	r.nSteps = nSteps;
	sim.Reset();
	sim.tEngine = sim.tExtract = 0.;
	sim.skipped = 0;
	const int n3 = 3 * c.nNode;
	std::vector<double> dv(n3, 0.), prevDv(n3, 0.), prevX(n3);
	PcgWork work;
	const double gScale = h * kGrav * c.massFreeNorm;   // ‖h M g‖ (자유 자유도)
	double strainMax = 0.;
	const double sqrtA = kEll / std::sqrt(2.);          // 삼각형마다 같은 uv 넓이 ℓ²/2
	for (int step = 0; step <= nSteps; ++step) {
		const Clock::time_point t0 = Clock::now();
		Assemble(sim, vw, tm);
		const Clock::time_point t1 = Clock::now();
		if (step > 0) {   // 앞 단계의 backward Euler 잔차: M Δv − h F(x₊, v₊)
			double r2 = 0.;
			for (int n = 0; n < c.nNode; ++n) {
				if (c.pinned[n]) continue;
				for (int t = 0; t < 3; ++t) { const double e = c.mass[n] * prevDv[3 * n + t] - h * sim.F[3 * n + t]; r2 += e * e; }
			}
			r.beRes.Add(std::sqrt(r2) / gScale);
			strainMax = std::max(strainMax, sim.maxStrainC / sqrtA);
		}
		if (step == nSteps) break;   // 마지막 조립은 잔차만 — 시간에 넣지 않는다
		r.tAsm += Sec(t0, t1);
		if (keep && step == keepAt) keep->assign(sim.X, sim.X + n3);
		const Clock::time_point t2 = Clock::now();
		const PcgOut po = Pcg(c, sim.M, sim.b, dv, work);
		const Clock::time_point t3 = Clock::now();
		r.tPcg += Sec(t2, t3);
		r.iters += po.iters;
		r.maxIters = std::max(r.maxIters, po.iters);
		if (po.breakdown) ++r.breakdownSteps;
		else if (!po.converged) ++r.capSteps;
		if (po.precond) ++r.precondSteps;
		if (keep) std::copy(sim.X, sim.X + n3, prevX.begin());
		double vmax = 0.;
		bool finite = true;
		for (int n = 0; n < c.nNode; ++n) {
			if (c.pinned[n]) { sim.V[3 * n] = sim.V[3 * n + 1] = sim.V[3 * n + 2] = 0.; continue; }
			double v2 = 0.;
			for (int t = 0; t < 3; ++t) {
				const int d = 3 * n + t;
				sim.V[d] += dv[d];
				sim.X[d] += h * sim.V[d];
				v2 += sim.V[d] * sim.V[d];
				if (!std::isfinite(sim.X[d])) finite = false;
			}
			vmax = std::max(vmax, std::sqrt(v2));
		}
		std::copy(dv.begin(), dv.end(), prevDv.begin());
		r.steps = step + 1;
		if (!finite || !(vmax <= kVBlow)) {
			r.diverged = true;
			if (keep && step < keepAt) *keep = prevX;
			break;
		}
		r.maxSpeed = std::max(r.maxSpeed, vmax);
		if (keepSamples && (step + 1) % perSample == 0) { r.samples.insert(r.samples.end(), sim.X, sim.X + n3); ++r.nSamples; }
	}
	if (keep && keep->empty()) keep->assign(sim.X, sim.X + n3);
	double ke = 0.;
	for (int n = 0; n < c.nNode; ++n) for (int t = 0; t < 3; ++t) ke += 0.5 * c.mass[n] * sim.V[3 * n + t] * sim.V[3 * n + t];
	r.finalKE = r.diverged ? NAN : ke;
	r.maxStrain = strainMax;
	r.tEngine = sim.tEngine;
	r.tExtract = sim.tExtract;
	r.skipped = sim.skipped;
	r.stable = !r.diverged && r.breakdownSteps == 0 && r.capSteps == 0 && r.precondSteps == 0;
	return r;
}

// D(P, Q) = max_t RMS_노드 |x_P − x_Q| (공통 표본 시각)
static double Dist(const RunResult& p, const RunResult& q, int nNode)
{
	const int ns = std::min(p.nSamples, q.nSamples);
	double worst = 0.;
	for (int k = 0; k < ns; ++k) {
		const double* a = &p.samples[static_cast<std::size_t>(k) * 3 * nNode];
		const double* b = &q.samples[static_cast<std::size_t>(k) * 3 * nNode];
		double s = 0.;
		for (int d = 0; d < 3 * nNode; ++d) s += (a[d] - b[d]) * (a[d] - b[d]);
		worst = std::max(worst, std::sqrt(s / nNode));
	}
	return worst;
}

// ---- I0: 전역 조립 대 유한차분 ------------------------------------------------------------
static double CheckI0(Sim& sim, const Views& vw, const std::vector<double>& x05, Terms tm, const char* name)
{
	const Cloth& c = sim.c;
	const int n3 = 3 * c.nNode;
	std::copy(x05.begin(), x05.end(), sim.X);
	std::fill(sim.V, sim.V + n3, 0.);
	tm.mass = false;
	Assemble(sim, vw, tm);
	const BlockMat A = sim.M;
	double amax = 0.;
	for (std::size_t i = 0; i < A.D.size(); ++i) amax = std::max(amax, std::fabs(A.D[i]));
	for (std::size_t i = 0; i < A.O.size(); ++i) amax = std::max(amax, std::fabs(A.O[i]));
	std::vector<double> e(n3, 0.), col(n3), Fp(n3), Fm(n3);
	const double h = tm.h, dx = 1e-6 * kEll, dvs = 1e-3;
	double worst = 0.;
	int nCol = 0, free = 0;
	for (int j = 0; j < n3; ++j) {
		if (c.pinned[j / 3]) continue;
		if (free++ % 7 != 0) continue;
		std::fill(e.begin(), e.end(), 0.);
		e[j] = 1.;
		MatVec(c, A, &e[0], &col[0]);   // 대칭이라 열 = 행
		const double x0 = sim.X[j];
		sim.X[j] = x0 + dx; Assemble(sim, vw, tm); Fp = sim.F;
		const double xp = sim.X[j];
		sim.X[j] = x0 - dx; Assemble(sim, vw, tm); Fm = sim.F;
		const double hx = xp - sim.X[j];
		sim.X[j] = x0;
		std::vector<double> fdx(n3);
		for (int d = 0; d < n3; ++d) fdx[d] = (Fp[d] - Fm[d]) / hx;
		sim.V[j] = dvs; Assemble(sim, vw, tm); Fp = sim.F;
		const double vp = sim.V[j];
		sim.V[j] = -dvs; Assemble(sim, vw, tm); Fm = sim.F;
		const double hv = vp - sim.V[j];
		sim.V[j] = 0.;
		for (int d = 0; d < n3; ++d) {
			const double fd = -h * (Fp[d] - Fm[d]) / hv - h * h * fdx[d];
			worst = std::max(worst, std::fabs(col[d] - fd));
		}
		++nCol;
	}
	const double rel = Rel(worst, amax);
	std::printf("  I0 %-8s 열 %d 개, max|A − A_FD| / max|A| = %.2e\n", name, nCol, rel);
	return rel;
}

// ---- 출력 ---------------------------------------------------------------------------
static std::string Key(const Config& cf)
{
	char s[64];
	std::snprintf(s, sizeof s, "%s%s_%s", kModeName[cf.stretch], kModeName[cf.shear], kBendName[cf.bend]);
	return s;
}

static void PrintRun(const RunResult& r)
{
	const double perS = 1e3 / kT;   // s 합 → ms / 모의 1 초
	std::printf("  %-14s k%-5g h1/%-3d %s%s 단계 %3d/%3d  PCG 평균 %5.1f 최대 %4d 깨짐 %3d 상한 %3d  조립 %6.1f · 풀이 %7.1f ms/s  "
		"BE 잔차 중앙 %.1e 최대 %.1e  최대 속력 %.2f  변형률 %.3f\n",
		Key(r.cfg).c_str(), kBendK[r.cfg.stiff], Sps(r.cfg), r.stable ? "안정" : (r.diverged ? "발산" : "실패"),
		r.precondSteps ? "*" : " ", r.steps, r.nSteps, r.steps ? static_cast<double>(r.iters) / r.steps : 0., r.maxIters,
		r.breakdownSteps, r.capSteps, r.tAsm * perS, r.tPcg * perS, r.beRes.Median(), r.beRes.Max(), r.maxSpeed, r.maxStrain);
}

static void JsonRun(std::FILE* js, const RunResult& r, const char* tail)
{
	const double perS = 1e3 / kT;
	std::fprintf(js, "    {\"stretch\": \"%s\", \"shear\": \"%s\", \"bend\": \"%s\", \"kb\": %g, \"steps_per_s\": %d, \"stable\": %s, \"diverged\": %s, "
		"\"steps\": %d, \"n_steps\": %d, \"pcg_iters\": %lld, \"pcg_mean\": %.2f, \"pcg_max\": %d, \"breakdown_steps\": %d, \"cap_steps\": %d, "
		"\"precond_steps\": %d, \"asm_ms_per_s\": %.2f, \"pcg_ms_per_s\": %.2f, \"engine_ms_per_s\": %.2f, \"extract_ms_per_s\": %.2f, "
		"\"be_res_median\": %.3e, \"be_res_max\": %.3e, \"max_speed\": %.4f, \"final_ke\": %.4e, \"max_strain\": %.4f, \"engine_skipped\": %d}%s\n",
		kModeName[r.cfg.stretch], kModeName[r.cfg.shear], kBendName[r.cfg.bend], kBendK[r.cfg.stiff], Sps(r.cfg),
		r.stable ? "true" : "false", r.diverged ? "true" : "false", r.steps, r.nSteps, r.iters,
		r.steps ? static_cast<double>(r.iters) / r.steps : 0., r.maxIters, r.breakdownSteps, r.capSteps, r.precondSteps,
		r.tAsm * perS, r.tPcg * perS, r.tEngine * perS, r.tExtract * perS, r.beRes.Median(), r.beRes.Max(), r.maxSpeed,
		std::isfinite(r.finalKE) ? r.finalKE : -1., r.maxStrain, r.skipped, tail);
}

// ---- 권고 규칙 (가) · (나) ----------------------------------------------------------------
static int Index(int stretch, int shear, int bend, int stiff, int ih) { return (((stretch * 2 + shear) * 5 + bend) * 2 + stiff) * kNH + ih; }

struct Rule {
	const char* name;
	int pairs, violations, diffPairs, diffFails;
	double worstRatio;
	bool refE, refG;
	std::vector<std::string> notes;
	Rule() : name(""), pairs(0), violations(0), diffPairs(0), diffFails(0), worstRatio(0.), refE(false), refG(false) {}
	bool A() const { return violations == 0; }
	bool B() const { return diffPairs > 0 && diffFails == 0; }
};

static std::string Desc(const Config& cf)
{
	char s[96];
	std::snprintf(s, sizeof s, "%s k=%g h=1/%d", Key(cf).c_str(), kBendK[cf.stiff], Sps(cf));
	return s;
}

// 힘 which(0 인장, 1 전단, 2 Baraff 굽힘, 3 Bridson 굽힘)의 판을 g(kE · kGN)로 둔 회의 번호.
// 나머지: which < 2 이면 o1 = 다른 삼각형 힘의 판, o2 = 굽힘 변형. which ≥ 2 이면 o1 = 인장 판, o2 = 전단 판
static int PairIndex(int which, int g, int o1, int o2, int stiff, int ih)
{
	if (which == 0) return Index(g, o1, o2, stiff, ih);
	if (which == 1) return Index(o1, g, o2, stiff, ih);
	if (which == 2) return Index(o1, o2, g == kE ? kBaraffE : kBaraffG, stiff, ih);
	return Index(o1, o2, g == kE ? kBridsonE : kBridsonG, stiff, ih);
}

static Rule MakeRule(int which, const std::vector<RunResult>& runs, int nNode, bool engine)
{
	static const char* names[4] = { "stretch", "shear", "baraff", "bridson" };
	Rule r;
	r.name = names[which];
	const int n2 = which < 2 ? (engine ? 5 : 4) : 2;
	for (int o1 = 0; o1 < 2; ++o1)
		for (int o2 = 0; o2 < n2; ++o2)
			for (int stiff = 0; stiff < 2; ++stiff)
				for (int ih = 0; ih < kNH; ++ih) {
					const int iE = PairIndex(which, kE, o1, o2, stiff, ih), iG = PairIndex(which, kGN, o1, o2, stiff, ih);
					const RunResult& E = runs[iE];
					const RunResult& G = runs[iG];
					if (!E.ran || !G.ran) continue;
					++r.pairs;
					if (E.stable && !G.stable) { ++r.violations; r.notes.push_back("(가) E 안정 · G 불안정: " + Desc(G.cfg)); }
					const bool isRef = stiff == 0 && ih == 1 && (which < 2 ? (o1 == kE && o2 == kBridsonG) : (o1 == kE && o2 == kE));
					if (isRef) { r.refE = E.stable; r.refG = G.stable; }
					if (ih == 0) continue;
					const RunResult& H = runs[PairIndex(which, kE, o1, o2, stiff, ih - 1)];
					if (!(E.stable && G.stable && H.stable)) continue;
					const double dGE = Dist(G, E, nNode), dEH = Dist(E, H, nNode);
					++r.diffPairs;
					r.worstRatio = std::max(r.worstRatio, dGE / dEH);
					if (dGE > dEH) {
						++r.diffFails;
						char s[160];
						std::snprintf(s, sizeof s, "(나) D(G,E) %.2e > D(E,E/2) %.2e: %s", dGE, dEH, Desc(G.cfg).c_str());
						r.notes.push_back(s);
					}
				}
	return r;
}

static void JsonRule(std::FILE* js, const Rule& r, const char* tail)
{
	std::fprintf(js, "    \"%s\": {\"pairs\": %d, \"stable_violations\": %d, \"A_ok\": %s, \"diff_pairs\": %d, \"diff_fails\": %d, "
		"\"diff_worst_ratio\": %.4f, \"B_ok\": %s, \"ref_E_stable\": %s, \"ref_G_stable\": %s, \"notes\": [",
		r.name, r.pairs, r.violations, r.A() ? "true" : "false", r.diffPairs, r.diffFails, r.worstRatio, r.B() ? "true" : "false",
		r.refE ? "true" : "false", r.refG ? "true" : "false");
	for (std::size_t i = 0; i < r.notes.size(); ++i) std::fprintf(js, "%s\"%s\"", i ? ", " : "", r.notes[i].c_str());
	std::fprintf(js, "]}%s\n", tail);
}

// ---- 시간 모드: 기준 짝 (+ 보고: 뻣뻣한 굽힘의 같은 짝) ------------------------------------------
// 기준 짝은 판정문대로 부드러운 굽힘 · h = 1/60. 첫 실행에서 그 짝이 정확판 · GN 모두 불안정이라 (다) 를 잴 수 없어,
// 같은 짝을 뻣뻣한 굽힘(모두 안정)으로도 잰다 — 판정 밖 보고(판정문 「측정 방법 변경」). 열쇠 끝에 _stiff.
static int TimingMode(Sim& sim, const Views& vw)
{
	const Config refs[] = {
		{ kE, kE, kBridsonG, 0, 1 }, { kE, kE, kBridsonE, 0, 1 }, { kE, kE, kBaraffE, 0, 1 }, { kE, kE, kBaraffG, 0, 1 },
		{ kE, kE, kCompactE, 0, 1 }, { kE, kE, kCompactG, 0, 1 }, { kGN, kE, kBridsonG, 0, 1 }, { kE, kGN, kBridsonG, 0, 1 },
#ifdef SX_WITH_ENGINE
		{ kE, kE, kEngine, 0, 1 },
#endif
		{ kE, kE, kBridsonG, 1, 1 }, { kE, kE, kBridsonE, 1, 1 }, { kE, kE, kBaraffE, 1, 1 }, { kE, kE, kBaraffG, 1, 1 },
		{ kE, kE, kCompactE, 1, 1 }, { kE, kE, kCompactG, 1, 1 }, { kGN, kE, kBridsonG, 1, 1 }, { kE, kGN, kBridsonG, 1, 1 },
#ifdef SX_WITH_ENGINE
		{ kE, kE, kEngine, 1, 1 },
#endif
	};
	const int n = sizeof refs / sizeof refs[0];
	std::vector<RunResult> rs;
	for (int i = 0; i < n; ++i) { rs.push_back(Run(sim, vw, refs[i], false, -1, 0)); PrintRun(rs.back()); }
	std::FILE* js = std::fopen("results/gate_integration_time.json", "w");
	if (!js) return 2;
	const double perS = 1e3 / kT;
	const char* blocks[] = { "ms_per_sim_second", "asm_ms_per_sim_second", "pcg_iters", "stable" };
	std::fprintf(js, "{\n");
	for (int k = 0; k < 4; ++k) {
		std::fprintf(js, "  \"%s\": {", blocks[k]);
		for (int i = 0; i < n; ++i) {
			const RunResult& r = rs[i];
			const double v = k == 0 ? (r.tAsm + r.tPcg) * perS : k == 1 ? r.tAsm * perS : k == 2 ? static_cast<double>(r.iters) : (r.stable ? 1. : 0.);
			std::fprintf(js, "%s\"%s%s\": %.3f", i ? ", " : "", Key(r.cfg).c_str(), r.cfg.stiff ? "_stiff" : "", v);
		}
		std::fprintf(js, "}%s\n", k < 3 ? "," : "");
	}
	std::fprintf(js, "}\n");
	std::fclose(js);
	return 0;
}

// ---- 인장 C > 0 판 게이트 (docs/GATE_TENSION.md) ---------------------------------------------
static bool Usable(const RunResult& r) { return r.ran && r.stable && r.maxStrain <= kUsableStrain; }
static char Code(const RunResult& r) { return Usable(r) ? 'O' : r.stable ? 'X' : r.diverged ? 'D' : 'F'; }   // X = 안정이지만 변형률 > 0.25

static int TensionMode(Sim& sim, const Views& vw, int nNode, bool engine)
{
	std::vector<RunResult> runs(3 * 2 * 5 * 2 * kNH);
	const Clock::time_point tStart = Clock::now();
	for (int bend = 0; bend < 5; ++bend) {
		if (bend == kEngine && !engine) continue;
		for (int st = 0; st < 3; ++st)
			for (int sh = 0; sh < 2; ++sh)
				for (int stiff = 0; stiff < 2; ++stiff)
					for (int ih = 0; ih < kNH; ++ih) {
						const Config cf = { st, sh, bend, stiff, ih, 0 };
						runs[Index(st, sh, bend, stiff, ih)] = Run(sim, vw, cf, true, -1, 0);
						PrintRun(runs[Index(st, sh, bend, stiff, ih)]);
					}
	}
	std::printf("  (적분 %.1f s)\n", Sec(tStart, Clock::now()));

	// 표: 행 = 인장 · 전단 판, 칸 = 굽힘 k · 굽힘 변형, 글자 = h 1/120 · 1/60 · 1/30 (O 쓸 수 있음 · X 안정이지만 변형률 > 0.25 · D 발산 · F 풀이 실패)
	std::printf("쓸 수 있음 표 (O · X 변형률 > %.2f · D 발산 · F 풀이 실패), 칸 = h 1/120 · 1/60 · 1/30\n", kUsableStrain);
	for (int st = 0; st < 3; ++st)
		for (int sh = 0; sh < 2; ++sh) {
			std::printf("  %s%s ", kModeName[st], kModeName[sh]);
			for (int stiff = 0; stiff < 2; ++stiff) {
				std::printf(" | k%-5g", kBendK[stiff]);
				for (int bend = 0; bend < 5; ++bend) {
					if (bend == kEngine && !engine) continue;
					std::printf(" %s ", kBendName[bend]);
					for (int ih = 0; ih < kNH; ++ih) std::printf("%c", Code(runs[Index(st, sh, bend, stiff, ih)]));
				}
			}
			std::printf("\n");
		}

	// (가′) · (나′) · (라′) 와 보고
	int pairs = 0, viol = 0, diffPairs = 0, diffFails = 0, nE = 0, nT = 0, nG = 0, gNotT = 0;
	double worst = 0.;
	std::vector<std::string> notes;
	for (int sh = 0; sh < 2; ++sh)
		for (int bend = 0; bend < 5; ++bend) {
			if (bend == kEngine && !engine) continue;
			for (int stiff = 0; stiff < 2; ++stiff)
				for (int ih = 0; ih < kNH; ++ih) {
					const RunResult& E = runs[Index(kE, sh, bend, stiff, ih)];
					const RunResult& G = runs[Index(kGN, sh, bend, stiff, ih)];
					const RunResult& T = runs[Index(kTen, sh, bend, stiff, ih)];
					++pairs;
					const bool uE = Usable(E), uG = Usable(G), uT = Usable(T);
					nE += uE; nG += uG; nT += uT;
					if (uE && !uT) { ++viol; notes.push_back("(가′) E 쓸 수 있음 · T 못 씀: " + Desc(T.cfg)); }
					if (uG && !uT) { ++gNotT; notes.push_back("보고 — G 쓸 수 있음 · T 못 씀: " + Desc(T.cfg)); }
					if (ih == 0) continue;
					const RunResult& H = runs[Index(kE, sh, bend, stiff, ih - 1)];
					if (!(uT && uE && Usable(H))) continue;
					const double dTE = Dist(T, E, nNode), dEH = Dist(E, H, nNode);
					++diffPairs;
					worst = std::max(worst, dTE / dEH);
					if (dTE > dEH) {
						++diffFails;
						char s[160];
						std::snprintf(s, sizeof s, "(나′) D(T,E) %.2e > D(E,E/2) %.2e: %s", dTE, dEH, Desc(T.cfg).c_str());
						notes.push_back(s);
					}
				}
		}
	const bool aOk = viol == 0, bOk = diffPairs > 0 && diffFails == 0, dOk = nT > nE;
	const RunResult& refE = runs[Index(kE, kE, kBridsonG, 1, 1)];
	const RunResult& refT = runs[Index(kTen, kE, kBridsonG, 1, 1)];
	const RunResult& defT = runs[Index(kTen, kE, kBridsonG, 0, 1)];
	std::printf("인장 C > 0 판 — 짝 %d (인장 E 대 T)\n", pairs);
	std::printf("  (가′) %s  위반 %d\n", aOk ? "넘음" : "못 넘음", viol);
	std::printf("  (나′) %s  견준 짝 %d, 넘은 짝 %d, 최대 D(T,E)/D(E,E/2) %.4f\n", bOk ? "넘음" : "못 넘음", diffPairs, diffFails, worst);
	std::printf("  (라′) %s  쓸 수 있는 회 T %d · E %d (보고: G %d)\n", dOk ? "넘음" : "못 넘음", nT, nE, nG);
	std::printf("  (다′) 기준 짝(전단 E · Bridson G · k 1e-3 · 1/60) 쓸 수 있음 E %s · T %s — 시간은 tests/timing_bases.py tension\n",
		Usable(refE) ? "예" : "아니오", Usable(refT) ? "예" : "아니오");
	std::printf("  보고 — 기본 조합(전단 E · Bridson G · 부드러움 · 1/60) T %s, G 는 쓸 수 있는데 T 는 못 쓰는 짝 %d\n",
		Usable(defT) ? "쓸 수 있음" : "못 씀", gNotT);
	for (std::size_t i = 0; i < notes.size() && i < 20; ++i) std::printf("      %s\n", notes[i].c_str());
	if (notes.size() > 20) std::printf("      … 외 %zu 개\n", notes.size() - 20);

	std::FILE* js = std::fopen("results/gate_tension_integration.json", "w");
	if (!js) return 2;
	std::fprintf(js, "{\n  \"scenario\": {\"G\": %d, \"L\": %g, \"rho\": %g, \"g\": %g, \"stretch_k\": %g, \"stretch_kd\": %g, \"shear_k\": %g, "
		"\"shear_kd\": %g, \"bend_k\": [%g, %g], \"T\": %g, \"steps_per_s\": [%d, %d, %d], \"pcg_eps\": %g, \"pcg_cap\": %d, \"v_blow\": %g, "
		"\"usable_strain\": %g},\n",
		kG, kL, kRho, kGrav, kStretchK, kStretchKd, kShearK, kShearKd, kBendK[0], kBendK[1], kT, kStepsPerSec[0], kStepsPerSec[1],
		kStepsPerSec[2], kEps, kCap, kVBlow, kUsableStrain);
	std::fprintf(js, "  \"engine\": %s,\n  \"pairs\": %d,\n", engine ? "true" : "false", pairs);
	std::fprintf(js, "  \"rules\": {\"A\": {\"violations\": %d, \"ok\": %s}, \"B\": {\"diff_pairs\": %d, \"diff_fails\": %d, \"worst_ratio\": %.4f, \"ok\": %s}, "
		"\"D\": {\"usable_T\": %d, \"usable_E\": %d, \"usable_G\": %d, \"ok\": %s}},\n",
		viol, aOk ? "true" : "false", diffPairs, diffFails, worst, bOk ? "true" : "false", nT, nE, nG, dOk ? "true" : "false");
	std::fprintf(js, "  \"ref_stiff\": {\"E_usable\": %s, \"T_usable\": %s},\n  \"default_soft_T_usable\": %s,\n  \"g_usable_t_not\": %d,\n",
		Usable(refE) ? "true" : "false", Usable(refT) ? "true" : "false", Usable(defT) ? "true" : "false", gNotT);
	std::fprintf(js, "  \"notes\": [");
	for (std::size_t i = 0; i < notes.size(); ++i) std::fprintf(js, "%s\"%s\"", i ? ", " : "", notes[i].c_str());
	std::fprintf(js, "],\n  \"distances_over_L\": [\n");
	{
		bool first = true;
		for (int sh = 0; sh < 2; ++sh)
			for (int bend = 0; bend < 5; ++bend) {
				if (bend == kEngine && !engine) continue;
				for (int stiff = 0; stiff < 2; ++stiff)
					for (int ih = 0; ih < kNH; ++ih) {
						const RunResult& T = runs[Index(kTen, sh, bend, stiff, ih)];
						const double dE = Dist(T, runs[Index(kE, sh, bend, stiff, ih)], nNode) / kL;
						const double dG = Dist(T, runs[Index(kGN, sh, bend, stiff, ih)], nNode) / kL;
						const double dH = ih > 0 ? Dist(runs[Index(kE, sh, bend, stiff, ih)], runs[Index(kE, sh, bend, stiff, ih - 1)], nNode) / kL : -1.;
						std::fprintf(js, "%s    {\"run\": \"%s\", \"kb\": %g, \"steps_per_s\": %d, \"D_T_vs_E\": %.4e, \"D_T_vs_G\": %.4e, \"D_E_vs_half_h\": %.4e}",
							first ? "" : ",\n", Key(T.cfg).c_str(), kBendK[stiff], kStepsPerSec[ih], dE, dG, dH);
						first = false;
					}
			}
		std::fprintf(js, "\n  ],\n");
	}
	std::fprintf(js, "  \"runs\": [\n");
	{
		std::vector<const RunResult*> list;
		for (std::size_t i = 0; i < runs.size(); ++i) if (runs[i].ran) list.push_back(&runs[i]);
		for (std::size_t i = 0; i < list.size(); ++i) JsonRun(js, *list[i], i + 1 < list.size() ? "," : "");
	}
	std::fprintf(js, "  ]\n}\n");
	std::fclose(js);
	std::printf("\n인장 C > 0 판 적분 → results/gate_tension_integration.json (권고는 tests/timing_bases.py tension)\n");
	return 0;
}

// 기준 짝만 돌려 시간 · 조립 시간 · PCG 반복 · 쓸 수 있음을 낸다(판정은 tests/timing_bases.py)
static int TimingRefs(Sim& sim, const Views& vw, const Config* refs, int n, const char* path)
{
	std::vector<RunResult> rs;
	for (int i = 0; i < n; ++i) { rs.push_back(Run(sim, vw, refs[i], false, -1, 0)); PrintRun(rs.back()); }
	std::FILE* js = std::fopen(path, "w");
	if (!js) return 2;
	const double perS = 1e3 / kT;
	const char* blocks[] = { "ms_per_sim_second", "asm_ms_per_sim_second", "pcg_iters", "usable" };
	std::fprintf(js, "{\n");
	for (int k = 0; k < 4; ++k) {
		std::fprintf(js, "  \"%s\": {", blocks[k]);
		for (int i = 0; i < n; ++i) {
			const RunResult& r = rs[i];
			const double v = k == 0 ? (r.tAsm + r.tPcg) * perS : k == 1 ? r.tAsm * perS : k == 2 ? static_cast<double>(r.iters) : (Usable(r) ? 1. : 0.);
			std::fprintf(js, "%s\"%s%s\": %.3f", i ? ", " : "", Key(r.cfg).c_str(), r.cfg.stiff ? "_stiff" : "", v);
		}
		std::fprintf(js, "}%s\n", k < 3 ? "," : "");
	}
	std::fprintf(js, "}\n");
	std::fclose(js);
	return 0;
}

// (다′): 뻣뻣한 굽힘의 기준 짝(판정) + 부드러운 굽힘의 같은 짝(보고)
static int TimingTension(Sim& sim, const Views& vw)
{
	const Config refs[] = { { kE, kE, kBridsonG, 1, 1, 0 }, { kTen, kE, kBridsonG, 1, 1, 0 }, { kE, kE, kBridsonG, 0, 1, 0 }, { kTen, kE, kBridsonG, 0, 1, 0 } };
	return TimingRefs(sim, vw, refs, sizeof refs / sizeof refs[0], "results/gate_tension_time.json");
}

// ---- 전단 GN 재비교 (docs/GATE_SHEAR_GN.md) — 인장은 모두 T --------------------------------------
static int ShearGnMode(Sim& sim, const Views& vw, int nNode, bool engine)
{
	std::vector<RunResult> runs(3 * 2 * 5 * 2 * kNH);   // Index 는 인장 셋을 다 담는다. 여기서는 인장 T 만 채운다
	const Clock::time_point tStart = Clock::now();
	for (int bend = 0; bend < 5; ++bend) {
		if (bend == kEngine && !engine) continue;
		for (int sh = 0; sh < 2; ++sh)
			for (int stiff = 0; stiff < 2; ++stiff)
				for (int ih = 0; ih < kNH; ++ih) {
					const Config cf = { kTen, sh, bend, stiff, ih, 0 };
					runs[Index(kTen, sh, bend, stiff, ih)] = Run(sim, vw, cf, true, -1, 0);
					PrintRun(runs[Index(kTen, sh, bend, stiff, ih)]);
				}
	}
	std::printf("  (적분 %.1f s)\n", Sec(tStart, Clock::now()));

	std::printf("쓸 수 있음 표 (O · X 변형률 > %.2f · D 발산 · F 풀이 실패), 칸 = h 1/120 · 1/60 · 1/30\n", kUsableStrain);
	for (int sh = 0; sh < 2; ++sh) {
		std::printf("  T%s ", kModeName[sh]);
		for (int stiff = 0; stiff < 2; ++stiff) {
			std::printf(" | k%-5g", kBendK[stiff]);
			for (int bend = 0; bend < 5; ++bend) {
				if (bend == kEngine && !engine) continue;
				std::printf(" %s ", kBendName[bend]);
				for (int ih = 0; ih < kNH; ++ih) std::printf("%c", Code(runs[Index(kTen, sh, bend, stiff, ih)]));
			}
		}
		std::printf("\n");
	}

	// (가″) · (나″) · (라″) 와 보고
	int pairs = 0, viol = 0, diffPairs = 0, diffFails = 0, nE = 0, nG = 0, gNotE = 0;
	double worst = 0.;
	std::vector<std::string> notes;
	for (int bend = 0; bend < 5; ++bend) {
		if (bend == kEngine && !engine) continue;
		for (int stiff = 0; stiff < 2; ++stiff)
			for (int ih = 0; ih < kNH; ++ih) {
				const RunResult& E = runs[Index(kTen, kE, bend, stiff, ih)];
				const RunResult& G = runs[Index(kTen, kGN, bend, stiff, ih)];
				++pairs;
				const bool uE = Usable(E), uG = Usable(G);
				nE += uE; nG += uG;
				if (uE && !uG) { ++viol; notes.push_back("(가″) E 쓸 수 있음 · G 못 씀: " + Desc(G.cfg)); }
				if (uG && !uE) { ++gNotE; notes.push_back("보고 — G 쓸 수 있음 · E 못 씀: " + Desc(G.cfg)); }
				if (ih == 0) continue;
				const RunResult& H = runs[Index(kTen, kE, bend, stiff, ih - 1)];
				if (!(uG && uE && Usable(H))) continue;
				const double dGE = Dist(G, E, nNode), dEH = Dist(E, H, nNode);
				++diffPairs;
				worst = std::max(worst, dGE / dEH);
				if (dGE > dEH) {
					++diffFails;
					char s[160];
					std::snprintf(s, sizeof s, "(나″) D(G,E) %.2e > D(E,E/2) %.2e: %s", dGE, dEH, Desc(G.cfg).c_str());
					notes.push_back(s);
				}
			}
	}
	const bool aOk = viol == 0, bOk = diffPairs > 0 && diffFails == 0, dOk = nG > nE;
	const RunResult& softE = runs[Index(kTen, kE, kBridsonG, 0, 1)];
	const RunResult& softG = runs[Index(kTen, kGN, kBridsonG, 0, 1)];
	const RunResult& stiffE = runs[Index(kTen, kE, kBridsonG, 1, 1)];
	const RunResult& stiffG = runs[Index(kTen, kGN, kBridsonG, 1, 1)];
	std::printf("전단 GN 재비교 (인장 T) — 짝 %d (전단 E 대 G)\n", pairs);
	std::printf("  (가″) %s  위반 %d\n", aOk ? "넘음" : "못 넘음", viol);
	std::printf("  (나″) %s  견준 짝 %d, 넘은 짝 %d, 최대 D(G,E)/D(E,E/2) %.4f\n", bOk ? "넘음" : "못 넘음", diffPairs, diffFails, worst);
	std::printf("  (라″) %s  쓸 수 있는 회 G %d · E %d (G 만 쓸 수 있는 짝 %d)\n", dOk ? "넘음" : "못 넘음", nG, nE, gNotE);
	std::printf("  (다″) 기준 짝 쓸 수 있음 — ① 부드러움 E %s · G %s, ② 뻣뻣함 E %s · G %s — 시간은 tests/timing_bases.py shear_gn\n",
		Usable(softE) ? "예" : "아니오", Usable(softG) ? "예" : "아니오", Usable(stiffE) ? "예" : "아니오", Usable(stiffG) ? "예" : "아니오");
	for (std::size_t i = 0; i < notes.size() && i < 20; ++i) std::printf("      %s\n", notes[i].c_str());
	if (notes.size() > 20) std::printf("      … 외 %zu 개\n", notes.size() - 20);

	std::FILE* js = std::fopen("results/gate_shear_gn_integration.json", "w");
	if (!js) return 2;
	std::fprintf(js, "{\n  \"scenario\": {\"G\": %d, \"L\": %g, \"rho\": %g, \"g\": %g, \"stretch\": \"T\", \"stretch_k\": %g, \"stretch_kd\": %g, \"shear_k\": %g, "
		"\"shear_kd\": %g, \"bend_k\": [%g, %g], \"T\": %g, \"steps_per_s\": [%d, %d, %d], \"pcg_eps\": %g, \"pcg_cap\": %d, \"v_blow\": %g, "
		"\"usable_strain\": %g},\n",
		kG, kL, kRho, kGrav, kStretchK, kStretchKd, kShearK, kShearKd, kBendK[0], kBendK[1], kT, kStepsPerSec[0], kStepsPerSec[1],
		kStepsPerSec[2], kEps, kCap, kVBlow, kUsableStrain);
	std::fprintf(js, "  \"engine\": %s,\n  \"pairs\": %d,\n", engine ? "true" : "false", pairs);
	std::fprintf(js, "  \"rules\": {\"A\": {\"violations\": %d, \"ok\": %s}, \"B\": {\"diff_pairs\": %d, \"diff_fails\": %d, \"worst_ratio\": %.4f, \"ok\": %s}, "
		"\"D\": {\"usable_G\": %d, \"usable_E\": %d, \"ok\": %s}},\n",
		viol, aOk ? "true" : "false", diffPairs, diffFails, worst, bOk ? "true" : "false", nG, nE, dOk ? "true" : "false");
	std::fprintf(js, "  \"ref_soft\": {\"E_usable\": %s, \"G_usable\": %s},\n  \"ref_stiff\": {\"E_usable\": %s, \"G_usable\": %s},\n  \"g_usable_e_not\": %d,\n",
		Usable(softE) ? "true" : "false", Usable(softG) ? "true" : "false", Usable(stiffE) ? "true" : "false", Usable(stiffG) ? "true" : "false", gNotE);
	std::fprintf(js, "  \"notes\": [");
	for (std::size_t i = 0; i < notes.size(); ++i) std::fprintf(js, "%s\"%s\"", i ? ", " : "", notes[i].c_str());
	std::fprintf(js, "],\n  \"distances_over_L\": [\n");
	{
		bool first = true;
		for (int bend = 0; bend < 5; ++bend) {
			if (bend == kEngine && !engine) continue;
			for (int stiff = 0; stiff < 2; ++stiff)
				for (int ih = 0; ih < kNH; ++ih) {
					const RunResult& G = runs[Index(kTen, kGN, bend, stiff, ih)];
					const RunResult& E = runs[Index(kTen, kE, bend, stiff, ih)];
					const double dGE = Dist(G, E, nNode) / kL;
					const double dH = ih > 0 ? Dist(E, runs[Index(kTen, kE, bend, stiff, ih - 1)], nNode) / kL : -1.;
					std::fprintf(js, "%s    {\"bend\": \"%s\", \"kb\": %g, \"steps_per_s\": %d, \"G\": \"%c\", \"E\": \"%c\", \"D_G_vs_E\": %.4e, \"D_E_vs_half_h\": %.4e}",
						first ? "" : ",\n", kBendName[bend], kBendK[stiff], kStepsPerSec[ih], Code(G), Code(E), dGE, dH);
					first = false;
				}
		}
		std::fprintf(js, "\n  ],\n");
	}
	std::fprintf(js, "  \"runs\": [\n");
	{
		std::vector<const RunResult*> list;
		for (std::size_t i = 0; i < runs.size(); ++i) if (runs[i].ran) list.push_back(&runs[i]);
		for (std::size_t i = 0; i < list.size(); ++i) JsonRun(js, *list[i], i + 1 < list.size() ? "," : "");
	}
	std::fprintf(js, "  ]\n}\n");
	std::fclose(js);
	std::printf("\n전단 GN 재비교 적분 → results/gate_shear_gn_integration.json (권고는 tests/timing_bases.py shear_gn)\n");
	return 0;
}

// (다″): 기준 짝 둘 — ① 부드러운 굽힘 · Bridson G · 1/60, ② 뻣뻣한 굽힘 · Bridson G · 1/60, 인장은 T
static int TimingShearGn(Sim& sim, const Views& vw)
{
	const Config refs[] = { { kTen, kE, kBridsonG, 0, 1, 0 }, { kTen, kGN, kBridsonG, 0, 1, 0 }, { kTen, kE, kBridsonG, 1, 1, 0 }, { kTen, kGN, kBridsonG, 1, 1, 0 } };
	return TimingRefs(sim, vw, refs, sizeof refs / sizeof refs[0], "results/gate_shear_gn_time.json");
}

int main(int argc, char** argv)
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const Cloth cloth = MakeCloth();
	Sim sim(cloth);
	const Views vw(sim.mesh);
#ifdef SX_WITH_ENGINE
	const bool engine = true;
#else
	const bool engine = false;
#endif
	std::printf("시간 적분 비교 — 노드 %d (고정 2), 삼각형 %d, 힌지 %d, 노드쌍 블록 %d, T = %g s%s\n",
		cloth.nNode, cloth.nTri, cloth.nHinge, cloth.nPair, kT, engine ? ", 엔진 Bridson 포함" : ", 엔진 없음");
	if (argc > 1 && std::strcmp(argv[1], "timing") == 0) return TimingMode(sim, vw);
	if (argc > 1 && std::strcmp(argv[1], "tension") == 0) return TensionMode(sim, vw, cloth.nNode, engine);
	if (argc > 1 && std::strcmp(argv[1], "timing_tension") == 0) return TimingTension(sim, vw);
	if (argc > 1 && std::strcmp(argv[1], "shear_gn") == 0) return ShearGnMode(sim, vw, cloth.nNode, engine);
	if (argc > 1 && std::strcmp(argv[1], "timing_shear_gn") == 0) return TimingShearGn(sim, vw);

	// ---- 120 회 (엔진 없으면 96) ---------------------------------------------------------
	std::vector<RunResult> runs(2 * 2 * 5 * 2 * kNH);
	std::vector<double> x05;
	const Clock::time_point tStart = Clock::now();
	for (int bend = 0; bend < 5; ++bend) {
		if (bend == kEngine && !engine) continue;
		for (int st = 0; st < 2; ++st)
			for (int sh = 0; sh < 2; ++sh)
				for (int stiff = 0; stiff < 2; ++stiff)
					for (int ih = 0; ih < kNH; ++ih) {
						const Config cf = { st, sh, bend, stiff, ih, 0 };
						const bool d0 = st == kE && sh == kE && bend == kBridsonG && stiff == 0 && ih == 1;
						const int keepAt = static_cast<int>(0.5 * kStepsPerSec[ih] + 0.5);
						runs[Index(st, sh, bend, stiff, ih)] = Run(sim, vw, cf, true, d0 ? keepAt : -1, d0 ? &x05 : 0);
						PrintRun(runs[Index(st, sh, bend, stiff, ih)]);
					}
	}
	std::printf("  (적분 %.1f s)\n", Sec(tStart, Clock::now()));

	// ---- I0 -------------------------------------------------------------------------
	Report rpt("results/gate_integration.json");
	std::FILE* js = rpt.json;
	std::printf("I0 전역 조립 대 유한차분 (기본 조합의 t = 0.5 s, v = 0, 정확판, 힘 종류 하나씩)\n");
	Terms t0;
	t0.kb = kBendK[0]; t0.h = StepH(1); t0.mass = false;
	double i0[4];
	t0.stretch = kE; t0.shear = -1; t0.bend = -1; i0[0] = CheckI0(sim, vw, x05, t0, "인장");
	t0.stretch = -1; t0.shear = kE; i0[1] = CheckI0(sim, vw, x05, t0, "전단");
	t0.shear = -1; t0.bend = kBaraffE; i0[2] = CheckI0(sim, vw, x05, t0, "Baraff");
	t0.bend = kBridsonE; i0[3] = CheckI0(sim, vw, x05, t0, "Bridson");
	const double i0max = std::max(std::max(i0[0], i0[1]), std::max(i0[2], i0[3]));
	rpt.Gate("I0", i0max <= kI0Limit, "max %.2e", i0max, kI0Limit);

	// ---- 권고 규칙 (가) · (나) -------------------------------------------------------------
	const char* label[4] = { "인장", "전단", "Baraff 굽힘", "Bridson 굽힘" };
	Rule rules[4];
	std::printf("권고 규칙 (가) 안정 · (나) 차이 — (다) 시간은 tests/timing_bases.py integ\n");
	for (int w = 0; w < 4; ++w) {
		rules[w] = MakeRule(w, runs, cloth.nNode, engine);
		const Rule& r = rules[w];
		std::printf("  %-12s 짝 %d: (가) %s (위반 %d) · (나) %s (견준 짝 %d, 넘은 짝 %d, 최대 D(G,E)/D(E,E/2) %.3f) · 기준 짝 안정 E %s G %s\n",
			label[w], r.pairs, r.A() ? "넘음" : "못 넘음", r.violations, r.B() ? "넘음" : "못 넘음", r.diffPairs, r.diffFails, r.worstRatio,
			r.refE ? "예" : "아니오", r.refG ? "예" : "아니오");
		for (std::size_t i = 0; i < r.notes.size() && i < 12; ++i) std::printf("      %s\n", r.notes[i].c_str());
		if (r.notes.size() > 12) std::printf("      … 외 %zu 개\n", r.notes.size() - 12);
	}

	// ---- 보고: 엔진 대 생성 Bridson ---------------------------------------------------------
	if (engine) {
		std::printf("보고 — 엔진 Bridson 대 생성 Bridson (D / L)\n");
		for (int st = 0; st < 2; ++st) for (int sh = 0; sh < 2; ++sh) for (int stiff = 0; stiff < 2; ++stiff) for (int ih = 0; ih < kNH; ++ih) {
			const RunResult& en = runs[Index(st, sh, kEngine, stiff, ih)];
			const RunResult& be = runs[Index(st, sh, kBridsonE, stiff, ih)];
			const RunResult& bg = runs[Index(st, sh, kBridsonG, stiff, ih)];
			std::printf("  %s%s k%-5g h1/%-3d 엔진 %s · E %s · G %s  D(엔진,E) %.2e  D(엔진,G) %.2e  D(G,E) %.2e\n",
				kModeName[st], kModeName[sh], kBendK[stiff], kStepsPerSec[ih], en.stable ? "안정" : "불안정", be.stable ? "안정" : "불안정",
				bg.stable ? "안정" : "불안정", Dist(en, be, cloth.nNode) / kL, Dist(en, bg, cloth.nNode) / kL, Dist(bg, be, cloth.nNode) / kL);
		}
	}

	// ---- JSON -----------------------------------------------------------------------
	std::fprintf(js, "{\n  \"scenario\": {\"G\": %d, \"L\": %g, \"rho\": %g, \"g\": %g, \"stretch_k\": %g, \"stretch_kd\": %g, \"shear_k\": %g, "
		"\"shear_kd\": %g, \"bend_k\": [%g, %g], \"T\": %g, \"steps_per_s\": [%d, %d, %d], \"pcg_eps\": %g, \"pcg_cap\": %d, \"v_blow\": %g},\n",
		kG, kL, kRho, kGrav, kStretchK, kStretchKd, kShearK, kShearKd, kBendK[0], kBendK[1], kT, kStepsPerSec[0], kStepsPerSec[1],
		kStepsPerSec[2], kEps, kCap, kVBlow);
	std::fprintf(js, "  \"mesh\": {\"nodes\": %d, \"triangles\": %d, \"hinges\": %d, \"pairs\": %d},\n", cloth.nNode, cloth.nTri, cloth.nHinge, cloth.nPair);
	std::fprintf(js, "  \"engine\": %s,\n", engine ? "true" : "false");
	std::fprintf(js, "  \"I0\": {\"stretch\": %.3e, \"shear\": %.3e, \"baraff\": %.3e, \"bridson\": %.3e, \"limit\": %g, \"pass\": %s},\n",
		i0[0], i0[1], i0[2], i0[3], kI0Limit, i0max <= kI0Limit ? "true" : "false");
	std::fprintf(js, "  \"rules\": {\n");
	for (int w = 0; w < 4; ++w) JsonRule(js, rules[w], w < 3 ? "," : "");
	std::fprintf(js, "  },\n  \"distances_over_L\": [\n");
	{   // 짝의 거리와 시간 간격 척도(보고)
		bool first = true;
		for (std::size_t i = 0; i < runs.size(); ++i) {
			const RunResult& r = runs[i];
			if (!r.ran) continue;
			const Config& cf = r.cfg;
			double dE = -1., dHalf = -1., dEng = -1.;
			int iE = -1;
			if (cf.bend == kBaraffG) iE = Index(cf.stretch, cf.shear, kBaraffE, cf.stiff, cf.ih);
			else if (cf.bend == kBridsonG || cf.bend == kEngine) iE = Index(cf.stretch, cf.shear, kBridsonE, cf.stiff, cf.ih);
			if (iE >= 0 && runs[iE].ran) dE = Dist(r, runs[iE], cloth.nNode) / kL;
			if (cf.ih > 0) dHalf = Dist(r, runs[Index(cf.stretch, cf.shear, cf.bend, cf.stiff, cf.ih - 1)], cloth.nNode) / kL;
			if (cf.stretch == kGN) dEng = Dist(r, runs[Index(kE, cf.shear, cf.bend, cf.stiff, cf.ih)], cloth.nNode) / kL;
			std::fprintf(js, "%s    {\"run\": \"%s\", \"kb\": %g, \"steps_per_s\": %d, \"D_bend_vs_E\": %.4e, \"D_vs_half_h\": %.4e, \"D_stretchG_vs_E\": %.4e",
				first ? "" : ",\n", Key(cf).c_str(), kBendK[cf.stiff], Sps(cf), dE, dHalf, dEng);
			const double dSh = cf.shear == kGN ? Dist(r, runs[Index(cf.stretch, kE, cf.bend, cf.stiff, cf.ih)], cloth.nNode) / kL : -1.;
			std::fprintf(js, ", \"D_shearG_vs_E\": %.4e}", dSh);
			first = false;
		}
		std::fprintf(js, "\n  ],\n");
	}
	std::fprintf(js, "  \"runs\": [\n");
	{
		std::vector<const RunResult*> list;
		for (std::size_t i = 0; i < runs.size(); ++i) if (runs[i].ran) list.push_back(&runs[i]);
		for (std::size_t i = 0; i < list.size(); ++i) JsonRun(js, *list[i], i + 1 < list.size() ? "," : "");
	}
	std::fprintf(js, "  ],\n  \"gates_failed\": %d\n}\n", rpt.fail);
	std::printf("\nI0 %s (실패 %d) → results/gate_integration.json\n", rpt.fail == 0 ? "PASS — 하네스가 맞다" : "FAIL — 하네스를 고친다", rpt.fail);
	return rpt.fail == 0 ? 0 : 1;
}
