//---------------------------------------------------------------------------
// test_link_l1.cpp — 엔진 연결 L1(docs/GATE_LINK_L1.md): 엔진용 생성 변형 EvalEngine* 의 정의 · 행 합 복원 · 전역 곱.
// L1-1: 노드쌍 J · V 블록과 f + d, C 를 평탄판 Eval 의 C · ∇C · ∇²C 로 만든 기준 요소 행렬과 견준다.
// L1-2: 블록과 행 합 규약으로 되살린 A = −hV − h²J, b = h f + h² J v 를 짝이 되는 EvalImplicit* 와 견준다.
// L1-3: 격자 천에서 블록을 노드쌍 슬롯에 모아 엔진 곱 연산(SDESparseMatrix.cpp:4256-4308, 전치 처리 켬)을 흉내 내고,
//       기준 요소 행렬을 직접 곱해 흩뿌린 값과 견준다. 전치 처리 끔은 보고만 한다.
// L1-4: 감쇠 없는 GN EvalEngineGNUndamped(x, ev, {k}) 의 C · f · J 가 EvalEngineGN(x, ev, v = 0, {k, 0}) 과 같다
//       (2026-11-13 추가, docs/GATE_LINK_L4.md 「2 차 (다)」 — 엔진 굽힘 갈래가 이것을 부른다).
// L1-5: 엔진 처방 T 판 EvalEngineTensionRate(인장 · Provot 전단) — J 에 −kd Σ max(Ċ, 0) ∇²C 가 더 든다. 기준 요소 행렬 · 전역 곱과
//       견주고(L1-1 · L1-3 꼴), v = 0 에서 EvalEngineTension 과 같은지 본다(L1-5b). 짝이 되는 EvalImplicit* 가 없어 L1-2 는 없다.
//       (2026-11-13 추가 — 엔진 연결 L3 의 준비. Provot 전단 명세 ProvotShear 도 여섯째 명세로 넣었다.)
// L1-0(기존 출력 불변)은 tests/check_link_l1_p0.py.
// 실행 디렉터리 = 리포 루트. results/hinges.txt · triangles.txt 를 읽고 results/gate_link_l1.json 을 쓴다.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffBend.h"
#include "../generated/BaraffBendCompact.h"
#include "../generated/BaraffStretch.h"
#include "../generated/BaraffShear.h"
#include "../generated/BridsonBend.h"
#include "../generated/ProvotShear.h"
#include "../generated/BaraffBend_Engine.h"
#include "../generated/BaraffBendCompact_Engine.h"
#include "../generated/BaraffStretch_Engine.h"
#include "../generated/BaraffShear_Engine.h"
#include "../generated/BridsonBend_Engine.h"
#include "../generated/ProvotShear_Engine.h"
#include "gate_common.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;

static const double kPrmK = 1.0, kPrmKd = 0.1, kPrmH = 0.01;   // Gauss–Newton 게이트 N1 과 같다
static const double kPrm3[3] = { kPrmK, kPrmKd, kPrmH };        // EvalImplicit*
static const double kPrm2[2] = { kPrmK, kPrmKd };               // EvalEngine*
static const double kLimit = 1e-12;
static const int kGridG = 12;
static const int kGridVectors = 3;

static const int kMinBranch = 100;                             // L1-5 갈래

enum Kind { kExact, kGN, kTension, kTensionRate };
static const char* KindName(Kind k) { return k == kExact ? "exact" : k == kGN ? "gn" : k == kTension ? "tension" : "tension_rate"; }

// EvalEngine* 와 EvalImplicit* 는 인자 꼴이 같다: (x, ev, v, prm, C, f, J|A, V|b)
typedef void (*Fn8)(const double*, const double*, const double*, const double*, double*, double*, double*, double*);
struct Variant { Kind kind; Fn8 eng; Fn8 imp; };   // imp = 0 이면 짝이 없다(L1-2 를 하지 않는다)

struct Sample { std::vector<double> x, ev, v; double cscale; };

// Bridson 표본 가중치 w̄ = 3|e|² / (½|nA| + ½|nB|) (Gauss–Newton 게이트와 같은 식)
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

// 기준 요소 행렬(전체 N × N, 행 우선): 평탄판 Eval 의 C · ∇C · ∇²C 에서.
// 정확판 = AssembleConditionForce 의 K + ∂d/∂x(비대칭 항 버림) · ∂d/∂v, GN = 같은 함수에 ∇²C = 0,
// T = GN 에 −k Σ C⁺ ∇²C 를 더한 것, 엔진 처방 T = 거기에 −kd Σ Ċ⁺ ∇²C 를 더 더한 것. f = f + d. Cdot 이 있으면 Ċ = ∇C · v 를 낸다.
template <class S>
static void RefMatrices(Kind kind, const double* x, const double* ev, const double* v, double* C, double* Jr, double* Vr, double* fr,
	double* Cdot = 0)
{
	enum { N = S::kDof, NH = S::kHess, NC = S::kComp };
	double g[NC * N], H[NC * NH], H0[NC * NH];
	S::Eval(x, ev, C, g, H);
	for (int i = 0; i < NC * NH; ++i) H0[i] = 0.;
	sx::ConditionForce<N> cf;
	sx::AssembleConditionForce<N>(NC, C, g, kind == kExact ? H : H0, v, kPrmK, kPrmKd, false, cf);
	for (int i = 0; i < N * N; ++i) { Jr[i] = cf.K[i] + cf.Dx[i]; Vr[i] = cf.Dv[i]; }
	for (int a = 0; a < N; ++a) fr[a] = cf.f[a] + cf.d[a];
	for (int c = 0; c < NC; ++c) {
		double cd = 0.;
		for (int a = 0; a < N; ++a) cd += g[c * N + a] * v[a];
		if (Cdot) Cdot[c] = cd;
		if (kind != kTension && kind != kTensionRate) continue;
		const double cp = C[c] > 0. ? C[c] : 0.;
		const double cdp = (kind == kTensionRate && cd > 0.) ? cd : 0.;
		const double beta = kind == kTensionRate ? -kPrmK * cp - kPrmKd * cdp : -kPrmK * cp;
		for (int a = 0; a < N; ++a)
			for (int b = 0; b < N; ++b) Jr[a * N + b] += beta * H[c * NH + sx::UpperIndex<N>(a, b)];
	}
}

// 노드쌍 블록 → 전체 요소 행렬(행 합 규약): (b, a) = (a, b)ᵀ, 대각 블록 = −(그 행의 이웃 블록 합)
template <class E>
static void FullFromPairs(const double* B, double* F)
{
	enum { N = E::kDof, NN = E::kNodes, NP = E::kPairs };
	for (int i = 0; i < N * N; ++i) F[i] = 0.;
	for (int p = 0; p < NP; ++p) {
		const int a = E::PairA(p), b = E::PairB(p);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j) {
				F[(3 * a + i) * N + 3 * b + j] = B[9 * p + 3 * i + j];
				F[(3 * b + j) * N + 3 * a + i] = B[9 * p + 3 * i + j];
			}
	}
	for (int a = 0; a < NN; ++a)
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j) {
				double s = 0.;
				for (int b = 0; b < NN; ++b)
					if (b != a) s += F[(3 * a + i) * N + 3 * b + j];
				F[(3 * a + i) * N + 3 * a + j] = -s;
			}
}

struct Acc {
	Stat c, J, V, f, A, b;     // L1-1 (c · J · V · f), L1-2 (A · b)
	Stat gJ, gV, gJoff, gVoff; // L1-3 (전치 처리 켬), 보고(끔)
	int n, nPos, nNeg, nRatePos, nRateNeg;   // 성분 C > 0 · C ≤ 0, Ċ > 0 · Ċ ≤ 0
	Acc() : n(0), nPos(0), nNeg(0), nRatePos(0), nRateNeg(0) {}
};

template <class S, class E>
static void CheckSample(const Variant& var, const Sample& sm, Acc& acc)
{
	enum { N = S::kDof, NH = S::kHess, NC = S::kComp, NP = E::kPairs };
	double C[NC], Jr[N * N], Vr[N * N], fr[N], Cdot[NC];
	RefMatrices<S>(var.kind, &sm.x[0], &sm.ev[0], &sm.v[0], C, Jr, Vr, fr, Cdot);
	double Ce[NC], fe[N], Je[9 * NP], Ve[9 * NP];
	var.eng(&sm.x[0], &sm.ev[0], &sm.v[0], kPrm2, Ce, fe, Je, Ve);
	// L1-1
	double dJ = 0., dV = 0.;
	for (int p = 0; p < NP; ++p) {
		const int a = E::PairA(p), b = E::PairB(p);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j) {
				dJ = std::max(dJ, std::fabs(Je[9 * p + 3 * i + j] - Jr[(3 * a + i) * N + 3 * b + j]));
				dV = std::max(dV, std::fabs(Ve[9 * p + 3 * i + j] - Vr[(3 * a + i) * N + 3 * b + j]));
			}
	}
	acc.c.Add(MaxAbsDiff(Ce, C, NC) / sm.cscale);
	acc.J.Add(Rel(dJ, MaxAbs(Jr, N * N)));
	acc.V.Add(Rel(dV, MaxAbs(Vr, N * N)));
	acc.f.Add(Rel(MaxAbsDiff(fe, fr, N), MaxAbs(fr, N)));
	for (int c = 0; c < NC; ++c) {
		if (C[c] > 0.) ++acc.nPos; else ++acc.nNeg;
		if (Cdot[c] > 0.) ++acc.nRatePos; else ++acc.nRateNeg;
	}
	++acc.n;
	if (!var.imp) return;   // 짝이 없는 변형(엔진 처방 T)
	// L1-2
	double Jf[N * N], Vf[N * N];
	FullFromPairs<E>(Je, Jf);
	FullFromPairs<E>(Ve, Vf);
	double Arec[NH], brec[N];
	for (int r = 0; r < N; ++r) {
		for (int c = r; c < N; ++c) Arec[sx::UpperIndex<N>(r, c)] = -kPrmH * Vf[r * N + c] - kPrmH * kPrmH * Jf[r * N + c];
		double s = 0.;
		for (int c = 0; c < N; ++c) s += Jf[r * N + c] * sm.v[c];
		brec[r] = kPrmH * fe[r] + kPrmH * kPrmH * s;
	}
	double Ci[NC], fi[N], Ai[NH], bi[N];
	var.imp(&sm.x[0], &sm.ev[0], &sm.v[0], kPrm3, Ci, fi, Ai, bi);
	acc.A.Add(Rel(MaxAbsDiff(Arec, Ai, NH), MaxAbs(Ai, NH)));
	acc.b.Add(Rel(MaxAbsDiff(brec, bi, N), std::max(MaxAbs(bi, N), kPrmH * MaxAbs(fe, N))));
}

// ---- L1-3: 격자 천의 전역 곱 ---------------------------------------------------
struct Stencil { std::vector<int> nodes; std::vector<double> ev; };

// 엔진 곱 연산의 흉내(SDESparseMatrix.cpp:4256-4308): 블록 B(행 노드 i, 열 노드 j)마다
// out_i += B (r_j − r_i),  out_j −= B' (r_j − r_i),  B' = 전치 처리 켬이면 Bᵀ, 끔이면 B
static void EngineProduct(const std::vector<int>& rows, const std::vector<int>& cols, const std::vector<double>& blocks,
	bool transposeLower, const std::vector<double>& r, std::vector<double>& out)
{
	std::fill(out.begin(), out.end(), 0.);
	for (std::size_t s = 0; s < rows.size(); ++s) {
		const int i = rows[s], j = cols[s];
		const double* B = &blocks[9 * s];
		const double d[3] = { r[3 * j] - r[3 * i], r[3 * j + 1] - r[3 * i + 1], r[3 * j + 2] - r[3 * i + 2] };
		for (int a = 0; a < 3; ++a) {
			const double pij = B[3 * a] * d[0] + B[3 * a + 1] * d[1] + B[3 * a + 2] * d[2];
			const double pji = transposeLower ? B[a] * d[0] + B[3 + a] * d[1] + B[6 + a] * d[2] : pij;
			out[3 * i + a] += pij;
			out[3 * j + a] -= pji;
		}
	}
}

template <class S, class E>
static void CheckGrid(const Variant& var, const std::vector<Stencil>& st, const std::vector<double>& X, const std::vector<double>& Vel,
	std::mt19937_64& rng, Acc& acc)
{
	enum { N = S::kDof, NN = S::kNodes, NC = S::kComp, NP = E::kPairs };
	const int nDof = static_cast<int>(X.size());
	std::unordered_map<long long, int> slot;
	std::vector<int> rows, cols;
	std::vector<double> bJ, bV;
	std::vector<std::vector<double> > refJ(st.size()), refV(st.size());
	for (std::size_t e = 0; e < st.size(); ++e) {
		double x[N], v[N];
		for (int a = 0; a < NN; ++a)
			for (int c = 0; c < 3; ++c) { x[3 * a + c] = X[3 * st[e].nodes[a] + c]; v[3 * a + c] = Vel[3 * st[e].nodes[a] + c]; }
		double C[NC], f[N], J[9 * NP], V[9 * NP];
		var.eng(x, &st[e].ev[0], v, kPrm2, C, f, J, V);
		for (int p = 0; p < NP; ++p) {
			const int ga = st[e].nodes[E::PairA(p)], gb = st[e].nodes[E::PairB(p)];
			const bool flip = ga > gb;
			const int lo = flip ? gb : ga, hi = flip ? ga : gb;
			const long long key = static_cast<long long>(lo) * (nDof / 3) + hi;
			std::unordered_map<long long, int>::iterator it = slot.find(key);
			int s;
			if (it == slot.end()) {
				s = static_cast<int>(rows.size());
				slot[key] = s;
				rows.push_back(lo); cols.push_back(hi);
				bJ.resize(bJ.size() + 9, 0.); bV.resize(bV.size() + 9, 0.);
			} else s = it->second;
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j) {
					const int k = flip ? 9 * p + 3 * j + i : 9 * p + 3 * i + j;   // (gb, ga) 슬롯에는 블록의 전치
					bJ[9 * s + 3 * i + j] += J[k];
					bV[9 * s + 3 * i + j] += V[k];
				}
		}
		double Cr[NC], Jr[N * N], Vr[N * N], fr[N];
		RefMatrices<S>(var.kind, x, &st[e].ev[0], v, Cr, Jr, Vr, fr);
		refJ[e].assign(Jr, Jr + N * N);
		refV[e].assign(Vr, Vr + N * N);
	}
	std::normal_distribution<double> nd(0., 1.);
	std::vector<double> r(nDof), outE(nDof), outD(nDof);
	for (int t = 0; t < kGridVectors; ++t) {
		for (int i = 0; i < nDof; ++i) r[i] = nd(rng);
		for (int m = 0; m < 2; ++m) {   // m = 0: J, 1: V
			const std::vector<double>& blocks = m == 0 ? bJ : bV;
			std::fill(outD.begin(), outD.end(), 0.);
			for (std::size_t e = 0; e < st.size(); ++e) {
				const std::vector<double>& M = m == 0 ? refJ[e] : refV[e];
				for (int a = 0; a < N; ++a) {
					double s = 0.;
					for (int b = 0; b < N; ++b) s += M[a * N + b] * r[3 * st[e].nodes[b / 3] + b % 3];
					outD[3 * st[e].nodes[a / 3] + a % 3] += s;
				}
			}
			const double scale = MaxAbs(&outD[0], nDof);
			EngineProduct(rows, cols, blocks, true, r, outE);
			(m == 0 ? acc.gJ : acc.gV).Add(Rel(MaxAbsDiff(&outE[0], &outD[0], nDof), scale));
			EngineProduct(rows, cols, blocks, false, r, outE);
			(m == 0 ? acc.gJoff : acc.gVoff).Add(Rel(MaxAbsDiff(&outE[0], &outD[0], nDof), scale));
		}
	}
}

struct SpecResult { std::string spec; Kind kind; Acc acc; };

// ---- L1-4: 감쇠 없는 GN ---------------------------------------------------------
typedef void (*FnU)(const double*, const double*, const double*, double*, double*, double*);
struct UndampedResult { std::string spec; int n, bad; };

// 같다 = 값이 같다(== 로 견줘 −0 과 +0 은 같게 본다; NaN 은 다르게 센다).
template <class E>
static UndampedResult CheckUndamped(const char* name, Fn8 gn, FnU un, const std::vector<Sample>& samples)
{
	UndampedResult r = { name, 0, 0 };
	const double ks[2] = { kPrmK, 3.7 };
	for (std::size_t s = 0; s < samples.size(); ++s)
		for (int q = 0; q < 2; ++q) {
			const std::vector<double> v0(samples[s].v.size(), 0.);
			const double prm2[2] = { ks[q], 0. }, prm1[1] = { ks[q] };
			double C1[E::kComp], f1[E::kDof], J1[9 * E::kPairs], V1[9 * E::kPairs];
			double C2[E::kComp], f2[E::kDof], J2[9 * E::kPairs];
			gn(samples[s].x.data(), samples[s].ev.data(), v0.data(), prm2, C1, f1, J1, V1);
			un(samples[s].x.data(), samples[s].ev.data(), prm1, C2, f2, J2);
			bool same = true;
			for (int i = 0; i < E::kComp; ++i) same = same && C1[i] == C2[i];
			for (int i = 0; i < E::kDof; ++i) same = same && f1[i] == f2[i];
			for (int i = 0; i < 9 * E::kPairs; ++i) same = same && J1[i] == J2[i];
			r.n++;
			if (!same) r.bad++;
		}
	return r;
}

// ---- L1-5b: 엔진 처방 T 판은 v = 0 에서 T 판과 같다 ------------------------------
template <class E>
static UndampedResult CheckRateAtRest(const char* name, Fn8 ten, Fn8 rate, const std::vector<Sample>& samples)
{
	UndampedResult r = { name, 0, 0 };
	for (std::size_t s = 0; s < samples.size(); ++s) {
		const std::vector<double> v0(samples[s].v.size(), 0.);
		double C1[E::kComp], f1[E::kDof], J1[9 * E::kPairs], V1[9 * E::kPairs];
		double C2[E::kComp], f2[E::kDof], J2[9 * E::kPairs], V2[9 * E::kPairs];
		ten(samples[s].x.data(), samples[s].ev.data(), v0.data(), kPrm2, C1, f1, J1, V1);
		rate(samples[s].x.data(), samples[s].ev.data(), v0.data(), kPrm2, C2, f2, J2, V2);
		bool same = true;
		for (int i = 0; i < E::kComp; ++i) same = same && C1[i] == C2[i];
		for (int i = 0; i < E::kDof; ++i) same = same && f1[i] == f2[i];
		for (int i = 0; i < 9 * E::kPairs; ++i) same = same && J1[i] == J2[i] && V1[i] == V2[i];
		r.n++;
		if (!same) r.bad++;
	}
	return r;
}

template <class S, class E>
static void RunSpec(const char* name, const std::vector<Variant>& vars, const std::vector<Sample>& samples,
	const std::vector<Stencil>& grid, const std::vector<double>& X, const std::vector<double>& Vel, std::mt19937_64& rng,
	std::vector<SpecResult>& out)
{
	for (std::size_t k = 0; k < vars.size(); ++k) {
		SpecResult r;
		r.spec = name;
		r.kind = vars[k].kind;
		for (std::size_t s = 0; s < samples.size(); ++s) CheckSample<S, E>(vars[k], samples[s], r.acc);
		CheckGrid<S, E>(vars[k], grid, X, Vel, rng, r.acc);
		out.push_back(r);
	}
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Hinge> hs = ReadHinges("results/hinges.txt");
	const std::vector<Tri> ts = ReadTriangles("results/triangles.txt");
	const std::vector<std::vector<double> > dg = ReadRows("results/provot_diag.txt");
	if (hs.empty() || ts.empty()) { std::printf("results/hinges.txt · triangles.txt 없음\n"); return 2; }
	if (dg.size() < ts.size()) { std::printf("results/provot_diag.txt 없음 — model/make_provot_diag.py 를 먼저 돌린다\n"); return 2; }
	Report rpt("results/gate_link_l1.json");
	std::FILE* js = rpt.json;
	std::printf("엔진 연결 L1 — 힌지 %zu 개, 삼각형 %zu 개, k = %g, kd = %g, h = %g, 격자 G = %d\n",
		hs.size(), ts.size(), kPrmK, kPrmKd, kPrmH, kGridG);

	// 표본(Gauss–Newton 게이트 N1 과 같은 스텐실 값과 C 척도)
	std::vector<Sample> sTheta, sBridson, sTri, sProvot;
	for (std::size_t s = 0; s < hs.size(); ++s) {
		Sample a;
		a.x.assign(hs[s].x, hs[s].x + 12); a.v.assign(hs[s].v, hs[s].v + 12);
		a.ev.assign(1, hs[s].theta0); a.cscale = 1.;
		sTheta.push_back(a);
		const double w = BendWeight(hs[s].x);
		a.ev.assign(1, w); a.cscale = std::sqrt(w);
		sBridson.push_back(a);
	}
	for (std::size_t s = 0; s < ts.size(); ++s) {
		Sample a;
		a.x.assign(ts[s].x, ts[s].x + 9); a.v.assign(ts[s].v, ts[s].v + 9);
		a.ev.assign(ts[s].uv, ts[s].uv + 6); a.cscale = std::sqrt(UvArea(ts[s].uv));
		sTri.push_back(a);
		if (dg[s].size() != 2) { std::printf("results/provot_diag.txt %zu 번째 줄이 두 칸이 아니다\n", s + 1); return 2; }
		a.ev.push_back(dg[s][0]); a.ev.push_back(dg[s][1]);   // Provot 전단: uv 6 칸 + 대각선 길이 2 칸
		sProvot.push_back(a);
	}

	// 격자 천(MakeGrid, 시간 측정용과 같은 만드는 법). θ0 = 0, Bridson w̄ = 그 상태의 값, 쉬는 uv = 평평한 격자(1 cm)
	std::mt19937_64 rng(20261113);
	const Grid grid = MakeGrid(kGridG, rng);
	std::vector<Stencil> gTheta, gBridson, gTri, gProvot;
	for (int s = 0; s < grid.nHinge; ++s) {
		Stencil a;
		a.nodes.assign(&grid.hinges[4 * s], &grid.hinges[4 * s] + 4);
		a.ev.assign(1, 0.);
		gTheta.push_back(a);
		double x[12];
		for (int k = 0; k < 4; ++k)
			for (int c = 0; c < 3; ++c) x[3 * k + c] = grid.X[3 * a.nodes[k] + c];
		a.ev.assign(1, BendWeight(x));
		gBridson.push_back(a);
	}
	for (int s = 0; s < grid.nTri; ++s) {
		Stencil a;
		a.nodes.assign(&grid.tris[3 * s], &grid.tris[3 * s] + 3);
		for (int k = 0; k < 3; ++k) { a.ev.push_back(0.01 * (a.nodes[k] % kGridG)); a.ev.push_back(0.01 * (a.nodes[k] / kGridG)); }
		gTri.push_back(a);
		a.ev.push_back(std::sqrt(2.0)); a.ev.push_back(std::sqrt(2.0));
		gProvot.push_back(a);
	}
	std::printf("격자: 노드 %d, 삼각형 %d, 힌지 %d\n", grid.nNode, grid.nTri, grid.nHinge);

	using namespace sxgen;
	std::vector<SpecResult> res;
	{
		std::vector<Variant> v;
		Variant e = { kExact, BaraffBendEngine::EvalEngine, BaraffBend::EvalImplicit }; v.push_back(e);
		Variant g = { kGN, BaraffBendEngine::EvalEngineGN, BaraffBend::EvalImplicitGN }; v.push_back(g);
		RunSpec<BaraffBend, BaraffBendEngine>("BaraffBend", v, sTheta, gTheta, grid.X, grid.V, rng, res);
	}
	{
		std::vector<Variant> v;
		Variant e = { kExact, BaraffBendCompactEngine::EvalEngine, BaraffBendCompact::EvalImplicit }; v.push_back(e);
		Variant g = { kGN, BaraffBendCompactEngine::EvalEngineGN, BaraffBendCompact::EvalImplicitGN }; v.push_back(g);
		RunSpec<BaraffBendCompact, BaraffBendCompactEngine>("BaraffBendCompact", v, sTheta, gTheta, grid.X, grid.V, rng, res);
	}
	{
		std::vector<Variant> v;
		Variant e = { kExact, BridsonBendEngine::EvalEngine, BridsonBend::EvalImplicit }; v.push_back(e);
		Variant g = { kGN, BridsonBendEngine::EvalEngineGN, BridsonBend::EvalImplicitGN }; v.push_back(g);
		RunSpec<BridsonBend, BridsonBendEngine>("BridsonBend", v, sBridson, gBridson, grid.X, grid.V, rng, res);
	}
	{
		std::vector<Variant> v;
		Variant e = { kExact, BaraffStretchEngine::EvalEngine, BaraffStretch::EvalImplicit }; v.push_back(e);
		Variant g = { kGN, BaraffStretchEngine::EvalEngineGN, BaraffStretch::EvalImplicitGN }; v.push_back(g);
		Variant t = { kTension, BaraffStretchEngine::EvalEngineTension, BaraffStretch::EvalImplicitTension }; v.push_back(t);
		Variant r = { kTensionRate, BaraffStretchEngine::EvalEngineTensionRate, 0 }; v.push_back(r);
		RunSpec<BaraffStretch, BaraffStretchEngine>("BaraffStretch", v, sTri, gTri, grid.X, grid.V, rng, res);
	}
	{
		std::vector<Variant> v;
		Variant e = { kExact, BaraffShearEngine::EvalEngine, BaraffShear::EvalImplicit }; v.push_back(e);
		Variant g = { kGN, BaraffShearEngine::EvalEngineGN, BaraffShear::EvalImplicitGN }; v.push_back(g);
		RunSpec<BaraffShear, BaraffShearEngine>("BaraffShear", v, sTri, gTri, grid.X, grid.V, rng, res);
	}
	{
		std::vector<Variant> v;
		Variant e = { kExact, ProvotShearEngine::EvalEngine, ProvotShear::EvalImplicit }; v.push_back(e);
		Variant g = { kGN, ProvotShearEngine::EvalEngineGN, ProvotShear::EvalImplicitGN }; v.push_back(g);
		Variant t = { kTension, ProvotShearEngine::EvalEngineTension, ProvotShear::EvalImplicitTension }; v.push_back(t);
		Variant r = { kTensionRate, ProvotShearEngine::EvalEngineTensionRate, 0 }; v.push_back(r);
		RunSpec<ProvotShear, ProvotShearEngine>("ProvotShear", v, sProvot, gProvot, grid.X, grid.V, rng, res);
	}

	std::fprintf(js, "{\n  \"hinges\": %zu, \"triangles\": %zu, \"k\": %g, \"kd\": %g, \"h\": %g, \"limit\": %g,\n"
		"  \"grid\": {\"G\": %d, \"nodes\": %d, \"triangles\": %d, \"hinges\": %d, \"vectors\": %d},\n  \"variants\": [\n",
		hs.size(), ts.size(), kPrmK, kPrmKd, kPrmH, kLimit, kGridG, grid.nNode, grid.nTri, grid.nHinge, kGridVectors);
	for (std::size_t k = 0; k < res.size(); ++k) {
		const SpecResult& r = res[k];
		const Acc& a = r.acc;
		std::printf("\n %s · %s (표본 %d 개)\n", r.spec.c_str(), KindName(r.kind), a.n);
		rpt.Gate("C", a.n > 0 && a.c.Max() <= kLimit, "L1-1 |ΔC|/척도 %.2e", a.c.Max(), kLimit);
		rpt.Gate("J", a.n > 0 && a.J.Max() <= kLimit, "L1-1 블록 rel max %.2e", a.J.Max(), kLimit);
		rpt.Gate("V", a.n > 0 && a.V.Max() <= kLimit, "L1-1 블록 rel max %.2e", a.V.Max(), kLimit);
		rpt.Gate("f+d", a.n > 0 && a.f.Max() <= kLimit, "L1-1 rel max %.2e", a.f.Max(), kLimit);
		if (r.kind != kTensionRate) {
			rpt.Gate("A", a.n > 0 && a.A.Max() <= kLimit, "L1-2 행 합 복원 rel max %.2e", a.A.Max(), kLimit);
			rpt.Gate("b", a.n > 0 && a.b.Max() <= kLimit, "L1-2 행 합 복원 rel max %.2e", a.b.Max(), kLimit);
		} else {   // L1-5: 짝이 없어 L1-2 는 없다. 대신 Ċ 의 두 갈래가 표본에 있어야 한다
			rpt.Gate("Ċ>0", a.nRatePos >= kMinBranch, "L1-5 갈래 %.0f 개", a.nRatePos, kMinBranch);
			rpt.Gate("Ċ≤0", a.nRateNeg >= kMinBranch, "L1-5 갈래 %.0f 개", a.nRateNeg, kMinBranch);
		}
		rpt.Gate("gJ", !a.gJ.v.empty() && a.gJ.Max() <= kLimit, "L1-3 전역 곱 J rel max %.2e", a.gJ.Max(), kLimit);
		rpt.Gate("gV", !a.gV.v.empty() && a.gV.Max() <= kLimit, "L1-3 전역 곱 V rel max %.2e", a.gV.Max(), kLimit);
		std::printf("  보고 — 전치 처리 끔: J rel max %.2e, V rel max %.2e;  성분 C > 0 %d 개 · C ≤ 0 %d 개\n",
			a.gJoff.Max(), a.gVoff.Max(), a.nPos, a.nNeg);
		std::fprintf(js, "    {\"spec\": \"%s\", \"kind\": \"%s\", \"samples\": %d, "
			"\"L1_1\": {\"C_abs_over_scale\": %.3e, \"J_rel\": %.3e, \"V_rel\": %.3e, \"f_rel\": %.3e, "
			"\"J_rel_median\": %.3e, \"V_rel_median\": %.3e, \"f_rel_median\": %.3e}, "
			"\"L1_2\": {\"A_rel\": %.3e, \"b_rel\": %.3e, \"A_rel_median\": %.3e, \"b_rel_median\": %.3e}, "
			"\"L1_3\": {\"J_rel\": %.3e, \"V_rel\": %.3e}, "
			"\"report\": {\"transpose_off_J_rel\": %.3e, \"transpose_off_V_rel\": %.3e, \"components_C_pos\": %d, \"components_C_nonpos\": %d, "
			"\"components_Cdot_pos\": %d, \"components_Cdot_nonpos\": %d}}%s\n",
			r.spec.c_str(), KindName(r.kind), a.n, a.c.Max(), a.J.Max(), a.V.Max(), a.f.Max(), a.J.Median(), a.V.Median(), a.f.Median(),
			a.A.Max(), a.b.Max(), a.A.Median(), a.b.Median(), a.gJ.Max(), a.gV.Max(), a.gJoff.Max(), a.gVoff.Max(), a.nPos, a.nNeg,
			a.nRatePos, a.nRateNeg,
			k + 1 < res.size() ? "," : "");
	}
	// L1-4
	std::vector<UndampedResult> und;
	und.push_back(CheckUndamped<BaraffBendEngine>("BaraffBend", BaraffBendEngine::EvalEngineGN, BaraffBendEngine::EvalEngineGNUndamped, sTheta));
	und.push_back(CheckUndamped<BaraffBendCompactEngine>("BaraffBendCompact", BaraffBendCompactEngine::EvalEngineGN, BaraffBendCompactEngine::EvalEngineGNUndamped, sTheta));
	und.push_back(CheckUndamped<BridsonBendEngine>("BridsonBend", BridsonBendEngine::EvalEngineGN, BridsonBendEngine::EvalEngineGNUndamped, sBridson));
	und.push_back(CheckUndamped<BaraffStretchEngine>("BaraffStretch", BaraffStretchEngine::EvalEngineGN, BaraffStretchEngine::EvalEngineGNUndamped, sTri));
	und.push_back(CheckUndamped<BaraffShearEngine>("BaraffShear", BaraffShearEngine::EvalEngineGN, BaraffShearEngine::EvalEngineGNUndamped, sTri));
	und.push_back(CheckUndamped<ProvotShearEngine>("ProvotShear", ProvotShearEngine::EvalEngineGN, ProvotShearEngine::EvalEngineGNUndamped, sProvot));
	std::printf("\n L1-4 감쇠 없는 GN = EvalEngineGN(v = 0, kd = 0) — C · f · J, k = %g · 3.7\n", kPrmK);
	std::fprintf(js, "  ],\n  \"L1_4\": [\n");
	for (std::size_t k = 0; k < und.size(); ++k) {
		const UndampedResult& u = und[k];
		const bool ok = u.n > 0 && u.bad == 0;
		std::printf("  %-18s %s  다른 표본 %d / %d  (기준 0)\n", u.spec.c_str(), ok ? "PASS" : "FAIL", u.bad, u.n);
		if (!ok) rpt.fail++;
		std::fprintf(js, "    {\"spec\": \"%s\", \"samples\": %d, \"different\": %d}%s\n", u.spec.c_str(), u.n, u.bad,
			k + 1 < und.size() ? "," : "");
	}
	// L1-5b
	std::vector<UndampedResult> rest;
	rest.push_back(CheckRateAtRest<BaraffStretchEngine>("BaraffStretch", BaraffStretchEngine::EvalEngineTension, BaraffStretchEngine::EvalEngineTensionRate, sTri));
	rest.push_back(CheckRateAtRest<ProvotShearEngine>("ProvotShear", ProvotShearEngine::EvalEngineTension, ProvotShearEngine::EvalEngineTensionRate, sProvot));
	std::printf("\n L1-5b 엔진 처방 T 판 = EvalEngineTension (v = 0) — C · f · J · V\n");
	std::fprintf(js, "  ],\n  \"L1_5b\": [\n");
	for (std::size_t k = 0; k < rest.size(); ++k) {
		const UndampedResult& u = rest[k];
		const bool ok = u.n > 0 && u.bad == 0;
		std::printf("  %-18s %s  다른 표본 %d / %d  (기준 0)\n", u.spec.c_str(), ok ? "PASS" : "FAIL", u.bad, u.n);
		if (!ok) rpt.fail++;
		std::fprintf(js, "    {\"spec\": \"%s\", \"samples\": %d, \"different\": %d}%s\n", u.spec.c_str(), u.n, u.bad,
			k + 1 < rest.size() ? "," : "");
	}
	std::fprintf(js, "  ],\n  \"gates_failed\": %d\n}\n", rpt.fail);
	std::printf("\n엔진 연결 L1 %s (실패 %d) → results/gate_link_l1.json\n", rpt.fail == 0 ? "PASS" : "FAIL", rpt.fail);
	return rpt.fail == 0 ? 0 : 1;
}
