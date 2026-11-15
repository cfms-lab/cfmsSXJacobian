//---------------------------------------------------------------------------
// measure_frozen_force.cpp — 길이 고정 근사의 가짜 힘을 주름진 천에서 잰다 (docs/MEASURE_FROZEN_FORCE.md, 판정 밖 측정).
//   M1 힘(정적)  ·  M2 변위, 굽힘만(양해법 이완)  ·  M3 변위, 인장 · 전단을 켠 정적 평형(z 고정, 면 안 자유도)
//   세 팔: 정확(BaraffBend::Eval / EvalStaged 의 ∇θ) · 길이 고정(EvalFrozen 의 ∇θ) · 보정(길이 고정 ∇θ 에서 변 방향 성분을 뺀 것)
//   → results/measure_frozen_force.json
// 천 · 삼각형 나눔 · 힌지는 tests/test_integration.cpp 의 MakeCloth 와 같다. 실행 디렉터리 = 리포 루트.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffBend.h"
#include "../generated/BaraffShear.h"
#include "../generated/BaraffStretch.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

const double kPi = 3.14159265358979323846;
const double kL = 0.5;            // 한 변 (m)
const double kRho = 0.2;          // 면밀도 (kg/m²)
const double kStretchK = 1000.;
const double kShearK = 100.;
const int kWaves = 2;             // 물결 수 → 파장 L / 2
const double kDtFactor = 0.005;   // M2 시간 간격 / τ
const double kCollapse = 0.2;     // M2: 힌지 변이 쉬는 길이의 이 배 아래면 무너짐
enum { kExact = 0, kFrozen = 1, kCorrected = 2 };
const char* const kArmName[3] = { "exact", "frozen", "corrected" };

struct Mesh {
	int G, nNode, nTri, nHinge;
	double ell;
	std::vector<double> X0, mass, uv;
	std::vector<int> tris, hinges;                   // 삼각형(반시계) 3 개 · 힌지 (h, i, j, k) 4 개
	std::unordered_map<long long, int> edgeHinge;    // 변 (작은 번호, 큰 번호) → 힌지
};

long long EdgeKey(const Mesh& m, int p, int q)
{
	if (p > q) std::swap(p, q);
	return static_cast<long long>(p) * m.nNode + q;
}

Mesh MakeMesh(int G)
{
	Mesh c;
	c.G = G;
	c.nNode = G * G;
	c.ell = kL / (G - 1);
	c.X0.assign(3 * static_cast<std::size_t>(c.nNode), 0.);
	for (int j = 0; j < G; ++j)
		for (int i = 0; i < G; ++i) { const int n = j * G + i; c.X0[3 * n] = c.ell * i; c.X0[3 * n + 1] = c.ell * j; }
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
		for (int e = 0; e < 3; ++e)
			third[static_cast<long long>(c.tris[3 * t + e]) * c.nNode + c.tris[3 * t + (e + 1) % 3]] = c.tris[3 * t + (e + 2) % 3];
	for (int t = 0; t < c.nTri; ++t)
		for (int e = 0; e < 3; ++e) {
			const int p = c.tris[3 * t + e], q = c.tris[3 * t + (e + 1) % 3];
			if (p > q) continue;
			std::unordered_map<long long, int>::const_iterator it = third.find(static_cast<long long>(q) * c.nNode + p);
			if (it == third.end()) continue;
			const int hh[4] = { c.tris[3 * t + (e + 2) % 3], p, q, it->second };
			c.edgeHinge[EdgeKey(c, p, q)] = static_cast<int>(c.hinges.size() / 4);
			c.hinges.insert(c.hinges.end(), hh, hh + 4);
		}
	c.nHinge = static_cast<int>(c.hinges.size() / 4);
	c.mass.assign(c.nNode, 0.);
	c.uv.resize(6 * static_cast<std::size_t>(c.nTri));
	for (int t = 0; t < c.nTri; ++t) {
		for (int a = 0; a < 3; ++a) { const int n = c.tris[3 * t + a]; c.uv[6 * t + 2 * a] = c.X0[3 * n]; c.uv[6 * t + 2 * a + 1] = c.X0[3 * n + 1]; }
		const double area = 0.5 * c.ell * c.ell;
		for (int a = 0; a < 3; ++a) c.mass[c.tris[3 * t + a]] += kRho * area / 3.;
	}
	return c;
}

std::vector<double> Corrugate(const Mesh& m, double A, double phi)
{
	std::vector<double> x(m.X0);
	const double cs = std::cos(phi), sn = std::sin(phi);
	for (int n = 0; n < m.nNode; ++n) x[3 * n + 2] = A * std::sin(2. * kPi * kWaves * (x[3 * n] * cs + x[3 * n + 1] * sn) / kL);
	return x;
}

inline void Gather(const std::vector<int>& idx, int s, int nn, const std::vector<double>& x, double* xl)
{
	for (int a = 0; a < nn; ++a) { const int n = idx[nn * s + a]; for (int t = 0; t < 3; ++t) xl[3 * a + t] = x[3 * n + t]; }
}

// 힌지의 ∇θ 를 팔마다. 돌려주는 값은 θ. staged = 정확 팔에서 단계판(EvalStaged — 같은 값, 빠르다)을 쓴다.
inline double HingeGrad(const double* x, int arm, bool staged, double* g)
{
	const double ev[1] = { 0. };
	double C = 0., H[sxgen::BaraffBend::kHess];
	if (arm == kExact) {
		if (staged) sxgen::BaraffBend::EvalStaged(x, ev, &C, g, H);
		else sxgen::BaraffBend::Eval(x, ev, &C, g, H);
		return C;
	}
	sxgen::BaraffBend::EvalFrozen(x, ev, &C, g, H);
	if (arm == kCorrected) {   // 변의 두 끝 i · j 에서 변 방향 성분을 뺀다
		double e[3] = { x[6] - x[3], x[7] - x[4], x[8] - x[5] };
		const double len = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
		for (int t = 0; t < 3; ++t) e[t] /= len;
		for (int a = 1; a <= 2; ++a) {
			const double d = g[3 * a] * e[0] + g[3 * a + 1] * e[1] + g[3 * a + 2] * e[2];
			for (int t = 0; t < 3; ++t) g[3 * a + t] -= d * e[t];
		}
	}
	return C;
}

void AddBend(const Mesh& m, const std::vector<double>& x, double k, int arm, bool staged, std::vector<double>& f)
{
	for (int s = 0; s < m.nHinge; ++s) {
		double xl[12], g[12];
		Gather(m.hinges, s, 4, x, xl);
		const double w = -k * HingeGrad(xl, arm, staged, g);
		for (int a = 0; a < 4; ++a) { const int n = m.hinges[4 * s + a]; for (int t = 0; t < 3; ++t) f[3 * n + t] += w * g[3 * a + t]; }
	}
}

inline int NodeClass(const Mesh& m, int n)   // 0 안쪽 · 1 y 가장자리(주름이 끝나는 곳, 모서리 포함) · 2 x 가장자리
{
	const int i = n % m.G, j = n / m.G;
	if (j == 0 || j == m.G - 1) return 1;
	if (i == 0 || i == m.G - 1) return 2;
	return 0;
}

inline double Norm3(const double* a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }

struct Stat {
	double mx, ss; int n;
	Stat() : mx(0.), ss(0.), n(0) {}
	void Add(double v) { mx = std::max(mx, std::fabs(v)); ss += v * v; ++n; }
	double Rms() const { return n ? std::sqrt(ss / n) : 0.; }
};

// ---- M1 힘 ---------------------------------------------------------------------------------
struct M1Out {
	double thetaMaxDeg, identityErr, correctionErr, stagedErr, fExactMax, fExactRms, edgeFormulaErr, tMax;
	Stat df[3];
};

M1Out RunM1(const Mesh& m, double A, double phi, double k)
{
	M1Out o;
	const std::vector<double> x = Corrugate(m, A, phi);
	double thetaMax = 0., gMax = 0., idErr = 0., corErr = 0., stErr = 0.;
	for (int s = 0; s < m.nHinge; ++s) {
		double xl[12], gE[12], gS[12], gF[12], gC[12];
		Gather(m.hinges, s, 4, x, xl);
		const double th = HingeGrad(xl, kExact, false, gE);
		HingeGrad(xl, kExact, true, gS);
		HingeGrad(xl, kFrozen, false, gF);
		HingeGrad(xl, kCorrected, false, gC);
		thetaMax = std::max(thetaMax, std::fabs(th));
		double e[3] = { xl[6] - xl[3], xl[7] - xl[4], xl[8] - xl[5] };
		const double len = Norm3(e);
		const double sc = std::sin(th) * std::cos(th) / len;
		for (int a = 0; a < 4; ++a)
			for (int t = 0; t < 3; ++t) {
				const double expect = (a == 2 ? sc : a == 1 ? -sc : 0.) * e[t] / len;
				idErr = std::max(idErr, std::fabs(gF[3 * a + t] - gE[3 * a + t] - expect));
				corErr = std::max(corErr, std::fabs(gC[3 * a + t] - gE[3 * a + t]));
				stErr = std::max(stErr, std::fabs(gS[3 * a + t] - gE[3 * a + t]));
				gMax = std::max(gMax, std::fabs(gE[3 * a + t]));
			}
	}
	o.thetaMaxDeg = thetaMax * 180. / kPi;
	o.identityErr = idErr / gMax;
	o.correctionErr = corErr / gMax;
	o.stagedErr = stErr / gMax;
	std::vector<double> fE(3 * static_cast<std::size_t>(m.nNode), 0.), fF(fE);
	AddBend(m, x, k, kExact, false, fE);
	AddBend(m, x, k, kFrozen, false, fF);
	Stat fe;
	for (int n = 0; n < m.nNode; ++n) {
		fe.Add(Norm3(&fE[3 * n]));
		const double d[3] = { fF[3 * n] - fE[3 * n], fF[3 * n + 1] - fE[3 * n + 1], fF[3 * n + 2] - fE[3 * n + 2] };
		o.df[NodeClass(m, n)].Add(Norm3(d));
	}
	o.fExactMax = fe.mx;
	o.fExactRms = fe.Rms();
	o.edgeFormulaErr = -1.;
	o.tMax = 0.;
	if (phi == 0.) {   // y 가장자리 노드의 가짜 힘 = T_e · (주름을 따라 안쪽)
		double err = 0.;
		for (int side = 0; side < 2; ++side)
			for (int i = 1; i + 1 < m.G; ++i) {
				const int n = (side ? (m.G - 1) * m.G : 0) + i, n2 = side ? n - m.G : n + m.G;
				const int s = m.edgeHinge.find(EdgeKey(m, n, n2))->second;
				double xl[12], g[12];
				Gather(m.hinges, s, 4, x, xl);
				const double th = HingeGrad(xl, kExact, false, g);
				double e[3] = { x[3 * n2] - x[3 * n], x[3 * n2 + 1] - x[3 * n + 1], x[3 * n2 + 2] - x[3 * n + 2] };
				const double len = Norm3(e);
				const double T = k * th * std::sin(th) * std::cos(th) / len;
				o.tMax = std::max(o.tMax, std::fabs(T));
				for (int t = 0; t < 3; ++t) err = std::max(err, std::fabs(fF[3 * n + t] - fE[3 * n + t] - T * e[t] / len));
			}
		o.edgeFormulaErr = o.tMax > 0. ? err / o.tMax : err;
	}
	return o;
}

// ---- M2 변위, 굽힘만 -------------------------------------------------------------------------
struct M2Arm {
	bool collapsed;
	double tStop, tEnd, minEdgeRatio, maxSpeedEnd;
	int steps;
	std::vector<double> x;
};

// maxSteps > 0 이면 그 걸음 수에서 멈춘다(무너진 길이 고정 팔과 같은 시각에서 견주려고)
M2Arm RunM2Arm(const Mesh& m, const std::vector<double>& x0, double k, int arm, double cFactor, int maxSteps)
{
	M2Arm o;
	const double q = 2. * kPi * kWaves / kL;
	const double omega = q * q * std::sqrt(k / kRho);
	const double c = cFactor * 2. * omega;
	const double dt = kDtFactor * m.ell * m.ell * std::sqrt(kRho / k);
	o.tEnd = 20. / omega;
	int nSteps = static_cast<int>(std::ceil(o.tEnd / dt));
	if (maxSteps > 0) nSteps = std::min(nSteps, maxSteps);
	o.x = x0;
	o.collapsed = false;
	o.minEdgeRatio = 1e300;
	o.tStop = o.tEnd;
	std::vector<double> v(o.x.size(), 0.), f(o.x.size());
	int step = 0;
	for (; step < nSteps; ++step) {
		std::fill(f.begin(), f.end(), 0.);
		AddBend(m, o.x, k, arm, true, f);
		for (int n = 0; n < m.nNode; ++n)
			for (int t = 0; t < 3; ++t) {
				const int a = 3 * n + t;
				v[a] = (v[a] + dt * f[a] / m.mass[n]) / (1. + c * dt);
				o.x[a] += dt * v[a];
			}
		if (step % 10 == 9 || step + 1 == nSteps) {
			double mn = 1e300;
			bool finite = true;
			for (int s = 0; s < m.nHinge; ++s) {
				const int i = m.hinges[4 * s + 1], j = m.hinges[4 * s + 2];
				const double e[3] = { o.x[3 * j] - o.x[3 * i], o.x[3 * j + 1] - o.x[3 * i + 1], o.x[3 * j + 2] - o.x[3 * i + 2] };
				const double e0[3] = { m.X0[3 * j] - m.X0[3 * i], m.X0[3 * j + 1] - m.X0[3 * i + 1], 0. };
				const double r = Norm3(e) / Norm3(e0);
				if (!(r == r) || r > 1e6) finite = false;
				mn = std::min(mn, r);
			}
			o.minEdgeRatio = std::min(o.minEdgeRatio, mn);
			if (!finite || mn < kCollapse) { o.collapsed = true; o.tStop = (step + 1) * dt; ++step; break; }
		}
	}
	o.steps = step;
	double vm = 0.;
	for (std::size_t a = 0; a < v.size(); ++a) vm = std::max(vm, std::fabs(v[a]));
	o.maxSpeedEnd = vm;
	return o;
}

// 가장자리에서 센 줄 d 마다 안쪽 방향 y 변위(기준 위치 ref 에 대해)의 평균 · 최대
void RowProfile(const Mesh& m, const std::vector<double>& x, const std::vector<double>& ref, std::vector<double>& mean, std::vector<double>& mx)
{
	const int half = (m.G - 1) / 2;
	mean.assign(half + 1, 0.);
	mx.assign(half + 1, 0.);
	std::vector<int> cnt(half + 1, 0);
	for (int j = 0; j < m.G; ++j) {
		const int d = std::min(j, m.G - 1 - j);
		const double sign = j < m.G - 1 - j ? 1. : j > m.G - 1 - j ? -1. : 1.;
		for (int i = 0; i < m.G; ++i) {
			const int n = j * m.G + i;
			const double u = sign * (x[3 * n + 1] - ref[3 * n + 1]);
			mean[d] += u;
			++cnt[d];
			if (std::fabs(u) > std::fabs(mx[d])) mx[d] = u;
		}
	}
	for (int d = 0; d <= half; ++d) mean[d] /= cnt[d];
}

// ---- M3 변위, 인장 · 전단을 켠 정적 평형 ---------------------------------------------------------
struct Band {   // 대칭 띠 행렬의 아래쪽: (i, j), i − bw ≤ j ≤ i
	int n, bw;
	std::vector<double> a;
	bool overflow;
	void Reset(int n_, int bw_) { n = n_; bw = bw_; overflow = false; a.assign(static_cast<std::size_t>(n) * (bw + 1), 0.); }
	double& At(int i, int j) { return a[static_cast<std::size_t>(i) * (bw + 1) + (j - i + bw)]; }
	void Add(int i, int j, double v) { if (i < j) return; if (i - j > bw) { overflow = true; return; } At(i, j) += v; }
};

bool Cholesky(Band& B)
{
	for (int i = 0; i < B.n; ++i) {
		const int j0 = std::max(0, i - B.bw);
		for (int j = j0; j <= i; ++j) {
			double s = B.At(i, j);
			for (int k = std::max(j0, j - B.bw); k < j; ++k) s -= B.At(i, k) * B.At(j, k);
			if (j < i) B.At(i, j) = s / B.At(j, j);
			else { if (!(s > 0.)) return false; B.At(i, i) = std::sqrt(s); }
		}
	}
	return true;
}

void Solve(Band& B, std::vector<double>& b)
{
	for (int i = 0; i < B.n; ++i) {
		double s = b[i];
		for (int k = std::max(0, i - B.bw); k < i; ++k) s -= B.At(i, k) * b[k];
		b[i] = s / B.At(i, i);
	}
	for (int i = B.n - 1; i >= 0; --i) {
		double s = b[i];
		for (int k = i + 1; k <= std::min(B.n - 1, i + B.bw); ++k) s -= B.At(k, i) * b[k];
		b[i] = s / B.At(i, i);
	}
}

template <int N>
inline int Upper(int a, int b) { if (a > b) std::swap(a, b); return a * N - a * (a - 1) / 2 + (b - a); }

// 면 안 힘 f(2 × 노드)와, K 가 있으면 준 Newton 행렬(인장 C > 0 판 + 전단 · 굽힘 Gauss–Newton)
void InPlane(const Mesh& m, const std::vector<double>& x, double kb, int arm, std::vector<double>& f, Band* K)
{
	std::fill(f.begin(), f.end(), 0.);
	for (int t = 0; t < m.nTri; ++t) {
		double xl[9], C[2], g[18], H[90];
		Gather(m.tris, t, 3, x, xl);
		sxgen::BaraffStretch::Eval(xl, &m.uv[6 * t], C, g, H);
		for (int c = 0; c < 2; ++c) {
			const double* gc = g + 9 * c;
			const double* Hc = H + 45 * c;
			for (int a = 0; a < 3; ++a)
				for (int p = 0; p < 2; ++p) {
					const int ra = 2 * m.tris[3 * t + a] + p;
					f[ra] += -kStretchK * C[c] * gc[3 * a + p];
					if (!K) continue;
					for (int b = 0; b < 3; ++b)
						for (int r = 0; r < 2; ++r)
							K->Add(ra, 2 * m.tris[3 * t + b] + r,
							       kStretchK * (gc[3 * a + p] * gc[3 * b + r] + std::max(C[c], 0.) * Hc[Upper<9>(3 * a + p, 3 * b + r)]));
				}
		}
		double Cs[1], gs[9], Hs[45];
		sxgen::BaraffShear::Eval(xl, &m.uv[6 * t], Cs, gs, Hs);
		for (int a = 0; a < 3; ++a)
			for (int p = 0; p < 2; ++p) {
				const int ra = 2 * m.tris[3 * t + a] + p;
				f[ra] += -kShearK * Cs[0] * gs[3 * a + p];
				if (!K) continue;
				for (int b = 0; b < 3; ++b)
					for (int r = 0; r < 2; ++r) K->Add(ra, 2 * m.tris[3 * t + b] + r, kShearK * gs[3 * a + p] * gs[3 * b + r]);
			}
	}
	for (int s = 0; s < m.nHinge; ++s) {
		double xl[12], g[12];
		Gather(m.hinges, s, 4, x, xl);
		const double C = HingeGrad(xl, arm, true, g);
		for (int a = 0; a < 4; ++a)
			for (int p = 0; p < 2; ++p) {
				const int ra = 2 * m.hinges[4 * s + a] + p;
				f[ra] += -kb * C * g[3 * a + p];
				if (!K) continue;
				for (int b = 0; b < 4; ++b)
					for (int r = 0; r < 2; ++r) K->Add(ra, 2 * m.hinges[4 * s + b] + r, kb * g[3 * a + p] * g[3 * b + r]);
			}
	}
	// 가운데 노드의 x, y 와 가운데 줄 오른쪽 끝 노드의 y 를 고정한다 — z 를 붙든 채 면 안으로 도는 것은 z 축 둘레의 강체 회전이라
	// 가운데 한 점만 붙들면 행렬이 특이하다(처음 실행에서 이 때문에 수렴하지 않았다).
	const int hj = m.G / 2;
	const int pins[3] = { 2 * (hj * m.G + hj), 2 * (hj * m.G + hj) + 1, 2 * (hj * m.G + m.G - 1) + 1 };
	for (int p = 0; p < 3; ++p) {
		const int d = pins[p];
		f[d] = 0.;
		if (!K) continue;
		for (int j = std::max(0, d - K->bw); j < d; ++j) K->At(d, j) = 0.;
		for (int i = d + 1; i <= std::min(K->n - 1, d + K->bw); ++i) K->At(i, d) = 0.;
		K->At(d, d) = 1.;
	}
}

inline double MaxAbs(const std::vector<double>& a) { double v = 0.; for (std::size_t i = 0; i < a.size(); ++i) v = std::max(v, std::fabs(a[i])); return v; }

struct M3Arm { bool converged; int iters; double residual; std::vector<double> x; };

M3Arm SolveM3(const Mesh& m, const std::vector<double>& x0, double kb, int arm)
{
	M3Arm o;
	o.x = x0;
	o.converged = false;
	const int n2 = 2 * m.nNode;
	const double tol = 1e-12 * kStretchK * m.ell;
	std::vector<double> f(n2), ft(n2), d;
	std::vector<double> xt;
	// 띠 너비는 스텐실에서 잰다 — 힌지 스텐실은 노드 번호가 2G + 1 까지 벌어진다(처음에 삼각형만 보고 2(G + 1) + 1 로 두어 굽힘 블록이
	// 띠 밖에 쓰였다. 행렬은 수렴에만 쓰여 수렴한 값은 같았지만 두 경우가 수렴하지 않았다).
	int span = 0;
	for (int t = 0; t < m.nTri; ++t)
		span = std::max(span, *std::max_element(&m.tris[3 * t], &m.tris[3 * t] + 3) - *std::min_element(&m.tris[3 * t], &m.tris[3 * t] + 3));
	for (int s = 0; s < m.nHinge; ++s)
		span = std::max(span, *std::max_element(&m.hinges[4 * s], &m.hinges[4 * s] + 4) - *std::min_element(&m.hinges[4 * s], &m.hinges[4 * s] + 4));
	Band K;
	int it = 0;
	for (; it < 2000; ++it) {
		K.Reset(n2, 2 * span + 1);
		InPlane(m, o.x, kb, arm, f, &K);
		if (K.overflow) break;
		const double res = MaxAbs(f);
		if (res <= tol) { o.converged = true; break; }
		if (!Cholesky(K)) break;
		d = f;
		Solve(K, d);
		double step = 1.;
		bool accepted = false;
		for (int h = 0; h < 30 && !accepted; ++h, step *= 0.5) {
			xt = o.x;
			for (int n = 0; n < m.nNode; ++n) { xt[3 * n] += step * d[2 * n]; xt[3 * n + 1] += step * d[2 * n + 1]; }
			InPlane(m, xt, kb, arm, ft, 0);
			if (MaxAbs(ft) < res) { o.x = xt; accepted = true; }
		}
		if (!accepted) break;
	}
	o.iters = it;
	InPlane(m, o.x, kb, arm, f, 0);
	o.residual = MaxAbs(f);
	return o;
}

// ---- 출력 -----------------------------------------------------------------------------------
void JsonArray(std::FILE* js, const char* key, const std::vector<double>& a, double scale)
{
	std::fprintf(js, "\"%s\": [", key);
	for (std::size_t i = 0; i < a.size(); ++i) std::fprintf(js, "%s%.9e", i ? ", " : "", a[i] * scale);
	std::fprintf(js, "]");
}

double MaxDiff(const std::vector<double>& a, const std::vector<double>& b)
{
	double v = 0.;
	for (std::size_t i = 0; i < a.size(); ++i) v = std::max(v, std::fabs(a[i] - b[i]));
	return v;
}

struct Clock {
	std::chrono::steady_clock::time_point t0;
	Clock() : t0(std::chrono::steady_clock::now()) {}
	double Sec() const { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }
};

} // namespace

int main()
{
	const Clock clock;
	const int kGs[3] = { 11, 21, 41 };
	std::vector<Mesh> meshes;
	for (int g = 0; g < 3; ++g) meshes.push_back(MakeMesh(kGs[g]));
	const double kAmps[4] = { 0.005, 0.01, 0.02, 0.04 };
	const double kPhi = 30. * kPi / 180.;
	std::FILE* js = std::fopen("results/measure_frozen_force.json", "w");
	if (!js) { std::printf("results/measure_frozen_force.json 을 열 수 없다\n"); return 1; }
	std::fprintf(js, "{\n \"L\": %g, \"rho\": %g, \"waves\": %d, \"stretch_k\": %g, \"shear_k\": %g, \"dt_over_tau\": %g, \"collapse_ratio\": %g,\n",
	             kL, kRho, kWaves, kStretchK, kShearK, kDtFactor, kCollapse);

	// ---- M1 ----
	struct C1 { int g; double A, phi, k; };
	const C1 c1[] = { { 1, 0.005, 0., 1e-3 }, { 1, 0.01, 0., 1e-3 }, { 1, 0.02, 0., 1e-3 }, { 1, 0.04, 0., 1e-3 }, { 0, 0.02, 0., 1e-3 },
	                  { 2, 0.02, 0., 1e-3 }, { 1, 0.02, 0., 1e-5 }, { 1, 0.02, kPhi, 1e-3 }, { 2, 0.02, kPhi, 1e-3 }, { 1, 0.04, kPhi, 1e-3 } };
	const int n1 = static_cast<int>(sizeof(c1) / sizeof(c1[0]));
	std::printf("M1 힘 (정적) — 가짜 힘 Δf = f_길이고정 − f_정확, 척도 = 정확한 굽힘 힘의 노드 최대\n");
	std::printf("   G   A(cm) 주름    k      θ최대   식 오차   보정 오차  |  Δf/척도: 안쪽 최대(rms) · y 가장자리 최대(rms) · x 가장자리 최대  |  가장자리 식 오차\n");
	std::fprintf(js, " \"M1\": [\n");
	for (int q = 0; q < n1; ++q) {
		const Mesh& m = meshes[c1[q].g];
		const M1Out o = RunM1(m, c1[q].A, c1[q].phi, c1[q].k);
		char ef[24];
		if (o.edgeFormulaErr < 0.) std::snprintf(ef, sizeof(ef), "-");
		else std::snprintf(ef, sizeof(ef), "%.1e", o.edgeFormulaErr);
		std::printf("  %3d  %4.1f  %-6s %.0e  %5.1f°  %.1e  %.1e  |  %.2e (%.2e) · %.2e (%.2e) · %.2e  |  %s\n",
		            m.G, c1[q].A * 100., c1[q].phi == 0. ? "나란" : "비스듬", c1[q].k, o.thetaMaxDeg, o.identityErr, o.correctionErr,
		            o.df[0].mx / o.fExactMax, o.df[0].Rms() / o.fExactMax, o.df[1].mx / o.fExactMax, o.df[1].Rms() / o.fExactMax,
		            o.df[2].mx / o.fExactMax, ef);
		std::fprintf(js, "  {\"G\": %d, \"A\": %g, \"phi_deg\": %g, \"k\": %g, \"theta_max_deg\": %.6f, \"hinge_identity_err\": %.3e, "
		             "\"correction_err\": %.3e, \"staged_vs_eval\": %.3e, \"f_exact_max\": %.9e, \"f_exact_rms\": %.9e, "
		             "\"df_interior_max\": %.9e, \"df_interior_rms\": %.9e, \"df_yedge_max\": %.9e, \"df_yedge_rms\": %.9e, "
		             "\"df_xedge_max\": %.9e, \"df_xedge_rms\": %.9e, \"yedge_formula_err\": %.3e, \"T_max\": %.9e}%s\n",
		             m.G, c1[q].A, c1[q].phi * 180. / kPi, c1[q].k, o.thetaMaxDeg, o.identityErr, o.correctionErr, o.stagedErr, o.fExactMax,
		             o.fExactRms, o.df[0].mx, o.df[0].Rms(), o.df[1].mx, o.df[1].Rms(), o.df[2].mx, o.df[2].Rms(), o.edgeFormulaErr, o.tMax,
		             q + 1 < n1 ? "," : "");
	}
	std::fprintf(js, " ],\n");

	// ---- M2 ----
	struct C2 { int g; double A, phi, k, cF; };
	const C2 c2[] = { { 1, 0.005, 0., 1e-3, 1. }, { 1, 0.01, 0., 1e-3, 1. }, { 1, 0.02, 0., 1e-3, 1. }, { 1, 0.04, 0., 1e-3, 1. },
	                  { 0, 0.02, 0., 1e-3, 1. }, { 2, 0.02, 0., 1e-3, 1. }, { 1, 0.02, 0., 1e-5, 1. }, { 1, 0.02, 0., 1e-3, 0.5 },
	                  { 1, 0.02, 0., 1e-3, 2. }, { 1, 0.02, kPhi, 1e-3, 1. } };
	const int n2 = static_cast<int>(sizeof(c2) / sizeof(c2[0]));
	std::printf("\nM2 변위, 굽힘만 — y 변위를 ℓ 로 나눈 값(안쪽이 +), 가장자리에서 센 줄 d = 0 · 1 · 2 … 가운데. 무너지면 그 시각에서 견준다\n");
	std::fprintf(js, " \"M2\": [\n");
	for (int q = 0; q < n2; ++q) {
		const Mesh& m = meshes[c2[q].g];
		const std::vector<double> x0 = Corrugate(m, c2[q].A, c2[q].phi);
		M2Arm arm[3];
		for (int a = 0; a < 3; ++a) arm[a] = RunM2Arm(m, x0, c2[q].k, a, c2[q].cF, 0);
		// 길이 고정 팔이 무너졌으면 정확 팔을 같은 걸음 수에서 멈춘 것과 견준다
		const M2Arm same = arm[kFrozen].collapsed ? RunM2Arm(m, x0, c2[q].k, kExact, c2[q].cF, arm[kFrozen].steps) : arm[kExact];
		std::vector<double> meanF, maxF, meanE, maxE, meanD, maxD;
		RowProfile(m, arm[kFrozen].x, m.X0, meanF, maxF);     // 길이 고정, 쉬는 위치 기준
		RowProfile(m, same.x, m.X0, meanE, maxE);             // 정확(같은 시각), 쉬는 위치 기준
		RowProfile(m, arm[kFrozen].x, same.x, meanD, maxD);   // 길이 고정 − 정확
		double exactDy = 0.;
		for (int n = 0; n < m.nNode; ++n) exactDy = std::max(exactDy, std::fabs(arm[kExact].x[3 * n + 1] - m.X0[3 * n + 1]));
		Stat dEdge, dIn;   // 길이 고정 팔과 정확 팔(같은 시각)의 위치 차 (벡터 크기)
		for (int n = 0; n < m.nNode; ++n) {
			const double d[3] = { arm[kFrozen].x[3 * n] - same.x[3 * n], arm[kFrozen].x[3 * n + 1] - same.x[3 * n + 1],
			                      arm[kFrozen].x[3 * n + 2] - same.x[3 * n + 2] };
			(NodeClass(m, n) ? dEdge : dIn).Add(Norm3(d));
		}
		const double corr = MaxDiff(arm[kCorrected].x, arm[kExact].x);
		const int half = (m.G - 1) / 2;
		std::printf("  G %2d  A %3.1f cm  %-6s k %.0e  감쇠 ×%.1f | 길이 고정: %s t/t끝 %.3f, 변 최소 %.2f | 길이 고정 − 정확, 줄 평균 %+.4f %+.4f %+.4f … %+.4f, "
		            "줄 0 최대 %+.4f | 줄 0 평균: 길이 고정 %+.4f · 정확 %+.4f | 정확 |Δy| 최대 %.1e · 보정 − 정확 %.1e | 팔 차이 가장자리 %.3f(rms %.3f) · "
		            "안쪽 %.3f(rms %.3f)\n",
		            m.G, c2[q].A * 100., c2[q].phi == 0. ? "나란" : "비스듬", c2[q].k, c2[q].cF, arm[kFrozen].collapsed ? "무너짐" : "끝까지",
		            arm[kFrozen].tStop / arm[kFrozen].tEnd, arm[kFrozen].minEdgeRatio, meanD[0] / m.ell, meanD[1] / m.ell, meanD[2] / m.ell,
		            meanD[half] / m.ell, maxD[0] / m.ell, meanF[0] / m.ell, meanE[0] / m.ell, exactDy / m.ell, corr / m.ell, dEdge.mx / m.ell,
		            dEdge.Rms() / m.ell, dIn.mx / m.ell, dIn.Rms() / m.ell);
		std::fprintf(js, "  {\"G\": %d, \"A\": %g, \"phi_deg\": %g, \"k\": %g, \"damping_factor\": %g, \"ell\": %.9e,\n", m.G, c2[q].A,
		             c2[q].phi * 180. / kPi, c2[q].k, c2[q].cF, m.ell);
		for (int a = 0; a < 3; ++a)
			std::fprintf(js, "   \"%s\": {\"collapsed\": %s, \"t_stop_over_t_end\": %.6f, \"steps\": %d, \"min_edge_ratio\": %.6f, \"max_speed_end\": %.3e},\n",
			             kArmName[a], arm[a].collapsed ? "true" : "false", arm[a].tStop / arm[a].tEnd, arm[a].steps, arm[a].minEdgeRatio,
			             arm[a].maxSpeedEnd);
		std::fprintf(js, "   \"exact_max_abs_dy_over_ell\": %.3e, \"corrected_minus_exact_max_over_ell\": %.3e,\n   ", exactDy / m.ell, corr / m.ell);
		JsonArray(js, "frozen_row_mean_inward_dy_over_ell", meanF, 1. / m.ell);
		std::fprintf(js, ",\n   ");
		JsonArray(js, "exact_same_time_row_mean_inward_dy_over_ell", meanE, 1. / m.ell);
		std::fprintf(js, ",\n   ");
		JsonArray(js, "frozen_minus_exact_row_mean_inward_dy_over_ell", meanD, 1. / m.ell);
		std::fprintf(js, ",\n   ");
		JsonArray(js, "frozen_minus_exact_row_max_inward_dy_over_ell", maxD, 1. / m.ell);
		std::fprintf(js, ",\n   \"frozen_minus_exact_over_ell\": {\"edge_max\": %.6e, \"edge_rms\": %.6e, \"interior_max\": %.6e, \"interior_rms\": %.6e}}%s\n",
		             dEdge.mx / m.ell, dEdge.Rms() / m.ell, dIn.mx / m.ell, dIn.Rms() / m.ell, q + 1 < n2 ? "," : "");
	}
	std::fprintf(js, " ],\n");

	// ---- M3 ----
	struct C3 { int g; double A, k; };
	const C3 c3[] = { { 1, 0.005, 1e-3 }, { 1, 0.01, 1e-3 }, { 1, 0.02, 1e-3 }, { 1, 0.04, 1e-3 }, { 0, 0.02, 1e-3 }, { 2, 0.02, 1e-3 }, { 1, 0.02, 1e-5 } };
	const int n3 = static_cast<int>(sizeof(c3) / sizeof(c3[0]));
	std::printf("\nM3 변위, 인장 · 전단 켬 (z 고정 정적 평형) — 길이 고정 팔 − 정확 팔의 y 변위, 가장자리의 안쪽 변위(사다리꼴 평균)\n");
	std::fprintf(js, " \"M3\": [\n");
	for (int q = 0; q < n3; ++q) {
		const Mesh& m = meshes[c3[q].g];
		const std::vector<double> x0 = Corrugate(m, c3[q].A, 0.);
		M3Arm arm[3];
		for (int a = 0; a < 3; ++a) arm[a] = SolveM3(m, x0, c3[q].k, a);
		double exactDy = 0.;
		for (int n = 0; n < m.nNode; ++n) exactDy = std::max(exactDy, std::fabs(arm[kExact].x[3 * n + 1] - m.X0[3 * n + 1]));
		const double corr = MaxDiff(arm[kCorrected].x, arm[kExact].x);
		const int G = m.G, half = (G - 1) / 2;
		// 가운데 줄 기준, 줄마다 안쪽 변위의 사다리꼴 평균 · 최대
		std::vector<double> prof(half + 1, 0.), profMax(half + 1, 0.);
		for (int d = 0; d <= half; ++d)
			for (int side = 0; side < 2; ++side) {
				const int j = side ? G - 1 - d : d;
				const double sign = side ? -1. : 1.;
				for (int i = 0; i < G; ++i) {
					const int n = j * G + i, nc = half * G + i;
					const double u = sign * ((arm[kFrozen].x[3 * n + 1] - arm[kExact].x[3 * n + 1]) - (arm[kFrozen].x[3 * nc + 1] - arm[kExact].x[3 * nc + 1]));
					prof[d] += 0.5 * u * ((i == 0 || i == G - 1) ? 1. : 2.) / (2. * (G - 1));
					profMax[d] = std::max(profMax[d], std::fabs(u));
				}
			}
		// 힘의 평형에서: Σ_줄 Σ_i T_i / (k_인장 (G − 1)), 두 반쪽의 평균
		double sumT = 0., thetaMax = 0.;
		for (int j = 0; j + 1 < G; ++j)
			for (int i = 1; i + 1 < G; ++i) {
				const int n = j * G + i, nb = n + G;
				const int s = m.edgeHinge.find(EdgeKey(m, n, nb))->second;
				double xl[12], g[12];
				Gather(m.hinges, s, 4, arm[kFrozen].x, xl);
				const double th = HingeGrad(xl, kExact, true, g);
				const double e[3] = { xl[6] - xl[3], xl[7] - xl[4], xl[8] - xl[5] };
				sumT += c3[q].k * th * std::sin(th) * std::cos(th) / Norm3(e);
				thetaMax = std::max(thetaMax, std::fabs(th));
			}
		const double predicted = 0.5 * sumT / (kStretchK * (G - 1));
		std::printf("  G %2d  A %3.1f cm  k %.0e | θ최대 %5.1f° | 수렴 %d%d%d (반복 %d · %d · %d, 잔차 %.0e · %.0e · %.0e) | 가장자리 평균 %.4e m (최대 %.4e) = "
		            "변형률 %.3e | 힘의 평형 식 %.4e m, 비 %.4f | 정확 |Δy| %.1e · 보정 − 정확 %.1e m\n",
		            G, c3[q].A * 100., c3[q].k, thetaMax * 180. / kPi, arm[0].converged, arm[1].converged, arm[2].converged, arm[0].iters, arm[1].iters,
		            arm[2].iters, arm[0].residual, arm[1].residual, arm[2].residual, prof[0], profMax[0], prof[0] / (0.5 * kL), predicted,
		            predicted != 0. ? prof[0] / predicted : 0., exactDy, corr);
		std::fprintf(js, "  {\"G\": %d, \"A\": %g, \"k\": %g, \"ell\": %.9e, \"theta_max_deg\": %.6f,\n", G, c3[q].A, c3[q].k, m.ell, thetaMax * 180. / kPi);
		for (int a = 0; a < 3; ++a)
			std::fprintf(js, "   \"%s\": {\"converged\": %s, \"iterations\": %d, \"residual\": %.3e},\n", kArmName[a], arm[a].converged ? "true" : "false",
			             arm[a].iters, arm[a].residual);
		std::fprintf(js, "   \"exact_max_abs_dy\": %.3e, \"corrected_minus_exact_max\": %.3e, \"edge_inward_mean\": %.9e, \"edge_inward_max\": %.9e, "
		             "\"mean_strain\": %.9e, \"force_balance_prediction\": %.9e, \"measured_over_prediction\": %.6f,\n   ",
		             exactDy, corr, prof[0], profMax[0], prof[0] / (0.5 * kL), predicted, predicted != 0. ? prof[0] / predicted : 0.);
		JsonArray(js, "row_inward_mean", prof, 1.);
		std::fprintf(js, "}%s\n", q + 1 < n3 ? "," : "");
	}
	std::fprintf(js, " ],\n \"seconds\": %.1f\n}\n", clock.Sec());
	std::fclose(js);
	std::printf("\n→ results/measure_frozen_force.json (%.0f s)\n", clock.Sec());
	return 0;
}
