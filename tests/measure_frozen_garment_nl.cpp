//---------------------------------------------------------------------------
// measure_frozen_garment_nl.cpp — 옷깃 패널의 비선형 평형: 정확한 굽힘 힘과 길이 고정 굽힘 힘으로 따로 평형시켜 그 차를 잰다
// (docs/MEASURE_FROZEN_FORCE.md 「M4c」, 판정 밖 측정). tests/measure_frozen_garment_nl.py 가 패널을 뽑아 부른다.
//   입력(글): measure_frozen_garment.cpp 와 같다 — E SE / nNode nTri / 노드마다 x y z u v 고정여부 / 삼각형마다 a b c
//   힘: 인장(Baraff–Witkin, k = E) + Provot 전단(k = SE) + 굽힘 C = θ(θ₀ = 0, 힌지 강성 k = BE · w̄ / 2) — measure_frozen_garment.cpp 와 같다
//   A. 옷깃만 다시 드레이프: 봉제 노드 고정 + 중력(−y). 정확한 힘의 평형 x_E → 거기서 길이 고정 힘의 평형 x_F, u_NL = x_F − x_E.
//      x_E 에서의 선형 응답 둘: u_L = K̃⁻¹Δf(준정부호 행렬), u_T = K_T⁻¹Δf(참 접선).
//   B. CLO 모양 x₀ 에서의 죽은 하중 문제: 하중 = −(정확한 모형의 내력)(x₀). 길이 고정 팔이 평형에 닿는지 본다.
//   풀이: Newton 꼴 반복 + 질량 비례 감쇠(Levenberg–Marquardt). 접선은 참 접선, 양의 정부호가 아니면 준정부호 행렬.
//         선형 풀이는 노드를 패턴의 긴 쪽 좌표로 줄 세운 띠 행렬 Cholesky.
//   사용:  measure_frozen_garment_nl <입력> <출력 json> <면밀도 g/cm²> <중력 cm/s²> <BE> [<BE> ...]
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffBend.h"
#include "../generated/BaraffStretch.h"
#include "../generated/ProvotShear.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace {

const double kPi = 3.14159265358979323846;
const double kDiag = static_cast<double>(1.41421356237f);   // 엔진 Provot 전단의 b (float 리터럴)
const int kExact = 0, kFrozen = 1, kCorrected = 2;

struct Panel {
	int nNode, nTri, nHinge, hb;
	double E, SE;
	std::vector<double> x0, uv, mass, wbar;
	std::vector<char> pinned;
	std::vector<int> tris, hinges, rank;
};

inline double Norm3(const double* a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }

bool Read(const char* path, Panel& p)
{
	std::FILE* f = std::fopen(path, "r");
	if (!f) return false;
	bool ok = std::fscanf(f, "%lf %lf %d %d", &p.E, &p.SE, &p.nNode, &p.nTri) == 4;
	p.x0.resize(3 * static_cast<std::size_t>(p.nNode));
	p.uv.resize(2 * static_cast<std::size_t>(p.nNode));
	p.pinned.assign(p.nNode, 0);
	for (int n = 0; ok && n < p.nNode; ++n) {
		int pin = 0;
		ok = std::fscanf(f, "%lf %lf %lf %lf %lf %d", &p.x0[3 * n], &p.x0[3 * n + 1], &p.x0[3 * n + 2], &p.uv[2 * n], &p.uv[2 * n + 1], &pin) == 6;
		p.pinned[n] = static_cast<char>(pin);
	}
	p.tris.resize(3 * static_cast<std::size_t>(p.nTri));
	for (int t = 0; ok && t < p.nTri; ++t) ok = std::fscanf(f, "%d %d %d", &p.tris[3 * t], &p.tris[3 * t + 1], &p.tris[3 * t + 2]) == 3;
	std::fclose(f);
	// 좌표를 패널의 무게중심 기준으로 옮긴다(병진 — 힘 · 평형은 같다). 옷 좌표는 160 cm 근처라 그대로 두면 위치의 반올림(3e-14 cm)에
	// 인장 강성(1e6 dyn/cm)을 곱한 3e-7 dyn 이 잔차의 바닥이 된다.
	double c[3] = { 0., 0., 0. };
	for (int n = 0; ok && n < p.nNode; ++n) for (int t = 0; t < 3; ++t) c[t] += p.x0[3 * n + t] / p.nNode;
	for (int n = 0; ok && n < p.nNode; ++n) for (int t = 0; t < 3; ++t) p.x0[3 * n + t] -= c[t];
	return ok;
}

// 힌지 (h, i, j, k) · 쉬는 가중치 w̄ · 노드 질량 · 띠 행렬의 줄 세우기 — 힌지와 w̄ 는 measure_frozen_garment.cpp 와 같은 짜임
void Build(Panel& p, double density)
{
	std::unordered_map<long long, int> third;
	for (int t = 0; t < p.nTri; ++t)
		for (int e = 0; e < 3; ++e)
			third[static_cast<long long>(p.tris[3 * t + e]) * p.nNode + p.tris[3 * t + (e + 1) % 3]] = p.tris[3 * t + (e + 2) % 3];
	for (int t = 0; t < p.nTri; ++t)
		for (int e = 0; e < 3; ++e) {
			const int a = p.tris[3 * t + e], b = p.tris[3 * t + (e + 1) % 3];
			std::unordered_map<long long, int>::const_iterator it = third.find(static_cast<long long>(b) * p.nNode + a);
			if (it == third.end() || a > b) continue;
			const int hh[4] = { p.tris[3 * t + (e + 2) % 3], a, b, it->second };
			p.hinges.insert(p.hinges.end(), hh, hh + 4);
		}
	p.nHinge = static_cast<int>(p.hinges.size() / 4);
	p.wbar.resize(p.nHinge);
	for (int s = 0; s < p.nHinge; ++s) {   // 쉬는 uv 에서 w̄ = 3|e|² / (A₁ + A₂)
		const int* h = &p.hinges[4 * s];
		const double* u0 = &p.uv[2 * h[0]]; const double* u1 = &p.uv[2 * h[1]]; const double* u2 = &p.uv[2 * h[2]]; const double* u3 = &p.uv[2 * h[3]];
		const double a1 = 0.5 * std::fabs((u1[0] - u0[0]) * (u2[1] - u0[1]) - (u2[0] - u0[0]) * (u1[1] - u0[1]));
		const double a2 = 0.5 * std::fabs((u2[0] - u3[0]) * (u1[1] - u3[1]) - (u1[0] - u3[0]) * (u2[1] - u3[1]));
		const double e2 = (u2[0] - u1[0]) * (u2[0] - u1[0]) + (u2[1] - u1[1]) * (u2[1] - u1[1]);
		p.wbar[s] = 3. * e2 / (a1 + a2);
	}
	p.mass.assign(p.nNode, 0.);
	for (int t = 0; t < p.nTri; ++t) {   // 삼각형의 쉬는 넓이를 세 꼭짓점에 똑같이
		const double* u0 = &p.uv[2 * p.tris[3 * t]]; const double* u1 = &p.uv[2 * p.tris[3 * t + 1]]; const double* u2 = &p.uv[2 * p.tris[3 * t + 2]];
		const double area = 0.5 * std::fabs((u1[0] - u0[0]) * (u2[1] - u0[1]) - (u2[0] - u0[0]) * (u1[1] - u0[1]));
		for (int a = 0; a < 3; ++a) p.mass[p.tris[3 * t + a]] += density * area / 3.;
	}
	// 줄 세우기: 패턴의 긴 쪽 좌표 순서 → 요소가 닿는 자유도의 차(반띠폭)가 작다
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
	for (int n = 0; n < p.nNode; ++n)
		for (int c = 0; c < 2; ++c) { lo[c] = std::min(lo[c], p.uv[2 * n + c]); hi[c] = std::max(hi[c], p.uv[2 * n + c]); }
	const int axis = (hi[0] - lo[0]) >= (hi[1] - lo[1]) ? 0 : 1;
	std::vector<int> order(p.nNode);
	std::iota(order.begin(), order.end(), 0);
	const std::vector<double>& uv = p.uv;
	std::stable_sort(order.begin(), order.end(), [&uv, axis](int a, int b) { return uv[2 * a + axis] < uv[2 * b + axis]; });
	p.rank.resize(p.nNode);
	for (int i = 0; i < p.nNode; ++i) p.rank[order[i]] = i;
	int span = 0;
	for (int t = 0; t < p.nTri; ++t) {
		const int r[3] = { p.rank[p.tris[3 * t]], p.rank[p.tris[3 * t + 1]], p.rank[p.tris[3 * t + 2]] };
		span = std::max(span, *std::max_element(r, r + 3) - *std::min_element(r, r + 3));
	}
	for (int s = 0; s < p.nHinge; ++s) {
		const int r[4] = { p.rank[p.hinges[4 * s]], p.rank[p.hinges[4 * s + 1]], p.rank[p.hinges[4 * s + 2]], p.rank[p.hinges[4 * s + 3]] };
		span = std::max(span, *std::max_element(r, r + 4) - *std::min_element(r, r + 4));
	}
	p.hb = 3 * span + 2;
}

// 대칭 띠 행렬의 아래 삼각(반띠폭 hb): A(i, j), i − hb ≤ j ≤ i
struct Band {
	int n, hb;
	std::vector<double> a;
	void Init(int n_, int hb_) { n = n_; hb = hb_; a.assign(static_cast<std::size_t>(n) * (hb + 1), 0.); }
	void Zero() { std::fill(a.begin(), a.end(), 0.); }
	inline double& At(int i, int j) { return a[static_cast<std::size_t>(i) * (hb + 1) + (j - i + hb)]; }
	inline double Get(int i, int j) const { return a[static_cast<std::size_t>(i) * (hb + 1) + (j - i + hb)]; }
	bool Cholesky()
	{
		for (int i = 0; i < n; ++i) {
			const int j0 = std::max(0, i - hb);
			for (int j = j0; j <= i; ++j) {
				double s = Get(i, j);
				const int k0 = std::max(j0, j - hb);
				for (int k = k0; k < j; ++k) s -= Get(i, k) * Get(j, k);
				if (j < i) At(i, j) = s / Get(j, j);
				else { if (!(s > 0.)) return false; At(i, i) = std::sqrt(s); }
			}
		}
		return true;
	}
	void Solve(std::vector<double>& b) const
	{
		for (int i = 0; i < n; ++i) {
			double s = b[i];
			for (int k = std::max(0, i - hb); k < i; ++k) s -= Get(i, k) * b[k];
			b[i] = s / Get(i, i);
		}
		for (int i = n - 1; i >= 0; --i) {
			double s = b[i];
			const int k1 = std::min(n - 1, i + hb);
			for (int k = i + 1; k <= k1; ++k) s -= Get(k, i) * b[k];
			b[i] = s / Get(i, i);
		}
	}
};

inline void Gather(const std::vector<int>& idx, int s, int nn, const std::vector<double>& x, double* xl)
{
	for (int a = 0; a < nn; ++a) { const int n = idx[nn * s + a]; for (int t = 0; t < 3; ++t) xl[3 * a + t] = x[3 * n + t]; }
}

template <int N>
inline int Upper(int a, int b) { if (a > b) std::swap(a, b); return a * N - a * (a - 1) / 2 + (b - a); }

inline void TriValues(const Panel& p, int t, double* ev)
{
	for (int a = 0; a < 3; ++a) { ev[2 * a] = p.uv[2 * p.tris[3 * t + a]]; ev[2 * a + 1] = p.uv[2 * p.tris[3 * t + a] + 1]; }
	ev[6] = kDiag;
	ev[7] = kDiag;
}

// 힌지의 굽힘 각 θ · 정확한 기울기 gE · 변형 기울기 gV(variant) · ∇²θ
inline void HingeEval(const double* xl, int variant, double& th, double* gE, double* gV, double* H)
{
	const double ev[1] = { 0. };
	sxgen::BaraffBend::Eval(xl, ev, &th, gE, H);
	if (variant == kExact) { for (int a = 0; a < 12; ++a) gV[a] = gE[a]; return; }
	double c2 = 0., H2[sxgen::BaraffBend::kHess];
	sxgen::BaraffBend::EvalFrozen(xl, ev, &c2, gV, H2);
	if (variant == kCorrected) {   // 보정: 변의 두 끝에서 변 방향 성분을 뺀다
		const double e[3] = { xl[6] - xl[3], xl[7] - xl[4], xl[8] - xl[5] };
		const double len = Norm3(e);
		for (int a = 1; a <= 2; ++a) {
			const double d = (gV[3 * a] * e[0] + gV[3 * a + 1] * e[1] + gV[3 * a + 2] * e[2]) / len;
			for (int t = 0; t < 3; ++t) gV[3 * a + t] -= d * e[t] / len;
		}
	}
}

// 내력(인장 + 전단 + 굽힘[variant]) + 하중. 고정 노드의 힘은 0. energy = 정확한 모형의 위치 에너지 − 하중 · x.
// scale = 자유 노드에서 서로 지우는 요소 힘 · 하중의 크기 합의 최대(잔차의 반올림 바닥을 재는 자).
void Forces(const Panel& p, const std::vector<double>& x, double BE, int variant, const std::vector<double>& load, std::vector<double>& f, double* energy,
            double* scale = 0)
{
	const int n3 = 3 * p.nNode;
	f.assign(n3, 0.);
	std::vector<double> fa(scale ? n3 : 0, 0.);
	double en = 0.;
	for (int t = 0; t < p.nTri; ++t) {
		double xl[9], ev[8], C[2], g[18], H[90];
		Gather(p.tris, t, 3, x, xl);
		TriValues(p, t, ev);
		for (int kind = 0; kind < 2; ++kind) {
			if (kind == 0) sxgen::BaraffStretch::Eval(xl, ev, C, g, H);
			else sxgen::ProvotShear::Eval(xl, ev, C, g, H);
			const double k = kind == 0 ? p.E : p.SE;
			for (int c = 0; c < 2; ++c) {
				en += 0.5 * k * C[c] * C[c];
				for (int a = 0; a < 9; ++a) {
					const int d = 3 * p.tris[3 * t + a / 3] + a % 3;
					f[d] -= k * C[c] * g[9 * c + a];
					if (scale) fa[d] += std::fabs(k * C[c] * g[9 * c + a]);
				}
			}
		}
	}
	for (int s = 0; s < p.nHinge; ++s) {
		double xl[12], th = 0., gE[12], gV[12], H[sxgen::BaraffBend::kHess];
		Gather(p.hinges, s, 4, x, xl);
		HingeEval(xl, variant, th, gE, gV, H);
		const double k = BE * p.wbar[s] / 2.;
		en += 0.5 * k * th * th;
		for (int a = 0; a < 12; ++a) {
			const int d = 3 * p.hinges[4 * s + a / 3] + a % 3;
			f[d] -= k * th * gV[a];
			if (scale) fa[d] += std::fabs(k * th * gV[a]);
		}
	}
	for (int i = 0; i < n3; ++i) {
		f[i] += load[i];
		en -= load[i] * x[i];
		if (scale) fa[i] += std::fabs(load[i]);
	}
	for (int n = 0; n < p.nNode; ++n)
		if (p.pinned[n]) for (int t = 0; t < 3; ++t) f[3 * n + t] = 0.;
	if (energy) *energy = en;
	if (scale) {
		*scale = 0.;
		for (int n = 0; n < p.nNode; ++n)
			if (!p.pinned[n]) *scale = std::max(*scale, Norm3(&fa[3 * n]));
	}
}

// 접선 −∂f/∂x 의 띠 행렬. kind 0 = 정확한 모형의 참 접선 Σ k(∇C∇Cᵀ + C∇²C), kind 1 = 준정부호(막은 C > 0 판, 굽힘은 Gauss–Newton).
// shift 는 질량에 비례해 대각에 더한다(Levenberg–Marquardt). 고정 자유도는 단위 행.
void Tangent(const Panel& p, const std::vector<double>& x, double BE, int kind, double shift, Band& K)
{
	K.Zero();
	for (int t = 0; t < p.nTri; ++t) {
		double xl[9], ev[8], C[2], g[18], H[90];
		Gather(p.tris, t, 3, x, xl);
		TriValues(p, t, ev);
		for (int force = 0; force < 2; ++force) {
			if (force == 0) sxgen::BaraffStretch::Eval(xl, ev, C, g, H);
			else sxgen::ProvotShear::Eval(xl, ev, C, g, H);
			const double k = force == 0 ? p.E : p.SE;
			for (int c = 0; c < 2; ++c) {
				const double ch = kind == 0 ? C[c] : std::max(C[c], 0.);
				for (int a = 0; a < 9; ++a) {
					const int ra = 3 * p.rank[p.tris[3 * t + a / 3]] + a % 3;
					for (int b = 0; b < 9; ++b) {
						const int rb = 3 * p.rank[p.tris[3 * t + b / 3]] + b % 3;
						if (ra >= rb) K.At(ra, rb) += k * (g[9 * c + a] * g[9 * c + b] + ch * H[45 * c + Upper<9>(a, b)]);
					}
				}
			}
		}
	}
	for (int s = 0; s < p.nHinge; ++s) {
		double xl[12], th = 0., gE[12], gV[12], H[sxgen::BaraffBend::kHess];
		Gather(p.hinges, s, 4, x, xl);
		HingeEval(xl, kExact, th, gE, gV, H);
		const double k = BE * p.wbar[s] / 2.;
		const double ch = kind == 0 ? th : 0.;
		for (int a = 0; a < 12; ++a) {
			const int ra = 3 * p.rank[p.hinges[4 * s + a / 3]] + a % 3;
			for (int b = 0; b < 12; ++b) {
				const int rb = 3 * p.rank[p.hinges[4 * s + b / 3]] + b % 3;
				if (ra >= rb) K.At(ra, rb) += k * (gE[a] * gE[b] + ch * H[Upper<12>(a, b)]);
			}
		}
	}
	for (int n = 0; n < p.nNode; ++n)
		for (int t = 0; t < 3; ++t) {
			const int d = 3 * p.rank[n] + t;
			if (p.pinned[n]) {
				for (int j = std::max(0, d - K.hb); j < d; ++j) K.At(d, j) = 0.;
				for (int i = d + 1; i <= std::min(K.n - 1, d + K.hb); ++i) K.At(i, d) = 0.;
				K.At(d, d) = 1.;
			}
			else K.At(d, d) += shift * p.mass[n];
		}
}

// b(노드 순서) → K⁻¹ b(노드 순서). K 는 Cholesky 가 끝난 띠 행렬.
void SolveNodes(const Panel& p, const Band& K, const std::vector<double>& b, std::vector<double>& out)
{
	std::vector<double> bb(b.size());
	for (int n = 0; n < p.nNode; ++n) for (int t = 0; t < 3; ++t) bb[3 * p.rank[n] + t] = b[3 * n + t];
	K.Solve(bb);
	out.resize(b.size());
	for (int n = 0; n < p.nNode; ++n) for (int t = 0; t < 3; ++t) out[3 * n + t] = bb[3 * p.rank[n] + t];
}

double MaxNodal(const Panel& p, const std::vector<double>& v, int* node)
{
	double mx = 0.;
	if (node) *node = -1;
	for (int n = 0; n < p.nNode; ++n) {
		const double a = Norm3(&v[3 * n]);
		if (a > mx) { mx = a; if (node) *node = n; }
	}
	return mx;
}

double Norm2(const std::vector<double>& v) { double s = 0.; for (std::size_t i = 0; i < v.size(); ++i) s += v[i] * v[i]; return std::sqrt(s); }
double Rms(const Panel& p, const std::vector<double>& v) { return Norm2(v) / std::sqrt(static_cast<double>(p.nNode)); }

struct Info {
	int iters, nTrue, nPsd, factorizations;
	bool converged;
	double res, mu, tol;
};

const double kFloor = 1e-11;   // 수렴 기준의 반올림 몫: 서로 지우는 힘의 크기(scale)의 이 배

// f(x) = 0 을 Newton 꼴 반복으로 푼다. 걸음 d = (K + μM)⁻¹ f, 받아들이면 μ 를 줄이고 못 받아들이면 키운다.
// 정확한 팔(variant 0)은 에너지가 줄 때(또는 잔차가 줄고 에너지가 반올림 안에서 늘지 않을 때), 나머지는 잔차 노름이 줄 때 받아들인다.
Info Equilibrate(const Panel& p, std::vector<double>& x, double BE, int variant, const std::vector<double>& load, double tol, int maxIter, Band& K)
{
	Info info = { 0, 0, 0, 0, false, 0., 0., tol };
	const std::size_t n3 = x.size();
	std::vector<double> f, ft, d, xt(n3);
	double e0 = 0., et = 0., mu = 0., sc = 0., sct = 0.;
	Forces(p, x, BE, variant, load, f, &e0, &sc);
	for (int it = 0; it < maxIter; ++it) {
		info.res = MaxNodal(p, f, 0);
		info.tol = std::max(tol, kFloor * sc);
		if (info.res <= info.tol) { info.converged = true; break; }
		const double fn = Norm2(f);
		bool accepted = false;
		for (int attempt = 0; attempt < 80 && !accepted; ++attempt) {
			int used = -1;
			for (int kind = 0; kind < 2 && used < 0; ++kind) {
				Tangent(p, x, BE, kind, mu, K);
				++info.factorizations;
				if (K.Cholesky()) used = kind;
			}
			if (used < 0) { mu = mu > 0. ? 4. * mu : 1.; continue; }
			SolveNodes(p, K, f, d);
			double slope = 0.;
			for (std::size_t i = 0; i < n3; ++i) { slope += f[i] * d[i]; xt[i] = x[i] + d[i]; }
			Forces(p, xt, BE, variant, load, ft, &et, &sct);
			const double ftn = Norm2(ft);
			bool ok;
			if (variant == kExact) ok = slope > 0. && (et <= e0 - 1e-4 * slope || (ftn < fn && et <= e0 + 1e-12 * (std::fabs(e0) + 1.)));
			else ok = ftn < fn;
			if (ok) {
				x.swap(xt);
				f.swap(ft);
				e0 = et;
				sc = sct;
				if (used == 0) ++info.nTrue; else ++info.nPsd;
				mu = mu > 1e-6 ? 0.25 * mu : 0.;
				accepted = true;
			}
			else mu = mu > 0. ? 4. * mu : 1.;
		}
		++info.iters;
		if (!accepted) break;
	}
	info.res = MaxNodal(p, f, 0);
	info.tol = std::max(tol, kFloor * sc);
	info.converged = info.res <= info.tol;
	info.mu = mu;
	return info;
}

struct Disp { double mx, rms; int node; };

Disp Measure(const Panel& p, const std::vector<double>& u)
{
	Disp d;
	d.mx = MaxNodal(p, u, &d.node);
	d.rms = Rms(p, u);
	return d;
}

double Cosine(const std::vector<double>& a, const std::vector<double>& b)
{
	double ab = 0.;
	for (std::size_t i = 0; i < a.size(); ++i) ab += a[i] * b[i];
	const double na = Norm2(a), nb = Norm2(b);
	return na > 0. && nb > 0. ? ab / (na * nb) : 0.;
}

// 모양의 접힌 각(도) rms · 최대, 인장 변형률 rms · 최대
void ShapeStats(const Panel& p, const std::vector<double>& x, double out[4])
{
	double s2 = 0., mx = 0.;
	for (int s = 0; s < p.nHinge; ++s) {
		double xl[12], th = 0., gE[12], gV[12], H[sxgen::BaraffBend::kHess];
		Gather(p.hinges, s, 4, x, xl);
		HingeEval(xl, kExact, th, gE, gV, H);
		s2 += th * th;
		mx = std::max(mx, std::fabs(th));
	}
	out[0] = std::sqrt(s2 / std::max(1, p.nHinge)) * 180. / kPi;
	out[1] = mx * 180. / kPi;
	double e2 = 0., emx = 0.;
	for (int t = 0; t < p.nTri; ++t) {
		double xl[9], ev[8], C[2], g[18], H[90];
		Gather(p.tris, t, 3, x, xl);
		TriValues(p, t, ev);
		sxgen::BaraffStretch::Eval(xl, ev, C, g, H);
		const double sa = std::sqrt(0.5 * std::fabs((ev[2] - ev[0]) * (ev[5] - ev[1]) - (ev[4] - ev[0]) * (ev[3] - ev[1])));
		for (int c = 0; c < 2; ++c) { e2 += (C[c] / sa) * (C[c] / sa); emx = std::max(emx, std::fabs(C[c] / sa)); }
	}
	out[2] = std::sqrt(e2 / std::max(1, 2 * p.nTri));
	out[3] = emx;
}

double EdgeChange(const Panel& p, const std::vector<double>& xa, const std::vector<double>& xb)
{
	double mx = 0.;
	for (int t = 0; t < p.nTri; ++t)
		for (int e = 0; e < 3; ++e) {
			const int a = p.tris[3 * t + e], b = p.tris[3 * t + (e + 1) % 3];
			double da[3], db[3];
			for (int q = 0; q < 3; ++q) { da[q] = xa[3 * b + q] - xa[3 * a + q]; db[q] = xb[3 * b + q] - xb[3 * a + q]; }
			mx = std::max(mx, std::fabs(Norm3(db) / Norm3(da) - 1.));
		}
	return mx;
}

std::vector<double> Diff(const std::vector<double>& a, const std::vector<double>& b)
{
	std::vector<double> d(a.size());
	for (std::size_t i = 0; i < a.size(); ++i) d[i] = a[i] - b[i];
	return d;
}

} // namespace

int main(int argc, char** argv)
{
	if (argc < 6) { std::printf("사용: measure_frozen_garment_nl <입력> <출력 json> <면밀도> <중력> <BE> [<BE> ...]\n"); return 2; }
	Panel p;
	if (!Read(argv[1], p)) { std::printf("입력을 읽을 수 없다: %s\n", argv[1]); return 1; }
	const double density = std::atof(argv[3]), grav = std::atof(argv[4]);
	Build(p, density);
	const int n3 = 3 * p.nNode;
	Band K;
	K.Init(n3, p.hb);
	int nPinned = 0;
	double wmax = 0.;
	for (int n = 0; n < p.nNode; ++n) { nPinned += p.pinned[n]; wmax = std::max(wmax, p.mass[n] * grav); }
	const double tol = 1e-9 * wmax;
	std::vector<double> gravity(n3, 0.), zero(n3, 0.);
	for (int n = 0; n < p.nNode; ++n) gravity[3 * n + 1] = -p.mass[n] * grav;   // 중력은 −y

	std::FILE* js = std::fopen(argv[2], "w");
	if (!js) { std::printf("출력을 열 수 없다: %s\n", argv[2]); return 1; }
	std::fprintf(js, "{\"nodes\": %d, \"triangles\": %d, \"hinges\": %d, \"pinned_nodes\": %d, \"half_bandwidth\": %d, \"E\": %g, \"SE\": %g,\n"
	             " \"density\": %g, \"gravity\": %g, \"tolerance_dyn\": %.6e, \"max_node_weight_dyn\": %.6e,\n \"cases\": [\n",
	             p.nNode, p.nTri, p.nHinge, nPinned, p.hb, p.E, p.SE, density, grav, tol, wmax);
	std::printf("옷깃: 노드 %d(고정 %d) · 삼각형 %d · 힌지 %d · 반띠폭 %d | 수렴 기준 %.2e dyn (가장 무거운 노드 무게 %.3f dyn 의 1e-9)\n",
	            p.nNode, nPinned, p.nTri, p.nHinge, p.hb, tol, wmax);

	for (int ci = 5; ci < argc; ++ci) {
		const double BE = std::atof(argv[ci]);
		std::vector<double> f, fF, fC, u;

		// ---- A. 옷깃만 다시 드레이프
		std::vector<double> xE(p.x0);
		const Info ie = Equilibrate(p, xE, BE, kExact, gravity, tol, 2000, K);
		double shape[4];
		ShapeStats(p, xE, shape);
		const Disp moved = Measure(p, Diff(xE, p.x0));
		// x_E 에서의 가짜 힘과 선형 응답 둘
		Forces(p, xE, BE, kExact, gravity, f, 0);
		Forces(p, xE, BE, kFrozen, gravity, fF, 0);
		Forces(p, xE, BE, kCorrected, gravity, fC, 0);
		const std::vector<double> df = Diff(fF, f);
		const double dfMax = MaxNodal(p, df, 0), dfcMax = MaxNodal(p, Diff(fC, f), 0);
		Disp uL = { 0., 0., -1 }, uT = { 0., 0., -1 };
		std::vector<double> uLv(n3, 0.), uTv(n3, 0.);
		Tangent(p, xE, BE, 1, 0., K);
		const bool psdOk = K.Cholesky();
		if (psdOk) { SolveNodes(p, K, df, uLv); uL = Measure(p, uLv); }
		Tangent(p, xE, BE, 0, 0., K);
		const bool truePd = K.Cholesky();
		if (truePd) { SolveNodes(p, K, df, uTv); uT = Measure(p, uTv); }
		// 길이 고정 팔 · 보정 팔
		std::vector<double> xF(xE), xC(xE);
		const Info ifz = Equilibrate(p, xF, BE, kFrozen, gravity, tol, 2000, K);
		const Info ic = Equilibrate(p, xC, BE, kCorrected, gravity, tol, 2000, K);
		const std::vector<double> uNLv = Diff(xF, xE);
		const Disp uNL = Measure(p, uNLv), uC = Measure(p, Diff(xC, xE));
		const double edge = EdgeChange(p, xE, xF);
		const double cosT = truePd ? Cosine(uNLv, uTv) : 0., cosL = psdOk ? Cosine(uNLv, uLv) : 0.;

		// ---- B. CLO 모양에서의 죽은 하중 문제
		std::vector<double> dead;
		Forces(p, p.x0, BE, kExact, zero, dead, 0);
		for (int i = 0; i < n3; ++i) dead[i] = -dead[i];
		Forces(p, p.x0, BE, kFrozen, dead, fF, 0);                 // = 가짜 힘 Δf(x₀)
		std::vector<double> uL0v(n3, 0.);
		Disp uL0 = { 0., 0., -1 };
		Tangent(p, p.x0, BE, 1, 0., K);
		if (K.Cholesky()) { SolveNodes(p, K, fF, uL0v); uL0 = Measure(p, uL0v); }   // M4b 의 선형 응답과 같아야 한다
		Tangent(p, p.x0, BE, 0, 0., K);
		const bool truePd0 = K.Cholesky();
		std::vector<double> xB(p.x0);
		const Info ib = Equilibrate(p, xB, BE, kFrozen, dead, tol, 500, K);
		const Disp uB = Measure(p, Diff(xB, p.x0));

		std::printf("  굽힘 %-6g | A: x_E 반복 %d(참 %d · 준 %d)%s 잔차 %.1e, CLO 에서 %.2f cm, 각 rms %.1f° 최대 %.1f°, 변형률 rms %.4f, 참 접선 %s | "
		            "가짜 힘 최대 %.3e dyn | u_NL %.4e cm (반복 %d%s) · u_T %.4e · u_L %.4e → NL/T %.4f · NL/L %.4f · cos(NL,T) %.5f | 보정 %.1e cm | "
		            "B: x0 참 접선 %s, u_L(x0) %.4e, 길이 고정 %s(반복 %d, 잔차 %.1e) u %.4e cm\n",
		            BE, ie.iters, ie.nTrue, ie.nPsd, ie.converged ? "" : " [미수렴]", ie.res, moved.mx, shape[0], shape[1], shape[2], truePd ? "정부호" : "부정치",
		            dfMax, uNL.mx, ifz.iters, ifz.converged ? "" : " [미수렴]", uT.mx, uL.mx, uT.mx > 0. ? uNL.mx / uT.mx : 0., uL.mx > 0. ? uNL.mx / uL.mx : 0., cosT,
		            uC.mx, truePd0 ? "정부호" : "부정치", uL0.mx, ib.converged ? "수렴" : "미수렴", ib.iters, ib.res, uB.mx);
		std::fflush(stdout);
		std::fprintf(js,
		             "  {\"BE\": %g,\n"
		             "   \"A\": {\"exact\": {\"converged\": %s, \"iterations\": %d, \"true_tangent_steps\": %d, \"psd_steps\": %d, \"factorizations\": %d, \"residual_dyn\": %.6e, \"tolerance_dyn\": %.6e},\n"
		             "         \"moved_from_clo_max_cm\": %.9e, \"theta_rms_deg\": %.6f, \"theta_max_deg\": %.6f, \"strain_rms\": %.6e, \"strain_max\": %.6e,\n"
		             "         \"true_tangent_pd_at_xE\": %s, \"df_max_dyn\": %.9e, \"df_corrected_max_dyn\": %.3e,\n"
		             "         \"frozen\": {\"converged\": %s, \"iterations\": %d, \"residual_dyn\": %.6e, \"tolerance_dyn\": %.6e},\n"
		             "         \"u_nl_max_cm\": %.9e, \"u_nl_rms_cm\": %.9e, \"u_nl_max_node\": %d, \"edge_length_change_max\": %.9e,\n"
		             "         \"u_true_tangent_max_cm\": %.9e, \"u_true_tangent_rms_cm\": %.9e, \"u_psd_max_cm\": %.9e, \"u_psd_rms_cm\": %.9e,\n"
		             "         \"cos_nl_true\": %.9f, \"cos_nl_psd\": %.9f,\n"
		             "         \"corrected\": {\"converged\": %s, \"iterations\": %d, \"u_max_cm\": %.3e}},\n"
		             "   \"B\": {\"true_tangent_pd_at_x0\": %s, \"u_psd_at_x0_max_cm\": %.9e, \"frozen\": {\"converged\": %s, \"iterations\": %d, \"residual_dyn\": %.6e, \"tolerance_dyn\": %.6e, "
		             "\"u_max_cm\": %.9e, \"u_rms_cm\": %.9e}}}%s\n",
		             BE, ie.converged ? "true" : "false", ie.iters, ie.nTrue, ie.nPsd, ie.factorizations, ie.res, ie.tol,
		             moved.mx, shape[0], shape[1], shape[2], shape[3],
		             truePd ? "true" : "false", dfMax, dfcMax,
		             ifz.converged ? "true" : "false", ifz.iters, ifz.res, ifz.tol,
		             uNL.mx, uNL.rms, uNL.node, edge,
		             uT.mx, uT.rms, uL.mx, uL.rms, cosT, cosL,
		             ic.converged ? "true" : "false", ic.iters, uC.mx,
		             truePd0 ? "true" : "false", uL0.mx, ib.converged ? "true" : "false", ib.iters, ib.res, ib.tol, uB.mx, uB.rms, ci + 1 < argc ? "," : "");
		std::fflush(js);
	}
	std::fprintf(js, " ]}\n");
	std::fclose(js);
	return 0;
}
