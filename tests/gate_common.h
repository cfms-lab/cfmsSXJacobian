//---------------------------------------------------------------------------
// gate_common.h — 게이트 시험들이 함께 쓰는 것: 표본 읽기 · 오차 · 유한차분 · 판정 출력 · 격자 메시.
//---------------------------------------------------------------------------
#pragma once

#include "../src/sx_condition_force.h"
#include "../src/sx_stencil.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace gate {

struct Hinge { double x[12]; double theta0; double v[12]; };

// results/hinges.txt 한 줄 = x(12) · θ0 · v(12)
inline std::vector<Hinge> ReadHinges(const char* path)
{
	std::vector<Hinge> r;
	std::ifstream f(path);
	std::string line;
	while (std::getline(f, line)) {
		if (line.empty()) continue;
		std::istringstream is(line);
		Hinge h;
		for (int i = 0; i < 12; ++i) is >> h.x[i];
		is >> h.theta0;
		for (int i = 0; i < 12; ++i) is >> h.v[i];
		if (is) r.push_back(h);
	}
	return r;
}

// 삼각형 표본(인장 · 전단 게이트, Gauss–Newton 게이트)
struct Tri { double x[9]; double uv[6]; double v[9]; };

// results/triangles.txt 한 줄 = x(9) · uv(6) · v(9)
inline std::vector<Tri> ReadTriangles(const char* path)
{
	std::vector<Tri> r;
	std::ifstream f(path);
	std::string line;
	while (std::getline(f, line)) {
		if (line.empty()) continue;
		std::istringstream is(line);
		Tri t;
		for (int i = 0; i < 9; ++i) is >> t.x[i];
		for (int i = 0; i < 6; ++i) is >> t.uv[i];
		for (int i = 0; i < 9; ++i) is >> t.v[i];
		if (is) r.push_back(t);
	}
	return r;
}

inline double UvArea(const double* uv)
{
	const double du1 = uv[2] - uv[0], dv1 = uv[3] - uv[1], du2 = uv[4] - uv[0], dv2 = uv[5] - uv[1];
	return 0.5 * std::fabs(du1 * dv2 - du2 * dv1);
}

// 빈칸으로 나뉜 수의 줄들(results/mma_ref.txt)
inline std::vector<std::vector<double> > ReadRows(const char* path)
{
	std::vector<std::vector<double> > r;
	std::ifstream f(path);
	std::string line;
	while (std::getline(f, line)) {
		if (line.empty()) continue;
		std::vector<double> row;
		const char* p = line.c_str();
		char* end = 0;
		for (;;) {
			const double v = std::strtod(p, &end);
			if (end == p) break;
			row.push_back(v);
			p = end;
		}
		r.push_back(row);
	}
	return r;
}

inline double MaxAbs(const double* a, int n) { double m = 0.; for (int i = 0; i < n; ++i) m = std::max(m, std::fabs(a[i])); return m; }
inline double MaxAbsDiff(const double* a, const double* b, int n) { double m = 0.; for (int i = 0; i < n; ++i) m = std::max(m, std::fabs(a[i] - b[i])); return m; }
inline double Rel(double diff, double scale) { return scale > 0. ? diff / scale : diff; }
inline double EdgeLength(const double* x) { const double dx = x[6] - x[3], dy = x[7] - x[4], dz = x[8] - x[5]; return std::sqrt(dx * dx + dy * dy + dz * dz); }

struct Stat {   // 최댓값과 중앙값
	std::vector<double> v;
	void Add(double x) { v.push_back(x); }
	double Max() const { return v.empty() ? 0. : *std::max_element(v.begin(), v.end()); }
	double Median() const { if (v.empty()) return 0.; std::vector<double> c(v); std::nth_element(c.begin(), c.begin() + c.size() / 2, c.end()); return c[c.size() / 2]; }
};

typedef void (*EvalFn)(const double*, const double*, double*, double*, double*);

// 유한차분: g ≈ ∂C/∂x, H ≈ ∂g/∂x (위 삼각). 분모는 실제로 표현된 간격(좌표 ~1 m 에서 2·step 과 다르다).
template <int N>
void FiniteDiff(EvalFn f, const double* x, const double* ev, double step, double* gFD, double* HFD)
{
	const int kH = N * (N + 1) / 2;
	double xp[N], xm[N], Cp, Cm, gp[N], gm[N], Hp[kH], Hm[kH];
	for (int b = 0; b < N; ++b) {
		std::memcpy(xp, x, sizeof xp); std::memcpy(xm, x, sizeof xm);
		xp[b] += step; xm[b] -= step;
		const double h2 = xp[b] - xm[b];
		f(xp, ev, &Cp, gp, Hp);
		f(xm, ev, &Cm, gm, Hm);
		gFD[b] = (Cp - Cm) / h2;
		for (int a = 0; a <= b; ++a) HFD[sx::UpperIndex<N>(a, b)] = (gp[a] - gm[a]) / h2;
	}
}

// 판정 출력과 실패 수
struct Report {
	std::FILE* json;
	int fail;
	explicit Report(const char* path) : json(std::fopen(path, "w")), fail(0) {}
	~Report() { if (json) std::fclose(json); }
	void Gate(const char* name, bool ok, const char* fmt, double value, double limit)
	{
		std::printf("  %-4s %s  ", name, ok ? "PASS" : "FAIL");
		std::printf(fmt, value);
		std::printf("  (기준 %.0e)\n", limit);
		if (!ok) ++fail;
	}
};

// 시간 측정용 격자 천: G × G 노드, 간격 1 cm, 살짝 물결치는 면. 삼각형(반시계) → 방향 있는 변 → 힌지 (h, i, j, k)
struct Grid {
	int nNode, nHinge, nTri;
	std::vector<double> X, V;
	std::vector<int> hinges;   // 힌지마다 4 개
	std::vector<int> tris;     // 삼각형마다 3 개(반시계). 인장 · 전단 게이트가 쓴다
};

inline Grid MakeGrid(int G, std::mt19937_64& rng)
{
	Grid r;
	r.nNode = G * G;
	const int nNode = r.nNode;
	r.X.resize(3 * nNode); r.V.resize(3 * nNode);
	std::uniform_real_distribution<double> noise(-1., 1.);
	for (int j = 0; j < G; ++j)
		for (int i = 0; i < G; ++i) {
			const int n = j * G + i;
			r.X[3 * n] = 0.01 * i + 0.3; r.X[3 * n + 1] = 0.01 * j + 0.8;
			r.X[3 * n + 2] = 0.02 * std::sin(0.2 * i) * std::cos(0.15 * j) + 0.001 * noise(rng);
			for (int c = 0; c < 3; ++c) r.V[3 * n + c] = 0.1 * noise(rng);
		}
	std::vector<int> tri;
	for (int j = 0; j + 1 < G; ++j)
		for (int i = 0; i + 1 < G; ++i) {
			const int a = j * G + i, b = a + 1, c = a + G + 1, d = a + G;
			const int t1[3] = { a, b, c }, t2[3] = { a, c, d };
			tri.insert(tri.end(), t1, t1 + 3);
			tri.insert(tri.end(), t2, t2 + 3);
		}
	std::unordered_map<long long, int> third;   // 방향 변 p→q → 세 번째 꼭짓점
	for (std::size_t t = 0; t < tri.size(); t += 3)
		for (int e = 0; e < 3; ++e) third[static_cast<long long>(tri[t + e]) * nNode + tri[t + (e + 1) % 3]] = tri[t + (e + 2) % 3];
	for (std::size_t t = 0; t < tri.size(); t += 3)
		for (int e = 0; e < 3; ++e) {
			const int p = tri[t + e], q = tri[t + (e + 1) % 3];
			if (p > q) continue;
			std::unordered_map<long long, int>::const_iterator it = third.find(static_cast<long long>(q) * nNode + p);
			if (it == third.end()) continue;
			const int hh[4] = { tri[t + (e + 2) % 3], p, q, it->second };
			r.hinges.insert(r.hinges.end(), hh, hh + 4);
		}
	r.nHinge = static_cast<int>(r.hinges.size() / 4);
	r.tris = tri;
	r.nTri = static_cast<int>(tri.size() / 3);
	return r;
}

// 명세의 선언대로 SXMesh 변수를 만들고 격자 · 쉬는 각(theta0[s])을 채운다
template <class Spec>
void FillMesh(SXMesh& mesh, const Grid& g, const std::vector<double>& theta0)
{
	sx::CreateVariables<Spec>(mesh, g.nNode, g.nHinge);
	SXNodeNamedContainer& px = mesh.GetNodeV(Spec::PositionVar());
	SXNodeNamedContainer& pv = mesh.GetNodeV(Spec::VelocityVar());
	SXElementNamedContainer& ph = mesh.GetElementV(Spec::StencilVar());
	SXNodeNamedContainer& pt0 = mesh.GetNodeV(Spec::ElemVar(0));
	for (int n = 0; n < g.nNode; ++n) { px.Add(n, &g.X[3 * n], 3); pv.Add(n, &g.V[3 * n], 3); }
	for (int s = 0; s < g.nHinge; ++s) { ph.Add(s, &g.hinges[4 * s], 4); pt0[s][0] = theta0[s]; }
}

} // namespace gate
