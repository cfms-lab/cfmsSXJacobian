//---------------------------------------------------------------------------
// measure_frozen_garment.cpp — 드레이프된 실제 옷의 한 패널(옷깃)에서 길이 고정 근사의 가짜 힘과 그 변위를 푼다
// (docs/MEASURE_FROZEN_FORCE.md 「M4」, 판정 밖 측정). tests/measure_frozen_garment.py 가 패널을 뽑아 부른다.
//   입력(글): E SE / nNode nTri / 노드마다 x y z u v 고정여부 / 삼각형마다 a b c
//   굽힘: Baraff C = θ(θ₀ = 0), 힌지 강성 k = BE · w̄ / 2 (w̄ = 3|e|² / (A₁ + A₂), 쉬는 uv) — 엔진 Bridson 굽힘과 작은 각에서 같다
//   Δf = f_길이고정 − f_정확,  K u = Δf,  K = 인장(C > 0 판, k = E) + Provot 전단(C > 0 판, k = SE) + 굽힘 Gauss–Newton
//   사용:  measure_frozen_garment <입력> <출력 json> [--dump <접두사>] <BE> [<BE> ...]
//          --dump 를 주면 굽힘 값마다 <접두사>_<순번>.txt 에 노드별 변위(3) · 변 길이 변화의 |평균| · 부호 있는 평균 · 최대를 쓴다(그림용)
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/BaraffBend.h"
#include "../generated/BaraffStretch.h"
#include "../generated/ProvotShear.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

const double kPi = 3.14159265358979323846;
const double kDiag = static_cast<double>(1.41421356237f);   // 엔진 Provot 전단의 b (float 리터럴)

struct Panel {
	int nNode, nTri, nHinge;
	double E, SE;
	std::vector<double> x, uv;
	std::vector<char> pinned, boundary;
	std::vector<int> tris, hinges;
	std::vector<double> wbar;
	std::unordered_map<long long, int> third;   // 방향 변 (a → b) → 맞은편 꼭짓점
};

inline double Norm3(const double* a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }

bool Read(const char* path, Panel& p)
{
	std::FILE* f = std::fopen(path, "r");
	if (!f) return false;
	bool ok = std::fscanf(f, "%lf %lf %d %d", &p.E, &p.SE, &p.nNode, &p.nTri) == 4;
	p.x.resize(3 * static_cast<std::size_t>(p.nNode));
	p.uv.resize(2 * static_cast<std::size_t>(p.nNode));
	p.pinned.assign(p.nNode, 0);
	for (int n = 0; ok && n < p.nNode; ++n) {
		int pin = 0;
		ok = std::fscanf(f, "%lf %lf %lf %lf %lf %d", &p.x[3 * n], &p.x[3 * n + 1], &p.x[3 * n + 2], &p.uv[2 * n], &p.uv[2 * n + 1], &pin) == 6;
		p.pinned[n] = static_cast<char>(pin);
	}
	p.tris.resize(3 * static_cast<std::size_t>(p.nTri));
	for (int t = 0; ok && t < p.nTri; ++t) ok = std::fscanf(f, "%d %d %d", &p.tris[3 * t], &p.tris[3 * t + 1], &p.tris[3 * t + 2]) == 3;
	std::fclose(f);
	return ok;
}

// 힌지 (h, i, j, k): 삼각형 A = (h, i, j), B = (k, j, i) — tests/test_integration.cpp 의 MakeCloth 와 같은 짜임
void BuildHinges(Panel& p, int& nBoundaryEdge)
{
	std::unordered_map<long long, int>& third = p.third;
	for (int t = 0; t < p.nTri; ++t)
		for (int e = 0; e < 3; ++e)
			third[static_cast<long long>(p.tris[3 * t + e]) * p.nNode + p.tris[3 * t + (e + 1) % 3]] = p.tris[3 * t + (e + 2) % 3];
	p.boundary.assign(p.nNode, 0);
	nBoundaryEdge = 0;
	for (int t = 0; t < p.nTri; ++t)
		for (int e = 0; e < 3; ++e) {
			const int a = p.tris[3 * t + e], b = p.tris[3 * t + (e + 1) % 3];
			std::unordered_map<long long, int>::const_iterator it = third.find(static_cast<long long>(b) * p.nNode + a);
			if (it == third.end()) {
				++nBoundaryEdge;
				p.boundary[a] = 1;
				p.boundary[b] = 1;
				continue;
			}
			if (a > b) continue;
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
}

inline void Gather(const std::vector<int>& idx, int s, int nn, const std::vector<double>& x, double* xl)
{
	for (int a = 0; a < nn; ++a) { const int n = idx[nn * s + a]; for (int t = 0; t < 3; ++t) xl[3 * a + t] = x[3 * n + t]; }
}

template <int N>
inline int Upper(int a, int b) { if (a > b) std::swap(a, b); return a * N - a * (a - 1) / 2 + (b - a); }

bool Cholesky(std::vector<double>& A, int n)   // 아래 삼각, 행 우선 n × n
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j <= i; ++j) {
			double s = A[static_cast<std::size_t>(i) * n + j];
			for (int k = 0; k < j; ++k) s -= A[static_cast<std::size_t>(i) * n + k] * A[static_cast<std::size_t>(j) * n + k];
			if (j < i) A[static_cast<std::size_t>(i) * n + j] = s / A[static_cast<std::size_t>(j) * n + j];
			else { if (!(s > 0.)) return false; A[static_cast<std::size_t>(i) * n + i] = std::sqrt(s); }
		}
	return true;
}

void Solve(const std::vector<double>& L, int n, std::vector<double>& b)
{
	for (int i = 0; i < n; ++i) {
		double s = b[i];
		for (int k = 0; k < i; ++k) s -= L[static_cast<std::size_t>(i) * n + k] * b[k];
		b[i] = s / L[static_cast<std::size_t>(i) * n + i];
	}
	for (int i = n - 1; i >= 0; --i) {
		double s = b[i];
		for (int k = i + 1; k < n; ++k) s -= L[static_cast<std::size_t>(k) * n + i] * b[k];
		b[i] = s / L[static_cast<std::size_t>(i) * n + i];
	}
}

struct Stat {
	double mx, ss; int n;
	Stat() : mx(0.), ss(0.), n(0) {}
	void Add(double v) { mx = std::max(mx, std::fabs(v)); ss += v * v; ++n; }
	double Rms() const { return n ? std::sqrt(ss / n) : 0.; }
};

} // namespace

int main(int argc, char** argv)
{
	if (argc < 4) { std::printf("사용: measure_frozen_garment <입력> <출력 json> [--dump <접두사>] <BE> [<BE> ...]\n"); return 2; }
	const char* dump = 0;   // 있으면 굽힘 값마다 <접두사>_<순번>.txt 로 노드별 변위 · 변 길이 변화를 쓴다(그림용)
	int first = 3;
	if (argc >= 6 && std::string(argv[3]) == "--dump") { dump = argv[4]; first = 5; }
	Panel p;
	if (!Read(argv[1], p)) { std::printf("입력을 읽을 수 없다: %s\n", argv[1]); return 1; }
	int nBoundaryEdge = 0;
	BuildHinges(p, nBoundaryEdge);
	const int n3 = 3 * p.nNode;
	int nPinned = 0, nBoundaryNode = 0;
	for (int n = 0; n < p.nNode; ++n) { nPinned += p.pinned[n]; nBoundaryNode += p.boundary[n]; }

	// 드레이프된 모양의 접힌 각과 막 변형률
	Stat theta;
	int over10 = 0, over30 = 0;
	std::vector<double> th(p.nHinge), gE(12 * static_cast<std::size_t>(p.nHinge)), gF(gE.size()), gC(gE.size());
	double idErr = 0., gMax = 0.;
	for (int s = 0; s < p.nHinge; ++s) {
		double xl[12], H[sxgen::BaraffBend::kHess], C = 0.;
		const double ev[1] = { 0. };
		Gather(p.hinges, s, 4, p.x, xl);
		sxgen::BaraffBend::Eval(xl, ev, &C, &gE[12 * s], H);
		sxgen::BaraffBend::EvalFrozen(xl, ev, &C, &gF[12 * s], H);
		th[s] = C;
		theta.Add(C);
		over10 += std::fabs(C) > 10. * kPi / 180.;
		over30 += std::fabs(C) > 30. * kPi / 180.;
		double e[3] = { xl[6] - xl[3], xl[7] - xl[4], xl[8] - xl[5] };
		const double len = Norm3(e);
		const double sc = std::sin(C) * std::cos(C) / len;
		for (int a = 0; a < 4; ++a)
			for (int t = 0; t < 3; ++t) {
				const double expect = (a == 2 ? sc : a == 1 ? -sc : 0.) * e[t] / len;
				idErr = std::max(idErr, std::fabs(gF[12 * s + 3 * a + t] - gE[12 * s + 3 * a + t] - expect));
				gMax = std::max(gMax, std::fabs(gE[12 * s + 3 * a + t]));
				gC[12 * s + 3 * a + t] = gF[12 * s + 3 * a + t];
			}
		for (int a = 1; a <= 2; ++a) {   // 보정: 변 방향 성분을 뺀다
			const double d = (gC[12 * s + 3 * a] * e[0] + gC[12 * s + 3 * a + 1] * e[1] + gC[12 * s + 3 * a + 2] * e[2]) / len;
			for (int t = 0; t < 3; ++t) gC[12 * s + 3 * a + t] -= d * e[t] / len;
		}
	}
	Stat strain;
	for (int t = 0; t < p.nTri; ++t) {
		double xl[9], uvl[6], C[2], g[18], H[90];
		Gather(p.tris, t, 3, p.x, xl);
		for (int a = 0; a < 3; ++a) { uvl[2 * a] = p.uv[2 * p.tris[3 * t + a]]; uvl[2 * a + 1] = p.uv[2 * p.tris[3 * t + a] + 1]; }
		sxgen::BaraffStretch::Eval(xl, uvl, C, g, H);
		const double sa = std::sqrt(0.5 * std::fabs((uvl[2] - uvl[0]) * (uvl[5] - uvl[1]) - (uvl[4] - uvl[0]) * (uvl[3] - uvl[1])));
		strain.Add(C[0] / sa);
		strain.Add(C[1] / sa);
	}

	std::FILE* js = std::fopen(argv[2], "w");
	if (!js) { std::printf("출력을 열 수 없다: %s\n", argv[2]); return 1; }
	std::fprintf(js, "{\"nodes\": %d, \"triangles\": %d, \"hinges\": %d, \"boundary_edges\": %d, \"boundary_nodes\": %d, \"pinned_nodes\": %d,\n"
	             " \"E\": %g, \"SE\": %g, \"theta_max_deg\": %.6f, \"theta_rms_deg\": %.6f, \"hinges_over_10deg\": %d, \"hinges_over_30deg\": %d,\n"
	             " \"stretch_strain_max\": %.6e, \"stretch_strain_rms\": %.6e, \"hinge_identity_err\": %.3e,\n \"cases\": [\n",
	             p.nNode, p.nTri, p.nHinge, nBoundaryEdge, nBoundaryNode, nPinned, p.E, p.SE, theta.mx * 180. / kPi, theta.Rms() * 180. / kPi, over10,
	             over30, strain.mx, strain.Rms(), idErr / gMax);
	std::printf("옷깃: 노드 %d(가장자리 %d, 고정 %d) · 삼각형 %d · 힌지 %d · 가장자리 변 %d | 접힌 각 최대 %.1f° rms %.1f° (10° 넘는 힌지 %d, 30° 넘는 힌지 %d) | "
	            "인장 변형률 최대 %.3f rms %.4f | 힌지 식 오차 %.1e\n",
	            p.nNode, nBoundaryNode, nPinned, p.nTri, p.nHinge, nBoundaryEdge, theta.mx * 180. / kPi, theta.Rms() * 180. / kPi, over10, over30, strain.mx,
	            strain.Rms(), idErr / gMax);

	for (int ci = first; ci < argc; ++ci) {
		const double BE = std::atof(argv[ci]);
		std::vector<double> fE(n3, 0.), dF(n3, 0.), dC(n3, 0.), K(static_cast<std::size_t>(n3) * n3, 0.);
		// 막: 인장 · Provot 전단, C > 0 판
		for (int t = 0; t < p.nTri; ++t) {
			double xl[9], ev[8], C[2], g[18], H[90];
			Gather(p.tris, t, 3, p.x, xl);
			for (int a = 0; a < 3; ++a) { ev[2 * a] = p.uv[2 * p.tris[3 * t + a]]; ev[2 * a + 1] = p.uv[2 * p.tris[3 * t + a] + 1]; }
			ev[6] = kDiag;
			ev[7] = kDiag;
			for (int kind = 0; kind < 2; ++kind) {
				if (kind == 0) sxgen::BaraffStretch::Eval(xl, ev, C, g, H);
				else sxgen::ProvotShear::Eval(xl, ev, C, g, H);
				const double k = kind == 0 ? p.E : p.SE;
				for (int c = 0; c < 2; ++c)
					for (int a = 0; a < 9; ++a) {
						const int ra = 3 * p.tris[3 * t + a / 3] + a % 3;
						for (int b = 0; b < 9; ++b) {
							const int rb = 3 * p.tris[3 * t + b / 3] + b % 3;
							if (ra >= rb) K[static_cast<std::size_t>(ra) * n3 + rb] += k * (g[9 * c + a] * g[9 * c + b] + std::max(C[c], 0.) * H[45 * c + Upper<9>(a, b)]);
						}
					}
			}
		}
		// 굽힘: 정확한 힘, 가짜 힘, Gauss–Newton 행렬
		for (int s = 0; s < p.nHinge; ++s) {
			const double k = BE * p.wbar[s] / 2.;
			for (int a = 0; a < 12; ++a) {
				const int ra = 3 * p.hinges[4 * s + a / 3] + a % 3;
				fE[ra] += -k * th[s] * gE[12 * s + a];
				dF[ra] += -k * th[s] * (gF[12 * s + a] - gE[12 * s + a]);
				dC[ra] += -k * th[s] * (gC[12 * s + a] - gE[12 * s + a]);
				for (int b = 0; b < 12; ++b) {
					const int rb = 3 * p.hinges[4 * s + b / 3] + b % 3;
					if (ra >= rb) K[static_cast<std::size_t>(ra) * n3 + rb] += k * gE[12 * s + a] * gE[12 * s + b];
				}
			}
		}
		for (int n = 0; n < p.nNode; ++n) {
			if (!p.pinned[n]) continue;
			for (int t = 0; t < 3; ++t) {
				const int d = 3 * n + t;
				for (int j = 0; j < d; ++j) K[static_cast<std::size_t>(d) * n3 + j] = 0.;
				for (int i = d + 1; i < n3; ++i) K[static_cast<std::size_t>(i) * n3 + d] = 0.;
				K[static_cast<std::size_t>(d) * n3 + d] = 1.;
			}
		}
		Stat fe, dfAll, dfBnd, dfIn;
		for (int n = 0; n < p.nNode; ++n) {
			fe.Add(Norm3(&fE[3 * n]));
			const double v = Norm3(&dF[3 * n]);
			dfAll.Add(v);
			(p.boundary[n] ? dfBnd : dfIn).Add(v);
		}
		std::vector<double> u(dF), uc(dC);
		for (int n = 0; n < p.nNode; ++n) if (p.pinned[n]) for (int t = 0; t < 3; ++t) { u[3 * n + t] = 0.; uc[3 * n + t] = 0.; }
		const bool ok = Cholesky(K, n3);
		Stat uAll, uBnd, uIn, ucAll, dStrain;
		int nMax = -1;
		if (ok) {
			Solve(K, n3, u);
			Solve(K, n3, uc);
			for (int n = 0; n < p.nNode; ++n) {
				const double v = Norm3(&u[3 * n]);
				if (v > uAll.mx) nMax = n;
				uAll.Add(v);
				(p.boundary[n] ? uBnd : uIn).Add(v);
				ucAll.Add(Norm3(&uc[3 * n]));
			}
			// 변위가 만드는 변 길이의 상대 변화. 노드마다는 그 노드에 닿는 변(한 번씩)의 |변화| 평균 · 부호 있는 평균 · 최대
			std::vector<double> sAbs(p.nNode, 0.), sSigned(p.nNode, 0.), sMax(p.nNode, 0.);
			std::vector<int> sCount(p.nNode, 0);
			for (int t = 0; t < p.nTri; ++t)
				for (int e = 0; e < 3; ++e) {
					const int a = p.tris[3 * t + e], b = p.tris[3 * t + (e + 1) % 3];
					double d0[3], d1[3];
					for (int q = 0; q < 3; ++q) { d0[q] = p.x[3 * b + q] - p.x[3 * a + q]; d1[q] = d0[q] + u[3 * b + q] - u[3 * a + q]; }
					const double ds = Norm3(d1) / Norm3(d0) - 1.;
					dStrain.Add(ds);
					if (a > b && p.third.count(static_cast<long long>(b) * p.nNode + a)) continue;   // 안쪽 변은 한 번만
					const int ends[2] = { a, b };
					for (int q = 0; q < 2; ++q) {
						sAbs[ends[q]] += std::fabs(ds);
						sSigned[ends[q]] += ds;
						sMax[ends[q]] = std::max(sMax[ends[q]], std::fabs(ds));
						++sCount[ends[q]];
					}
				}
			if (dump) {   // 그림용: 노드마다 변위 3 · 변 길이 변화의 |평균| · 부호 있는 평균 · 최대
				char path[1024];
				std::snprintf(path, sizeof(path), "%s_%d.txt", dump, ci - first);
				std::FILE* df = std::fopen(path, "w");
				if (df) {
					for (int n = 0; n < p.nNode; ++n)
						std::fprintf(df, "%.9e %.9e %.9e %.9e %.9e %.9e\n", u[3 * n], u[3 * n + 1], u[3 * n + 2], sCount[n] ? sAbs[n] / sCount[n] : 0.,
						             sCount[n] ? sSigned[n] / sCount[n] : 0., sMax[n]);
					std::fclose(df);
				}
			}
		}
		std::printf("  굽힘 %-7g | 정확한 굽힘 힘 최대 %.3e dyn · 가짜 힘 최대 %.3e (가장자리 %.3e · 안쪽 %.3e) = 정확한 힘의 %.3f 배 | %s 변위 최대 %.3e cm = %.4f µm "
		            "(노드 %d, %s; 가장자리 최대 %.3e · 안쪽 최대 %.3e · rms %.3e) | 변 길이 변화 최대 %.2e | 보정 팔 변위 최대 %.1e cm\n",
		            BE, fe.mx, dfAll.mx, dfBnd.mx, dfIn.mx, fe.mx > 0. ? dfAll.mx / fe.mx : 0., ok ? "풀림," : "행렬이 양의 정부호가 아님,", uAll.mx,
		            uAll.mx * 1e4, nMax, nMax >= 0 && p.boundary[nMax] ? "가장자리" : "안쪽", uBnd.mx, uIn.mx, uAll.Rms(), dStrain.mx, ucAll.mx);
		std::fprintf(js, "  {\"BE\": %g, \"solved\": %s, \"f_exact_bend_max\": %.9e, \"df_max\": %.9e, \"df_boundary_max\": %.9e, \"df_interior_max\": %.9e, "
		             "\"df_rms\": %.9e, \"u_max_cm\": %.9e, \"u_max_node\": %d, \"u_max_on_boundary\": %s, \"u_boundary_max_cm\": %.9e, "
		             "\"u_interior_max_cm\": %.9e, \"u_rms_cm\": %.9e, \"edge_length_change_max\": %.9e, \"u_corrected_max_cm\": %.3e}%s\n",
		             BE, ok ? "true" : "false", fe.mx, dfAll.mx, dfBnd.mx, dfIn.mx, dfAll.Rms(), uAll.mx, nMax,
		             nMax >= 0 && p.boundary[nMax] ? "true" : "false", uBnd.mx, uIn.mx, uAll.Rms(), dStrain.mx, ucAll.mx, ci + 1 < argc ? "," : "");
	}
	std::fprintf(js, " ]}\n");
	std::fclose(js);
	return 0;
}
