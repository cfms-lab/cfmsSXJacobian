//---------------------------------------------------------------------------
// sx_condition_force.h — 조건 함수에서 힘 · Jacobian · 감쇠를 조립한다 (손으로 한 번만 쓰는 공통 코드).
//
// Baraff & Witkin 1998
//   (7)  f  = −k ∇C C
//   (8)  K  = ∂f/∂x = −k (∇C ∇Cᵀ + C ∇²C)                      대칭
//   (11) d  = −kd ∇C Ċ,   Ċ = ∇C · v
//   (12) ∂d/∂x = −kd (∇C (∇²C v)ᵀ + Ċ ∇²C)                      첫 항이 비대칭 — 원 논문은 버린다
//        ∂d/∂v = −kd ∇C ∇Cᵀ                                     대칭
// 벡터 조건(성분 nComp 개)은 성분마다 더한다.
// 입력 C · ∇C · ∇²C 는 생성된 코드(generated/<명세>.h)가 낸다. ∇²C 는 위 삼각 행 우선.
//
// 기본 조립 경로는 이 파일이 아니라 생성된 EvalImplicit(평가 + 식 (6) 조립을 한 함수로)다(사용자 결정 2026-11-13,
// docs/GATE_ASSEMBLY.md). 여기 두 함수는 각 행렬을 따로 볼 때(AssembleConditionForce)와 검산(AssembleImplicit)에 쓴다.
//---------------------------------------------------------------------------
#pragma once

namespace sx {

// 위 삼각(대각 포함, 행 우선) 저장에서 (a, b) 의 위치
template <int N>
inline int UpperIndex(int a, int b)
{
	if (a > b) { const int t = a; a = b; b = t; }
	return a * N - a * (a - 1) / 2 + (b - a);
}

template <int N>
struct ConditionForce {
	double f[N];        // 탄성력
	double K[N * N];    // ∂f/∂x (행 우선)
	double d[N];        // 감쇠력
	double Dx[N * N];   // ∂d/∂x
	double Dv[N * N];   // ∂d/∂v
};

// keepAsymDamping = false 이면 ∂d/∂x 의 비대칭 항 ∇C (∇²C v)ᵀ 를 버린다(원 논문 방식).
template <int N>
inline void AssembleConditionForce(int nComp, const double* C, const double* g, const double* H,
                                   const double* v, double k, double kd, bool keepAsymDamping,
                                   ConditionForce<N>& out)
{
	const int nHess = N * (N + 1) / 2;
	for (int a = 0; a < N; ++a) { out.f[a] = 0.; out.d[a] = 0.; }
	for (int a = 0; a < N * N; ++a) { out.K[a] = 0.; out.Dx[a] = 0.; out.Dv[a] = 0.; }

	for (int c = 0; c < nComp; ++c) {
		const double* gc = g + c * N;
		const double* Hc = H + c * nHess;
		double Cdot = 0.;
		for (int a = 0; a < N; ++a) Cdot += gc[a] * v[a];
		double Hv[N];
		for (int a = 0; a < N; ++a) {
			double s = 0.;
			for (int b = 0; b < N; ++b) s += Hc[UpperIndex<N>(a, b)] * v[b];
			Hv[a] = s;
		}
		for (int a = 0; a < N; ++a) {
			out.f[a] += -k * C[c] * gc[a];
			out.d[a] += -kd * Cdot * gc[a];
			for (int b = 0; b < N; ++b) {
				const double Hab = Hc[UpperIndex<N>(a, b)];
				const double gg  = gc[a] * gc[b];
				out.K[a * N + b]  += -k * (gg + C[c] * Hab);
				out.Dv[a * N + b] += -kd * gg;
				out.Dx[a * N + b] += -kd * (Cdot * Hab + (keepAsymDamping ? gc[a] * Hv[b] : 0.));
			}
		}
	}
}

// 원문 식 (6) 에 M 을 곱한 backward Euler 선형계 (M − h ∂f/∂v − h² ∂f/∂x) Δv = h (f + h ∂f/∂x v) 의
// 요소 기여만 바로 모은다. 감쇠는 비대칭 항을 버린 원문 꼴이라 A 가 대칭이다(위 삼각만).
//   A = (h·kd + h²·k) ∇C∇Cᵀ + h²(k·C + kd·Ċ) ∇²C
//   b = −h(k·C + kd·Ċ) ∇C − h²[ k·Ċ ∇C + (k·C + kd·Ċ) ∇²C v ]
// K · ∂d/∂x · ∂d/∂v 를 따로 만들지 않는다. docs/GATE_ASSEMBLY.md.
template <int N>
struct ImplicitBlock {
	double f[N];                  // 전체 힘 f + d
	double A[N * (N + 1) / 2];    // −h ∂f/∂v − h² ∂f/∂x, 위 삼각 행 우선
	double b[N];                  // h (f + h ∂f/∂x v)
};

template <int N>
inline void AssembleImplicit(int nComp, const double* C, const double* g, const double* H,
                             const double* v, double k, double kd, double h, ImplicitBlock<N>& out)
{
	const int nHess = N * (N + 1) / 2;
	for (int a = 0; a < N; ++a) { out.f[a] = 0.; out.b[a] = 0.; }
	for (int i = 0; i < nHess; ++i) out.A[i] = 0.;
	const double alpha = h * kd + h * h * k;   // ∇C∇Cᵀ 의 계수(성분과 무관)

	for (int c = 0; c < nComp; ++c) {
		const double* gc = g + c * N;
		const double* Hc = H + c * nHess;
		double Cdot = 0.;
		for (int a = 0; a < N; ++a) Cdot += gc[a] * v[a];
		double Hv[N];   // ∇²C v, 위 삼각 저장에서 대칭으로
		for (int a = 0; a < N; ++a) Hv[a] = 0.;
		int i = 0;
		for (int a = 0; a < N; ++a) {
			Hv[a] += Hc[i] * v[a];
			++i;
			for (int j = a + 1; j < N; ++j, ++i) {
				Hv[a] += Hc[i] * v[j];
				Hv[j] += Hc[i] * v[a];
			}
		}
		const double s = k * C[c] + kd * Cdot;   // 힘 = −s ∇C
		const double beta = h * h * s;           // ∇²C 의 계수
		for (int a = 0; a < N; ++a) {
			out.f[a] += -s * gc[a];
			out.b[a] += -h * s * gc[a] - h * h * (k * Cdot * gc[a] + s * Hv[a]);
		}
		i = 0;
		for (int a = 0; a < N; ++a) {
			const double ga = alpha * gc[a];
			for (int j = a; j < N; ++j, ++i) out.A[i] += ga * gc[j] + beta * Hc[i];
		}
	}
}

} // namespace sx
