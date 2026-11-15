//---------------------------------------------------------------------------
// diag_engine_terms.cpp — 판정 밖 진단(docs/GATE_ENGINE.md 「판정 뒤 진단」). 엔진 Macri 굽힘 코드의 중간 배열이
// 무엇의 미분인지 기호 미분과 견준다. E3 가 「어느 쪽도 아님」이라 판정문에 미리 적은 대로 한다.
//
//   P = nA · nB, Q = (nA × nB) · e (엔진의 정규화하지 않은 벡터) — tests/diag_engine_terms_spec.py 를 생성기로 미분
//   (가) 엔진 _dcosth_dxm · _dsinth_dxm 이 ∇P · ∇Q 인가
//   (나) 엔진 _d2costh_dxmdxn · _d2sinth_dxmdxn 이 ∇²P · ∇²Q 인가(그대로 / 블록 안 전치)
//   (다) 엔진 ∇C · ∇²C 가 엔진 자신의 중간 배열에 읽은 가지 공식을 건 것과 같은가(공식을 바르게 읽었나)
//   (라) 엔진 ∇C · ∇²C 가 ∇P · ∇Q · ∇²P · ∇²Q(기호 미분)에 그 공식을 건 것과 같은가 — 같으면 엔진 값이 설명된다
//   (마) 꼭짓점 블록마다 엔진 dcos · dsin 이 곱의 법칙 항(diag_engine_split_spec.py)의 어떤 부호 조합인가
// 실행 디렉터리 = 리포 루트. results/diag_engine_terms.json 을 쓴다.
//---------------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS
#include "../generated/EngineTerms.h"
#include "../generated/EngineTermsSplit.h"
#include "gate_common.h"
#include "engine_ref/engine_macri.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

using namespace gate;

typedef sxgen::EngineTerms T;
typedef sxgen::EngineTermsSplit TS;

static void Cross(const double* a, const double* b, double* c)
{
	c[0] = a[1] * b[2] - a[2] * b[1]; c[1] = a[2] * b[0] - a[0] * b[2]; c[2] = a[0] * b[1] - a[1] * b[0];
}

// (마) 조합: 계수 −1 · 0 · +1. P 는 항 2 개(9 조합), Q 는 항 3 개(27 조합)
static int Coef(int combo, int k) { for (int i = 0; i < k; ++i) combo /= 3; return combo % 3 - 1; }
static const int N = 12;
static const int NH = 78;

static void Full(const double* Hu, double* Hf)
{
	for (int a = 0; a < N; ++a)
		for (int b = 0; b < N; ++b) Hf[a * N + b] = Hu[a <= b ? sx::UpperIndex<N>(a, b) : sx::UpperIndex<N>(b, a)];
}

// 블록 안 전치: [m][n][s][t] → [m][n][t][s]
static void BlockT(const double* H, double* Ht)
{
	for (int m = 0; m < 4; ++m) for (int n = 0; n < 4; ++n) for (int s = 0; s < 3; ++s) for (int t = 0; t < 3; ++t)
		Ht[(3 * m + s) * N + 3 * n + t] = H[(3 * m + t) * N + 3 * n + s];
}

// 엔진 CheckBendingChange_Macri_SOA #b.9 · #b.10 의 가지 공식(SDECore.cpp 7878 – 7975 을 읽은 그대로)
static void BranchFormula(bool sinBigger, double sn, double cs, const double* dcos, const double* dsin,
	const double* d2cos, const double* d2sin, double* g, double* H)
{
	for (int a = 0; a < N; ++a) g[a] = sinBigger ? dcos[a] / -sn : dsin[a] / cs;
	for (int m = 0; m < 4; ++m) for (int n = 0; n < 4; ++n) for (int s = 0; s < 3; ++s) for (int t = 0; t < 3; ++t) {
		const int r = 3 * m + s, c = 3 * n + t;
		if (sinBigger)   // _TEMP[0..8] = Outer(dsin[n], dcos[m])
			H[r * N + c] = 1. / (sn * sn) * (dsin[3 * n + s] * dcos[3 * m + t]) - 1. / sn * d2cos[r * N + c];
		else             // _TEMP[9..17] = Outer(dcos[n], dsin[m])
			H[r * N + c] = -1. / (cs * cs) * (dcos[3 * n + s] * dsin[3 * m + t]) + 1. / cs * d2sin[r * N + c];
	}
}

int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif
	const std::vector<Hinge> hs = ReadHinges("results/hinges.txt");
	if (hs.empty()) { std::printf("results/hinges.txt 없음\n"); return 2; }
	const int nS = static_cast<int>(hs.size());
	std::vector<double> X(12 * nS), V(12 * nS);
	std::vector<int> vid(4 * nS);
	for (int s = 0; s < nS; ++s) {
		for (int i = 0; i < 12; ++i) { X[12 * s + i] = hs[s].x[i]; V[12 * s + i] = hs[s].v[i]; }
		for (int m = 0; m < 4; ++m) vid[4 * s + m] = 4 * s + m;
	}
	engine_ref::Macri* em = engine_ref::Create(nS, &vid[0], 4 * nS, &X[0], &V[0]);
	engine_ref::Run(em);

	Stat dP, dQ, d2P, d2PT, d2Q, d2QT, selfG, selfH, symG, symH, symHT;
	double errP[4][9] = {}, errQ[4][27] = {};   // 블록 · 조합마다 표본 최대 상대오차
	for (int s = 0; s < nS; ++s) {
		engine_ref::Raw r;
		engine_ref::GetRaw(em, s, r);
		const bool sb = engine_ref::SinBigger(em, s);
		double ge[N], He[N * N];
		engine_ref::Grad(em, s, ge);
		engine_ref::Hess(em, s, He);

		double PQ[2], gPQ[2 * N], HPQ[2 * NH], HP[N * N], HQ[N * N];
		T::Eval(hs[s].x, PQ, gPQ, HPQ);
		Full(HPQ, HP);
		Full(HPQ + NH, HQ);
		const double* gP = gPQ;
		const double* gQ = gPQ + N;

		dP.Add(Rel(MaxAbsDiff(r.dcos, gP, N), MaxAbs(gP, N)));
		dQ.Add(Rel(MaxAbsDiff(r.dsin, gQ, N), MaxAbs(gQ, N)));
		double t[N * N];
		d2P.Add(Rel(MaxAbsDiff(r.d2cos, HP, N * N), MaxAbs(HP, N * N)));
		BlockT(r.d2cos, t);
		d2PT.Add(Rel(MaxAbsDiff(t, HP, N * N), MaxAbs(HP, N * N)));
		d2Q.Add(Rel(MaxAbsDiff(r.d2sin, HQ, N * N), MaxAbs(HQ, N * N)));
		BlockT(r.d2sin, t);
		d2QT.Add(Rel(MaxAbsDiff(t, HQ, N * N), MaxAbs(HQ, N * N)));

		double g[N], H[N * N];
		BranchFormula(sb, r.sin, r.cos, r.dcos, r.dsin, r.d2cos, r.d2sin, g, H);
		selfG.Add(Rel(MaxAbsDiff(ge, g, N), MaxAbs(ge, N)));
		selfH.Add(Rel(MaxAbsDiff(He, H, N * N), MaxAbs(He, N * N)));
		// (라) 기호 미분에 공식을 건다. ∇²P · ∇²Q 는 엔진 배열이 그대로 맞으면 그대로, 블록 전치로 맞으면 전치해서
		BranchFormula(sb, r.sin, r.cos, gP, gQ, HP, HQ, g, H);
		symG.Add(Rel(MaxAbsDiff(ge, g, N), MaxAbs(ge, N)));
		symH.Add(Rel(MaxAbsDiff(He, H, N * N), MaxAbs(He, N * N)));
		double HPt[N * N], HQt[N * N];
		BlockT(HP, HPt);
		BlockT(HQ, HQt);
		BranchFormula(sb, r.sin, r.cos, gP, gQ, HPt, HQt, g, H);
		symHT.Add(Rel(MaxAbsDiff(He, H, N * N), MaxAbs(He, N * N)));

		// (마) 엔진 벡터를 고정값으로 넣고 항별 기울기
		const double* x = hs[s].x;
		double a1[3], a2[3], b1[3], b2[3], ev[9];
		for (int c = 0; c < 3; ++c) {
			a1[c] = x[6 + c] - x[3 + c]; a2[c] = x[3 + c] - x[c];          // qA[0] = xJ − xI, qA[2] = xI − xH
			b1[c] = x[6 + c] - x[9 + c]; b2[c] = x[9 + c] - x[3 + c];      // qB[1] = xJ − xK, qB[2] = xK − xI
			ev[6 + c] = x[3 + c] - x[6 + c];                               // e = qB[3] = xI − xJ
		}
		Cross(a1, a2, ev);
		Cross(b1, b2, ev + 3);
		double C5[5], g5[5 * N], H5[5 * NH];
		TS::Eval(x, ev, C5, g5, H5);
		for (int m = 0; m < 4; ++m) {
			for (int combo = 0; combo < 9; ++combo) {
				double v[3], sc = MaxAbs(r.dcos + 3 * m, 3);
				for (int c = 0; c < 3; ++c) {
					v[c] = Coef(combo, 0) * g5[0 * N + 3 * m + c] + Coef(combo, 1) * g5[1 * N + 3 * m + c];
					sc = std::max(sc, std::max(std::fabs(g5[0 * N + 3 * m + c]), std::fabs(g5[1 * N + 3 * m + c])));
				}
				errP[m][combo] = std::max(errP[m][combo], Rel(MaxAbsDiff(r.dcos + 3 * m, v, 3), sc));
			}
			for (int combo = 0; combo < 27; ++combo) {
				double v[3], sc = MaxAbs(r.dsin + 3 * m, 3);
				for (int c = 0; c < 3; ++c) {
					v[c] = 0.;
					for (int k = 0; k < 3; ++k) {
						v[c] += Coef(combo, k) * g5[(2 + k) * N + 3 * m + c];
						sc = std::max(sc, std::fabs(g5[(2 + k) * N + 3 * m + c]));
					}
				}
				errQ[m][combo] = std::max(errQ[m][combo], Rel(MaxAbsDiff(r.dsin + 3 * m, v, 3), sc));
			}
		}
	}
	engine_ref::Destroy(em);

	std::printf("진단 — 엔진 Macri 중간 배열 대 기호 미분 (%d 힌지, 상대오차 max)\n", nS);
	std::printf("  (가) _dcosth_dxm 대 ∇P %.2e,  _dsinth_dxm 대 ∇Q %.2e\n", dP.Max(), dQ.Max());
	std::printf("  (나) _d2costh 대 ∇²P: 그대로 %.2e · 블록 전치 %.2e,  _d2sinth 대 ∇²Q: 그대로 %.2e · 블록 전치 %.2e\n",
		d2P.Max(), d2PT.Max(), d2Q.Max(), d2QT.Max());
	std::printf("  (다) 엔진 중간 배열 + 가지 공식 대 엔진 ∇C %.2e · ∇²C %.2e\n", selfG.Max(), selfH.Max());
	std::printf("  (라) 기호 미분 + 가지 공식 대 엔진 ∇C %.2e · ∇²C %.2e (∇²P · ∇²Q 블록 전치 %.2e)\n",
		symG.Max(), symH.Max(), symHT.Max());
	static const char* vName[4] = { "H", "I", "J", "K" };
	int bestP[4], bestQ[4];
	std::printf("  (마) 꼭짓점 블록마다 가장 맞는 부호 조합(표본 최대 상대오차) — dcos = a·∇PA + b·∇PB, dsin = a·∇QA + b·∇QB + c·∇QE\n");
	for (int m = 0; m < 4; ++m) {
		bestP[m] = 0; bestQ[m] = 0;
		for (int k = 1; k < 9; ++k) if (errP[m][k] < errP[m][bestP[m]]) bestP[m] = k;
		for (int k = 1; k < 27; ++k) if (errQ[m][k] < errQ[m][bestQ[m]]) bestQ[m] = k;
		std::printf("      %s: dcos (a, b) = (%+d, %+d) %.1e  [정답 (+1, +1) %.1e],  dsin (a, b, c) = (%+d, %+d, %+d) %.1e  [정답 (+1, +1, +1) %.1e]\n",
			vName[m], Coef(bestP[m], 0), Coef(bestP[m], 1), errP[m][bestP[m]], errP[m][8],
			Coef(bestQ[m], 0), Coef(bestQ[m], 1), Coef(bestQ[m], 2), errQ[m][bestQ[m]], errQ[m][26]);
	}
	std::FILE* js = std::fopen("results/diag_engine_terms.json", "w");
	if (js) {
		std::fprintf(js, "{\n  \"engine_commit\": \"378d20c\",\n  \"samples\": %d,\n", nS);
		std::fprintf(js, "  \"dcos_vs_gradP\": %.3e, \"dsin_vs_gradQ\": %.3e,\n", dP.Max(), dQ.Max());
		std::fprintf(js, "  \"d2cos_vs_hessP\": %.3e, \"d2cos_blockT_vs_hessP\": %.3e, \"d2sin_vs_hessQ\": %.3e, \"d2sin_blockT_vs_hessQ\": %.3e,\n",
			d2P.Max(), d2PT.Max(), d2Q.Max(), d2QT.Max());
		std::fprintf(js, "  \"formula_on_engine_terms\": {\"grad\": %.3e, \"hess\": %.3e},\n", selfG.Max(), selfH.Max());
		std::fprintf(js, "  \"formula_on_symbolic\": {\"grad\": %.3e, \"hess\": %.3e, \"hess_with_blockT_second\": %.3e},\n",
			symG.Max(), symH.Max(), symHT.Max());
		std::fprintf(js, "  \"best_sign_combo_by_vertex\": {");
		for (int m = 0; m < 4; ++m)
			std::fprintf(js, "%s\"%s\": {\"dcos_PA_PB\": [%d, %d], \"dcos_err\": %.3e, \"dcos_err_if_correct\": %.3e, \"dsin_QA_QB_QE\": [%d, %d, %d], \"dsin_err\": %.3e, \"dsin_err_if_correct\": %.3e}",
				m ? ", " : "", vName[m], Coef(bestP[m], 0), Coef(bestP[m], 1), errP[m][bestP[m]], errP[m][8],
				Coef(bestQ[m], 0), Coef(bestQ[m], 1), Coef(bestQ[m], 2), errQ[m][bestQ[m]], errQ[m][26]);
		std::fprintf(js, "}\n}\n");
		std::fclose(js);
	}
	return 0;
}
