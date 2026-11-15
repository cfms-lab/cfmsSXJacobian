//---------------------------------------------------------------------------
// engine_macri.cpp — 엔진 손 코드 비교 게이트의 엔진 쪽 번역 단위(docs/GATE_ENGINE.md · docs/GATE_BRIDSON.md).
//
// generated/engine_ref/*.inc 는 tests/engine_ref/extract_engine.py 가 엔진 cfmsDrape_dev 378d20c 에서 바이트 그대로
// 뽑은 것이다(E1 · B5). 여기서는 그 파일을 고치지 않고 #include 하고, 그 코드가 기대하는 주변만 최소로 채운다:
//   - SDECore 껍데기: 엔진 선언(SDECore_decl.inc)과 이 함수들이 읽는 멤버만
//   - SDEVector3D · SDEMatrix3f · SDEHalfEdge · 재질 · 천 · 풀이기의 최소형(엔진과 같은 멤버 이름 · 같은 연산)
//   - sbv 할당: 엔진 SetMem 을 옮기지 않고 같은 배치로 새로 짰다(배열마다 new, _TEMP333 = _TEMP[10 …] 별칭, 0 으로 시작)
//   - _Input_F_and_J: 빈 함수(Baraff 경로의 조립은 옮기지 않는다)
// 두 경로: Macri(Baraff, bUseC2 = 1) 는 CheckBendingChange_Macri_SOA 만, Bridson(엔진 기본값)은 __CalculateBending 의
// Bridson 갈래 그대로 bUseC2 = 0 → CheckBendingChange_Macri_SOA(θ 만) → GetImpBending_Oh_SOA(힘 · Jacobian 블록).
// 엔진 코드 번역 단위라 이 리포의 /W4 /WX 를 걸지 않는다(run_gate.bat). 최적화 스위치는 같다(/O2).
//---------------------------------------------------------------------------
#define _USE_MATH_DEFINES   // 엔진 SDEPreHeaders.h:142 — Oh 경로가 M_PI 를 쓴다
#include "engine_macri.h"

#include <math.h>
#include <string.h>
#include <unordered_map>
#include <vector>

#include "SDEPreHeaders.h"
#include "../../generated/engine_ref/SDEmysse.h"

// 엔진 SDEMatrix3f(_USE_DYNAMIC_SMatrix3f 꺼짐 — 엔진 기본)에서 Oh 경로가 부르는 것만. 정의는 엔진 것(SDEMatrix3f_oh.inc)
class SDEMatrix3f {
public:
	dFloat Data[3][3];
	void AXPY(double a, double x[3][3], double b);
	void Set(double* d, double alpha);
};
#include "../../generated/engine_ref/SDEMatrix3f_oh.inc"

// 엔진 SDEVector3D 와 같은 꼴(x, y, z 와 Data[3] 공용체, 성분마다 같은 연산). Macri 경로는 꺼진 if(0) 블록에서만 쓴다
struct SDEVector3D {
	union {
		struct { dFloat x, y, z; };
		dFloat Data[3];
	};
	SDEVector3D(dFloat a = 0, dFloat b = 0, dFloat c = 0) { x = a; y = b; z = c; }
	explicit SDEVector3D(const dFloat* p) { x = p[0]; y = p[1]; z = p[2]; }
	void Set(dFloat a, dFloat b, dFloat c) { x = a; y = b; z = c; }
	void operator+=(const SDEVector3D& a) { x += a.x; y += a.y; z += a.z; }
	SDEVector3D operator*(dFloat s) const { SDEVector3D t; t.Set(x * s, y * s, z * s); return t; }
	SDEVector3D operator+(const SDEVector3D& o) const { return SDEVector3D(x + o.x, y + o.y, z + o.z); }
	SDEVector3D operator-(const SDEVector3D& o) const { return SDEVector3D(x - o.x, y - o.y, z - o.z); }
	void Normalize() { const dFloat s = sqrt(x * x + y * y + z * z); if (s > 0.) { x /= s; y /= s; z /= s; } }
};
inline SDEVector3D operator*(dFloat s, const SDEVector3D& a) { SDEVector3D t; t.Set(a.x * s, a.y * s, a.z * s); return t; }

struct SDEMaterialPropertyShim { dFloat Theta0, BE_U, BV; };
struct SDEHalfEdge { SDEMaterialPropertyShim* pMaterialProperty; dFloat Length0; };
struct SDESewingShim { int nStitch; };
struct SDEClothShim { dFloat** ppNodes; dFloat* p2DNodeMass; SDESewingShim* pSewing; };
struct SDEPCGSolverShim {   // 엔진 SDEPCGSolver 에서 Oh 경로가 읽고 쓰는 멤버만
	bool bUseTransposeAwareUpperJacobian, bUseUpperJacobianMethod, bUseUDJSerialJacobian;
	dFloat* pBendingF;
	SDEMatrix3f *pBendingJ, *pBendingV, *pSerialJacobianJ, *pSerialJacobianV;
	int *piBendingJ, *piBendingV;
};

#include "../../generated/engine_ref/SDEImplicitVariables_sbv.inc"
#include "../../generated/engine_ref/SDECore_globals.inc"

class SDECore {
public:
#include "../../generated/engine_ref/SDECore_decl.inc"

	_SerialBendingVariables sbv;
	int I_Bend_0;
	dFloat bUseC2;
	SDEImplicitBendingForceType iImplicitBendingForceType;
	dFloat DimensionFactor;
	SDEHalfEdge* _pBendingSpring;
	SDEClothShim* pCloth;
	dFloatVec _v0;
	// Oh 경로(GetImpBending_Oh_SOA)가 읽는 것
	SDEPCGSolverShim* pCGSolver;
	bool bCalculateSewing, bSeamClosureMeasured, bBendDamping, bBendDampingImplicit;
	dFloat SeamClosureRatio, GlobalBendingFactor, _h, BendDampCapFactor;
	dFloatVec _f0;
	int** BMPointer_Oh;
	dFloat** _df_dx_2d;
};

#include "../../generated/engine_ref/SDECore_macri.inc"
#include "../../generated/engine_ref/SDECore_oh.inc"

// 조립(힘을 _f0 에, Jacobian 을 풀이기 구조에)은 옮기지 않는다 — 이 게이트는 미분까지만 견준다
void SDECore::_Input_F_and_J(int, int) {}

// 엔진 __PrepareImplicit 의 _DO_NOT_USE_SSE 경로 대입문(바이트 그대로)
static void AssignBestSSE()
{
#include "../../generated/engine_ref/SDECore_assign.inc"
}

// ---- sbv 할당 (엔진 SetMem 과 같은 배치, 새로 짬) ------------------------------
static double g_bytes = 0.;

static void Alloc(dFloat*& p, int n) { p = new dFloat[n](); g_bytes += 8. * n; }

_SerialBendingVariables::_SerialBendingVariables() { ResetPointers(); }
_SerialBendingVariables::~_SerialBendingVariables() {}   // 하네스는 프로세스 끝까지 쓴다(Destroy 가 지운다)
void _SerialBendingVariables::Reset(void) {}
void _SerialBendingVariables::ResetPointers(void) { memset(this, 0, sizeof *this); }

void _SerialBendingVariables::SetMem(int size)
{
	const int n = size;
	nData = n;
	bSSESupported = false;
	_vid = new int[4 * n]();
	g_bytes += 16. * n;
	for (int i = 0; i < NUM_TEMP; ++i) { Alloc(_TEMP[i], n); b_TEMP[i] = false; }
	int r = 10;   // 엔진 SetMem 과 같은 별칭
	for (int s = 0; s < 3; ++s) for (int t = 0; t < 3; ++t) for (int q = 0; q < 3; ++q) _TEMP333[s][t][q] = _TEMP[r++];
	for (int m = 0; m < 4; ++m) for (int s = 0; s < 3; ++s) {
		Alloc(_qA[m][s], n); Alloc(_qB[m][s], n);
		Alloc(_x0[m][s], n); Alloc(_v0[m][s], n);
		Alloc(_dC_dxm[m][s], n);
		Alloc(_dcosth_dxm[m][s], n); Alloc(_dsinth_dxm[m][s], n);
		for (int t = 0; t < 3; ++t) {
			Alloc(_dnA_dxm[m][s][t], n); Alloc(_dnB_dxm[m][s][t], n);
			for (int p = 0; p < 3; ++p) { Alloc(_dnA_tilde_dxm[m][s][t][p], n); Alloc(_dnB_tilde_dxm[m][s][t][p], n); }
		}
	}
	for (int s = 0; s < 3; ++s) {
		Alloc(_nhatA[s], n); Alloc(_nA[s], n); Alloc(_nhatB[s], n); Alloc(_nB[s], n); Alloc(_ehat[s], n);
		Alloc(N12[s], n); Alloc(_nA_cross_nB[s], n);
		for (int t = 0; t < 3; ++t) {
			Alloc(_nA_tilde[s][t], n); Alloc(_nB_tilde[s][t], n); Alloc(_e_tilde[s][t], n);
			Alloc(_I_tildetilde_nA_e[s][t], n); Alloc(_I_tildetilde_nB_e[s][t], n);
		}
	}
	Alloc(_nAim, n); Alloc(_nBim, n); Alloc(_eim, n); Alloc(_costh, n); Alloc(_sinth, n); Alloc(_C, n);
	Alloc(_k_bend, n); Alloc(_k_bend_damp, n); Alloc(_dC_dt, n);
	Alloc(dot0332, n); Alloc(dot1332, n); Alloc(dot2032, n); Alloc(dot2132, n); Alloc(x32_dot, n); Alloc(x32_size, n); Alloc(_cb, n);
	for (int m = 0; m < 4; ++m) for (int nn = 0; nn < 4; ++nn) for (int s = 0; s < 3; ++s) for (int t = 0; t < 3; ++t) {
		Alloc(_d2costh_dxmdxn[m][nn][s][t], n); Alloc(_d2sinth_dxmdxn[m][nn][s][t], n); Alloc(_d2C_dxmdxn[m][nn][s][t], n);
	}
	pbSinBigger = new bool[n]();
	g_bytes += n;
	// 쓰지 않는 것: _d2nA_dxmdxn · _d2nB_dxmdxn(엔진은 864 개 배열을 할당하지만 이 함수들이 건드리지 않는다), _Jacobian, _f
}

// ---- 창구 -------------------------------------------------------------------
namespace engine_ref {

struct Macri {
	SDECore core;
	SDEClothShim cloth;
	std::vector<dFloat*> nodes;
	SDEMaterialPropertyShim mat;
	std::vector<SDEHalfEdge> edges;
	double bytes;
};

// 두 경로가 같이 쓰는 준비: 노드 배열 · 재질(Theta0 = 0, BE_U = 1) · 반변 · sbv. length0 이 0 이면 Length0 = 0
static void SetupCore(SDECore& c, SDEClothShim& cloth, std::vector<dFloat*>& nodes, SDEMaterialPropertyShim& mat,
                      std::vector<SDEHalfEdge>& edges, int nHinge, const int* hinges, int nNode, double* X, double* V,
                      const double* length0, double& bytes)
{
	static bool assigned = false;
	if (!assigned) { AssignBestSSE(); assigned = true; }
	nodes.resize(nNode);
	for (int i = 0; i < nNode; ++i) nodes[i] = X + 3 * i;
	cloth.ppNodes = nodes.empty() ? 0 : &nodes[0];
	cloth.p2DNodeMass = 0;
	cloth.pSewing = 0;
	mat.Theta0 = 0.;
	mat.BE_U = 1.;
	mat.BV = 0.;
	edges.resize(nHinge);
	for (int i = 0; i < nHinge; ++i) { edges[i].pMaterialProperty = &mat; edges[i].Length0 = length0 ? length0[i] : 0.; }
	c.pCloth = &cloth;
	c._v0 = V;
	c.DimensionFactor = 1.;
	c.I_Bend_0 = nHinge;
	c._pBendingSpring = edges.empty() ? 0 : &edges[0];
	c.pCGSolver = 0;
	c.bCalculateSewing = false;      // 봉제 없음 → _bendDampingScale = 1 (감쇠는 아래에서 꺼진다)
	c.bSeamClosureMeasured = false;
	c.bBendDamping = false;          // 엔진 기본값(SDECore.h bBendDamping)
	c.bBendDampingImplicit = false;
	c.SeamClosureRatio = 0.;
	c.GlobalBendingFactor = 1.;
	c._h = 0.;
	c.BendDampCapFactor = 0.;
	c._f0 = 0;                       // 옛(비 upper) 경로에서만 쓴다
	c.BMPointer_Oh = 0;
	c._df_dx_2d = 0;
	g_bytes = 0.;
	c.sbv.SetMem(nHinge);
	bytes = g_bytes;
	for (int i = 0; i < 4 * nHinge; ++i) c.sbv._vid[i] = hinges[i];
}

Macri* Create(int nHinge, const int* hinges, int nNode, double* X, double* V)
{
	Macri* m = new Macri;
	SetupCore(m->core, m->cloth, m->nodes, m->mat, m->edges, nHinge, hinges, nNode, X, V, 0, m->bytes);
	m->core.bUseC2 = 1.;
	m->core.iImplicitBendingForceType = SDEibfImplicitBaraff;
	return m;
}

void Destroy(Macri* m) { delete m; }   // 배열은 프로세스가 끝날 때 돌려준다(시험 exe)

void Run(Macri* m) { m->core.CheckBendingChange_Macri_SOA(); }

int Size(const Macri* m) { return m->core.I_Bend_0; }
double Theta(const Macri* m, int i) { return m->core.sbv._C[i]; }

void Grad(const Macri* m, int i, double g[12])
{
	for (int a = 0; a < 4; ++a) for (int s = 0; s < 3; ++s) g[3 * a + s] = m->core.sbv._dC_dxm[a][s][i];
}

void Hess(const Macri* m, int i, double H[144])
{
	for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) for (int s = 0; s < 3; ++s) for (int t = 0; t < 3; ++t)
		H[(3 * a + s) * 12 + 3 * b + t] = m->core.sbv._d2C_dxmdxn[a][b][s][t][i];
}

void GetRaw(const Macri* m, int i, Raw& r)
{
	const _SerialBendingVariables& v = m->core.sbv;
	r.cos = v._costh[i];
	r.sin = v._sinth[i];
	for (int a = 0; a < 4; ++a) for (int s = 0; s < 3; ++s) { r.dcos[3 * a + s] = v._dcosth_dxm[a][s][i]; r.dsin[3 * a + s] = v._dsinth_dxm[a][s][i]; }
	for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) for (int s = 0; s < 3; ++s) for (int t = 0; t < 3; ++t) {
		r.d2cos[(3 * a + s) * 12 + 3 * b + t] = v._d2costh_dxmdxn[a][b][s][t][i];
		r.d2sin[(3 * a + s) * 12 + 3 * b + t] = v._d2sinth_dxmdxn[a][b][s][t][i];
	}
}

bool SinBigger(const Macri* m, int i) { return m->core.sbv.pbSinBigger[i]; }
double Bytes(const Macri* m) { return m->bytes; }

double Checksum(const Macri* m)
{
	double s = 0.;
	const int n = m->core.I_Bend_0;
	for (int i = 0; i < n; i += 97) s += m->core.sbv._C[i] + m->core.sbv._dC_dxm[1][2][i] + m->core.sbv._d2C_dxmdxn[1][2][0][1][i];
	return s;
}

// ---- Bridson(Oh) 경로 — docs/GATE_BRIDSON.md -----------------------------------
struct Bridson {
	SDECore core;
	SDEClothShim cloth;
	std::vector<dFloat*> nodes;
	SDEMaterialPropertyShim mat;
	std::vector<SDEHalfEdge> edges;
	SDEPCGSolverShim solver;
	std::vector<dFloat> F;                 // pBendingF  [12 × 힌지]
	std::vector<SDEMatrix3f> J, SJ;        // pBendingJ [6 × 힌지](힌지 소유), pSerialJacobianJ [노드쌍]
	std::vector<int> iJ;                   // piBendingJ [6 × 힌지] → 노드쌍 슬롯
	int nPair;
	double bytes;
};

// Oh 함수의 지역 번호 nn = (vid0, vid3, vid2, vid1) — 지역 번호 → 명세 꼭짓점 (h, i, j, k) = (0, 1, 2, 3)
static const int kOhToSpec[4] = { 0, 3, 2, 1 };

Bridson* CreateBridson(int nHinge, const int* hinges, int nNode, double* X, double* V, const double* length0)
{
	Bridson* b = new Bridson;
	SDECore& c = b->core;
	SetupCore(c, b->cloth, b->nodes, b->mat, b->edges, nHinge, hinges, nNode, X, V, length0, b->bytes);
	c.bUseC2 = 0.;
	c.iImplicitBendingForceType = SDEibfImplicitBridson;
	b->F.assign(12 * static_cast<std::size_t>(nHinge), 0.);
	SDEMatrix3f zero;
	memset(&zero, 0, sizeof zero);
	b->J.assign(6 * static_cast<std::size_t>(nHinge), zero);
	// 노드쌍 슬롯: 엔진처럼 힌지의 여섯 쌍(지역 번호 (0,1)(0,2)(0,3)(1,2)(1,3)(2,3))을 공유 슬롯에 잇는다
	static const int prs[6][2] = { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} };
	std::unordered_map<long long, int> slot;
	b->iJ.resize(6 * static_cast<std::size_t>(nHinge));
	for (int i = 0; i < nHinge; ++i) {
		const int nn[4] = { hinges[4 * i], hinges[4 * i + 3], hinges[4 * i + 2], hinges[4 * i + 1] };
		for (int k = 0; k < 6; ++k) {
			int p = nn[prs[k][0]], q = nn[prs[k][1]];
			if (p > q) { const int t = p; p = q; q = t; }
			const long long key = static_cast<long long>(p) * nNode + q;
			std::unordered_map<long long, int>::const_iterator it = slot.find(key);
			const int s = it == slot.end() ? static_cast<int>(slot.size()) : it->second;
			if (it == slot.end()) slot[key] = s;
			b->iJ[6 * i + k] = s;
		}
	}
	b->nPair = static_cast<int>(slot.size());
	b->SJ.assign(b->nPair, zero);
	b->bytes += 8. * b->F.size() + sizeof(SDEMatrix3f) * (b->J.size() + b->SJ.size()) + 4. * b->iJ.size();
	SDEPCGSolverShim& s = b->solver;
	s.bUseTransposeAwareUpperJacobian = false;
	s.bUseUpperJacobianMethod = true;     // 엔진 기본값
	s.bUseUDJSerialJacobian = true;       // 엔진 기본값(SetSerialBridson 로 바꾼다)
	s.pBendingF = b->F.empty() ? 0 : &b->F[0];
	s.pBendingJ = b->J.empty() ? 0 : &b->J[0];
	s.pSerialJacobianJ = b->SJ.empty() ? 0 : &b->SJ[0];
	s.piBendingJ = b->iJ.empty() ? 0 : &b->iJ[0];
	s.pBendingV = 0;                      // 감쇠가 꺼져 있어 쓰지 않는다
	s.pSerialJacobianV = 0;
	s.piBendingV = 0;
	c.pCGSolver = &s;
	return b;
}

void DestroyBridson(Bridson* b) { delete b; }
void SetSerialBridson(Bridson* b, bool serial) { b->solver.bUseUDJSerialJacobian = serial; }
void SetStiffnessBridson(Bridson* b, double k) { b->mat.BE_U = k; }   // _kb = BE_U × GlobalBendingFactor, 힌지마다 읽는다

void RunBridson(Bridson* b)
{
	SDECore& c = b->core;   // __CalculateBending 의 Bridson 갈래(SDECore.cpp:3577)
	c.bUseC2 = 0;
	c.CheckBendingChange_Macri_SOA();
	c.GetImpBending_Oh_SOA();
}

int SizeBridson(const Bridson* b) { return b->core.I_Bend_0; }
int PairsBridson(const Bridson* b) { return b->nPair; }
double ThetaBridson(const Bridson* b, int i) { return b->core.sbv._C[i]; }
bool SkippedBridson(const Bridson* b, int i) { return s_BendSkip[i] != 0; }
double BytesBridson(const Bridson* b) { return b->bytes; }

double WeightBridson(const Bridson* b, int i)
{
	const _SerialBendingVariables& sbv = b->core.sbv;   // Oh 함수와 같은 식 · 같은 순서
	const dFloat _areaAB = 0.5/sbv._nAim[i] + 0.5/sbv._nBim[i];
	return 3. * b->edges[i].Length0 * b->edges[i].Length0 / _areaAB;
}

void ForceBridson(const Bridson* b, int i, double f[12])
{
	for (int a = 0; a < 4; ++a) for (int s = 0; s < 3; ++s) f[3 * kOhToSpec[a] + s] = b->F[12 * static_cast<std::size_t>(i) + 3 * a + s];
}

void JacobianBridson(const Bridson* b, int i, double J[144])
{
	static const int prs[6][2] = { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} };
	for (int k = 0; k < 144; ++k) J[k] = 0.;
	for (int k = 1; k <= 4; ++k) {   // 엔진이 쓰는 블록: (0,2)(0,3)(1,2)(1,3). (0,1) · (2,3) 은 쓰지 않는다
		const SDEMatrix3f& M = b->J[6 * static_cast<std::size_t>(i) + k];
		const int r = kOhToSpec[prs[k][0]], q = kOhToSpec[prs[k][1]];
		for (int s = 0; s < 3; ++s) for (int t = 0; t < 3; ++t) {
			J[(3 * r + s) * 12 + 3 * q + t] = M.Data[s][t];
			J[(3 * q + t) * 12 + 3 * r + s] = M.Data[s][t];   // 아래 삼각 = 전치
		}
	}
	for (int m = 0; m < 4; ++m)      // 대각 = −행 합(풀이기 규약)
		for (int s = 0; s < 3; ++s) for (int t = 0; t < 3; ++t) {
			double sum = 0.;
			for (int n = 0; n < 4; ++n) if (n != m) sum += J[(3 * m + s) * 12 + 3 * n + t];
			J[(3 * m + s) * 12 + 3 * m + t] = -sum;
		}
}

double ChecksumBridson(const Bridson* b)
{
	double s = 0.;
	const int n = b->core.I_Bend_0;
	for (int i = 0; i < n; i += 97) s += b->F[12 * static_cast<std::size_t>(i) + 7] + b->J[6 * static_cast<std::size_t>(i) + 2].Data[1][2];
	for (int p = 0; p < b->nPair; p += 97) s += b->SJ[p].Data[0][0];
	return s;
}

} // namespace engine_ref
