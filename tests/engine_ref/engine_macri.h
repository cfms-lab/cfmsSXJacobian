//---------------------------------------------------------------------------
// engine_macri.h — 엔진 Macri 굽힘 손 코드(CheckBendingChange_Macri_SOA)를 부르는 창구. docs/GATE_ENGINE.md.
// 구현은 engine_macri.cpp(엔진 코드를 바이트 그대로 #include 하는 번역 단위).
//---------------------------------------------------------------------------
#pragma once

namespace engine_ref {

struct Macri;

// 힌지 nHinge 개(hinges[4s .. 4s+3] = 엔진 꼭짓점 H, I, J, K = 명세 h, i, j, k), 노드 좌표 X(3 × nNode) · 속도 V.
// X · V 는 호출하는 쪽이 들고 있다(엔진도 노드 배열을 가리키기만 한다). sbv 배열을 모두 할당하고 0 으로 채운다.
Macri* Create(int nHinge, const int* hinges, int nNode, double* X, double* V);
void Destroy(Macri* m);

// CheckBendingChange_Macri_SOA() 한 번 — _Input_b012 부터 ∇²C 까지, 모든 힌지(_Input_F_and_J 는 빈 함수)
void Run(Macri* m);

int Size(const Macri* m);
double Theta(const Macri* m, int i);                 // sbv._C[i] = atan2(sin, cos), 쉬는 각 없음
void Grad(const Macri* m, int i, double g[12]);      // g[3m + s] = sbv._dC_dxm[m][s][i]
void Hess(const Macri* m, int i, double H[144]);     // H[(3m + s) · 12 + 3n + t] = sbv._d2C_dxmdxn[m][n][s][t][i]
bool SinBigger(const Macri* m, int i);               // sbv.pbSinBigger[i] — |sin| > |cos| 가지
double Bytes(const Macri* m);                        // 하네스가 할당한 sbv 배열 바이트
double Checksum(const Macri* m);                     // 시간 측정 뒤 결과를 읽어 최적화로 사라지지 않게

// 판정 밖 진단(tests/diag_engine_terms.cpp)용 중간 배열 — 배치는 Grad · Hess 와 같다
struct Raw {
	double cos, sin;                 // sbv._costh · _sinth (단위 법선)
	double dcos[12], dsin[12];       // sbv._dcosth_dxm · _dsinth_dxm
	double d2cos[144], d2sin[144];   // sbv._d2costh_dxmdxn · _d2sinth_dxmdxn
};
void GetRaw(const Macri* m, int i, Raw& r);

// ---- Bridson(Oh) 경로 — 엔진 기본 굽힘, docs/GATE_BRIDSON.md ----------------------------
struct Bridson;

// 힌지 · 노드는 Create 와 같다. length0[s] = 힌지 변의 Length0(쉬는 변 길이). 재질 BE_U = 1, Theta0 = 0,
// GlobalBendingFactor = 1, 감쇠 꺼짐, bUseUpperJacobianMethod = bUseUDJSerialJacobian = true(엔진 기본값).
// 노드쌍 슬롯(piBendingJ)은 엔진처럼 힌지의 여섯 쌍을 공유 슬롯에 잇는다.
Bridson* CreateBridson(int nHinge, const int* hinges, int nNode, double* X, double* V, const double* length0);
void DestroyBridson(Bridson* b);
void SetSerialBridson(Bridson* b, bool serial);       // bUseUDJSerialJacobian — false 면 쌍 슬롯 환원 없음
void SetStiffnessBridson(Bridson* b, double k);      // 재질 BE_U(기본 1) — 시간 적분 비교(docs/GATE_INTEGRATION.md)

// __CalculateBending 의 Bridson 갈래: bUseC2 = 0 → CheckBendingChange_Macri_SOA() → GetImpBending_Oh_SOA()
void RunBridson(Bridson* b);

int SizeBridson(const Bridson* b);
int PairsBridson(const Bridson* b);                   // 노드쌍 슬롯 수
double ThetaBridson(const Bridson* b, int i);         // sbv._C[i]
bool SkippedBridson(const Bridson* b, int i);         // 퇴화 가드로 건너뜀(s_BendSkip)
double WeightBridson(const Bridson* b, int i);        // 엔진이 쓴 _cb = 3 Length0² / (0.5/_nAim + 0.5/_nBim) (같은 식으로 다시 계산)
void ForceBridson(const Bridson* b, int i, double f[12]);       // pBendingF, 명세 순서 (h, i, j, k)
void JacobianBridson(const Bridson* b, int i, double J[144]);   // pBendingJ 블록 넷 + 전치 + 대각 = −행 합, 명세 순서, 행 우선
double BytesBridson(const Bridson* b);
double ChecksumBridson(const Bridson* b);

} // namespace engine_ref
