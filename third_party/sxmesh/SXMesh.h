//---------------------------------------------------------------------------
// SXMesh.h — 변수 이름을 실행 중에 정해 쓰는 메시(VC++ · C++11 판).
//
// 원본: cfmsCIPC_dev/drape/Sdll/Headers/Mesh/SXMesh.h (C++Builder). 2026-11-12 점검에서 나온
// 치명 · 높음 항목을 고친 판이다. 점검 노트:
//   Projects/Drape/현황과 기록/2026-11-12 SXMesh 코드 점검 — 이름으로 꺼내 쓰는 메시.md
//
// 바뀐 점(요약)
//  - 변수(SXNamedContainer)가 자기 저장소를 소유한다 → SXMesh 복사 · 대입 · vector<SXMesh> 가 안전하다.
//  - 변수 목록은 std::deque → 새 변수를 만들어도 앞서 받은 변수 참조가 무효가 되지 않는다.
//  - 같은 이름으로 다시 만들면 그 변수를 다시 잡는다(중복 변수가 생기지 않는다).
//  - 없는 이름: Find*V 는 NULL, Get*V 는 std::out_of_range, GetN* 는 -1.
//  - NMF 저장 · 불러오기는 표준 스트림. 저장은 UTF-8, 값은 왕복 정밀도(최대 유효숫자 17).
//    불러오기는 옛 C++Builder 판 파일도 읽는다(bData 키, CP949 문자열, 한 줄 1열 데이터, CRLF).
//  - 실패하면 false 와 오류 문구를 돌려주고 메시는 그대로 둔다(크래시 · 반쯤 읽힌 상태 없음).
//
// 그리기(옛 GLDraw · _BindMTLData · S_IO::ShowNormalVector)는 SXMeshGL.h 에 따로 있다.
// 이 파일과 SXMesh.cpp 는 OpenGL 없이 컴파일된다.
// 아직 옮기지 않은 것: SBaseClass(그래프 노드) 연결, 텍스처 이미지 읽기, 재질 경로의 상대 경로 변환.
//---------------------------------------------------------------------------
#ifndef SXMeshH
#define SXMeshH

#include <cstddef>
#include <deque>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

#include "SXNamedContainer.h"

// 옛 enumSGLDisplayOptionPropertyBit(Preheaders.h) 와 같은 값
enum SXDisplayBit {
	dopbVisible           = 1 << 0,
	dopbUseDisplayList    = 1 << 1,
	dopbGouraudNormal     = 1 << 2,
	dopbShowAsSolid       = 1 << 3,
	dopbShowAsWireframe   = 1 << 4,
	dopbShowAsSelected    = 1 << 5,
	dopbShowNormal        = 1 << 6,
	dopbShowGouraudNormal = 1 << 7,
	dopbShowAABB          = 1 << 8,
	dopbShowTexture       = 1 << 9,
	dopbColorMaterial     = 1 << 10,
	dopbUseMipMap         = 1 << 11,
	dopbBackFaceCulling   = 1 << 12,
	dopbShadowDonor       = 1 << 13,
	dopbShadowAcceptor    = 1 << 14
};

typedef std::pair<std::wstring, std::wstring> SXKeyValue;

// 표시 옵션(옛 SGLDisplayOption 중 NMF 에 저장되던 부분)
struct SXMeshDisplay {
	unsigned long PropertyBit;
	float PointSize;
	float LineWidth;
	float GLAlpha;

	SXMeshDisplay() : PropertyBit(0), PointSize(1.f), LineWidth(1.f), GLAlpha(1.f) {}

	bool GetPropertyBit(unsigned long bit) const { return (PropertyBit & bit) != 0; }
	void SetPropertyBit(unsigned long bit, bool on) { if (on) PropertyBit |= bit; else PropertyBit &= ~bit; }
	void TurnOnPropertyBit(unsigned long bit)  { PropertyBit |= bit; }
	void TurnOffPropertyBit(unsigned long bit) { PropertyBit &= ~bit; }
};

// 재질(옛 SMTLData 중 파일에 들어가는 부분 + 그리기용 GL 텍스처 번호). 이미지 읽기는 아직 없다.
struct SXMeshMaterial {
	int ID;
	std::wstring Name;
	std::wstring map_Kd, map_Ks, map_d, bump;   // 경로는 받은 그대로 저장한다
	float Ka[4], Kd[4], Ks[4];
	float Shininess;
	float GLAlpha;
	std::vector<SXKeyValue> Extra;              // 모르는 키는 그대로 보존

	// 그리기용 GL 텍스처 이름(실행 중에만 씀, 파일에 저장하지 않음). 0 = 없음.
	// SXMesh 는 텍스처를 만들거나 지우지 않는다 — 만든 쪽이 지운다. 메시를 복사하면 같은 이름을 함께 쓴다.
	unsigned int tex_Kd, tex_Ks, tex_d, tex_bump;

	SXMeshMaterial();
};

// 이 클래스가 해석하지 않는 섹션(예: [ManikinOBB]). 불러온 그대로 보존하고 다시 저장한다.
struct SXIniSection {
	std::wstring Title;
	std::vector<SXKeyValue> Pairs;
};

class SXMesh {
public:
	SXMesh();
	// 복사 · 대입 · 이동 · 소멸: 컴파일러 기본(모든 멤버가 값 의미) → 깊은 복사

	std::wstring   Name;          // 파일 이름(파일에는 저장하지 않는다)
	SXMeshDisplay  Display;
	std::vector<SXMeshMaterial> Materials;
	std::vector<SXIniSection>   ExtraSections;

	// ---- 이름 붙은 변수 ------------------------------------------------------
	// 새로 만들거나, 같은 이름이 있으면 그 변수를 nRow × nCol 로 다시 잡는다.
	// 이름이 규칙에 맞지 않거나 반대쪽(노드↔요소)에 같은 이름이 있으면 std::invalid_argument.
	SXNodeNamedContainer&    NewNodeV(   const std::wstring& name, int nRow, int nCol);
	SXElementNamedContainer& NewElementV(const std::wstring& name, int nRow, int nCol);

	SXNodeNamedContainer*          FindNodeV(const std::wstring& name);        // 없으면 NULL
	const SXNodeNamedContainer*    FindNodeV(const std::wstring& name) const;
	SXElementNamedContainer*       FindElementV(const std::wstring& name);
	const SXElementNamedContainer* FindElementV(const std::wstring& name) const;

	SXNodeNamedContainer&          GetNodeV(const std::wstring& name);         // 없으면 std::out_of_range
	const SXNodeNamedContainer&    GetNodeV(const std::wstring& name) const;
	SXElementNamedContainer&       GetElementV(const std::wstring& name);
	const SXElementNamedContainer& GetElementV(const std::wstring& name) const;

	int GetNNode(   const std::wstring& name) const;   // 행 개수, 없으면 -1 (옛 API 와 같음)
	int GetNElement(const std::wstring& name) const;

	bool RemoveNodeV(   const std::wstring& name);     // 지우면 다른 변수 참조도 무효가 될 수 있다
	bool RemoveElementV(const std::wstring& name);

	const std::deque<SXNodeNamedContainer>&    NodeVs()    const { return NodeVs_; }
	const std::deque<SXElementNamedContainer>& ElementVs() const { return ElementVs_; }

	void Reset();   // 변수 · 재질 · 보존 섹션 · 표시 옵션 · Name 을 모두 비운다

	// 변수 이름 규칙: 비어 있지 않고, 앞뒤 공백 · 제어 문자 · [ ] = , / 가 없고,
	// "Info" 가 아니며 "pMTLData(" 로 시작하지 않는다(NMF 섹션 이름과 겹치지 않게).
	static bool IsValidVariableName(const std::wstring& name);

	// ---- NMF 입출력 ----------------------------------------------------------
	bool SaveToStream(std::ostream& os, std::wstring* err = 0) const;
	bool SaveToFile(const std::wstring& path, std::wstring* err = 0) const;

	// 불러오기는 내용을 통째로 바꾼다(Name 은 그대로). 실패하면 false, 메시는 그대로.
	bool LoadFromBuffer(const char* buffer, std::size_t size, std::wstring* err = 0);
	bool LoadFromStream(std::istream& is, std::wstring* err = 0);
	bool LoadFromFile(const std::wstring& path, std::wstring* err = 0);

private:
	void CheckNewName(const std::wstring& name, bool forNode) const;

	std::deque<SXNodeNamedContainer>    NodeVs_;
	std::deque<SXElementNamedContainer> ElementVs_;
};

typedef std::vector<SXMesh> SXMeshVector;   // 이제 안전하다

// 옛 매크로 호환: 변수를 만들고 참조를 P_NAME 으로 받는다. P_STR 은 "ppNodes" 또는 L"ppNodes".
// 옛 판과 달리 다음 매크로를 불러도 앞서 받은 참조가 무효가 되지 않는다.
#define _SXMESH_NEW_NODEV(pXMESH, P_NAME, P_STR, P_I, P_J) \
	SXNodeNamedContainer& P_NAME = (pXMESH)->NewNodeV(std::wstring(L"" P_STR), (P_I), (P_J));
#define _SXMESH_NEW_ELEMENTV(pXMESH, P_NAME, P_STR, P_I, P_J) \
	SXElementNamedContainer& P_NAME = (pXMESH)->NewElementV(std::wstring(L"" P_STR), (P_I), (P_J));

#endif
