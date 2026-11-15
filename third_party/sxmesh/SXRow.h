//---------------------------------------------------------------------------
// SXRow.h — 이름 붙은 변수의 한 행(노드 좌표 하나, 요소 하나)을 가리키는 보기(view).
//
// 옛 SXVector3D · SXIntList 를 대신한다(VC++ · C++11).
//  - 저장소는 SXNamedContainer 가 std::vector 로 소유한다. SXRow 는 아무것도 지우지 않는다.
//  - 복사 생성 = 같은 곳을 가리키는 보기 복사. 대입 = 값 복사(옛 SXVector3D::operator= 와 같은 뜻).
//  - 새 값을 만들어 돌려주는 연산(+, -, 외적)은 두지 않는다. 옛 코드는 저장소 없는 임시 객체에
//    값을 써서 접근 위반이 났다. 제자리 연산(+=, -=, *=, /=)만 둔다.
//  - 행을 가리키는 포인터는 그 변수를 다시 Allocate 하거나 메시를 지우면 무효가 된다.
//---------------------------------------------------------------------------
#ifndef SXRowH
#define SXRowH

#include <cassert>

template <class T>
class SXRow {
public:
	T*  rpData;   // 행의 첫 값
	int npData;   // 이 행의 값 개수

	SXRow() : rpData(0), npData(0) {}
	SXRow(T* p, int n) : rpData(p), npData(n) {}
	SXRow(const SXRow& o) : rpData(o.rpData), npData(o.npData) {}

	// 값 복사: 두 행에 겹치는 개수만큼 복사한다. 행 길이(npData)는 바꾸지 않는다.
	SXRow& operator=(const SXRow& o)
	{
		if (rpData != o.rpData) {
			const int n = (npData < o.npData) ? npData : o.npData;
			for (int j = 0; j < n; ++j) rpData[j] = o.rpData[j];
		}
		return *this;
	}

	T& operator[](int j) const
	{
		assert(rpData != 0 && j >= 0 && j < npData);
		return rpData[j];
	}

	int  size()  const { return npData; }
	bool empty() const { return npData == 0; }
	T*   begin() const { return rpData; }
	T*   end()   const { return rpData + npData; }

	// 제자리 연산
	SXRow& operator+=(const SXRow& o) { const int n = Common(o); for (int j = 0; j < n; ++j) rpData[j] += o.rpData[j]; return *this; }
	SXRow& operator-=(const SXRow& o) { const int n = Common(o); for (int j = 0; j < n; ++j) rpData[j] -= o.rpData[j]; return *this; }
	SXRow& operator*=(T s) { for (int j = 0; j < npData; ++j) rpData[j] *= s; return *this; }
	SXRow& operator/=(T s) { for (int j = 0; j < npData; ++j) rpData[j] /= s; return *this; }

	void Set(T x, T y, T z)
	{
		assert(npData >= 3);
		rpData[0] = x; rpData[1] = y; rpData[2] = z;
	}

private:
	int Common(const SXRow& o) const { return (npData < o.npData) ? npData : o.npData; }
};

// 옛 이름 호환. 옛 SXVector3D 의 산술 연산은 없다(SXRow 머리말 참고).
typedef SXRow<double> SXVector3D;
typedef SXRow<int>    SXIntList;

#endif
