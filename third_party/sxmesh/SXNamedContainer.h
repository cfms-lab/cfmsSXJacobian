//---------------------------------------------------------------------------
// SXNamedContainer.h — 이름 붙은 변수 하나(ppNodes, pMTLID ...)를 담는 컨테이너.
//
// 옛 SXNodeNamedContainer · SXElementNamedContainer 를 대신한다(VC++ · C++11).
//  - 값은 nRow × nCol 연속 블록 하나(std::vector<T>)에 담는다. 컨테이너가 블록을 소유하므로
//    복사 · 대입 · 이동 · 소멸은 컴파일러 기본 동작으로 깊은 복사와 해제가 된다(rule of zero).
//    옛 InitSerialData / DeleteSerialData 는 필요 없다.
//  - 행마다 값 개수(RowSize)를 따로 둔다. nCol 은 행 길이의 상한이다.
//    삼각형 · 사각형이 섞인 요소(nonmanifold)도 nCol = 4 로 잡고 행마다 3 또는 4 를 둔다.
//  - 행 번호가 범위를 벗어나면 std::out_of_range 를 던진다(옛 vData.at(i) 와 같은 동작).
//  - 옛 operator double*() 암시 변환은 없앴다(delete c; 같은 실수가 컴파일된다). Data() 를 쓴다.
//---------------------------------------------------------------------------
#ifndef SXNamedContainerH
#define SXNamedContainerH

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "SXRow.h"

template <class T>
class SXNamedContainer {
public:
	typedef T value_type;

	std::wstring Name;

	SXNamedContainer() : nRow_(0), nCol_(0) {}
	explicit SXNamedContainer(const std::wstring& name) : Name(name), nRow_(0), nCol_(0) {}

	// nRow × nCol 을 0 으로 채워 잡는다. 행마다 값 개수 = nCol.
	// 이 변수의 행 보기 · Data() 포인터는 무효가 된다.
	void Allocate(int nRow, int nCol)
	{
		CheckShape(nRow, nCol);
		std::vector<T>   data(static_cast<std::size_t>(nRow) * static_cast<std::size_t>(nCol), T());
		std::vector<int> count(static_cast<std::size_t>(nRow), nCol);
		data_.swap(data);
		count_.swap(count);
		nRow_ = nRow;
		nCol_ = nCol;
	}

	void Clear()
	{
		std::vector<T>().swap(data_);
		std::vector<int>().swap(count_);
		nRow_ = nCol_ = 0;
	}

	int  size()  const { return nRow_; }
	int  nCol()  const { return nCol_; }
	bool empty() const { return nRow_ == 0; }

	SXRow<T> operator[](int i)
	{
		CheckRow(i);
		return SXRow<T>(data_.data() + static_cast<std::size_t>(i) * nCol_, count_[i]);
	}
	SXRow<const T> operator[](int i) const
	{
		CheckRow(i);
		return SXRow<const T>(data_.data() + static_cast<std::size_t>(i) * nCol_, count_[i]);
	}
	SXRow<T>       at(int i)       { return (*this)[i]; }
	SXRow<const T> at(int i) const { return (*this)[i]; }

	int RowSize(int i) const { CheckRow(i); return count_[i]; }

	// 행 i 의 값 개수를 n 으로 바꾼다(0 <= n <= nCol). 값은 그대로 둔다.
	void SetRowSize(int i, int n)
	{
		CheckRow(i);
		if (n < 0 || n > nCol_) throw std::out_of_range("SXNamedContainer::SetRowSize: n must be in [0, nCol]");
		count_[i] = n;
	}

	// 행 i 에 값 n 개를 복사하고 값 개수를 n 으로 둔다(옛 Add 와 같은 뜻, n <= nCol 확인 추가).
	void Add(int i, const T* p, int n)
	{
		CheckRow(i);
		if (n < 0 || n > nCol_) throw std::out_of_range("SXNamedContainer::Add: n must be in [0, nCol]");
		if (n > 0 && p == 0) throw std::invalid_argument("SXNamedContainer::Add: null data");
		T* dst = data_.data() + static_cast<std::size_t>(i) * nCol_;
		for (int j = 0; j < n; ++j) dst[j] = p[j];
		count_[i] = n;
	}

	// 연속 블록(nRow × nCol, 행 우선). 비어 있으면 NULL 일 수 있다.
	T*       Data()       { return data_.data(); }
	const T* Data() const { return data_.data(); }

	// 모든 행의 값 개수가 nCol 인가
	bool IsUniform() const
	{
		for (std::size_t i = 0; i < count_.size(); ++i)
			if (count_[i] != nCol_) return false;
		return true;
	}

	// 가장 긴 행의 값 개수
	int MaxRowSize() const
	{
		int m = 0;
		for (std::size_t i = 0; i < count_.size(); ++i)
			if (count_[i] > m) m = count_[i];
		return m;
	}

	static void CheckShape(int nRow, int nCol)
	{
		if (nRow < 0 || nCol < 0) throw std::invalid_argument("SXNamedContainer: nRow and nCol must be >= 0");
		if (nCol > 0 && static_cast<std::size_t>(nRow) > (std::numeric_limits<std::size_t>::max)() / sizeof(T) / static_cast<std::size_t>(nCol))
			throw std::length_error("SXNamedContainer: nRow * nCol too large");
	}

private:
	void CheckRow(int i) const
	{
		if (i < 0 || i >= nRow_) throw std::out_of_range("SXNamedContainer: row index out of range");
	}

	std::vector<T>   data_;
	std::vector<int> count_;
	int nRow_;
	int nCol_;
};

typedef SXNamedContainer<double> SXNodeNamedContainer;     // 노드 쪽 변수(좌표, 법선, UV ...)
typedef SXNamedContainer<int>    SXElementNamedContainer;  // 요소 쪽 변수(꼭짓점 번호, 재질 번호 ...)

#endif
