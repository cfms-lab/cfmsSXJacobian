//---------------------------------------------------------------------------
// sx_stencil.h — 명세가 선언한 변수를 SXMesh 에 이름으로 만들고, 스텐실마다 값을 모은다.
//
// SXMesh 의 두 칸은 값의 종류로 나뉜다: NodeV = double, ElementV = int.
// 그래서 스텐실(힌지 · 삼각형)마다 붙는 double 값(쉬는 각 · 면적 ...)도 NodeV 에 둔다(행 수 = 스텐실 수).
// 각 변수는 자기 배열 하나를 가진다(변수 단위 SoA). 한 변수 안에서는 행 우선(x y z x y z ...)이다.
//---------------------------------------------------------------------------
#pragma once

#include <stdexcept>
#include <string>

#include "../third_party/sxmesh/SXMesh.h"

namespace sx {

// 명세가 선언한 변수를 모두 만든다(0 으로 채움). 같은 이름이 있으면 다시 잡는다.
template <class Spec>
void CreateVariables(SXMesh& m, int nNode, int nStencil)
{
	for (int i = 0; i < Spec::kNodeVars; ++i) m.NewNodeV(Spec::NodeVar(i), nNode, Spec::NodeVarCols(i));
	m.NewElementV(Spec::StencilVar(), nStencil, Spec::kNodes);
	for (int i = 0; i < Spec::kElemVars; ++i) m.NewNodeV(Spec::ElemVar(i), nStencil, Spec::ElemVarCols(i));
}

// 스텐실 s 의 좌표 · 속도 · 스텐실 값을 모으는 읽기 전용 보기. 변수는 만들 때 한 번만 찾는다.
template <class Spec>
class StencilView {
public:
	explicit StencilView(const SXMesh& m)
		: x_(Need(m.FindNodeV(Spec::PositionVar()), 3)), v_(m.FindNodeV(Spec::VelocityVar())),
		  st_(m.FindElementV(Spec::StencilVar())), nNode_(x_->size())
	{
		if (st_ == 0 || st_->nCol() < Spec::kNodes) throw std::invalid_argument("StencilView: stencil variable missing or too narrow");
		if (v_ != 0 && (v_->nCol() < 3 || v_->size() != nNode_)) throw std::invalid_argument("StencilView: velocity shape");
		for (int i = 0; i < Spec::kElemVars; ++i) {
			ev_[i] = m.FindNodeV(Spec::ElemVar(i));
			if (ev_[i] == 0 || ev_[i]->size() != st_->size() || ev_[i]->nCol() < Spec::ElemVarCols(i))
				throw std::invalid_argument("StencilView: stencil value variable missing or wrong shape");
		}
	}

	int size() const { return st_->size(); }
	bool HasVelocity() const { return v_ != 0; }
	const int* Nodes(int s) const { return st_->Data() + static_cast<std::size_t>(s) * st_->nCol(); }

	// 번호가 범위 밖이면 false
	bool Gather(int s, double* x, double* v, double* ev) const
	{
		const int* n = Nodes(s);
		for (int a = 0; a < Spec::kNodes; ++a) {
			if (n[a] < 0 || n[a] >= nNode_) return false;
			const double* p = x_->Data() + static_cast<std::size_t>(n[a]) * x_->nCol();
			x[3 * a] = p[0]; x[3 * a + 1] = p[1]; x[3 * a + 2] = p[2];
			if (v != 0) {
				if (v_ != 0) {
					const double* q = v_->Data() + static_cast<std::size_t>(n[a]) * v_->nCol();
					v[3 * a] = q[0]; v[3 * a + 1] = q[1]; v[3 * a + 2] = q[2];
				}
				else v[3 * a] = v[3 * a + 1] = v[3 * a + 2] = 0.;
			}
		}
		int k = 0;
		for (int i = 0; i < Spec::kElemVars; ++i) {
			const double* p = ev_[i]->Data() + static_cast<std::size_t>(s) * ev_[i]->nCol();
			for (int c = 0; c < Spec::ElemVarCols(i); ++c) ev[k++] = p[c];
		}
		return true;
	}

private:
	static const SXNodeNamedContainer* Need(const SXNodeNamedContainer* c, int nCol)
	{
		if (c == 0 || c->nCol() < nCol) throw std::invalid_argument("StencilView: position variable missing or narrow");
		return c;
	}

	const SXNodeNamedContainer*    x_;
	const SXNodeNamedContainer*    v_;
	const SXElementNamedContainer* st_;
	const SXNodeNamedContainer*    ev_[Spec::kElemVars > 0 ? Spec::kElemVars : 1];
	int nNode_;
};

} // namespace sx
