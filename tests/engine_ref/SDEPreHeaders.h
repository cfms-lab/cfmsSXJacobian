//---------------------------------------------------------------------------
// SDEPreHeaders.h — 엔진 cfmsDrape_dev 의 같은 이름 헤더 대역(엔진 손 코드 비교 게이트, docs/GATE_ENGINE.md).
// generated/engine_ref/SDEmysse.h 가 이 이름으로 #include 한다. 뽑아 온 엔진 코드가 쓰는 것만, 엔진 값 그대로 둔다
// (378d20c drape/src/sde_core/SDEPreHeaders.h 의 typedef · enum · 매크로).
//---------------------------------------------------------------------------
#pragma once

#include <math.h>

typedef double dFloat;
typedef double *dFloatVec;

enum SDEImplicitBendingForceType {
	SDEibfNone = 0,
	SDEibfImplicitBaraff = 1,
	SDEibfImplicitBridson,
	SDEibfImplicitSyflex };

#define _DO_NOT_USE_SSE
#define _GLOBAL_IMPLICIT_BENDING_2ND_DERIVATIVE
// OMP 는 정의하지 않는다 — 단일 스레드로 견준다(엔진 빌드는 CMake 가 OMP 를 정의한다)
