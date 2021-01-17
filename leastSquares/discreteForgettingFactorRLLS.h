#pragma once
#include "forgettingFactorRLLS.h"
using namespace std;
namespace zlzLS {
template<typename _Scalar, int _identifiedNum, int _rank>
class discreteForgettingFactorRLLS:public forgettingFactorRLLS<_Scalar, _identifiedNum> {
public:
	discreteForgettingFactorRLLS() 
		:forgettingFactorRLLS<_Scalar, _identifiedNum>() {}
	discreteForgettingFactorRLLS(_Scalar lamda) 
		:forgettingFactorRLLS<_Scalar, _identifiedNum>(lamda) {}
	discreteForgettingFactorRLLS(_Scalar theta0, _Scalar P0, _Scalar lamda = -1) 
		:forgettingFactorRLLS<_Scalar, _identifiedNum>(theta0, P0,lamda) {}
	discreteForgettingFactorRLLS(_Scalar theta0[_identifiedNum], _Scalar P0 = 10000, _Scalar lamda = -1)
		:forgettingFactorRLLS<_Scalar, _identifiedNum>(theta0[_identifiedNum],P0,lamda) {}
	virtual void step(_Scalar Y, _Scalar u);
private:
	Matrix<_Scalar, _identifiedNum, 1> _phi;
	Index _step;
};
template<typename _Scalar, int _identifiedNum, int _rank>
void discreteForgettingFactorRLLS<_Scalar, _identifiedNum, _rank>::step(_Scalar Y, _Scalar u) {
	if (_step > _identifiedNum) {
		forgettingFactorRLLS<_Scalar, _identifiedNum>::step(_phi, Y);
	}
	_step++;
	for (Index i = _identifiedNum - 1; i > _rank; --i) {
		_phi(i, 0) = _phi(i - 1, 0);
	}
	_phi(_rank, 0) = u;
	for (Index i = _rank - 1; i > 0; --i) {
		_phi(i, 0) = _phi(i - 1, 0);
	}
	_phi(0, 0) = -Y;
}
}
