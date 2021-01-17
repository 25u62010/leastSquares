#pragma once
#include "RLLS_Std.h"
namespace zlzLS {
template<typename _Scalar, int _identifiedNum, int _rank>
class discreteRLLS:public RLLS_Std<_Scalar, _identifiedNum>
{
public:
	discreteRLLS();
	discreteRLLS(_Scalar theta0, _Scalar P0);
	discreteRLLS(_Scalar theta0[_identifiedNum], _Scalar P0 = 10000);
	~discreteRLLS();
	virtual void step(_Scalar Y,_Scalar u);
private:
	Matrix<_Scalar, _identifiedNum, 1> _phi;
	Index _step;
};

template<typename _Scalar, int _identifiedNum, int _rank>
discreteRLLS<_Scalar, _identifiedNum, _rank>::discreteRLLS()
	:RLLS_Std<_Scalar, _identifiedNum>(){
}

template<typename _Scalar, int _identifiedNum, int _rank>
discreteRLLS<_Scalar, _identifiedNum, _rank>::discreteRLLS(_Scalar theta0, _Scalar P0)
	: RLLS_Std<_Scalar, _identifiedNum>(theta0,P0) {
}

template<typename _Scalar, int _identifiedNum, int _rank>
discreteRLLS<_Scalar, _identifiedNum, _rank>::discreteRLLS(_Scalar theta0[_identifiedNum], _Scalar P0)
	: RLLS_Std<_Scalar, _identifiedNum>(theta0[_identifiedNum], P0) {
}

template<typename _Scalar, int _identifiedNum, int _rank>
void discreteRLLS<_Scalar, _identifiedNum, _rank>::step(_Scalar Y, _Scalar u){
	if (_step > _identifiedNum) {
		RLLS_Std<_Scalar, _identifiedNum>::step(_phi,Y);
	}
	_step++;
	for (Index i = _identifiedNum-1; i > _rank; --i) {
		_phi(i, 0) = _phi(i - 1, 0);
	}
	_phi(_rank, 0) = u;
	for (Index i = _rank - 1; i > 0; --i) {
		_phi(i, 0) = _phi(i - 1, 0);
	}
	_phi(0, 0) = -Y;
}

template<typename _Scalar, int _identifiedNum, int _rank>
discreteRLLS<_Scalar, _identifiedNum, _rank>::~discreteRLLS(){

}
}

