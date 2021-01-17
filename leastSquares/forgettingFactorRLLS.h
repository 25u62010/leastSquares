#pragma once
#include "iostream"
#include "RLLS_Std.h"
using namespace std;
namespace zlzLS {
template<typename _Scalar, int _identifiedNum>
class forgettingFactorRLLS :public RLLS_Std <_Scalar, _identifiedNum> {
public:
	forgettingFactorRLLS();
	forgettingFactorRLLS(_Scalar lamda);
	forgettingFactorRLLS(_Scalar theta0, _Scalar P0, _Scalar lamda=-1);
	forgettingFactorRLLS(_Scalar theta0[_identifiedNum], _Scalar P0 = 10000, _Scalar lamda=-1);
	
	virtual void step(Matrix<_Scalar, _identifiedNum, 1> phi, _Scalar Y);//读取数据并迭代
	virtual void step(vector<_Scalar> phi, _Scalar Y);
protected:
	_Scalar _lamda=1;
};
template<typename _Scalar, int _identifiedNum>
forgettingFactorRLLS<_Scalar, _identifiedNum>::forgettingFactorRLLS()
								:RLLS_Std <_Scalar, _identifiedNum>(){
}
template<typename _Scalar, int _identifiedNum>
forgettingFactorRLLS<_Scalar, _identifiedNum>::forgettingFactorRLLS(_Scalar lamda)
	: RLLS_Std <_Scalar, _identifiedNum>(),_lamda(lamda) {
}
template<typename _Scalar, int _identifiedNum>
forgettingFactorRLLS<_Scalar, _identifiedNum>::forgettingFactorRLLS(_Scalar theta0, _Scalar P0, _Scalar lamda)
	: RLLS_Std <_Scalar, _identifiedNum>(theta0,P0), _lamda(lamda) {
}
template<typename _Scalar, int _identifiedNum>
forgettingFactorRLLS<_Scalar, _identifiedNum>::forgettingFactorRLLS(_Scalar theta0[_identifiedNum], _Scalar P0, _Scalar lamda)
	: RLLS_Std <_Scalar, _identifiedNum>(theta0, P0), _lamda(lamda) {
}

template<typename _Scalar, int _identifiedNum>
void forgettingFactorRLLS<_Scalar, _identifiedNum>::step(Matrix<_Scalar, _identifiedNum, 1> phi, _Scalar Y) {
	Matrix<_Scalar, 1, _identifiedNum> phiTranspose = phi.transpose();
	Matrix<_Scalar, 1, 1> Phi_P_PhiTemp = phiTranspose * (this->P)*phi;
	_Scalar Phi_P_Phi = Phi_P_PhiTemp(0, 0);
	this->K = this->P / (_lamda + Phi_P_Phi) * phi;
	Matrix<_Scalar, 1, 1> phiT_theta = phiTranspose * this->theta;
	this->theta += this->K * (Y - phiT_theta(0, 0));
	this->P += -this->K * this->K.transpose()*(1 + Phi_P_Phi);
	this->P *= 1 / _lamda;
	if (this->recordTheta) {
		(this->thetaBuffer).push_back(this->theta);
	}
}
template<typename _Scalar, int _identifiedNum>
void forgettingFactorRLLS<_Scalar, _identifiedNum>::step(vector<_Scalar> phi, _Scalar Y) {
	if (phi.size() != _identifiedNum) {
		cout << "wrong with the size of phi" << endl;
		return;
	}
	Matrix<_Scalar, _identifiedNum, 1> phiM;
	for (int i = 0; i < _identifiedNum; i++) {
		phiM(i, 1) = phi[i];
	}
	step(phiM, Y);
}
}