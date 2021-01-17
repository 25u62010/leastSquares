#pragma once
#include <Eigen/core>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

using namespace Eigen;
using namespace std;

namespace zlzLS {

template<typename _Scalar, int _identifiedNum>
class RLLS_Std
{
public:
	RLLS_Std(Matrix<_Scalar, _identifiedNum, 1> theta0, Matrix<_Scalar, _identifiedNum, _identifiedNum> P0) :theta(theta0), P(P0) {

	}
	RLLS_Std(_Scalar theta0, _Scalar P0) {
		for (int i = 0; i < _identifiedNum;i++) {
			theta(i,0) = theta0;
			
			for (int j = 0; j < _identifiedNum; j++) {
				P(i,j) = 0;
			}
			P(i, i) = P0;
		}
	}
	RLLS_Std(_Scalar theta0[_identifiedNum], _Scalar P0=INT_MAX) {
		for (int i = 0; i < _identifiedNum; i++) {
			theta(i, 0) = theta0[i];

			for (int j = 0; j < _identifiedNum; j++) {
				P(i, j) = 0;
			}
			P(i, i) = P0;
		}
	}
	RLLS_Std() {
		for (int i = 0; i < _identifiedNum; i++) {
			theta(i, 0) = 0;

			for (int j = 0; j < _identifiedNum; j++) {
				P(i, j) = 0;
			}
			P(i, i) = 1000;
		}
	}
	virtual void step(Matrix<_Scalar, _identifiedNum, 1> phi, _Scalar Y);//读取数据并迭代
	virtual void step(vector<_Scalar> phi, _Scalar Y);
	virtual _Scalar weightFcn();
	void setRecordTheta(bool flag) {//设置是否记录数据
		recordTheta = flag;
	};
	//读取数据
	void readTheta(Matrix<_Scalar, _identifiedNum, 1>& _theta) {
		_theta = theta;
	}
	void readTheta(vector<_Scalar>& _theta) {
		_theta = vector<_Scalar>(_identifiedNum);
		for (int i = 0; i < _identifiedNum; i++) {
			_theta[i] = theta(i,0);
		}
	}
	void writeThetaBuffer(ostream &out);
protected:
	Matrix<_Scalar, _identifiedNum, 1> theta;
	//数据缓存区
	vector<Matrix<_Scalar, _identifiedNum, 1>> thetaBuffer;
	Matrix<_Scalar, _identifiedNum, _identifiedNum> P;
	Matrix<_Scalar, _identifiedNum, 1> K;
	//记录数据标志
	bool recordTheta;

};

template<typename _Scalar, int _identifiedNum>
void RLLS_Std<_Scalar, _identifiedNum>::step(Matrix<_Scalar, _identifiedNum,1> phi, _Scalar Y) {
	_Scalar lamda = weightFcn();
	Matrix<_Scalar,1, _identifiedNum> phiTranspose=phi.transpose();
	Matrix<_Scalar, 1, 1> Phi_P_PhiTemp = phiTranspose *P*phi;
	_Scalar Phi_P_Phi = Phi_P_PhiTemp(0, 0);
	K = P / (1/ lamda + Phi_P_Phi) * phi;
	Matrix<_Scalar, 1, 1> phiT_theta = phiTranspose * theta;
	theta += K*(Y- phiT_theta(0,0));
	P += -K* K.transpose()*(1+ Phi_P_Phi);
	if (recordTheta) {
		thetaBuffer.push_back(theta);
	}
}
template<typename _Scalar, int _identifiedNum>
void RLLS_Std<_Scalar, _identifiedNum>::step(vector<_Scalar> phi, _Scalar Y) {
	if (phi.size() != _identifiedNum) {
		cout << "wrong with the size of phi" << endl;
		return;
	}
	Matrix<_Scalar, _identifiedNum, 1> phiM;
	for (int i = 0; i < _identifiedNum;i++) {
		phiM(i, 1) = phi[i];
	}
	step(phiM, Y);
}
template<typename _Scalar, int _identifiedNum>
_Scalar RLLS_Std<_Scalar, _identifiedNum>::weightFcn() {
	return _Scalar(1);
}
template<typename _Scalar, int _identifiedNum>
void RLLS_Std<_Scalar, _identifiedNum>::writeThetaBuffer(ostream &out) {
	if (thetaBuffer.size() == 0) {
		return;
	}
	for (Matrix<_Scalar, _identifiedNum, 1> thetaI : thetaBuffer) {
		out << thetaI.transpose()<<endl;
	}
	out << endl;
}
}
