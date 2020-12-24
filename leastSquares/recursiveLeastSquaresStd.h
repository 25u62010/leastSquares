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
class recursiveLeastSquaresStd
{
public:
	recursiveLeastSquaresStd(Matrix<_Scalar, _identifiedNum, 1> theta0, Matrix<_Scalar, _identifiedNum, _identifiedNum> P0) :theta(theta0), P(P0) {

	}
	recursiveLeastSquaresStd(_Scalar theta0, _Scalar P0) {
		for (int i = 0; i < _identifiedNum;i++) {
			theta(i,0) = theta0;
			
			for (int j = 0; j < _identifiedNum; j++) {
				P(i,j) = 0;
			}
			P(i, i) = P0;
		}
	}
	recursiveLeastSquaresStd() {
		for (int i = 0; i < _identifiedNum; i++) {
			theta(i, 0) = 0;

			for (int j = 0; j < _identifiedNum; j++) {
				P(i, j) = 0;
			}
			P(i, i) = 1000;
		}
	}
	void step(Matrix<_Scalar, _identifiedNum, 1> phi, _Scalar Y);
	void step(vector<_Scalar> phi, _Scalar Y);
	void setRecordTheta(bool flag) {
		recordTheta = flag;
	};
	void readTheta(Matrix<_Scalar, _identifiedNum, 1> _theta) {
		_theta = theta;
	}
	void readTheta(vector<_Scalar> _theta) {
		_theta = vector<_Scalar>(_identifiedNum);
		for (int i = 0; i < _identifiedNum; i++) {
			_theta[i] = theta(i,0);
		}
	}
	void writeThetaBuffer(ostream &out);
private:
	Matrix<_Scalar, _identifiedNum, 1> theta;
	vector<Matrix<_Scalar, _identifiedNum, 1>> thetaBuffer;
	Matrix<_Scalar, _identifiedNum, _identifiedNum> P;
	Matrix<_Scalar, _identifiedNum, 1> K;

	bool recordTheta;

};

template<typename _Scalar, int _identifiedNum>
void recursiveLeastSquaresStd<_Scalar, _identifiedNum>::step(Matrix<_Scalar, _identifiedNum,1> phi, _Scalar Y) {
	Matrix<_Scalar,1, _identifiedNum> phiTranspose=phi.transpose();
	Matrix<_Scalar, 1, 1> Phi_P_PhiTemp = phiTranspose *P*phi;
	_Scalar Phi_P_Phi = Phi_P_PhiTemp(0, 0);
	K = P / (1 + Phi_P_Phi) * phi;
	Matrix<_Scalar, 1, 1> phiT_theta = phiTranspose * theta;
	theta += K*(Y- phiT_theta(0,0));
	P += -K* phiTranspose*P;
	if (recordTheta) {
		thetaBuffer.push_back(theta);
	}
}
template<typename _Scalar, int _identifiedNum>
void recursiveLeastSquaresStd<_Scalar, _identifiedNum>::step(vector<_Scalar> phi, _Scalar Y) {
	if (phi.size() != _identifiedNum || Y.size() != _identifiedNum) {
		cout << "wrong with the size of phi or Y" << endl;
		return;
	}
	Matrix<_Scalar, _identifiedNum, 1> phiM;
	for (int i = 0; i < _identifiedNum;i++) {
		phiM(i, 1) = phi[i];
	}
	step(phiM, Y);
}
template<typename _Scalar, int _identifiedNum>
void recursiveLeastSquaresStd<_Scalar, _identifiedNum>::writeThetaBuffer(ostream &out) {
	if (thetaBuffer.size() == 0) {
		return;
	}
	for (Matrix<_Scalar, _identifiedNum, 1> thetaI : thetaBuffer) {
		out << thetaI.transpose()<<endl;
	}
	out << endl;
}
}