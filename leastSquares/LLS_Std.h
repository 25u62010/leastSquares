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
class LLS_Std
{
public:
	LLS_Std() {
	}
	Matrix<_Scalar, _identifiedNum, 1> optimize();

	void addInputData(Matrix<_Scalar, 1, _identifiedNum> datas);
	void addInputData(std::vector<_Scalar> datas);
	void addInputData(Matrix<_Scalar, Dynamic, _identifiedNum> datas);
	void addInputData(std::vector<std::vector<_Scalar>> datas);

	void addOutputData(Matrix<_Scalar, Dynamic, 1> datas);
	void addOutputData(std::vector<_Scalar> datas);

	void readTheta(Matrix< _Scalar, _identifiedNum, 1>& result);
	void readTheta(vector< _Scalar>& result);

	void printPHI();
	void printTheta();
	friend ostream & operator<<(ostream &out, LLS_Std<_Scalar, _identifiedNum> &obj) {
		out << "PHI="<<endl<<obj.PHI << endl<<"Y="<<endl << obj.Y << endl<<"theta="<<endl<<obj.theta.transpose()<<endl;
		return out;
	}
protected:
	Index N=0;
	Matrix< _Scalar, Dynamic, _identifiedNum > PHI;
	Matrix< _Scalar, Dynamic, 1> Y;
	Matrix< _Scalar, _identifiedNum, 1> theta;
};

template<typename _Scalar, int _identifiedNum>
Matrix<_Scalar, _identifiedNum, 1> LLS_Std<_Scalar, _identifiedNum>::optimize() {
	Index rows = PHI.rows();
	Matrix<_Scalar, _identifiedNum, Dynamic> phiTranspose;
	phiTranspose.resize(_identifiedNum,rows);
	phiTranspose = PHI.transpose();
	Matrix<_Scalar, _identifiedNum, _identifiedNum> phiTphi = phiTranspose *PHI;

	Matrix<_Scalar, _identifiedNum, _identifiedNum> phiTphiInverse= phiTphi.inverse();
	if (Y.rows() != PHI.rows()) {
		std::cout << "It has a problem of the rows of the Y and the PHI" <<std::endl;
		return theta;
	}
	theta = phiTphiInverse * phiTranspose*Y;
	return theta;
}

template<typename _Scalar, int _identifiedNum>
void LLS_Std<_Scalar, _identifiedNum>::addInputData(Matrix<_Scalar, 1, _identifiedNum> datas) {
	N++;
	PHI.conservativeResize(N, _identifiedNum);
	
	for (int i = 0; i < _identifiedNum;i++) {
		PHI(N - 1, i) = datas(0,i);
	}
}

template<typename _Scalar, int _identifiedNum>
void LLS_Std<_Scalar, _identifiedNum>::addInputData(std::vector<_Scalar> datas) {
	if (datas.size() != _identifiedNum) {
		std::cout << "It has a problem of the size of insered data" << std::endl;
		return;
	}
	N++;
	PHI.conservativeResize(N, _identifiedNum);

	for (Index i = 0; i < _identifiedNum; i++) {
		PHI(N - 1, i) = datas[i];
	}
}
template<typename _Scalar, int _identifiedNum>
void LLS_Std<_Scalar, _identifiedNum>::addInputData(Matrix<_Scalar, Dynamic, _identifiedNum> datas) {
	Index rows = datas.rows();
	N += rows;
	PHI.conservativeResize(N, _identifiedNum);
	for (int row = 0; row <rows; row++) {
		for (int i = 0; i < _identifiedNum; i++) {
			PHI(N - rows+row, i) = datas(row,i);
		}
	}
}
template<typename _Scalar, int _identifiedNum>
void LLS_Std<_Scalar, _identifiedNum>::addInputData(std::vector<std::vector<_Scalar>> datas) {
	Index rows = datas.size();
	for (Index row = 0; row < rows; row++) {
		std::vector<_Scalar> data = datas[row];
		if (data.size()!=_identifiedNum) {
			std::cout << "It has a problem of the size of insered data" << std::endl;
			return;
		}
		addInputData(data);
	}
}

template<typename _Scalar, int _identifiedNum>
void LLS_Std<_Scalar, _identifiedNum>::addOutputData(Matrix<_Scalar, Dynamic, 1> datas) {
	Y.resize(datas.rows(),1);
	Y << datas;
}
template<typename _Scalar, int _identifiedNum>
void LLS_Std<_Scalar, _identifiedNum>::addOutputData(std::vector<_Scalar> datas) {
	Y.resize(datas.size(), 1);
	for (unsigned int i = 0; i < datas.size(); i++) {
		Y(i,0) = datas[i];
	}
}

template<typename _Scalar, int _identifiedNum>
void LLS_Std<_Scalar, _identifiedNum>::readTheta(Matrix< _Scalar, _identifiedNum, 1>& result) {
	result<<theta;
}
template<typename _Scalar, int _identifiedNum>
void LLS_Std<_Scalar, _identifiedNum>::readTheta(vector< _Scalar>& result) {
	for (int i = 0; i < _identifiedNum; i++) {
		result.push_back(theta(i, 0));
	}
}

template<typename _Scalar, int _identifiedNum>
void LLS_Std<_Scalar, _identifiedNum>::printPHI() {
	std::cout << "PHI=" << std::endl << PHI << std::endl;
}

template<typename _Scalar, int _identifiedNum>
void LLS_Std<_Scalar, _identifiedNum>::printTheta() {
	std::cout << "theta=" << std::endl << theta << std::endl;
}
}