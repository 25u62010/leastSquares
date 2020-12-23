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
class linearLeastSquaresStd
{
public:
	linearLeastSquaresStd() {
	}
	Matrix<_Scalar, _identifiedNum, 1> optimize();
	void addInputData(Matrix<_Scalar, 1, _identifiedNum> datas);
	void addInputData(std::vector<_Scalar> datas);
	void addInputData(Matrix<_Scalar, Dynamic, _identifiedNum> datas);
	void addInputData(std::vector<std::vector<_Scalar>> datas);

	void addOutputData(Matrix<_Scalar, Dynamic, 1> datas);
	void addOutputData(std::vector<_Scalar> datas);

	static vector<vector<_Scalar>>  readFromFile(const string &filename,int startCol,int endCol,int startRow=1,int endRow=-1);

	void readTheta(Matrix< _Scalar, _identifiedNum, 1>& result);
	void readTheta(vector< _Scalar>& result);

	void printPHI();
	void printTheta();
	friend ostream & operator<<(ostream &out, linearLeastSquaresStd<_Scalar, _identifiedNum> &obj) {
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
Matrix<_Scalar, _identifiedNum, 1> linearLeastSquaresStd<_Scalar, _identifiedNum>::optimize() {
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
void linearLeastSquaresStd<_Scalar, _identifiedNum>::addInputData(Matrix<_Scalar, 1, _identifiedNum> datas) {
	N++;
	PHI.conservativeResize(N, _identifiedNum);
	
	for (int i = 0; i < _identifiedNum;i++) {
		PHI(N - 1, i) = datas(0,i);
	}
}

template<typename _Scalar, int _identifiedNum>
void linearLeastSquaresStd<_Scalar, _identifiedNum>::addInputData(std::vector<_Scalar> datas) {
	if (datas.size() != _identifiedNum) {
		std::cout << "It has a problem of the size of the insered data" << std::endl;
		return;
	}
	N++;
	PHI.conservativeResize(N, _identifiedNum);

	for (Index i = 0; i < _identifiedNum; i++) {
		PHI(N - 1, i) = datas[i];
	}
}
template<typename _Scalar, int _identifiedNum>
void linearLeastSquaresStd<_Scalar, _identifiedNum>::addInputData(Matrix<_Scalar, Dynamic, _identifiedNum> datas) {
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
void linearLeastSquaresStd<_Scalar, _identifiedNum>::addInputData(std::vector<std::vector<_Scalar>> datas) {
	Index rows = datas.size();
	for (Index row = 0; row < rows; row++) {
		std::vector<_Scalar> data = datas[row];
		if (data.size()!=_identifiedNum) {
			std::cout << "It has a problem of the size of the insered data" << std::endl;
			return;
		}
		addInputData(data);
	}
}

template<typename _Scalar, int _identifiedNum>
void linearLeastSquaresStd<_Scalar, _identifiedNum>::addOutputData(Matrix<_Scalar, Dynamic, 1> datas) {
	Y.resize(datas.rows(),1);
	Y << datas;
}
template<typename _Scalar, int _identifiedNum>
void linearLeastSquaresStd<_Scalar, _identifiedNum>::addOutputData(std::vector<_Scalar> datas) {
	Y.resize(datas.size(), 1);
	for (unsigned int i = 0; i < datas.size(); i++) {
		Y(i,0) = datas[i];
	}
}

template<typename _Scalar, int _identifiedNum>
vector<vector<_Scalar>> linearLeastSquaresStd<_Scalar, _identifiedNum>::readFromFile(const string &filename,int startCol,int endCol,int startRow,int endRow) {
	vector<vector<_Scalar>> result;
	//输入文件流
	ifstream fAssociation;
	//打开关联文件
	fAssociation.open(filename.c_str());
	int row = 1;
	int size = 0;
	string s;
	if (!fAssociation.is_open()) {
		cout << "cann't open file"<<filename << endl;
		return result;
	}
	//一直读取,知道文件结束
	while (!fAssociation.eof()){
		//读取一行的内容到字符串s中
		getline(fAssociation, s);
		//如果不是空行就可以分析数据了
		if (!s.empty() && row >= startRow&&(row<=endRow|| endRow ==-1)) {
			//字符串流
			stringstream ss;
			ss << s;
			_Scalar val;
			int col = 1;
			result.push_back(vector<_Scalar>());

			while (ss >> val) {
				if (col >= startCol&&col <= endCol) {
					result[size].push_back(val);
				}
				col++;
			}
			size++;
		}	
		row++;
	}
	fAssociation.close();
	return result;
}
template<typename _Scalar, int _identifiedNum>
void linearLeastSquaresStd<_Scalar, _identifiedNum>::readTheta(Matrix< _Scalar, _identifiedNum, 1>& result) {
	result<<theta;
}
template<typename _Scalar, int _identifiedNum>
void linearLeastSquaresStd<_Scalar, _identifiedNum>::readTheta(vector< _Scalar>& result) {
	for (int i = 0; i < _identifiedNum; i++) {
		result.push_back(theta(i, 0));
	}
}

template<typename _Scalar, int _identifiedNum>
void linearLeastSquaresStd<_Scalar, _identifiedNum>::printPHI() {
	std::cout << "PHI=" << std::endl << PHI << std::endl;
}

template<typename _Scalar, int _identifiedNum>
void linearLeastSquaresStd<_Scalar, _identifiedNum>::printTheta() {
	std::cout << "theta=" << std::endl << theta << std::endl;
}


}