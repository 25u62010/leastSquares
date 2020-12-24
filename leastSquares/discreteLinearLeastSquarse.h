#pragma once
#include "iostream"
#include "linearLeastSquaresStd.h"
using namespace std;
namespace zlzLS {
template<typename _Scalar, int _identifiedNum,int _rank>
class discreteLinearSquarse{
	public:
		discreteLinearSquarse() {

		}
		~discreteLinearSquarse() {
			
		}
		void addInputData(vector<_Scalar> datas);
		void addOutputData(vector<_Scalar> datas);
		vector<_Scalar> optimize();
		static vector<_Scalar>  readFromFile(const string &filename, int col, int startRow = 1, int endRow = -1);
		void readTheta(Matrix< _Scalar, _identifiedNum, 1>& result) {
			core.readTheta(result);
		}
		void readTheta(vector< _Scalar>& result) {
			core.readTheta(result);
		}
		friend ostream & operator<<(ostream &out, discreteLinearSquarse<_Scalar, _identifiedNum,_rank> &obj) {
			Matrix< _Scalar, _identifiedNum, 1> theta;
			obj.core.readTheta(theta);
			out<<"theta=" << endl << theta.transpose() << endl;
			return out;
		}
	private:
		linearLeastSquaresStd<_Scalar, _identifiedNum> core;
		vector<_Scalar> inputBuffer;
		vector<_Scalar> outputBuffer;

};

template<typename _Scalar, int _identifiedNum, int _rank>
void discreteLinearSquarse<_Scalar, _identifiedNum, _rank>::addInputData(vector<_Scalar> datas) {
	for (_Scalar data : datas) {
		inputBuffer.push_back(data);
	}
}
template<typename _Scalar, int _identifiedNum, int _rank>
void discreteLinearSquarse<_Scalar, _identifiedNum, _rank>::addOutputData(vector<_Scalar> datas) {
	for (_Scalar data : datas) {
		outputBuffer.push_back(data);
	}
}
template<typename _Scalar, int _identifiedNum, int _rank>
vector<_Scalar> discreteLinearSquarse<_Scalar, _identifiedNum, _rank>::optimize() {
	vector<_Scalar> result;
	if (inputBuffer.size() != outputBuffer.size()) {
		cout << "The number of inputs and outputs varies" << endl;
		return result;
	}
	if (inputBuffer.size()<_identifiedNum) {
		cout << "Data volume is too small" << endl;
		return result;
	}
	int identifieNum = _identifiedNum;
	int rank = _rank;
	int outputNum = identifieNum - rank;	
	size_t bufferSize = outputBuffer.size();
	vector<_Scalar> Y(outputBuffer.begin()+ identifieNum,outputBuffer.end());
	core.addOutputData(Y);
	size_t N = bufferSize - identifieNum;
	Matrix<_Scalar,Dynamic,_identifiedNum> PHI;
	PHI.resize(N,_identifiedNum);
	for (int i = identifieNum; i < bufferSize; i++) {
		for (int j = 0; j < outputNum; j++) {
			PHI(i - identifieNum,j) = -outputBuffer[i-j-1];
		}
		for (int j = 0; j < rank; j++) {
			PHI(i - identifieNum,j+ outputNum) = inputBuffer[i - j - 1];
		}
	}
	core.addInputData(PHI);
	core.optimize();
	
	core.readTheta(result);
	return result;
}
template<typename _Scalar, int _identifiedNum, int _rank>
vector<_Scalar> discreteLinearSquarse<_Scalar, _identifiedNum, _rank>::readFromFile(const string &filename, int col, int startRow, int endRow) {
	vector<_Scalar> result;
	ifstream fAssociation;
	fAssociation.open(filename.c_str());
	int row = 1;
	string s;
	if (!fAssociation.is_open()) {
		cout << "cann't open file" << filename << endl;
		return result;
	}
	while (!fAssociation.eof()) {
		getline(fAssociation, s);
		if (!s.empty() && row >= startRow && (row <= endRow || endRow == -1)) {
			stringstream ss;
			ss << s;
			_Scalar val;
			int nowCol = 1;

			while (ss >> val) {
				if (col == nowCol) {
					result.push_back(val);
				}
				nowCol++;
			}
		}
		row++;
	}
	fAssociation.close();
	return result;
}
}

