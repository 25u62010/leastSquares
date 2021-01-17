#pragma once
#include "iostream"
#include "LLS_Std.h"
using namespace std;
namespace zlzLS {
template<typename _Scalar, int _identifiedNum,int _rank>
class discreteLLS{
public:
	discreteLLS() {

	}
	~discreteLLS() {
			
	}
	void addInputData(vector<_Scalar> datas);
	void addOutputData(vector<_Scalar> datas);
	vector<_Scalar> optimize();
	void readTheta(Matrix< _Scalar, _identifiedNum, 1>& result) {
		core.readTheta(result);
	}
	void readTheta(vector< _Scalar>& result) {
		core.readTheta(result);
	}
	friend ostream & operator<<(ostream &out, discreteLLS<_Scalar, _identifiedNum,_rank> &obj) {
		Matrix< _Scalar, _identifiedNum, 1> theta;
		obj.core.readTheta(theta);
		out<<"theta=" << endl << theta.transpose() << endl;
		return out;
	}
private:
	LLS_Std<_Scalar, _identifiedNum> core;
	vector<_Scalar> inputBuffer;
	vector<_Scalar> outputBuffer;

};

template<typename _Scalar, int _identifiedNum, int _rank>
void discreteLLS<_Scalar, _identifiedNum, _rank>::addInputData(vector<_Scalar> datas) {
	for (_Scalar data : datas) {
		inputBuffer.push_back(data);
	}
}
template<typename _Scalar, int _identifiedNum, int _rank>
void discreteLLS<_Scalar, _identifiedNum, _rank>::addOutputData(vector<_Scalar> datas) {
	for (_Scalar data : datas) {
		outputBuffer.push_back(data);
	}
}
template<typename _Scalar, int _identifiedNum, int _rank>
vector<_Scalar> discreteLLS<_Scalar, _identifiedNum, _rank>::optimize() {
	vector<_Scalar> result;
	if (inputBuffer.size() != outputBuffer.size()) {
		cout << "error: The number of inputs and outputs different" << endl;
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
}

