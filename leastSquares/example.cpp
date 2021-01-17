#include "iostream"
#include <Eigen/core>
#include "LLS_Std.h"
#include "discreteLLS.h"
#include "RLLS_Std.h"
#include "discreteRLLS.h"
#include "forgettingFactorRLLS.h"
#include "discreteForgettingFactorRLLS.h"
#include "fileOperations.h"
#include <fstream>
#include <string>
#include <sstream>
using namespace std;
using namespace Eigen;
using namespace zlzLS;
void LLS_StdTest(double a, double b);//普通最小二乘
void discreteLLSTest(string filePath);//用于控制系统的最小二乘
void RLLS_StdTest(double a, double b);//递推最小二乘
void discreteRLLSTest(string filePath);//用于控制系统的递推最小二乘
void forgettingFactorRLLSTest(double a, double b);//渐消记忆最小二乘
void discreteForgettingFactorRLLSTest(string filePath);//用于控制系统的渐消记忆二乘
int main(int argc,char* argv[]) {
	
	LLS_StdTest(2,3);
	discreteLLSTest("..\\dataGenerator.txt");
	RLLS_StdTest(2,3);
	discreteRLLSTest("..\\dataGenerator.txt");
	forgettingFactorRLLSTest(2, 3);
	discreteForgettingFactorRLLSTest("..\\dataGenerator.txt");
}
void LLS_StdTest(double a, double b) {//普通最小二乘
	LLS_Std<double, 2> test;
	vector<vector<double>> PHITest;
	zlzLS::fileOperations<double>::readFromFile(".\\test.txt", PHITest,1, 2);
	test.addInputData(PHITest);
	std::vector<double> YTest;
	for (std::vector<double> phi : PHITest) {
		YTest.push_back(a*phi[0]+b*phi[1]);
	}
	test.addOutputData(YTest);
	test.optimize();
	Matrix<double,2,1> theta;
	test.readTheta(theta);
	cout << "LLS_StdTest:" << theta.transpose() <<endl;
}
void discreteLLSTest(string filePath) {//用于控制系统的最小二乘解
	discreteLLS<double, 4, 2> test;
	vector<double> inputTest;
	zlzLS::fileOperations<double>::readFromFile(filePath, inputTest, 1);
	vector<double> outputTest;
	zlzLS::fileOperations<double>::readFromFile(filePath, outputTest, 2);
	test.addInputData(inputTest);
	test.addOutputData(outputTest);
	test.optimize();
	Matrix<double, 4, 1> theta;
	test.readTheta(theta);
	cout << "discreteLLSTest:" <<theta.transpose()<< endl;
}
void RLLS_StdTest(double a,double b) {//递推最小二乘
	RLLS_Std<double, 2> test;
	test.setRecordTheta(true);
	double YTest;
	Matrix<double, 2, 1> phi(10000, 0);
	for (double x = 15; x < 1000; x += 0.1) {
		YTest = a * x + b;
		phi(0, 0) = x;
		phi(1, 0) = 1;
		test.step(phi, YTest);
	}
	Matrix<double, 2, 1> theta;
	test.readTheta(theta);
	cout << "RLLS_StdTest:" << theta.transpose() << endl;
	ofstream fAssociation;
	fAssociation.open("RLLS_StdTest.txt");
	test.writeThetaBuffer(fAssociation);
	fAssociation.close();
}
void discreteRLLSTest(string filePath) {//用于控制系统的递推最小二乘
	discreteRLLS<double,4, 2> test;

	fileOperations<double> dataFile(filePath);
	if (!dataFile.isOpen()) {
		cout << "cannot open file" << endl;
	}
	test.setRecordTheta(true);
	while (!dataFile.isEnd()) {
		vector<double> data = dataFile.readLine(1, 2);
		if (data.size() == 2) {
			double u = data[0];
			double Y = data[1];
			test.step(Y, u);
		}
	}
	Matrix<double, 4, 1> theta;
	test.readTheta(theta);
	cout << "discreteRLLSTest:" << theta.transpose() << endl;
	ofstream fAssociation;
	fAssociation.open("RLLS_StdTest.txt");
	test.writeThetaBuffer(fAssociation);
	fAssociation.close();
}
void forgettingFactorRLLSTest(double a, double b) {//渐消记忆最小二乘
	forgettingFactorRLLS<double, 2> test(1.12);
	test.setRecordTheta(true);
	double YTest;
	Matrix<double, 2, 1> phi(10000, 0);
	for (double x = 15; x < 1000; x += 0.1) {
		YTest = a * x + b;
		phi(0, 0) = x;
		phi(1, 0) = 1;
		test.step(phi, YTest);
	}
	Matrix<double, 2, 1> theta;
	test.readTheta(theta);
	cout << "forgettingFactorRLLSTest:" << theta.transpose() << endl;
	ofstream fAssociation;
	fAssociation.open("forgettingFactorRLLS.txt");
	test.writeThetaBuffer(fAssociation);
	fAssociation.close();
}
void discreteForgettingFactorRLLSTest(string filePath) {//用于控制系统的减小基于最小二乘
	discreteForgettingFactorRLLS<double, 4, 2> test(1.0001);

	fileOperations<double> dataFile(filePath);
	if (!dataFile.isOpen()) {
		cout << "cannot open file" << endl;
	}
	test.setRecordTheta(true);
	while (!dataFile.isEnd()) {
		vector<double> data = dataFile.readLine(1, 2);
		if (data.size() == 2) {
			double u = data[0];
			double Y = data[1];
			test.step(Y, u);
		}
	}
	Matrix<double, 4, 1> theta;
	test.readTheta(theta);
	cout << "discreteForgettingFactorRLLS:" << theta.transpose() << endl;
	ofstream fAssociation;
	fAssociation.open("RLLS_StdTest.txt");
	test.writeThetaBuffer(fAssociation);
	fAssociation.close();
}