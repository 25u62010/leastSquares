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
void LLS_StdTest(double a, double b);//��ͨ��С����
void discreteLLSTest(string filePath);//���ڿ���ϵͳ����С����
void RLLS_StdTest(double a, double b);//������С����
void discreteRLLSTest(string filePath);//���ڿ���ϵͳ�ĵ�����С����
void forgettingFactorRLLSTest(double a, double b);//����������С����
void discreteForgettingFactorRLLSTest(string filePath);//���ڿ���ϵͳ�Ľ����������
int main(int argc,char* argv[]) {
	
	LLS_StdTest(2,3);
	discreteLLSTest("..\\dataGenerator.txt");
	RLLS_StdTest(2,3);
	discreteRLLSTest("..\\dataGenerator.txt");
	forgettingFactorRLLSTest(2, 3);
	discreteForgettingFactorRLLSTest("..\\dataGenerator.txt");
}
void LLS_StdTest(double a, double b) {//��ͨ��С����
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
void discreteLLSTest(string filePath) {//���ڿ���ϵͳ����С���˽�
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
void RLLS_StdTest(double a,double b) {//������С����
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
void discreteRLLSTest(string filePath) {//���ڿ���ϵͳ�ĵ�����С����
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
void forgettingFactorRLLSTest(double a, double b) {//����������С����
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
void discreteForgettingFactorRLLSTest(string filePath) {//���ڿ���ϵͳ�ļ�С������С����
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