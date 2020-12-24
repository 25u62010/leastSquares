#include "iostream"
#include <Eigen/core>
#include "linearLeastSquaresStd.h"
#include "discreteLinearLeastSquarse.h"
#include "recursiveLeastSquaresStd.h"
#include <fstream>
#include <string>
#include <sstream>
using namespace std;
using namespace Eigen;
using namespace zlzLS;
int main(int argc,char* argv[]) {
	//Matrix<double, Dynamic,5> data;
	//Matrix<double, Dynamic, 1> dataY;
	//dataY.resize(5,1);
	//dataY << 1, 2, 3, 4, 5;
	//data.resize(2,5);
	//data << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
	//double a[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ,11,12,13,14,15,16,17,18,19,20};
	//vector<double> testVec1(a,a+3);
	//vector<double> testVec2(a+5, a + 10);
	//vector<double> testVec3(a + 10, a + 15);
	//vector<double> testVec4(a + 15, a + 20);
	//vector<vector<double>> testVecs;
	//testVecs.push_back(testVec2);
	//testVecs.push_back(testVec3);
	//testVecs.push_back(testVec4);
	//linearLeastSquaresStd<double, 5> test;
	//test.addInputData(data);
	//test.printPHI();
	//test.addInputData(testVecs);
	//test.printPHI();
	//test.addOutputData(dataY);
	//test.optimize();
	//test.printTheta();
	//test.addOutputData(testVec1);
	//test.optimize();
	//test.printTheta();
	//linearLeastSquaresStd<double, 2> test2;

	//vector<double> YTest;
	//for (double x = 0; x < 10;x+=0.1) {
	//	YTest.push_back(3.5*x+1.2);
	//}
	//vector<vector<double>> PHITest = zlzLS::linearLeastSquaresStd<double,5>::readFromFile(".\\test.txt", 1, 2);
	//test2.addInputData(PHITest);
	////test2.printPHI();
	//test2.addOutputData(YTest);
	//test2.optimize();
	//
	//cout << test2<<endl;
	//discreteLinearSquarse<double, 4, 2> test1;
	//vector<double> inputTest = discreteLinearSquarse<double, 4, 2>::readFromFile(".\\dataGenerator.txt",1);
	//vector<double> outputTest = discreteLinearSquarse<double, 4, 2>::readFromFile(".\\dataGenerator.txt", 2);
	//test1.addInputData(inputTest);
	//test1.addOutputData(outputTest);
	//test1.optimize();
	//cout << test1 << endl;
	recursiveLeastSquaresStd<double, 2> test3;
	test3.setRecordTheta(true);
	double YTest3;
	Matrix<double, 2, 1> phi3(10000,0);
	for (double x = 15; x < 1000;x+=0.1) {
		YTest3 = 3.5 * x + 1.234;
		phi3(0, 0) = x;
		phi3(1, 0) = 1;
		test3.step(phi3,YTest3);
	}
	ofstream fAssociation;
	fAssociation.open("outTest3.txt");
	test3.writeThetaBuffer(fAssociation);
	fAssociation.close();
}

