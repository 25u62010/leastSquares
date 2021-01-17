#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;
namespace zlzLS {
template<typename _Scalar>
class fileOperations
{
public:
	fileOperations(string filePath);
	~fileOperations();
	inline bool isOpen(){
		return _fileIsOpen;
	}
	inline bool isEnd() {
		return _fAssociation.eof();
	}
	vector<_Scalar> readLine(int startCol=1, int endCol=-1);
	
	static bool readFromFile(const string &filename, vector<vector<_Scalar>>& result,int startCol=1, int endCol=-1, int startRow = 1, int endRow = -1);
	static bool readFromFile(const string &filename, vector<_Scalar>& result, int col, int startRow = 1, int endRow = -1);
	
private:
	bool _fileIsOpen;
	string _filePath;
	ifstream _fAssociation;
};
template<typename _Scalar>
fileOperations<_Scalar>::fileOperations(string filePath):_filePath(filePath){
	_fAssociation.open(_filePath.c_str());
	_fileIsOpen = _fAssociation.is_open();
	if (!_fileIsOpen) {
		cout << "cann't open file" << _filePath << endl;
	}
}
template<typename _Scalar>
fileOperations<_Scalar>::~fileOperations(){
	_fAssociation.close();
}
template<typename _Scalar>
vector<_Scalar> fileOperations<_Scalar>::readLine(int startCol , int endCol) {
	vector<_Scalar> result;
	string s;
	if (!_fAssociation.eof()) {
		getline(_fAssociation, s);
		stringstream ss;
		ss << s;
		_Scalar val;
		int col = 1;

		while (ss >> val) {
			if (col >= startCol && (col <= endCol || col == -1)) {
				result.push_back(val);
			}
			col++;
		}
	}
	return result;
}
template<typename _Scalar>
bool fileOperations<_Scalar>::readFromFile(const string &filename, vector<vector<_Scalar>>& result , int startCol, int endCol, int startRow, int endRow) {
	ifstream fAssociation;
	fAssociation.open(filename.c_str());
	int row = 1;
	int size = 0;
	string s;
	if (!fAssociation.is_open()) {
		cout << "cann't open file" << filename << endl;
		return false;
	}
	while (!fAssociation.eof()) {
		getline(fAssociation, s);
		if (!s.empty() && row >= startRow && (row <= endRow || endRow == -1)) {
			stringstream ss;
			ss << s;
			_Scalar val;
			int col = 1;
			result.push_back(vector<_Scalar>());

			while (ss >> val) {
				if (col >= startCol && (col <= endCol||col==-1)) {
					result[size].push_back(val);
				}
				col++;
			}
			size++;
		}
		row++;
	}
	fAssociation.close();
	return true;
}
template<typename _Scalar>
bool fileOperations<_Scalar>::readFromFile(const string &filename, vector<_Scalar>& result, int col, int startRow, int endRow) {
	ifstream fAssociation;
	fAssociation.open(filename.c_str());
	int row = 1;
	string s;
	if (!fAssociation.is_open()) {
		cout << "cann't open file" << filename << endl;
		return false;
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
					break;
				}
				nowCol++;
			}
		}
		row++;
	}
	fAssociation.close();
	return true;
}

}
