#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;

class Matrix{
private:
	double **data;
	
public:
	int n_row, n_column;
	
	Matrix(const int n_row, const int n_column);
		
	double& operator() (const int n_row, const int n_column);
	
	friend ostream& operator << (ostream &o, Matrix&m);
	
	
};

ostream& operator << (ostream &o, Matrix&m);

#endif