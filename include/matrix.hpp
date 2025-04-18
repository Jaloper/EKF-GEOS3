#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix {
public:
    int n_row, n_column;
	double **data;

    // Parameterized constructor
    Matrix(const int n_row, const int n_column);
	Matrix(const int n);
	
	// Member operators
	double& operator () (const int row, const int column);
	double& operator () (const int n);
	Matrix& operator + (Matrix &m);
	Matrix& operator - (Matrix &m);
	Matrix& operator * (Matrix &m);
	Matrix& operator / (Matrix &m);
	Matrix& operator = (Matrix &m);
	Matrix& operator + (double d);
	Matrix& operator - (double d);
	Matrix& operator * (double d);
	Matrix& operator / (double d);
	
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int n_row, const int n_column);
Matrix& inv(Matrix &m);
Matrix& eye(int n);
Matrix& transpose(Matrix &m);
Matrix& zeros(const int n);
double norm(Matrix &m);
double dot(Matrix& v1, Matrix& v2);
Matrix& cross(Matrix &v1, Matrix &v2);
Matrix& extract_vector(Matrix &v, int init, int fin);
Matrix& union_vector(Matrix &v1, Matrix &v2);
Matrix& extract_row(Matrix &m, int n);
Matrix& extract_column(Matrix &m, int n);
Matrix& assign_row(Matrix &m, Matrix &row, int n);
Matrix& assign_column(Matrix &m, Matrix &row, int n);


#endif