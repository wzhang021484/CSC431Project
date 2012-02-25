#pragma once

#include <iostream>
#include <vector>
using namespace std;

#define datatype float // Please change this line if you want to use other data type, such as int, double, etc.

const double PRECISION=1e-6;

class matrix
{
public:
	int rows;
	int cols;
	vector<datatype> data;
	matrix(void);
	matrix(int rows);
	matrix(int rows, int cols);
	~matrix(void);

	datatype operator()(int i, int j=0) const 
	{
		return data[i*cols+j];
	}
	datatype &operator()(int i, int j=0) 
	{
		return data[i*cols+j];
	}
	void swap_rows(int i, int j) 
	{
		if(i<0 || i>=rows || j<0 || j>=rows) 
		{
			cout << "OutOfBounds" << endl;
			exit(-1);
		}
		datatype tmp;
		for(int c=0; c<cols; c++) 
		{
			tmp=(*this)(i,c); 
			(*this)(i,c)=(*this)(j,c); 
			(*this)(j,c)=tmp;
		}
	}
};

// operators on matrix:
ostream &operator<<(ostream &out, const matrix &A);
matrix operator+(const matrix &A, const matrix &B);
matrix operator-(const matrix &A, const matrix &B);
matrix operator*(datatype a, const matrix &B);
matrix operator*(const matrix &A, const matrix &B);
matrix operator/(const matrix &A, double a);
bool operator!=(const matrix& A, const matrix& B);
bool operator==(const matrix& A, const matrix& B);

// function on matrix
matrix inv(matrix A);

double norm_1(const matrix& A);
double norm_2(const matrix& A);
double norm_infinite(const matrix& A);

bool is_almost_symmetric(matrix A, double ap = 0.000001, double rp = 0.0001);
bool is_almost_zero(matrix A, double ap = 0.000001, double rp = 0.0001);

double condition_number(const matrix& A);
matrix trans(const matrix& A);
matrix my_identity(int n);

matrix exp(matrix x, double ap=1e-6, double rp=1e-4, int ns=40);

matrix Cholesky(const matrix& A);

bool is_positive(const matrix& A);

// numeric functions
matrix Markoviz(matrix mu, const matrix& A, double r_free);
matrix fit_least_squares(const matrix& points, int n);