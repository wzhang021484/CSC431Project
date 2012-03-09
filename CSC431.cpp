// CSC431.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "matrix.h"

void test_matrix();
void test_Cholesky_and_is_almost_zero();
void test_Markoviz();

int _tmain(int argc, _TCHAR* argv[])
{
	//test_matrix();
	//test_Cholesky_and_is_almost_zero();
	test_Markoviz();

	return 0;
}

void test_Markoviz() 
{
    matrix A(3,3);
    A(0,0)=pow(0.1,2);
    A(1,1)=pow(0.2,2);
    A(2,2)=pow(0.3,2);
    A(0,1)=A(1,0)=0.10*0.20*0.30;
    A(0,2)=A(2,0)=0.25*0.30*0.45;
    A(1,2)=A(2,1)=0.30*0.40*0.20;
    matrix mu(3,1);
    mu(0,0)=0.10;
    mu(1,0)=0.13;
    mu(2,0)=0.16;
    float r_bar=0.05;

    matrix x=Markoviz(mu,A,r_bar);

	cout << "A:" << endl << A << endl;
	cout << "mu:" << endl << mu << endl;
	cout << "r_bar:" << endl << r_bar << endl;
    cout << "x:" << endl << x << endl;
    cout << "return:" << trans(mu)*x << endl;
    cout << "risk:" << sqrt((trans(x)*A*x)(0,0)) << endl;
}

void test_Cholesky_and_is_almost_zero() // also test trans()
{
    matrix A(2,2);
    A(0,0)=0.001; 
    A(1,1)=0.002; 
    A(0,1)=A(1,0)=0.001; 
	cout << "A:" << endl << A << endl;

    matrix L=Cholesky(A);
	cout << "Cholesky(A):" << endl << L << endl;
	cout << "is_almost_zero(A - L*L'), expect true(1): " << is_almost_zero(A - L * trans(L)) << endl;   
}

void test_matrix() 
{
	matrix A(3,3);
	matrix B(3,3);
	A(0,0)=1; A(0,1)=2; A(0,2)=3;
	A(1,0)=-1; A(1,1)=0; A(1,2)=1;
	A(2,0)=4; A(2,1)=2; A(2,2)=7;
	B = inv(A);

	// MAGIC FORMULA
	cout.setf(ios::fixed);
	cout.setf(ios::showpoint);
	cout.precision(2);

	cout << "A:" << endl << A << endl;
	cout << "A*A:" << endl << A*A << endl;
	cout << "inv(A):" << endl << B << endl;
	cout << "inv(A)*A:" << endl << B*A << endl;
	return ;
}

