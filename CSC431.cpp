// CSC431.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "matrix.h"

void test_matrix();

int _tmain(int argc, _TCHAR* argv[])
{
	test_matrix();

	return 0;
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

