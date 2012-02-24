#include "matrix.h"


matrix::matrix(void)
{
	rows=0;
	cols=0;
	data.resize(0);
}

matrix::~matrix(void)
{
	data.clear();
}

matrix::matrix(int _rows, int _cols)
{
	rows=_rows;
	cols=_cols;
	data.resize(rows*cols);
	for(int r=0; r<rows; r++)
		for(int c=0; c<cols; c++)
			data[r*cols+c]=0;
}

ostream &operator<<(ostream &out, const matrix& A)
{
	out << "[";
	for(int r=0; r<A.rows; r++) 
	{
		if(r>0) out << ",";
		out << "[";
		for(int c=0; c<A.cols; c++) 
		{
			if(c>0) out << ",";
			out << A(r,c);    
		}
		out << "]";
	}
	out << "]";    
	return out;
}

matrix operator+(const matrix &A, const matrix &B) 
{
	if(A.cols!=B.rows)
	{
		cout << "error! dimension is not the same!\n";
		exit(-1);
	}
	matrix C(A.rows,A.cols);
	for(int r=0; r<A.rows; r++)
		for(int c=0; c<A.cols; c++)
			C(r,c)=A(r,c)+B(r,c);
	return C;
}

matrix operator-(const matrix &A, const matrix &B)
{
	if(A.cols!=B.rows)
	{
		cout << "error! dimension is not the same!\n";
		exit(-1);
	}
	matrix C(A.rows,A.cols);
	for(int r=0; r<A.rows; r++)
		for(int c=0; c<A.cols; c++)
			C(r,c)=A(r,c)-B(r,c);
	return C;
}

matrix operator*(datatype a, const matrix &B)
{
	matrix C(B.rows,B.cols);
	for(int r=0; r<B.rows; r++)
		for(int c=0; c<B.cols; c++)
			C(r,c)=a*B(r,c);
	return C;
}

matrix operator*(const matrix &A, const matrix &B)
{
	if(A.cols!=B.rows)
	{
		cout << "error! dimension is not the same!\n";
		exit(-1);
	}
	matrix C(A.rows,B.cols);
	for(int r=0; r<A.rows; r++)
		for(int c=0; c<B.cols; c++)
			for(int k=0; k<A.cols; k++)
				C(r,c)+=A(r,k)*B(k,c);
	return C;
}

matrix inv(matrix A) 
{
	if(A.cols!=A.rows)
	{
		cout << "error! dim1 and dim2 is not the same!\n";
		exit(-1);
	}
	matrix B(A.rows,A.cols);
	datatype p;
	datatype q;
	int m;
	for(int r=0; r<B.cols;r++) B(r,r)=1;
	for(int c=0; c<A.cols;c++) 
	{    
		m=c; 
		p=A(c,c);
		for(int i=c+1; i<A.rows; i++)
		{
			if(abs(A(i,c)) > abs(p)) 
			{
				m=i; 
				p=A(i,c);
			}
		}

		/*for(int i=0; i<A.cols; i++) 
		{
		    swap(A(m,i),A(c,i));
		    swap(B(m,i),B(c,i));
		}*/
		A.swap_rows(m, c);
		B.swap_rows(m, c);

		for(int i=0; i<A.cols; i++) 
		{
			A(c,i) /= p; 
			B(c,i) /= p;
		}
		for(int r=0; r<A.rows; r++) 
		{
			if(r!=c)
			{
				q = A(r,c);
				for(int i=0; i<A.cols; i++)
				{
					A(r,i)-=q*A(c,i);
					B(r,i)-=q*B(c,i);
				}
			}
		}
	}
	return B;
}

datatype norm_1(const matrix& A) 
{
	datatype z,m=0;
	for(int j=0; j<A.cols; j++) 
	{
		z=0;
		for(int i=0; i<A.rows; i++) 
		{
			z+=abs(A(i,j));
		}
		if(z>m) 
		{
			m=z;
		}
	}
	return m;
}

datatype norm_2(const matrix& A)
{
	if(A.cols!=1) 
	{
		cout << "error! function norm_2 is only for vectors!\n";
		exit(-1);			
	}
	datatype m=0;
	for(int i=0; i<A.rows; i++)
	{
		m+=A(i)*A(i);
	}
	return sqrt(m);
}

datatype norm_infinite(const matrix& A)
{
	datatype z,m=0;
	for(int j=0; j<A.rows; j++) 
	{
		z=0;
		for(int i=0; i<A.cols; i++) 
		{
			z+=abs(A(j,i));
		}
		if(z>m)
		{
			m=z;
		}
	}
	return m;
}

bool is_almost_symmetric(matrix A, double ap, double rp)       
{
	if (A.rows != A.cols)
	{
		return false;
	}
	else
	{
		for (int r = 0; r < A.rows; r++)
		{
			for (int c = 0; c < r; c++)
			{                       
				double delta = abs(A(r, c) - A(c, r));                      
				double abs_arc = abs(A(r, c));
				double abs_acr = abs(A(c, r));                       
				if ((delta > ap) && (delta > max(abs_arc, abs_acr) * rp))
				{
					return false;
				}
			}
		}
	}

	return true;
}

bool is_almost_zero(matrix A, double ap, double rp)
{
	bool result = true;

	for (int r = 0; r < A.rows; r++)
	{
		for (int c = 0; c < A.rows; c++)
		{
			double delta = abs(A(r, c) - A(c, r));
			double abs_arc = abs(A(r, c));
			double abs_acr = abs(A(c, r));
			if ((delta > ap) && (delta > max(abs_arc, abs_acr) * rp))
			{
				result = false;
			}
		}
	}

	return result;
}