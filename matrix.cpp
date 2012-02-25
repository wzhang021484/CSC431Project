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

matrix::matrix(int _rows)
{
	matrix(_rows, _rows);
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

matrix operator/(const matrix &A, double a)
{
	matrix C(A.rows,A.cols);
	for(int r=0; r<A.rows; r++)
		for(int c=0; c<A.cols; c++)
			C(r,c)=A(r,c)/a;
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

bool operator!=(const matrix& A, const matrix& B) 
{
    return (norm_1(A-B)>PRECISION);
}

bool operator==(const matrix& A, const matrix& B) 
{
    return (norm_1(A-B)<PRECISION);
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

double norm_1(const matrix& A) 
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

double norm_2(const matrix& A)
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

double norm_infinite(const matrix& A)
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

double condition_number(const matrix& A) 
{
	return norm_1(A)/norm_1(inv(A));
}

matrix trans(const matrix& A) 
{
	matrix B(A.cols,A.rows);
	for(int i=0; i<A.rows; i++)
	{
		for(int j=0; j<A.cols; j++)
		{
			B(j,i)=A(i,j);
		}
	}
	return B;
}

matrix my_identity(int n) 
{
	matrix A(n, n);
	for (int i = 0; i < n; i++) 
	{
		A(i, i) = 1;
	}

	return A;
}

matrix exp(matrix x, double ap, double rp,int ns)
{
	matrix t = my_identity(x.cols);
	matrix s = my_identity(x.cols);

	for (int k=0; k<ns; k++)
	{
		t = t*x/k;   // next term
		s = s + t;   // add next term
		if (norm_1(t)<max(ap,norm_1(s)*rp))
		{
			return s;
		}
	}

	cout << "error! exp: no convergence!\n";
	exit(-1);
}

matrix Cholesky(const matrix& A) 
{
	if(A!=trans(A)) 
	{
		cout << "error! Cholesky: not symmetric!\n";
		exit(-1);
	}
	matrix L=A;
	for(int k=0; k<L.cols; k++) 
	{
		if(L(k,k)<=0)
		{
			cout << "error! Cholesky: not positive!\n";
			exit(-1);
		}
		L(k,k)=sqrt(L(k,k));
		for(int i=k+1; i<L.rows; i++)
		{
			L(i,k)/=L(k,k);
		}
		for(int j=k+1; j<L.rows; j++)
		{
			for(int i=k+1; i<L.rows; i++)
			{
				L(i,j)-=L(i,k)*L(j,k);
			}
		}

	}
	for(int i=0; i<L.rows; i++) 
	{
		for(int j=i+1; j<L.cols; j++)
		{
			L(i,j)=0;
		}
	}
	return L;	
}

bool is_positive(const matrix& A)
{
	Cholesky(A);
	return true;
}

matrix Markoviz(matrix mu, const matrix& A, double r_free) 
{
	matrix x(A.rows,1);
	for(int r=0; r<mu.rows; r++)
	{
		mu(r,0)-=r_free;
	}
	x=inv(A)*mu;
	float x_norm=0;
	for(int r=0; r<mu.rows; r++)
	{
		x_norm+=x(r,0);
	}
	for(int r=0; r<mu.rows; r++)
	{
		x(r,0)/=x_norm;
	}

	return x;
}

matrix fit_least_squares(const matrix& points, int n)
{
	matrix b(points.rows);
	matrix A(points.rows,n);
	for(int i=0; i<points.rows; i++) 
	{
		for(int j=0; j<n; j++)
		{
			A(i,j)=pow(points(i,0),j);
		}
		b(i)=points(i,1);
	}
	matrix x=inv(trans(A)*A)*trans(A)*b;

	return x;
}