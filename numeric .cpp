#include "numeric.h"

double numeric::solve_fixed_point(double x_guess) 
{
	double x=x_guess;
	double x_old=x+2.0*PRECISION;
	while(abs(x_old-x)>=PRECISION) 
	{           
		if(abs(d1g(x))>=1) 
		{
			cout << "error! solve_fixed_point: no convergence since D(g)(x)>=1!\n";
			exit(-1);
		}
		x_old=x;
		x=g(x);
	}
	return x;		
}