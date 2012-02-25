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

double numeric::solve_bisection(double a, double b)
{
	double fa=f(a);
	double fb=f(b);
	double x, fx;
	if(fa==0) return a;
	if(fb==0) return b;
	if(fa*fb>0) 
	{
		cout << "error! solve_bisection: f(a) and f(b) must have opposite sign!\n";
		exit(-1);
	}

	for(int k=0; k<NS; k++) 
	{          
		x=(a+b)/2;
		fx=f(x);
		if(abs(fx)<PRECISION)
		{
			return x;
		}
		else if(fx*fa<0)
		{
			b=x; fb=fx;
		} 
		else if(fx*fb<0) 
		{
			a=x; fa=fx;
		}
	}
	return x;
}

double numeric::solve_newton(double x)
{
	double x_old=x+PRECISION;
	double f1x;
	for(int k=0; k<NS; k++) 
	{
		cout << x << endl;
		f1x=d1f(x);	
		if(abs(f1x)<PRECISION) 
		{
			cout << "error! solve_newton: Instability!\n";
			exit(-1);
		}
		x_old=x;
		x=x-f(x)/f1x;
		if(abs(x-x_old)<PRECISION) return x;
	}

	cout << "error! exp: no convergence!\n";
	exit(-1);
}

double numeric::solve_secant(double x) 
{
	double x_old;
	double fx, f1x, f_old;
	x_old=x-0.0001;
	f_old=f(x_old);
	for(int k=0; k<NS; k++)
	{
		fx=f(x);
		f1x=(fx-f_old)/(x-x_old);	
		if(abs(f1x)<PRECISION) 
		{
			cout << "error! solve_newton: Instability!\n";
			exit(-1);
		}
		f_old=fx;
		x_old=x;
		x=x-fx/f1x;
		if((k>1) && (abs(x-x_old)<PRECISION))
		{
			return x;
		}
	}
	cout << "error! exp: no convergence!\n";
	exit(-1);
}

double numeric::solve_newton_stabalized(double a, double b)
{
	double fa=f(a);
	double fb=f(b);
	double x, fx;
	if(fa==0) return a;
	if(fb==0) return b;
	if(fa*fb>0)
	{
		cout << "error! solve_bisection: f(a) and f(b) must have opposite sign!\n";
		exit(-1);
	}
	double f1x=0;
	for(int k=0; k<NS; k++) 
	{		
		if(abs(f1x)>PRECISION) 
		{
			x=x-fx/f1x;
		}
		if((abs(f1x)<=PRECISION) || (x<=a) || (x>=b)) 
		{
			x=(a+b)/2;
		}
		fx=f(x);
		f1x=d1f(x);
		if(abs(fx)<PRECISION)
		{
			return x;
		} else if(fx*fa<0)
		{
			b=x; fb=fx;
		} else if(fx*fb<0)
		{
			a=x; fa=fx;
		}
	}
	return x;
}