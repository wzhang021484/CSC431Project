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
			cout << "error! solve_secant: Instability!\n";
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
	cout << "error! solve_secant: no convergence!\n";
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
		cout << "error! solve_newton_stabalized: f(a) and f(b) must have opposite sign!\n";
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

double numeric::optimize_bisection(double a, double b)
{
	double fa=d1f(a);
	double fb=d1f(b);
	double x, fx;
	if(fa==0)
	{
		return a;
	}
	if(fb==0)
	{
		return b;
	}
	if(fa*fb>0)
	{
		cout << "error! optimize_bisection: f(a) and f(b) must have opposite sign!\n";
		exit(-1);
	}

	for(int k=0; k<NS; k++)
	{           
		x=(a+b)/2;
		fx=d1f(x);
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

double numeric::optimize_newton(double x)
{
	double x_old=x+PRECISION;
	double f2x;
	for(int k=0; k<NS; k++)
	{            
		f2x=d2f(x);	
		if(abs(f2x)<PRECISION)
		{
			cout << "error! optimize_newton: Instability!\n";
			exit(-1);
		}
		x_old=x;
		x=x-d1f(x)/f2x;
		if(abs(x-x_old)<PRECISION)
		{
			return x;
		}
	}


	cout << "error! optimize_newton: no convergence!\n";
	exit(-1);
}

double numeric::optimize_secant(double x) 
{
	double x_old;
	double f1x, f2x, f1_old;
	x_old=x-0.0001;
	f1_old=d1f(x_old);
	for(int k=0; k<NS; k++)
	{          
		f1x=d1f(x);
		f2x=(f1x-f1_old)/(x-x_old);	
		if(abs(f2x)<PRECISION) 
		{
			cout << "error! optimize_secant: Instability!\n";
			exit(-1);
		}
		f1_old=f1x;
		x_old=x;
		x=x-f1x/f2x;
		if((k>1) && (abs(x-x_old)<PRECISION))
		{
			return x;
		}
	}


	cout << "error! optimize_secant: no convergence!\n";
	exit(-1);
}

double numeric::optimize_newton_stabalized(double a, double b)
{
	double f1a=d1f(a);
	double f1b=d1f(b);
	double x, f1x;
	if(f1a==0) 
	{
		return a;
	}
	if(f1b==0)
	{
		return b;
	}
	if(f1a*f1b>0)
	{
		cout << "error! optimize_newton_stabalized: f(a) and f(b) must have opposite sign!\n";
		exit(-1);
	}
	double f2x=0;
	for(int k=0; k<NS; k++)
	{		
		if(abs(f1x)>PRECISION) 
		{
			x=x-f1x/f2x;
		}
		if((abs(f2x)<=PRECISION) || (x<=a) || (x>=b)) 
		{
			x=(a+b)/2;
		}
		f1x=d1f(x);
		f2x=d2f(x);
		if(abs(f1x)<PRECISION) 
		{
			return x;
		}
		else if(f1x*f1a<0)
		{
			b=x; f1b=f1x;
		}
		else if(f1x*f1b<0) 
		{
			a=x; f1a=f1x;
		}
	}
	return x;
}

double numeric::optimize_golden_search(double a, double b) 
{
	double t=(sqrt(5.0)-1)/2;
	double x1=a+(1.0-t)*(b-a);
	double x2=a+(t)*(b-a);
	double fa=f(a);
	double fb=f(b);
	double f1=f(x1);
	double f2=f(x2);
	while(abs(b-a)>PRECISION) 
	{           
		if(f1>f2)
		{
			a=x1;
			fa=f1;	
			x1=x2;
			f1=f2;
			x2=a+(t)*(b-a);
			f2=f(x2);
		} 
		else 
		{
			b=x2;
			fb=f2;
			x2=x1;
			f2=f1;
			x1=a+(1.0-t)*(b-a);
			f1=f(x1);
		}	
	}

	return b;
}	