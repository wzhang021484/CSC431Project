#pragma once

#include "matrix.h"

const double EPSILON=1e-3;

class numeric 
{
public:
    virtual double f(double)=NULL;
    virtual double g(double x)
	{
        return f(x)+x;
    }
    double d1f(double x) 
	{
        double h=EPSILON;
        return (f(x+h)-f(x))/h;
    }
    double d2f(double x) 
	{
        double h=EPSILON;
        return (f(x+h)-2.0*f(x)+f(x-h))/(h*h);
    }
    double d1g(double x) 
	{
        double h=EPSILON;
        return (g(x+h)-g(x))/h;
    }   

	double solve_fixed_point(double x_guess);
};