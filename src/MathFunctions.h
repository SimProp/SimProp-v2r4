// File   : MathFunctions.h
//
// Purpose: Mathematics propagate nuclei
//
// Authors: Roberto Aloisio (fortran), Sergio Petrera (C++), F.Salamida 
// 
// date: 19/02/2010
// Last update: 19/02/2013


#ifndef _MathFunctions_hh_
#define _MathFunctions_hh_

//double DivDif (const double*, const double*,const int&, const double&, const int&);   
double qgaus (double (*)(double), double, double); 
void spline (const double*, const double*, int, double , double , double*);
double splint (const double*, const double*, double*, int, double );
void splie2(const double*, const double*, const double*, int, int, double*);
double splin2 (const double*, const double*, const double*, double*, int, int, double, double);
double odeint (double &, double, double, double, double, double, int &, int &, double (*)(double, double));
void rkqs(double &, double &, double & , double & , double, double &, double &, double &, double (*)(double, double));
void rkck(double &, double &, double & , double & , double &, double &, double (*)(double, double)); 
void splie2_inv(const double x1a[], const double x2a[],const double ya[], int m, int n, double *y2a);
double splin2_inv(const double x1a[], const double x2a[], const double ya[], double y2a[], int m, int n, double x1, double y);


#endif

