//
// File   : MathFunctions.cc
//
// Purpose: Mathematics for propagation.
//

#include <iomanip>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <list>
#include <cmath>

#include <TGraph.h>

#include <MathFunctions.h>
#include <Constants.h>

using namespace std;

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


//double DivDif(const double f[], const double a[], const int &nn,
//                          const double &x, const int &mm)
//{
//  TGraph * dd = new TGraph (nn);
//  for (int i=0; i<nn; i++) dd->SetPoint(i, a[i], f[i]);
//  double value = dd->Eval (x, 0, "S");
//  delete dd;
//  return value;
//}
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
double qgaus (double (*func) (double), double a, double b)  
{   
  const double x[5] = {0.1488743389, 0.4333953941, 0.6794095682,
		       0.8650633666, 0.9739065285};
  const double w[5] = {0.2955242247, 0.2692667193, 0.2190863625,
		       0.1494513491, 0.0666713443};   
  double xm = 0.5 * (a + b);   
  double xr = 0.5 * (b - a);   
  double ss = 0.0;
  double dx;

  for (int j = 0; j<5; j++) {  
    dx = xr * x[j];   
    ss = ss + w[j] * (func(xm + dx) + func(xm - dx));   
  }   
  ss = xr * ss;
  return ss;
}




void spline(const double x[], const double y[], int n, double  yp1, double  ypn, double y2[])
{
  double p,qn,sig,un;
  double* u = new double[n];
                            
  if (yp1 > 0.99e30)
    
    y2[0]=u[0]=0.0;
                                
  else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
                                
  for (int i=1; i<=n-2; i++) {
                                                 
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
                                                
    p=sig*y2[i-1]+2.0;
                                                 
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  
  if (ypn > 0.99e30)
                                    
    qn=un=0.0;
                                
  else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  
  for (int k=n-2; k>=0; k--){
                                   
    y2[k]=y2[k]*y2[k+1]+u[k];
  }
 
  delete [] u;
  return;  
}

 


double splint(const double xa[], const double ya[], double y2a[], int n, double  x)

{
  int k,khi,klo;
  double a,b,h;
  double y;
  
  klo=0;
  khi=n-1;
  while (khi-klo>1) {
    k=(khi+klo) >> 1;
    //k=(khi+klo)/2;
    if(xa[k]>x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h==0.0) cout<< " Bad xa input in routine splint "<<endl;
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  
  return y;
}



void splie2(const double x1a[], const double x2a[],const double ya[], int m, int n, double *y2a)
{
  for (int j=0; j<m; j++)
    spline(x2a,ya+n*j,n,1.0e30,1.0e30,y2a+n*j);;
  return;
}   

                       
double splin2(const double x1a[], const double x2a[], const double ya[], double y2a[], int m, int n, double x1, double x2)                                                                     
{            
  double * yytmp = new double[m];
  double * ytmp = new double[m];

  for (int j=0; j<m; j++) {
    yytmp[j] = splint(x2a,ya+n*j,y2a+n*j,n,x2);
  }
  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
  double y = splint(x1a,yytmp,ytmp,m,x1);

  delete [] yytmp;
  delete [] ytmp;
  return y; 
}

void splie2_inv(const double x1a[], const double x2a[],const double ya[], int m, int n, double *y2a)
// initializes table of 2nd derivatives for splin2_inv
// Armando di Matteo, 30 Apr 2013
{
  for (int j=0; j<m; j++) 
    spline(ya+n*j,x2a,n,1.0e30,1.0e30,y2a+n*j);
  return;
}   

double splin2_inv(const double x1a[], const double x2a[], const double ya[], double y2a[], int m, int n, double x1, double y)                                                                     
// computes x2 as a function of x1 and y
// Armando di Matteo, 30 Apr 2013
{            
  double * x2tmp = new double[m];
  for (int j=0; j<m; j++)
    x2tmp[j] = splint(ya+n*j,x2a,y2a+n*j,n,y);
  double * x22tmp = new double[m];
  spline(x1a,x2tmp,m,1.0e30,1.0e30,x22tmp);
  double x2 = splint(x1a,x2tmp,x22tmp,m,x1);
  delete [] x22tmp;
  delete [] x2tmp;
  return x2; 
}

double odeint(double &ystart, double x1, double x2, double eps, double h1, double hmin, 
		 int &nok, int &nbad, double (*func) (double, double))
// Adapted from Numerical Recipes (Set nvar = 1)
// Runge-Kutta driver with adaptive stepsize control. 
// The routine integrates starting values ystart[0,.nvar-1] from xl to x2 
// with accuracy eps, storing intermediate results in global variables. 
// h1 should be set as a guessedf irst stepsize, hmin as the minimum 
// allowed stepsize (can be zero). On output nok and nbad are the number 
// of good and bad (but retried and fixed) steps taken, and ystart is 
// replaced by values at the end of the integration interval.
// derivs is the user-supplied routine for calculating the right-hand side 
// derivative, while rkqs is the name of the stepper routine to be used.
{
  const int    MAXSTP = 10000;
  const double TINY   = 1.0e-30;
  int nstp;
  double xsav, x, hnext, hdid, h;
  
  double yscal, y, dydx;
  double dxsav;
  int kount;
  const int kmax = 0;
  double xp[1];
  double yp[1];

  x = x1;
  h = SIGN(h1, x2-x1) ;
  
  nok = nbad = kount = 0 ;

  y = ystart ;
  if (kmax > 0) xsav = x - dxsav*2.0;  // Assures storage of first step.
  for (nstp=0; nstp<MAXSTP; nstp++) {
    dydx = func(x, y) ;
      // Scaling used to monitor accuracy. This general-purpose choice can be modified
    // if need be.
    yscal = fabs(y) + fabs (dydx * h) + TINY;
    if ((kmax > 0) && (kount < kmax-1) && (fabs(x-xsav) > fabs(dxsav))) {
      yp[kount] = y;
      xp[kount++] = x; // Store intermediate results.
      xsav = x;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) h= x2-x; // If stepsize can overshoot,decrease.
    rkqs (y, dydx, x, h, eps, yscal, hdid, hnext, func) ;
    if (hdid == h) ++nok; else ++nbad;
    if ((x-x2)*(x2-x1) >= 0.0) { // Are we done?
      ystart = y ;
      if (kmax != 0) {
	yp[kount]= y ;
	xp[kount++] = x; 
      }  // Save final step.
      return y; // Normal exit.
    }

    if (fabs(hnext) <= hmin) cout<<" odeint ==> Step size too small in odeint"<<endl;
    h = hnext;
  }
  cout<<" odeint ==> Too many steps in routine odeint"<<endl;
    
  return 0.;
}

 

void rkqs(double & y, double & dydx, double & x, double & htry, double eps, double  & yscal, double & hdid, double & hnext, double (*func) (double, double))
{
  double h, htemp,xnew;
  double errmax, yerr, ytmp;
  double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;

  h = htry;
 
  for(;;){
    rkck(y,dydx,x,h,ytmp,yerr,func);
    errmax=0.0;
    errmax=max(errmax,fabs(yerr/yscal));
    errmax=errmax/eps;

    if (errmax<=1.0) break;
    
    htemp=SAFETY*h*pow(errmax,PSHRNK);
    h=(h>=0.0 ? max(htemp,0.1*h) : min(htemp,0.1*h));
    xnew=x+h;
    if(xnew==x) cout<<" stepsize underflow in rkqs "<<endl;
  }
 
  if(errmax>ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
  else hnext=5.0*h;

  x+=(hdid=h);
  y=ytmp; 
}

 

void rkck(double & y, double & dydx, double & x, double & h, double & yout, double & yerr, double (*func) (double, double))
{
  
  const double A2=0.2;
  const double A3=0.3;
  const double A4=0.6;
  const double A5=1.0;
  const double A6=0.875;
  const double B21=0.2;
  const double B31=3.0/40.0;
  const double B32=9.0/40.0;
  const double B41=0.3;
  const double B42=-0.9;
  const double B43=1.2;
  const double B51=-11.0/54.0;
  const double B52=2.5;
  const double B53=-70.0/27.0;
  const double B54=35.0/27.0;
  const double B61=1631.0/55296.0;
  const double B62=175.0/512.0;
  const double B63=575.0/13824.0;
  const double B64=44275.0/110592.0;
  const double B65=253.0/4096.0;
  const double C1=37.0/378.0;
  const double C3=250.0/621.0;
  const double C4=125.0/594.0;
  const double C6=512.0/1771.0;
  const double DC1=C1-2825.0/27648.0;
  const double DC3=C3-18575.0/48384.0;
  const double DC4=C4-13525.0/55296.0;
  const double DC5=-277.0/14336.0;
  const double DC6=C6-0.25;

  double ak2,ak3,ak4,ak5,ak6;
  double t21,t22,t31,t32,t41,t42,t51,t52,t61,t62;

  t21=y+B21*h*dydx;
  t22=x+A2*h;
  ak2=func(t22,t21);

  t31=y+h*(B31*dydx+B32*ak2);
  t32=x+A3*h;
  ak3=func(t32,t31);
 
  t41=y+h*(B41*dydx+B42*ak2+B43*ak3);
  t42=x+A4*h;
  ak4=func(t42,t41);
 
  t51=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4);
  t52=x+A5*h;
  ak5=func(t52,t51);
 
  t61=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+
	   B65*ak5);
  t62=x+A6*h;
  ak6=func(t62,t61);
  
 
  yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6);
  
 
  yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*
	  ak6);
}

