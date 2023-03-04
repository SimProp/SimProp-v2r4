
// File   : PhysicsFunctions.h
//
// Purpose: Physics propagate nuclei
//
// Authors: Roberto Aloisio (fortran), 
//          Denise Boncioli, Sergio Petrera (C++), 
//          Francesco Salamida, Armando di Matteo
// 
// date: 21/07/2011
// Last update: 17/02/2016


#ifndef _PhysicsFunctions_hh_
#define _PhysicsFunctions_hh_

double derivsN(double, double);
double evolveN (double, double, double);
double rr(double );
double GetR(double, double);
double W (double, double, double);
//double GetTtot(double&, double&,double&);
//double Getdtdzeff(double&, double&, int&, int&, double&, double&, double&, double&);
//double GetTheta(double&, double&,int&, int&, double&, double&, double&);
double phi(double );
double phi_inv(double );
double sample_s(double , double , double& );
double sample_eps(double );
//double OptDepth(double );
double integCMB(double, double);

#endif
