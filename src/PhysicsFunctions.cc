#include <iomanip>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <list>
#include <cmath>
#include <TRandom3.h>

#include <MathFunctions.h>
#include <PhysicsFunctions.h>
#include <Constants.h>
#include "lossesN.h"
#include <sigmaA_disi.h>
#include <NucModel.h>

using namespace std;

// -------------------------------------------------------------------------

double derivsN (double z, double LnEg)
{
  double lbeta_pg = 0.;
  double lEg = LnEg * log10(exp(1.)) + log10(1.+z); 
  // double * y2a = new double[nE_p];
  double beta_pg = 0.;
  double dtdz = 1. / (H0 * (1.+z) * sqrt(om_m * pow(1+z,3) + om_l));
  // double dtdz = dtdz_com;
  double dLnEgdz = 1. / (1.+z); 
  if (lEg >= 17.) {
    static double *y2a = NULL;
    if (!y2a) {
      y2a = new double[nE_p];
      spline(LogE,lbeta0,nE_p,1.e30,1.e30,y2a);
    }
    lbeta_pg = splint(LogE,lbeta0,y2a,nE_p,lEg);
    beta_pg = pow(1+z,3) * pow(10, lbeta_pg);
    dLnEgdz += dtdz * beta_pg;
  }
  return dLnEgdz;
}


double evolveN (double z_in, double Gam_in, double z_out)
{
  // .............................................................................
  //  Proton evolution according to the energy losses of BGG 2002, the generation 
  //  energy of protons is determined for any intitial condition (z_in,Gam_in).
  // .............................................................................

  static const double eps=1.e-9;
  double h1 =1.e-4;
  double hmin = 0.;

  int nok, nbad;
  
  double LnEg0 = log(Gam_in * mN);
  double LnEg = odeint(LnEg0, z_in, z_out, eps, h1, hmin, nok, nbad, derivsN);
  double Gam_out = exp(LnEg)/mN;
  return Gam_out;
}

double rr(double z)
{
  double ddz = 1./(H0 * sqrt(om_m * pow(1+z,3) + om_l));
  return SpeedOfLight*m2cm*SecInY*ddz/Mpc2cm;
}


double GetR(double zmin, double z)
{	
  return qgaus(rr,zmin,z);
}

//double GetTheta(double &zi, double &zf, int &A, int &Ze, double & G, double & b, double &lc)
//{
//  double zmin = zi;
//  double zmax = zf;
//  double L = qgaus(rr,zmin,zmax);
//  double Rl=108*(mN*A*G/1e20)/Ze/b;
//  //  double al=1./3.;
//  double al=1.;
////  double Ld=lc*(pow((Rl/lc),2)+pow((Rl/lc),al));
//  double Ld=lc*(pow((Rl/lc),2));
//  double angle = sqrt(L/Ld);
//
//  return angle;
//}

//double Getdtdzeff(double &zi, double &zf, int &A, int &Ze, double & G, double &b, double &lc, double &t)
//{
//  double zmin = zi;
//  double zmax = zf;
//  double L = qgaus(rr,zmin,zmax);
//  double Lnew;
//  double Ld=L/t/t;
//  if(b<1e-10) {Lnew=L;}
//  else
//      {      double  Ln=Ld*((exp(1*t*t)-1));
//        Lnew=Ln;}
//
//  double diff=0;
//   diff = Lnew*Mpc2cm/((1.+zmin)*(zmax-zmin)*SpeedOfLight*SecInY*m2cm);
//  return diff;
//}




// -------------------------------------------------------------------------

double W (double eth, double e0, double d0)
{
     
  double pi = 3.1416;

  // here we corrected a typo in eq. (3) of the Stecker & Salamon paper

  double x1 = (e_1-e0)/(d0/sqrt(2.));     
  double x2 = (e0-eth)/(d0/sqrt(2.));     

  double s_erfx1 = erf(x1);
  double s_erfx2 = erf(x2);

  double Wval = d0 * sqrt(pi/8.) * (s_erfx1 + s_erfx2);

  return Wval;
}

const int n_phi = 54;
const double x_phi[] = { //s in GeV^2
1.165, 1.17, 1.175, 1.18, 1.185, 1.19, 1.195, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3., 3.5, 4., 4.5, 5.
};
const double y_phi[] = { //int (s-m**2)*sigma(s)ds in GeV^4 mubarn
0., 0.00421300986, 0.020118723, 0.0492934866, 0.0923954114, 0.14975072, 0.221511203, 0.307728178, 0.954606438, 1.95783016, 3.32413408, 5.08585576, 7.31149413, 10.1106731, 13.6257831, 18.0064959, 23.3748365, 29.791079, 37.2159201, 45.3659542, 53.5104124, 61.2190907, 68.2301781, 74.4674277, 79.9939617, 84.9744172, 89.6221849, 94.0209122, 98.2327261, 102.309672, 106.297433, 110.236213, 114.167759, 118.267121, 122.62728, 127.260338, 132.173961, 137.396617, 142.957783, 148.888999, 175.334789, 207.105032, 246.15604, 286.323605, 321.939067, 356.319849, 392.997575, 433.403795, 475.513311, 516.204338, 711.311197, 931.899478, 1173.89448, 1439.60659
}; //data obtained by integrating SOPHIA output

double phi(double s)
// Integral from s_th to s_max of (s-m**2)*sigma(s) ds
// s in eV^2, sigma in mbarn, phi in eV^4 mbarn
// (eq. 9 in Protheroe and Johnson (1995), Astropart. Phys. 4 253-269)
{
      static double *y2a = NULL;
      double retval;
      if (!y2a) {
        y2a = new double[n_phi];
        spline(x_phi,y_phi,n_phi,1.e30,1.e30,y2a);
      }
      s /= 1.e18; // convert to GeV^2
      if (s < 1.165)
          retval = 0.;
      else  if (s < 5.)
          retval = splint(x_phi,y_phi,y2a,n_phi,s);
      else
          retval = 56.6*s*s;
      return retval*1.e33; // convert to eV^4 mbarn
}

double phi_inv(double p)
// Inverse function of the above; p in eV^4 mbarn, phi_inv in eV^2
{
      static double *y2a = NULL;
      double retval;
      if (!y2a) {
        y2a = new double[n_phi];
        spline(y_phi,x_phi,n_phi,1.e30,1.e30,y2a);
      }
      p /= 1.e33; // convert to GeV^4 mubarn
      if (p < 1439.)
          retval = splint(y_phi,x_phi,y2a,n_phi,p);
      else
          retval = sqrt(p/56.6);
      //cout << "# phi()_inv: p = " << p << ", phi_inv = " << retval << endl;
      return retval * 1.e18; // convert to eV^2
}

double sample_s(double eps, double E, double& angle)
// Given photon energy eps and nucleon energy E, samples centre-of-mass energy squared s
// See eq. 19 in Muecke et al. (1999), arXiv:astro-ph/9903478v1
// or eq. 10 in Protheroe and Johnson (1995), Astropart. Phys. 4 253-269.
// Energies are in eV, output is in eV^2, angle is in radians
{
        double r = gRandom->Uniform(0.,1.);
        double s = phi_inv(r*phi(mN*mN + 4.*eps*E));
        // phi(s)/phi(s_max) is the cumulative probability distribution
        // so the sampled s is such that r = phi(s)/phi(s_max), i.e. phi(s) = r*phi(s_max)
        double cosangle = (mN*mN - s)/(2.*E*eps)+1.;
        if (cosangle < -1.) // can happen due to rounding errors
                cosangle = -1.;
        if (cosangle > 1.)
                cosangle = 1.;
        angle = acos(cosangle);
        // the angle between the nucleon and the photon in the lab frame        
        return s;
}

double sample_eps(double E)
// samples photon energy given nucleon energy
// E and eps in eV
{
        double u1, u2;
        double kT = .2348e-3;
        double eps_min = 6.788e16/E; // 6.788e16 = (2*m_pi*m_N + m_pi**2)/4 
        double eps_max = 0.019; // 0.007*(T/K) eV, see SOPHIA
        double eps_pmax = (3.e-3*pow(E*kT*1e-18,-0.97)+0.047)/3.9e2*2.725; //see SOPHIA
        double p_max = phi(mN*mN + 4.*eps_pmax*E)/(exp(eps_pmax/kT)-1.);
        do {
                u1 = gRandom->Uniform(eps_min,eps_max);
                u2 = gRandom->Uniform(0.,p_max);
        //cout << "sample_eps " << eps_min << " " << u1 << " " << eps_max << " " << eps_pmax << " " << u2 << " " << p_max << endl;        
        } while ( u2 > phi(mN*mN + 4.*u1*E)/(exp(u1/kT)-1.) );
        return u1;
}


//double OptDepth(double E)
//// Armando di Matteo, 27 Feb 2013
//// Returns the redshift from which a photon with energy E
//// has probability exp(-1) of reaching Earth
//{
//  static const int nOpt =  27;
//  static const double EOpt[nOpt] = {
//    10.77, 11.12, 11.60, 12.03, 12.35, 12.81, 13.14, 13.59, 14.04, 14.53,
//    15.00, 15.36, 15.75, 16.41, 17.45, 18.26, 18.77, 19.10, 19.32, 19.60,
//    19.74, 19.95, 20.15, 20.60, 21.01, 21.50, 21.78
//  };
//  static const double zOpt[nOpt] = {
//    +0.36, -0.05, -0.62, -1.02, -1.16, -1.37, -1.74, -2.65, -3.28, -4.92,
//    -5.66, -5.76, -5.67, -5.31, -4.53, -3.86, -3.44, -3.21, -3.11, -3.04,
//    -3.06, -3.06, -3.04, -2.92, -2.72, -2.43, -2.26    
//  }; // data from arXiv:1302.6460v1 [astro-ph.HE]
//  static double *y2a = NULL;
//  if (!y2a) {
//    y2a = new double[nOpt];
//    spline(EOpt,zOpt,nOpt,1.e30,1.e30,y2a);
//  }
//  return pow(10., splint(EOpt,zOpt,y2a,nOpt,log10(E)));
//}


double integCMB(double eps, double z)
// integral from eps to +inf of n_CMB(e)/(2e^2) de
{
  const double kT = 2.327e-4 * (1+z);    // in eV
  const double norm = 1.452e12 * (1+z);  // kT/(2 pi^2) 
  return -norm*log(1-exp(-eps/kT));       // in 1/(y MeV^2 mb)
}

