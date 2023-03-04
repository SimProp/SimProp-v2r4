// File   : Constants.h
//
// Define cosmology & injection constants


#ifndef _Constants_h_
#define _Constants_h_
//.........................................................
  const double om_m=0.3;            // Omega matter
  const double om_l=0.7;            // Omega lambda
  const double H0=7.17e-11;         // Hubble constant yr^{-1} (H0=70 Km/Mpc/s)
  const double mN=938.27e6;         // Proton mass (eV)
  const double mpi=134.977e6;       // Pion mass (eV)
  const double mmu=105.658e6;       // Muon mass (eV)
  const double me=.511e6;           // Electron mass (eV)
  const double epsthr = (mpi + (mpi*mpi)/(2*mN))/1e6; // Pion production threshold (MeV)
  const double Mpc2cm = 3.086e24;   // Mpc -> cm
  const double SpeedOfLight = 3.e8; // speed of light m s^{-1}
  const double m2cm = 1.e2;         // m -> cm
  const double SecInY = 3.1557e7;   // year -> seconds
//.........................................................
#endif
