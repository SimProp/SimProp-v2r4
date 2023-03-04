// File   : Particle.hh
// Particle class
// Authors: F. Salamida, A. di Matteo
// 
// date: 05/07/2011
// Last update: 17/02/2016


#ifndef _Particle_hh_
#define _Particle_hh_

#include <cstdio>
#include <list>
#include <vector>
#include <fstream>

#include "Constants.h"

extern double gEinject;
extern int gAinject;
extern double gzinject;

 
enum eParticleType {
   eNull,
   eProton,
   eNeutron,
   eNucleus,
   ePion,
   eMuon,
   ePhoton,
   eNeutrino_e,
   eAntineutrino_e,
   eNeutrino_mu,
   eAntineutrino_mu,
   eElectron
};

class Particle {

 public:
  
  //Constructors
  Particle();
  Particle(int massnum, int charge, int branch, double Eprod, double zprod, eParticleType type);
  Particle(int massnum, int charge, int branch, double Eprod, double Zprod, double Eint, double Zint, eParticleType type);
  ~Particle();
 
  
  //setters
  void SetMassNum(int val)  { fmassnum = val; }
  void SetCharge(int val)   { fcharge = val; }
  void SetBranch(int val)   { fbranch = val; }
  void SetIntMult(int val)  { fIntMult = val; }
  void SetEprod(double val) { fEprod = val; }
  void SetZprod(double val) { fZprod= val;}
  void SetEint(double val)  { fEint = val; }
  void SetZint(double val)  { fZint = val; }
  void SetGprod(double val) { fEprod = val * mN * fmassnum; }
  void SetGint(double val)  { fEint  = val * mN * fmassnum; }
  void SetType(eParticleType val)  { fType = val; }
//  void SetTtot(double val)  { fTtot = val; }
//  void Setdelta_fin(double val)  { fdelta_fin = val; }
//  void SetB(double val)  { fB = val; }
//  void SetLcoh(double val)  { fLcoh = val; }

 //getters
  int GetNeutrons() const { return fmassnum-fcharge; }
  int GetMassNum() const  { return fmassnum; }
  int GetCharge() const   { return fcharge; }
  int GetBranch() const   { return fbranch; }
  int GetIntMult() const  { return fIntMult; }
  double GetEprod() const { return fEprod; }
  double GetZprod() const { return fZprod;}
  double GetEint() const  { return fEint; }
  double GetZint() const  { return fZint; }
  double GetGprod() const { return fEprod/(mN * fmassnum); }
  double GetGint() const  { return fEint/(mN * fmassnum); }
//  double GetTtot()  {return fTtot; }
//  double Getdelta_fin()  { return fdelta_fin ; }
//  double GetB()  { return fB ; }
//  double GetLcoh()  { return fLcoh ; }
  eParticleType GetType() const { return fType; }
  int GetFlavor() const;
  const char* GetName() const;
  int GetID() const;
  
 private:  

 int fmassnum; 
 int fcharge;
 int fbranch;
 int fIntMult;
 double fEprod;
 double fZprod;
 double fEint; 
 double fZint;
// double fTtot;
// double fdelta_fin;
// double fB;
// double fLcoh;
 eParticleType fType;
 
};

#endif  // _Particle_hh_
