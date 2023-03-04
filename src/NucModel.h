// File   : NucModel.hh
// Nuclear Model class
// Authors: F. Salamida, A. di Matteo
// 
// date: 05/07/2011
// Last update: 14/11/2016


#ifndef _NucModel_hh_
#define _NucModel_hh_

#include <cstdio>
#include <list>
#include <vector>
#include <fstream>
#include <TRandom3.h>
#include <Particle.h>
#include <math.h>

extern int fLosses;
extern double epsIBL_min, epsIBL_max;

extern std::vector<int> fZ;
extern std::vector<int> fA;
extern int k_com;
extern double Gam_com;
extern double z_com;
extern double s_com;
//extern double dtdz_com;

double derivsA0 (double, double);
double integIBL(double, double);
double inv_integIBL(double, double);
double fdisi (double );
double fdisi_IR (double);
double losseA0pair(double, double);
double losseA0(double, double);
void sigma_pion (double, double &);

enum eModelType {
   eNone,
   eOnlyProton,
   eSalomon,
   eGauss,
   eBreitWigner,
   eBreitWigner2,
   eGauss2
 };

enum eProcType {
  eEarth,
  eHadron,
  eHadron_IR,
  eDisi,
  eNucleonProd,
  eAlphaProd,
//  eDecay,
};

class NucModel {

 public:
  
  //Constructors
  NucModel(eModelType model, int LossType, int Stoch, int Bdecay, int PairProd);
  ~NucModel();
   
  //getters
  int GetNElements()    {return fnelements;}
  int GetLossType()    {return fLosses;}
//  std::vector<double> GetMag()  {return fB;}

  //functions
  //void Validate(Particle* );
  int TestPrimary(int );
  Particle GetRandomParticle(double, double);
  std::vector<Particle> PropagateParticle(Particle* );
  Particle* PropagateProton(Particle*);
  std::vector<Particle> PropagateNucleus(Particle*);
  std::vector<Particle> PropagatePion(Particle*);
  std::vector<Particle> PropagatePhoton(Particle*);
  std::vector<Particle> PropagateMuon(Particle*);
  std::vector<Particle> PropagateNeutrino(Particle*);
  std::vector<Particle> PropagateElectron(Particle*);

  eProcType GetProcess(Particle *, std::vector<Particle>*);
  double GetEdetN(Particle *);
  int GetJNuc(int);
  double losseAdisi(int, double, double, int, int & );
  double HadrRate(double, double); 
  double HadrRate_IR(double, double); 
  double NuclRate(int, double, double); 
  double AlphRate(int, double, double); 
  double TotalRate(int, int, int, double, double, eProcType*);
  double evolveA0 (int, double, double, double);
  double BRdisi(int , int );
  void sigma_disi (int, double, double &, double &, double &);
//  void SetMag(std::vector<double> val) {fB = val;}
  double Photopion(double z, double Gam, eProcType proc);
  double sample_sIBL(double, double, double& );
  void DetermProton(Particle* input);
  double GetDecay(int, int);
  int GetBetaDecayStableZ(int);
  double GetDecay(const Particle* p) {return GetDecay( p->GetMassNum(), p->GetCharge() );}
  std::vector<Particle> BetaDecay(Particle*);
private:  

  int binary(int );
  int fnelements;
  eModelType fType;
  int fStoch;
  int fBdecay;
  int fPairProd;
  
  double fe_1;
  double fe_max;
  std::vector<double> feth_1;
  std::vector<double> feth_2;
  std::vector<double> fe0_1;
  std::vector<double> fcsi_1;
  std::vector<double> fdelta_1;
  std::vector<double> fe0_2;
  std::vector<double> fcsi_2;
  std::vector<double> fdelta_2;
  std::vector<double> fzita;
  std::vector<double> fNuclt;
  std::vector<double> fNuclh1;
  std::vector<double> fNuclx1;
  std::vector<double> fNuclw1;
  std::vector<double> fNuclh2;
  std::vector<double> fNuclx2;
  std::vector<double> fNuclw2;
  std::vector<double> fNuclc;
  std::vector<double> fAlpht;
  std::vector<double> fAlphh1;
  std::vector<double> fAlphx1;
  std::vector<double> fAlphw1;
  std::vector<double> fAlphh2;
  std::vector<double> fAlphx2;
  std::vector<double> fAlphw2;
  std::vector<double> fAlphc;
  
//  std::vector<double> fB;
  
};

extern double qval[32][37];

#endif  // _NucModel_hh_
