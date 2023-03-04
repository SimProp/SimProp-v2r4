// File   : Output.h
// Particle class
// Authors: F. Salamida, A. di Matteo
// 
// date: 05/07/2011
// Last update: 14/11/2016


#ifndef _Output_hh_
#define _Output_hh_

#include <cstdio>
#include <list>
#include <vector>
#include <fstream>

#include "Constants.h"
#include "Particle.h" 

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class Output {

 public:
  
  //Constructors
  Output(char *, int, int);
  ~Output();
  void Stream(const Particle&);
  void Stream(double, int, unsigned long int);
  void Injec(int, const Particle&);
  void Earth(const Particle&);
  void Prod(const Particle&);
  void Write() { if (fsumm) fsumm->Fill(); }
 
 private:  

  TFile* fout;

  TTree* fopt;

  TTree* fnuc;
  int fevt; 
  int fbranch;
  int fintmult;
  int fAcurr;
  int fZecurr;
  int fFlav;
  double fzOri; 
  double fzEnd;
  double fEOri;
  double fEEnd;
  double fDist;
//  double fdelta_final; 
//  double fttot;

  TTree* fev;
  int fbranxev;
  double ftimexev;
  unsigned long int fRandomSeed;
  
  TTree* fsumm;
//  int fID_i;
  double fEnergy_i;   
  double fRedshift_i;
  double fDist_i;
  int fA_i;
  int fZ_i;
//  int fFlav_i;

  unsigned fN_nuc;
  static const unsigned nnuc_max = 100;
  double fE_nuc[nnuc_max];
  int fA_nuc[nnuc_max];
  int fZ_nuc[nnuc_max];

  unsigned fN_pho;
  static const unsigned npho_max = 10000;
  double *fE_pho;
  double *fz_pho;

  unsigned fN_neu;
  static const unsigned nneu_max = 10000;
  double *fE_neu;
  int *fF_neu;

  int fpairs;
  unsigned fN_ele;
  static const unsigned nele_max = 1000000;
  int *fZ_ele;
  double *fE_ele;
  double *fz_ele;
};
  static const double log10z_min = -4., log10z_max = +1., dlog10z = 0.2;
  static const int Nlog10z = 25; //int((log10z_max-log10z_min)/dlog10z + 0.5);

#endif  // _Output_hh_
