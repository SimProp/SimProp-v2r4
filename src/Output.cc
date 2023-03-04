#include <stdlib.h>
#include <Output.h>
#include <assert.h>

#include <sstream>
#include <iostream>

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "Particle.h"
#include <PhysicsFunctions.h>
#include "NucModel.h"

using namespace std;
extern int Nevts;
extern unsigned int Seed;
extern int Areq;
extern double Emin;
extern double Emax;
extern double gamma_inj;
extern double zmin;
extern double zmax;
extern double Ls;
extern int LossType;
extern eModelType ModelType;
extern double nFactor;
extern double aFactor;
extern int Bdecay;
extern int Stoch;
extern int PairProd;
//extern double B;
//extern double Lcoh;
//extern double Bs, Rsources;
extern int OutType;

int nucl_model;

Output::Output(char *filename, int type, int pairs) 
    : fopt(NULL), fnuc(NULL), fev(NULL), fsumm(NULL), fpairs(pairs)
{

  fout = new TFile(filename,"RECREATE");
  int comp = 1;
  fout->SetCompressionLevel(comp);

  TDirectory* const save = gDirectory;
  // Spectrum output tree                                                       
  fout->cd();
  
  fopt = new TTree("opt", "command-line options");
  fopt->Branch("Nevts", &Nevts, "Nevts/I");
  fopt->Branch("rand_seed", &Seed, "rand_seed/i");
  fopt->Branch("Ainj", &Areq, "Ainj/I");
  fopt->Branch("Emin", &Emin, "Emin/D");
  fopt->Branch("Emax", &Emax, "Emax/D");
  fopt->Branch("gam_inj", &gamma_inj, "gam_inj/D");
  fopt->Branch("zmin", &zmin, "zmin/D");
  fopt->Branch("zmax", &zmax, "zmax/D");
  fopt->Branch("src_distance", &Ls, "src_distance/D");
  fopt->Branch("EBL_model", &LossType, "EBLmodel/I");
  nucl_model = int(ModelType) - 2;
  fopt->Branch("nucl_model", &nucl_model, "nucl_model/I");
  fopt->Branch("n_factor", &nFactor, "n_factor/D");
  fopt->Branch("a_factor", &aFactor, "a_factor/D");
  fopt->Branch("beta_decay", &Bdecay, "beta_decay/I");
  fopt->Branch("stoch_pion", &Stoch, "stoch_pion/I");
  fopt->Branch("pair_prod", &PairProd, "pair_prod/I");
  fopt->Branch("output_type", &OutType, "output_type/I");
  fopt->Fill();

  if (type == 0 || type == 2) {
    fnuc = new TTree("nuc","Prop Nuclei"); // also contains other particles, the name is for historical reasons
    fnuc->Branch("evt",&fevt,"evt/I");
    fnuc->Branch("branch",&fbranch,"branch/I");
    fnuc->Branch("intmult",&fintmult,"intmult/I");
    fnuc->Branch("Acurr",&fAcurr,"Acurr/I"); 
    fnuc->Branch("Zecurr",&fZecurr,"Zecurr/I");
    fnuc->Branch("Flav",&fFlav,"Flav/I"); 
    fnuc->Branch("zOri",&fzOri,"zOri/D"); 
    fnuc->Branch("zEnd",&fzEnd,"zEnd/D"); 
    fnuc->Branch("EOri",&fEOri,"EOri/D"); 
    fnuc->Branch("EEnd",&fEEnd,"EEnd/D"); 
    fnuc->Branch("Dist",&fDist,"Dist/D"); 
  //  fnuc->Branch("delta_fin",&fdelta_final,"delta_fin/D");
  //  fnuc->Branch("Ttot",&fttot,"Ttot/D");

    fev = new TTree("ev","Events"); 
    fev->Branch("timexev",&ftimexev,"timexev/D");
    fev->Branch("branxev",&fbranxev,"branxev/I");
    fev->Branch("seed",&fRandomSeed,"RandomSeed/i");
  }
  
  if (type == 1 || type == 2) {
    fsumm = new TTree("summary", "Simulated events and products at Earth");
    fsumm->Branch("event",&fevt,"event/I");
//    fsumm->Branch("injID",&fID_i,"injID/I");
    fsumm->Branch("injEnergy",&fEnergy_i,"injEnergy/D");
    fsumm->Branch("injRedshift",&fRedshift_i,"injRedshift/D");
    fsumm->Branch("injDist",&fDist_i,"injDist/D");
    fsumm->Branch("injA",&fA_i,"injA/I");
    fsumm->Branch("injZ",&fZ_i,"injZ/I");
//    fsumm->Branch("injFlav",&fFlav_i,"injFlav/I");

    fsumm->Branch("nNuc",     &fN_nuc,"nNuc/i");
    fsumm->Branch("nucEnergy",fE_nuc,"nucEnergy[nNuc]/D");
    fsumm->Branch("nucA",     fA_nuc,"nucA[nNuc]/I");
    fsumm->Branch("nucZ",     fZ_nuc,"nucZ[nNuc]/I");

    fsumm->Branch("nPho",     &fN_pho,"nPho/i");
    fE_pho = new double[npho_max];
    fsumm->Branch("phoEProd", fE_pho,"phoEProd[nPho]/D");
    fz_pho = new double[npho_max];
    fsumm->Branch("phozProd", fz_pho,"phozProd[nPho]/D");

    fsumm->Branch("nNeu",     &fN_neu,"nNeu/i");
    fE_neu = new double[nneu_max];
    fsumm->Branch("neuEnergy",fE_neu,"neuEnergy[nNeu]/D");
    fF_neu = new int[nneu_max];
    fsumm->Branch("neuFlav",  fF_neu,"neuFlav[nNeu]/I");

    fsumm->Branch("nEle",     &fN_ele,"nEle/i");
    fZ_ele = new int[nele_max];
    fsumm->Branch("eleZ",     fZ_ele,"eleZ[nEle]/I");
    fE_ele = new double[nele_max];
    fsumm->Branch("eleEProd", fE_ele,"eleEProd[nEle]/D");
    fz_ele = new double[nele_max];
    fsumm->Branch("elezProd", fz_ele,"elezProd[nEle]/D");
  }
  
  save->cd();

}


Output::~Output() {

  TDirectory* const save = gDirectory;

  fout->cd();
  fopt->Write();
  if (fnuc) fnuc->Write();
  if (fev) fev->Write();
  if (fsumm) {
    fsumm->Write();
    delete[] fE_pho;
    delete[] fz_pho;
    delete[] fE_neu;
    delete[] fF_neu;
    delete[] fZ_ele;
    delete[] fE_ele;
    delete[] fz_ele;
  }
  fout->Close();
  save->cd();

}


void Output::Stream(const Particle &mypart) {
 if (fnuc) {
  fbranch = mypart.GetBranch();
  fintmult = mypart.GetIntMult();
  fAcurr = mypart.GetMassNum();
  fZecurr = mypart.GetCharge();
  fFlav = mypart.GetFlavor();
  fzOri = mypart.GetZprod();
  fzEnd = mypart.GetZint(); 
  fEOri = mypart.GetEprod(); 
  fEEnd = mypart.GetEint();
  fDist = GetR(fzEnd, fzOri);
//  fttot = mypart.GetTtot(); 
//  fdelta_final = mypart.Getdelta_fin(); 
  fnuc->Fill();
 }
}

void Output::Stream(double timexev, int branxev, unsigned long int cseed) {
 if (fev) {
  fbranxev = branxev;
  ftimexev = timexev;
  fRandomSeed = cseed;
  fev->Fill();
 }
}

void Output::Injec(int evt, const Particle& primary)
{
 fevt = evt;
 if (fsumm) {
//  fID_i = primary.GetID();
  fEnergy_i = primary.GetEprod();
  fRedshift_i = primary.GetZprod();
  fDist_i = GetR(0., fRedshift_i);
  fA_i = primary.GetMassNum();
  fZ_i = primary.GetCharge();
//  fFlav_i = primary.GetFlavor();
  fN_nuc = fN_pho = fN_neu = 0;
  switch (fpairs) {
   case 0: // e+ and e- not written
   case 1: // e+ and e- written individually
    fN_ele = 0;
    break;
   case 2: // e+ and e- binned in redshift
    fN_ele = Nlog10z;
    for (int i=0; i<Nlog10z; i++) {
      fZ_ele[i] = 0;
      fE_ele[i] = 0.;
      fz_ele[i] = pow(10, log10z_min + (i+0.5)*dlog10z);
    }
    break;
   case 3: // only total energy written
    fN_ele = 1;
    *fZ_ele = 0;
    *fE_ele = 0.;
    *fz_ele = 0.;
    break;
  }
 }
}

void Output::Earth(const Particle& product)
{
  if (fsumm) switch (product.GetType()) {
   case eProton:
   case eNeutron:
   case eNucleus:
    fE_nuc[fN_nuc] = product.GetEint();
    fA_nuc[fN_nuc] = product.GetMassNum();
    fZ_nuc[fN_nuc] = product.GetCharge();
    fN_nuc++;
    assert(fN_nuc < nnuc_max);
    break;
   case eNeutrino_e:
   case eAntineutrino_e:
   case eNeutrino_mu:
   case eAntineutrino_mu:
    fE_neu[fN_neu] = product.GetEint();
    fF_neu[fN_neu] = product.GetFlavor();
    fN_neu++;
    assert(fN_neu < nneu_max);
    break;
   default:
    break;
  }
}

void Output::Prod(const Particle& product)
// for particles produced but not propagated (electrons and photons),
// to be propagated with external codes e.g. Elmag
{
  if (fsumm) switch (product.GetType()) {
   case ePhoton:
    fE_pho[fN_pho] = product.GetEprod();
    fz_pho[fN_pho] = product.GetZprod();
    fN_pho++;
    assert(fN_pho < npho_max);
    break;     
   case eElectron:
    if (fpairs == 1) { // individually written
      fZ_ele[fN_ele] = product.GetCharge();
      fE_ele[fN_ele] = product.GetEprod();
      fz_ele[fN_ele] = product.GetZprod();
      fN_ele++;
      assert(fN_ele < nele_max);
    } else if (fpairs == 2) { // binned in z
      int i = int((log10(product.GetZprod()) - log10z_min)/dlog10z);
      if (i < 0 || i >= Nlog10z) return;
      fZ_ele[i] += product.GetCharge();
      fE_ele[i] += product.GetEprod();
    } else if (fpairs == 3) { // only total energy written
      *fZ_ele += product.GetCharge();
      *fE_ele += product.GetEprod()/(1.+product.GetZprod());
    }
    break;
   default:
    break;
  }
}
