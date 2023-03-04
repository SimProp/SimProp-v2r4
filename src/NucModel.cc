#include <NucModel.h>
#include <Particle.h>
#include <PhysicsFunctions.h>
#include <MathFunctions.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>
#include <sstream>
#include <iostream>
#include <sigmaA_disi.h>
#include <lossesN.h>
#include <IR_fast.h>
#include <Dominguez11EBL.h>
#include <Gilmore12EBL.h>
#include <Output.h>
using namespace std;
NucModel* gModel;
extern double nFactor, aFactor;

int fLosses;
double epsIBL_min, epsIBL_max;

std::vector<int> fZ;
std::vector<int> fA;
int k_com;
double Gam_com;
double z_com;
double s_com;
//double dtdz_com;

NucModel::NucModel(eModelType model, int LossType, int Stoch, int Bdecay,
                   int PairProd) :
  fnelements(0),
  fType(model),
  fStoch(Stoch),
  fBdecay(Bdecay),
  fPairProd(PairProd)
{
  fLosses = LossType;
  switch (fLosses) {
   case 1:
    epsIBL_min = 3.3e-3; // (2.d-3) minimum contributing IR/V/UV energy (eV)
    epsIBL_max = 1.e1; // (1.6d1) maximum contributing IR/V/UV energy (eV)
    break;
   case 2:
    epsIBL_min = 2.e-3; // (3.e-3) minimum contributing IR/V/UV energy (eV)
    epsIBL_max = 1.; // (13.6)  maximum contributing IR/V/UV energy (eV)
    break;
   case 3:
    epsIBL_min = 1.26e-3;
    epsIBL_max = 11.5;
    break;
   case 4: case 5: case 6:
    epsIBL_min = pow(10,logeps_EBL[0]);
    epsIBL_max = pow(10,logeps_EBL[neps_EBL - 1]);
    break;
   case 7:
    epsIBL_min = pow(10,logeps_EBLg[0]);
    epsIBL_max = pow(10,logeps_EBLg[neps_EBLg - 1]);
    break;
  }
  if(model==eSalomon) {
    fnelements = 51;
    fe_1  =  e_1; // Photodisintegration minimum energy (MeV) 
    fe_max =  e_max; // Photodisintegration maximum energy (MeV)     
    for (int i=0; i<fnelements; i++) {
      fZ.push_back(ZZ[i]);
      fA.push_back(AA[i]);
      feth_1.push_back(eth_1[i]);
      feth_2.push_back(eth_2[i]);
      fe0_1.push_back(e0_1[i]);
      fcsi_1.push_back(csi_1[i]);
      fdelta_1.push_back(delta_1[i]);
      fe0_2.push_back(e0_2[i]);
      fcsi_2.push_back(csi_2[i]);
      fdelta_2.push_back(delta_2[i]);
      fzita.push_back(zita[i]);
    }
  } else if (model == eGauss || model == eBreitWigner) {
    std::cin >> fnelements >> fe_1 >> fe_max; // Photodisintegration maximum energy (MeV)     
    for (int i=0; i<fnelements; i++) {
      int ZZ, AA;
      double eth_1, eth_2, e0_1, csi_1, delta_1, e0_2, csi_2, delta_2, zita;
      std::cin >> ZZ >> AA >> eth_1 >> eth_2 >> e0_1 >> csi_1 >> delta_1 >> e0_2 >> csi_2 >> delta_2 >> zita;
      fZ.push_back(ZZ);
      fA.push_back(AA);
      feth_1.push_back(eth_1);
      feth_2.push_back(eth_2);
      fe0_1.push_back(e0_1);
      fcsi_1.push_back(csi_1);
      fdelta_1.push_back(delta_1);
      fe0_2.push_back(e0_2);
      fcsi_2.push_back(csi_2);
      fdelta_2.push_back(delta_2);
      fzita.push_back(zita);
    }
  } else if (model == eBreitWigner2) {
    std::cin >> fnelements >> fe_1 >> fe_max; // Photodisintegration maximum energy (MeV)  
    for (int i = 0; i < fnelements; i++) {
      int ZZ, AA;
      double nuclt, nuclh1, nuclx1, nuclw1, nuclh2, nuclx2, nuclw2, nuclc;
      double alpht, alphh1, alphx1, alphw1, alphh2, alphx2, alphw2, alphc; 
      std::cin >> ZZ >> AA 
        >> nuclt >> nuclh1 >> nuclx1 >> nuclw1 >> nuclh2 >> nuclx2 >> nuclw2 >> nuclc
        >> alpht >> alphh1 >> alphx1 >> alphw1 >> alphh2 >> alphx2 >> alphw2 >> alphc;
      fZ.push_back(ZZ);
      fA.push_back(AA);
      fNuclt.push_back(nuclt);
      fNuclh1.push_back(nuclh1);
      fNuclx1.push_back(nuclx1);
      fNuclw1.push_back(nuclw1);
      fNuclh2.push_back(nuclh2);
      fNuclx2.push_back(nuclx2);
      fNuclw2.push_back(nuclw2);
      fNuclc.push_back(nuclc);
      fAlpht.push_back(alpht);
      fAlphh1.push_back(alphh1);
      fAlphx1.push_back(alphx1);
      fAlphw1.push_back(alphw1);
      fAlphh2.push_back(alphh2);
      fAlphx2.push_back(alphx2);
      fAlphw2.push_back(alphw2);
      fAlphc.push_back(alphc);
    }   
  } else if (model == eGauss2) {
    std::cin >> fnelements >> fe_1 >> fe_max;
    for (int i = 0; i < fnelements; i++) {
      int ZZ, AA;
      double nuclt, nuclh1, nuclx1, nuclw1, nuclc;
      double alpht, alphh1, alphx1, alphw1, alphc; 
      std::cin >> AA >> ZZ 
        >> nuclt >> nuclh1 >> nuclx1 >> nuclw1 >> nuclc
        >> alpht >> alphh1 >> alphx1 >> alphw1 >> alphc;
      fZ.push_back(ZZ);
      fA.push_back(AA);
      fNuclt.push_back(nuclt);
      fNuclh1.push_back(nuclh1);
      fNuclx1.push_back(nuclx1);
      fNuclw1.push_back(nuclw1);
      fNuclc.push_back(nuclc);
      fAlpht.push_back(alpht);
      fAlphh1.push_back(alphh1);
      fAlphx1.push_back(alphx1);
      fAlphw1.push_back(alphw1);
      fAlphc.push_back(alphc);
    }   
  }
}

NucModel::~NucModel() {}

//void NucModel::Validate(Particle* input) {}


int NucModel::TestPrimary(int mass) {

  if(fType==eOnlyProton && mass==1) { //each particle is transformed into proton
    return 1;
  }

  int found = binary(mass); 
  return found;

}


int NucModel::binary(int mass) { //return the charge of the element ... if found
    
  int first = 0;
  int last = fnelements;

  while(first<=last) {
    int mid = (first + last)/2;
    if(mass > fA[mid])
      last = mid - 1;
    else if (mass < fA[mid])
      first = mid + 1;
    else
      return fZ[mid];
  }

  return -1;
}


Particle NucModel::GetRandomParticle(double E, double z) {
  
  const int index = gRandom->Integer(fnelements);
  const int A = fA[index];
  const int Z = fZ[index];

  return Particle(A, Z, 0, E, z, (A == 1) ? eProton : eNucleus);
}


vector<Particle> NucModel::PropagateParticle(Particle* input)
{
  cout<<"... propagating " << input->GetName() << " (" << input->GetMassNum() 
      << ", " << input->GetCharge() << ") " << " Ecurr " << input->GetEprod()
      << " from z " << input->GetZprod() << endl;

  vector<Particle> output;
  switch (input->GetType()) {
   case eProton:
    if (fStoch <= 0)
      DetermProton(input);
    else
      output = PropagateNucleus(input);
    break;      
   case eNeutron:
   case eNucleus:
    output = PropagateNucleus(input);
    break;      
   case ePion:
    output = PropagatePion(input);
    break;
   case ePhoton:
    output = PropagatePhoton(input);
    break;
   case eMuon:
    output = PropagateMuon(input);
    break;
   case eNeutrino_e:
   case eAntineutrino_e:
   case eNeutrino_mu:
   case eAntineutrino_mu:
    output = PropagateNeutrino(input);
    break;
   case eElectron:
    output = PropagateElectron(input);
   case eNull:
    break;
  }

  // Now, input->GetIntMult is: 
  // <= 0   if the particle reaches Earth
  // n      if the particle interacts (n = number of particles after the interaction)
  // 1000   for electrons and photons (not propagated by SimProp; use external tools e.g. Elmag)
  // 1000+n if the particle decays (n = number of decay products)  

  if (input->GetIntMult() > 1000) {
    cout << "........ " << input->GetName() << " decays into " << input->GetIntMult()-1000
         << " particles at z=" << input->GetZint() << endl;
  } else if (input->GetIntMult() == 1000) {
    cout << "........ " << input->GetName() << " produced with E=" << input->GetEprod()
         << " at z=" << input->GetZprod() << endl;
  } else if (input->GetIntMult()>0) {
    assert(output.size()>0);
    cout << "...... int. type " << input->GetIntMult() << " (" << input->GetMassNum()
         << ", " << input->GetCharge() << ") -> (" << output[0].GetMassNum()
         << ", " << output[0].GetCharge() << ") zfin=" << input->GetZint()
         << " Efin=" << input->GetEint() << endl; 
  } else if (input->GetZint() <= 0.) {
      input->SetZint(0);
      cout << "........... " << input->GetName() << " (" << input->GetMassNum() <<", "
           << input->GetCharge() << ") reaches Earth with E = " << input->GetEint() 
	   << endl;
  } 
  return output;
}


void NucModel::DetermProton(Particle* input)
{
  double Gamdet = evolveN(input->GetZprod(),input->GetGprod(),0.); 
  input->SetZint(0.);
  input->SetGint(Gamdet);
  input->SetIntMult(0);
}

int NucModel::GetBetaDecayStableZ(int A)
{
  int Zstable[] = {0, 1,1,2,2,2, 3,3,4,4,5, 5,6,6,7,7, 8,8,8,9,10, 
  10,10,11,12,12, 12,13,14,14,14, 15,16,16,16,17, 16,17,18,19,18,
  19,20,20,20,21, 22,22,22,22,22, 23,24,24,24,25, 26,26,26,27,28};
  int Z = binary(A);
  if (Z > 0)
    return Z;
  else
    return Zstable[A];
}

vector<Particle> NucModel::PropagateNucleus(Particle * input)
{
    const int Acurr = input->GetMassNum();
    const int Zcurr = input->GetCharge();
    if (Zcurr != GetBetaDecayStableZ(Acurr)) {
      if (fBdecay)
        return BetaDecay(input);
      else {
        input->SetCharge(GetBetaDecayStableZ(Acurr));
        if (input->GetType() == eNeutron) input->SetType(eProton);      
      }
    }
    vector<Particle> output;
    if (input->GetType() == eProton && fStoch <= 0) {
      DetermProton(input);
      return output;
    }
    const int nBr = input->GetBranch() + 1; 
    eProcType proc;
    vector<Particle> electrons;
    if (Acurr == 5 || Acurr == 8) {
      input->SetZint(input->GetZprod());
      input->SetGint(input->GetGprod());
      input->SetIntMult(1002);
      if (Acurr == 8)
        proc = eAlphaProd;
      if (Acurr == 5)
        proc = eNucleonProd;
    } else
      proc = GetProcess(input, fPairProd ? &electrons : 0);
    const double zfin = input->GetZint(); 
    const double Gfin = input->GetGint();
    int Afin, Zfin, intType, protons, neutrons;
    switch (proc) {
     case eEarth:
      return electrons;
     case eAlphaProd:
      Afin = Acurr - 4;
      Zfin = Zcurr - 2;
      intType = 4;
      protons = 0;
      neutrons = 0;
      break;
     case eNucleonProd:
     case eHadron:
     case eHadron_IR:
      Afin = Acurr - 1;
      if (Acurr == 3 || (Acurr != 5 && gRandom->Uniform(0.,Acurr) < Zcurr)) {
        protons = 1;
        neutrons = 0;
        Zfin = Zcurr - 1;      
      }
      else {
        protons = 0;
        neutrons = 1;
        Zfin = Zcurr;
      }
      intType = 1;
      break;
     case eDisi:
      intType = input->GetIntMult();
      Afin = Acurr - intType;
      if (Afin > 1 && GetJNuc(Afin) < 0) { // e.g. Afin == 8, Afin == 7...
        int jnuc = -GetJNuc(Afin); 
        Afin = fA[jnuc];            // i.e. Afin == 4
        intType = Acurr - Afin;
      }
      if (Acurr == 9 && Afin == 4) {
           proc = eAlphaProd;
           intType = 1;
           Zfin = 2;
           protons = 0;
           neutrons = 1;
      } else {
          Zfin = Zcurr;
          for (int i = 0; i < intType; i++)
              if (gRandom->Uniform(0.,Acurr-i) < Zfin)
                  Zfin--; //loses a proton; else: loses a neutron
          if (Afin > 1) {
            if (Zfin == 0) Zfin++;
            if (Zfin == Afin) Zfin--; // avoid stuff like 2He ...
          }
          protons = Zcurr - Zfin;
          neutrons = intType - protons;
      }
      break;
    }
    double Epion = 0;
    int pioncharge = 0;
    if (proc == eHadron || proc == eHadron_IR) {
        Epion = Photopion(zfin,input->GetGint(),proc);
        if (gRandom->Uniform(0.,1.) < 1./3.) { // assuming isospin invariance
            if (protons) {
                protons--;
                neutrons++;
                pioncharge = +1;
                cout << "....... stacking pion+, z=" << zfin << ", E=" << Epion;
            } else {
                protons++;
                neutrons--;
                pioncharge = -1;
                cout << "....... stacking pion-, z=" << zfin << ", E=" << Epion;
            }           
        } else
            cout << "....... stacking pion0, z=" << zfin << ", E=" << Epion;
        if (proc == eHadron)    cout << "\n";
        if (proc == eHadron_IR) cout << " (from IBL)\n";
    }
    if (Afin == 1) {
      if (Zfin == 1)
        protons++;
      else
        neutrons++;
      Afin = 0;
    }
    if (!fBdecay) {
        if (Afin > 1) Zfin = GetBetaDecayStableZ(Afin);
        protons = intType;
        neutrons = 0;
    }
    if (Afin != 0)
        output.push_back(Particle(Afin,Zfin,nBr,Afin*mN*Gfin,zfin,eNucleus));
    else
        if (Acurr > 1) cout << "Nucleus fully fragmented\n";
    if (proc == eAlphaProd) {
        cout << "....... stacking 4He, z=" << zfin << ", E=" << 4*mN*Gfin << '\n';    
        output.push_back(Particle(4,2,nBr,4*mN*Gfin,zfin,eNucleus));
    }
    if (protons || neutrons)
        cout << "....... stacking " << protons << " p + " << neutrons 
             << " n, z=" << zfin << ", E=" << mN*Gfin-Epion << '\n';
    for (int i = 0; i < protons; i++)
        output.push_back(Particle(1,1,nBr,mN*Gfin-Epion,zfin,eProton));
    for (int i = 0; i < neutrons; i++)
        output.push_back(Particle(1,0,nBr,mN*Gfin-Epion,zfin,eNeutron));
    if (proc == eHadron || proc == eHadron_IR)
        output.push_back(Particle(0,pioncharge,nBr,Epion,zfin,ePion));
    output.insert(output.end(), electrons.begin(), electrons.end());
    input->SetIntMult(intType);
    return output; 
}

// -------------------------------------------------------------------------

int NucModel::GetJNuc(int ANuc)
{
  int jNuc;
  int i;
  if(ANuc > fA[0] || ANuc < 0) return 0;
  if(ANuc == 1) return -1;
  for(i = 0; i < fnelements; i++) if(fA[i] <= ANuc) break;
  if(fA[i] == ANuc) jNuc = i;
  else jNuc = -i;
  return jNuc;
}
// -------------------------------------------------------------------------

eProcType NucModel::GetProcess(Particle *input, vector<Particle>* electrons=0)
{
  const int A =  input->GetMassNum();
  const int Z =  input->GetCharge();
  const int J = GetJNuc(A);
  double Gam  = input->GetGprod();
  double z = input->GetZprod();
// double B = input->GetB();
// double Lcoh = input->GetLcoh();
// double B = fB[0];
// double Lcoh = fB[1];

  const int nstep = 200*z+20;
  const double zOri = z;
  const double zEnd = 0;
  double rate=TotalRate(A, Z, J, Gam, z, NULL)/(H0*(1+z)*sqrt(om_m*pow(1+z,3)+om_l)); // int. rate
  double logpSurv=0;
//  double thetasumq=0.,deltatheta=0.;
//  double Ttot=0;
//  const double rad = 180./3.14159;
//  const double pi = 3.14159;
//  double delta_fin=0;

  const double logu = log(gRandom->Rndm());
  double z_old, Gam_old, rate_old, logpSurv_old;
  eProcType proc = eEarth;
  for (int j = 1; j <= nstep; j++) {
    z_old = z;
    Gam_old = Gam;
    rate_old = rate;
    logpSurv_old = logpSurv;
    double t = (double)j/nstep;
    z = (1-t)*zOri + t*zEnd;
    if (j == nstep) z = 0;  // just in case there are rounding errors
//    double theta=GetTheta(z,z_old,A,Z,G_old,B,Lcoh);
//    double th = gRandom->Gaus(0., sqrt(theta));
//    double th = gRandom->Uniform(0., theta);  //    ****** test ****
//    thetasumq += theta*rad*theta*rad;
//    double dtdz= Getdtdzeff(z,z_old,A,Z,G_old,B,Lcoh,theta);
//    dtdz_com=dtdz;
//    Ttot-=((dtdz*(1+z)*(z-z_old)));
    Gam = evolveA0(A,z_old,Gam_old,z);
    rate = TotalRate(A, Z, J, Gam, z, NULL)/(H0 *(1+z) * sqrt(om_m *pow(1+z,3) + om_l));
    logpSurv += (z-z_old)*(rate+rate_old)/2.;	    
    if (logpSurv <= logu) {
      z -= (logu-logpSurv)/(logpSurv_old-logpSurv)*(z-z_old);
      Gam = evolveA0(A,z_old,Gam_old,z);
      (void)TotalRate(A, Z, J, Gam, z, &proc);
    }
    if (electrons) {
      // We are assuming that the fractional energy loss in each interaction
      // is that at threshold; this is inaccurate at large nucleus energies,
      // but the spectrum of low-energy gamma rays at Earth doesn't strongly
      // depend on the spectrum of electrons/positrons but mostly only on their
      // total energy, so it doesn't matter very much
      // (AdM, 18 Feb 2016)
      const double K = 2*me/(A*mN+2*me);
      double npairs_avg = log(Gam_old/(1+z_old)/Gam*(1+z))/K;
      int npairs = gRandom->Poisson(npairs_avg);
      const int nBr = input->GetBranch() + 1;
      for (int i = 0; i < npairs; i++) {
        double Ee = K*A*exp(gRandom->Uniform(log(Gam),log(Gam_old)))*mN/2;
        double zprod = gRandom->Uniform(z,z_old);
        electrons->push_back(Particle(0, -1, nBr, Ee, zprod, eElectron));
        electrons->push_back(Particle(0, +1, nBr, Ee, zprod, eElectron));
      }
    }
    if (proc != eEarth)
      break;
  }
//  delta_fin=sqrt(thetasumq);
//  distInt = qgaus(rr,z_int,z_ori); 
//  input->SetTtot(Ttot);
//  input->Setdelta_fin(delta_fin);
//  const double dfin = GetR(0., z_int);
  int ndis = 0;
  switch (proc) {
   case eEarth:
    assert (z == 0);
    ndis = 0; break;
   case eHadron:
   case eHadron_IR:
   case eNucleonProd:
    assert (z > 0);
    ndis = 1; break;
   case eAlphaProd:
    assert (z > 0);
    ndis = 4; break;
   case eDisi:
    assert (z > 0);
    double bdis = SecInY*SpeedOfLight*m2cm/(Mpc2cm*losseAdisi(A,z,Gam,J,ndis)); 
    break;
//   case eDecay:
//    ndis = 1003; break;
  }
  input->SetZint(z);
  input->SetGint(Gam);
  input->SetIntMult(ndis);
  return proc; 
}

// -------------------------------------------------------------------------


double NucModel::evolveA0 (int A_0, double z_in, double Gam_in,
		 double z_out)
{
//.............................................................................
// Evolution with red-shift of Gamma and A for nuclei,
// .............................................................................
  k_com = GetJNuc(A_0);
  //cout<<k_com<<"  "<<A_0<<"  "<<Gam_in<<endl;

  const double eps=1.e-9; //adm
  double h1 =1.e-4;  //adm
  double hmin = 0.;

  int nok, nbad;

  double Lny0 = log(Gam_in);

  double Lny = odeint(Lny0,  z_in, z_out, eps, h1, hmin, nok , nbad, derivsA0);

  double Gam_out = exp(Lny);
  return Gam_out;
}  



// -------------------------------------------------------------------------


double NucModel::losseAdisi(int Acurr, double z, double Gam, int k0, int &multi)
{
  //............................................................................. 
  // Subroutine that computes the photo-disintegration energy losses for nuclei 
  // expressed in 1/y. All backgrounds (CMB+IR/V/UV) are taken into account
  //.............................................................................
  if (Acurr == 1)
    return 0.;

  double beta_disi=0.;

  Gam_com = Gam;
  z_com = z;
  k_com = k0;
  s_com = 0;
  gModel = this;

  double sum1=0;
  double sum2=0;
  double sum3=0;
  
  //.... photodisintegration energy losses (CMB + IRV/UV) ....
  
  double lxmin1 = log(eth_1[k0]);
  double lxmax1 = log(eth_2[k0]);
  sum1 = qgaus (fdisi, lxmin1, lxmax1);
  
  double lxmin2 = log(eth_2[k0]);
  double lxmax2 = log(e_1);
  sum2 = qgaus (fdisi, lxmin2, lxmax2);
  
  double lxmin3 = log(e_1);
  double lxmax3 = log(e_max);
  sum3 = qgaus (fdisi, lxmin3, lxmax3);

  double sum=sum1+sum2+sum3;
 
  // sp
  if (sum<=0) sum = 0;
  if (multi < 0)
      return sum;
  
  multi = 0; 
  // sp
  
  double urange = gRandom->Rndm();
  double umult  = gRandom->Rndm();

//  int erange = 1;
  if (urange <= sum1/sum) 
    multi = 1;
  else if (urange>(sum1/sum) && urange<=(sum1+sum2)/sum) {
//    erange = 2;
    s_com = 1;
    double sum21 = qgaus (fdisi, lxmin2, lxmax2);
    multi = 2;
    if (umult < sum21/sum2) multi = 1;
  } else {
//    erange = 3;
    int dim = (int)BRdisi(Acurr, 0);
    double cum = 0;
    for (int n = 1; n <= dim; n++ ) {
      cum += BRdisi(Acurr, n);
      if (umult < cum) {
        multi = n;
        break;
      }
    }
    if (Acurr == 9) multi = 5;
  }
  
  beta_disi = sum; 
  return beta_disi;
}

// -------------------------------------------------------------------------

double derivsA0 (double z, double Lny)
{
// Right hand side of the differential equations for energy losses
 extern int Stoch;
 double Gam = exp(Lny);      
 double dtdz = 1./ (H0 * (1+z) * sqrt(om_m * pow((1+z),3) + om_l)); 
 // double dtdz=dtdz_com;
 double beta;
 if (Stoch) beta = losseA0pair(z,Gam);
 else beta = losseA0(z,Gam);
 if(k_com >= 0 && fA[k_com]>1 && fA[k_com]<56) beta *= (fA[k_com]/4.);
 // beta_pair *= fZ[k_com]*fZ[k_com]/fA[k_com];

 return dtdz * beta + 1/(1+z);
}

// -------------------------------------------------------------------------

double integIBL(double eps, double z)
// integral from e to +inf of n_IBL(e)/(2e^2) de
{
  double In_IR = 0.;
  if (fLosses == 3) eps /= (1+z);
  if (eps < epsIBL_max) {
    if (z > 5.)
      z = 5.;
    if (eps < epsIBL_min)
      eps = epsIBL_min;
    if (fLosses == 1) { // assuming the spectrum of Malkan & Stecker (2006)
      static double *y2a = NULL;
      if (!y2a) {
        y2a = new  double[nE_IR*nz_IR];
        splie2(z_IR,lE_IR,lint_IR,nz_IR,nE_IR,y2a);
      }
      double logeps = log10(eps);
      double logIn_IR = -splin2(z_IR,lE_IR,lint_IR,y2a,nz_IR,nE_IR,z,logeps);
      In_IR = pow(10.,logIn_IR);
    } else if (fLosses == 2) { // using a power-law approximation
      static const double gam_IR = 2.5; // IR/V/UV background power law index 
      static const double N0_IR = 2.*5.207e-2/(gam_IR+1); // (min) normalization 1/y/MeV^2/mbarn
      double sum = pow(eps, -(gam_IR+1)) - pow(epsIBL_max, -(gam_IR+1));
      if (z <= 1.4) In_IR = N0_IR * pow(1+z, 3.1) * sum;
      else In_IR = N0_IR * pow(1+1.4, 3.1) * sum;
    } else if (fLosses == 3) { // Kneise (2004), with scaling in z as in CRPropa (?)
      static const int nz = 8;
      static const double z_[nz] = { 0.0, 0.2, 0.4, 0.6, 1.0, 2.0, 3.0, 4.0, };
      static const double scale_[nz] = {
        1., 1.69373789, 2.58847457, 3.61780921, 5.19804508, 7.38709744,
        8.54708679,  7.86045266, 
      };
      static const int nE = 64;
      static const double logE_[nE] = {
        -2.9375, -2.8750, -2.8125, -2.7500, -2.6875, -2.6250, -2.5625, -2.5000, 
        -2.4375, -2.3750, -2.3125, -2.2500, -2.1875, -2.1250, -2.0625, -2.0000, 
        -1.9375, -1.8750, -1.8125, -1.7500, -1.6875, -1.6250, -1.5625, -1.5000, 
        -1.4375, -1.3750, -1.3125, -1.2500, -1.1875, -1.1250, -1.0625, -1.0000, 
        -0.9375, -0.8750, -0.8125, -0.7500, -0.6875, -0.6250, -0.5625, -0.5000, 
        -0.4375, -0.3750, -0.3125, -0.2500, -0.1875, -0.1250, -0.0625,  0.0000, 
         0.0625,  0.1250,  0.1875,  0.2500,  0.3125,  0.3750,  0.4375,  0.5000, 
         0.5625,  0.6250,  0.6875,  0.7500,  0.8125,  0.8750,  0.9375,  1.0000, 
      };
      const double logI_[nE] = {
         7.6425,  7.6275,  7.5860,  7.5319,  7.4591,  7.3881,  7.3012,  7.1957, 
         7.0854,  6.9568,  6.8105,  6.6459,  6.4594,  6.2580,  6.0260,  5.7639, 
         5.4873,  5.1737,  4.8592,  4.5428,  4.2482,  4.0003,  3.7928,  3.6128, 
         3.4472,  3.2884,  3.1248,  2.9538,  2.7780,  2.5914,  2.4035,  2.2252, 
         2.0528,  1.9009,  1.7787,  1.6628,  1.5511,  1.4426,  1.3305,  1.2056, 
         1.0659,  0.9215,  0.7609,  0.5834,  0.3927,  0.1881, -0.0266, -0.2510, 
        -0.4875, -0.7322, -0.9862, -1.2402, -1.4892, -1.7349, -1.9741, -2.2072, 
        -2.4470, -2.7033, -2.9959, -3.3409, -3.7292, -4.1682, -4.6959, -5.3645, 
      };
      static double *scale2a, *logI2a;
      if (!scale2a) {
        scale2a = new double[nz];
        spline(z_,scale_,nz,1.e30,1.e30,scale2a);
      }
      if (!logI2a) {
        logI2a = new double[nE];
        spline(logE_,logI_,nE,1.e30,1.e30,logI2a);
      }
      double scale = splint(z_,scale_,scale2a,nz,z)/pow(1+z,3);
      double I = pow(10,splint(logE_,logI_,logI2a,nE,log10(eps)));
      In_IR = (1+z)*I*scale;
    } else if (fLosses == 4 || fLosses == 5 || fLosses == 6) { //Dominguez (2011)
      static const double *logI_EBL = 
      	fLosses == 4 ? logImid_EBL :
      	fLosses == 5 ? logIlow_EBL :
      	logIupp_EBL;
      static double *y2a = NULL;
      if (!y2a) {
        y2a = new  double[neps_EBL*nz_EBL];
        splie2(z_EBL,logeps_EBL,logI_EBL,nz_EBL,neps_EBL,y2a);
      }
      double logeps = log10(eps);
      double logIn_IR = -splin2(z_EBL,logeps_EBL,logI_EBL,y2a,nz_EBL,neps_EBL,z,logeps);
      In_IR = pow(10.,logIn_IR) * pow(1+z, 3);
    } else if (fLosses == 7) { // Gilmore (2012)
      static double *y2a = NULL;
      if (!y2a) {
        y2a = new  double[neps_EBLg*nz_EBLg];
        splie2(z_EBLg,logeps_EBLg,logI_EBLg,nz_EBLg,neps_EBLg,y2a);
      }
      double logeps = log10(eps);
      double logIn_IR = -splin2(z_EBLg,logeps_EBLg,logI_EBLg,y2a,nz_EBLg,neps_EBLg,z,logeps);
      In_IR = pow(10.,logIn_IR) * pow(1+z, 3);
    }
  }
  return In_IR;
}


double inv_integIBL(double In_IR, double z)
// inverse function of the above
// Armando di Matteo, 30 Apr 2013
{
  double eps=0.;
  if (z > 5.)
    z = 5.;
  if (fLosses == 1) {
    static double *y2a = NULL;
    if (!y2a) {
      y2a = new double[nE_IR*nz_IR];
      splie2_inv(z_IR,lE_IR,lint_IR,nz_IR,nE_IR,y2a);
    }
    if (In_IR <= 0.) return epsIBL_max;
    double logIn_IR = log10(In_IR);
    double logeps = splin2_inv(z_IR,lE_IR,lint_IR,y2a,nz_IR,nE_IR,z,-logIn_IR);
    eps = pow(10.,logeps);
  } else if (fLosses == 2) {
    static const double gam_IR = 2.5; // IR/V/UV background power law index 
    static const double N0_IR = 2.*5.207e-2/(gam_IR+1); // (min) normalization 1/y/MeV^2/mbarn
    double sum;
    if (z <= 1.4) sum = In_IR/N0_IR*pow(1+z,-3.1) + pow(epsIBL_max, -(gam_IR+1));
    else sum = In_IR/N0_IR*pow(1+1.4,-3.1) + pow(epsIBL_max, -(gam_IR+1));
    eps = pow(sum, -1./(gam_IR+1));
  } else if (fLosses == 3) {
    static const int nz = 8;
    static const double z_[nz] = { 0.0, 0.2, 0.4, 0.6, 1.0, 2.0, 3.0, 4.0, };
    static const double scale_[nz] = {
      1., 1.69373789, 2.58847457, 3.61780921, 5.19804508, 7.38709744,
      8.54708679,  7.86045266, 
    };
    static const int nE = 64;
    const double mlogI_[nE] = {
      -7.6425, -7.6275, -7.5860, -7.5319, -7.4591, -7.3881, -7.3012, -7.1957, 
      -7.0854, -6.9568, -6.8105, -6.6459, -6.4594, -6.2580, -6.0260, -5.7639, 
      -5.4873, -5.1737, -4.8592, -4.5428, -4.2482, -4.0003, -3.7928, -3.6128, 
      -3.4472, -3.2884, -3.1248, -2.9538, -2.7780, -2.5914, -2.4035, -2.2252, 
      -2.0528, -1.9009, -1.7787, -1.6628, -1.5511, -1.4426, -1.3305, -1.2056, 
      -1.0659, -0.9215, -0.7609, -0.5834, -0.3927, -0.1881,  0.0266,  0.2510, 
       0.4875,  0.7322,  0.9862,  1.2402,  1.4892,  1.7349,  1.9741,  2.2072, 
       2.4470,  2.7033,  2.9959,  3.3409,  3.7292,  4.1682,  4.6959,  5.3645, 
    };
    static const double logE_[nE] = {
      -2.9375, -2.8750, -2.8125, -2.7500, -2.6875, -2.6250, -2.5625, -2.5000, 
      -2.4375, -2.3750, -2.3125, -2.2500, -2.1875, -2.1250, -2.0625, -2.0000, 
      -1.9375, -1.8750, -1.8125, -1.7500, -1.6875, -1.6250, -1.5625, -1.5000, 
      -1.4375, -1.3750, -1.3125, -1.2500, -1.1875, -1.1250, -1.0625, -1.0000, 
      -0.9375, -0.8750, -0.8125, -0.7500, -0.6875, -0.6250, -0.5625, -0.5000, 
      -0.4375, -0.3750, -0.3125, -0.2500, -0.1875, -0.1250, -0.0625,  0.0000, 
       0.0625,  0.1250,  0.1875,  0.2500,  0.3125,  0.3750,  0.4375,  0.5000, 
       0.5625,  0.6250,  0.6875,  0.7500,  0.8125,  0.8750,  0.9375,  1.0000, 
    };
    static double *scale2a, *logE2a;
    if (!scale2a) {
      scale2a = new double[nz];
      spline(z_,scale_,nz,1.e30,1.e30,scale2a);
    }
    if (!logE2a) {
      logE2a = new double[nE];
      spline(mlogI_,logE_,nE,1.e30,1.e30,logE2a);
    }
    double scale = splint(z_,scale_,scale2a,nz,z)/pow(1+z,3);
    double I = In_IR/(1+z)/scale;
    if (I <= 0) return (1+z)*epsIBL_max;
    eps = (1+z)*pow(10,splint(mlogI_,logE_,logE2a,nE,-log10(I)));
  } else if (fLosses == 4 || fLosses == 5 || fLosses == 6) { //Dominguez (2011)
    static double *y2a = NULL;
    static const double *logI_EBL = 
    	fLosses == 4 ? logImid_EBL :
    	fLosses == 5 ? logIlow_EBL :
    	logIupp_EBL;
    if (!y2a) {
      y2a = new double[neps_EBL*nz_EBL];
      splie2_inv(z_EBL,logeps_EBL,logI_EBL,nz_EBL,neps_EBL,y2a);
    }
    if (In_IR <= 0.) return epsIBL_max;
    double logIn_IR = log10(In_IR / pow(1+z, 3));
    double logeps = splin2_inv(z_EBL,logeps_EBL,logI_EBL,y2a,nz_EBL,neps_EBL,z,-logIn_IR);
    eps = pow(10.,logeps);
  } else if (fLosses == 7) { // Gilmore (2012)
    static double *y2a = NULL;
    if (!y2a) {
      y2a = new double[neps_EBLg*nz_EBLg];
      splie2_inv(z_EBLg,logeps_EBLg,logI_EBLg,nz_EBLg,neps_EBLg,y2a);
    }
    if (In_IR <= 0.) return epsIBL_max;
    double logIn_IR = log10(In_IR / pow(1+z, 3));
    double logeps = splin2_inv(z_EBLg,logeps_EBLg,logI_EBLg,y2a,nz_EBLg,neps_EBLg,z,-logIn_IR);
    eps = pow(10.,logeps);
  }
  //assert(epsIBL_min <= eps && eps <= epsIBL_max);
  //if (eps < epsIBL_min) eps = epsIBL_min;
  //if (eps > epsIBL_max) eps = epsIBL_max; 
  return eps;
}


double fdisi(double lx)
// lx = ln(eps'/MeV)
{
  double epsprime = exp(lx); // in MeV
  double sigma = 0.; // in mb
  double sigma_1 = 0., sigma_2 = 0., sigma_3 = 0.;
  gModel->sigma_disi (k_com, epsprime, sigma_1, sigma_2, sigma_3);  
  if(s_com == 0) sigma = sigma_1 + sigma_2 + sigma_3;
  if(s_com == 1) sigma = sigma_1;
  if(s_com == 2) sigma = sigma_2;
  if(s_com == 3) sigma = sigma_3;
  double eps = epsprime / (2*Gam_com) * 1.e6; // in eV
  double integ = integCMB(eps, z_com) + integIBL(eps, z_com); // in 1/(y Mev^2 mb)
  return epsprime*epsprime/(Gam_com*Gam_com)*sigma*integ; //in y^-1
}

// double fpion(double lx)
// {
// // integral function of the pion production energy losses

//   double Ene = exp(lx);

//   double sigma = 0;

//   sigma_pion (Ene, sigma);

//   double N0_CMB = 1.452e12;  // normalization CMB 1/y/MeV^2/mbarn
//   double kT     = 2.327e-10; // KT (MeV)
//   double xx1    = Ene / (2.*(1.+z_com) * Gam_com * kT);
//   double In_CMB = pow(1+z_com,3) * N0_CMB * log(1./(1-exp(-xx1)));

//   double ln = In_CMB / ((1.+z_com) * (1.+z_com));
//   double retvalue = Ene * sigma * Ene * ln / (1. * Gam_com * Gam_com); // 1/y

//   return retvalue;
// }

double fpion_IR(double lx)
// lx = ln(eps'/MeV)
{
  double epsprime = exp(lx); // in MeV
  double sigma = 0.; // in mb
  sigma_pion(epsprime, sigma);  
  double eps = epsprime / (2*Gam_com) * 1.e6; // in eV
  double integ = integIBL(eps, z_com);
  return epsprime*epsprime/(Gam_com*Gam_com)*sigma*integ; //in y^-1
}

// -------------------------------------------------------------------------
double NucModel::BRdisi(int A, int imult)
{
  const int Alow[4] = {1, 9, 10, 23};
  const int Aup[4]  = {4, 9, 22, 56};
  const int dim[4]  = {2, 1,  6, 15};

  const double fBR[] = {0.8, 0.2,                          // He
			1.0,                               // Be
			0.1, 0.3, 0.1, 0.1, 0.2, 0.2,      // B,C,N,O,F,Ne
			0.1, 0.35, 0.1, 0.05, 0.15,        // Na,...,Fe
			0.045, 0.040, 0.035, 0.030, 0.025,
			0.020, 0.018, 0.015, 0.012, 0.010};
  int indx = -1;
  int ptr  = 0;
  for (int i = 0; i < 4; i++) {
    if (A>=Alow[i] && A<=Aup[i]) {
      indx = i;
      break;
    }
    ptr += dim[i];
  }

  int mult = 0;
  if (indx >= 0) mult = dim[indx];

  if (imult <= 0) return (double)mult;
  else return fBR[ptr+imult-1];
}

// -------------------------------------------------------------------------


void NucModel::sigma_disi (int k, double Ene, 
		 double &sigma_1, double &sigma_2, double &sigma_3)
{
// ............................................................................. 
//  Subroutine that computes the photodisintegration cross-section,
//  requires (Z,A) and photon energy (MeV) in input. The output is in mbarn 
//  and is divided into sigma_1: single nucleon emission (g,n) and (g,p); 
//  sigma_2: double nucleon emission (g,2n), (g,2p) and (g,np); 
//  sigma_3: constant normalized cross-section. We have followed the recipes 
//  of Stecker and Salamon ApJ 512 (1999) 521.
// .............................................................................
  sigma_1 = 0.;
  sigma_2 = 0.;
  sigma_3 = 0.;
  if (Ene >= e_max) return;
  assert(k >= 0);
  int A  = fA[k];
  int Ze = fZ[k];
  //...... Normalization factor (mb*MeV) ..................
  
  if (A <= 4) {
    assert(A == 2 || A == 3 || A == 4); k = nA + 1 - A;  
    //...... Normalization factor (mb*MeV) ..................
    double Sigma_d = 60.* (double)((A-Ze)*Ze)/ (double)(A); // mb*MeV
    
    if ((Ene >= eth_1[k]) && (Ene <= e_1)) {
      double W1 = W(eth_1[k], e0_1[k], delta_1[k]);
      sigma_1 = (csi_1[k] * Sigma_d/W1) *
        exp(-2.*((Ene-e0_1[k])/delta_1[k])*((Ene-e0_1[k])/delta_1[k]));
    } else sigma_1 = 0.;

    if ((Ene >= eth_2[k]) && (Ene <= e_1)) {
      double W2 = W(eth_2[k], e0_2[k], delta_2[k]);
      sigma_2 = (csi_2[k] * Sigma_d/W2 )*
        exp(-2.*((Ene-e0_2[k])/delta_2[k])*((Ene-e0_2[k])/delta_2[k]));
    } else sigma_2=0.;

    if ((Ene >= e_1) && (Ene <=e_max))
      {
      sigma_3 = zita[k] * Sigma_d/(e_max-e_1);
      }
    else sigma_3=0.;
    
    if (fType == eBreitWigner2 || fType == eGauss2) {
      sigma_1 += 2 * sigma_2 + 1.2 * sigma_3;
      sigma_2 = 0; sigma_3 = 0;
    }
    return;
  }
  
  if (fType == eBreitWigner2) {
    sigma_2 = 0;
    if (Ene <= fe_1) {
      sigma_1 = fNuclh1[k]/(1.+pow((Ene-fNuclx1[k])/fNuclw1[k],2))
              + fNuclh2[k]/(1.+pow((Ene-fNuclx2[k])/fNuclw2[k],2));
      sigma_3 = fAlphh1[k]/(1.+pow((Ene-fAlphx1[k])/fAlphw1[k],2))
              + fAlphh2[k]/(1.+pow((Ene-fAlphx2[k])/fAlphw2[k],2));
    } else if (Ene <= fe_max) {
      sigma_1 = fNuclc[k];
      sigma_3 = fAlphc[k];
    } else {
      sigma_1 = 0;
      sigma_3 = 0;
    }
    return;
  }

  if (fType == eGauss2) {
    sigma_2 = 0;
    if (Ene <= fe_1) {
      sigma_1 = fNuclh1[k]*exp(-pow(Ene-fNuclx1[k],2)/fNuclw1[k]);
      sigma_3 = fAlphh1[k]*exp(-pow(Ene-fAlphx1[k],2)/fAlphw1[k]);
    } else if (Ene <= fe_max) {
      sigma_1 = fNuclc[k];
      sigma_3 = fAlphc[k];
    } else {
      sigma_1 = 0;
      sigma_3 = 0;
    }
    return;
  }

  double Sigma_d = 60.* (double)((A-Ze)*Ze)/ (double)(A); // mb*MeV
  
  if ((Ene >= feth_1[k]) && (Ene <= fe_1)) {
    if (fType == eBreitWigner && A > 4) {
      sigma_1 = fcsi_1[k] / 
        (1. + ((Ene-fe0_1[k])/fdelta_1[k])*((Ene-fe0_1[k])/fdelta_1[k]));
    } else {
      double W1 = W(feth_1[k], fe0_1[k], fdelta_1[k]);
      sigma_1 = (fcsi_1[k] * Sigma_d/W1) *
        exp(-2.*((Ene-fe0_1[k])/fdelta_1[k])*((Ene-fe0_1[k])/fdelta_1[k]));
    }
  } else sigma_1 = 0.;

  if ((Ene >= feth_2[k]) && (Ene <= fe_1)) {
    if (fType == eBreitWigner && A > 4) {
      sigma_2 = fcsi_2[k] / 
        (1. + ((Ene-fe0_2[k])/fdelta_2[k])*((Ene-fe0_2[k])/fdelta_2[k]));
    } else {
      double W2 = W(feth_2[k], fe0_2[k], fdelta_2[k]);
      sigma_2 = (fcsi_2[k] * Sigma_d/W2 )*
        exp(-2.*((Ene-fe0_2[k])/fdelta_2[k])*((Ene-fe0_2[k])/fdelta_2[k]));
    }
  } else sigma_2=0.;

  if ((Ene >= fe_1) && (Ene <= fe_max))
    {
    sigma_3 = fzita[k] * Sigma_d/(fe_max-fe_1);
    }
  else sigma_3=0.;
  
}
// -------------------------------------------------------------------------

void sigma_pion(double Ene, double& sigma)
// Ene = photon energy in NRF (MeV)
// sigma in mb
{
  if (Ene >= epsthr) {
    double s = (mN/1e9)*(mN/1e9)  + 2*(mN/1e9)*(Ene/1e3); // in GeV^2
    if (s <= 5.) {
      static const int n = 38;
      static const double s_[n] = {
        1.160, 1.165, 1.195, 1.250, 1.305, 1.375, 1.455, 1.470, 1.485, 1.500,
        1.625, 1.680, 1.750, 1.795, 1.820, 1.875, 1.950, 1.980, 2.000, 2.050,
        2.100, 2.150, 2.200, 2.250, 2.300, 2.350, 2.400, 2.450, 2.500, 2.625,
        2.750, 2.875, 3.000, 3.250, 3.500, 4.000, 4.500, 5.000,
      };
      static const double sigma_[n] = {
        .0000, .0003, .0502, .1279, .1952, .3173, .4970, .5229, .5414, .5317,
        .2799, .2233, .1851, .1719, .1676, .1811, .2014, .2115, .2189, .2253,
        .2343, .2491, .2690, .2874, .2915, .2752, .2490, .2257, .2112, .2051,
        .2166, .2109, .1873, .1620, .1564, .1493, .1395, .1359
      }; // data computed using SOPHIA (arXiv:astro-ph/9903478v1)
      static double *y2a = NULL;
      if (!y2a) {
        y2a = new double[n];
        spline(s_,sigma_,n,1.e30,1.e30,y2a);
      }
      sigma = splint(s_,sigma_,y2a,n,s);
    } else {
      static const int n = 7;
      static const double logsqrts[n] = {
        .3495,
        .5044, 1., 2., 3., 4., 4.548
      };
      static const double logsigma[n] = { 
        -.867, // same as the last point above
        -.915, -.941, -.851, -.684, -.514, -.423
      }; // high-energy extrapolation (green line) in Fig. 46.9 of RPP 2012
         //  http://pdg.web.cern.ch/pdg/2012/reviews/rpp2012-rev-cross-section-plots.pdf
      static double *y2a = NULL;
      if (!y2a) {
        y2a = new double[n];
        spline(logsqrts,logsigma,n,1.e30,1.e30,y2a);
      }
      //assert(s <= 1.294e9);
      sigma = pow(10.,splint(logsqrts,logsigma,y2a,n,0.5*log10(s)));
    }
    if (sigma < 0.) sigma = 0.;  
  } else 
    sigma = 0.;
}

double losseA0pair(double z, double Gam)
{
//............................................................................. 
// Subroutine that computes the pair production energy losses for nuclei 
// expressed in 1/y. Only the CMB background is relevant
//.............................................................................

//........ pair production energy losses (only CMB) ......................

  double E = mN * Gam * (1+z); // corresponding proton energy (eV)
  double beta_pair = 0.;
  double lbeta;
  if (E <= 1.e17)
    beta_pair = 0.;
  else {
    if (E <= 1.e24) {
      static double *y2a = NULL;
      if (!y2a) {
        y2a = new double[nE_ee];
        spline(lEne_ee,lbeta0_ee,nE_ee,1.e30,1.e30,y2a);
      }
      lbeta = splint(lEne_ee,lbeta0_ee,y2a,nE_ee,log10(E));
    } else
      lbeta = -11.91234 - 0.61714286 * (log10(E)-24.);
    beta_pair = pow(1+z,3) * pow(10, lbeta);
  }
  return beta_pair;
}

double losseA0(double z, double Gam)
// pair production + pion production (approximated as continuous)
{
  double E = mN * Gam * (1+z);
  double beta = 0.;
  if (E > 1.e17) {
    static double *y2a = NULL;
    if (!y2a) {
      y2a = new double[nE_p];
      spline(LogE,lbeta0,nE_p,1.e30,1.e30,y2a);
    }
    double lbeta = splint(LogE,lbeta0,y2a,nE_p,log10(E));
    beta = pow(1+z,3) * pow(10, lbeta);
  }
  return beta;
}

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

//double fdisi_IR (double lx)
//{
//  // integral function for the IR energy losses
//  
//  double x = exp(lx);
//  double lxx = log10(x);
//  static double *y2a = NULL;
//  if (!y2a) {
//    y2a = new  double[nE_IR*nz_IR];
//    splie2(lE_IR,z_IR,ln_IR,nE_IR,nz_IR,y2a);
//  }
//  double lnn = splin2(lE_IR,z_IR,ln_IR,y2a,nE_IR,nz_IR,lxx,z_com);
//  
//  double ln=pow(10.0,lnn)/x;
//  double f_IR=x*(ln)/x/x;
//  
//  return f_IR;
//
//}

double NucModel::HadrRate(double Gam, double z)
// Interaction rate for pion photoproduction by protons, in yr^{-1}
// (that for nuclei is A times this)
// CMB only
{
  double rate = 0.;
  if (fStoch > 0 && Gam*mN*(1.+z) > 1.e19) {
    static double *y2a = NULL;
    if (!y2a) {
      y2a = new double[nE_p];
      spline(LogE_pi,lbeta_pi,nE_pi,1.e30,1.e30,y2a);
    }
    double lrate = splint(LogE_pi,lbeta_pi,y2a,nE_pi,log10(Gam*mN*(1.+z)));
    rate = pow(1.+z,3.)*pow(10.,lrate);
  }
  return rate;
}

double NucModel::HadrRate_IR(double Gam, double z)
// Interaction rate for pion photoproduction by protons, in yr^{-1}
// (that for nuclei is A times this)
// IBL only
{ 
  double rate = 0.;
  double epsmin = epsthr;
  double epsmax = 2.e-6*Gam*epsIBL_max;
  if (fLosses == 3) epsmax *= 1 + z;
  if (fStoch == 2 && epsmax > epsmin) {
    if (fLosses == 1) {
      if (z > 5.) z = 5.;
      static double *y2a = NULL;
      if (!y2a) {
        y2a = new double[nG_IR*nz_IR];
        splie2(lG_IR,z_IR,hrate_IR,nG_IR,nz_IR,y2a);
      }
      rate = pow(10.,splin2(lG_IR,z_IR,hrate_IR,y2a,nG_IR,nz_IR,log10(Gam),z));
    } else if (fLosses == 2 || fLosses == 3 || fLosses == 4 || fLosses == 5 || fLosses == 6 || fLosses == 7) {
      z_com = z;
      Gam_com = Gam;
      rate = qgaus(fpion_IR, log(epsmin), log(epsmax));
    }
  }
  return rate;
}

double NucModel::NuclRate(int j, double Gam, double z)
{
  gModel = this;
  s_com = 1.;
  k_com = j;
  z_com = z;
  Gam_com = Gam;
  return qgaus (fdisi, log(fNuclt[j]), log(fe_1)) + qgaus(fdisi,log(fe_1),log(fe_max));
}

double NucModel::AlphRate(int j, double Gam, double z)
{
  gModel = this;
  s_com = 3.;
  k_com = j;
  z_com = z;
  Gam_com = Gam;
  return qgaus (fdisi, log(fAlpht[j]), log(fe_1)) + qgaus(fdisi,log(fe_1),log(fe_max));
}


double NucModel::TotalRate(int A, int Z, int j, double Gam, double z, eProcType *type)
{
  int ndisDummy = -1;
  double rate_had = A*HadrRate(Gam,z);
  double rate_had_IR = A*HadrRate_IR(Gam,z);
  double rate_dis = 0, rate_nuc = 0, rate_alp = 0;
  if (fType == eBreitWigner2 || fType == eGauss2) {
    if (A > 1) rate_nuc = nFactor * NuclRate(j,Gam,z);
    if (A > 4) rate_alp = aFactor * AlphRate(j,Gam,z);
  } else
    rate_dis = losseAdisi(A,z,Gam,j,ndisDummy);
  double rate_tot = rate_had + rate_had_IR + rate_dis + rate_nuc + rate_alp;
  if (type != NULL) {
    double u = gRandom->Uniform(rate_tot);
    if (u < rate_had)
      *type = eHadron;
    else if (u < rate_had + rate_had_IR)
      *type = eHadron_IR;
    else if (u < rate_had + rate_had_IR + rate_dis)
      *type = eDisi;
    else if (u < rate_had + rate_had_IR + rate_dis + rate_nuc)
      *type = eNucleonProd;
    else
      *type = eAlphaProd;
  }
  return rate_tot;
}

 
vector<Particle> NucModel::PropagatePion(Particle * input)
{
  vector<Particle> output;
  int charge = input->GetCharge();
  double E = input->GetEprod();
  double z = input->GetZprod();
  int nBr = input->GetBranch() + 1;
  assert(E >= 0. && z >= 0.);
  double r = gRandom->Rndm();
  if (charge == 0) {
    output.push_back(Particle(0,0,nBr,r*E,z,ePhoton));
    output.push_back(Particle(0,0,nBr,(1.-r)*E,z,ePhoton));
    cout << "photon " << r*E <<"; photon " << (1.-r)*E << endl;
  } else {
    double E_nu = (1.-(mmu*mmu)/(mpi*mpi))*E*r;
    //nu_from_pion++;
    output.push_back(Particle(0,charge,nBr,E-E_nu,z,eMuon));
    output.push_back(Particle(0,0,nBr,E_nu,z,(charge > 0) ? eNeutrino_mu : eAntineutrino_mu));
    cout << "muon " << E-E_nu <<"; neutrino " << E_nu << endl;
  }
  input->SetEint(E);
  input->SetZint(z);
  input->SetIntMult(1002);
  return output;
}

vector<Particle> NucModel::PropagateMuon(Particle * input)
{
  int charge = input->GetCharge();
  double E = input->GetEprod();
  double z = input->GetZprod();
  int nBr = input->GetBranch() + 1;
  assert(E >= 0. && z >= 0.);
  double gamma = E/mmu;
  double E1, E2, E3, p1, p2, p3; //quantities in CoM frame
  // 1 = muon neutrino, 2 = electron neutrino, 3 = electron
  const double Enumax = (mmu - me*me/mmu)/2.;
  do {
    E1 = gRandom->Uniform(Enumax);
    E2 = gRandom->Uniform(Enumax);
    if (E1 + E2 < Enumax) {
      E1 = Enumax - E1;
      E2 = Enumax - E2;
    } // (E1, E2) uniform in triangle (0, Enumax), (Enumax, Enumax), (Enumax, 0)
    E3 = mmu - E1 - E2;
    p1 = E1;
    p2 = E2;
    p3 = sqrt(E3*E3 - me*me);
  } while (E3 < me || p3 > p1 + p2 || p1 > p2 + p3 || p2 > p3 + p1);

  double c12 = (p3*p3 - p1*p1 - p2*p2)/(2.*p1*p2); //angle btw neutrinos
  double s12 = sqrt(1. - c12*c12);
  double c1 = gRandom->Uniform(-1.,1.); //angle btw neutrino 1 and line of sight
  double s1 = sqrt(1. - c1*c1);
  double ca = cos(gRandom->Uniform(6.283185307)); //angle bw nu2 and (nu1, line of sight) plane
  double c2 = c12*c1 - s12*s1*ca;

  double Enu1 = gamma*(E1 + p1*c1); // in lab frame
  double Enu2 = gamma*(E2 + p2*c2);
  double Ee = E - Enu1 - Enu2;

  input->SetEint(E);
  input->SetZint(z);
  input->SetIntMult(1003);
  vector<Particle> output;
  output.push_back(Particle(0,0,nBr,Enu1,z,(charge > 0) ? eAntineutrino_mu : eNeutrino_mu));
  output.push_back(Particle(0,0,nBr,Enu2,z,(charge > 0) ? eNeutrino_e : eAntineutrino_e));
  output.push_back(Particle(0, charge, nBr, Ee, z, eElectron));
  cout << "electron " << Ee << "; neutrino " << Enu1 << "; neutrino " << Enu2 << endl;
  return output;
}

vector<Particle> NucModel::PropagateNeutrino(Particle * input)
{
  vector<Particle> output;
  double E = input->GetEprod();
  double z = input->GetZprod();
  input->SetEint(E/(1.+z));
  input->SetZint(0.);
  input->SetIntMult(0);
  return output;
}

vector<Particle> NucModel::PropagatePhoton(Particle *input)
{
  vector<Particle> output;
  double E = input->GetEprod();
  double zOri = input->GetZprod();
  input->SetEint(E);
  input->SetZint(zOri);
  input->SetIntMult(1000);
  return output;
}

vector<Particle> NucModel::PropagateElectron(Particle *input)
{
  vector<Particle> output;
  double E = input->GetEprod();
  double zOri = input->GetZprod();
  input->SetEint(E);
  input->SetZint(zOri);
  input->SetIntMult(1000);
  return output;
}

double NucModel::Photopion(double z, double Gam, eProcType proc)
{
// Restituisce l'energia di un pione fotoprodotto da un protone a energia Gam*M e redshift z
// (L'energia del protone sarà quella iniziale meno quella del pione.)
  double angle;
  double E = Gam*mN;
  double s;
  if (proc == eHadron_IR)
    s = sample_sIBL(E,z,angle);
  else {
    double eps = (1.+z)*sample_eps(E*(1.+z));
    s = sample_s(eps,E,angle);
  }
  double r_s = sqrt(s);
  double p_star = sqrt((s-(mpi+mN)*(mpi+mN))*(s-(mpi-mN)*(mpi-mN)))/(2.*r_s);
  double E_star = (s-mN*mN+mpi*mpi)/(2.*r_s);
  double cos_theta = gRandom->Uniform(-1.,1.); // supp. distribuz. isotropa nel CM
  double gamm = E/r_s; // da CM a Lab
  assert(gamm*(E_star+p_star*cos_theta) <= mN*Gam);
  return gamm*(E_star+p_star*cos_theta);
}

double NucModel::sample_sIBL(double E, double z, double& angle)
// samples CoM energy of photohadronic interaction on the IBL
// (complicated, maybe I'll write something about how it works later)
// Armando di Matteo,  2 May 2013
{
  double I, ph, eps, s;
  const double Gam = E/mN;
  double eps_min = std::max(epsthr/(2.e-6*Gam), epsIBL_min);
  //double s_thr = (mpi + mN)*(mpi + mN);
  const double I_max = integIBL(eps_min,z);
  const double I_min = std::max(integIBL(epsIBL_max,z), 1e-4*I_max);
  static const int N = 8;
  double I_[N + 1], phi_[N + 1];
  for (int i = 0; i <= N; i++) {
    I_[i] = ((N-i)*I_min + i*I_max)/N;
    phi_[i] = phi(mN*mN+4*E*inv_integIBL(I_[i],z));
  }
  double area_cum[N + 1];
  area_cum[0] = 0.;
  for (int i = 1; i <= N; i++) {
    area_cum[i] = area_cum[i - 1] + phi_[i - 1] * (I_[i] - I_[i - 1]);
  }
  do {
    double u = gRandom->Uniform(0., area_cum[N]);
    int i_min = 0, i_max = N;
    while (i_max - i_min > 1) {
      int i_mid = (i_min + i_max)/2;
      if (u > area_cum[i_mid]) i_min = i_mid;
      else i_max = i_mid;
    }
    I = gRandom->Uniform(I_[i_min], I_[i_max]);
    ph = gRandom->Uniform(0., phi_[i_min]);
    eps = inv_integIBL(I,z);
    s = phi_inv(ph);
  } while (s > mN*mN + 4*E*eps);
  angle = acos(1-(s-mN*mN)/(2*E*eps));
  return s;
}

double NucModel::GetDecay(int A, int Z)
{
    int N = A-Z;
    double Q = NAN;
    if (Z<32 && N<37)
      Q = qval[Z][N];
    if (isnan(Q)) { // estimate using semi-empirical mass formula
       int Zstable=GetBetaDecayStableZ(A);
       const double mp = 938.2720e6, mn = 939.5653e6, me = .5110e6;
       const double aC = 0.697e6, aA = 23.285e6, aP = 12.0e6;
       if (Z == Zstable) {
           Q = 0.;
       } else if (Z < Zstable) { // beta-minus decay
           Q = mn-mp-me-(2*Z+1)*pow((double)A,-1./3.)*aC+4.*aA*(double)(A-2*Z-1)/A;
           if (A % 2) {
               if (Z % 2) // odd-odd to even-even
                   Q += 2*aP*pow((double)A,-0.5);
               else // even-even to odd-odd
                   Q -= 2*aP*pow((double)A,-0.5);
           }
       } else { // beta-plus decay
           Q = mp-me-mn+(2*Z-1)*pow(A,-1./3.)*aC-4.*aA*(double)(A-2*Z+1)/A;
           if (A % 2) {
               if (Z % 2) // odd-odd to even-even
                   Q += 2*aP*pow((double)A,-0.5);
               else // even-even to odd-odd
                   Q -= 2*aP*pow((double)A,-0.5);
           }
           Q = -Q;
       }
    }
    return Q;
}

vector<Particle> NucModel::BetaDecay(Particle * input)
{
  int A = input->GetMassNum();
  int Z = input->GetCharge();
  double Gam = input->GetGprod();
  double z = input->GetZprod();
  int nBr = input->GetBranch() + 1;
  double Q = GetDecay(A, Z);
  int charge = -1; // beta-minus decay
  if (Q < 0.) {    // beta-plus decay
      Q = -Q;
      charge = +1;
  }
  double Ee, u;
  double max = sqrt((me+Q)*(me+Q) - me*me)*(me+Q)*Q*Q;
  do {
    Ee = gRandom->Uniform(me,me+Q);
    u = gRandom->Uniform(0.,max);
  } while (u > sqrt(Ee*Ee - me*me)*Ee*(me+Q-Ee)*(me+Q-Ee)); // neglecting nucleus recoil and EM effects
  double Enu = me+Q - Ee; // in CoM frame

  Ee = Gam*Ee*gRandom->Uniform(0.,2.); // (1+cos(theta))
  Enu = Gam*Enu*gRandom->Uniform(0.,2.); // in lab frame
  double En = input->GetEprod() - Ee - Enu;

  input->SetGint(Gam);
  input->SetZint(z);
  input->SetIntMult(1003);
  vector<Particle> output;
  output.push_back(Particle(A,Z-charge,nBr,En,z,(A > 1) ? eNucleus : eProton));
  output.push_back(Particle(0,0,nBr,Enu,z,(charge > 0) ? eNeutrino_e : eAntineutrino_e));
  output.push_back(Particle(0, charge, nBr, Ee, z, eElectron));
  cout << "nucleus " << En << "; electron " << Ee << "; neutrino " << Enu << endl;
  return output;
}

double qval[32][37] = {{NAN, 782347, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {-0, -0, 18591, 2.348e+07, 2.151e+07, 2.427e+07, 2.303e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, -0, -0, -0, 3.5083e+06, 1.1193e+07, 1.0651e+07, 1.5985e+07, 1.576e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {-1.374e+07, -2.29e+07, -290000, -0, -0, 1.60052e+07, 1.36066e+07, 2.0444e+07, 2.0623e+07, 2.502e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, -2.632e+07, -4.288e+06, -861890, -0, -0, 555900, 1.1506e+07, 1.1708e+07, 1.669e+07, 1.629e+07, 2.083e+07, 2.06e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, -2.523e+07, -1.21e+07, -1.79798e+07, -1.068e+06, -0, -0, 1.33689e+07, 1.34372e+07, 2.0644e+07, 1.9099e+07, 2.339e+07, 2.273e+07, 2.74e+07, 2.694e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, -1.2173e+07, -1.64948e+07, -3.64795e+06, -1.9824e+06, -0, -0, 156476, 9.7717e+06, 8.01e+06, 1.3167e+07, 1.181e+07, 1.656e+07, 1.579e+07, 2.071e+07, 2.124e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, -2.31e+07, -1.365e+07, -1.73381e+07, -2.22047e+06, -0, -0, 1.04207e+07, 8.68e+06, 1.3896e+07, 1.2527e+07, 1.797e+07, 1.719e+07, 2.275e+07, 2.378e+07, 2.847e+07, 2.906e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, -1.471e+07, -1.7767e+07, -5.14394e+06, -2.7542e+06, -0, -0, -0, 4.8223e+06, 3.8149e+06, 8.11e+06, 6.49e+06, 1.128e+07, 1.151e+07, 1.617e+07, 1.744e+07, 2.003e+07, 2.062e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, -2.465e+07, -1.392e+07, -1.5417e+07, -2.76051e+06, -1.6552e+06, -0, 7.02453e+06, 5.6842e+06, 1.0818e+07, 8.48e+06, 1.351e+07, 1.338e+07, 1.784e+07, 1.786e+07, 2.198e+07, 2.224e+07, 2.58e+07, 2.545e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, -1.3316e+07, -1.4509e+07, -4.4435e+06, -3.23883e+06, -0, -0, -0, 4.37581e+06, 2.4666e+06, 7.25e+06, 7.292e+06, 1.259e+07, 1.223e+07, 1.539e+07, 1.474e+07, 1.819e+07, 1.821e+07, 2.111e+07, 2.036e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.887e+07, -1.1175e+07, -1.389e+07, -3.5476e+06, -2.8423e+06, -0, 5.51545e+06, 3.835e+06, 9.352e+06, 9.069e+06, 1.4029e+07, 1.3284e+07, 1.7272e+07, 1.587e+07, 2.002e+07, 2e+07, 2.395e+07, 2.343e+07, 2.653e+07, 2.603e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, -2.011e+07, -1.0723e+07, -1.3095e+07, -4.7855e+06, -4.0561e+06, -0, -0, -0, 2.61001e+06, 1.8318e+06, 7.596e+06, 6.962e+06, 1.1736e+07, 1.011e+07, 1.342e+07, 1.174e+07, 1.628e+07, 1.564e+07, 1.93e+07, 1.895e+07, 2.217e+07, 2.094e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.521e+07, -1.858e+07, -1.2243e+07, -1.38766e+07, -4.2767e+06, -4.00427e+06, -0, 4.64236e+06, 3.6797e+06, 8.561e+06, 7.995e+06, 1.302e+07, 1.196e+07, 1.702e+07, 1.423e+07, 1.826e+07, 1.653e+07, 2.012e+07, 1.947e+07, 2.383e+07, 2.214e+07, 2.524e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.398e+07, -1.7e+07, -1.0812e+07, -1.274e+07, -5.066e+06, -4.81236e+06, -0, -0, -0, 1.49188e+06, 224310, 5.845e+06, 4.601e+06, 1.05e+07, 7.77e+06, 1.241e+07, 1.069e+07, 1.48e+07, 1.357e+07, 1.884e+07, 1.75e+07, 2.093e+07, 2.074e+07, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -2.124e+07, -1.505e+07, -1.812e+07, -1.1667e+07, -1.4334e+07, -4.9424e+06, -4.2324e+06, -0, 1.71048e+06, 248500, 5.374e+06, 3.9886e+06, 1.0413e+07, 7.9e+06, 1.21e+07, 1.029e+07, 1.476e+07, 1.374e+07, 1.862e+07, 1.773e+07, 2.122e+07, 2.116e+07, 2.481e+07, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.5e+07, -1.826e+07, -1.123e+07, -1.379e+07, -6.138e+06, -5.3962e+06, -0, -0, -0, 167180, 1.14222e+06, 4.86517e+06, 2.937e+06, 6.64e+06, 4.69e+06, 8.29e+06, 7.24e+06, 1.22e+07, 1.111e+07, 1.511e+07, 1.541e+07, 1.852e+07, 1.79e+07, 2.17e+07, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -2.248e+07, -1.63e+07, -1.851e+07, -1.198e+07, -1.2686e+07, -5.5826e+06, -5.49201e+06, -0, 709680, -0, 4.9165e+06, 3.442e+06, 7.48e+06, 5.76e+06, 9.51e+06, 7.84e+06, 1.244e+07, 1.141e+07, 1.501e+07, 1.539e+07, 1.901e+07, 1.844e+07, 2.181e+07, 2.129e+07, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.564e+07, -1.836e+07, -1.113e+07, -1.16193e+07, -6.0626e+06, -5.9661e+06, -0, -813870, -0, 565000, 1.50469e+06, 2.4916e+06, 599000, 4.583e+06, 3.14e+06, 6.838e+06, 5.7e+06, 9.79e+06, 8.41e+06, 1.217e+07, 1.085e+07, 1.421e+07, 1.32e+07, 1.66e+07, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -2.262e+07, -1.615e+07, -1.69e+07, -1.1879e+07, -1.2805e+07, -6.14746e+06, -5.91386e+06, -0, 1.31107e+06, -0, 3.52552e+06, 1.815e+06, 5.66e+06, 4.204e+06, 7.717e+06, 6.644e+06, 1.209e+07, 1.097e+07, 1.422e+07, 1.386e+07, 1.631e+07, 1.59e+07, 1.849e+07, 1.785e+07}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.463e+07, -1.577e+07, -1.099e+07, -1.1638e+07, -6.741e+06, -6.5326e+06, -0, -421310, -0, -0, -0, 255800, 1.378e+06, 1.992e+06, 282000, 5.2631e+06, 4.966e+06, 7.35e+06, 7.85e+06, 9.73e+06, 1.033e+07, 1.146e+07, 1.183e+07}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -2.034e+07, -1.6e+07, -1.712e+07, -1.3106e+07, -1.4323e+07, -6.49537e+06, -6.42583e+06, -2.2207e+06, -3.6524e+06, -0, 2.3663e+06, 600300, 3.992e+06, 2.006e+06, 6.89e+06, 6.51e+06, 9.11e+06, 9.21e+06, 1.138e+07, 1.209e+07, 1.367e+07, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.404e+07, -1.567e+07, -1.167e+07, -1.294e+07, -7e+06, -6.867e+06, -267600, -2.0621e+06, -0, -0, -0, -0, 2.2051e+06, 2.4735e+06, 1.976e+06, 5.02e+06, 4.3e+06, 7.48e+06, 7.14e+06, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.918e+07, -1.55e+07, -1.695e+07, -1.13e+07, -1.343e+07, -7.126e+06, -7.0504e+06, -2.93034e+06, -4.0123e+06, -601900, 1.0379e+06, -0, 3.9756e+06, 3.436e+06, 7.042e+06, 5.96e+06, 9.2e+06, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.416e+07, -1.589e+07, -1.066e+07, -1.291e+07, -7.599e+06, -7.444e+06, -1.656e+06, -2.6265e+06, -0, -752580, -0, -0, 1.3772e+06, 2.6031e+06, 1.6285e+06, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.986e+07, -1.385e+07, -1.71e+07, -1.23e+07, -1.35e+07, -7.715e+06, -7.63269e+06, -3.2075e+06, -4.7115e+06, -596800, 697100, -0, 3.69564e+06, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.869e+07, -1.313e+07, -1.564e+07, -1.116e+07, -1.303e+07, -8.15e+06, -8.019e+06, -2.374e+06, -3.7426e+06, -0, -231210, -0, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.733e+07, -1.98e+07, -1.501e+07, -1.728e+07, -1.295e+07, -1.442e+07, -8.3e+06, -8.24292e+06, -3.4518e+06, -4.566e+06, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.676e+07, -1.857e+07, -1.34e+07, -1.584e+07, -1.126e+07, -1.328e+07, -8.8e+06, -8.692e+06, -2.136e+06, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -2.003e+07, -1.591e+07, -1.752e+07, -1.371e+07, -1.53e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -1.513e+07, -1.67e+07, -1.287e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}, {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, -2.099e+07, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN}};
