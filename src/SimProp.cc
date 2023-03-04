//
// File   : SimProp.cc
//
// Purpose: Program to simulate propagation of nuclei
//
// Authors: D. Boncioli, A. di Matteo, A. Grillo, S. Petrera, F. Salamida
//
// Based on the analytical calculation by R. Aloisio 
//      (arXiv:1006.2484 [astro-ph.CO] and arXiv:0802.4452 [astro-ph])
// 
// date: 10/09/2010
// Last update: 14/11/2016
 
int VERSION = 2;
int REVISION = 4;

#include <iomanip>  
#include <vector>
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <list>

#include <cmath>
#include <time.h>

//ROOT headers
#include <TROOT.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include <Particle.h>
#include <NucModel.h>
#include <Output.h> 

//double delta_fin,Ttot; 

int Nevts = 100;
unsigned int Seed = 65539;
int Areq = 56;
double Emin = 1.e17;
double Emax = 1.e21;
double gamma_inj = 1.;
double zmin = 0.0;
double zmax = 1.0;
double Ls=0.; // distance between sources
int LossType = 1; // EBL model
eModelType ModelType = eSalomon; // nuclear model
double nFactor = 1.;
double aFactor = 1.;
int Bdecay = 1;
int Stoch = 1; //stochastic pion prod.
int PairProd = 0;
//double B = 1e-20;
//double Lcoh = 1;
//double Bs=1e-20, Rsources=1;
int OutType = 0;

int NatEarth[56] = {0};
int Photons = 0, Electrons = 0, Positrons = 0;
int NuatEarth[4] = {0};
using namespace std;

void Follow(vector<Particle>*, NucModel*, Output*);

bool GetOptions(int, char **);

void PrintSummary(int);


int main(int argc,char*argv[])
{
  if (!GetOptions(argc, argv))
    return 0;
 
  cout << '\n';
  cout << "=====================================================================\n";
  cout << "SimProp v" << VERSION << "r" << REVISION ;

  cout << '\n';
 
  bool random = (Areq == 0);
  bool proton = (Areq == 1);
  
  gRandom->SetSeed(Seed);
  
//  std::vector<double> magfield;
//  magfield.push_back(B);
//  magfield.push_back(Lcoh);
  NucModel* model = new NucModel(ModelType,LossType,Stoch,Bdecay,PairProd);
//  model->SetMag(magfield);
  int Zreq = 1;
  if (!proton && Areq > 0) {
    Zreq = model->TestPrimary(Areq);
    if(Zreq<0) {
      cout<<"** Initial nucleus not in the list"<<endl;
      exit(EXIT_FAILURE);
    }
  }
  if (Areq < 0) {
    cout<<"** Wrong choice for requested A"<<endl;
    exit(EXIT_FAILURE);
  }
  cout << "** Initial nucleus is: " << Areq << endl;
  cout << "** Losses are: ";
  cout << " CMB";
  if (LossType == 0) cout << " only";
  else if (LossType == 1) cout << " + EBL (Stecker et al. 2006)"; 
  else if (LossType == 2) cout << " + EBL (power law)"; 
  else if (LossType == 3) cout << " + EBL (Kneiske et al. 2004)";
  else if (LossType == 4) cout << " + EBL (Dominguez et al. 2011 best)";
  else if (LossType == 5) cout << " + EBL (Dominguez et al. 2011 lower)";
  else if (LossType == 6) cout << " + EBL (Dominguez et al. 2011 upper)";
  else if (LossType == 7) cout << " + EBL (Gilmore et al. 2012 fiducial)";
  else {
    cout<<"** Wrong choice for EBL model"<<endl;
    exit(EXIT_FAILURE);
  }
  cout << endl;
  
  cout<<"** Limits on redshift are: "<<zmin<<" "<<zmax<<endl; 
  cout<<"** Limits on energy   are: "<<Emin<<" eV "<<Emax<<" eV"<<endl;
  
  char FileOut[150];
  sprintf(FileOut,"SimProp-v%dr%d_N%d_A%d_e%4.1f_E%4.1f_g%1.2f_z%4.2f_Z%4.2f_r%4.2f_L%d_M%d_n%1.2f_a%1.2f_D%d_S%d_p%d_o%d_s%d.root",
    VERSION, REVISION, Nevts, Areq, log10(Emin), log10(Emax), gamma_inj, zmin,
    zmax, Ls, LossType, (int)ModelType - 2, nFactor, aFactor, Bdecay, Stoch,
    PairProd, OutType, Seed);
  Output *myout = new Output(FileOut, OutType, PairProd);
  cout << "** Output files is: " << FileOut << endl;
  cout << endl;
  cout << "** Events to be generated: " << Nevts << endl;


  
  //...... Loop over events ........
  double t_tot = 0.;
  clock_t start, end;
  for(int event = 0; event < Nevts; event++) { 
    int branchxev = 0;
    start = clock();
    cout << "** Event " << event << " ****************************************" << endl;
    //cout << "===> Seed: " <<   gRandom->GetSeed() << endl;
    double Zini = gRandom->Uniform(zmin, zmax);
    // mod by afg 23/04/12
    // to allow for sources at discrete positions with relative distance Ls
    // up to z=0.1
    if(Ls > 0. && Zini < 0.1) {
        const double deltaz = Ls/4285.; 
    	const int ddz = (int)(Zini/deltaz) + 1; 
    	Zini = (double)ddz * deltaz;
    }
    // end mod

    const double Eini = (gamma_inj == 1)
     ? exp(gRandom->Uniform(log(Emin), log(Emax)))
     : pow(gRandom->Uniform(pow(Emin, 1.-gamma_inj), pow(Emax, 1.-gamma_inj)),
           1./(1.-gamma_inj));

    Particle primary;
    if (proton)
        primary = Particle(1, 1, 0, Eini, Zini, eProton);
    else if (random)
        primary = model->GetRandomParticle(Eini, Zini);
    else
        primary = Particle(Areq, Zreq, 0, Eini, Zini, eNucleus);
    cout << "==> Primary nucleus (A,Z,E,z): " << primary.GetMassNum() << " " 
	 << primary.GetCharge() << " " << primary.GetEprod() << " " 
	 << primary.GetZprod() << endl;
    myout->Injec(event, primary);
    vector<Particle> stack;
    stack.push_back(primary);
    while (!stack.empty()) {
      Follow(&stack, model, myout);
      branchxev++;
      //cout << "---------------------------- size of stack " << stack.size() << endl;
    }

    end = clock();
    const double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout << "** event ended - elapsed time (s) "<<cpu_time_used<<endl;
    t_tot += cpu_time_used;
    unsigned long int CurrentSeed = (unsigned long int)gRandom->GetSeed();
    myout->Stream(cpu_time_used, branchxev, CurrentSeed);
    myout->Write();
  }
 
  delete myout;

  cout << endl;
  cout << "  ++++  Succesfully processed " << Nevts << " events ++++" << endl;
  cout << " --------------------------------------------------" << endl; 
  PrintSummary(Areq);
  cout << endl << "Total time elapsed: " << t_tot << endl;
  return 0;
} 



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Follow(vector<Particle> *mystack, NucModel* mymodel, Output* myout) 
{
  
  Particle it = mystack->front();
  mystack->erase(mystack->begin());
  // if(it.GetType() == eProton) {
  //   mymodel->PropagateProton(&it);
  // }
  // else {
  vector<Particle> products = mymodel->PropagateParticle(&it);
  if(!products.empty()) { 
    mystack->insert(mystack->end(),products.begin(),products.end());
  }
      //}
  myout->Stream(it);
  //cout << " it.GetZint() = " << it.GetZint() << " it.GetMassNum() = " << it.GetMassNum() << endl;  
  if(it.GetZint()==0) {
    switch (it.GetType()) {
      case eProton: case eNeutron: case eNucleus:
        NatEarth[it.GetMassNum()-1]++; break;
      case eNeutrino_e:
        NuatEarth[0]++; break;
      case eAntineutrino_e:
        NuatEarth[1]++; break;
      case eNeutrino_mu:
        NuatEarth[2]++; break;
      case eAntineutrino_mu:
        NuatEarth[3]++; break;
      default:
        break;
    }
    myout->Earth(it);
  }
  if (it.GetIntMult() == 1000) {
    switch (it.GetType()) {
     case ePhoton:
      Photons++;
      break;
     case eElectron:
      if (it.GetCharge() > 0) Positrons++;
      else Electrons++;
      break;
     default:
      break;
    }
    myout->Prod(it);
  }
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool GetOptions(int argc, char **argv)
{
  int c;
  ostringstream help;
  bool PrintHelp = false;
  help << " SimProp -N [number of events, D=100]  -s [random seed, D=65539]\n"
       << "         -A [A (0 for all), D=56]  -e [min energy log10(Emin/eV), D=17.]\n"
       << "         -E [max energy log10(Emax/eV), D=21.]\n"
       << "         -g [injection spectral index, D=1.]\n"
       << "         -z [min redshift, D=0.] -Z [max redshift, D=1.]\n"
       << "         -r [source relative distance (Mpc), D=0.]\n"
       << "         -L [EBL, 0=none, 1=Stecker+ '06, 2=power law, 3=Kneiske+ '04,\n"
       << "             4=Dominguez+ '11 (best), 5=Dominguez (lower), 6=Dominguez (upper),\n"
       << "             7=Gilmore+ '12 (fiducial), D=1]\n"
       << "         -M [nuclear model type, 0=Stecker-Salamon, 1=arb. Gaussians, \n"
       << "             2=arb. Breit-Wigner, 3=arb. Breit-Wigner with alpha, \n"
       << "             4=arb. Gausssians with alpha, D=0]\n"
       << "         -n [scales nucleon photodis. cross sections (only with M>=3), D=1.] \n"
       << "         -a [scales alpha photodis. cross sections (only with M>=3), D=1.] \n"
       << "         -D [beta decay, 0=none, 1=istantaneous, D=1]\n"
       << "         -S [pion photoproduction, -1=continuous for protons only, \n"
       << "             0=continuous for all nuclei, 1=stochastic on CMB only, \n"
       << "             2=stochastic on CMB+EBL, D=1]\n"
//       << "         -p [1=e+e- pairs written in output (WARNING: very large files!), \n"
//       << "             2=only total energy written, 0=neglected, D=0]\n"
       << "         -p [e+e-, 0=neglected, 1=individually written to output (WARNING: very\n"
       << "             large files!), 2=binned in z, 3=only total energy written, D=0]\n"
//       << "         -B [vacuum mag field (nG), D=0] -c [coherence lenght (Mpc), D=1]\n"
//       << "         -b [source mag field (nG), D=0] -d [magnetic dimension of the source (Mpc), D=10]\n"
       << "         -o [output types, 0=all branches, 1=summary tree, 2=both, D=0]\n";
  while ((c = getopt (argc, argv, "N:s:A:e:E:g:z:Z:r:L:M:n:a:D:S:p:o:h")) != -1)
//  while ((c = getopt (argc, argv, "N:s:A:e:E:g:z:Z:r:L:M:n:a:D:S:p:B:c:b:d:o:h")) != -1)

    switch (c)
      {
      case 'N':
	Nevts=atoi(optarg);
	break;
      case 's':
	Seed=atoi(optarg);
	break;
      case 'A':
	Areq=atoi(optarg);
	break;
      case 'e':
	Emin=pow(10,atof(optarg));
	break;
      case 'E':
	Emax=pow(10,atof(optarg));
	break;
      case 'g':
        gamma_inj=atof(optarg);
        break;
      case 'z':
	zmin=atof(optarg);
	break;
      case 'Z':
	zmax=atof(optarg);
	break;
      case 'r':
	Ls=atof(optarg);
	break;
      case 'L':
	LossType=atoi(optarg);
	break;
      case 'M':
        ModelType=eModelType(atoi(optarg)+2);
        break;
      case 'n':
        nFactor=atof(optarg);
        break;
      case 'a':
        aFactor=atof(optarg);
        break;
      case 'D':
	Bdecay=atoi(optarg);
	break;
      case 'S':
	Stoch=atoi(optarg);
	break;
      case 'p':
        PairProd=atoi(optarg);
        break;
//      case 'B':
//	B=atof(optarg);
//	break;
//      case 'c':
//	Lcoh=atof(optarg);
//	break;
//      case 'b':
//	Bs=atof(optarg);
//	break;
//      case 'd':
//	Rsources=atof(optarg);
//	break;
      case 'o':
        OutType=atoi(optarg);
        break;
      case 'h':
	PrintHelp=true;
	break;
      }

  if(PrintHelp) {
    cout << help.str() << endl;
    return false;
  }

  cout << "\n          >>  SimProp v" << VERSION << "r" << REVISION << "  <<\n\n"
       <<"-------------------------------------------------------------------------------\n"
       <<"  Authors: D. Boncioli, A. di Matteo, A.F. Grillo, S. Petrera and F. Salamida  \n" 
       <<"-------------------------------------------------------------------------------\n"
       << "  events             : " << Nevts << "\n"
       << "  random seed        : " << Seed << "\n" 
       << "  A                  : " << Areq << "\n"
       << "  Emin               : " << Emin << " eV\n"
       << "  Emax               : " << Emax << " eV\n"
       << "  inj. spectr. index : " << gamma_inj << "\n"
       << "  zmin               : " << zmin << "\n"
       << "  zmax               : " << zmax << "\n"
       << "  dist. btw. sources : " << Ls << "\n"
       << "  EBL model          : " << LossType << "\n"
       << "  nuclear model      : " << (int)ModelType - 2 << "\n"
       << "  nucleon ej. scaling: " << nFactor << "\n"
       << "  alpha ej. scaling  : " << aFactor << "\n"
       << "  beta decay         : " << Bdecay << "\n"
       << "  pion photoprod.    : " << Stoch << "\n"
       << "  pair photoprod.    : " << PairProd << "\n"
//       << "  B                  : " << B << "\n"
//       << "  Lcoh               : " << Lcoh << "\n"
//       << "  Bs                 : " << Bs << "\n" 
       << "  output type        : " << OutType << "\n";
  
  return true;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PrintSummary(int A)
{
 
  printf("       Nucleus Type   Number reaching Earth\n");
  printf(" --------------------------------------------------\n"); 
  int Ntot = 0;
  for (int i=0; i<A; i++) {
    Ntot += NatEarth[i];
    if(NatEarth[i]>0) 
      printf("      %8d     %14d\n",i+1,NatEarth[i]);
      //cout << "Nucleus n." << i+1 << " --> " << NatEarth[i] << endl;
  }
  printf(" --------------------------------------------------\n"); 
  printf("         Total     %14d\n\n",Ntot);
  printf(" --------------------------------------------------\n");
  printf(" photons from pi0 decay: %14d\n",Photons);
  printf(" electrons produced:     %14d\n",Electrons);
  printf(" positrons produced:     %14d\n",Positrons);
  printf(" --------------------------------------------------\n");
  printf("       Neutrino Type     Number reaching Earth\n");
  printf(" --------------------------------------------------\n");
  printf(" electron neutrinos:     %14d\n",NuatEarth[0]);
  printf(" electron antineutrinos: %14d\n",NuatEarth[1]);
  printf(" muon neutrinos:         %14d\n",NuatEarth[2]);
  printf(" muon antineutrinos:     %14d\n",NuatEarth[3]);
  printf(" --------------------------------------------------\n"); 
  printf("         Total           %14d\n\n",NuatEarth[0]+NuatEarth[1]+NuatEarth[2]+NuatEarth[3]);
}
