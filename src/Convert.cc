#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

TFile * f;
TTree * nuc;
 
int currev = 0;
double z1[200],z2[200];
double E1[200],E2[200];
int A[200], Z[200], intm[200];

int GetEvent(int &, int &);
 
int main(int argc,char*argv[])
{

  if(argc<2) {
    cerr << "   USAGE: Convert SimPropOutput.root" << endl;
    exit( 1 );
  }
  char filename[130];
  strcpy(filename, argv[1]);
  gSystem->Load("libTree.so");
  f = new TFile(filename);
  nuc = (TTree*)f->Get("nuc");

  cout<<"Root file name is: "<< filename << endl;

  string outname(filename);
  size_t found = outname.find_last_not_of("root");
  outname = outname.substr(0,found);
  outname.append(".out");

  ofstream outFile(outname.c_str()); 
  if ( !outFile ) {
    cerr << "Output File could not be opened." << endl;
    exit( 1 );
  }
  cout<<"Output Ascii file name is: "<< outname <<endl;
  int nentries = (int)(nuc->GetEntries());
  cout << "Entries are = " << nentries <<endl;

  int kontev = 0;
  int maxev = nentries;
 
  int nall  = 0;
  int numev = 0;
  int i     = 0;

  while (i < nentries) {

    nall = GetEvent(i, numev);
    kontev++;

    for (int nb = 0; nb < nall; nb++) {
      if (intm[nb] <= 0) 
	outFile<<z1[0]<<" "<<E1[0]<<" "<<A[0]<<" "<<Z[0]<<" "
	       <<E2[nb]<<" "<<A[nb]<<" "<<Z[nb]<<endl;
    }

    if (kontev >= maxev) break;

  }
  cout << "Events are = " << kontev <<endl;
  outFile.close();
  return 0;
}
 
//////////////////////////////////////////////////////////////////////////////////

int GetEvent(int &i, int &numev) {


  int evt;
  int branch;
  int intmult;
  int Acurr, Zecurr;
  double zOri, zEnd;
  double EOri, EEnd;
  double Dist; 

  int nb = 0;

  while (1) {
    
    if (i == 0) {
  
      nuc->SetBranchAddress("evt",&evt);
      nuc->SetBranchAddress("branch",&branch);
      nuc->SetBranchAddress("intmult",&intmult);
      nuc->SetBranchAddress("Acurr",&Acurr); 
      nuc->SetBranchAddress("Zecurr",&Zecurr); 
      nuc->SetBranchAddress("zOri",&zOri);
      nuc->SetBranchAddress("zEnd",&zEnd);
      nuc->SetBranchAddress("EOri",&EOri);
      nuc->SetBranchAddress("EEnd",&EEnd);
      nuc->SetBranchAddress("Dist",&Dist);

//       nuc->GetEntry(nentries-1);
//       cout << "Entries are = " << nentries <<endl;
    } 
    nuc->GetEntry(i);
    if (evt != numev) {
      numev = evt;
      break;
    }
    currev = evt;
    A[nb] = Acurr;
    Z[nb] = Zecurr;
    intm[nb] = intmult;
    z1[nb] = zOri;
    z2[nb] = zEnd;
    E1[nb] = EOri;
    E2[nb] = EEnd;
    nb++;
    i++;
    if (i == (int)(nuc->GetEntries())) break;
    
  }
  

  return nb;
}
