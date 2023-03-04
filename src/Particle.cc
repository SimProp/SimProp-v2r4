#include "Particle.h"
#include <sstream>
#include <iostream>

using namespace std;

Particle::Particle() :
  fmassnum(-1), 
  fcharge(-1),
  fbranch(-1),
  fIntMult(0),
  fEprod(-1.),
  fZprod(-1.),
  fEint(-1.), 
  fZint(-1.),
//  fTtot(0.),
//  fdelta_fin(0.),
  fType(eNull)
{

}

Particle::~Particle()
{

}

Particle::Particle(int massnum, int charge, int branch, double Eprod, double Zprod, eParticleType type) : 
  fmassnum(massnum),
  fcharge(charge),
  fbranch(branch),
  fIntMult(0),
  fEprod(Eprod),
  fZprod(Zprod),
  fEint(-1.), 
  fZint(-1.),
//  fTtot(0.),
//  fdelta_fin(0.),
  fType(type)
{  

}

Particle::Particle(int massnum, int charge, int branch, double Eprod, double Zprod, double Eint, double Zint, eParticleType type) :
  fmassnum(massnum),
  fcharge(charge),
  fbranch(branch),
  fIntMult(0),
  fEprod(Eprod),
  fZprod(Zprod),
  fEint(Eint), 
  fZint(Zint),
//  fTtot(0.),
//  fdelta_fin(0.),
  fType(type)
{  

}

const char* Particle::GetName() const
{
  switch (fType) {
  case eNull: return "null";
  case eProton: return "proton";
  case eNeutron: return "neutron";
  case eNucleus: return "nucleus";
  case ePion: return "pion";
  case eMuon: return "muon";
  case ePhoton: return "photon";
  case eNeutrino_e: return "electron neutrino";
  case eAntineutrino_e: return "electron antineutrino";
  case eNeutrino_mu: return "muon neutrino";
  case eAntineutrino_mu: return "muon antineutrino";
  case eElectron: return fcharge > 0 ? "positron" : "electron";
  }
  return "null";
}

int Particle::GetFlavor() const
{
  switch (fType) {
  case eNeutrino_e: return 1;
  case eAntineutrino_e: return -1;
  case eNeutrino_mu: return 2;
  case eAntineutrino_mu: return -2;
  case eMuon: return fcharge > 0 ? -2 : +2;
  case eElectron: return fcharge > 0 ? -1 : +1;
  default: return 0;
  }
}

int Particle::GetID() const
// The numbering scheme used by CRPropa
// see https://crpropa.desy.de/Basic_Concepts#Particle_ID
{
  switch (fType) {
  case eNull:
   return 0;
  case eProton:
   return 2212;
  case eNeutron:
   return 2112;
  case eNucleus:
   return 1000000000 + 10000*fcharge + 10*fmassnum;
  case ePion:
   return fcharge ? 211*fcharge : 111;
  case eMuon:
   return 13 * -fcharge;
  case ePhoton:
   return 22;
  case eNeutrino_e:
   return 12;
  case eAntineutrino_e:
   return -12;
  case eNeutrino_mu:
   return 14;
  case eAntineutrino_mu:
   return -14;
  case eElectron:
   return fcharge > 0 ? -11 : 11;
  }
  return 0;
}
