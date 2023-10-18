#include "TOMPrimaryGeneratorAction.hh"
#include "TOMRunAction.hh"
#include "TOMAnalysisManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4AutoLock.hh"

G4Mutex PrimGenMutex = G4MUTEX_INITIALIZER;

TOMPrimaryGeneratorAction::TOMPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{  
  //gun setup
  eParticleGun  = new G4ParticleGun(1); 	
  particle = particleTable->FindParticle(particleName="e-");
  eParticleGun->SetParticleDefinition(particle);
}

TOMPrimaryGeneratorAction::~TOMPrimaryGeneratorAction()
{
  G4AutoLock lock(&PrimGenMutex);
  delete eParticleGun;
}

void TOMPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  //this function is called at the begining of each event
  runAction = static_cast<const TOMRunAction*>(G4RunManager::GetRunManager() -> GetUserRunAction());
  eParticleGun->SetParticlePosition(G4ThreeVector(-10*mm,0,0));
  G4double Energy = 120*keV;
  eParticleGun->SetParticleEnergy(Energy);
  eParticleGun->GeneratePrimaryVertex(anEvent);
}

