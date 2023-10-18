#include "TOMPhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Decay.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4Proton.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ProductionCutsTable.hh"

TOMPhysicsList::TOMPhysicsList() : G4VModularPhysicsList()
{   
  // EM physics
  G4LossTableManager::Instance()->SetVerbose(0);
      
  G4double defaultCutValue1 = 0.000001*mm;
  G4double defaultCutValue2 = 0.0001*mm;

  cutForGamma     = defaultCutValue1;
  cutForElectron  = defaultCutValue1;
  cutForPositron  = defaultCutValue1;
  cutForProton    = defaultCutValue2;

  SetVerboseLevel(0);

  emPhysicsList = new G4EmStandardPhysics_option4;
}

TOMPhysicsList::~TOMPhysicsList()
{
  delete emPhysicsList;
}

void TOMPhysicsList::ConstructParticle()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();
  
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();  

  // mesons
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

  // barions
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();

  // ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
}

void TOMPhysicsList::ConstructProcess()
{
  AddTransportation();
  emPhysicsList->ConstructProcess();
  //G4EmParameters* param = G4EmParameters::Instance();
  //param->SetFluo(true); 
  //param->SetAuger(true); 
  //param->SetPixe(true); 
  SetCuts();
  
  AddStepMax();
}

void TOMPhysicsList::SetCuts()
{
  SetCutsWithDefault();
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(1*keV,10*GeV);
}

void TOMPhysicsList::AddStepMax()
{  
  G4StepLimiter* stepLimiter =new G4StepLimiter;
  auto theParticleIterator = GetParticleIterator();  
  theParticleIterator->reset();
  while ((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    pmanager ->AddDiscreteProcess(stepLimiter);
  }
}

