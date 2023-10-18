#include "TOMSteppingAction.hh"
#include "TOMEventAction.hh"
#include "TOMDetectorConstruction.hh"
#include "TOMAnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "TOMPrimaryGeneratorAction.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"
#include <variant>

TOMSteppingAction::TOMSteppingAction(TOMEventAction* eventAction)
: G4UserSteppingAction()
{}

TOMSteppingAction::~TOMSteppingAction()
{}

void TOMSteppingAction::UserSteppingAction(const G4Step* step)
{
  TOMAnalysisManager* analysis = TOMAnalysisManager::getInstance();
  
  // End step
  if((step -> GetPostStepPoint() -> GetPhysicalVolume() == NULL) )
  {
    step->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }
  
  // Save data
  std::vector<G4int> out_int;
  std::vector<G4double> out_dou;
  std::vector<G4String> out_str;    
  if ((step->GetTrack()->GetDynamicParticle()->GetDefinition()) == G4Gamma::Definition() )
  { 
    if (step->GetPostStepPoint()->GetPosition()[2]<-27*mm &&\
        step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "World" &&\
        step->GetPostStepPoint()->GetMomentum()[2]<0)//
    {      
      out_int.push_back(G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
      out_int.push_back(step->GetTrack()->GetTrackID());
      out_int.push_back(step->GetTrack()->GetParentID());
      out_int.push_back(step->GetTrack()->GetCurrentStepNumber());
      out_str.push_back(step->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName());   
      if (step->GetTrack()->GetTrackID() != 1 ) 
      {
        out_str.push_back(step->GetTrack()->GetCreatorProcess()->GetProcessName());
      }
      else
      {
        out_str.push_back("primary");
      }
      if (step->GetTrack()->GetCurrentStepNumber() != 1 ) 
      {
        out_str.push_back(step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName());
      }
      else
      {
        out_str.push_back("start");
      }
      out_str.push_back(step->GetPreStepPoint()->GetPhysicalVolume()->GetName());
      out_str.push_back(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
      out_str.push_back(step->GetPostStepPoint()->GetPhysicalVolume()->GetName());      
      out_dou.push_back(step->GetTrack()->GetGlobalTime());
      out_dou.push_back(step->GetTrack()->GetLocalTime());
      out_dou.push_back(step->GetTrack()->GetTrackLength());
      out_dou.push_back(step->GetTrack()->GetStepLength());
      out_dou.push_back(step->GetPreStepPoint()->GetKineticEnergy()); 
      out_dou.push_back(step->GetPostStepPoint()->GetKineticEnergy());        
      out_dou.push_back(step->GetPreStepPoint()->GetPosition()[0]);
      out_dou.push_back(step->GetPreStepPoint()->GetPosition()[1]);
      out_dou.push_back(step->GetPreStepPoint()->GetPosition()[2]);
      out_dou.push_back(step->GetPreStepPoint()->GetMomentum()[0]);
      out_dou.push_back(step->GetPreStepPoint()->GetMomentum()[1]);
      out_dou.push_back(step->GetPreStepPoint()->GetMomentum()[2]);
      out_dou.push_back(step->GetPostStepPoint()->GetPosition()[0]);
      out_dou.push_back(step->GetPostStepPoint()->GetPosition()[1]);
      out_dou.push_back(step->GetPostStepPoint()->GetPosition()[2]);
      out_dou.push_back(step->GetPostStepPoint()->GetMomentum()[0]);
      out_dou.push_back(step->GetPostStepPoint()->GetMomentum()[1]);
      out_dou.push_back(step->GetPostStepPoint()->GetMomentum()[2]);
      analysis->save_pho(out_int, out_str, out_dou);
    }
  }
}

