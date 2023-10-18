#include "TOMRunAction.hh"
#include "TOMPrimaryGeneratorAction.hh"
#include "TOMAnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

TOMRunAction::TOMRunAction()
 : G4UserRunAction()
{
  TOMAnalysisManager* analysis = TOMAnalysisManager::getInstance();
  analysis -> book();
}

TOMRunAction::~TOMRunAction()
{}

void TOMRunAction::BeginOfRunAction(const G4Run*)
{}

void TOMRunAction::EndOfRunAction(const G4Run* run)
{
  TOMAnalysisManager* analysis = TOMAnalysisManager::getInstance();  
  analysis -> finish();

}

