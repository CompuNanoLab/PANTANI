#include "G4AnalysisManager.hh"
#include "G4VProcess.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VDataSetAlgorithm.hh"
#include "TOMAnalysisManager.hh"
#include "G4WorkerThread.hh"
#include <variant>

TOMAnalysisManager* TOMAnalysisManager::instance = 0;

TOMAnalysisManager::TOMAnalysisManager()
  :outputFileName("Data.root")
{
  G4AnalysisManager::Instance();
}

TOMAnalysisManager::~TOMAnalysisManager() 
{
  delete instance;
  delete G4AnalysisManager::Instance();  
}

TOMAnalysisManager* TOMAnalysisManager::getInstance()
{
  if (instance == 0) 
  {
    instance = new TOMAnalysisManager;
  }
  return instance;
}

void TOMAnalysisManager::book()
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->SetVerboseLevel(0);
  man -> OpenFile(outputFileName);   
  std::vector<G4String> to_save_int={"EventID","TrackID","ParentID","Step_Number"};
  std::vector<G4String> to_save_str={"Name","Creator_Process","Process_Prestep","Volume_Prestep","Process_Poststep","Volume_Poststep"};
  std::vector<G4String> to_save_dou={  "Global_Time","Local_Time","Track_Length","Step_Length","Kinetic_Energy_Prestep","Kinetic_Energy_Poststep",\
  "Position_x_Prestep","Position_y_Prestep","Position_z_Prestep",\
  "Momentum_px_Prestep","Momentum_py_Prestep","Momentum_pz_Prestep",\
  "Position_x_Poststep","Position_y_Poststep","Position_z_Poststep",\
  "Momentum_px_Poststep","Momentum_py_Poststep"  ,"Momentum_pz_Poststep"};
  
  //PHOTONS
  man->CreateNtuple("PHO", "Data");
  for (G4int i=0;i<to_save_int.size();i++)
  {
    man->CreateNtupleIColumn(0,to_save_int[i]);
  }
  for (G4int i=0;i<to_save_str.size();i++)
  {
    man->CreateNtupleSColumn(0,to_save_str[i]);
  }
  for (G4int i=0;i<to_save_dou.size();i++)
  {
    man->CreateNtupleDColumn(0,to_save_dou[i]);
  }
  man->FinishNtuple(0);
}

void TOMAnalysisManager::save_pho(std::vector<G4int> out_int, std::vector<G4String> out_str, std::vector<G4double> out_dou)
{
  G4int i,j,k;
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  for (i = 0 ; i < out_int.size(); i++)
  {
    man->FillNtupleIColumn(0, i, out_int[i]);
  }
  for (j = 0 ; j < out_str.size(); j++)
  {
    man->FillNtupleSColumn(0, j+out_int.size(), out_str[j]);
  }  
  for (k = 0 ; k < out_dou.size(); k++)
  {
    man->FillNtupleDColumn(0, k+out_str.size()+out_int.size(), out_dou[k]);
  }
  man->AddNtupleRow(0);
}

void TOMAnalysisManager::finish()
{ 
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
}

