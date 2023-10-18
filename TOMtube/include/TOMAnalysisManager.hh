#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4VProcess.hh"
#include <map>
#include <variant>

class G4Step;

class TOMAnalysisManager
{
  public:
    virtual ~TOMAnalysisManager();
    static TOMAnalysisManager* getInstance();
    void book();
    void save_pho(std::vector<G4int> out_int, std::vector<G4String> out_str, std::vector<G4double> out_dou);
    void CountProcesses(const G4VProcess* process);    
    void count();
    void finish();
  private:
    TOMAnalysisManager();
    G4String outputFileName;
    static TOMAnalysisManager* instance;
};

#endif

