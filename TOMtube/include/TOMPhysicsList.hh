#ifndef TOMPhysicsList_h
#define TOMPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;

class TOMPhysicsList: public G4VModularPhysicsList
{
  public:
    TOMPhysicsList();
    virtual ~TOMPhysicsList();
    void ConstructParticle(); 
    void ConstructProcess();
    void AddStepMax();      
    void SetCuts(); 
  private:
    G4VPhysicsConstructor*  emPhysicsList;   
    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;    
    G4double cutForProton;
};
#endif

