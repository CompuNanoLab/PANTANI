#ifndef TOMPrimaryGeneratorAction_h
#define TOMPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <fstream>
#include "G4ParticleTable.hh"

class G4ParticleGun;
class G4Event;
class G4Box;
class TOMFileReader;
class TOMRunAction;

class TOMPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    TOMPrimaryGeneratorAction();
    virtual~TOMPrimaryGeneratorAction();
    virtual void GeneratePrimaries(G4Event*);
    const G4ParticleGun* GeteParticleGun() const { return eParticleGun; }

  private:
    const TOMRunAction* runAction; 
    G4ParticleGun* eParticleGun = nullptr; // pointer a to G4 gun class  
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle;
};
#endif
