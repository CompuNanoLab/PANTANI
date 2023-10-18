#ifndef TOMRunAction_h
#define TOMRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

class G4Run;

class TOMRunAction : public G4UserRunAction
{
  public:
    TOMRunAction();
    ~TOMRunAction();
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
};
#endif

