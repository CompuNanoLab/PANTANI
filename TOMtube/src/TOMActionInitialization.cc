#include "TOMActionInitialization.hh"
#include "TOMPrimaryGeneratorAction.hh"
#include "TOMRunAction.hh"
#include "TOMEventAction.hh"
#include "TOMSteppingAction.hh"

TOMActionInitialization::TOMActionInitialization() : G4VUserActionInitialization()
{}

TOMActionInitialization::~TOMActionInitialization()
{}

void TOMActionInitialization::BuildForMaster() const
{
  SetUserAction(new TOMRunAction);
}

void TOMActionInitialization::Build() const
{
  SetUserAction(new TOMPrimaryGeneratorAction);
  SetUserAction(new TOMRunAction);
  TOMEventAction* eventAction = new TOMEventAction;
  SetUserAction(eventAction);
  SetUserAction(new TOMSteppingAction(eventAction));
}
