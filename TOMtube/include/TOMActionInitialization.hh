#ifndef TOMActionInitialization_h
#define TOMActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class TOMActionInitialization : public G4VUserActionInitialization
{
  public:
    TOMActionInitialization();
    ~TOMActionInitialization() override;
    void BuildForMaster() const override;
    void Build() const override;
};
#endif
