#ifndef LSExpActionInitialization_hh
#define LSExpActionInitialization_hh

#include "G4VUserActionInitialization.hh"


class LSExpActionInitialization : public G4VUserActionInitialization
{
public:
    LSExpActionInitialization();
    ~LSExpActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

    virtual G4VSteppingVerbose* InitializeSteppingVerbose() const;

};

#endif
