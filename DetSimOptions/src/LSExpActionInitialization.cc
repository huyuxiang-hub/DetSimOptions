#include "LSExpActionInitialization.hh"
#include "G4Event.hh"
#include "LSExpPrimaryGeneratorAction.hh"


LSExpActionInitialization::LSExpActionInitialization()
    : G4VUserActionInitialization()
{}

LSExpActionInitialization::~LSExpActionInitialization()
{}

void LSExpActionInitialization::BuildForMaster() const
{

}

void LSExpActionInitialization::Build() const
{
    // Thread Local
    // the registered actions are created during the initialization of slave run manager
    LSExpPrimaryGeneratorAction* pga = new LSExpPrimaryGeneratorAction();
    SetUserAction(pga);
}

G4VSteppingVerbose* LSExpActionInitialization::InitializeSteppingVerbose() const
{
  return 0;
}
