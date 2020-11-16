#ifndef DetPrototypeSvc_hh
#define DetPrototypeSvc_hh

#include <SniperKernel/SvcBase.h>
#include <DetSimAlg/IDetSimFactory.h>
#include <vector>
#include <string>

class DetPrototypeSvc: public SvcBase, public IDetSimFactory {
public:
    DetPrototypeSvc(const std::string& name);
    ~DetPrototypeSvc();

    bool initialize();
    bool finalize();

    G4VUserDetectorConstruction* createDetectorConstruction();
    G4VUserPhysicsList* createPhysicsList();
    G4VUserPrimaryGeneratorAction* createPrimaryGenerator();

    G4UserRunAction*      createRunAction(); 
    G4UserEventAction*    createEventAction(); 
    G4UserStackingAction* createStackingAction();
    G4UserTrackingAction* createTrackingAction();
    G4UserSteppingAction* createSteppingAction();
private:
    std::vector<std::string> m_ana_list;

    std::string m_pmt_sd_type;
};

#endif
