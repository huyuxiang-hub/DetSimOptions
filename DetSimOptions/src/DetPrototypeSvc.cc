#include "DetPrototypeSvc.hh"
#include "SniperKernel/SvcFactory.h"
#include "SniperKernel/ToolBase.h"

#include "LSExpDetectorConstruction.hh"
#include "LSExpPhysicsList.hh"
#include "LSExpEventAction.hh"
#include "LSExpRunAction.hh"
#include "LSExpSteppingAction.hh"
#include "LSExpTrackingAction.hh"
#include "LSExpPrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "DetSimAlg/MgrOfAnaElem.h"
#include "DetSimAlg/DetSimAlg.h"
#include "DetSimAlg/IAnalysisElement.h"

DECLARE_SERVICE(DetPrototypeSvc);

DetPrototypeSvc::DetPrototypeSvc(const std::string& name)
    : SvcBase(name)
{
    declProp("AnaMgrList", m_ana_list);
}

DetPrototypeSvc::~DetPrototypeSvc()
{

}

bool
DetPrototypeSvc::initialize()
{
    MgrOfAnaElem& mgr = MgrOfAnaElem::instance();
    // register
    if (m_ana_list.size() == 0) {
        m_ana_list.push_back("NormalAnaMgr");
    }
    // TODO
    // user can create the Ana Element more Early
    SniperPtr<DetSimAlg> detsimalg(*getParent(), "DetSimAlg");
    if (detsimalg.invalid()) {
        LogError << "Can't Load DetSimAlg" << std::endl;
        return false;
    }
    for (std::vector<std::string>::iterator it = m_ana_list.begin();
            it != m_ana_list.end();
            ++it) {
        ToolBase* t = 0;
        // find the tool first
        t = detsimalg->findTool(*it);
        //
        // create the tool if not exist
        if (not t) {
            t = detsimalg->createTool(*it);
        }
        //
        // register the tool into MgrOfAnaElem
        if (t) {
            IAnalysisElement* ae = dynamic_cast<IAnalysisElement*>(t);
            if (ae) {
                mgr.reg(*it, ae);
                continue;
            }
        } 
        LogError << "Can't Load Tool " << *it << std::endl;
        return false;
    }

    return true;
}

bool
DetPrototypeSvc::finalize()
{
    return true;
}

G4VUserDetectorConstruction* 
DetPrototypeSvc::createDetectorConstruction()
{
    LSExpDetectorConstruction* dc = new LSExpDetectorConstruction;
    dc->setScope(getParent());
    dc -> setCDName("Prototype");
    dc -> setPMTMother("lBuffer");

    return dc;
}

G4VUserPhysicsList* 
DetPrototypeSvc::createPhysicsList()
{
    return new LSExpPhysicsList;
}

G4VUserPrimaryGeneratorAction* 
DetPrototypeSvc::createPrimaryGenerator()
{
    LSExpPrimaryGeneratorAction* pga = new LSExpPrimaryGeneratorAction;
    pga->setScope(getParent());
    return pga;
}

G4UserRunAction*  
DetPrototypeSvc::createRunAction() 
{
    return new LSExpRunAction;
}

G4UserEventAction*  
DetPrototypeSvc::createEventAction() 
{
    return new LSExpEventAction;
}

G4UserStackingAction*  
DetPrototypeSvc::createStackingAction() 
{
    return NULL;
} 

G4UserTrackingAction*  
DetPrototypeSvc::createTrackingAction() 
{
    return new LSExpTrackingAction;
} 

G4UserSteppingAction*  
DetPrototypeSvc::createSteppingAction() 
{
    return new LSExpSteppingAction;
} 
