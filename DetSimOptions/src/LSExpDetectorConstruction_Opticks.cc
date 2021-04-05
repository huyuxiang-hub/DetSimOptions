#include <vector>
#include <string>

#include "LSExpDetectorConstruction_Opticks.hh"

#ifdef WITH_G4OPTICKS
#include "G4PVPlacement.hh"

#include "junoSD_PMT_v2.hh"

#include "G4Opticks.hh"
#include "PMTEfficiencyTable.hh"
#include "PLOG.hh"
#include "SSys.hh"
#endif

#ifdef WITH_G4OPTICKS
/**
LSExpDetectorConstruction_Opticks::Setup
------------------------------------------

1. pass geometry to Opticks, translate it to GPU and return sensor placements 
2. use the placements to pass sensor data : efficiencies, categories, identifiers
3. pass theta dependent efficiency tables for all sensor categories




                             |--------- 2230 ----------------|-- 120--|
                             20050                           17820    17700
                          / /                               /         /
                         / /                               /         /
                        / pInnerWater                     /         /
                       / /                               /         /
                      / /                  (0)          /         /
                     pTyvek                  \         pAcylic   /
                    / /                       \       /         /
                   / /                         \     /         pTarget:LS
                  / /                           \   /         /
                 / /                             \ /         /
                / /                              (1)        /
               / /                               / \       /
              / /                               /   \     /
             / /                               /     \   /         
            / /                               /       \ /
           / /                          Wa   /  Ac    (2)             
          / /                               /         / \
         / /                               /         /   \
        / /                               /         /     \        LS    

**/


G4Opticks* LSExpDetectorConstruction_Opticks::Setup(const G4VPhysicalVolume* world, const G4VSensitiveDetector* sd_, int opticksMode )  // static
{
    if( opticksMode == 0 ) return nullptr ; 
    LOG(info) << "[ WITH_G4OPTICKS opticksMode " << opticksMode  ; 

    assert(world); 

    // 1. pass geometry to Opticks, translate it to GPU and return sensor placements  

    G4Opticks* g4ok = new G4Opticks ; 

    bool outer_volume = true ; 
    bool profile = true ; 

    // waymask selects whether to match on pvname and/or boundary for the recorder boundary position
    const char* geospecific_default = "--way --pvname pAcylic  --boundary Water///Acrylic --waymask 3" ; // (1): gives radius 17820
    //const char* geospecific_default = "--way --pvname pTarget  --boundary Acrylic///LS --waymask 3" ;      // (2): gives radius 17700
    const char* embedded_commandline_extra = SSys::getenvvar("LSXDC_GEOSPECIFIC", geospecific_default ) ; // see OGeo::initWayControl

    LOG(info) << " embedded_commandline_extra " << embedded_commandline_extra ;
    g4ok->setPlacementOuterVolume(outer_volume); 
    g4ok->setProfile(profile); 
    g4ok->setEmbeddedCommandLineExtra(embedded_commandline_extra);
    g4ok->setGeometry(world); 

    const std::vector<G4PVPlacement*>& sensor_placements = g4ok->getSensorPlacements() ;       
    unsigned num_sensor = sensor_placements.size(); 

    // 2. use the placements to pass sensor data : efficiencies, categories, identifiers  

    const junoSD_PMT_v2* sd = dynamic_cast<const junoSD_PMT_v2*>(sd_) ;  
    assert(sd) ; 

    LOG(info) << "[ setSensorData num_sensor " << num_sensor ; 
    for(unsigned i=0 ; i < num_sensor ; i++)
    {
        const G4PVPlacement* pv = sensor_placements[i] ; // i is 0-based unlike sensor_index
        unsigned sensor_index = 1 + i ; // 1-based 
        assert(pv); 
        G4int copyNo = pv->GetCopyNo();  
        int pmtid = copyNo ; 
        int pmtcat = sd->getPMTCategory(pmtid); 
        float efficiency_1 = sd->getQuantumEfficiency(pmtid); 
        float efficiency_2 = sd->getEfficiencyScale() ; 

        g4ok->setSensorData( sensor_index, efficiency_1, efficiency_2, pmtcat, pmtid ); 
    }
    LOG(info) << "] setSensorData num_sensor " << num_sensor ; 

    // 3. pass theta dependent efficiency tables for all sensor categories 

    PMTEfficiencyTable* pt = sd->getPMTEfficiencyTable(); 
    assert(pt);

    const std::vector<int>& shape = pt->getShape(); 
    const std::vector<float>& data = pt->getData(); 

    int   theta_steps = pt->getThetaSteps();  
    float theta_min = pt->getThetaMin(); 
    float theta_max = pt->getThetaMax();
    LOG(info) 
         << "[ setSensorAngularEfficiency "
         << " theta_steps " << theta_steps 
         << " theta_min " << theta_min
         << " theta_max " << theta_max
         ; 
  
    g4ok->setSensorAngularEfficiency(shape, data, theta_steps, theta_min, theta_max); 
    LOG(info) << "] setSensorAngularEfficiency " ;

    g4ok->saveSensorLib("$TMP/LSExpDetectorConstruction__SetupOpticks/SensorLib") ; // just for debug 

    LOG(info) << "] WITH_G4OPTICKS " ; 
    return g4ok ; 
}
#else
G4Opticks* LSExpDetectorConstruction_Opticks::Setup(const G4VPhysicalVolume*, const G4VSensitiveDetector*, int )  // static
{
    return nullptr ; 
}
#endif

