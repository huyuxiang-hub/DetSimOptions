#ifndef DetSim0Svc_hh
#define DetSim0Svc_hh

#include <SniperKernel/SvcBase.h>
#include <DetSimAlg/IDetSimFactory.h>
#include <vector>
#include <string>

class DetSim0Svc: public SvcBase, public IDetSimFactory {
public:
    DetSim0Svc(const std::string& name);
    ~DetSim0Svc();

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
    std::string m_muon_output_file;
    bool m_muon_merge_all;

    bool m_hits_merge;
    double m_hits_timewindow;

    std::string m_pmt_sd_type;
    double m_pmt_qe_scale_for_elec;

    std::string m_cd_name;
    std::string m_pmt_name;
    std::string m_extra_lpmt;
    std::string m_3inchpmt_name;
    std::string m_pmt_mother;
    std::string m_pmt_pos_mode;
    std::string m_pmt_pos_file;
    std::string m_3inch_pmt_pos_file;
    double m_3inch_pmt_pos_offset;

    std::string m_strut_name;
    std::string m_strut_mother;
    std::string m_strut_pos_file;
    
    std::string m_fastener_name;
    std::string m_fastener_mother;
    std::string m_fastener_pos_file;

    std::string m_xjanchor_name;
    std::string m_xjanchor_mother;
    std::string m_xjanchor_pos_file;

    std::string m_xjfixture_name;
    std::string m_xjfixture_mother;
    std::string m_xjfixture_pos_file;

    std::string m_sjclsanchor_name;
    std::string m_sjclsanchor_mother;
    std::string m_sjclsanchor_pos_file;

    std::string m_sjreceiver_name;
    std::string m_sjreceiver_mother;
    std::string m_sjreceiver_pos_file;

    std::string m_sjfixture_name;
    std::string m_sjfixture_mother;
    std::string m_sjfixture_pos_file;


    // enable CD
    bool m_cd_enabled;

    //veto
    std::string m_veto_pmt_name;
    std::string m_veto_pmt_pos_mode;
    // enable WP 
    bool m_wp_enabled;

    //define the name of the top tracker tool
    //(to enable TT sim at run time)
    std::string m_tt_name;
    bool m_tt_enabled;

    // Calib Unit related
    std::string m_cu_name; // calib unit name
    std::vector<std::string> m_cu_extras; // CalibUnitExtras
    bool m_cu_enable;      // enable or disable CU

    //(to enable Chimney sim at run time)/////////////////////////////////////
    std::string m_chimney_top_name;
    std::string m_chimney_lower_name;
    bool m_top_chimney_enable;      // enable or disable top Chimney 
    bool m_lower_chimney_enable;      // enable or disable lower Chimney 

    // == physics lists related
    // Default is JUNO. it is migrated from dyb.
    // If you want to use the geant4 pre-packaged physics lists,
    // just set its name.
    std::string m_physics_lists_name;

private:
    bool initializeOpticks() ; 
    bool finalizeOpticks() ;
    void dumpOpticks(const char* msg) ;  
    int m_opticksMode ; 


private:
    int m_GdLSAbsLengthMode;
    int m_turn_off;

};

#endif
