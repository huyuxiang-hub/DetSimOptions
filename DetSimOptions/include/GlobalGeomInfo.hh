#ifndef GlobalGeomInfo_hh
#define GlobalGeomInfo_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IDetElement.h"
#include "globals.hh"

class GlobalGeomInfo: public IDetElement,
                      public ToolBase {
public:
    G4LogicalVolume* getLV() {return 0;};
    bool inject(std::string /* motherName */, IDetElement* /* other */, IDetElementPos* /* pos */) {return false;}

    // FIXME: 
    // how the user check the info is right or not???
    double geom_info(const std::string& /* param */ );
    double geom_info(const std::string& /* param */ , int /* idx */);

    GlobalGeomInfo(const std::string& /* name */);
    ~GlobalGeomInfo();

private:
    void init_variables();

private:
    // = World =
    // == Experiment Hall ==
    double m_expHallXLen;
    double m_expHallYLen;
    double m_expHallZLen;
    double m_expHallXmax;
    double m_expHallYmax;
    double m_expHallZmax;
    double m_expHallXmin;
    double m_expHallYmin;
    double m_expHallZmin;
    //
    // == Water Pool ==
    double m_waterPoolRadius;
    double m_waterPoolHeight;
    // == Central Detector ==

    // == material ==
    double m_LS_abslen_at430;
    double m_LS_raylen_at430;
    double m_LS_light_yield;

    // = load from other service =
    bool m_use_param_svc;
};

#endif
