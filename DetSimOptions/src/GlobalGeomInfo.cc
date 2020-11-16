#include <boost/python.hpp>
#include "G4SystemOfUnits.hh"
#include "GlobalGeomInfo.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperException.h"
#include "SniperKernel/ToolFactory.h"

#include "globals.hh"

DECLARE_TOOL(GlobalGeomInfo);

GlobalGeomInfo::GlobalGeomInfo(const std::string& name)
    : ToolBase(name)
{
    init_variables();

    declProp("LS.AbsLen", m_LS_abslen_at430);
    declProp("LS.RayleighLen", m_LS_raylen_at430);
    declProp("LS.LightYield", m_LS_light_yield=-1); // if LY<0, use default one.

    declProp("UseParamSvc", m_use_param_svc=false); // Using IMCParamsSvc
}

GlobalGeomInfo::~GlobalGeomInfo()
{

}

double
GlobalGeomInfo::geom_info(const std::string& param)
{
    // TODO use a map to replace these if else.
    if (param == "") {
        
    } else if (param == "ExpHall.XLen") {
        return m_expHallXLen;
    } else if (param == "ExpHall.YLen") {
        return m_expHallYLen;
    } else if (param == "ExpHall.ZLen") {
        return m_expHallZLen;
    } else if (param == "ExpHall.+X") {
        return m_expHallXmax;
    } else if (param == "ExpHall.-X") {
        return m_expHallXmin;
    } else if (param == "ExpHall.+Y") {
        return m_expHallYmax;
    } else if (param == "ExpHall.-Y") {
        return m_expHallYmin;
    } else if (param == "ExpHall.+Z") {
        return m_expHallZmax;
    } else if (param == "ExpHall.-Z") {
        return m_expHallZmin;

    } else if (param == "WaterPool.R") {
        return m_waterPoolRadius;
    } else if (param == "WaterPool.H") {
        return m_waterPoolHeight;
    } else if (param == "LS.AbsLen") {
        return m_LS_abslen_at430;
    } else if (param == "LS.RayleighLen") {
        return m_LS_raylen_at430;
    } else if (param == "LS.LightYield") {
        return m_LS_light_yield;
    } else if (param == "UseParamSvc") {
        return m_use_param_svc;
    } else {

    }
    throw SniperException(std::string("Can't find param: ") + param);
    return 0;
}

double
GlobalGeomInfo::geom_info(const std::string& param, int /* idx */)
{
    throw SniperException(std::string("Can't find param: ") + param);
    return 0;
}

void
GlobalGeomInfo::init_variables()
{
    // == Experiment Hall ==
    m_expHallXLen = 48*m;
    m_expHallYLen = 48*m; // 55.25m in the future.
    m_expHallZLen = 18.6*m; // height

    m_expHallXmax = m_expHallXLen/2.;
    m_expHallYmax = m_expHallYLen/2.;
    m_expHallZmax = m_expHallZLen/2.;

    m_expHallXmin = -m_expHallXLen/2.;
    m_expHallYmin = -m_expHallYLen/2.;
    m_expHallZmin = -m_expHallZLen/2.;
    // == Water Pool ==
    m_waterPoolRadius = 43.5*m/2;
    m_waterPoolHeight = 43.5*m;
    // == Central Detector ==
    //
    // == material ==
    m_LS_abslen_at430 = 77.;
    m_LS_raylen_at430 = 27.;
}
