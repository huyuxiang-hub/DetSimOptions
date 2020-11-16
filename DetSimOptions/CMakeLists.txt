include (PKG)

PKG(DetSimOptions
    DEPENDS
        AnalysisCode

        CentralDetector
        TopTracker
        Chimney
        CalibUnit

        PhysiSim
        PMTSim
        SimUtil
        MCParamsSvc

        DetSimAlg

        ${GEANT4_LIBRARIES}
)
