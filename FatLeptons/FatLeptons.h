#ifndef FatLeptons_FatLeptons_H
#define FatLeptons_FatLeptons_H

#include <EventLoop/Algorithm.h>

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

#include "xAODTruth/TruthParticleContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODCaloEvent/CaloCluster.h"

#include <FatLeptons/HistoManager.h>

namespace CP {
    class EgammaCalibrationAndSmearingTool;
}

namespace TrigConf {
      class xAODConfigTool;
}

namespace Trig {
      class TrigDecisionTool;
}

class AsgElectronLikelihoodTool;
class JetCalibrationTool;
class GoodRunsListSelectionTool;

class Minitree;

class FatLeptons : public EL::Algorithm {
    public:
        FatLeptons();

        virtual EL::StatusCode setupJob(EL::Job& job);
        virtual EL::StatusCode fileExecute(){return EL::StatusCode::SUCCESS;}
        virtual EL::StatusCode histInitialize(){return EL::StatusCode::SUCCESS;}
        virtual EL::StatusCode changeInput(bool){return EL::StatusCode::SUCCESS;}
        virtual EL::StatusCode initialize();
        virtual EL::StatusCode execute();
        virtual EL::StatusCode postExecute(){return EL::StatusCode::SUCCESS;}
        virtual EL::StatusCode finalize(){return EL::StatusCode::SUCCESS;}
        virtual EL::StatusCode histFinalize(){return EL::StatusCode::SUCCESS;}

        void prepareHistograms();
        const xAOD::TruthParticle * getTruthZ(xAOD::Jet *j);
        int getNumberOfTruthMatchedElectrons(const xAOD::TruthParticleContainer *tps, const xAOD::ElectronContainer *electrons, std::vector<const xAOD::Electron*> &goodElectrons);
        int getNMatched(const xAOD::Electron *el, const xAOD::CaloClusterContainer *cls);
        std::vector<const xAOD::CaloCluster*> getMatchedCaloClusters(const xAOD::Electron *el, const xAOD::CaloClusterContainer *cls);
        int getNMatchedEls(const xAOD::Jet *jet, const xAOD::ElectronContainer *els);
        void doMatchingStudy(const xAOD::Jet *jet, const xAOD::ElectronContainer *els, const xAOD::CaloClusterContainer *cls);
        void doRecoStudy(const xAOD::Jet *jet);

        ClassDef(FatLeptons, 1);

        // executable options
        std::string outputName;

        // local variables
        xAOD::TEvent *m_event; //!
        xAOD::TStore *m_store; //!
        int m_eventCounter; //!
        bool m_isMC; //!
        int m_eventNumber; //!
        int m_runNumber; //!
        std::string m_derivation; //!

        // local tools
#ifndef __CINT__
        CP::EgammaCalibrationAndSmearingTool *m_egammaCalibrationTool; //!
        AsgElectronLikelihoodTool *m_electronLikelihoodTool; //!
        JetCalibrationTool *m_jetCalibrationTool; //!
        GoodRunsListSelectionTool *m_grl; //!
        TrigConf::xAODConfigTool *m_trigConfigTool; //!
        Trig::TrigDecisionTool *m_trigDecisionTool; //!
        HistoManager *m_hm; //!
        Minitree *mt; //!
#endif

};

#endif
