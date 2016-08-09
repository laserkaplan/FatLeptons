#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>

#include <FatLeptons/FatLeptons.h>
#include <FatLeptons/Minitree.h>

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/xAODTruthHelpers.h"
#include "xAODJet/JetConstituentVector.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

#include "ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "JetCalibTools/JetCalibrationTool.h"
#include "JetMomentTools/JetVertexTaggerTool.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"

const static float ZMASSPDG = 91187.6;

struct caloClusterComp {
    const xAOD::CaloCluster *refCC;
    bool operator() (const xAOD::CaloCluster *c1, const xAOD::CaloCluster *c2) {return (c1->p4().DeltaR(refCC->p4()) < c2->p4().DeltaR(refCC->p4()));}
};

bool compPt(const xAOD::IParticle *p1, const xAOD::IParticle *p2) {return (p1->pt() > p2->pt());}

ClassImp(FatLeptons)

FatLeptons::FatLeptons() {
}

EL::StatusCode FatLeptons::setupJob(EL::Job &job) {
    job.useXAOD();
    xAOD::Init();
    
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode FatLeptons::initialize() {
    m_event = wk()->xaodEvent();
    m_store = wk()->xaodStore();

    m_eventCounter = 0;

    // set derivation name
    TTree *md = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
    md->LoadTree(0);
    if (!md->GetBranch("StreamAOD")) {
        for (auto br : *md->GetListOfBranches()) {
            std::string brName = br->GetName();
            if (brName.find("StreamDAOD") != std::string::npos) {
                m_derivation = brName.substr(brName.find("_") + 1);
                break;
            }
        }
    }
    std::cout << "Derivation = " << m_derivation << std::endl;

    // EgammaCalibrationAndSmearingTool
    /*
    m_egammaCalibrationTool = new CP::EgammaCalibrationAndSmearingTool("EgammaCalibrationAndSmearingTool");
    m_egammaCalibrationTool->setProperty("ESModel", "es2015PRE");
    m_egammaCalibrationTool->setProperty("decorrelationModel", "1NP_v1");
    m_egammaCalibrationTool->initialize();
    */

    // AsgElectronLikelihoodTool
    /*
    m_electronLikelihoodTool = new AsgElectronLikelihoodTool("AsgElectronLikelihoodTool");
    m_electronLikelihoodTool->setProperty("primaryVertexContainer", "PrimaryVertices");
    m_electronLikelihoodTool->setProperty("ConfigFile", "ElectronPhotonSelectorTools/offline/mc15_20150712/ElectronLikelihoodMediumOfflineConfig2015.conf");
    m_electronLikelihoodTool->initialize();
    */

    // JetCalibrationTool
    // LASER: note that this is only for data!
    /*
    m_jetCalibrationTool = new JetCalibrationTool("JetCalibrationTool", "AntiKt4EMTopo", "JES_MC15Prerecommendation_April2015.config", "JetArea_Residual_Origin_EtaJES_GSC_Insitu", true);
    m_jetCalibrationTool->initialize();
    */

    // GoodRunsListTool
    m_grl = new GoodRunsListSelectionTool("GoodRunsListTool");
    std::vector<std::string> grlvec;
    grlvec.push_back("$ROOTCOREBIN/data/FatLeptons/data15_13TeV.periodAllYear_DetStatus-v79-repro20-02_DQDefects-00-02-02_PHYS_StandardGRL_All_Good_25ns.xml");
    grlvec.push_back("$ROOTCOREBIN/data/FatLeptons/data16_13TeV.periodAllYear_DetStatus-v80-pro20-08_DQDefects-00-02-02_PHYS_StandardGRL_All_Good_25ns.xml");
    m_grl->setProperty("GoodRunsListVec", grlvec);
    m_grl->initialize();

    // xAODConfigTool
    /*
    m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool");
    m_trigConfigTool->initialize();
    ToolHandle<TrigConf::ITrigConfigTool> handle(m_trigConfigTool);
    handle->initialize();
    */
    
    // TrigDecisionTool
    /*
    m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
    m_trigDecisionTool->setProperty("ConfigTool", handle);
    m_trigDecisionTool->setProperty("TrigDecisionKey", "xTrigDecision");
    m_trigDecisionTool->initialize();
    */
    
    // HistoManager
    m_hm = new HistoManager();
    m_hm->setWeight(1.0);
    prepareHistograms();
    for (auto h : m_hm->getHistos()) wk()->addOutput(h);

    // Minitree
    TFile *outFile = wk()->getOutputFile(outputName);
    mt = new Minitree();
    mt->setDirectory(outFile);

    return EL::StatusCode::SUCCESS;
}

EL::StatusCode FatLeptons::execute() {
    m_store->clear();

    m_eventCounter++;
    if (m_eventCounter % 1000 == 0) std::cout << "Event: " << m_eventCounter << std::endl;

    // retrieve containers
    const xAOD::EventInfo *event_info = 0;
    if (!m_event->retrieve(event_info, "EventInfo").isSuccess()) return EL::StatusCode::FAILURE;
    if (event_info->eventType(xAOD::EventInfo::IS_SIMULATION)) {
        m_runNumber = event_info->mcChannelNumber();
        m_isMC = true;
    }
    else {
        m_runNumber = event_info->runNumber();
        m_isMC = false;
    }
    m_eventNumber = event_info->eventNumber();

    const xAOD::ElectronContainer *electrons = 0;
    if (!m_event->retrieve(electrons, "Electrons").isSuccess()) return EL::StatusCode::FAILURE;
    
    const xAOD::JetContainer *jets = 0;
    if (!m_event->retrieve(jets, "AntiKt4EMTopoJets").isSuccess()) return EL::StatusCode::FAILURE;
   
    // LASER - CANNOT DO IN JETM1!!!
    const xAOD::CaloClusterContainer *cls = 0;
    if (!m_event->retrieve(cls, "CaloCalTopoClusters").isSuccess()) return EL::StatusCode::FAILURE;

    // MC ONLY: set weight and fill sumOfWeights
    if (m_isMC) {
        if (m_eventCounter == 1 && m_derivation != "") {
            const xAOD::CutBookkeeperContainer *cbcs = 0;
            if (!m_event->retrieveMetaInput(cbcs, "CutBookkeepers").isSuccess()) return EL::StatusCode::FAILURE;

            const xAOD::CutBookkeeper *cbc = 0;
            for (const auto &cbk : *cbcs) {
                if (cbk->name() == m_derivation + "Kernel") cbc = cbk;
            }
            if (!cbc) return EL::StatusCode::FAILURE;

            m_hm->setWeight(cbc->sumOfEventWeights());
            m_hm->fill("sumOfWeights", 0);
        }

        m_hm->setWeight(event_info->mcEventWeight());
    }

    // LASER: get number of truth matched reco electrons
    // MC ONLY
    /*
    const xAOD::TruthParticleContainer *tps = 0;
    if (!m_event->retrieve(tps, "TruthParticles").isSuccess()) return EL::StatusCode::FAILURE;
    */
    
    /* 
    // basic electron selection, pt > 25 GeV, |eta| < 2.47, medium quality
    // LASER: for now, remove quality cut and reduce pt cut to 7 GeV
    xAOD::ElectronContainer *goodElectrons = new xAOD::ElectronContainer();
    xAOD::ElectronAuxContainer *goodElectronsAux = new xAOD::ElectronAuxContainer();
    goodElectrons->setStore(goodElectronsAux);
    int nElectrons = 0;
    for (const xAOD::Electron *el : *electrons) {
        xAOD::Electron *calibEl = new xAOD::Electron();
        calibEl->makePrivateStore(*el);
        m_egammaCalibrationTool->applyCorrection(*calibEl);
        //if ((calibEl->pt() / 1000. > 25.) && (fabs(calibEl->caloCluster()->etaBE(2)) < 2.47) && (m_electronLikelihoodTool->accept(calibEl))) {
        if ((calibEl->pt() / 1000. > 7.) && (fabs(calibEl->caloCluster()->etaBE(2)) < 2.47)) {
            goodElectrons->push_back(calibEl);
            nElectrons++;
        }
        else delete calibEl;
    }

    m_store->record(goodElectrons, "goodElectrons");
    m_store->record(goodElectronsAux, "goodElectronsAux.");
    */

    // basic jet selection, pt > 500 GeV, |eta| < 2.4
    // NB: no calibration applied!
    std::vector<const xAOD::Jet*> goodJets;
    int nJets = 0;
    for (const xAOD::Jet *j : *jets) {
        if ((j->jetP4("JetEMScaleMomentum").Pt() / 1000. > 500.) && (fabs(j->jetP4("JetEMScaleMomentum").Eta()) < 2.4)) {
            goodJets.push_back(j);
            nJets++;
            j->auxdecor<int>("emfraccut") = 0;
            j->auxdecor<int>("ntrackscut") = 0;
            float emfrac;
            std::vector<int> ntracks;
            j->getAttribute(xAOD::JetAttribute::EMFrac, emfrac);
            j->getAttribute(xAOD::JetAttribute::NumTrkPt1000, ntracks);
            if (emfrac > 0.97) j->auxdecor<int>("emfraccut") = 1;
            if (ntracks[0] < 5) j->auxdecor<int>("ntrackscut") = 1;
        }
    }

    // LASER: truth electron selection
    // MC ONLY
    /*
    std::vector<const xAOD::Electron*> goodElectrons;
    int nElectrons = getNumberOfTruthMatchedElectrons(tps, electrons, goodElectrons);
    if (nElectrons == -1) return EL::StatusCode::SUCCESS;
    */

    // LASER: calo cluster study
    // MC ONLY
    /*
    if (nElectrons == 1) {
        std::vector<const xAOD::CaloCluster*> ccs = getMatchedCaloClusters(goodElectrons[0], cls);
        m_hm->fill("nMatchedTopo2", ccs.size());
        //TLorentzVector vec = goodElectrons[0]->caloCluster()->p4();
        TLorentzVector vec = ccs[0]->p4();
        if (ccs.size() > 1) {vec += ccs[1]->p4(); m_hm->fill("mcc_2", vec.M() / 1000.);}
        if (ccs.size() > 2) {vec += ccs[2]->p4(); m_hm->fill("mcc_3", vec.M() / 1000.);}
        if (ccs.size() > 3) {vec += ccs[3]->p4(); m_hm->fill("mcc_4", vec.M() / 1000.);}
        if (ccs.size() > 4) {vec += ccs[4]->p4(); m_hm->fill("mcc_5", vec.M() / 1000.);}
    }
    */

    // LASER: truth studies for Christos (kill me now)
    // MC ONLY
    /*
    if (nElectrons == 1) {
        m_hm->fill("Christos_nTracks1TruthElectron", goodElectrons[0]->nTrackParticles());
        m_hm->fill("Christos_nMatchedTopo", getNMatched(goodElectrons[0], cls));
    }
    */

    // GRL cut
    if (!(m_isMC)) {
        if (!(m_grl->passRunLB(*event_info))) return EL::StatusCode::SUCCESS;
    }

    // Trigger cut
    // LASER: remove for now
    /*
    if (!(m_isMC)) {
        if (!(m_trigDecisionTool->isPassed("HLT_j360"))) return EL::StatusCode::SUCCESS;
    }
    */

    // fill histograms pre cuts
    // m_hm->fill("nElectrons", nElectrons);
    m_hm->fill("nJets", nJets);
    for (auto j : goodJets) {
        float mj = j->jetP4("JetEMScaleMomentum").M() / 1000.;
        float ptj = j->jetP4("JetEMScaleMomentum").Pt() / 1000.;
        float emfrac;
        std::vector<int> ntracks;
        j->getAttribute(xAOD::JetAttribute::EMFrac, emfrac);
        j->getAttribute(xAOD::JetAttribute::NumTrkPt1000, ntracks);
        int ntopoc = j->getConstituents().size();
        m_hm->fill("mj", mj);
        m_hm->fill("ptj", ptj);
        m_hm->fill("ntracksj", ntracks[0]);
        m_hm->fill("emfracj", emfrac);
        m_hm->fill("ntopocj", ntopoc);

        // fill minitree
        fatLeptonsEvent evt;
        evt.mj = mj;
        evt.ptj = ptj;
        evt.emfrac = emfrac;
        evt.ntracksj = ntracks[0];
        evt.nmatchedels = getNMatchedEls(j, electrons);
        evt.weight = (m_isMC) ? event_info->mcEventWeight() : 1.0;
        mt->fillTree(evt);
    }

    /*
    // electrons
    if (nElectrons >= 2) {
        float mdiff = 100000.;
        float curMass = 0, curPt = 0;
        for (auto e1 : goodElectrons) {
            for (auto e2 : goodElectrons) {
                if (e1 == e2) continue;
                TLorentzVector cur = e1->p4() + e2->p4();
                if (fabs(cur.M() - ZMASSPDG) < mdiff) {
                    curMass = cur.M() / 1000.;
                    curPt = cur.Pt() / 1000.;
                    mdiff = fabs(cur.M() - ZMASSPDG);
                }
            }
        }
        if (curMass != 0) {
            m_hm->fill("mee", curMass);
            m_hm->fill("ptee", curPt);
        }
    }
    */

    // jets - inclusive
    if (nJets != 0) {
        for (auto j : goodJets) {
            // variables to fill
            float mj = j->jetP4("JetEMScaleMomentum").M() / 1000.;
            float ptj = j->jetP4("JetEMScaleMomentum").Pt() / 1000.;
            int ntopoc = j->getConstituents().size();
            float emfrac;
            std::vector<int> ntracks;
            j->getAttribute(xAOD::JetAttribute::EMFrac, emfrac);
            j->getAttribute(xAOD::JetAttribute::NumTrkPt1000, ntracks);

            // postemfrac
            if (j->auxdecor<int>("emfraccut")) {
                m_hm->fill("mj_postemfrac", mj);
                m_hm->fill("ptj_postemfrac", ptj);
                m_hm->fill("ntracksj_postemfrac", ntracks[0]);
                m_hm->fill("ntopocj_postemfrac", ntopoc);
            }
            
            // postntracks
            if (j->auxdecor<int>("ntrackscut")) {
                m_hm->fill("mj_postntracks", mj);
                m_hm->fill("ptj_postntracks", ptj);
                m_hm->fill("emfracj_postntracks", emfrac);
                m_hm->fill("ntopocj_postntracks", ntopoc);
            }
            
            // postemfracntracks
            if (j->auxdecor<int>("emfraccut") && j->auxdecor<int>("ntrackscut")) {
                m_hm->fill("mj_postemfracntracks", mj);
                m_hm->fill("ptj_postemfracntracks", ptj);
                m_hm->fill("ntopocj_postemfracntracks", ntopoc);
            }
            
            // failemfrac
            if (!(j->auxdecor<int>("emfraccut"))) {
                m_hm->fill("mj_failemfrac", mj);
                m_hm->fill("ptj_failemfrac", ptj);
                m_hm->fill("ntracksj_failemfrac", ntracks[0]);
                m_hm->fill("ntopocj_failemfrac", ntopoc);
            }
            
            // failntracks
            if (!(j->auxdecor<int>("ntrackscut"))) {
                m_hm->fill("mj_failntracks", mj);
                m_hm->fill("ptj_failntracks", ptj);
                m_hm->fill("emfracj_failntracks", emfrac);
                m_hm->fill("ntopocj_failntracks", ntopoc);
            }
            
            // failemfracntracks
            if (!(j->auxdecor<int>("emfraccut") && j->auxdecor<int>("ntrackscut"))) {
                m_hm->fill("mj_failemfracntracks", mj);
                m_hm->fill("ptj_failemfracntracks", ptj);
                m_hm->fill("ntopocj_failemfracntracks", ntopoc);
            }
        }
    }

    // LASER - data study on matching electrons
    // CANNOT DO ON JETM1!!
    if (nJets != 0) {
        for (auto j : goodJets) {
            if (j->auxdecor<int>("emfraccut")) {
                doMatchingStudy(j, electrons, cls);
                doRecoStudy(j);
            }
        }
    }

    /*     
    // jets - 2lep
    if (nJets != 0 && nElectrons >= 2) {
        m_hm->fill("nJets_2lep", nJets);
        for (auto j : goodJets) {
            float mj = j->m() / 1000.;
            float ptj = j->pt() / 1000.;
            float emfrac;
            std::vector<int> ntracks;
            j->getAttribute(xAOD::JetAttribute::EMFrac, emfrac);
            j->getAttribute(xAOD::JetAttribute::NumTrkPt1000, ntracks);
            float mindR = 100000.;
            for (auto el : goodElectrons) {
                float curdR = j->p4().DeltaR(el->p4());
                if (curdR < mindR) mindR = curdR;
            }
            m_hm->fill("mj_2lep", mj);
            m_hm->fill("ptj_2lep", ptj);
            m_hm->fill("ntracksj_2lep", ntracks[0]);
            m_hm->fill("emfracj_2lep", emfrac);
            m_hm->fill("dRej_2lep", mindR);
        }
        
        int nJets_emfrac = 0;
        for (auto j : goodJets) {
            j->auxdata<int>("emfraccut") = 0;
            float emfrac;
            j->getAttribute(xAOD::JetAttribute::EMFrac, emfrac);
            if (emfrac > 0.97) {
                j->auxdata<int>("emfraccut") = 1;
                nJets_emfrac++;
            }
        }

        if (nJets_emfrac != 0) {
            m_hm->fill("nJets_postemfrac_2lep", nJets_emfrac);
            for (auto j : goodJets) {
                if (!(j->auxdata<int>("emfraccut"))) continue;
                float mj = j->m() / 1000.;
                float ptj = j->pt() / 1000.;
                std::vector<int> ntracks;
                j->getAttribute(xAOD::JetAttribute::NumTrkPt1000, ntracks);
                float mindR = 100000.;
                for (auto *el : goodElectrons) {
                    float curdR = j->p4().DeltaR(el->p4());
                    if (curdR < mindR) mindR = curdR;
                }
                m_hm->fill("mj_postemfrac_2lep", mj);
                m_hm->fill("ptj_postemfrac_2lep", ptj);
                m_hm->fill("ntracksj_postemfrac_2lep", ntracks[0]);
                m_hm->fill("dRej_postemfrac_2lep", mindR);
            }
        }
    }
    
    // jets - 1lep
    if (nJets != 0 && nElectrons == 1) {
        m_hm->fill("nJets_1lep", nJets);
        for (auto j : goodJets) {
            float mj = j->m() / 1000.;
            float ptj = j->pt() / 1000.;
            float emfrac;
            std::vector<int> ntracks;
            j->getAttribute(xAOD::JetAttribute::EMFrac, emfrac);
            j->getAttribute(xAOD::JetAttribute::NumTrkPt1000, ntracks);
            float mindR = j->p4().DeltaR(goodElectrons.at(0)->p4());
            m_hm->fill("mj_1lep", mj);
            m_hm->fill("ptj_1lep", ptj);
            m_hm->fill("ntracksj_1lep", ntracks[0]);
            m_hm->fill("emfracj_1lep", emfrac);
            m_hm->fill("dRej_1lep", mindR);
        }
        
        int nJets_emfrac = 0;
        for (auto j : goodJets) {
            j->auxdata<int>("emfraccut") = 0;
            float emfrac;
            j->getAttribute(xAOD::JetAttribute::EMFrac, emfrac);
            if (emfrac > 0.97) {
                j->auxdata<int>("emfraccut") = 1;
                nJets_emfrac++;
            }
        }

        if (nJets_emfrac != 0) {
            m_hm->fill("nJets_postemfrac_1lep", nJets_emfrac);
            for (auto j : goodJets) {
                if (!(j->auxdata<int>("emfraccut"))) continue;
                float mj = j->m() / 1000.;
                float ptj = j->pt() / 1000.;
                std::vector<int> ntracks;
                j->getAttribute(xAOD::JetAttribute::NumTrkPt1000, ntracks);
                float mindR = j->p4().DeltaR(goodElectrons.at(0)->p4());
                    m_hm->fill("mj_postemfrac_1lep", mj);
                m_hm->fill("ptj_postemfrac_1lep", ptj);
                m_hm->fill("ntracksj_postemfrac_1lep", ntracks[0]);
                m_hm->fill("dRej_postemfrac_1lep", mindR);
            }
        }
    }

    // jets - 0lep
    if (nJets != 0 && nElectrons == 0) {
        m_hm->fill("nJets_0lep", nJets);
        for (auto j : goodJets) {
            float mj = j->m() / 1000.;
            float ptj = j->pt() / 1000.;
            float emfrac;
            std::vector<int> ntracks;
            j->getAttribute(xAOD::JetAttribute::EMFrac, emfrac);
            j->getAttribute(xAOD::JetAttribute::NumTrkPt1000, ntracks);
            m_hm->fill("mj_0lep", mj);
            m_hm->fill("ptj_0lep", ptj);
            m_hm->fill("ntracksj_0lep", ntracks[0]);
            m_hm->fill("emfracj_0lep", emfrac);
        }
        
        int nJets_emfrac = 0;
        for (auto j : goodJets) {
            j->auxdata<int>("emfraccut") = 0;
            float emfrac;
            j->getAttribute(xAOD::JetAttribute::EMFrac, emfrac);
            if (emfrac > 0.97) {
                j->auxdata<int>("emfraccut") = 1;
                nJets_emfrac++;
            }
        }

        if (nJets_emfrac != 0) {
            m_hm->fill("nJets_postemfrac_0lep", nJets_emfrac);
            for (auto j : goodJets) {
                if (!(j->auxdata<int>("emfraccut"))) continue;
                float mj = j->m() / 1000.;
                float ptj = j->pt() / 1000.;
                std::vector<int> ntracks;
                j->getAttribute(xAOD::JetAttribute::NumTrkPt1000, ntracks);
                m_hm->fill("mj_postemfrac_0lep", mj);
                m_hm->fill("ptj_postemfrac_0lep", ptj);
                m_hm->fill("ntracksj_postemfrac_0lep", ntracks[0]);
            }
        }
    }
    */

    // just in case
    return EL::StatusCode::SUCCESS;
}

void FatLeptons::prepareHistograms() {
    // LASER: CaloCluster study
    /*
    m_hm->book(new TH1F("nMatchedTopo2", "nMatchedTopo2", 20, -0.5, 19.5), "nMatchedTopo2");
    m_hm->book(new TH1F("mcc_2", "mcc_2", 200, 0, 200), "mcc_2");
    m_hm->book(new TH1F("mcc_3", "mcc_3", 200, 0, 200), "mcc_3");
    m_hm->book(new TH1F("mcc_4", "mcc_4", 200, 0, 200), "mcc_4");
    m_hm->book(new TH1F("mcc_5", "mcc_5", 200, 0, 200), "mcc_5");
    */
    
    // LASER: Christos' shit
    /*
    m_hm->book(new TH1F("nTracks1TruthElectron", "nTracks1TruthElectron", 6, -0.5, 5.5), "Christos_nTracks1TruthElectron");
    m_hm->book(new TH1F("nMatchedTopo", "nMatchedTopo", 20, -0.5, 19.5), "Christos_nMatchedTopo");
    */

    // event weights
    m_hm->book(new TH1F("sumOfWeights", "sumOfWeights", 1, -0.5, 0.5), "sumOfWeights");

    // pre event selection
    m_hm->book(new TH1F("nElectrons", "nElectrons", 6, -0.5, 5.5), "nElectrons");
    m_hm->book(new TH1F("nJets", "nJets", 6, -0.5, 5.5), "nJets");
    m_hm->book(new TH1F("mj", "mj", 200, 0, 200), "mj");
    m_hm->book(new TH1F("ptj", "ptj", 200, 0, 2000), "ptj");
    m_hm->book(new TH1F("ntracksj", "ntracksj", 50, -0.5, 49.5), "ntracksj");
    m_hm->book(new TH1F("emfracj", "emfracj", 100, 0, 1), "emfracj");
    m_hm->book(new TH1F("ntopocj", "ntopocj", 30, -0.5, 29.5), "ntopocj");

    /*      
    // 2lep plots - electrons
    m_hm->book(new TH1F("mee", "mee", 200, 0, 200), "mee");
    m_hm->book(new TH1F("ptee", "ptee", 1000, 0, 1000), "ptee");
    */

    // inclusive plots - jets
    // postemfrac
    m_hm->book(new TH1F("mj_postemfrac", "mj_postemfrac", 200, 0, 200), "mj_postemfrac");
    m_hm->book(new TH1F("ptj_postemfrac", "ptj_postemfrac", 200, 0, 2000), "ptj_postemfrac");
    m_hm->book(new TH1F("ntracksj_postemfrac", "ntracksj_postemfrac", 50, -0.5, 49.5), "ntracksj_postemfrac");
    m_hm->book(new TH1F("ntopocj_postemfrac", "ntopocj_postemfrac", 30, -0.5, 29.5), "ntopocj_postemfrac");
    // postntracks
    m_hm->book(new TH1F("mj_postntracks", "mj_postntracks", 200, 0, 200), "mj_postntracks");
    m_hm->book(new TH1F("ptj_postntracks", "ptj_postntracks", 200, 0, 2000), "ptj_postntracks");
    m_hm->book(new TH1F("emfracj_postntracks", "emfracj_postntracks", 100, 0, 1), "emfracj_postntracks");
    m_hm->book(new TH1F("ntopocj_postntracks", "ntopocj_postntracks", 30, -0.5, 29.5), "ntopocj_postntracks");
    // postemfracntracks
    m_hm->book(new TH1F("mj_postemfracntracks", "mj_postemfracntracks", 200, 0, 200), "mj_postemfracntracks");
    m_hm->book(new TH1F("ptj_postemfracntracks", "ptj_postemfracntracks", 200, 0, 2000), "ptj_postemfracntracks");
    m_hm->book(new TH1F("ntopocj_postemfracntracks", "ntopocj_postemfracntracks", 30, -0.5, 29.5), "ntopocj_postemfracntracks");
    // failemfrac
    m_hm->book(new TH1F("mj_failemfrac", "mj_failemfrac", 200, 0, 200), "mj_failemfrac");
    m_hm->book(new TH1F("ptj_failemfrac", "ptj_failemfrac", 200, 0, 2000), "ptj_failemfrac");
    m_hm->book(new TH1F("ntracksj_failemfrac", "ntracksj_failemfrac", 50, -0.5, 49.5), "ntracksj_failemfrac");
    m_hm->book(new TH1F("ntopocj_failemfrac", "ntopocj_failemfrac", 30, -0.5, 29.5), "ntopocj_failemfrac");
    // failntracks
    m_hm->book(new TH1F("mj_failntracks", "mj_failntracks", 200, 0, 200), "mj_failntracks");
    m_hm->book(new TH1F("ptj_failntracks", "ptj_failntracks", 200, 0, 2000), "ptj_failntracks");
    m_hm->book(new TH1F("emfracj_failntracks", "emfracj_failntracks", 100, 0, 1), "emfracj_failntracks");
    m_hm->book(new TH1F("ntopocj_failntracks", "ntopocj_failntracks", 30, -0.5, 29.5), "ntopocj_failntracks");
    // failemfracntracks
    m_hm->book(new TH1F("mj_failemfracntracks", "mj_failemfracntracks", 200, 0, 200), "mj_failemfracntracks");
    m_hm->book(new TH1F("ptj_failemfracntracks", "ptj_failemfracntracks", 200, 0, 2000), "ptj_failemfracntracks");
    m_hm->book(new TH1F("ntopocj_failemfracntracks", "ntopocj_failemfracntracks", 30, -0.5, 29.5), "ntopocj_failemfracntracks");

    // LASER - matching study
    m_hm->book(new TH1F("matching_nMatchedElectrons", "matching_nMatchedElectrons", 4, -0.5, 3.5), "matching_nMatchedElectrons");
    m_hm->book(new TH1F("matching_mj", "matching_mj", 200, 0, 200), "matching_mj");
    m_hm->book(new TH1F("matching_nMCC", "matching_nMCC", 8, -0.5, 7.5), "matching_nMCC");
    m_hm->book(new TH1F("matching_mcc2", "matching_mcc2", 200, 0, 200), "matching_mcc2");
    m_hm->book(new TH1F("matching_mcc3", "matching_mcc3", 200, 0, 200), "matching_mcc3");
    m_hm->book(new TH1F("matching_mcc4", "matching_mcc4", 200, 0, 200), "matching_mcc4");
    m_hm->book(new TH1F("matching_mcc5", "matching_mcc5", 200, 0, 200), "matching_mcc5");

    // LASER - reco study
    m_hm->book(new TH1F("reco_mcc2", "reco_mcc2", 200, 0, 200), "reco_mcc2");
    m_hm->book(new TH1F("reco_mcc3", "reco_mcc3", 200, 0, 200), "reco_mcc3");
    m_hm->book(new TH1F("reco_mcc4", "reco_mcc4", 200, 0, 200), "reco_mcc4");
    m_hm->book(new TH1F("reco_mcc5", "reco_mcc5", 200, 0, 200), "reco_mcc5");

    /*      
    // 2lep plots - jets
    m_hm->book(new TH1F("nJets_2lep", "nJets_2lep", 6, -0.5, 5.5), "nJets_2lep");
    m_hm->book(new TH1F("mj_2lep", "mj_2lep", 200, 0, 200), "mj_2lep");
    m_hm->book(new TH1F("ptj_2lep", "ptj_2lep", 1000, 0, 1000), "ptj_2lep");
    m_hm->book(new TH1F("ntracksj_2lep", "ntracksj_2lep", 50, 0, 50), "ntracksj_2lep");
    m_hm->book(new TH1F("emfracj_2lep", "emfracj_2lep", 100, 0, 1), "emfracj_2lep");
    m_hm->book(new TH1F("dRej_2lep", "dRej_2lep", 40, 0, 4), "dRej_2lep");
    m_hm->book(new TH1F("nJets_postemfrac_2lep", "nJets_postemfrac_2lep", 6, -0.5, 5.5), "nJets_postemfrac_2lep");
    m_hm->book(new TH1F("mj_postemfrac_2lep", "mj_postemfrac_2lep", 200, 0, 200), "mj_postemfrac_2lep");
    m_hm->book(new TH1F("ptj_postemfrac_2lep", "ptj_postemfrac_2lep", 1000, 0, 1000), "ptj_postemfrac_2lep");
    m_hm->book(new TH1F("ntracksj_postemfrac_2lep", "ntracksj_postemfrac_2lep", 50, -0.5, 49.5), "ntracksj_postemfrac_2lep");
    m_hm->book(new TH1F("dRej_postemfrac_2lep", "dRej_postemfrac_2lep", 40, 0, 4), "dRej_postemfrac_2lep");

    // 1lep plots - jets
    m_hm->book(new TH1F("nJets_1lep", "nJets_1lep", 6, -0.5, 5.5), "nJets_1lep");
    m_hm->book(new TH1F("mj_1lep", "mj_1lep", 200, 0, 200), "mj_1lep");
    m_hm->book(new TH1F("ptj_1lep", "ptj_1lep", 1000, 0, 1000), "ptj_1lep");
    m_hm->book(new TH1F("ntracksj_1lep", "ntracksj_1lep", 50, 0, 50), "ntracksj_1lep");
    m_hm->book(new TH1F("emfracj_1lep", "emfracj_1lep", 100, 0, 1), "emfracj_1lep");
    m_hm->book(new TH1F("dRej_1lep", "dRej_1lep", 40, 0, 4), "dRej_1lep");
    m_hm->book(new TH1F("nJets_postemfrac_1lep", "nJets_postemfrac_1lep", 6, -0.5, 5.5), "nJets_postemfrac_1lep");
    m_hm->book(new TH1F("mj_postemfrac_1lep", "mj_postemfrac_1lep", 200, 0, 200), "mj_postemfrac_1lep");
    m_hm->book(new TH1F("ptj_postemfrac_1lep", "ptj_postemfrac_1lep", 1000, 0, 1000), "ptj_postemfrac_1lep");
    m_hm->book(new TH1F("ntracksj_postemfrac_1lep", "ntracksj_postemfrac_1lep", 50, -0.5, 49.5), "ntracksj_postemfrac_1lep");
    m_hm->book(new TH1F("dRej_postemfrac_1lep", "dRej_postemfrac_1lep", 40, 0, 4), "dRej_postemfrac_1lep");

    // 0lep plots - jets
    m_hm->book(new TH1F("nJets_0lep", "nJets_0lep", 6, -0.5, 5.5), "nJets_0lep");
    m_hm->book(new TH1F("mj_0lep", "mj_0lep", 200, 0, 200), "mj_0lep");
    m_hm->book(new TH1F("ptj_0lep", "ptj_0lep", 1000, 0, 1000), "ptj_0lep");
    m_hm->book(new TH1F("ntracksj_0lep", "ntracksj_0lep", 50, 0, 50), "ntracksj_0lep");
    m_hm->book(new TH1F("emfracj_0lep", "emfracj_0lep", 100, 0, 1), "emfracj_0lep");
    m_hm->book(new TH1F("nJets_postemfrac_0lep", "nJets_postemfrac_0lep", 6, -0.5, 5.5), "nJets_postemfrac_0lep");
    m_hm->book(new TH1F("mj_postemfrac_0lep", "mj_postemfrac_0lep", 200, 0, 200), "mj_postemfrac_0lep");
    m_hm->book(new TH1F("ptj_postemfrac_0lep", "ptj_postemfrac_0lep", 1000, 0, 1000), "ptj_postemfrac_0lep");
    m_hm->book(new TH1F("ntracksj_postemfrac_0lep", "ntracksj_postemfrac_0lep", 50, -0.5, 49.5), "ntracksj_postemfrac_0lep");
    */
}

const xAOD::TruthParticle * FatLeptons::getTruthZ(xAOD::Jet *j) {
    std::vector<const xAOD::TruthParticle*> tps;
    j->getAssociatedObjects("GhostZBosons", tps);
    if (tps.size() != 0) return tps[0];
    else return nullptr;
}

int FatLeptons::getNumberOfTruthMatchedElectrons(const xAOD::TruthParticleContainer *tps, const xAOD::ElectronContainer *electrons, std::vector<const xAOD::Electron*> &goodElectrons) {
    const xAOD::TruthParticle *tel1 = 0, *tel2 = 0;
    for (const xAOD::TruthParticle *tp : *tps) {
        if (tp->pdgId() == 23 && tp->status() == 22) {
            if (tp->nChildren() == 2) {
                if (abs(tp->child(0)->pdgId()) != 11 || abs(tp->child(1)->pdgId()) != 11) continue;
                tel1 = tp->child(0);
                tel2 = tp->child(1);
                break;
            }
        }
    }

    if (!tel1 || !tel2) return -1;

    int ntruth = 0;
    for (const xAOD::Electron *el : *electrons) {
        el->auxdecor<int>("matched") = 0;
        const xAOD::TruthParticle *tp = xAOD::TruthHelpers::getTruthParticle(*el);
        if (tp) {
            if (abs(tp->pdgId()) == 11) {
                const xAOD::TruthParticle *cur = tp;
                while (cur->nParents() > 0) {
                    if (cur == tel1 || cur == tel2) {
                        el->auxdecor<int>("matched") = 1;
                        goodElectrons.push_back(el);
                        ntruth++;
                        break;
                    }
                    cur = cur->parent();
                }
            }
        }
    }

    return ntruth;
}

int FatLeptons::getNMatched(const xAOD::Electron *el, const xAOD::CaloClusterContainer *cls) {
    int nmatchc = 0;
    for (const xAOD::CaloCluster *cl : *cls) {
        if(cl == el->caloCluster()) continue;
        if (cl->p4().DeltaR(el->caloCluster()->p4()) < 0.1) nmatchc++;
    }
    return nmatchc;
}
        
std::vector<const xAOD::CaloCluster*> FatLeptons::getMatchedCaloClusters(const xAOD::Electron *el, const xAOD::CaloClusterContainer *cls) {
    caloClusterComp comparator;
    comparator.refCC = el->caloCluster();
    std::vector<const xAOD::CaloCluster*> ccs;
    for (const xAOD::CaloCluster *cl : *cls) {
        if (cl == el->caloCluster()) continue;
        if (cl->p4().DeltaR(el->caloCluster()->p4()) < 0.2) ccs.push_back(cl);
    }
    std::sort(ccs.begin(), ccs.end(), comparator);
    return ccs;
}

int FatLeptons::getNMatchedEls(const xAOD::Jet *jet, const xAOD::ElectronContainer *els) {
    // get matched electrons dR < 0.2
    int nmatchedels = 0;
    for (const xAOD::Electron *el : *els) {
        TLorentzVector jtlv;
        jtlv.SetPtEtaPhiM(jet->jetP4("JetEMScaleMomentum").Pt(), jet->jetP4("JetEMScaleMomentum").Eta(), jet->jetP4("JetEMScaleMomentum").Phi(), jet->jetP4("JetEMScaleMomentum").M());
        if (el->p4().DeltaR(jtlv) < 0.4) nmatchedels++;
    }

    return nmatchedels;
}

void FatLeptons::doMatchingStudy(const xAOD::Jet *jet, const xAOD::ElectronContainer *els, const xAOD::CaloClusterContainer *cls) {
    // get matched electrons dR < 0.2
    std::vector<const xAOD::Electron*> matchedEls;
    for (const xAOD::Electron *el : *els) {
        TLorentzVector jtlv;
        jtlv.SetPtEtaPhiM(jet->jetP4("JetEMScaleMomentum").Pt(), jet->jetP4("JetEMScaleMomentum").Eta(), jet->jetP4("JetEMScaleMomentum").Phi(), jet->jetP4("JetEMScaleMomentum").M());
        if (el->p4().DeltaR(jtlv) < 0.4) {
            matchedEls.push_back(el);
        }
    }

    m_hm->fill("matching_nMatchedElectrons", matchedEls.size());

    // if there are none, return
    if (matchedEls.size() == 0) return;

    m_hm->fill("matching_mj", jet->jetP4("JetEMScaleMomentum").M() / 1000.);

    // sort by pt
    std::sort(matchedEls.begin(), matchedEls.end(), compPt);

    // and take the highest one
    const xAOD::Electron *el = matchedEls[0];

    // get matched calo clusters to electron
    caloClusterComp comparator;
    comparator.refCC = el->caloCluster();
    std::vector<const xAOD::CaloCluster*> ccs;
    for (const xAOD::CaloCluster *cl : *cls) {
        if (cl == el->caloCluster()) continue;
        if (cl->p4().DeltaR(el->caloCluster()->p4()) < 0.2) ccs.push_back(cl);
    }
    std::sort(ccs.begin(), ccs.end(), comparator);

    m_hm->fill("matching_nMCC", ccs.size());

    // if there are none, return
    if (ccs.size() == 0) return;

    // add and plot
    TLorentzVector vec = ccs[0]->p4();
    if (ccs.size() > 1) {vec += ccs[1]->p4(); m_hm->fill("matching_mcc2", vec.M() / 1000.);}
    if (ccs.size() > 2) {vec += ccs[2]->p4(); m_hm->fill("matching_mcc3", vec.M() / 1000.);}
    if (ccs.size() > 3) {vec += ccs[3]->p4(); m_hm->fill("matching_mcc4", vec.M() / 1000.);}
    if (ccs.size() > 4) {vec += ccs[4]->p4(); m_hm->fill("matching_mcc5", vec.M() / 1000.);}
}

void FatLeptons::doRecoStudy(const xAOD::Jet *jet) {
    // get consituents
    std::vector<const xAOD::IParticle*> ccs = jet->getConstituents().asIParticleVector();

    // sort by pt
    std::sort(ccs.begin(), ccs.end(), compPt);

    // add and plot
    TLorentzVector vec = ccs[0]->p4();
    if (ccs.size() > 1) {vec += ccs[1]->p4(); m_hm->fill("reco_mcc2", vec.M() / 1000.);}
    if (ccs.size() > 2) {vec += ccs[2]->p4(); m_hm->fill("reco_mcc3", vec.M() / 1000.);}
    if (ccs.size() > 3) {vec += ccs[3]->p4(); m_hm->fill("reco_mcc4", vec.M() / 1000.);}
    if (ccs.size() > 4) {vec += ccs[4]->p4(); m_hm->fill("reco_mcc5", vec.M() / 1000.);}
}
