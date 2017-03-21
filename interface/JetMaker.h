#ifndef ANALYSIS_VLQANA_JETMAKER_H
#define ANALYSIS_VLQANA_JETMAKER_H

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "Analysis/VLQAna/interface/JetID.h"
#include "AnalysisDataFormats/BoostedObjects/interface/Jet.h"

#include <boost/algorithm/string.hpp>
#include <string>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TF1.h>

using namespace std;
using namespace edm ; 

class JetMaker { 
  public:
    enum JETTYPES_t { AK4JET, AK8JET, N_JETTYPES} ; 
    JetMaker () ; 
    ~JetMaker () ; 
    JetMaker (edm::ParameterSet const&, edm::ConsumesCollector && iC) ;
    void operator()(edm::Event& evt, vlq::JetCollection& jetOut) ;
    void operator()(vlq::JetCollection jetsIn, vlq::JetCollection& jetsOut) ;

    double idxjetPtMin_          ;
    double idxjetPtMax_          ;
    double idxjetAbsEtaMax_      ;
    double idxjetCSVDiscMin_     ; 
    double idxjetCSVDiscMax_     ; 
    double idxjettau2Bytau1Min_  ; 
    double idxjettau2Bytau1Max_  ; 
    double idxjettau3Bytau2Min_  ; 
    double idxjettau3Bytau2Max_  ; 
    double idxjetMassMin_        ; 
    double idxjetMassMax_        ; 
    double idxjetPrunedMassMin_  ; 
    double idxjetPrunedMassMax_  ; 
    double idxjetTrimmedMassMin_ ; 
    double idxjetTrimmedMassMax_ ; 
    double idxjetSoftDropMassMin_; 
    double idxjetSoftDropMassMax_; 
    double idxjetnSubJetsMin_    ; 
    double idxjetnSubJetsMax_    ; 
    double idxsjMassMin_         ; 
    double idxsjCSVMin_          ; 
    double idxsjCSVMax_          ;
    double idxjetDoubleBMin_     ; //EJS
    double idxjetDoubleBMax_     ; //EJS
    double idxsjHighestCSVMin_   ; 
    double idxsjHighestCSVMax_   ; 
  private:
    edm::ParameterSet JetIDParams_ ; 
    edm::ParameterSet JetSubstrParams_ ; 
    edm::ParameterSet SubjetParams_ ; 

    bool         scaleJetP4_     ;
    double       scaledJetMass_  ; 

    double btaggedcsvlOP_        ; 
    double btaggedcsvmOP_        ; 
    double btaggedcsvtOP_        ; 

    edm::EDGetTokenT<int>                t_npv                ; 
    edm::EDGetTokenT<double>             t_rho                ; 
    edm::EDGetTokenT<std::vector<float>> t_jetPt              ; 
    edm::EDGetTokenT<std::vector<float>> t_jetEta             ; 
    edm::EDGetTokenT<std::vector<float>> t_jetPhi             ; 
    edm::EDGetTokenT<std::vector<float>> t_jetMass            ; 
    edm::EDGetTokenT<std::vector<float>> t_jetEnergy          ; 
    edm::EDGetTokenT<std::vector<float>> t_jetPartonFlavour   ; 
    edm::EDGetTokenT<std::vector<float>> t_jetHadronFlavour   ; 
    edm::EDGetTokenT<std::vector<float>> t_jetCSV             ; 
    edm::EDGetTokenT<std::vector<float>> t_jetCMVA            ; 
    edm::EDGetTokenT<std::vector<float>> t_jetCvsB            ; 
    edm::EDGetTokenT<std::vector<float>> t_jetCvsL            ; 
    edm::EDGetTokenT<std::vector<float>> t_jetDoubleBAK8      ; 
    edm::EDGetTokenT<std::vector<float>> t_jetDoubleBCA15     ; 
    edm::EDGetTokenT<std::vector<float>> t_jetJEC             ; 
    edm::EDGetTokenT<std::vector<float>> t_jetnHadEnergy      ;
    edm::EDGetTokenT<std::vector<float>> t_jetnEMEnergy       ;
    edm::EDGetTokenT<std::vector<float>> t_jetHFHadronEnergy  ;
    edm::EDGetTokenT<std::vector<float>> t_jetcHadEnergy      ;
    edm::EDGetTokenT<std::vector<float>> t_jetcEMEnergy       ;
    edm::EDGetTokenT<std::vector<float>> t_jetnumDaughters    ;
    edm::EDGetTokenT<std::vector<float>> t_jetcMultip         ;
    edm::EDGetTokenT<std::vector<float>> t_jetnMultip         ;
    edm::EDGetTokenT<std::vector<float>> t_jetY               ;
    edm::EDGetTokenT<std::vector<float>> t_jetArea            ; 
    edm::EDGetTokenT<std::vector<float>> t_jetGenJetPt        ; 
    edm::EDGetTokenT<std::vector<float>> t_jetGenJetEta       ; 
    edm::EDGetTokenT<std::vector<float>> t_jetGenJetPhi       ; 
    edm::EDGetTokenT<std::vector<float>> t_jetGenJetE         ; 
    edm::EDGetTokenT<std::vector<float>> t_jetGenJetCharge    ; 
    edm::EDGetTokenT<std::vector<float>> t_jettau1            ; 
    edm::EDGetTokenT<std::vector<float>> t_jettau2            ; 
    edm::EDGetTokenT<std::vector<float>> t_jettau3            ; 
    edm::EDGetTokenT<std::vector<float>> t_jetPrunedMass      ; 
    edm::EDGetTokenT<std::vector<float>> t_jetTrimmedMass     ; 
    edm::EDGetTokenT<std::vector<float>> t_jetFilteredMass    ; 
    edm::EDGetTokenT<std::vector<float>> t_jetSoftDropMass    ; 
    edm::EDGetTokenT<std::vector<float>> t_jetnSubJets        ; 
    edm::EDGetTokenT<std::vector<float>> t_vjetsjIdx0         ; 
    edm::EDGetTokenT<std::vector<float>> t_vjetsjIdx1         ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjHadFl         ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjPt            ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjEta           ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjPhi           ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjMass          ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjEnergy        ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjJEC           ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjCSV           ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjCMVA          ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjGenJetPt      ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjGenJetEta     ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjGenJetPhi     ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjGenJetE       ; 
    edm::EDGetTokenT<std::vector<float>> t_ak8sjGenJetCharge  ; 
    double                    jecShift_ ; 
    int                       jerShift_ ; 
    std::vector<std::string>  newJECPayloadNames_;
    std::string               jecUncPayloadName_;
    std::vector<std::string>  jecAK8GroomedPayloadNames_;
    bool                      doGroomedMassCorr_;
    std::string               SDMassCorrWt_   ;
    bool                      doNewJEC_;
    boost::shared_ptr<FactorizedJetCorrector>   ptr_jecAK8_  ;
    boost::shared_ptr<JetCorrectionUncertainty> ptr_jecUnc ; 
    boost::shared_ptr<FactorizedJetCorrector>   ptr_newJEC_  ;
    std::string jecAK8UncName_  ;
    std::unique_ptr<TFile> puppiSDMassCorrWtFile_ ; 
    std::unique_ptr<TF1> puppiSDMassCorrWtFunGen_ ; 
    std::unique_ptr<TF1> puppiSDMassCorrWtFunRecoCen_ ; 
    std::unique_ptr<TF1> puppiSDMassCorrWtFunRecoFor_ ; 

    JETTYPES_t type_ ; 
    JetID jetID_ ; 
};

#endif 
