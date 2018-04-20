// -*- C++ -*-
//
// Package:    Analysis/VLQAna
// Class:      OS2LAna
// 
/**\class VLQAna OS2LAna.cc Analysis/VLQAna/plugins/OS2LAna.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Devdatta Majumder
//         Created:  Fri, 27 Feb 2015 16:09:10 GMT
// Modified: Sadia Khalil
//           25 Mar 2016 17:11 CDT
//

#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <math.h>
#include <unordered_set>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "AnalysisDataFormats/BoostedObjects/interface/GenParticleWithDaughters.h"
#include "AnalysisDataFormats/BoostedObjects/interface/ResolvedVjj.h"

#include "Analysis/VLQAna/interface/Utilities.h"
#include "Analysis/VLQAna/interface/DileptonCandsProducer.h"
#include "Analysis/VLQAna/interface/CandidateFilter.h"
#include "Analysis/VLQAna/interface/MuonMaker.h"
#include "Analysis/VLQAna/interface/ElectronMaker.h"
#include "Analysis/VLQAna/interface/JetMaker.h"
#include "Analysis/VLQAna/interface/HT.h"
#include "Analysis/VLQAna/interface/ApplyLeptonIDSFs.h"
#include "Analysis/VLQAna/interface/CandidateCleaner.h"
#include "Analysis/VLQAna/interface/METMaker.h"
#include "Analysis/VLQAna/interface/PickGenPart.h"
#include "Analysis/VLQAna/interface/JetID.h"
#include "Analysis/VLQAna/interface/BTagSFUtils.h"
#include "Analysis/VLQAna/interface/ApplyLeptonTrigSFs.h"
#include "Analysis/VLQAna/interface/DYNLOEwkKfact.h"
#include "Analysis/VLQAna/interface/OS2LTree.h"
#include "Analysis/VLQAna/interface/HCandsProducer.h"
#include "Analysis/VLQAna/interface/ZCandsProducer.h"
#include "Analysis/VLQAna/interface/TopCandsProducer.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>

//
// class declaration
//

class OS2LAna : public edm::EDFilter {
  public:
    explicit OS2LAna(const edm::ParameterSet&);
    ~OS2LAna();

  private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    void fillAdditionalPlots( vlq::ElectronCollection goodElectrons,double evtwt);
    pair<double, double> vector_eval(vector<pair<double, double> >);
    double resolvedChi2(vector<TLorentzVector>, TLorentzVector, double, double);
    double boostedChi2(vector<TLorentzVector>, TLorentzVector, TLorentzVector, double, double, double);
    pair<double, double> doBoostedReco(vlq::JetCollection&, TLorentzVector, double, TLorentzVector, double);
    pair<double, double> doResolvedReco(vlq::JetCollection&, double, TLorentzVector);
    void fillPdfHistos(string name, double value, double evtwt, vector<pair<int, double>> &lhe); 
    // ----------member data ---------------------------
    edm::EDGetTokenT<string>   t_evttype         ;
    edm::EDGetTokenT<double>   t_evtwtGen        ;
    edm::EDGetTokenT<double>   t_evtwtPV         ;
    edm::EDGetTokenT<double>   t_evtwtPVLow      ;
    edm::EDGetTokenT<double>   t_evtwtPVHigh     ;
    edm::EDGetTokenT<unsigned> t_npv             ;
    edm::EDGetTokenT<bool>     t_hltdecision     ;
    edm::EDGetTokenT<bool>     t_hltreject       ;
    edm::EDGetTokenT<int>      t_evtno           ; 
    edm::EDGetTokenT<int>      t_runno           ;
    edm::EDGetTokenT<int>      t_lumisec         ;
    edm::EDGetTokenT<vector<int> > t_lhewtids    ;
    edm::EDGetTokenT<vector<double> > t_lhewts   ;
    edm::ParameterSet DilepCandParams_           ; 
    edm::ParameterSet ZCandParams_               ; 
    //edm::ParameterSet BoostedZCandParams_        ; 
    //edm::ParameterSet GenHSelParams_             ;
    edm::ParameterSet genParams_                 ;
    const unsigned int NAK4Min_                  ;
    const double HTMin_                          ;
    const double STMin_                          ; 
    const double STMaxControl_                   ;
    const bool skim_                             ;
    const bool isData_                           ;
    const bool run_photons_                      ;
    const bool filterSignal_                     ;
    const bool additionalPlots_                  ;
    const std::string signalType_                ;
    const std::string zdecayMode_                ;
    const bool btagsf_bcUp_                      ;
    const bool btagsf_bcDown_                    ;
    const bool btagsf_lUp_                       ;
    const bool btagsf_lDown_                     ;
    const bool sbtagsf_bcUp_                     ;
    const bool sbtagsf_bcDown_                   ;
    const bool sbtagsf_lUp_                      ;
    const bool sbtagsf_lDown_                    ;
    const bool PileupUp_                         ;
    const bool PileupDown_                       ;   
    const int lheId_                             ;
    const bool applyLeptonIDSFs_                 ;
    const bool applyLeptonTrigSFs_               ;
    const bool applyBTagSFs_                     ;
    const bool applyDYNLOCorr_                   ;
    const bool DYSFSyst_                         ;
    const int  tau21Shift_                       ;
    const int  tau32Shift_                       ;
    const int  WZShift_                          ;
    const int  higgsShift_                       ;
    const int  topShift_                         ;
    const int pdfID_offset_                      ;
    const int scale_offset_                      ;
    const bool syst_                             ;
    const bool vv_                               ;
    const int elSyst_                            ;  
    const std::string fname_DYNLOCorr_           ; 
    const std::string funname_DYNLOCorr_         ; 
    DYNLOEwkKfact dynloewkkfact                  ;
    ApplyLeptonIDSFs lepIdSFs                    ;
    ApplyLeptonTrigSFs lepTrigSFs                ;
    METMaker metmaker                            ;
    MuonMaker muonmaker                          ; 
    ElectronMaker electronmaker                  ; 
    JetMaker jetAK4maker                         ; 
    JetMaker jetAK4BTaggedmaker                  ; 
    JetMaker jetAK8maker                         ; 
    JetMaker jetHTaggedmaker                     ; 
    JetMaker jetWTaggedmaker                     ; 
    JetMaker jetTopTaggedmaker                   ; 
    edm::Service<TFileService> fs                ; 
    std::map<std::string, TH1D*> h1_             ; 
    std::map<std::string, TH2D*> h2_             ; 
    PickGenPart genpart                          ;
    const std::string fnamebtagSF_               ;
    const std::string fnameSJbtagSF_             ;
    const std::string btageffmap_                ;
    std::unique_ptr<BTagSFUtils> btagsfutils_    ; 
    std::unique_ptr<BTagSFUtils> sjbtagsfutils_  ; 

    const bool maketree_ ;
    TTree* tree_ ; 
    os2l::OS2LAnaTree os2ltree_ ; 
};

using namespace std;

// static data member definitions
void OS2LAna::fillAdditionalPlots( vlq::ElectronCollection goodElectrons,double evtwt){

  for  (unsigned int iele=0; iele<goodElectrons.size(); ++iele){
    float scEta = goodElectrons.at(iele).getscEta();
    if(fabs(scEta) <= 1.479){
      h1_["Eta_EB_el_pre"]-> Fill(goodElectrons.at(iele).getEta(), evtwt);
      h1_["Iso03_EB_el_pre"]->Fill(goodElectrons.at(iele).getIso03(), evtwt);
      h1_["dEtaInSeed_EB_el_pre"]->Fill(goodElectrons.at(iele).getdEtaInSeed(), evtwt);
      h1_["dPhiIn_EB_el_pre"]->Fill(goodElectrons.at(iele).getdPhiIn(), evtwt);
      h1_["Dz_EB_el_pre"]->Fill(goodElectrons.at(iele).getDz(), evtwt);
      h1_["Dxy_EB_el_pre"]->Fill(goodElectrons.at(iele).getDxy(), evtwt);
    }
    else if  (fabs(scEta) > 1.479 && fabs(scEta) < 2.4){
      h1_["Eta_EE_el_pre"]->Fill(goodElectrons.at(iele).getEta(), evtwt);
      h1_["Iso03_EE_el_pre"]->Fill(goodElectrons.at(iele).getIso03(), evtwt);
      h1_["dEtaInSeed_EE_el_pre"]->Fill(goodElectrons.at(iele).getdEtaInSeed(), evtwt);
      h1_["dPhiIn_EE_el_pre"]->Fill(goodElectrons.at(iele).getdPhiIn(), evtwt);
      h1_["Dz_EE_el_pre"]->Fill(goodElectrons.at(iele).getDz(), evtwt);
      h1_["Dxy_EE_el_pre"]->Fill(goodElectrons.at(iele).getDxy(), evtwt);
    }
  }
}

// constructors and destructor
OS2LAna::OS2LAna(const edm::ParameterSet& iConfig) :
  t_evttype               (consumes<string>  (iConfig.getParameter<edm::InputTag>("evttype"))),
  t_evtwtGen              (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtGen"))),
  t_evtwtPV               (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPV"))),
  t_evtwtPVLow            (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPVLow"))),
  t_evtwtPVHigh           (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPVHigh"))),
  t_npv                   (consumes<unsigned>(iConfig.getParameter<edm::InputTag>("npv"))),
  t_hltdecision           (consumes<bool>    (iConfig.getParameter<edm::InputTag>("hltdecision"))),
  t_hltreject             (consumes<bool>    (iConfig.getParameter<edm::InputTag>("hltreject"))),
  t_evtno                 (consumes<int>     (iConfig.getParameter<edm::InputTag>("evtno"))),
  t_runno                 (consumes<int>     (iConfig.getParameter<edm::InputTag>("runno"))),
  t_lumisec               (consumes<int>     (iConfig.getParameter<edm::InputTag>("lumisec"))),
  t_lhewtids              (consumes<vector<int> >     (iConfig.getParameter<edm::InputTag>("lhewtids"))),
  t_lhewts (consumes<vector<double> > (iConfig.getParameter<edm::InputTag>("lhewts"))),
  DilepCandParams_        (iConfig.getParameter<edm::ParameterSet> ("DilepCandParams")),
  ZCandParams_            (iConfig.getParameter<edm::ParameterSet> ("ZCandParams")),
  //BoostedZCandParams_     (iConfig.getParameter<edm::ParameterSet> ("BoostedZCandParams")),
  //GenHSelParams_          (iConfig.getParameter<edm::ParameterSet> ("GenHSelParams")),
  genParams_              (iConfig.getParameter<edm::ParameterSet> ("genParams")),
  NAK4Min_                (iConfig.getParameter<unsigned int>      ("NAK4Min")),
  HTMin_                  (iConfig.getParameter<double>            ("HTMin")),
  STMin_                  (iConfig.getParameter<double>            ("STMin")), 
  STMaxControl_           (iConfig.getParameter<double>            ("STMaxControl")), 
  skim_                   (iConfig.getParameter<bool>              ("skim")), 
  isData_                 (iConfig.getParameter<bool>              ("isData")),
  run_photons_            (iConfig.getParameter<bool>              ("run_photons")),
  filterSignal_           (iConfig.getParameter<bool>              ("filterSignal")), 
  additionalPlots_        (iConfig.getParameter<bool>              ("additionalPlots")), 
  signalType_             (iConfig.getParameter<std::string>       ("signalType")), 
  zdecayMode_             (iConfig.getParameter<std::string>       ("zdecayMode")),
  btagsf_bcUp_            (iConfig.getParameter<bool>              ("btagsf_bcUp")),
  btagsf_bcDown_          (iConfig.getParameter<bool>              ("btagsf_bcDown")),
  btagsf_lUp_             (iConfig.getParameter<bool>              ("btagsf_lUp")),
  btagsf_lDown_           (iConfig.getParameter<bool>              ("btagsf_lDown")),
  sbtagsf_bcUp_           (iConfig.getParameter<bool>              ("sbtagsf_bcUp")),
  sbtagsf_bcDown_         (iConfig.getParameter<bool>              ("sbtagsf_bcDown")),
  sbtagsf_lUp_            (iConfig.getParameter<bool>              ("sbtagsf_lUp")),
  sbtagsf_lDown_          (iConfig.getParameter<bool>              ("sbtagsf_lDown")),
  PileupUp_               (iConfig.getParameter<bool>              ("PileupUp")),
  PileupDown_             (iConfig.getParameter<bool>              ("PileupDown")),
  lheId_                  (iConfig.getParameter<int>               ("lheId")),
  applyLeptonIDSFs_       (iConfig.getParameter<bool>              ("applyLeptonIDSFs")), 
  applyLeptonTrigSFs_     (iConfig.getParameter<bool>              ("applyLeptonTrigSFs")),
  applyBTagSFs_           (iConfig.getParameter<bool>              ("applyBTagSFs")),
  applyDYNLOCorr_         (iConfig.getParameter<bool>              ("applyDYNLOCorr")),
  DYSFSyst_               (iConfig.getParameter<bool>              ("DYSFSyst")),
  tau21Shift_             (iConfig.getParameter<int>               ("tau21Shift")),
  tau32Shift_             (iConfig.getParameter<int>               ("tau32Shift")),
  WZShift_                (iConfig.getParameter<int>               ("WZShift")),
  higgsShift_             (iConfig.getParameter<int>               ("higgsShift")),
  topShift_             (iConfig.getParameter<int>               ("topShift")),
  pdfID_offset_           (iConfig.getParameter<int>               ("pdfID_offset")),
  scale_offset_           (iConfig.getParameter<int>               ("scale_offset")),
  syst_                   (iConfig.getParameter<bool>               ("syst")),
  vv_                     (iConfig.getParameter<bool>              ("vv")),
  elSyst_                 (iConfig.getParameter<int>               ("elSyst")),
  fname_DYNLOCorr_        (iConfig.getParameter<std::string>       ("File_DYNLOCorr")),
  funname_DYNLOCorr_      (iConfig.getParameter<std::string>       ("Fun_DYNLOCorr")),
  dynloewkkfact           (DYNLOEwkKfact(fname_DYNLOCorr_,funname_DYNLOCorr_)),
  lepIdSFs                (iConfig.getParameter<edm::ParameterSet> ("lepIdSFsParams")),
  lepTrigSFs              (iConfig.getParameter<edm::ParameterSet> ("lepTrigSFsParams")), 
  metmaker                (iConfig.getParameter<edm::ParameterSet> ("metselParams"),consumesCollector()),
  muonmaker               (iConfig.getParameter<edm::ParameterSet> ("muselParams"),consumesCollector()),
  electronmaker           (iConfig.getParameter<edm::ParameterSet> ("elselParams"),consumesCollector()),
  jetAK4maker             (iConfig.getParameter<edm::ParameterSet> ("jetAK4selParams"),consumesCollector()),
  jetAK4BTaggedmaker      (iConfig.getParameter<edm::ParameterSet> ("jetAK4BTaggedselParams"),consumesCollector()),
  jetAK8maker             (iConfig.getParameter<edm::ParameterSet> ("jetAK8selParams"),consumesCollector()),
  jetHTaggedmaker         (iConfig.getParameter<edm::ParameterSet> ("jetHTaggedselParams"),consumesCollector()),
  jetWTaggedmaker         (iConfig.getParameter<edm::ParameterSet> ("jetWTaggedselParams"),consumesCollector()),
  jetTopTaggedmaker       (iConfig.getParameter<edm::ParameterSet> ("jetTopTaggedselParams"),consumesCollector()),   
  genpart                 (genParams_, consumesCollector()),
  fnamebtagSF_            (iConfig.getParameter<std::string>       ("fnamebtagSF")),
  fnameSJbtagSF_          (iConfig.getParameter<std::string>       ("fnameSJbtagSF")),
  btageffmap_             (iConfig.getParameter<std::string>       ("btageffmap")),
  btagsfutils_            (new BTagSFUtils(fnamebtagSF_,BTagEntry::OP_MEDIUM,20., 1000., 20., 1000., 20., 1000.,btageffmap_)),
  sjbtagsfutils_          (new BTagSFUtils(fnameSJbtagSF_,BTagEntry::OP_LOOSE,30., 450., 30., 450., 20., 1000.,btageffmap_)),
  maketree_               (iConfig.getParameter<bool>("maketree"))

{
  produces<vlq::JetCollection>("ak4jets") ;
  produces<vlq::JetCollection>("ak8jets") ;
  produces<vlq::JetCollection>("tjets") ; 
  produces<vlq::JetCollection>("wjets") ; 
  produces<vlq::JetCollection>("hjets") ;
  produces<vlq::JetCollection>("bjets") ; 
  produces<vlq::JetCollection>("jets") ; 
  produces<vlq::CandidateCollection>("zllcands") ; 
  produces<double>("PreWeight");
  produces<double>("btagsf");
  produces<double>("btagsfbcUp");
  produces<double>("btagsfbcDown");
  produces<double>("btagsflUp");
  produces<double>("btagsflDown");
  produces<double>("sjbtagsf");
  produces<double>("sjbtagsfbcUp");
  produces<double>("sjbtagsfbcDown");
  produces<double>("sjbtagsflUp");
  produces<double>("sjbtagsflDown");
  produces<double>("finalWeight");
}


OS2LAna::~OS2LAna() {}

bool OS2LAna::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<string>   h_evttype     ; evt.getByToken(t_evttype     ,h_evttype    ) ; 
  Handle<double>   h_evtwtGen    ; evt.getByToken(t_evtwtGen    ,h_evtwtGen   ) ; 
  Handle<double>   h_evtwtPV     ; evt.getByToken(t_evtwtPV     ,h_evtwtPV    ) ; 
  Handle<double>   h_evtwtPVLow  ; evt.getByToken(t_evtwtPVLow  ,h_evtwtPVLow ) ;
  Handle<double>   h_evtwtPVHigh ; evt.getByToken(t_evtwtPVHigh ,h_evtwtPVHigh) ;  
  Handle<unsigned> h_npv         ; evt.getByToken(t_npv         ,h_npv        ) ; 
  Handle<bool>     h_hltdecision ; evt.getByToken(t_hltdecision ,h_hltdecision) ; 
  Handle<bool>     h_hltreject   ; evt.getByToken(t_hltreject   ,h_hltreject  ) ;
  Handle<int>      h_evtno       ; evt.getByToken(t_evtno       ,h_evtno      ) ; 
  Handle<int>      h_runno       ; evt.getByToken(t_runno       ,h_runno      ) ; 
  Handle<int>      h_lumisec     ; evt.getByToken(t_lumisec     ,h_lumisec    ) ; 
  Handle<vector<int> > h_lhewtids; evt.getByToken(t_lhewtids    ,h_lhewtids   ) ;
  Handle<vector<double> > h_lhewts;evt.getByToken(t_lhewts ,h_lhewts ) ;

  const int npv(*h_npv.product());
  const int evtno(*h_evtno.product()) ;
  const int runno(*h_runno.product()) ;
  const int lumisec(*h_lumisec.product()) ;
  const bool isData(evtno > 0 ? true : false) ; 
  vector<pair<int, double> > lhe_id_wts;

  
  if (!isData_){
    for (unsigned i=0; i<(*h_lhewtids.product()).size(); i++){
      int id = (*h_lhewtids.product()).at(i);
      double wt = (*h_lhewts.product()).at(i);
      lhe_id_wts.push_back(make_pair(id, wt));
    }
  }


  int signalType(-1);
  if ( !isData && filterSignal_ ) {
    if( (skim_ || maketree_ ) && signalType_.empty() ){ 
      if      (*h_evttype.product() == "EvtType_MC_bZbZ") signalType = 1; 
      else if (*h_evttype.product() == "EvtType_MC_bZbH") signalType = 2; 
      else if (*h_evttype.product() == "EvtType_MC_bZtW") signalType = 3; 
      else if (*h_evttype.product() == "EvtType_MC_bHbH") signalType = 4; 
      else if (*h_evttype.product() == "EvtType_MC_bHtW") signalType = 5; 
      else if (*h_evttype.product() == "EvtType_MC_tWtW") signalType = 6; 
      else if (*h_evttype.product() == "EvtType_MC_tZtZ") signalType = 7; 
      else if (*h_evttype.product() == "EvtType_MC_tZtZ") signalType = 7; 
      else if (*h_evttype.product() == "EvtType_MC_tZtH") signalType = 8; 
      else if (*h_evttype.product() == "EvtType_MC_tZbW") signalType = 9; 
      else if (*h_evttype.product() == "EvtType_MC_tHtH") signalType = 10; 
      else if (*h_evttype.product() == "EvtType_MC_tHbW") signalType = 11; 
      else if (*h_evttype.product() == "EvtType_MC_bWbW") signalType = 12; 
      h1_["signalEvts_all"] -> Fill(signalType);
    }
    else{
      if (*h_evttype.product()!=signalType_) return false ;
      else  h1_["signalEvts"] -> Fill(1) ;
    }
  }

  double evtwt(1.0);
  if (PileupUp_)         evtwt = (*h_evtwtGen.product()) * (*h_evtwtPVHigh.product()) ; 
  else if (PileupDown_)  evtwt = (*h_evtwtGen.product()) * (*h_evtwtPVLow.product()) ;
  else                   evtwt = (*h_evtwtGen.product()) * (*h_evtwtPV.product()) ;
  h1_["cutflow"] -> Fill(1, evtwt) ;

  if (!isData_ && !syst_ && !vv_) {
    for (unsigned i = 0; i < 9; i++){
      h1_[Form("pre_scale%d", i+1)] ->Fill(1, lhe_id_wts.at(i+scale_offset_).second);
    }
    for (unsigned i = 0; i < 100; i++) {
      h1_[Form("pre_pdf%d", i+1)] -> Fill(1, lhe_id_wts.at(i+pdfID_offset_).second);
    }
  }

  bool hltdecision_initializer(*h_hltdecision.product());
  if (run_photons_)
    hltdecision_initializer *= !(*h_hltreject.product());

  const bool hltdecision(hltdecision_initializer) ; 
  if ( hltdecision ) h1_["cutflow"] -> Fill(2, evtwt) ;
  else return false; //// Presel: HLT  

  vlq::MuonCollection goodMuons; 
  muonmaker(evt, goodMuons) ; 

  vlq::ElectronCollection goodElectrons; 
  electronmaker(evt, goodElectrons) ;

  vlq::MetCollection goodMet;
  metmaker(evt, goodMet) ;

  vlq::CandidateCollection dimuons, dielectrons, dileptons;   
  vlq::CandidateCollection zll; //generic collection

  // dilepton properties: M > 50, lead pt > 45, second pt > 25
  DileptonCandsProducer dileptonsprod(DilepCandParams_) ; 
  dileptonsprod.operator()<vlq::MuonCollection>(dimuons, goodMuons); 
  dileptonsprod.operator()<vlq::ElectronCollection>(dielectrons, goodElectrons) ; 

  //================================================================
  //First pre-selection: 1) 2 OS dileptons from boosted Z, >=3 jets
  //================================================================

  //// Dilepton candidate
  if (zdecayMode_ == "zmumu") {dileptons = dimuons; }
  else if (zdecayMode_ == "zelel") {dileptons = dielectrons;}

  if (dileptons.size() > 0) h1_["cutflow"] -> Fill(3, evtwt) ;
  else return false ; //// Presel: dileptons

  //// jets
  vlq::JetCollection goodAK4Jets;
  jetAK4maker(evt, goodAK4Jets) ;
  
  //// Get Dy EWK correction
  if ( applyDYNLOCorr_ ) {
    double EWKNLOkfact( dynloewkkfact(dileptons.at(0).getPt()) ) ; ////GetDYNLOCorr(dileptons.at(0).getPt())) ; 
    evtwt *= EWKNLOkfact ;
    
    if (DYSFSyst_) {
      if ( zdecayMode_ == "zmumu"){
        if (goodAK4Jets.size() == 3) { evtwt *= 0.94;}
        else if (goodAK4Jets.size() == 4) { evtwt*= 1.05;}
        else if (goodAK4Jets.size() == 5) { evtwt*= 1.18;}
        else if (goodAK4Jets.size() == 6) { evtwt*= 1.39;}
        else if (goodAK4Jets.size() >= 7) { evtwt*= 1.77;}
      }
      else if ( zdecayMode_ == "zelel"){
        if (goodAK4Jets.size() == 3) { evtwt *= 0.87;}
        else if (goodAK4Jets.size() == 4) { evtwt*= 1.02;}
        else if (goodAK4Jets.size() == 5) { evtwt*= 1.05;}
        else if (goodAK4Jets.size() == 6) { evtwt*= 1.24;}
        else if (goodAK4Jets.size() >= 7) { evtwt*= 1.60;}
      }
    }
  }

  //// Get lepton ID and Iso SF
  if  ( !isData_ ) {
    if (applyLeptonIDSFs_) {
      if ( zdecayMode_ == "zmumu" ) {
        evtwt *= lepIdSFs.IDSF(goodMuons.at(0).getPt(),goodMuons.at(0).getEta()) ;
        evtwt *= lepIdSFs.IDSF(goodMuons.at(1).getPt(),goodMuons.at(1).getEta()) ;
        evtwt *= lepIdSFs.IsoSF(goodMuons.at(0).getPt(),goodMuons.at(0).getEta()) ;
        evtwt *= lepIdSFs.IsoSF(goodMuons.at(1).getPt(), goodMuons.at(1).getEta()) ; 
      }
      else if ( zdecayMode_ == "zelel" ) {
        evtwt *= lepIdSFs.IDSF(goodElectrons.at(0).getPt(),goodElectrons.at(0).getEta()) ;
        evtwt *= lepIdSFs.IDSF(goodElectrons.at(1).getPt(), goodElectrons.at(1).getEta()) ;
      }
    }
    
    if (applyLeptonTrigSFs_) {
      if ( zdecayMode_ == "zmumu" ) evtwt *= lepTrigSFs(goodMuons.at(0).getPt(),goodMuons.at(0).getEta()) ; 
      else if ( zdecayMode_ == "zelel" ) {
        double leg1 = lepTrigSFs(goodElectrons.at(0).getPt(),goodElectrons.at(0).getEta()) ;
        double leg2 = lepTrigSFs(goodElectrons.at(1).getPt(),goodElectrons.at(1).getEta()) ;
        double elTrigSF = leg1 + leg2 - leg1*leg2 ;
        evtwt *= elTrigSF;
        if (goodElectrons.at(0).getPt() > 300)
          evtwt += evtwt * elSyst_ * 0.04;
        else
          evtwt += evtwt * elSyst_ * 0.02;
      }
    }
  } 

  //// Z mass candidate filter: 75 < M < 105, lead pt > 45, 2nd pt > 25, Z pt > 100
  CandidateFilter zllfilter(ZCandParams_) ; 
  zllfilter(dileptons, zll);

  // b-tagging:
  vlq::JetCollection goodBTaggedAK4Jets;
  jetAK4BTaggedmaker(evt, goodBTaggedAK4Jets) ; 

  //// jet cleaning w.r.t dileptons
  CandidateCleaner cleanjets(0.4, -1); //// The second argument is for lepton 2D iso, setting to -1 disables it
  CandidateCleaner cleanak8jets(0.8, -1); //// The second argument is for lepton 2D iso, setting to -1 disables it
  if (zdecayMode_ == "zmumu") {
    cleanjets(goodAK4Jets, goodMuons);
    cleanjets(goodBTaggedAK4Jets, goodMuons);
  }
  else if (zdecayMode_ == "zelel") {
    cleanjets(goodAK4Jets, goodElectrons);
    cleanjets(goodBTaggedAK4Jets, goodElectrons);
  } 

  //// Exactly one Z cand in event
  if(zll.size() == 1) {h1_["cutflow"] -> Fill(4, evtwt) ;}
  else return false ; //// Presel : Z->ll

  //// HT selection
  HT htak4(goodAK4Jets) ; 
  if ( htak4.getHT() > HTMin_ ) h1_["cutflow"] -> Fill(5, evtwt) ;  
  else return false ; //// Presel: HT cut

  //// At least 3 AK4 jets in event
  if (goodAK4Jets.size()  >= NAK4Min_ ) {h1_["cutflow"] -> Fill(6, evtwt) ;} 
  else return false; //// Presel N(AK4) 

  vlq::JetCollection goodAK8Jets ;
  jetAK8maker(evt, goodAK8Jets); 
  cleanak8jets(goodAK8Jets, goodMuons); 
  cleanak8jets(goodAK8Jets, goodElectrons); 

  vlq::JetCollection  goodWTaggedJets;
  jetWTaggedmaker(evt, goodWTaggedJets);
  cleanak8jets(goodWTaggedJets, goodMuons); 
  cleanak8jets(goodWTaggedJets, goodElectrons); 

  vlq::JetCollection goodHTaggedJets; 
  jetHTaggedmaker(evt, goodHTaggedJets);
  cleanak8jets(goodHTaggedJets, goodMuons); 
  cleanak8jets(goodHTaggedJets, goodElectrons); 

  double presel_wt(evtwt);
  double btagsf(1) ;
  double btagsf_bcUp(1) ; 
  double btagsf_bcDown(1) ; 
  double btagsf_lUp(1) ; 
  double btagsf_lDown(1) ; 
  if ( applyBTagSFs_ ) {
    std::vector<double>csvs;
    std::vector<double>pts;
    std::vector<double>etas;
    std::vector<int>   flhads;

    for (vlq::Jet jet : goodAK4Jets) {
      csvs.push_back(jet.getCSV()) ; 
      pts.push_back(jet.getPt()) ; 
      etas.push_back(jet.getEta()) ; 
      flhads.push_back(jet.getHadronFlavour()) ; 
    }

    btagsfutils_->getBTagSFs (csvs, pts, etas, flhads, jetAK4BTaggedmaker.idxjetCSVDiscMin_, btagsf, btagsf_bcUp, btagsf_bcDown, btagsf_lUp, btagsf_lDown) ;

    //// bTag SF along with sys. unc. options
    if (btagsf_bcUp_)
      evtwt *= btagsf_bcUp;
    else if (btagsf_bcDown_)
      evtwt *= btagsf_bcDown;
    else  if (btagsf_lUp_)
      evtwt *= btagsf_lUp;
    else if (btagsf_lDown_)
      evtwt *= btagsf_lDown;
    else
      evtwt *= btagsf;
  }

  double ST(htak4.getHT() + zll.at(0).getPt() + goodMet.at(0).getFullPt()); 

  if (goodAK4Jets.at(0).getPt() > 100 ) {
    h1_["cutflow"] -> Fill(7, presel_wt) ; 
    if (goodAK4Jets.at(1).getPt() > 50){
      h1_["cutflow"] -> Fill(8, presel_wt) ; 
      if ( goodBTaggedAK4Jets.size() > 0 ) { 
        h1_["cutflow"] -> Fill(9, evtwt) ;
        if ( ST > STMin_ ) h1_["cutflow"] -> Fill(10, evtwt) ;  
      }
    }
  } //// Completing the cutflow

  //// http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2016_245_v3.pdf
  //// SF of 0.93+/-0.09 required for top tag WP with mistag rate 1% (no subjet b tag): AN2016-245v3
  vlq::JetCollection goodTopTaggedJets;
  jetTopTaggedmaker(evt, goodTopTaggedJets);
  cleanak8jets(goodTopTaggedJets, goodMuons); 
  cleanak8jets(goodTopTaggedJets, goodElectrons);

  std::set<int> t_higgs_overlap, iso_top;
  for (unsigned i = 0; i < goodTopTaggedJets.size(); i++) {
    for (auto &higgs : goodHTaggedJets) {
      if (goodTopTaggedJets.at(i).getIndex() == higgs.getIndex())
        t_higgs_overlap.insert(i);
    }
    if (t_higgs_overlap.find(i) == t_higgs_overlap.end())
      iso_top.insert(i);
  }

  for (auto idx : t_higgs_overlap)
    h1_["overlapped_top_higgs"] -> Fill(ST, evtwt);

  double mistag_higgs(0.98), mistag_top(1.04), mistag_WZ(1.08);
  double mis_higgs_err(0.14), mis_top_err(0.07), mis_WZ_err(0.05);

  if (!isData_){                                                                                                                          
    GenParticleCollection genPartsInfo;                                                                                                                                 
    genPartsInfo = genpart(evt) ; 

    for (auto i : t_higgs_overlap) { 
      bool realTag = false;
      for (auto &gen : genPartsInfo) {
        if (abs(gen.getPdgID()) == 6 && (abs(gen.getMom0PdgID()) == 8000001 || abs(gen.getMom1PdgID()) == 8000001)) {
          if (goodTopTaggedJets.at(i).getP4().DeltaR(gen.getP4()) < 0.8) {
            realTag = true;
            break;
          }
        }
      }
      if (!isData_ && realTag) 
          evtwt *= ( 1.06 + (tau32Shift_ * .09));
      else 
          evtwt *= ( mistag_top + (topShift_ * mis_top_err));
    }

    for (auto i : iso_top) {
      bool realTag = false;
      for (auto &gen : genPartsInfo) {
        if (abs(gen.getPdgID()) == 6 && (abs(gen.getMom0PdgID()) == 8000001 || abs(gen.getMom1PdgID()) == 8000001)) {
          if (goodTopTaggedJets.at(i).getP4().DeltaR(gen.getP4()) < 0.8) {
            realTag = true;
            break;
          }
        }
      }
      if (!isData_ && realTag) 
          evtwt *= ( 1.06 + (tau32Shift_ * .09));
      else 
          evtwt *= ( mistag_top + (topShift_ * mis_top_err));
    }

    // check for W/Z & Higgs overlap
    for (auto &WZtag : goodWTaggedJets) {
      bool overlap_higgs_WZ(false), realW(false);
      for (auto &gen : genPartsInfo) {

        if ( abs(gen.getPdgID()) == 25 && (abs(gen.getMom0PdgID()) == 8000001 || abs(gen.getMom1PdgID()) == 8000001) ) {
          if (WZtag.getP4().DeltaR(gen.getP4()) < 0.8)
            overlap_higgs_WZ = true;
        }
        if ( (abs(gen.getPdgID()) == 24 || abs(gen.getPdgID()) == 23) && (abs(gen.getMom0PdgID()) == 8000001 || abs(gen.getMom1PdgID()) == 8000001) ) {
          if (WZtag.getP4().DeltaR(gen.getP4()) < 0.8)
            realW = true;
        }

      }
      if (!overlap_higgs_WZ){
        if (!isData_ && realW)
            evtwt *= ( 1.11 + (tau21Shift_ * .08) + (tau21Shift_ * 0.041 * log(WZtag.getPt() / 200)));
        else
           evtwt *= ( mistag_WZ + (WZShift_ * mis_WZ_err) ); 
      }      
      else 
        h1_["higgs_WZ_overlap"] -> Fill(ST, evtwt);
    } 
    for (auto& jet : goodHTaggedJets) {
      bool realH = false;
      for (auto& gen : genPartsInfo) {
        if ( abs(gen.getPdgID()) == 25 && (abs(gen.getMom0PdgID()) == 8000001 || abs(gen.getMom1PdgID()) == 8000001) ) {
          if (jet.getP4().DeltaR(gen.getP4()) < 0.8 && !isData_) {
            evtwt *= ( 1.11 + (tau21Shift_ * .08) + (tau21Shift_ * 0.041 * log(jet.getPt() / 200)));
            realH = true;
          }
        }
      }
      if (!realH)
        evtwt *= ( mistag_higgs + (higgsShift_ * mis_higgs_err) );
    }
  }


  double sjbtagsf(1) ;
  double sjbtagsf_bcUp(1) ; 
  double sjbtagsf_bcDown(1) ; 
  double sjbtagsf_lUp(1) ; 
  double sjbtagsf_lDown(1) ; 
  if ( applyBTagSFs_ ) {
    std::vector<double>csvs;
    std::vector<double>pts;
    std::vector<double>etas;
    std::vector<int>   flhads;

    for (vlq::Jet jet : goodHTaggedJets) {
      csvs.push_back(jet.getCSVSubjet0()) ; 
      csvs.push_back(jet.getCSVSubjet1()) ; 
      pts.push_back(jet.getPtSubjet0()) ; 
      pts.push_back(jet.getPtSubjet1()) ; 
      etas.push_back(jet.getEtaSubjet0()) ; 
      etas.push_back(jet.getEtaSubjet1()) ; 
      flhads.push_back(jet.getHadronFlavourSubjet0()) ; 
      flhads.push_back(jet.getHadronFlavourSubjet1()) ; 
    }

    sjbtagsfutils_->getBTagSFs (csvs, pts, etas, flhads, jetHTaggedmaker.idxsjCSVMin_, sjbtagsf, sjbtagsf_bcUp, sjbtagsf_bcDown, sjbtagsf_lUp, sjbtagsf_lDown) ;

    if (sjbtagsf_bcUp == 0) sjbtagsf_bcUp = 1.0;
    if (sjbtagsf_bcDown == 0) sjbtagsf_bcDown = 1.0;
    if (sjbtagsf_lUp == 0) sjbtagsf_lUp = 1.0;
    if (sjbtagsf_lDown == 0) sjbtagsf_lDown = 1.0;
    if (sjbtagsf == 0) sjbtagsf = 1.0;

    //// bTag SF along with sys. unc. options
    if (btagsf_bcUp_)
      evtwt *= sjbtagsf_bcUp;
    else if (btagsf_bcDown_)
      evtwt *= sjbtagsf_bcDown;
    else  if (btagsf_lUp_)
      evtwt *= sjbtagsf_lUp;
    else if (btagsf_lDown_)
      evtwt *= sjbtagsf_lDown;
    else
      evtwt *= sjbtagsf;

  }

  if ( skim_ ) {

    h1_["ht_preSel_bfFit"] -> Fill(htak4.getHT(), presel_wt);
    //fill jet flavor pt and eta for b-tagging efficiencies 
    if (goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50){// not sure if we should also cut on ST
      for (vlq::Jet jet : goodAK4Jets) {
        if ( abs(jet.getHadronFlavour()) == 5) h2_["pt_eta_b_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
        else if ( abs(jet.getHadronFlavour()) == 4) h2_["pt_eta_c_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
        else if ( abs(jet.getHadronFlavour()) == 0) h2_["pt_eta_l_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
      }
      if ( goodBTaggedAK4Jets.size() > 0 ){
        for (vlq::Jet jet : goodBTaggedAK4Jets) {
          if ( abs(jet.getHadronFlavour()) == 5) h2_["pt_eta_b_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
          else if ( abs(jet.getHadronFlavour()) == 4) h2_["pt_eta_c_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
          else if ( abs(jet.getHadronFlavour()) == 0) h2_["pt_eta_l_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
        }
      }
    }

    std::unique_ptr<double> ptr_fevtwt ( new double(evtwt) ) ;
    std::unique_ptr<double> ptr_prewt ( new double(presel_wt) ) ;
    std::unique_ptr<double> ptr_btagsf          ( new double(btagsf         ) ) ;
    std::unique_ptr<double> ptr_btagsf_bcUp     ( new double(btagsf_bcUp    ) ) ;
    std::unique_ptr<double> ptr_btagsf_bcDown   ( new double(btagsf_bcDown  ) ) ;
    std::unique_ptr<double> ptr_btagsf_lUp      ( new double(btagsf_lUp     ) ) ;
    std::unique_ptr<double> ptr_btagsf_lDown    ( new double(btagsf_lDown   ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf        ( new double(sjbtagsf       ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_bcUp   ( new double(sjbtagsf_bcUp  ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_bcDown ( new double(sjbtagsf_bcDown) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_lUp    ( new double(sjbtagsf_lUp   ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_lDown  ( new double(sjbtagsf_lDown ) ) ;
    evt.put(std::move(ptr_fevtwt), "finalWeight");
    evt.put(std::move(ptr_prewt), "PreWeight");
    evt.put(std::move(ptr_btagsf         ), "btagsf"         );
    evt.put(std::move(ptr_btagsf_bcUp    ), "btagsfbcUp"    );
    evt.put(std::move(ptr_btagsf_bcDown  ), "btagsfbcDown"  );
    evt.put(std::move(ptr_btagsf_lUp     ), "btagsflUp"     );
    evt.put(std::move(ptr_btagsf_lDown   ), "btagsflDown"   );
    evt.put(std::move(ptr_sjbtagsf       ), "sjbtagsf"       );
    evt.put(std::move(ptr_sjbtagsf_bcUp  ), "sjbtagsfbcUp"  );
    evt.put(std::move(ptr_sjbtagsf_bcDown), "sjbtagsfbcDown");
    evt.put(std::move(ptr_sjbtagsf_lUp   ), "sjbtagsflUp"   );
    evt.put(std::move(ptr_sjbtagsf_lDown ), "sjbtagsflDown" );

    if(goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50 && goodBTaggedAK4Jets.size() > 0 && ST > STMin_){
      std::unique_ptr<vlq::JetCollection> ptr_ak4jets( new vlq::JetCollection(goodAK4Jets) ) ;
      std::unique_ptr<vlq::JetCollection> ptr_ak8jets( new vlq::JetCollection(goodAK8Jets) ) ;
      std::unique_ptr<vlq::JetCollection> ptr_tjets( new vlq::JetCollection(goodTopTaggedJets) ) ; 
      std::unique_ptr<vlq::JetCollection> ptr_wjets( new vlq::JetCollection(goodWTaggedJets) ) ; 
      std::unique_ptr<vlq::JetCollection> ptr_hjets( new vlq::JetCollection(goodHTaggedJets) ) ;
      std::unique_ptr<vlq::JetCollection> ptr_bjets( new vlq::JetCollection(goodBTaggedAK4Jets ) ) ; 
      std::unique_ptr<vlq::JetCollection> ptr_jets ( new vlq::JetCollection(goodAK4Jets ) ) ; 
      std::unique_ptr<vlq::CandidateCollection> ptr_zllcands ( new vlq::CandidateCollection(zll) ) ; 

      evt.put(std::move(ptr_ak4jets), "ak4jets");
      evt.put(std::move(ptr_ak8jets), "ak8jets");
      evt.put(std::move(ptr_tjets), "tjets") ; 
      evt.put(std::move(ptr_wjets), "wjets") ; 
      evt.put(std::move(ptr_hjets), "hjets") ; 
      evt.put(std::move(ptr_bjets), "bjets") ; 
      evt.put(std::move(ptr_jets), "jets")  ; 
      evt.put(std::move(ptr_zllcands), "zllcands")  ;
    }   

  } //// if skim 
  else if ( !maketree_ ) { 
    
    std::string lep("");
    if(zdecayMode_ == "zmumu") {lep = "mu";}
    else if ( zdecayMode_ == "zelel") {lep = "el";}
    else edm::LogError("OS2LAna::filter") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;

    for (auto izll : zll) {
      h1_["mass_z"+lep+lep+"_pre"] -> Fill(izll.getMass(), presel_wt) ;
      h1_["pt_z"+lep+lep+"_pre"] -> Fill(izll.getPt(), presel_wt) ; 
    }
    h1_["nak4_pre"] -> Fill(goodAK4Jets.size(), presel_wt) ;
    h1_["ht_pre"] -> Fill(htak4.getHT(), presel_wt);
    h1_["st_pre"] -> Fill(ST, presel_wt) ;   

    //// Ak4 jet plots
    for(int j=0; j<3; ++j){
      h1_[Form("ptak4jet%d_pre", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), presel_wt) ; 
      h1_[Form("etaak4jet%d_pre", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), presel_wt) ;
      h1_[Form("cvsak4jet%d_pre", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), presel_wt) ;
    }
    h1_["phi_jet1MET_pre"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), presel_wt);

    if (goodAK8Jets.size() > 0) {
      h1_["ptak8jet1_pre"] -> Fill(goodAK8Jets.at(0).getPt(), evtwt);
      h1_["prunedMak8jet1_pre"] -> Fill(goodAK8Jets.at(0).getPrunedMass(), evtwt);
    }
    if (goodAK8Jets.size() > 1) {
      h1_["ptak8jet2_pre"] ->Fill(goodAK8Jets.at(1).getPt(), evtwt);
      h1_["prunedMak8jet2_pre"] -> Fill(goodAK8Jets.at(1).getPrunedMass(), evtwt);
    }
    for (auto& ak8 : goodAK8Jets) {
      h1_["ptak8_pre"] -> Fill(ak8.getPt(), evtwt);
      h1_["prunedMak8_pre"] -> Fill(ak8.getPrunedMass(), evtwt);
    }

    //// npv
    h1_["npv_noweight_pre"] -> Fill(npv, *h_evtwtGen.product()); 
    h1_["npv_pre"] -> Fill(npv, presel_wt);

    //// met
    h1_["met_pre"] -> Fill(goodMet.at(0).getFullPt(), presel_wt);
    h1_["metPhi_pre"] -> Fill(goodMet.at(0).getFullPhi(), presel_wt);

    //// Lepton specfic properties
    if ( zdecayMode_ == "zmumu" ){       
      for(int l=0; l<2; ++l){
        h1_["pt_"+lep+Form("%d_pre", l+1)]  -> Fill(goodMuons.at(l).getPt(), presel_wt) ; 
        h1_["eta_"+lep+Form("%d_pre", l+1)]  -> Fill(goodMuons.at(l).getEta(), presel_wt) ; 
      } 

      h1_["dr_mumu_pre"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), presel_wt);
    }
    else if (zdecayMode_ == "zelel" ) {
      for(int l=0; l<2; ++l){
        h1_["pt_"+lep+Form("%d_pre", l+1)]   -> Fill(goodElectrons.at(l).getPt(), presel_wt) ; 
        h1_["eta_"+lep+Form("%d_pre", l+1)]  -> Fill(goodElectrons.at(l).getEta(), presel_wt) ; 
      } 

      h1_["dr_elel_pre"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), presel_wt);
      if(additionalPlots_) fillAdditionalPlots(goodElectrons, presel_wt);
    }

    //// Z pt
    for (auto izll : zll) h1_["pt_z"+lep+lep+"_pre"] -> Fill(izll.getPt(), presel_wt) ;
    //========================================================
    // Preselection done, proceeding with control selections
    //========================================================
    
    // CR3 
    if (goodBTaggedAK4Jets.size() == 0 && ST > STMin_) {
      // CR3 control for BB
      if (goodWTaggedJets.size() > 0)
        h1_["st_bZ_boost_0b"] -> Fill(ST, evtwt);
      else {
        if (goodBTaggedAK4Jets.size() == 1)
          h1_["st_bZ_1b_0b"] -> Fill(ST, evtwt);
        else 
          h1_["st_bZ_2b_0b"] -> Fill(ST, evtwt);
      }
        
      if (goodHTaggedJets.size() > 0) 
        h1_["st_bH_boost_0b"] -> Fill(ST, evtwt);
      else {
        if (goodBTaggedAK4Jets.size() == 1) 
          h1_["st_bH_1b_0b"] -> Fill(ST, evtwt);
        else 
          h1_["st_bH_2b_0b"] -> Fill(ST, evtwt);
      }

      if (goodTopTaggedJets.size() > 0) 
        h1_["st_tW_boost_0b"] -> Fill(ST, evtwt);
      else {
        if (goodBTaggedAK4Jets.size() == 1) 
          h1_["st_tW_1b_0b"] -> Fill(ST, evtwt);
        else 
          h1_["st_tW_2b_0b"] -> Fill(ST, evtwt);
      }
    }

    //fill control plots
    if ( goodBTaggedAK4Jets.size() > 0 && ST < STMaxControl_  ) {
      for (auto izll : zll) {
        h1_["mass_z"+lep+lep+"_cnt"] -> Fill(izll.getMass(), evtwt) ;  
        h1_["pt_z"+lep+lep+"_cnt"] -> Fill(izll.getPt(), evtwt) ; 
      }
      h1_["nak4_cnt"] -> Fill(goodAK4Jets.size(), evtwt) ;
      h1_["ht_cnt"] -> Fill(htak4.getHT(), evtwt) ;
      h1_["st_cnt"] -> Fill(ST, evtwt) ;   
      h1_["npv_noweight_cnt"] -> Fill(npv, *h_evtwtGen.product()); 
      h1_["npv_cnt"] -> Fill(npv, evtwt);
      h1_["met_cnt"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
      h1_["metPhi_cnt"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);

      //lepton specfic properties
      if ( zdecayMode_ == "zmumu" ){       
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_mumu_cnt"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
      }
      else if (zdecayMode_ == "zelel" ) {
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d_cnt", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_elel_cnt"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
      }

      //ak4 jet plots
      for(int j=0; j<3; ++j){
        h1_[Form("ptak4jet%d_cnt", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
        h1_[Form("etaak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
        h1_[Form("cvsak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
      }
      h1_["phi_jet1MET_cnt"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);
      if (goodAK8Jets.size() > 0) {
        h1_["ptak8jet1_cnt"] -> Fill(goodAK8Jets.at(0).getPt(), evtwt);
        h1_["prunedMak8jet1_cnt"] -> Fill(goodAK8Jets.at(0).getPrunedMass(), evtwt);
      }
      if (goodAK8Jets.size() > 1) {
        h1_["ptak8jet2_cnt"] ->Fill(goodAK8Jets.at(1).getPt(), evtwt);
        h1_["prunedMak8jet2_cnt"] -> Fill(goodAK8Jets.at(1).getPrunedMass(), evtwt);
      }
      for (auto& ak8 : goodAK8Jets) {
        h1_["ptak8_cnt"] -> Fill(ak8.getPt(), evtwt);
        h1_["prunedMak8_cnt"] -> Fill(ak8.getPrunedMass(), evtwt);
      }

      // CR4 control for BB
      if (goodWTaggedJets.size() > 0)
        h1_["st_bZ_boost_cnt"] -> Fill(ST, evtwt);
      else {
        if (goodBTaggedAK4Jets.size() == 1)
          h1_["st_bZ_1b_cnt"] -> Fill(ST, evtwt);
        else 
          h1_["st_bZ_2b_cnt"] -> Fill(ST, evtwt);
      }
        
      if (goodHTaggedJets.size() > 0) 
        h1_["st_bH_boost_cnt"] -> Fill(ST, evtwt);
      else {
        if (goodBTaggedAK4Jets.size() == 1) 
          h1_["st_bH_1b_cnt"] -> Fill(ST, evtwt);
        else 
          h1_["st_bH_2b_cnt"] -> Fill(ST, evtwt);
      }

      if (goodTopTaggedJets.size() > 0) 
        h1_["st_tW_boost_cnt"] -> Fill(ST, evtwt);
      else {
        if (goodBTaggedAK4Jets.size() == 1) 
          h1_["st_tW_1b_cnt"] -> Fill(ST, evtwt);
        else 
          h1_["st_tW_2b_cnt"] -> Fill(ST, evtwt);
      }

      // CR4 + 1 Z-tag
      if ( goodWTaggedJets.size() > 0 ) {
        for (auto izll : zll) {
          h1_["mass_z"+lep+lep+"_cntZ"] -> Fill(izll.getMass(), evtwt) ;  
          h1_["pt_z"+lep+lep+"_cntZ"] -> Fill(izll.getPt(), evtwt) ; 
        }
        h1_["nak4_cntZ"] -> Fill(goodAK4Jets.size(), evtwt) ;
        h1_["ht_cntZ"] -> Fill(htak4.getHT(), evtwt) ;
        h1_["st_cntZ"] -> Fill(ST, evtwt) ;   
        h1_["npv_noweight_cntZ"] -> Fill(npv, *h_evtwtGen.product()); 
        h1_["npv_cntZ"] -> Fill(npv, evtwt);
        h1_["met_cntZ"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
        h1_["metPhi_cntZ"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
  
        //lepton specfic properties
        if ( zdecayMode_ == "zmumu" ){       
          for(int l=0; l<2; ++l){
            h1_["pt_"+lep+Form("%d_cntZ", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
            h1_["eta_"+lep+Form("%d_cntZ", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
          } 
          h1_["dr_mumu_cntZ"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
        }
        else if (zdecayMode_ == "zelel" ) {
          for(int l=0; l<2; ++l){
            h1_["pt_"+lep+Form("%d_cntZ", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
            h1_["eta_"+lep+Form("%d_cntZ", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
          } 
          h1_["dr_elel_cntZ"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
        }
  
        //ak4 jet plots
        for(int j=0; j<3; ++j){
          h1_[Form("ptak4jet%d_cntZ", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
          h1_[Form("etaak4jet%d_cntZ", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
          h1_[Form("cvsak4jet%d_cntZ", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
        }
        h1_["phi_jet1MET_cntZ"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);
        if (goodAK8Jets.size() > 0) {
          h1_["ptak8jet1_cntZ"] -> Fill(goodAK8Jets.at(0).getPt(), evtwt);
          h1_["prunedMak8jet1_cntZ"] -> Fill(goodAK8Jets.at(0).getPrunedMass(), evtwt);
        }
        if (goodAK8Jets.size() > 1) {
          h1_["ptak8jet2_cntZ"] ->Fill(goodAK8Jets.at(1).getPt(), evtwt);
          h1_["prunedMak8jet2_cntZ"] -> Fill(goodAK8Jets.at(1).getPrunedMass(), evtwt);
        }
        for (auto& ak8 : goodAK8Jets) {
          h1_["ptak8_cntZ"] -> Fill(ak8.getPt(), evtwt);
          h1_["prunedMak8_cntZ"] -> Fill(ak8.getPrunedMass(), evtwt);
        }
      }

      if ( goodHTaggedJets.size() > 0 ) {
        for (auto izll : zll) {
          h1_["mass_z"+lep+lep+"_cntH"] -> Fill(izll.getMass(), evtwt) ;  
          h1_["pt_z"+lep+lep+"_cntH"] -> Fill(izll.getPt(), evtwt) ; 
        }
        h1_["nak4_cntH"] -> Fill(goodAK4Jets.size(), evtwt) ;
        h1_["ht_cntH"] -> Fill(htak4.getHT(), evtwt) ;
        h1_["st_cntH"] -> Fill(ST, evtwt) ;   
        h1_["npv_noweight_cntH"] -> Fill(npv, *h_evtwtGen.product()); 
        h1_["npv_cntH"] -> Fill(npv, evtwt);
        h1_["met_cntH"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
        h1_["metPhi_cntH"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
  
        //lepton specfic properties
        if ( zdecayMode_ == "zmumu" ){       
          for(int l=0; l<2; ++l){
            h1_["pt_"+lep+Form("%d_cntH", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
            h1_["eta_"+lep+Form("%d_cntH", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
          } 
          h1_["dr_mumu_cntH"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
        }
        else if (zdecayMode_ == "zelel" ) {
          for(int l=0; l<2; ++l){
            h1_["pt_"+lep+Form("%d_cntH", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
            h1_["eta_"+lep+Form("%d_cntH", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
          } 
          h1_["dr_elel_cntH"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
        }
  
        //ak4 jet plots
        for(int j=0; j<3; ++j){
          h1_[Form("ptak4jet%d_cntH", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
          h1_[Form("etaak4jet%d_cntH", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
          h1_[Form("cvsak4jet%d_cntH", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
        }
        h1_["phi_jet1MET_cntH"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);
        if (goodAK8Jets.size() > 0) {
          h1_["ptak8jet1_cntH"] -> Fill(goodAK8Jets.at(0).getPt(), evtwt);
          h1_["prunedMak8jet1_cntH"] -> Fill(goodAK8Jets.at(0).getPrunedMass(), evtwt);
        }
        if (goodAK8Jets.size() > 1) {
          h1_["ptak8jet2_cntH"] ->Fill(goodAK8Jets.at(1).getPt(), evtwt);
          h1_["prunedMak8jet2_cntH"] -> Fill(goodAK8Jets.at(1).getPrunedMass(), evtwt);
        }
        for (auto& ak8 : goodAK8Jets) {
          h1_["ptak8_cntH"] -> Fill(ak8.getPt(), evtwt);
          h1_["prunedMak8_cntH"] -> Fill(ak8.getPrunedMass(), evtwt);
        }
      }

      if ( goodTopTaggedJets.size() > 0 ) {
        for (auto izll : zll) {
          h1_["mass_z"+lep+lep+"_cntT"] -> Fill(izll.getMass(), evtwt) ;  
          h1_["pt_z"+lep+lep+"_cntT"] -> Fill(izll.getPt(), evtwt) ; 
        }
        h1_["nak4_cntT"] -> Fill(goodAK4Jets.size(), evtwt) ;
        h1_["ht_cntT"] -> Fill(htak4.getHT(), evtwt) ;
        h1_["st_cntT"] -> Fill(ST, evtwt) ;   
        h1_["npv_noweight_cntT"] -> Fill(npv, *h_evtwtGen.product()); 
        h1_["npv_cntT"] -> Fill(npv, evtwt);
        h1_["met_cntT"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
        h1_["metPhi_cntT"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
  
        //lepton specfic properties
        if ( zdecayMode_ == "zmumu" ){       
          for(int l=0; l<2; ++l){
            h1_["pt_"+lep+Form("%d_cntT", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
            h1_["eta_"+lep+Form("%d_cntT", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
          } 
          h1_["dr_mumu_cntT"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
        }
        else if (zdecayMode_ == "zelel" ) {
          for(int l=0; l<2; ++l){
            h1_["pt_"+lep+Form("%d_cntT", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
            h1_["eta_"+lep+Form("%d_cntT", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
          } 
          h1_["dr_elel_cntT"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
        }
  
        //ak4 jet plots
        for(int j=0; j<3; ++j){
          h1_[Form("ptak4jet%d_cntT", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
          h1_[Form("etaak4jet%d_cntT", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
          h1_[Form("cvsak4jet%d_cntT", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
        }
        h1_["phi_jet1MET_cntT"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);
        if (goodAK8Jets.size() > 0) {
          h1_["ptak8jet1_cntT"] -> Fill(goodAK8Jets.at(0).getPt(), evtwt);
          h1_["prunedMak8jet1_cntT"] -> Fill(goodAK8Jets.at(0).getPrunedMass(), evtwt);
        }
        if (goodAK8Jets.size() > 1) {
          h1_["ptak8jet2_cntT"] ->Fill(goodAK8Jets.at(1).getPt(), evtwt);
          h1_["prunedMak8jet2_cntT"] -> Fill(goodAK8Jets.at(1).getPrunedMass(), evtwt);
        }
        for (auto& ak8 : goodAK8Jets) {
          h1_["ptak8_cntT"] -> Fill(ak8.getPt(), evtwt);
          h1_["prunedMak8_cntT"] -> Fill(ak8.getPrunedMass(), evtwt);
        }
      }


    } //// Control region 

    
    if ( goodAK4Jets.at(0).getPt() > 100 
        && goodAK4Jets.at(1).getPt() > 50
        && goodBTaggedAK4Jets.size() > 0
        && ST > STMin_ 
       ) { //// Signal region 

      //// fill all the plots in signal region
      for (auto izll : zll) {
        h1_["mass_z"+lep+lep] -> Fill(izll.getMass(), evtwt) ;  
        h1_["pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ; 
      }

      h1_["nak4"] -> Fill(goodAK4Jets.size(), evtwt) ;
      h1_["ht"] -> Fill(htak4.getHT(), evtwt) ;
      h1_["st"] -> Fill(ST, evtwt) ;   
      h1_["npv_noweight"] -> Fill(npv, *h_evtwtGen.product()); 
      h1_["npv"] -> Fill(npv, evtwt);
      h1_["met"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
      h1_["metPhi"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);

      //// lepton specfic properties
      if ( zdecayMode_ == "zmumu" ){       
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_mumu"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
      }
      else if (zdecayMode_ == "zelel" ) {
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_elel"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
      }
      //// ak4 jet plots
      for(int j=0; j<3; ++j){
        h1_[Form("ptak4jet%d", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
        h1_[Form("etaak4jet%d", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
        h1_[Form("cvsak4jet%d", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
      }
      h1_["phi_jet1MET"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);

      //// fill the b-tagging plots and efficiency maps
      h1_["nbjets"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
      h1_["ptbjetleading"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
      h1_["etabjetleading"] -> Fill(goodBTaggedAK4Jets.at(0).getEta(), evtwt) ;

      //// fill the additional plots
      h1_["nak8"] -> Fill(goodAK8Jets.size(), evtwt) ;
      h1_["nwjet"] -> Fill(goodWTaggedJets.size(), evtwt) ; 
      h1_["nhjet"] -> Fill(goodHTaggedJets.size(), evtwt) ; 
      h1_["ntjet"] -> Fill(goodTopTaggedJets.size(), evtwt) ; 

      if (goodAK8Jets.size() > 0) {
        h1_["ptak8leading"] -> Fill((goodAK8Jets.at(0)).getPt(), evtwt) ; 
        h1_["etaak8leading"] -> Fill((goodAK8Jets.at(0)).getEta(), evtwt) ;
        //h1_["mak8leading"] -> Fill((goodAK8Jets.at(0)).getMass(), evtwt) ; 
        //h1_["trimmedmak8leading"] -> Fill((goodAK8Jets.at(0)).getTrimmedMass(), evtwt) ;
        h1_["prunedmak8leading"] -> Fill((goodAK8Jets.at(0)).getPrunedMass(), evtwt) ;
        h1_["softdropmak8leading"] -> Fill((goodAK8Jets.at(0)).getSoftDropMass(), evtwt) ;
        h1_["subjetinessak8leading"] -> Fill((goodAK8Jets.at(0)).getTau2()/(goodAK8Jets.at(0)).getTau1(), evtwt);
      }
      if (goodAK8Jets.size() > 1) {
        h1_["ptak82nd"] -> Fill((goodAK8Jets.at(1)).getPt(), evtwt) ; 
        h1_["etaak82nd"] -> Fill((goodAK8Jets.at(1)).getEta(), evtwt) ;
        //h1_["mak82nd"] -> Fill((goodAK8Jets.at(1)).getMass(), evtwt) ; 
        //h1_["trimmedmak82nd"] -> Fill((goodAK8Jets.at(1)).getTrimmedMass(), evtwt) ;
        h1_["prunedmak82nd"] -> Fill((goodAK8Jets.at(1)).getPrunedMass(), evtwt) ;
        h1_["softdropmak82nd"] -> Fill((goodAK8Jets.at(1)).getSoftDropMass(), evtwt) ;
        h1_["subjetinessak82nd"] -> Fill((goodAK8Jets.at(1)).getTau2()/(goodAK8Jets.at(1)).getTau1(), evtwt);
      }
      
      // begin mass reconstruction secion
      fillPdfHistos("st", ST, evtwt, lhe_id_wts);
    
//      // pairs to hold <chi2, mass> 
//      pair<double, double> resReco_bZ(9999, -1), boostReco_bZ(9999, -1);
//      pair<double, double> resReco_bH(9999, -1), boostReco_bH(9999, -1);
//      pair<double, double> resReco_tW, boostReco_tW;

      if (ST > 4000) ST = 3999;    
        
      if (goodWTaggedJets.size() > 0) {
        h1_["st_bZ_boost"] -> Fill(ST, evtwt);
        fillPdfHistos("st_bZ_boost", ST, evtwt, lhe_id_wts); 
      }
      else {
        if (goodBTaggedAK4Jets.size() == 1) {
          h1_["st_bZ_1b"] -> Fill(ST, evtwt);
          fillPdfHistos("st_bZ_1b", ST, evtwt, lhe_id_wts);
        }
        else {
          h1_["st_bZ_2b"] -> Fill(ST, evtwt);
          fillPdfHistos("st_bZ_2b", ST, evtwt, lhe_id_wts);
        }
      }
        
      if (goodHTaggedJets.size() > 0) {
        h1_["st_bH_boost"] -> Fill(ST, evtwt);
        fillPdfHistos("st_bH_boost", ST, evtwt, lhe_id_wts);
      }
      else {
        if (goodBTaggedAK4Jets.size() == 1) {
          h1_["st_bH_1b"] -> Fill(ST, evtwt);
          fillPdfHistos("st_bH_1b", ST, evtwt, lhe_id_wts);
        }
        else {
          h1_["st_bH_2b"] -> Fill(ST, evtwt);
          fillPdfHistos("st_bH_2b", ST, evtwt, lhe_id_wts);
        }
      }

      if (goodTopTaggedJets.size() > 0) {
        h1_["st_tW_boost"] -> Fill(ST, evtwt);
        fillPdfHistos("st_tW_boost", ST, evtwt, lhe_id_wts);
      }
      else {
        if (goodBTaggedAK4Jets.size() == 1) {
          h1_["st_tW_1b"] -> Fill(ST, evtwt);
          fillPdfHistos("st_tW_1b", ST, evtwt, lhe_id_wts);
        }
        else {
          h1_["st_tW_2b"] -> Fill(ST, evtwt);
          fillPdfHistos("st_tW_2b", ST, evtwt, lhe_id_wts);
        }
      }

      // test area
//      std::vector<TLorentzVector> resolvedWJets, resolvedHJets;  
//      TLorentzVector test;
//      for (auto& jet1 : goodAK4Jets) {
//        for (auto& jet2 : goodAK4Jets) {
//          if (jet1.getP4() == jet2.getP4())
//            continue;
//          test = jet1.getP4() + jet2.getP4();
//          if (test.Pt() > 100 && test.M() > 80 && test.M() < 120) {
//            auto isbTag = std::find_if(goodBTaggedAK4Jets.begin(), goodBTaggedAK4Jets.end(), [jet2](vlq::Jet comp) { return comp.getP4() == jet2.getP4(); });
//            if (isbTag != goodBTaggedAK4Jets.end())
//              resolvedHJets.push_back(test);
//          }
//          else if (test.Pt() > 100 && test.M() > 70 && test.M() < 120)
//            resolvedWJets.push_back(test);
//        }
//      }
      
      if (goodTopTaggedJets.size() > 0) {
        h1_["t_category"] -> Fill(ST, evtwt);
        fillPdfHistos("t_category", ST, evtwt, lhe_id_wts);
      }
      else if (goodHTaggedJets.size() > 0) {
        h1_["H_category"] -> Fill(ST, evtwt);
        fillPdfHistos("H_category", ST, evtwt, lhe_id_wts);
      }
      else if (goodWTaggedJets.size() > 0) {
        h1_["WZ_category"] -> Fill(ST, evtwt);
        fillPdfHistos("WZ_category", ST, evtwt, lhe_id_wts);
      }
      else if (goodBTaggedAK4Jets.size() > 1) {
        h1_["2b_category"] -> Fill(ST, evtwt);
        fillPdfHistos("2b_category", ST, evtwt, lhe_id_wts);
      }
      else if (goodBTaggedAK4Jets.size() == 1) {
        h1_["1b_category"] -> Fill(ST, evtwt);
        fillPdfHistos("1b_category", ST, evtwt, lhe_id_wts);
      }
      else {
        h1_["uncategory"] -> Fill(ST, evtwt);
      }

//      if (goodAK4Jets.size() > 3) {
//        resReco_bH = doResolvedReco(goodAK4Jets, 125., zll.at(0).getP4());
//        resReco_bZ = doResolvedReco(goodAK4Jets, 91.2, zll.at(0).getP4());
//      }
//      if (goodWTaggedJets.size() > 0)
//        boostReco_bZ = doBoostedReco(goodAK4Jets, goodWTaggedJets.at(0).getP4(), 91.2, zll.at(0).getP4(), 150.);
//      if (goodHTaggedJets.size() > 0)
//        boostReco_bH = doBoostedReco(goodAK4Jets, goodHTaggedJets.at(0).getP4(), 125., zll.at(0).getP4(), 150.);
//
//      if (goodHTaggedJets.size() > 0){
//        h1_["H_reco"] -> Fill(boostReco_bH.second, evtwt);
//        fillPdfHistos("H_reco", boostReco_bH.second, evtwt, lhe_id_wts);
//      }
//      else if (goodWTaggedJets.size() > 0) {
//        h1_["Z_reco"] -> Fill(boostReco_bZ.second, evtwt);
//        fillPdfHistos("Z_reco", boostReco_bZ.second, evtwt, lhe_id_wts);
//      }
//      else if (goodBTaggedAK4Jets.size() == 1) { 
//        if (resReco_bZ.first < resReco_bH.first) {
//          h1_["b1_reco"] -> Fill(resReco_bZ.second, evtwt);
//        fillPdfHistos("b1_reco", resReco_bZ.second, evtwt, lhe_id_wts);
//      }
//        else {
//          h1_["b1_reco"] -> Fill(resReco_bH.second, evtwt);
//        fillPdfHistos("b1_reco", resReco_bH.second, evtwt, lhe_id_wts);
//      }
//      }
//      else if (goodBTaggedAK4Jets.size() > 1) {
//        if (resReco_bZ.first < resReco_bH.first) {
//          h1_["b2_reco"] -> Fill(resReco_bZ.second, evtwt);
//        fillPdfHistos("b2_reco", resReco_bZ.second, evtwt, lhe_id_wts);
//      }
//        else {
//          h1_["b2_reco"] -> Fill(resReco_bH.second, evtwt);
//        fillPdfHistos("b2_reco", resReco_bH.second, evtwt, lhe_id_wts);
//      }
//      }
//
//      if (boostReco_bZ.first < boostReco_bH.first) {
//        h1_["boost_reco"] -> Fill(boostReco_bZ.second, evtwt);
//        fillPdfHistos("boost_reco", boostReco_bZ.second, evtwt, lhe_id_wts);
//      }
//      else {
//        h1_["boost_reco"] -> Fill(boostReco_bH.second, evtwt);
//        fillPdfHistos("boost_reco", boostReco_bH.second, evtwt, lhe_id_wts);
//      }
//
//      if (goodWTaggedJets.size() > 0) {
//        h1_["boostReco_bZ"] -> Fill(boostReco_bZ.second, evtwt);
//        fillPdfHistos("boostReco_bZ", boostReco_bZ.second, evtwt, lhe_id_wts);
//      }
//      if (goodHTaggedJets.size() > 0) {
//        h1_["boostReco_bH"] -> Fill(boostReco_bH.second, evtwt);
//        fillPdfHistos("boostReco_bH", boostReco_bH.second, evtwt, lhe_id_wts);
//      }
//
//      if (goodAK4Jets.size() > 3) {
//         
//        if (goodWTaggedJets.size() == 0) {
//          resReco_bZ = doResolvedReco(goodAK4Jets, 91.2, zll.at(0).getP4());
//          if (goodBTaggedAK4Jets.size() == 1) {
//            h1_["resReco_bZ_1b"] -> Fill(resReco_bZ.second, evtwt);
//            fillPdfHistos("resReco_bZ_1b", resReco_bZ.second, evtwt, lhe_id_wts);
//          }
//          else {
//            h1_["resReco_bZ_2b"] -> Fill(resReco_bZ.second, evtwt);
//            fillPdfHistos("resReco_bZ_2b", resReco_bZ.second, evtwt, lhe_id_wts);
//          }
//        }    // close nZjet == 0
// 
//        if (goodHTaggedJets.size() == 0) {
//          resReco_bH = doResolvedReco(goodAK4Jets, 125., zll.at(0).getP4());
//          if (goodBTaggedAK4Jets.size() == 1) {
//            h1_["resReco_bH_1b"] -> Fill(resReco_bH.second, evtwt);
//            fillPdfHistos("resReco_bH_1b", resReco_bH.second, evtwt, lhe_id_wts);
//          }
//          else {
//            h1_["resReco_bH_2b"] -> Fill(resReco_bH.second, evtwt);
//            fillPdfHistos("resReco_bH_2b", resReco_bH.second, evtwt, lhe_id_wts);
//          }
//        }    // close nHjet == 0
//
//        if (goodTopTaggedJets.size() == 0) {
//          resReco_tW = doResolvedReco(goodAK4Jets, 180., zll.at(0).getP4());
//          if (goodBTaggedAK4Jets.size() == 1) {
//            h1_["resReco_tW_1b"] -> Fill(resReco_tW.second, evtwt);
//            fillPdfHistos("resReco_tW_1b", resReco_tW.second, evtwt, lhe_id_wts);
//          }
//          else {
//            h1_["resReco_tW_2b"] -> Fill(resReco_tW.second, evtwt);
//            fillPdfHistos("resReco_tW_2b", resReco_tW.second, evtwt, lhe_id_wts);
//          }
//        }    // close nTopjet == 0
//      }    // close nak4 > 3
//      else if (goodAK4Jets.size() == 3 && goodWTaggedJets.size() == 0 && goodHTaggedJets.size() == 0 && goodTopTaggedJets.size() == 0) {
//        h1_["st_residual"] -> Fill(ST, evtwt);
//        fillPdfHistos("st_residual", ST, evtwt, lhe_id_wts);
//      }

      // begin categorization
      vlq::JetCollection ak4matchedak8, ak4nonmatched1, ak4nonmatched2, ak4nonmatched3;
      vlq::JetCollection goodAK4Jetscleaned(goodAK4Jets), Hb(goodHTaggedJets), ZB(goodWTaggedJets), D(goodTopTaggedJets);
      vlq::JetCollection W,B;
      vlq::CandidateCollection tops, BC, Z,ZB1, H, ZH, ZHb;
      for (auto ak4 : goodAK4Jetscleaned) {
        
        for (auto& Htag : goodHTaggedJets) {
          if ((Htag.getP4()).DeltaR(ak4.getP4()) < 0.8)
            ak4matchedak8.push_back(ak4);
          else
            ak4nonmatched1.push_back(ak4);
        }
        for (auto& Wtag : goodWTaggedJets) {
          if ((Wtag.getP4()).DeltaR(ak4.getP4()) < 0.8)
            ak4matchedak8.push_back(ak4);
          else
            ak4nonmatched2.push_back(ak4);
        }
        for (auto& Ttag : goodTopTaggedJets) {
          if ((Ttag.getP4()).DeltaR(ak4.getP4()) < 0.8) 
            ak4matchedak8.push_back(ak4);
          else
            ak4nonmatched3.push_back(ak4);
        }
      }

      std::vector<decltype(ak4matchedak8.begin())> toerase;
      for (auto clean_it = goodAK4Jetscleaned.begin(); clean_it != goodAK4Jetscleaned.end(); clean_it++) {
        for (auto it = ak4matchedak8.begin(); it != ak4matchedak8.end(); it++) {
          if (it->getP4() == clean_it->getP4()) {
            toerase.push_back(clean_it);
            break;
          }
        }
      }
      
      for (auto it : toerase)
        goodAK4Jetscleaned.erase(it);

      HCandsProducer h;
      if (goodAK4Jetscleaned.size() > 1)
        h.operator()(goodAK4Jetscleaned.size(), 2, goodAK4Jetscleaned, H);

      ZCandsProducer z;
      if (goodAK4Jetscleaned.size() > 1)
        z.operator()(goodAK4Jetscleaned.size(), 2, goodAK4Jetscleaned, Z);

      TopCandsProducer top, w;
      // Category BC
      TLorentzVector bc1;
      for (auto wjet : goodWTaggedJets) {
        for (auto ak4 : goodAK4Jetscleaned) {
          bc1 = wjet.getP4() + ak4.getP4();
          if (bc1.Mag() >= 120 && bc1.Mag() <= 240 && bc1.Pt() >= 150) {
            vlq::Candidate bc2(bc1);
            
            W.push_back(wjet);
            B.push_back(ak4);
            BC.push_back(bc2);
          }
        }
      }

      if (goodAK4Jetscleaned.size() > 2)
        top.operator()(goodAK4Jetscleaned.size(), 3, goodAK4Jetscleaned, tops);

      double nHcandidates = Hb.size() + H.size();
      double nzcandidates = ZB.size() + Z.size();
      double ntopcandidates = D.size() + BC.size() + tops.size();

    	//n,Z,H,B
    	//(1)
    	if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
    	  if (nHcandidates >= 1.0){
    	    if( goodBTaggedAK4Jets.size() == 1 ){
    	      h1_["cutflow4"] -> Fill(1, evtwt) ;
    	      h1_["st_sigT1Z1H1b1"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT1Z1H1b1", ST, evtwt, lhe_id_wts);
    	    }
    	    else if( goodBTaggedAK4Jets.size() >= 2 ){
    	      h1_["cutflow4"] -> Fill(2, evtwt) ;
    	      h1_["st_sigT1Z1H1b2"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT1Z1H1b2", ST, evtwt, lhe_id_wts);
    	    }
    	  }
    	  else if (nHcandidates == 0.0){
    	    if( goodBTaggedAK4Jets.size() == 1 ){
    	      h1_["cutflow4"] -> Fill(3, evtwt) ;
    	      h1_["st_sigT1Z1H0b1"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT1Z1H0b1", ST, evtwt, lhe_id_wts);
    	    }
    	    else if( goodBTaggedAK4Jets.size() >= 2 ){
    	      h1_["cutflow4"] -> Fill(4, evtwt) ;
    	      h1_["st_sigT1Z1H0b2"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT1Z1H0b2", ST, evtwt, lhe_id_wts);
    	    }
    	  }  
    	}
     
    	//(2)
    	if (ntopcandidates ==0.0    &&   nzcandidates>=1.0){
    	  if (nHcandidates >= 1.0){
    	    if( goodBTaggedAK4Jets.size() == 1 ){
    	      h1_["cutflow4"] -> Fill(5, evtwt) ;
    	      h1_["st_sigT0Z1H1b1"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT0Z1H1b1", ST, evtwt, lhe_id_wts);
    	    }
    	    else if( goodBTaggedAK4Jets.size() >= 2 ){
    	      h1_["cutflow4"] -> Fill(6, evtwt) ;
    	      h1_["st_sigT0Z1H1b2"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT0Z1H1b2", ST, evtwt, lhe_id_wts);
    	    }
    	  }
    	  else if (nHcandidates == 0.0){
    	    if( goodBTaggedAK4Jets.size() == 1 ){
    	      h1_["cutflow4"] -> Fill(7, evtwt) ;
    	      h1_["st_sigT0Z1H0b1"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT0Z1H0b1", ST, evtwt, lhe_id_wts);
    	    }
    	    else if( goodBTaggedAK4Jets.size() >= 2 ){
    	      h1_["cutflow4"] -> Fill(8, evtwt) ;
    	      h1_["st_sigT0Z1H0b2"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT0Z1H0b2", ST, evtwt, lhe_id_wts);
    	    }
    	  }
    	}
    
    	//(3)
    	if (ntopcandidates >=1.0    &&   nzcandidates==0.0){
    	  if (nHcandidates >= 1.0){
    	    if( goodBTaggedAK4Jets.size() == 1 ){
    	      h1_["cutflow4"] -> Fill(9, evtwt) ;
    	      h1_["st_sigT1Z0H1b1"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT1Z0H1b1", ST, evtwt, lhe_id_wts);
    	    }
    	    else if( goodBTaggedAK4Jets.size() >= 2 ){
    	      h1_["cutflow4"] -> Fill(10, evtwt) ;
    	      h1_["st_sigT1Z0H1b2"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT1Z0H1b2", ST, evtwt, lhe_id_wts);
    	    }
    	  }
    	  else if (nHcandidates == 0.0){
    	    if( goodBTaggedAK4Jets.size() == 1 ){
    	      h1_["cutflow4"] -> Fill(11, evtwt) ;
    	      h1_["st_sigT1Z0H0b1"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT1Z0H0b1", ST, evtwt, lhe_id_wts);
    	    }
    	    else if( goodBTaggedAK4Jets.size() >= 2 ){
    	      h1_["cutflow4"] -> Fill(12, evtwt) ;
    	      h1_["st_sigT1Z0H0b2"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT1Z0H0b2", ST, evtwt, lhe_id_wts);
    	    }
    	  }
    	}
    
    	//(4)
    	if (ntopcandidates == 0.0    &&   nzcandidates ==0.0){
    	  if (nHcandidates >= 1.0){
    	    if( goodBTaggedAK4Jets.size() == 1 ){
    	      h1_["cutflow4"] -> Fill(13, evtwt) ;
    	      h1_["st_sigT0Z0H1b1"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT0Z0H1b1", ST, evtwt, lhe_id_wts);
    	    }
    	    else if( goodBTaggedAK4Jets.size() >= 2 ){
    	      h1_["cutflow4"] -> Fill(14, evtwt) ;
    	      h1_["st_sigT0Z0H1b2"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT0Z0H1b2", ST, evtwt, lhe_id_wts);
    	    }
    	  }
    	  else if (nHcandidates == 0.0){
    	    if( goodBTaggedAK4Jets.size() == 1 ){
    	      h1_["cutflow4"] -> Fill(15, evtwt) ;
    	      h1_["st_sigT0Z0H0b1"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT0Z0H0b1", ST, evtwt, lhe_id_wts);
    	    }
    	    else if( goodBTaggedAK4Jets.size() >= 2 ){
    	      h1_["cutflow4"] -> Fill(16, evtwt) ;
    	      h1_["st_sigT0Z0H0b2"] -> Fill(ST, evtwt) ;
            fillPdfHistos("st_sigT0Z0H0b2", ST, evtwt, lhe_id_wts);
    	    }
    	  }
    	}
    } //// Signal region 
    else return false;

  } //// if not skim  and not maketree

  if ( maketree_ ) {

    os2ltree_.clearTreeVectors();

    os2ltree_.t_signalType = signalType ; 

    os2ltree_.ta_npv = npv ; 
    os2ltree_.t_evtInfoEventNumber =  evtno ; 
    os2ltree_.t_evtInfoRunNumber = runno ; 
    os2ltree_.t_evtInfoLumiBlock = lumisec ; 

    for (auto& lhe : lhe_id_wts) {
      os2ltree_.t_lhewts      .push_back(lhe.second);
      os2ltree_.t_lhewtids    .push_back(lhe.first);
    }

    for (vlq::Electron e : goodElectrons) {
      os2ltree_.t_elPt        .push_back(e.getPt()) ; 
      os2ltree_.t_elPhi       .push_back(e.getPhi());
      os2ltree_.t_elEta       .push_back(e.getEta());
      os2ltree_.t_elE         .push_back(e.getE());
      os2ltree_.t_elCharge    .push_back(e.getCharge());
      os2ltree_.t_elIso03     .push_back(e.getIso03());
      os2ltree_.t_elecIsLoose .push_back(e.getisLoose()); 
      os2ltree_.t_elecIsMedium.push_back(e.getisMedium()); 
      os2ltree_.t_elecIsTight .push_back(e.getisTight()); 
      os2ltree_.t_elecIsVeto  .push_back(e.getisVeto());
    }

    for (vlq::Muon m : goodMuons) {
      os2ltree_.t_muPt         .push_back(m.getPt());
      os2ltree_.t_muPhi        .push_back(m.getPhi());
      os2ltree_.t_muEta        .push_back(m.getEta());
      os2ltree_.t_muE          .push_back(m.getE());
      os2ltree_.t_muCharge     .push_back(m.getCharge());
      os2ltree_.t_muIso04      .push_back(m.getIso04());
      os2ltree_.t_muonIsTight  .push_back(m.getIsTightMuon()); 
      os2ltree_.t_muonIsLoose  .push_back(m.getIsLooseMuon()); 
      os2ltree_.t_muonIsGlobal .push_back(m.getIsGlobalMuon()); 
      os2ltree_.t_muonIsPFMuon .push_back(m.getIsPFMuon()); 
      os2ltree_.t_muonIsTracker.push_back(m.getIsTrackerMuon()); 
    }

    os2ltree_.t_ZllPt   .push_back(zll.at(0).getPt())   ;
    os2ltree_.t_ZllEta  .push_back(zll.at(0).getEta())   ;
    os2ltree_.t_ZllPhi  .push_back(zll.at(0).getPhi())   ;
    os2ltree_.t_ZllE    .push_back(zll.at(0).getEnergy())   ;
    os2ltree_.t_ZllMass .push_back(zll.at(0).getMass())   ;

    for (vlq::Jet jet : goodAK4Jets) {
      os2ltree_.t_jetAK4Pt           .push_back(jet.getPt());
      os2ltree_.t_jetAK4Phi          .push_back(jet.getPhi());
      os2ltree_.t_jetAK4Eta          .push_back(jet.getEta());
      os2ltree_.t_jetAK4E            .push_back((jet.getP4()).E());
      os2ltree_.t_jetAK4CSV          .push_back(jet.getCSV());
      os2ltree_.t_jetAK4Mass         .push_back(jet.getMass());
      os2ltree_.t_jetAK4HadronFlavour.push_back(jet.getHadronFlavour());
      os2ltree_.t_jetAK4PartonFlavour.push_back(jet.getPartonFlavour());
    }

    for (vlq::Jet bjet : goodBTaggedAK4Jets) {
      os2ltree_.t_jetAK4BPt            .push_back(bjet.getPt());
      os2ltree_.t_jetAK4BPhi           .push_back(bjet.getPhi());
      os2ltree_.t_jetAK4BEta           .push_back(bjet.getEta());
      os2ltree_.t_jetAK4BE             .push_back((bjet.getP4()).E());
      os2ltree_.t_jetAK4BCSV           .push_back(bjet.getCSV());
      os2ltree_.t_jetAK4BMass          .push_back(bjet.getMass());
      os2ltree_.t_jetAK4BHadronFlavour .push_back(bjet.getHadronFlavour());
      os2ltree_.t_jetAK4BPartonFlavour .push_back(bjet.getPartonFlavour());
    }

    for (vlq::Jet ak8 : goodAK8Jets) {
      os2ltree_.t_jetAK8Pt            .push_back(ak8.getPt());
      os2ltree_.t_jetAK8Phi           .push_back(ak8.getPhi());
      os2ltree_.t_jetAK8Eta           .push_back(ak8.getEta());
      os2ltree_.t_jetAK8E             .push_back((ak8.getP4()).E());
      os2ltree_.t_jetAK8CSV           .push_back(ak8.getCSV());
      os2ltree_.t_jetAK8Mass          .push_back(ak8.getMass());
      os2ltree_.t_jetAK8HadronFlavour .push_back(ak8.getHadronFlavour());
      os2ltree_.t_jetAK8PartonFlavour .push_back(ak8.getPartonFlavour());
      os2ltree_.t_jetAK8_tau1         .push_back(ak8.getTau1());
      os2ltree_.t_jetAK8_tau2         .push_back(ak8.getTau2());
      os2ltree_.t_jetAK8_tau3         .push_back(ak8.getTau3());
      os2ltree_.t_jetAK8_MassPruned   .push_back(ak8.getPrunedMass());
      os2ltree_.t_jetAK8_SoftDropMass .push_back(ak8.getSoftDropMass());
      os2ltree_.t_jetAK8_NSubJets     .push_back(ak8.getNSubjets());
    }

    for (vlq::Jet wjet : goodWTaggedJets) {
      os2ltree_.t_jetWJetPt            .push_back(wjet.getPt());
      os2ltree_.t_jetWJetPhi           .push_back(wjet.getPhi());
      os2ltree_.t_jetWJetEta           .push_back(wjet.getEta());
      os2ltree_.t_jetWJetE             .push_back((wjet.getP4()).E());
      os2ltree_.t_jetWJetCSV           .push_back(wjet.getCSV());
      os2ltree_.t_jetWJetMass          .push_back(wjet.getMass());
      os2ltree_.t_jetWJetHadronFlavour .push_back(wjet.getHadronFlavour());
      os2ltree_.t_jetWJetPartonFlavour .push_back(wjet.getPartonFlavour());
      os2ltree_.t_jetWJet_tau1         .push_back(wjet.getTau1());
      os2ltree_.t_jetWJet_tau2         .push_back(wjet.getTau2());
      os2ltree_.t_jetWJet_tau3         .push_back(wjet.getTau3());
      os2ltree_.t_jetWJet_MassPruned   .push_back(wjet.getPrunedMass());
      os2ltree_.t_jetWJet_SoftDropMass .push_back(wjet.getSoftDropMass());
      os2ltree_.t_jetWJet_NSubJets     .push_back(wjet.getNSubjets());
    }

    for (vlq::Jet hjet : goodHTaggedJets) {
      os2ltree_.t_jetHJetPt            .push_back(hjet.getPt());
      os2ltree_.t_jetHJetPhi           .push_back(hjet.getPhi());
      os2ltree_.t_jetHJetEta           .push_back(hjet.getEta());
      os2ltree_.t_jetHJetE             .push_back((hjet.getP4()).E());
      os2ltree_.t_jetHJetCSV           .push_back(hjet.getCSV());
      os2ltree_.t_jetHJetMass          .push_back(hjet.getMass());
      os2ltree_.t_jetHJetHadronFlavour .push_back(hjet.getHadronFlavour());
      os2ltree_.t_jetHJetPartonFlavour .push_back(hjet.getPartonFlavour());
      os2ltree_.t_jetHJet_tau1         .push_back(hjet.getTau1());
      os2ltree_.t_jetHJet_tau2         .push_back(hjet.getTau2());
      os2ltree_.t_jetHJet_tau3         .push_back(hjet.getTau3());
      os2ltree_.t_jetHJet_MassPruned   .push_back(hjet.getPrunedMass());
      os2ltree_.t_jetHJet_SoftDropMass .push_back(hjet.getSoftDropMass());
      os2ltree_.t_jetHJet_NSubJets     .push_back(hjet.getNSubjets());
    }

    for (vlq::Jet tjet : goodTopTaggedJets) {
      os2ltree_.t_jetTopJetPt            .push_back(tjet.getPt());
      os2ltree_.t_jetTopJetPhi           .push_back(tjet.getPhi());
      os2ltree_.t_jetTopJetEta           .push_back(tjet.getEta());
      os2ltree_.t_jetTopJetE             .push_back((tjet.getP4()).E());
      os2ltree_.t_jetTopJetCSV           .push_back(tjet.getCSV());
      os2ltree_.t_jetTopJetMass          .push_back(tjet.getMass());
      os2ltree_.t_jetTopJetHadronFlavour .push_back(tjet.getHadronFlavour());
      os2ltree_.t_jetTopJetPartonFlavour .push_back(tjet.getPartonFlavour());
      os2ltree_.t_jetTopJet_tau1         .push_back(tjet.getTau1());
      os2ltree_.t_jetTopJet_tau2         .push_back(tjet.getTau2());
      os2ltree_.t_jetTopJet_tau3         .push_back(tjet.getTau3());
      os2ltree_.t_jetTopJet_MassPruned   .push_back(tjet.getPrunedMass());
      os2ltree_.t_jetTopJet_SoftDropMass .push_back(tjet.getSoftDropMass());
      os2ltree_.t_jetTopJet_NSubJets     .push_back(tjet.getNSubjets());
    }

    os2ltree_.t_HT = htak4.getHT();
    os2ltree_.t_ST = ST;

    os2ltree_.t_presel_wt = presel_wt;
    os2ltree_.t_evtwt = evtwt;
    os2ltree_.t_btagsf = btagsf;
    os2ltree_.t_btagsf_bcUp = btagsf_bcUp;
    os2ltree_.t_btagsf_bcDown = btagsf_bcDown;
    os2ltree_.t_btagsf_lUp = btagsf_lUp;
    os2ltree_.t_btagsf_lDown = btagsf_lDown;
    os2ltree_.t_sjbtagsf = sjbtagsf;
    os2ltree_.t_sjbtagsf_bcUp = sjbtagsf_bcUp;
    os2ltree_.t_sjbtagsf_bcDown = sjbtagsf_bcDown;
    os2ltree_.t_sjbtagsf_lUp = sjbtagsf_lUp;
    os2ltree_.t_sjbtagsf_lDown = sjbtagsf_lDown;

    if ( !isData_ && filterSignal_ ) {
      vlq::GenParticleCollection vlqGen = genpart(evt) ;
      for(vlq::GenParticle p : vlqGen) { 
        os2ltree_.t_genPartPt        .push_back(p.getP4().Pt());
        os2ltree_.t_genPartPhi       .push_back(p.getP4().Phi());
        os2ltree_.t_genPartEta       .push_back(p.getP4().Eta());
        os2ltree_.t_genPartE         .push_back(p.getP4().E());
        os2ltree_.t_genPartID        .push_back(p.getPdgID());
        os2ltree_.t_genPartStatus    .push_back(p.getStatus());
        os2ltree_.t_genPartMom1ID    .push_back(p.getMom0PdgID());
        os2ltree_.t_genPartMom2ID    .push_back(p.getMom1PdgID());
        os2ltree_.t_genPartDau1ID    .push_back(p.getDau0PdgID());
        os2ltree_.t_genPartDau2ID    .push_back(p.getDau1PdgID());
        os2ltree_.t_genPartMom1Status.push_back(p.getMom0Status());
        os2ltree_.t_genPartMom2Status.push_back(p.getMom1Status());
        os2ltree_.t_genPartDau1Status.push_back(p.getDau0Status());
        os2ltree_.t_genPartDau2Status.push_back(p.getDau1Status());
      }
    }

    for( vlq::Met m : goodMet ) { 
      os2ltree_.t_metPt .push_back(m.getPt());
      os2ltree_.t_metPhi.push_back(m.getPhi());
      os2ltree_.t_metEta.push_back(m.getEta());
      os2ltree_.t_metE  .push_back(m.getP4().E());
    }

    tree_->Fill(); 
  } //// maketree

  return true ; 
}

// ------------ method called once each job just before starting event loop  ------------
void OS2LAna::beginJob() {

  if (filterSignal_){
    if(skim_ || maketree_){
      const int nCh = 12;
      const char *channel[nCh] = {"bZbZ", "bZbH", "bZtW", "bHbH", "bHtW", "tWtW",
        "tZtZ", "tZtH", "tZbW", "tHtH", "tHbW", "bWbW"};
      h1_["signalEvts_all"] = fs->make<TH1D>("signalEvts_all", "All signal events", 12, 0.5, 12.5) ;
      for (int i=1;i<=nCh;i++) h1_["signalEvts_all"]->GetXaxis()->SetBinLabel(i,channel[i-1]);
    }
    else{
      h1_["signalEvts"] = fs->make<TH1D>("signalEvts", "signal events", 2, 0.5, 2.5) ;
    }
  }

  h1_["cutflow"] = fs->make<TH1D>("cutflow", "cut flow", 10, 0.5, 10.5) ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(1 , "All") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(2  , "Trig.") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(3  , "l^{+}l^{-}") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(4 , "75 #lt M(l^{+}l^{-}) #lt 105") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(5 , "H_{T} #geq 200") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(6 , "N(AK4) #geq 3") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(7 , "leading jet pt > 100") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(8 , "2nd jet pt > 50") ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(9 , "N(b jet) #geq 1") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(10, "S_{T} #geq 1000") ; 

  if (!isData_ && !vv_ && !syst_ && maketree_) {
    for (unsigned i = 0; i < 9; i++) {
      string preName_scale        = Form("pre_scale%d", i+1);
      h1_[preName_scale.c_str()]  = fs->make<TH1D>(preName_scale.c_str(), "preScale", 100, 0., 4000.);
    }
    for (unsigned i = 0; i < 101; i++) {
      string preName_pdf          = Form("pre_pdf%d", i+1);
      h1_[preName_pdf.c_str()]    = fs->make<TH1D>(preName_pdf.c_str(), "prePDF", 100, 0., 4000.);
    }
  }

  if(skim_){
    // b-tagging efficiency maps:
    h2_["pt_eta_b_all"] = fs->make<TH2D>("pt_eta_b_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_c_all"] = fs->make<TH2D>("pt_eta_c_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_l_all"] = fs->make<TH2D>("pt_eta_l_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 

    h2_["pt_eta_b_btagged"] = fs->make<TH2D>("pt_eta_b_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_c_btagged"] = fs->make<TH2D>("pt_eta_c_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_l_btagged"] = fs->make<TH2D>("pt_eta_l_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ;

    h1_["ht_preSel_bfFit"] = fs->make<TH1D>("ht_preSel_bfFit", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ; 
  }

  if (!skim_ and !maketree_) {
    TFileDirectory pre = fs->mkdir ("pre");
    TFileDirectory sig = fs->mkdir ("sig");
    TFileDirectory cnt = fs->mkdir ("cnt");
    TFileDirectory cat = fs->mkdir ("cat");


    TFileDirectory *bookDir[3]; bookDir[0] = &pre; bookDir[1] = &cnt; bookDir[2] = &sig; bookDir[3] = &cat;
    std::vector<string> suffix = {"_pre", "_cnt", "", "_cat"};


    if (!isData_ && !vv_ && !syst_) {
      for (unsigned i = 0; i < 9; i++) {
        string preName_scale      = Form("pre_scale%d", i+1);
        string STName_scale       = Form("st_scale%d", i+1);
        string st_residual_name   = Form("st_residual_scale%d", i+1);
        string st_bZ_boost_name   = Form("st_bZ_boost_scale%d", i+1);
        string st_bH_boost_name   = Form("st_bH_boost_scale%d", i+1);
        string st_tW_boost_name   = Form("st_tW_boost_scale%d", i+1);
        string st_bZ_1b_name      = Form("st_bZ_1b_scale%d", i+1);
        string st_bZ_2b_name      = Form("st_bZ_2b_scale%d", i+1);
        string st_bH_1b_name      = Form("st_bH_1b_scale%d", i+1);
        string st_bH_2b_name      = Form("st_bH_2b_scale%d", i+1);
        string st_tW_1b_name      = Form("st_tW_1b_scale%d", i+1);
        string st_tW_2b_name      = Form("st_tW_2b_scale%d", i+1);

        string t_category_name   = Form("t_category_scale%d", i+1);
        string H_category_name   = Form("H_category_scale%d", i+1);
        string WZ_category_name   = Form("WZ_category_scale%d", i+1);
        string b2_category_name   = Form("2b_category_scale%d", i+1);
        string b1_category_name   = Form("1b_category_scale%d", i+1);

        string ST_sigT1Z1H1b1_name = Form("st_sigT1Z1H1b1_scale%d", i+1); 
        string ST_sigT1Z1H1b2_name = Form("st_sigT1Z1H1b2_scale%d", i+1); 
        string ST_sigT1Z1H0b1_name = Form("st_sigT1Z1H0b1_scale%d", i+1); 
        string ST_sigT1Z1H0b2_name = Form("st_sigT1Z1H0b2_scale%d", i+1); 
        string ST_sigT0Z1H1b1_name = Form("st_sigT0Z1H1b1_scale%d", i+1); 
        string ST_sigT0Z1H1b2_name = Form("st_sigT0Z1H1b2_scale%d", i+1); 
        string ST_sigT0Z1H0b1_name = Form("st_sigT0Z1H0b1_scale%d", i+1); 
        string ST_sigT0Z1H0b2_name = Form("st_sigT0Z1H0b2_scale%d", i+1); 
        string ST_sigT1Z0H1b1_name = Form("st_sigT1Z0H1b1_scale%d", i+1); 
        string ST_sigT1Z0H1b2_name = Form("st_sigT1Z0H1b2_scale%d", i+1); 
        string ST_sigT1Z0H0b1_name = Form("st_sigT1Z0H0b1_scale%d", i+1); 
        string ST_sigT1Z0H0b2_name = Form("st_sigT1Z0H0b2_scale%d", i+1); 
        string ST_sigT0Z0H1b1_name = Form("st_sigT0Z0H1b1_scale%d", i+1); 
        string ST_sigT0Z0H1b2_name = Form("st_sigT0Z0H1b2_scale%d", i+1); 
        string ST_sigT0Z0H0b1_name = Form("st_sigT0Z0H0b1_scale%d", i+1); 
        string ST_sigT0Z0H0b2_name = Form("st_sigT0Z0H0b2_scale%d", i+1); 

//        string H_reco_name   = Form("H_reco_scale%d", i+1);
//        string Z_reco_name   = Form("Z_reco_scale%d", i+1);
//        string b2_reco_name   = Form("b2_reco_scale%d", i+1);
//        string b1_reco_name   = Form("b1_reco_scale%d", i+1);
//        string boost_reco_name = Form("boost_reco_scale%d", i+1);

//        string boostReco_bZ_name  = Form("boostReco_bZ_scale%d", i+1);
//        string boostReco_bH_name  = Form("boostReco_bH_scale%d", i+1);
//        string boostReco_name     = Form("boostReco_scale%d", i+1);
//        string resReco_1b_name    = Form("resReco_1b_scale%d", i+1);
//        string resReco_2b_name    = Form("resReco_2b_scale%d", i+1);

        h1_[ST_sigT1Z1H1b1_name.c_str()] =cat.make<TH1D>(ST_sigT1Z1H1b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z1H1b2_name.c_str()] =cat.make<TH1D>(ST_sigT1Z1H1b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z1H0b1_name.c_str()] =cat.make<TH1D>(ST_sigT1Z1H0b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z1H0b2_name.c_str()] =cat.make<TH1D>(ST_sigT1Z1H0b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.); 
        h1_[ST_sigT0Z1H1b1_name.c_str()] =cat.make<TH1D>(ST_sigT0Z1H1b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z1H1b2_name.c_str()] =cat.make<TH1D>(ST_sigT0Z1H1b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z1H0b1_name.c_str()] =cat.make<TH1D>(ST_sigT0Z1H0b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z1H0b2_name.c_str()] =cat.make<TH1D>(ST_sigT0Z1H0b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z0H1b1_name.c_str()] =cat.make<TH1D>(ST_sigT1Z0H1b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z0H1b2_name.c_str()] =cat.make<TH1D>(ST_sigT1Z0H1b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z0H0b1_name.c_str()] =cat.make<TH1D>(ST_sigT1Z0H0b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z0H0b2_name.c_str()] =cat.make<TH1D>(ST_sigT1Z0H0b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z0H1b1_name.c_str()] =cat.make<TH1D>(ST_sigT0Z0H1b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z0H1b2_name.c_str()] =cat.make<TH1D>(ST_sigT0Z0H1b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z0H0b1_name.c_str()] =cat.make<TH1D>(ST_sigT0Z0H0b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z0H0b2_name.c_str()] =cat.make<TH1D>(ST_sigT0Z0H0b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
 

        h1_[preName_scale.c_str()]      = fs->make<TH1D>(preName_scale.c_str(), "preScale", 100, 0., 4000.);
        h1_[STName_scale.c_str()]       = fs->make<TH1D>(STName_scale.c_str(), "scaleST", 100, 0., 4000.);
        h1_[st_residual_name.c_str()]   = fs->make<TH1D>(st_residual_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_bZ_boost_name.c_str()]   = fs->make<TH1D>(st_bZ_boost_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_bH_boost_name.c_str()]   = fs->make<TH1D>(st_bH_boost_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_tW_boost_name.c_str()]   = fs->make<TH1D>(st_tW_boost_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_bZ_1b_name.c_str()]      = fs->make<TH1D>(st_bZ_1b_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_bZ_2b_name.c_str()]      = fs->make<TH1D>(st_bZ_2b_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_bH_1b_name.c_str()]      = fs->make<TH1D>(st_bH_1b_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_bH_2b_name.c_str()]      = fs->make<TH1D>(st_bH_2b_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_tW_1b_name.c_str()]      = fs->make<TH1D>(st_tW_1b_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_tW_2b_name.c_str()]      = fs->make<TH1D>(st_tW_2b_name.c_str(), "reco", 100, 0., 4000.);

        h1_[t_category_name.c_str()]   = fs->make<TH1D>(t_category_name.c_str(), "ST", 100, 0., 4000.);
        h1_[H_category_name.c_str()]   = fs->make<TH1D>(H_category_name.c_str(), "ST", 100, 0., 4000.);
        h1_[WZ_category_name.c_str()]   = fs->make<TH1D>(WZ_category_name.c_str(), "ST", 100, 0., 4000.);
        h1_[b2_category_name.c_str()]   = fs->make<TH1D>(b2_category_name.c_str(), "ST", 100, 0., 4000.);
        h1_[b1_category_name.c_str()]   = fs->make<TH1D>(b1_category_name.c_str(), "ST", 100, 0., 4000.);

//        h1_[boost_reco_name.c_str()]   = fs->make<TH1D>(boost_reco_name.c_str(), "ST", 100, 0., 4000.);
//        h1_[H_reco_name.c_str()]   = fs->make<TH1D>(H_reco_name.c_str(), "ST", 100, 0., 4000.);
//        h1_[Z_reco_name.c_str()]   = fs->make<TH1D>(Z_reco_name.c_str(), "ST", 100, 0., 4000.);
//        h1_[b2_reco_name.c_str()]   = fs->make<TH1D>(b2_reco_name.c_str(), "ST", 100, 0., 4000.);
//        h1_[b1_reco_name.c_str()]   = fs->make<TH1D>(b1_reco_name.c_str(), "ST", 100, 0., 4000.);

//        h1_[boostReco_bZ_name.c_str()]  = fs->make<TH1D>(boostReco_bZ_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[boostReco_bH_name.c_str()]  = fs->make<TH1D>(boostReco_bH_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[boostReco_name.c_str()]     = fs->make<TH1D>(boostReco_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[resReco_1b_name.c_str()]    = fs->make<TH1D>(resReco_1b_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[resReco_2b_name.c_str()]    = fs->make<TH1D>(resReco_2b_name.c_str(), "reco", 1000, 0., 3000.);
      }

      for (unsigned i = 0; i < 101; i++) {
        string preName_pdf        = Form("pre_pdf%d", i+1);
        string STName_pdf         = Form("st_pdf%d", i+1);
        string st_residual_name   = Form("st_residual_pdf%d", i+1);
        string st_bZ_boost_name   = Form("st_bZ_boost_pdf%d", i+1);
        string st_bH_boost_name   = Form("st_bH_boost_pdf%d", i+1);
        string st_tW_boost_name   = Form("st_tW_boost_pdf%d", i+1);
        string st_bZ_1b_name      = Form("st_bZ_1b_pdf%d", i+1);
        string st_bZ_2b_name      = Form("st_bZ_2b_pdf%d", i+1);
        string st_bH_1b_name      = Form("st_bH_1b_pdf%d", i+1);
        string st_bH_2b_name      = Form("st_bH_2b_pdf%d", i+1);
        string st_tW_1b_name      = Form("st_tW_1b_pdf%d", i+1);
        string st_tW_2b_name      = Form("st_tW_2b_pdf%d", i+1);

        string t_category_name   = Form("t_category_pdf%d", i+1);
        string H_category_name   = Form("H_category_pdf%d", i+1);
        string WZ_category_name   = Form("WZ_category_pdf%d", i+1);
        string b2_category_name   = Form("2b_category_pdf%d", i+1);
        string b1_category_name   = Form("1b_category_pdf%d", i+1);

        string ST_sigT1Z1H1b1_name = Form("st_sigT1Z1H1b1_pdf%d", i+1); 
        string ST_sigT1Z1H1b2_name = Form("st_sigT1Z1H1b2_pdf%d", i+1); 
        string ST_sigT1Z1H0b1_name = Form("st_sigT1Z1H0b1_pdf%d", i+1); 
        string ST_sigT1Z1H0b2_name = Form("st_sigT1Z1H0b2_pdf%d", i+1); 
        string ST_sigT0Z1H1b1_name = Form("st_sigT0Z1H1b1_pdf%d", i+1); 
        string ST_sigT0Z1H1b2_name = Form("st_sigT0Z1H1b2_pdf%d", i+1); 
        string ST_sigT0Z1H0b1_name = Form("st_sigT0Z1H0b1_pdf%d", i+1); 
        string ST_sigT0Z1H0b2_name = Form("st_sigT0Z1H0b2_pdf%d", i+1); 
        string ST_sigT1Z0H1b1_name = Form("st_sigT1Z0H1b1_pdf%d", i+1); 
        string ST_sigT1Z0H1b2_name = Form("st_sigT1Z0H1b2_pdf%d", i+1); 
        string ST_sigT1Z0H0b1_name = Form("st_sigT1Z0H0b1_pdf%d", i+1); 
        string ST_sigT1Z0H0b2_name = Form("st_sigT1Z0H0b2_pdf%d", i+1); 
        string ST_sigT0Z0H1b1_name = Form("st_sigT0Z0H1b1_pdf%d", i+1); 
        string ST_sigT0Z0H1b2_name = Form("st_sigT0Z0H1b2_pdf%d", i+1); 
        string ST_sigT0Z0H0b1_name = Form("st_sigT0Z0H0b1_pdf%d", i+1); 
        string ST_sigT0Z0H0b2_name = Form("st_sigT0Z0H0b2_pdf%d", i+1); 

//        string H_reco_name   = Form("H_reco_pdf%d", i+1);
//        string Z_reco_name   = Form("Z_reco_pdf%d", i+1);
//        string b2_reco_name   = Form("b2_reco_pdf%d", i+1);
//        string b1_reco_name   = Form("b1_reco_pdf%d", i+1);
//        string boost_reco_name = Form("boost_reco_pdf%d", i+1);

//        string resReco_bZ_1b_name = Form("resReco_bZ_1b_pdf%d", i+1);
//        string resReco_bH_1b_name = Form("resReco_bH_1b_pdf%d", i+1);
//        string resReco_tW_1b_name = Form("resReco_tW_1b_pdf%d", i+1);
//        string resReco_bZ_2b_name = Form("resReco_bZ_2b_pdf%d", i+1);
//        string resReco_bH_2b_name = Form("resReco_bH_2b_pdf%d", i+1);
//        string resReco_tW_2b_name = Form("resReco_tW_2b_pdf%d", i+1);
//        string boostReco_bZ_name  = Form("boostReco_bZ_pdf%d", i+1);
//        string boostReco_bH_name  = Form("boostReco_bH_pdf%d", i+1);
//        string boostReco_tW_name  = Form("boostReco_tW_pdf%d", i+1);
//        string resReco_1b_name    = Form("resReco_1b_pdf%d", i+1);
//        string resReco_2b_name    = Form("resReco_2b_pdf%d", i+1);
//        string boostReco_name     = Form("boostReco_pdf%d", i+1);

        h1_[preName_pdf.c_str()]        = fs->make<TH1D>(preName_pdf.c_str(), "prePDF", 100, 0., 4000.);
        h1_[STName_pdf.c_str()]         = fs->make<TH1D>(STName_pdf.c_str(), "pdfST", 100, 0., 4000.);
        h1_[st_residual_name.c_str()]   = fs->make<TH1D>(st_residual_name.c_str(), "reco", 100, 0., 4000.);

        h1_[ST_sigT1Z1H1b1_name.c_str()] =cat.make<TH1D>(ST_sigT1Z1H1b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z1H1b2_name.c_str()] =cat.make<TH1D>(ST_sigT1Z1H1b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z1H0b1_name.c_str()] =cat.make<TH1D>(ST_sigT1Z1H0b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z1H0b2_name.c_str()] =cat.make<TH1D>(ST_sigT1Z1H0b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.); 
        h1_[ST_sigT0Z1H1b1_name.c_str()] =cat.make<TH1D>(ST_sigT0Z1H1b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z1H1b2_name.c_str()] =cat.make<TH1D>(ST_sigT0Z1H1b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z1H0b1_name.c_str()] =cat.make<TH1D>(ST_sigT0Z1H0b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z1H0b2_name.c_str()] =cat.make<TH1D>(ST_sigT0Z1H0b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z0H1b1_name.c_str()] =cat.make<TH1D>(ST_sigT1Z0H1b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z0H1b2_name.c_str()] =cat.make<TH1D>(ST_sigT1Z0H1b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z0H0b1_name.c_str()] =cat.make<TH1D>(ST_sigT1Z0H0b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT1Z0H0b2_name.c_str()] =cat.make<TH1D>(ST_sigT1Z0H0b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z0H1b1_name.c_str()] =cat.make<TH1D>(ST_sigT0Z0H1b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z0H1b2_name.c_str()] =cat.make<TH1D>(ST_sigT0Z0H1b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z0H0b1_name.c_str()] =cat.make<TH1D>(ST_sigT0Z0H0b1_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);
        h1_[ST_sigT0Z0H0b2_name.c_str()] =cat.make<TH1D>(ST_sigT0Z0H0b2_name.c_str(), ";S_{T} [Gev];;" , 50, 1000.,2500.);


        h1_[st_bZ_boost_name.c_str()]   = fs->make<TH1D>(st_bZ_boost_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_bH_boost_name.c_str()]   = fs->make<TH1D>(st_bH_boost_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_tW_boost_name.c_str()]   = fs->make<TH1D>(st_tW_boost_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_bZ_1b_name.c_str()]      = fs->make<TH1D>(st_bZ_1b_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_bZ_2b_name.c_str()]      = fs->make<TH1D>(st_bZ_2b_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_bH_1b_name.c_str()]      = fs->make<TH1D>(st_bH_1b_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_bH_2b_name.c_str()]      = fs->make<TH1D>(st_bH_2b_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_tW_1b_name.c_str()]      = fs->make<TH1D>(st_tW_1b_name.c_str(), "reco", 100, 0., 4000.);
        h1_[st_tW_2b_name.c_str()]      = fs->make<TH1D>(st_tW_2b_name.c_str(), "reco", 100, 0., 4000.);

        h1_[WZ_category_name.c_str()]   = fs->make<TH1D>(WZ_category_name.c_str(), "ST", 100, 0., 4000.);
        h1_[b2_category_name.c_str()]   = fs->make<TH1D>(b2_category_name.c_str(), "ST", 100, 0., 4000.);
        h1_[b1_category_name.c_str()]   = fs->make<TH1D>(b1_category_name.c_str(), "ST", 100, 0., 4000.);
        h1_[t_category_name.c_str()]   = fs->make<TH1D>(t_category_name.c_str(), "ST", 100, 0., 4000.);
        h1_[H_category_name.c_str()]   = fs->make<TH1D>(H_category_name.c_str(), "ST", 100, 0., 4000.);

//        h1_[Z_reco_name.c_str()]   = fs->make<TH1D>(Z_reco_name.c_str(), "ST", 100, 0., 4000.);
//        h1_[b2_reco_name.c_str()]   = fs->make<TH1D>(b2_reco_name.c_str(), "ST", 100, 0., 4000.);
//        h1_[b1_reco_name.c_str()]   = fs->make<TH1D>(b1_reco_name.c_str(), "ST", 100, 0., 4000.);
//        h1_[boost_reco_name.c_str()]   = fs->make<TH1D>(boost_reco_name.c_str(), "ST", 100, 0., 4000.);
//        h1_[H_reco_name.c_str()]   = fs->make<TH1D>(H_reco_name.c_str(), "ST", 100, 0., 4000.);

//        h1_[resReco_bZ_1b_name.c_str()] = fs->make<TH1D>(resReco_bZ_1b_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[resReco_bH_1b_name.c_str()] = fs->make<TH1D>(resReco_bH_1b_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[resReco_bZ_2b_name.c_str()] = fs->make<TH1D>(resReco_bZ_2b_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[resReco_bH_2b_name.c_str()] = fs->make<TH1D>(resReco_bH_2b_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[resReco_tW_1b_name.c_str()] = fs->make<TH1D>(resReco_tW_1b_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[resReco_tW_2b_name.c_str()] = fs->make<TH1D>(resReco_tW_2b_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[boostReco_bZ_name.c_str()]  = fs->make<TH1D>(boostReco_bZ_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[boostReco_bH_name.c_str()]  = fs->make<TH1D>(boostReco_bH_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[boostReco_tW_name.c_str()]  = fs->make<TH1D>(boostReco_tW_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[resReco_1b_name.c_str()]    = fs->make<TH1D>(resReco_1b_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[resReco_2b_name.c_str()]    = fs->make<TH1D>(resReco_2b_name.c_str(), "reco", 1000, 0., 3000.);
//        h1_[boostReco_name.c_str()]     = fs->make<TH1D>(boostReco_name.c_str(), "reco", 1000, 0., 3000.);
      }
    }

    for (int i=0; i<3; i++){
      h1_[("npv_noweight"+suffix[i]).c_str()] = bookDir[i]->make<TH1D>( ("npv_noweight"+suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
      h1_[("npv"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("npv"+suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
      h1_[("nak4"+suffix[i]).c_str()] =  bookDir[i]->make<TH1D>( ("nak4"+suffix[i]).c_str(), ";N(AK4 jets);;" , 21, -0.5, 20.5) ;
      h1_[("ht"+suffix[i]).c_str()]   =  bookDir[i]->make<TH1D>( ("ht"+suffix[i]).c_str(), ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
      h1_[("st"+suffix[i]).c_str()]   =  bookDir[i]->make<TH1D>( ("st"+suffix[i]).c_str() ,";S_{T} [GeV]", 100, 0., 4000.) ;
      h1_[("met"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("met"+suffix[i]).c_str(), "MET [GeV]", 100, 0., 1000.);
      h1_[("metPhi"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("metPhi"+suffix[i]).c_str(), "#Phi(MET)", 20, -5., 5.);

      //jets
      for(int j=1; j<4; ++j){
        string jetPtName = Form("ptak4jet%d", j)+suffix[i]; string jetPtTitle  = Form(";p_{T}(%d leading AK4 jet) [GeV];;",j);
        h1_[jetPtName.c_str()] = bookDir[i]->make<TH1D>(jetPtName.c_str(), jetPtTitle.c_str(), 50, 0., 1000.) ;
        string jetEtaName = Form("etaak4jet%d", j)+suffix[i]; string jetEtaTitle  = Form(";#eta(%d leading AK4 jet) ;;",j);
        h1_[jetEtaName.c_str()] = bookDir[i]->make<TH1D>(jetEtaName.c_str(), jetEtaTitle.c_str(), 80 ,-4. ,4.) ;
        string jetCVSName = Form("cvsak4jet%d", j)+suffix[i]; string jetCVSTitle  = Form(";CVS(%d leading AK4 jet) ;;",j); 
        h1_[jetCVSName.c_str()] = bookDir[i]->make<TH1D>(jetCVSName.c_str(), jetCVSTitle.c_str(), 50 ,0. ,1.) ;
      }
      string jet1METPhiName = "phi_jet1MET"+suffix[i];
      h1_[jet1METPhiName.c_str()] = bookDir[i]->make<TH1D>(jet1METPhiName.c_str(), ";#Phi(leading jet, MET)", 20, -5., 5.) ;

      //leptons
      std::string lep("");
      if(zdecayMode_ == "zmumu") {lep = "mu";}
      else if ( zdecayMode_ == "zelel") {lep = "el";}
      else edm::LogError("OS2LAna::beginJob") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;
      string mass_Z = "mass_z"+lep+lep+suffix[i];
      h1_[mass_Z.c_str()] = bookDir[i]->make<TH1D>(mass_Z.c_str(), ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 100, 20., 220.) ;
      string dr_ll = "dr_"+lep+lep+suffix[i];  
      h1_[dr_ll.c_str()] = bookDir[i]->make<TH1D>(dr_ll.c_str(), ";#DeltaR(l^{+}l^{-});;", 40, 0., 4.) ;
      string pt_Z = "pt_z"+lep+lep+suffix[i];
      h1_[pt_Z.c_str()] = bookDir[i]->make<TH1D>(pt_Z.c_str(), ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 50, 0., 1000.) ; 
      for(int l=1; l<3; ++l){
        string lepPtName = "pt_"+lep+Form("%d",l)+suffix[i]; string lepPtTitle = Form(";p_{T}(%d leading lepton) [GeV];;",l);
        h1_[lepPtName.c_str()] = bookDir[i]->make<TH1D>(lepPtName.c_str(), lepPtTitle.c_str(), 50, 0., 800.) ;
        string lepEtaName = "eta_"+lep+Form("%d",l)+suffix[i]; string lepEtaTitle  = Form(";#eta(%d leading lepton) ;;",l);
        h1_[lepEtaName.c_str()] = bookDir[i]->make<TH1D>(lepEtaName.c_str(), lepEtaTitle.c_str(), 80, -4., 4.) ;
      }

//      h1_[("pt_"+lep+lep+"-3j").c_str()]   = sig.make<TH1D>(("pt_"+lep+lep+"-3j").c_str()  , ("pt_"+lep+lep+"-3j").c_str()  , 50 , 0. , 1000.);
//      h1_[("pt_"+lep+"_lead-3j").c_str()]  = sig.make<TH1D>(("pt_"+lep+"_lead-3j").c_str() , ("pt_"+lep+"_lead-3j").c_str() , 50 , 0. , 500.);
//      h1_[("pt_"+lep+"_2nd-3j").c_str()]   = sig.make<TH1D>(("pt_"+lep+"_2nd-3j").c_str()  , ("pt_"+lep+"_2nd-3j").c_str()  , 50 , 0. , 500.);
//      h1_[("eta_"+lep+"_lead-3j").c_str()] = sig.make<TH1D>(("eta_"+lep+"_lead-3j").c_str(), ("eta_"+lep+"_lead-3j").c_str(), 80 , -4., 4.);
//      h1_[("eta_"+lep+"_2nd-3j").c_str()]  = sig.make<TH1D>(("eta_"+lep+"_2nd-3j").c_str() , ("eta_"+lep+"_2nd-3j").c_str() , 80 , -4., 4.);
//      h1_[("m_"+lep+lep+"-3j").c_str()]    = sig.make<TH1D>(("m_"+lep+lep+"-3j").c_str()   , ("m_"+lep+lep+"-3j").c_str()   , 100, 20., 220.);
//      h1_[("dr_"+lep+lep+"-3j").c_str()]   = sig.make<TH1D>(("dr_"+lep+lep+"-3j").c_str()  , ("dr_"+lep+lep+"-3j").c_str()  , 40 , 0. , 4.);
    }

    h1_["overlapped_top_higgs"] = sig.make<TH1D>("overlapped_top_higgs", "Overlapping Top and Higgs", 100, 0., 4000.);
    h1_["overlapped_top_WZ"] = sig.make<TH1D>("overlapped_top_WZ", "Overlapping Top and WZ", 100, 0., 4000.);
    h1_["overlapped_higgs_WZ"] = sig.make<TH1D>("overlapped_higgs_WZ", "Overlapping Higgs and WZ", 100, 0., 4000.);
    h1_["higgs_WZ_overlap"] = sig.make<TH1D>("higgs_WZ_overlap", "Overlapping Higgs and WZ", 100, 0., 4000.);
    h1_["top_higgs_overlap"] = sig.make<TH1D>("top_higgs_overlap", "Overlapping Top and Higgs", 100, 0., 4000.);

    h1_["ptak8jet1_pre"] = pre.make<TH1D>("ptak8jet1_pre", ";p_{T} leading AK8 jet;;", 100, 0., 1500.);
    h1_["prunedMak8jet1_pre"] = pre.make<TH1D>("prunedMak8jet1_pre", "Pruned Mass leading AK8 Jet;M [GeV];;", 100, 0., 200.);
    h1_["ptak8jet2_pre"] = pre.make<TH1D>("ptak8jet2_pre", ";p_{T} 2nd AK8 jet;;", 100, 0., 1500.);
    h1_["prunedMak8jet2_pre"] = pre.make<TH1D>("prunedMak8jet2_pre", "Pruned Mass 2nd AK8 Jet;M [GeV];;", 100, 0., 200.);
    h1_["ptak8_pre"] = pre.make<TH1D>("ptak8_pre", ";p_{T} all AK8 jets;;", 100, 0., 1500.);
    h1_["prunedMak8_pre"] = pre.make<TH1D>("prunedMak8_pre", "Pruned Mass all AK8 Jets;M [GeV];;", 100, 0., 200.);

    h1_["ptak8jet1_cnt"] = cnt.make<TH1D>("ptak8jet1_cnt", ";p_{T} leading AK8 jet;;", 100, 0., 1500.);
    h1_["prunedMak8jet1_cnt"] = cnt.make<TH1D>("prunedMak8jet1_cnt", "Pruned Mass leading AK8 Jet;M [GeV];;", 100, 0., 200.);
    h1_["ptak8jet2_cnt"] = cnt.make<TH1D>("ptak8jet2_cnt", ";p_{T} 2nd AK8 jet;;", 100, 0., 1500.);
    h1_["prunedMak8jet2_cnt"] = cnt.make<TH1D>("prunedMak8jet2_cnt", "Pruned Mass 2nd AK8 Jet;M [GeV];;", 100, 0., 200.);
    h1_["ptak8_cnt"] = cnt.make<TH1D>("ptak8_cnt", ";p_{T} all AK8 jets;;", 100, 0., 1500.);
    h1_["prunedMak8_cnt"] = cnt.make<TH1D>("prunedMak8_cnt", "Pruned Mass all AK8 Jets;M [GeV];;", 100, 0., 200.);

    h1_["ptak8jet1_cntT"] = cnt.make<TH1D>("ptak8jet1_cntT", ";p_{T} leading AK8 jet;;", 100, 0., 1500.);
    h1_["prunedMak8jet1_cntT"] = cnt.make<TH1D>("prunedMak8jet1_cntT", "Pruned Mass leading AK8 Jet;M [GeV];;", 100, 0., 200.);
    h1_["ptak8jet2_cntT"] = cnt.make<TH1D>("ptak8jet2_cntT", ";p_{T} 2nd AK8 jet;;", 100, 0., 1500.);
    h1_["prunedMak8jet2_cntT"] = cnt.make<TH1D>("prunedMak8jet2_cntT", "Pruned Mass 2nd AK8 Jet;M [GeV];;", 100, 0., 200.);
    h1_["ptak8_cntT"] = cnt.make<TH1D>("ptak8_cntT", ";p_{T} all AK8 jets;;", 100, 0., 1500.);
    h1_["prunedMak8_cntT"] = cnt.make<TH1D>("prunedMak8_cntT", "Pruned Mass all AK8 Jets;M [GeV];;", 100, 0., 200.);

    h1_["ptak8jet1_cntH"] = cnt.make<TH1D>("ptak8jet1_cntH", ";p_{T} leading AK8 jet;;", 100, 0., 1500.);
    h1_["prunedMak8jet1_cntH"] = cnt.make<TH1D>("prunedMak8jet1_cntH", "Pruned Mass leading AK8 Jet;M [GeV];;", 100, 0., 200.);
    h1_["ptak8jet2_cntH"] = cnt.make<TH1D>("ptak8jet2_cntH", ";p_{T} 2nd AK8 jet;;", 100, 0., 1500.);
    h1_["prunedMak8jet2_cntH"] = cnt.make<TH1D>("prunedMak8jet2_cntH", "Pruned Mass 2nd AK8 Jet;M [GeV];;", 100, 0., 200.);
    h1_["ptak8_cntH"] = cnt.make<TH1D>("ptak8_cntH", ";p_{T} all AK8 jets;;", 100, 0., 1500.);
    h1_["prunedMak8_cntH"] = cnt.make<TH1D>("prunedMak8_cntH", "Pruned Mass all AK8 Jets;M [GeV];;", 100, 0., 200.);

    h1_["ptak8jet1_cntZ"] = cnt.make<TH1D>("ptak8jet1_cntZ", ";p_{T} leading AK8 jet;;", 100, 0., 1500.);
    h1_["prunedMak8jet1_cntZ"] = cnt.make<TH1D>("prunedMak8jet1_cntZ", "Pruned Mass leading AK8 Jet;M [GeV];;", 100, 0., 200.);
    h1_["ptak8jet2_cntZ"] = cnt.make<TH1D>("ptak8jet2_cntZ", ";p_{T} 2nd AK8 jet;;", 100, 0., 1500.);
    h1_["prunedMak8jet2_cntZ"] = cnt.make<TH1D>("prunedMak8jet2_cntZ", "Pruned Mass 2nd AK8 Jet;M [GeV];;", 100, 0., 200.);
    h1_["ptak8_cntZ"] = cnt.make<TH1D>("ptak8_cntZ", ";p_{T} all AK8 jets;;", 100, 0., 1500.);
    h1_["prunedMak8_cntZ"] = cnt.make<TH1D>("prunedMak8_cntZ", "Pruned Mass all AK8 Jets;M [GeV];;", 100, 0., 200.);

    std::vector<string> cnt_suffix = {"_cntT", "_cntZ", "_cntH"};

    for (int i=0; i<3; i++){
      h1_[("npv_noweight"+cnt_suffix[i]).c_str()] = cnt.make<TH1D>( ("npv_noweight"+cnt_suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
      h1_[("npv"+cnt_suffix[i]).c_str()]  =  cnt.make<TH1D>( ("npv"+cnt_suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
      h1_[("nak4"+cnt_suffix[i]).c_str()] =  cnt.make<TH1D>( ("nak4"+cnt_suffix[i]).c_str(), ";N(AK4 jets);;" , 21, -0.5, 20.5) ;
      h1_[("ht"+cnt_suffix[i]).c_str()]   =  cnt.make<TH1D>( ("ht"+cnt_suffix[i]).c_str(), ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
      h1_[("st"+cnt_suffix[i]).c_str()]   =  cnt.make<TH1D>( ("st"+cnt_suffix[i]).c_str() ,";S_{T} [GeV]", 100, 0., 4000.) ;
      h1_[("met"+cnt_suffix[i]).c_str()]  =  cnt.make<TH1D>( ("met"+cnt_suffix[i]).c_str(), "MET [GeV]", 100, 0., 1000.);
      h1_[("metPhi"+cnt_suffix[i]).c_str()]  =  cnt.make<TH1D>( ("metPhi"+cnt_suffix[i]).c_str(), "#Phi(MET)", 20, -5., 5.);

      //jets
      for(int j=1; j<4; ++j){
        string jetPtName = Form("ptak4jet%d", j)+cnt_suffix[i]; string jetPtTitle  = Form(";p_{T}(%d leading AK4 jet) [GeV];;",j);
        h1_[jetPtName.c_str()] = cnt.make<TH1D>(jetPtName.c_str(), jetPtTitle.c_str(), 50, 0., 1000.) ;
        string jetEtaName = Form("etaak4jet%d", j)+cnt_suffix[i]; string jetEtaTitle  = Form(";#eta(%d leading AK4 jet) ;;",j);
        h1_[jetEtaName.c_str()] = cnt.make<TH1D>(jetEtaName.c_str(), jetEtaTitle.c_str(), 80 ,-4. ,4.) ;
        string jetCVSName = Form("cvsak4jet%d", j)+cnt_suffix[i]; string jetCVSTitle  = Form(";CVS(%d leading AK4 jet) ;;",j); 
        h1_[jetCVSName.c_str()] = cnt.make<TH1D>(jetCVSName.c_str(), jetCVSTitle.c_str(), 50 ,0. ,1.) ;
      }
      string jet1METPhiName = "phi_jet1MET"+cnt_suffix[i];
      h1_[jet1METPhiName.c_str()] = cnt.make<TH1D>(jet1METPhiName.c_str(), ";#Phi(leading jet, MET)", 20, -5., 5.) ;

      //leptons
      std::string lep("");
      if(zdecayMode_ == "zmumu") {lep = "mu";}
      else if ( zdecayMode_ == "zelel") {lep = "el";}
      else edm::LogError("OS2LAna::beginJob") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;
      string mass_Z = "mass_z"+lep+lep+cnt_suffix[i];
      h1_[mass_Z.c_str()] = cnt.make<TH1D>(mass_Z.c_str(), ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 100, 20., 220.) ;
      string dr_ll = "dr_"+lep+lep+cnt_suffix[i];  
      h1_[dr_ll.c_str()] = cnt.make<TH1D>(dr_ll.c_str(), ";#DeltaR(l^{+}l^{-});;", 40, 0., 4.) ;
      string pt_Z = "pt_z"+lep+lep+cnt_suffix[i];
      h1_[pt_Z.c_str()] = cnt.make<TH1D>(pt_Z.c_str(), ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 50, 0., 1000.) ; 
      for(int l=1; l<3; ++l){
        string lepPtName = "pt_"+lep+Form("%d",l)+cnt_suffix[i]; string lepPtTitle = Form(";p_{T}(%d leading lepton) [GeV];;",l);
        h1_[lepPtName.c_str()] = cnt.make<TH1D>(lepPtName.c_str(), lepPtTitle.c_str(), 50, 0., 800.) ;
        string lepEtaName = "eta_"+lep+Form("%d",l)+cnt_suffix[i]; string lepEtaTitle  = Form(";#eta(%d leading lepton) ;;",l);
        h1_[lepEtaName.c_str()] = cnt.make<TH1D>(lepEtaName.c_str(), lepEtaTitle.c_str(), 80, -4., 4.) ;
      }
    }


    //additional plots
    h1_["nbjets"] = sig.make<TH1D>("nbjets", ";N(b jets);;" , 11, -0.5,10.5) ; 
    h1_["ptbjetleading"]  = sig.make<TH1D>("ptbjetleading", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etabjetleading"] = sig.make<TH1D>("etabjetleading", ";#eta(leading b jet);;" , 80 ,-4. ,4.) ;

    h1_["nak8"] = sig.make<TH1D>("nak8", ";N(AK8 jets);;" , 11, -0.5,10.5) ; 
    h1_["nwjet"] = sig.make<TH1D>("nwjet", ";N(W jets );;" , 6, -0.5,5.5) ; 
    h1_["nhjet"] = sig.make<TH1D>("nhjet", ";N(H jets );;" , 6, -0.5,5.5) ; 
    h1_["ntjet"] = sig.make<TH1D>("ntjet", ";N(top jets);;" , 6, -0.5,5.5) ; 

    h1_["ptak8leading"]  = sig.make<TH1D>("ptak8leading", ";p_{T}(leading AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak8leading"] = sig.make<TH1D>("etaak8leading", ";#eta(leading AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["softdropmak8leading"] = sig.make<TH1D>("softdropmak8leading", ";M(leading AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["subjetinessak8leading"] = sig.make<TH1D>("subjetinessak8leading", ";#tau 2/1 (leading AK8 jet)", 20, 0., 1.);
    h1_["prunedmak8leading"] = sig.make<TH1D>("prunedmak8leading", ";Pruned M(leading AK8 jet) [GeV];;", 100, 0., 200.) ;
    h1_["ptak82nd"]  = sig.make<TH1D>("ptak82nd", ";p_{T}(2nd AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak82nd"] = sig.make<TH1D>("etaak82nd", ";#eta(2nd AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["softdropmak82nd"] = sig.make<TH1D>("softdropmak82nd", ";M(2nd AK8 jet) [GeV];;" ,100 ,0., 200.) ;
    h1_["subjetinessak82nd"]  = sig.make<TH1D>("subjetinessak82nd", "#tau 2/1 (2nd AK8 jet)", 20, 0., 1.);
    h1_["prunedmak82nd"] = sig.make<TH1D>("prunedmak82nd", ";Pruned M(2nd AK8 jet) [GeV];;", 100, 0., 200.) ;

    h1_["st_bZ_boost"] = sig.make<TH1D>("st_bZ_boost", "st-bZ-boost", 100, 0., 4000.);
    h1_["st_bH_boost"] = sig.make<TH1D>("st_bH_boost", "st-bH-boost", 100, 0., 4000.);
    h1_["st_tW_boost"] = sig.make<TH1D>("st_tW_boost", "st-tW-boost", 100, 0., 4000.);
    h1_["st_bZ_1b"] = sig.make<TH1D>("st_bZ_1b", "st-bZ-1b", 100, 0., 4000.);
    h1_["st_bH_1b"] = sig.make<TH1D>("st_bH_1b", "st-bH-1b", 100, 0., 4000.);
    h1_["st_tW_1b"] = sig.make<TH1D>("st_tW_1b", "st-tW-1b", 100, 0., 4000.);
    h1_["st_bZ_2b"] = sig.make<TH1D>("st_bZ_2b", "st-bZ-2b", 100, 0., 4000.);
    h1_["st_bH_2b"] = sig.make<TH1D>("st_bH_2b", "st-bH-2b", 100, 0., 4000.);
    h1_["st_tW_2b"] = sig.make<TH1D>("st_tW_2b", "st-tW-2b", 100, 0., 4000.);
    h1_["st_residual"] = sig.make<TH1D>("st_residual", "st-residual", 100, 0., 4000.);

    h1_["st_bZ_boost_cnt"] = cnt.make<TH1D>("st_bZ_boost_cnt", "st-bZ-boost", 100, 0., 4000.);
    h1_["st_bH_boost_cnt"] = cnt.make<TH1D>("st_bH_boost_cnt", "st-bH-boost", 100, 0., 4000.);
    h1_["st_tW_boost_cnt"] = cnt.make<TH1D>("st_tW_boost_cnt", "st-tW-boost", 100, 0., 4000.);
    h1_["st_bZ_1b_cnt"] = cnt.make<TH1D>("st_bZ_1b_cnt", "st-bZ-1b", 100, 0., 4000.);
    h1_["st_bH_1b_cnt"] = cnt.make<TH1D>("st_bH_1b_cnt", "st-bH-1b", 100, 0., 4000.);
    h1_["st_tW_1b_cnt"] = cnt.make<TH1D>("st_tW_1b_cnt", "st-tW-1b", 100, 0., 4000.);
    h1_["st_bZ_2b_cnt"] = cnt.make<TH1D>("st_bZ_2b_cnt", "st-bZ-2b", 100, 0., 4000.);
    h1_["st_bH_2b_cnt"] = cnt.make<TH1D>("st_bH_2b_cnt", "st-bH-2b", 100, 0., 4000.);
    h1_["st_tW_2b_cnt"] = cnt.make<TH1D>("st_tW_2b_cnt", "st-tW-2b", 100, 0., 4000.);

    h1_["st_bZ_boost_0b"] = cnt.make<TH1D>("st_bZ_boost_0b", "st-bZ-boost", 100, 0., 4000.);
    h1_["st_bH_boost_0b"] = cnt.make<TH1D>("st_bH_boost_0b", "st-bH-boost", 100, 0., 4000.);
    h1_["st_tW_boost_0b"] = cnt.make<TH1D>("st_tW_boost_0b", "st-tW-boost", 100, 0., 4000.);
    h1_["st_bZ_1b_0b"] = cnt.make<TH1D>("st_bZ_1b_0b", "st-bZ-1b", 100, 0., 4000.);
    h1_["st_bH_1b_0b"] = cnt.make<TH1D>("st_bH_1b_0b", "st-bH-1b", 100, 0., 4000.);
    h1_["st_tW_1b_0b"] = cnt.make<TH1D>("st_tW_1b_0b", "st-tW-1b", 100, 0., 4000.);
    h1_["st_bZ_2b_0b"] = cnt.make<TH1D>("st_bZ_2b_0b", "st-bZ-2b", 100, 0., 4000.);
    h1_["st_bH_2b_0b"] = cnt.make<TH1D>("st_bH_2b_0b", "st-bH-2b", 100, 0., 4000.);
    h1_["st_tW_2b_0b"] = cnt.make<TH1D>("st_tW_2b_0b", "st-tW-2b", 100, 0., 4000.);

    h1_["t_category"] = sig.make<TH1D>("t_category", "t_category", 100, 0., 4000.);
    h1_["H_category"] = sig.make<TH1D>("H_category", "H_category", 100, 0., 4000.);
    h1_["uncategory"] = sig.make<TH1D>("uncategory", "uncategory", 100, 0., 4000.);
    h1_["WZ_category"] = sig.make<TH1D>("WZ_category", "WZ_category", 100, 0., 4000.);
    h1_["2b_category"] = sig.make<TH1D>("2b_category", "2b_category", 100, 0., 4000.);
    h1_["1b_category"] = sig.make<TH1D>("1b_category", "1b_category", 100, 0., 4000.);

    h1_["boost_reco"] = sig.make<TH1D>("boost_reco", "boost_reco", 100, 0., 4000.);
    h1_["H_reco"] = sig.make<TH1D>("H_reco", "H_reco", 100, 0., 4000.);
    h1_["Z_reco"] = sig.make<TH1D>("Z_reco", "WZ_reco", 100, 0., 4000.);
    h1_["b2_reco"] = sig.make<TH1D>("b2_reco", "b2_reco", 100, 0., 4000.);
    h1_["b1_reco"] = sig.make<TH1D>("b1_reco", "b1_reco", 100, 0., 4000.);

    h1_["boostReco_bZ"] = sig.make<TH1D>("boostReco_bZ", "Boosted Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
    h1_["boostReco_bH"] = sig.make<TH1D>("boostReco_bH", "Boosted Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
    h1_["boostReco_tW"] = sig.make<TH1D>("boostReco_tW", "Boosted Reconstruction B->tW;M_{#chi^{2}}(B);;", 1000, 0., 3000);
    h1_["resReco_bZ_1b"] = sig.make<TH1D>("resReco_bZ_1b", "Resolved Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
    h1_["resReco_bH_1b"] = sig.make<TH1D>("resReco_bH_1b", "Resolved Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
    h1_["resReco_bZ_2b"] = sig.make<TH1D>("resReco_bZ_2b", "Resolved Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
    h1_["resReco_bH_2b"] = sig.make<TH1D>("resReco_bH_2b", "Resolved Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
    h1_["resReco_tW_1b"] = sig.make<TH1D>("resReco_tW_1b", "Resolved Reconstruction B->tW;M_{#chi^{2}}(B);;", 1000, 0., 3000);
    h1_["resReco_tW_2b"] = sig.make<TH1D>("resReco_tW_2b", "Resolved Reconstruction B->tW;M_{#chi^{2}}(B);;", 1000, 0., 3000);
    h1_["boostReco"] = sig.make<TH1D>("boostReco", "Boosted Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
    h1_["resReco_1b"] = sig.make<TH1D>("resReco_1b", "Resolved Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
    h1_["resReco_2b"] = sig.make<TH1D>("resReco_2b", "Resolved Reconstruction;M_{#chi^{2}}(B);;", 1000, 0., 3000);
//
//    h1_["nak4-3j"]             = sig.make<TH1D>("nak4-3j"       , "nak4-3j"       , 21, -0.5, 20.5);
//    h1_["nak8-3j"]             = sig.make<TH1D>("nak8-3j"       , "nak8-3j"       , 11, -0.5, 10.5);
//    h1_["nhjets-3j"]           = sig.make<TH1D>("nhjets-3j"     , "nhjets-3j"     , 6,  -0.5, 5.5);
//    h1_["nzjets-3j"]           = sig.make<TH1D>("nzjets-3j"     , "nzjets-3j"     , 6,  -0.5, 5.5);
//    h1_["nbjets-3j"]           = sig.make<TH1D>("nbjets-3j"     , "nbjets-3j"     , 11, -0.5, 10.5);
//    h1_["pt_bjet-3j"]          = sig.make<TH1D>("pt_bjet-3j"    , "pt_bjet-3j"    , 50, 0., 1000.);
//    h1_["pt_ak4_lead-3j"]      = sig.make<TH1D>("pt_ak4_lead-3j", "pt_ak4_lead-3j", 50, 0., 1000.);
//    h1_["pt_ak4_2nd-3j"]       = sig.make<TH1D>("pt_ak4_2nd-3j" , "pt_ak4_2nd-3j" , 50, 0., 1000.);
//    h1_["pt_ak4_3rd-3j"]       = sig.make<TH1D>("pt_ak4_3rd-3j" , "pt_ak4_3rd-3j" , 50, 0., 1000.);
//    h1_["pt_ak4_4th-3j"]       = sig.make<TH1D>("pt_ak4_4th-3j" , "pt_ak4_4th-3j" , 50, 0., 1000.);
//    h1_["eta_ak4_lead-3j"]     = sig.make<TH1D>("eta_ak4_lead-3", "eta_ak4_lead-3", 80, -4., 4.);
//    h1_["eta_ak4_2nd-3j"]      = sig.make<TH1D>("eta_ak4_2nd-3j", "eta_ak4_2nd-3j", 80, -4., 4.);
//    h1_["eta_ak4_3rd-3j"]      = sig.make<TH1D>("eta_ak4_3rd-3j", "eta_ak4_3rd-3j", 80, -4., 4.);
//    h1_["eta_ak4_4th-3j"]      = sig.make<TH1D>("eta_ak4_4th-3j", "eta_ak4_4th-3j", 80, -4., 4.);
//    h1_["st-3j"]               = sig.make<TH1D>("st-3j"         , "st-3j"         , 100, 0., 4000.);
//    h1_["ht-3j"]               = sig.make<TH1D>("ht-3j"         , "ht-3j"         , 100, 0., 4000.);
//    h1_["met-3j"]              = sig.make<TH1D>("met-3j"        , "met-3j"        , 100, 0., 1000.);
//    h1_["npv-3j"]              = sig.make<TH1D>("npv-3j"        , "npv-3j"        , 51, -0.5, 50.5);

    ////electrons specific varaibles in EE and EB at preselection level
    if (zdecayMode_ == "zelel" && additionalPlots_){
      h1_["Eta_EB_el_pre"] = pre.make<TH1D>("Eta_EB_el_pre", ";Eta (EB);;", 100,-4,4) ;
      h1_["Eta_EE_el_pre"] = pre.make<TH1D>("Eta_EE_el_pre", ";Eta (EE);;", 100,-4,4) ;
      h1_["Iso03_EB_el_pre"] = pre.make<TH1D>("Iso03_EB_el_pre", ";Iso03 (EB);;", 100,0,0.3) ;
      h1_["Iso03_EE_el_pre"] = pre.make<TH1D>("Iso03_EE_el_pre", ";Iso03 (EE);;", 100,0,0.3) ;
      h1_["dEtaInSeed_EB_el_pre"] = pre.make<TH1D>("dEtaInSeed_EB_el_pre", ";dEtaInSeed (EB);;", 200,-0.05,0.05) ;
      h1_["dEtaInSeed_EE_el_pre"] = pre.make<TH1D>("dEtaInSeed_EE_el_pre", ";dEtaInSeed (EE);;", 200,-0.05,0.05) ;
      h1_["dPhiIn_EB_el_pre"] = pre.make<TH1D>("dPhiIn_EB_el_pre", ";dPhiIn (EB);;", 100,-0.2,0.2) ;
      h1_["dPhiIn_EE_el_pre"] = pre.make<TH1D>("dPhiIn_EE_el_pre", ";dPhiIn (EE);;", 100,-0.2,0.2);
      h1_["Dz_EB_el_pre"] = pre.make<TH1D>("Dz_EB_el_pre",";dZ (EB);;", 200,-0.1,0.1) ;
      h1_["Dz_EE_el_pre"] = pre.make<TH1D>("Dz_EE_el_pre", ";dZ (EE);;", 200,-0.4,0.4) ;
      h1_["Dxy_EB_el_pre"] = pre.make<TH1D>("Dxy_EB_el_pre", ";d0 (EB);;", 100,-0.1,0.1) ;
      h1_["Dxy_EE_el_pre"] = pre.make<TH1D>("Dxy_EE_el_pre", ";d0 (EE);;", 100,-0.1,0.1) ;
    }
  
    h1_["ptbjetleading_pre"]  = pre.make<TH1D>("ptbjetleading_pre", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["etabjetleading_pre"] = pre.make<TH1D>("etabjetleading_pre", ";#eta(leading b jet);;" , 80 ,-4. ,4.) ;
    h1_["ptbjetsubleading_pre"]  = pre.make<TH1D>("ptbjetsubleading_pre", ";p_{T}(subleading b jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["etabjetsubleading_pre"] = pre.make<TH1D>("etabjetsubleading_pre", ";#eta(subleading b jet);;" , 80 ,-4. ,4.) ;

    h1_["ptak8leading"]  = sig.make<TH1D>("ptak8leading", ";p_{T}(leading AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak8leading"] = sig.make<TH1D>("etaak8leading", ";#eta(leading AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["softdropmak8leading"] = sig.make<TH1D>("softdropmak8leading", ";M(leading AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["ptak82nd"]  = sig.make<TH1D>("ptak82nd", ";p_{T}(2nd AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak82nd"] = sig.make<TH1D>("etaak82nd", ";#eta(2nd AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["softdropmak82nd"] = sig.make<TH1D>("softdropmak82nd", ";M(2nd AK8 jet) [GeV];;" ,100 ,0., 200.) ;

      
    //for ZH(combined category)
    h1_["cutflow3"] = cat.make<TH1D>("cutflow3", "cut flow", 15, 0.5, 15.5) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(1, "no t,Z,b ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(2, "b1  ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(3, "b2") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(4, "t1ZH1 ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(5, "t0ZH1") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(6, "t1ZH0 ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(7, "t0ZH0") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(8, "t1ZH1b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(9, "t1ZH1b2 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(10, "t0ZH1b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(11, "t0ZH0b2 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(12, "t1ZH0b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(13, "t1ZH0b2 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(14, "t0ZH0b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(15, "t0ZH0b2 " ) ;
      
    h1_["cutflow4"] = cat.make<TH1D>("cutflow4", "cut flow", 16, 0.5, 16.5) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(1, "t1Z1H1b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(2, "t1Z1H1b2") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(3, "t1Z1H0b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(4, "t1Z1H0b2  ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(5, "t0Z1H1b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(6, "t0Z1H1b2 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(7, "t0Z1H0b1") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(8, "t0Z1H0b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(9, "t1Z0H1b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(10, "t1Z0H1b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(11, "t1Z0H0b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(12, "t1Z0H0b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(13, "t0Z0H1b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(14, "t0Z0H1b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(15, "t0Z0H0b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(16, "t0Z0H0b2 " ) ;
    
    //signal region
      
    //ST plots
      
    //top,Z, Higgs
      
    h1_["st_sig"] =cat.make<TH1D>("ST_sig", ";S_{T} [Gev];;" , 50, 1000.,2500.);

    h1_["st_sigT1Z1H1b1"] =cat.make<TH1D>("ST_sigT1Z1H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H1b2"] =cat.make<TH1D>("ST_sigT1Z1H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H0b1"] =cat.make<TH1D>("ST_sigT1Z1H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H0b2"] =cat.make<TH1D>("ST_sigT1Z1H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    
    h1_["st_sigT0Z1H1b1"] =cat.make<TH1D>("ST_sigT0Z1H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H1b2"] =cat.make<TH1D>("ST_sigT0Z1H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H0b1"] =cat.make<TH1D>("ST_sigT0Z1H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H0b2"] =cat.make<TH1D>("ST_sigT0Z1H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    
    h1_["st_sigT1Z0H1b1"] =cat.make<TH1D>("ST_sigT1Z0H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H1b2"] =cat.make<TH1D>("ST_sigT1Z0H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H0b1"] =cat.make<TH1D>("ST_sigT1Z0H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H0b2"] =cat.make<TH1D>("ST_sigT1Z0H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
      
    h1_["st_sigT0Z0H1b1"] =cat.make<TH1D>("ST_sigT0Z0H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H1b2"] =cat.make<TH1D>("ST_sigT0Z0H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H0b1"] =cat.make<TH1D>("ST_sigT0Z0H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H0b2"] =cat.make<TH1D>("ST_sigT0Z0H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    
  }

  if (maketree_) {
    tree_ = fs->make<TTree>("tree", "HH4b") ; 
    os2ltree_.RegisterTree(tree_) ; 
  } 

}

void OS2LAna::endJob() {

  return ; 
}

pair<double, double> OS2LAna::vector_eval(vector<pair<double, double> > vec){
  double min_value = 9999.;
  double mass = -1;
  for( unsigned ind = 0; ind < vec.size(); ++ind) {
    if (vec[ind].first < min_value){
      min_value = vec[ind].first;
      mass = vec[ind].second;
    }
  }
  return std::make_pair(min_value, mass);
}

double OS2LAna::resolvedChi2(vector<TLorentzVector> jets, TLorentzVector Leptons, double bosMass, double mass){

  double Zup = abs((jets[2] + jets[3]).M() - bosMass);
  double Zup2 = Zup * Zup;
  double term1 = Zup2 / (13.5*13.5);

  double BHup = abs((jets[1] + jets[2] + jets[3]).M() - mass);
  double BHup2 = BHup * BHup;
  double term2 = BHup2 / (78.6888 * 78.6888);

  double BLup = abs((jets[0] + Leptons).M() - mass);
  double BLup2 = BLup * BLup;
  double term3 = BLup2 / (59.0947 * 59.0947);

  double result = term1 + term2 + term3;
  //  cout << term1 << " " << term2 << " " << term3 << endl;
  return(result);
}

double OS2LAna::boostedChi2(vector<TLorentzVector> ak4Jets, TLorentzVector ak8Jet, TLorentzVector Leptons, double bosMass, double mass, double pT){
  if (ak8Jet.Pt() > pT && ak8Jet.M() != ak4Jets.at(0).M() && ak8Jet.M() != ak4Jets.at(1).M() && ak4Jets.at(0).DeltaR(ak8Jet) > 1.){
    double Zup = abs(ak8Jet.M() - bosMass);
    double Zup2 = Zup * Zup;
    double term1 = Zup2 / (12 * 12);
    
    double BHup = abs((ak8Jet + ak4Jets[0]).M() - mass);
    double BHup2 = BHup * BHup;
    double term2 = BHup2 / (78 * 78);
    
    double BLup = abs((ak4Jets[1] + Leptons).M() - mass);
    double BLup2 = BLup * BLup;
    double term3 = BLup2 / (59 * 59);
    
    double result = term1 + term2 + term3;
    return(result);
  }
  else
    return(9999.);
}

pair<double, double> OS2LAna::doBoostedReco(vlq::JetCollection &ak4Jets, TLorentzVector fatJet, double bosMass, TLorentzVector Leptons, double pT){
  pair<double, double> chi2_result;
  double loop = 100000;
  vector<pair<double, double> > chi2s;
  pair<double, double> chi2_fill;
  chi2_fill.first = 10000;
  chi2_fill.second = 10000;
  TLorentzVector ak4[3];
  int index_array[] = {0, 1, 2};
  for (int mass=0; mass<=3000; mass+=5){
    double loop_check = 100000;
    int n = 0;
    do{
      int i0 = index_array[0];
      int i1 = index_array[1];
      int i2 = index_array[2];
      vector<TLorentzVector> passToChi2;
      if (!passToChi2.empty())
        passToChi2.clear();
      if (ak4Jets.size() > 2){
        ak4[0] = ak4Jets.at(i0).getP4();
        ak4[1] = ak4Jets.at(i1).getP4();
        ak4[2] = ak4Jets.at(i2).getP4();
        passToChi2.push_back(ak4[0]);
        passToChi2.push_back(ak4[1]);
        passToChi2.push_back(ak4[2]);
        n = 3;
      }
      else if (ak4Jets.size() == 2){
        ak4[0] = ak4Jets.at(i0).getP4();
        ak4[1] = ak4Jets.at(i1).getP4();
        passToChi2.push_back(ak4[0]);
        passToChi2.push_back(ak4[1]);
        n = 2;
      }
      else
        continue;

      loop = boostedChi2(passToChi2, fatJet, Leptons, bosMass,  mass, pT);

      if (loop < loop_check){
        loop_check = loop;
        chi2_result.first = loop_check;
        chi2_result.second = mass;
      }
    }
    while(std::next_permutation(index_array, index_array + n));

    chi2s.push_back(chi2_result);
  }

  chi2_fill = vector_eval(chi2s);
  return(chi2_fill);
}

pair<double, double> OS2LAna::doResolvedReco(vlq::JetCollection &collection, double bosMass, TLorentzVector Leptons){
  int next = 0;
  for(int i=0; i<4; i++){
    if (collection.at(i).getPt() > 0)
      ++next;
  }

  pair<double, double> chi2_result;
  double loop = 100000;
  vector<pair<double,double> > chi2s;

  pair<double, double> chi2_fill;
  chi2_fill.first = 10000;
  chi2_fill.second = 10000;

  TLorentzVector Jets[5];

  int index_array[] = {0, 1, 2, 3, 4};

  for (int mass = 0; mass <= 3000; mass+=5){

    double loop_check = 10000;
    int n = 0;

    do{
      int i0 = index_array[0];
      int i1 = index_array[1];
      int i2 = index_array[2];
      int i3 = index_array[3];
      int i4 = index_array[4];
      vector<TLorentzVector> passToChi2;
      if (!passToChi2.empty())
        passToChi2.clear();
      if (collection.size() > 4){

        Jets[0] = collection.at(i0).getP4();
        Jets[1] = collection.at(i1).getP4();
        Jets[2] = collection.at(i2).getP4();
        Jets[3] = collection.at(i3).getP4();
        Jets[4] = collection.at(i4).getP4();

        passToChi2.push_back(Jets[0]);
        passToChi2.push_back(Jets[1]);
        passToChi2.push_back(Jets[2]);
        passToChi2.push_back(Jets[3]);
        passToChi2.push_back(Jets[4]);
        n = 5;

      }
      else if (collection.size() == 4){

        Jets[0] = collection.at(i0).getP4();
        Jets[1] = collection.at(i1).getP4();
        Jets[2] = collection.at(i2).getP4();
        Jets[3] = collection.at(i3).getP4();

        passToChi2.push_back(Jets[0]);
        passToChi2.push_back(Jets[1]);
        passToChi2.push_back(Jets[2]);
        passToChi2.push_back(Jets[3]);
        Jets[2] = collection.at(i2).getP4();
        Jets[3] = collection.at(i3).getP4();

        passToChi2.push_back(Jets[0]);
        passToChi2.push_back(Jets[1]);
        passToChi2.push_back(Jets[2]);
        passToChi2.push_back(Jets[3]);
      
        n = 4;
      }
      else
        continue;

      loop = resolvedChi2(passToChi2, Leptons, bosMass,  mass);

      if (loop < loop_check){
        loop_check = loop;
        chi2_result.first = loop_check;
        chi2_result.second = mass;
      }
    }
    while(std::next_permutation(index_array, index_array + n));

    chi2s.push_back(chi2_result);
  }

  chi2_fill = vector_eval(chi2s);
return(chi2_fill);
}

void OS2LAna::fillPdfHistos(string name, double value, double evtwt, vector<pair<int,double>> &lhe_id_wts) {

  if (false) {
  if (!syst_ && !isData_ && !vv_) {

    for (unsigned i = 0; i < 101; i++) 
      h1_[Form((name+"_pdf%d").c_str(), i+1)]   -> Fill(value, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
    for (unsigned i = 0; i < 9; i++)
      h1_[Form((name+"_scale%d").c_str(), i+1)] -> Fill(value, evtwt*lhe_id_wts.at(i+scale_offset_).second);

  }
  }
}

DEFINE_FWK_MODULE(OS2LAna);
