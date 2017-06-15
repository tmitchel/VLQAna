// -*- C++ -*-
//
// Package:    Analysis/VLQAna
// Class:      VLQAna
// 
/**\class VLQAna VLQAna.cc Analysis/VLQAna/plugins/VLQAna.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Devdatta Majumder
//         Created:  Fri, 27 Feb 2015 16:09:10 GMT
//
//

#include <iostream>
#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisDataFormats/BoostedObjects/interface/GenParticleWithDaughters.h"
#include "AnalysisDataFormats/BoostedObjects/interface/Jet.h"

#include "Analysis/VLQAna/interface/HT.h"
#include "Analysis/VLQAna/interface/ElectronMaker.h"
#include "Analysis/VLQAna/interface/MuonMaker.h"
#include "Analysis/VLQAna/interface/JetMaker.h"
#include "Analysis/VLQAna/interface/Utilities.h"
#include "Analysis/VLQAna/interface/CandidateCleaner.h"
#include "Analysis/VLQAna/interface/TtHTree.h"
#include "Analysis/VLQAna/interface/BTagSFUtils.h"

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include <sstream>

class VLQAna : public edm::EDFilter {
  public:
    explicit VLQAna(const edm::ParameterSet&);
    ~VLQAna();

  private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    double getHTReweightingSF(double ht, double err); 

    // ----------member data ---------------------------
    edm::EDGetTokenT<int>            t_runno      ;
    edm::EDGetTokenT<int>            t_lumisec    ;
    edm::EDGetTokenT<int>            t_evtno      ;
    edm::EDGetTokenT<bool>           t_isData     ;
    edm::EDGetTokenT<bool>           t_hltdecision;
    edm::EDGetTokenT<string>         t_evttype    ;
    edm::EDGetTokenT<double>         t_evtwtGen   ;
    edm::EDGetTokenT<double>         t_evtwtPV    ;
    edm::EDGetTokenT<double>         t_evtwtPVLow ;
    edm::EDGetTokenT<double>         t_evtwtPVHigh;
    edm::EDGetTokenT<unsigned>       t_npv        ;
    edm::EDGetTokenT<int>            t_npuTrue    ;
    edm::EDGetTokenT<double>         t_htHat      ;
    edm::EDGetTokenT<vector<int>>    t_lhewtids   ;
    edm::EDGetTokenT<vector<double>> t_lhewts     ;

    ElectronMaker electronmaker                   ; 
    MuonMaker muonmaker                           ; 

    JetMaker jetAK4maker                          ; 
    JetMaker jetAK8maker                          ; 
    JetMaker jetHTaggedmaker                      ; 
    JetMaker jetTopTaggedmaker                    ; 
    JetMaker jetAntiTopTaggedmaker                    ; 
    JetMaker jetAntiHTaggedmaker                  ; 

    double leadingJetPtMin_                       ; 
    double leadingJetPrunedMassMin_               ; 
    double HTMin_                                 ; 
    bool   storePreselEvts_                       ; 
    bool   doPreselOnly_                          ; 

    edm::Service<TFileService> fs                 ; 
    std::map<std::string, TH1D*> h1_              ; 
    std::map<std::string, TH2D*> h2_              ; 

    TtHEventInfoBranches selectedevt_; 
    TtHJetInfoBranches jets_ ; 
    TTree* tree_ ; 

    std::unique_ptr<BTagSFUtils> btagsfutils_ ; 

};

using namespace std; 

VLQAna::VLQAna(const edm::ParameterSet& iConfig) :
  t_runno                 (consumes<int>            (iConfig.getParameter<edm::InputTag>("runno"))),   
  t_lumisec               (consumes<int>            (iConfig.getParameter<edm::InputTag>("lumisec"))),   
  t_evtno                 (consumes<int>            (iConfig.getParameter<edm::InputTag>("evtno"))),   
  t_isData                (consumes<bool>           (iConfig.getParameter<edm::InputTag>("isData"))),
  t_hltdecision           (consumes<bool>           (iConfig.getParameter<edm::InputTag>("hltdecision"))),
  t_evttype               (consumes<string>         (iConfig.getParameter<edm::InputTag>("evttype"))),
  t_evtwtGen              (consumes<double>         (iConfig.getParameter<edm::InputTag>("evtwtGen"))),
  t_evtwtPV               (consumes<double>         (iConfig.getParameter<edm::InputTag>("evtwtPV"))),
  t_evtwtPVLow            (consumes<double>         (iConfig.getParameter<edm::InputTag>("evtwtPVLow"))),
  t_evtwtPVHigh           (consumes<double>         (iConfig.getParameter<edm::InputTag>("evtwtPVHigh"))),
  t_npv                   (consumes<unsigned>       (iConfig.getParameter<edm::InputTag>("npv"))),
  t_npuTrue               (consumes<int>            (iConfig.getParameter<edm::InputTag>("npuTrue"))),
  t_htHat                 (consumes<double>         (iConfig.getParameter<edm::InputTag>("htHat"))),
  t_lhewtids              (consumes<vector<int>>    (iConfig.getParameter<edm::InputTag>("lhewtids"))),
  t_lhewts                (consumes<vector<double>> (iConfig.getParameter<edm::InputTag>("lhewts"))),
  electronmaker           (iConfig.getParameter<edm::ParameterSet>("elselParams"),consumesCollector()),
  muonmaker               (iConfig.getParameter<edm::ParameterSet>("muselParams"),consumesCollector()),
  jetAK4maker             (iConfig.getParameter<edm::ParameterSet>("jetAK4selParams"), consumesCollector()), 
  jetAK8maker             (iConfig.getParameter<edm::ParameterSet>("jetAK8selParams"), consumesCollector()), 
  jetHTaggedmaker         (iConfig.getParameter<edm::ParameterSet>("jetHTaggedselParams"), consumesCollector()), 
  jetTopTaggedmaker       (iConfig.getParameter<edm::ParameterSet>("jetTopTaggedselParams"), consumesCollector()),  
  jetAntiTopTaggedmaker   (iConfig.getParameter<edm::ParameterSet>("jetAntiTopTaggedselParams"), consumesCollector()),  
  jetAntiHTaggedmaker     (iConfig.getParameter<edm::ParameterSet>("jetAntiHTaggedselParams"), consumesCollector()), 
  leadingJetPtMin_        (iConfig.getParameter<double>           ("leadingJetPtMin")), 
  leadingJetPrunedMassMin_(iConfig.getParameter<double>           ("leadingJetPrunedMassMin")), 
  HTMin_                  (iConfig.getParameter<double>           ("HTMin")), 
  storePreselEvts_        (iConfig.getParameter<bool>             ("storePreselEvts")),
  doPreselOnly_           (iConfig.getParameter<bool>             ("doPreselOnly")), 
  btagsfutils_            (new BTagSFUtils())
{

}

VLQAna::~VLQAna() {
}

bool VLQAna::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<int>            h_runno         ; evt.getByToken(t_runno      , h_runno      ) ; 
  Handle<int>            h_lumisec       ; evt.getByToken(t_lumisec    , h_lumisec    ) ; 
  Handle<int>            h_evtno         ; evt.getByToken(t_evtno      , h_evtno      ) ; 
  Handle<bool>           h_isData        ; evt.getByToken(t_isData     , h_isData     ) ; 
  Handle<bool>           h_hltdecision   ; evt.getByToken(t_hltdecision, h_hltdecision) ; 
  Handle<string>         h_evttype       ; evt.getByToken(t_evttype    , h_evttype    ) ; 
  Handle<double>         h_evtwtGen      ; evt.getByToken(t_evtwtGen   , h_evtwtGen   ) ; 
  Handle<double>         h_evtwtPV       ; evt.getByToken(t_evtwtPV    , h_evtwtPV    ) ; 
  Handle<double>         h_evtwtPVLow    ; evt.getByToken(t_evtwtPVLow , h_evtwtPVLow ) ; 
  Handle<double>         h_evtwtPVHigh   ; evt.getByToken(t_evtwtPVHigh, h_evtwtPVHigh) ; 
  Handle<unsigned>       h_npv           ; evt.getByToken(t_npv        , h_npv        ) ; 
  Handle<int>            h_npuTrue       ; evt.getByToken(t_npuTrue    , h_npuTrue    ) ; 
  Handle<double>         h_htHat         ; evt.getByToken(t_htHat      , h_htHat      ) ; 
  Handle<vector<int>>    h_lhewtids      ; evt.getByToken(t_lhewtids   , h_lhewtids   ) ; 
  Handle<vector<double>> h_lhewts        ; evt.getByToken(t_lhewts     , h_lhewts     ) ; 

  const int runno(*h_runno.product()) ; 
  const int lumisec(*h_lumisec.product()) ; 
  const int evtno(*h_evtno.product()) ; 
  const bool isData(*h_isData.product()) ; 
  const bool hltdecision(*h_hltdecision.product()) ; 
  double evtwt((*h_evtwtGen.product()) * (*h_evtwtPV.product())) ; 
  int npv(*h_npv.product()) ; 
  const double htHat((*h_htHat.product())) ; 

  h1_["cutflow"] -> Fill(1, evtwt) ; 

  if ( isData && !hltdecision ) return false; 

  h1_["cutflow"] -> Fill(2, evtwt) ; 

  vlq::JetCollection goodAK4Jets;
  jetAK4maker(evt, goodAK4Jets) ;
  //// Event pre-selection
  if (goodAK4Jets.size() < 4) return false ; 
  h1_["cutflow"] -> Fill(3, evtwt) ; 

  HT htak4(goodAK4Jets) ; 
  //// Event pre-selection
  if (htak4.getHT() < HTMin_) return false ; 
  h1_["cutflow"] -> Fill(4, evtwt) ; 

  vlq::JetCollection goodAK8Jets; 
  jetAK8maker(evt, goodAK8Jets); 

  //// Event pre-selection
  if (goodAK8Jets.size() < 2) return false ; 
  h1_["cutflow"] -> Fill(5, evtwt) ; 

  //// Event pre-selection
  if (goodAK8Jets.at(0).getPt() < leadingJetPtMin_ 
      && goodAK8Jets.at(0).getPrunedMass() < leadingJetPrunedMassMin_) return false ; 
  h1_["cutflow"] -> Fill(6, evtwt) ; 

  h1_["npv_noreweight"] -> Fill(*h_npv.product(), *h_evtwtGen.product()); 
  h1_["npv"] -> Fill(*h_npv.product(), evtwt); 

  if ( storePreselEvts_ || doPreselOnly_ ) { 
    h1_["Presel_HT"] -> Fill(htak4.getHT(), evtwt) ; 

    if (goodAK8Jets.size() > 1) 
      h1_["Presel_Mjj"] -> Fill( (goodAK8Jets.at(0).getP4() + goodAK8Jets.at(1).getP4()).Mag() ) ; 
    else 
      h1_["Presel_Mjj"] -> Fill( -1 ); 

    h1_["Presel_ptAK4_0"] -> Fill (goodAK4Jets.at(0).getPt(),evtwt) ;  
    h1_["Presel_ptAK4_1"] -> Fill (goodAK4Jets.at(1).getPt(),evtwt) ;  
    h1_["Presel_ptAK4_2"] -> Fill (goodAK4Jets.at(2).getPt(),evtwt) ;  
    h1_["Presel_ptAK4_3"] -> Fill (goodAK4Jets.at(3).getPt(),evtwt) ;  
    h1_["Presel_ptAK8_0"] -> Fill (goodAK8Jets.at(0).getPt(),evtwt) ;  

    h1_["Presel_etaAK4_0"] -> Fill (goodAK4Jets.at(0).getEta(),evtwt) ;  
    h1_["Presel_etaAK4_1"] -> Fill (goodAK4Jets.at(1).getEta(),evtwt) ;  
    h1_["Presel_etaAK4_2"] -> Fill (goodAK4Jets.at(2).getEta(),evtwt) ;  
    h1_["Presel_etaAK4_3"] -> Fill (goodAK4Jets.at(3).getEta(),evtwt) ;  
    h1_["Presel_etaAK8_0"] -> Fill (goodAK8Jets.at(0).getEta(),evtwt) ;  

    h1_["Presel_NAK4"] -> Fill(goodAK4Jets.size(), evtwt) ; 
    h1_["Presel_NAK8"] -> Fill(goodAK8Jets.size(), evtwt) ; 

  }

  vlq::JetCollection goodHTaggedJets, goodTopTaggedJets, antiHTaggedJets, antiTopTaggedJets ; 

  jetHTaggedmaker(evt, goodHTaggedJets);
  jetTopTaggedmaker(evt, goodTopTaggedJets);
  jetAntiHTaggedmaker(evt, antiHTaggedJets);
  jetAntiTopTaggedmaker(evt, antiTopTaggedJets);

  unsigned nAK4      (goodAK4Jets.size());
  unsigned nAK8      (goodAK8Jets.size());
  unsigned nHiggs    (goodHTaggedJets.size()); 
  unsigned nTop      (goodTopTaggedJets.size()); 
  unsigned nAntiHiggs(antiHTaggedJets.size()); 
  unsigned nAntiTop(antiTopTaggedJets.size()); 

  //// Create 4 regions of the ABCD method according the the scheme below 
  //// | A: Anti-H Anti-top | B: Anti-H Good top | 
  //// | C: Good H Anti-top | D: Good H Good top | 

  bool isRegionA(false), isRegionB(false), isRegionC(false), isRegionD(false), isRegionNotABCD(false); 

  std::unique_ptr<vlq::Jet> theHiggs ;  
  std::unique_ptr<vlq::Jet> theTop ;
  std::unique_ptr<vlq::Jet> theAntiHiggs ; 
  std::unique_ptr<vlq::Jet> theAntiTop ;
  if (nTop > 0) {

    theTop = std::unique_ptr<vlq::Jet>(new vlq::Jet(goodTopTaggedJets.at(0))) ; 

    if (nHiggs > 0) {
      for (vlq::Jet& hjet : goodHTaggedJets) {
        double dphi = abs((theTop->getP4()).DeltaPhi(hjet.getP4())) ;  
        if (dphi > 2.0) {
          theHiggs = std::unique_ptr<vlq::Jet>(new vlq::Jet(hjet)) ; 
          isRegionD = true ; 
          break ; 
        }
      }
    } //// isRegionD 
    else isRegionD = false ; 

    if (isRegionD == false && nAntiHiggs > 0) {
      for (vlq::Jet& antihjet : antiHTaggedJets) {
        double dphi = abs((theTop->getP4()).DeltaPhi(antihjet.getP4())) ;  
        if (dphi > 2.0) {
          theAntiHiggs = std::unique_ptr<vlq::Jet>(new vlq::Jet(antihjet)) ; 
          isRegionB = true ; 
          break ; 
        }
      }
    } //// isRegionB
    else isRegionNotABCD = true ;

  }
  else if (nAntiTop > 0){

    theAntiTop = std::unique_ptr<vlq::Jet>(new vlq::Jet(antiTopTaggedJets.at(0))) ;

    if (nHiggs > 0) {
      for (vlq::Jet& hjet : goodHTaggedJets) {
        double dphi = abs((theAntiTop->getP4()).DeltaPhi(hjet.getP4())) ;  
        if (dphi > 2.0) {
          theHiggs = std::unique_ptr<vlq::Jet>(new vlq::Jet(hjet)) ; 
          isRegionC = true;
          break ;
        }
      }
    } //// isRegionC 
    else if (nAntiHiggs > 0) {
      for (vlq::Jet& antihjet : antiHTaggedJets) {
        double dphi = abs((theAntiTop->getP4()).DeltaPhi(antihjet.getP4())) ;  
        if (dphi > 2.0) {
          theAntiHiggs = std::unique_ptr<vlq::Jet>(new vlq::Jet(antihjet)) ; 
          isRegionA = true;
          break ;
        }
      }
    } //// isRegionA 
    else isRegionNotABCD = true ;

  }

  if ( !storePreselEvts_ && !isRegionA && !isRegionB && !isRegionC && !isRegionD ) return false ; 

  bool isOneOfABCD(isRegionA ^ isRegionB ^ isRegionC ^ isRegionD) ;
  if ( !storePreselEvts_ && isOneOfABCD == false) edm::LogInfo("ERROR ABCD") << ">>>> Check ABCD logic: Only one of A, B C, D should be true\n" ; 

  if ( theHiggs != nullptr ) h1_["cutflow"] -> Fill(7, evtwt) ; 
  if ( theTop != nullptr ) h1_["cutflow"] -> Fill(8, evtwt) ; 
  if ( isRegionA ) h1_["cutflow"] -> Fill(9, evtwt) ; 
  if ( isRegionB ) h1_["cutflow"] -> Fill(10, evtwt) ; 
  if ( isRegionC ) h1_["cutflow"] -> Fill(11, evtwt) ; 
  if ( isRegionD ) h1_["cutflow"] -> Fill(12, evtwt) ; 

  TLorentzVector p4_tprime, p4_TprimeDummy, p4_TprimeDummyA, p4_TprimeDummyC ; 
  double Mtprime(0), Mtprime_corr(0), MtprimeDummy(0), MtprimeDummy_corr(0), MtprimeDummyA(0), MtprimeDummyC(0) ; 

  if (isRegionD) {
    p4_tprime = theTop->getP4() + theHiggs->getP4() ; 
    Mtprime = p4_tprime.Mag();
    Mtprime_corr = Mtprime - theTop->getSoftDropMass() - theHiggs->getPrunedMass() + 172.5 + 125. ; 
  }
  else if (isRegionB) {
    p4_TprimeDummy = theTop->getP4() + theAntiHiggs->getP4() ; 
    MtprimeDummy = p4_TprimeDummy.Mag() ; 
    MtprimeDummy_corr = MtprimeDummy - theTop->getSoftDropMass() - theAntiHiggs->getPrunedMass() + 172.5 + 125. ; 
  }
  else if (isRegionA) {
    p4_TprimeDummyA = theAntiTop->getP4() + theAntiHiggs->getP4() ; 
    MtprimeDummyA = p4_TprimeDummyA.Mag() ; 
//    MtprimeDummyA_corr = MtprimeDummyA - theAntiTop->getSoftDropMass() - theAntiHiggs->getPrunedMass() + 172.5 + 125. ; 
  }
  else if (isRegionC) {
    p4_TprimeDummyC = theAntiTop->getP4() + theHiggs->getP4() ; 
    MtprimeDummyC = p4_TprimeDummyC.Mag() ; 
//    MtprimeDummyC_corr = MtprimeDummyC - theAntiTop->getSoftDropMass() - theHiggs->getPrunedMass() + 172.5 + 125. ; 
  }
 
  double evtwtHT(1), evtwtHTUp(1), evtwtHTDown(1);
  double toptagsf(1), toptagsfUp(1), toptagsfDown(1);
  double btagsf(1), btagsf_bcUp(1), btagsf_bcDown(1), btagsf_lUp(1), btagsf_lDown(1) ; 

  if ( !isData ) { 

    evtwtHT = getHTReweightingSF(htak4.getHT(), 0) ; 
    evtwtHTUp = getHTReweightingSF(htak4.getHT(), 1) ; 
    evtwtHTDown = getHTReweightingSF(htak4.getHT(), -1) ; 

    std::vector<double> sjcsvs(0) ; 
    std::vector<double> sjpts (0) ; 
    std::vector<double> sjetas(0) ; 
    std::vector<int> sjfls (0) ; 

    //// This needs to be revisied when the top tag WP is updated
    if ( theTop != nullptr ){
      if ( theTop->getPt() >= 400. && theTop->getPt() < 550.) {
        toptagsf = 0.88 ;
        toptagsfUp = 0.88+0.13 ; 
        toptagsfDown = 0.88-0.13;
      }
      else if ( theTop != nullptr && theTop->getPt() >= 550.)  {
        toptagsf = 0.94 ;
        toptagsfUp = 0.94+0.29 ; 
        toptagsfDown = 0.94-0.29;
      }
      sjcsvs.push_back( (theTop->getCSVSubjet0()) ); 
      sjpts .push_back( (theTop->getPtSubjet0()) ); 
      sjetas.push_back( (theTop->getEtaSubjet0()) ); 
      sjfls .push_back( (theTop->getHadronFlavourSubjet0()) ); 
      sjcsvs.push_back( (theTop->getCSVSubjet1()) ); 
      sjpts .push_back( (theTop->getPtSubjet1()) ); 
      sjetas.push_back( (theTop->getEtaSubjet1()) ); 
      sjfls .push_back( (theTop->getHadronFlavourSubjet1()) ); 
    }
    if ( theAntiTop != nullptr ){
      sjcsvs.push_back( (theAntiTop->getCSVSubjet0()) ); 
      sjpts .push_back( (theAntiTop->getPtSubjet0()) ); 
      sjetas.push_back( (theAntiTop->getEtaSubjet0()) ); 
      sjfls .push_back( (theAntiTop->getHadronFlavourSubjet0()) ); 
      sjcsvs.push_back( (theAntiTop->getCSVSubjet1()) ); 
      sjpts .push_back( (theAntiTop->getPtSubjet1()) ); 
      sjetas.push_back( (theAntiTop->getEtaSubjet1()) ); 
      sjfls .push_back( (theAntiTop->getHadronFlavourSubjet1()) ); 
    }

    if ( theHiggs != nullptr ) {
      sjcsvs.push_back( (theHiggs->getCSVSubjet0()) ); 
      sjpts .push_back( (theHiggs->getPtSubjet0()) ); 
      sjetas.push_back( (theHiggs->getEtaSubjet0()) ); 
      sjfls .push_back( (theHiggs->getHadronFlavourSubjet0()) ); 
      sjcsvs.push_back( (theHiggs->getCSVSubjet1()) ); 
      sjpts .push_back( (theHiggs->getPtSubjet1()) ); 
      sjetas.push_back( (theHiggs->getEtaSubjet1()) ); 
      sjfls .push_back( (theHiggs->getHadronFlavourSubjet1()) ); 
    }
    else if ( theAntiHiggs != nullptr ) {
      sjcsvs.push_back( (theAntiHiggs->getCSVSubjet0()) ); 
      sjpts .push_back( (theAntiHiggs->getPtSubjet0()) ); 
      sjetas.push_back( (theAntiHiggs->getPtSubjet0()) ); 
      sjfls .push_back( (theAntiHiggs->getHadronFlavourSubjet0()) ); 
      sjcsvs.push_back( (theAntiHiggs->getCSVSubjet1()) ); 
      sjpts .push_back( (theAntiHiggs->getPtSubjet1()) ); 
      sjetas.push_back( (theAntiHiggs->getPtSubjet1()) ) ; 
      sjfls .push_back( (theAntiHiggs->getHadronFlavourSubjet1()) ); 
    }

    ////// Get btag SFs
    btagsfutils_-> getBTagSFs(sjcsvs, sjpts, sjetas, sjfls, jetHTaggedmaker.idxsjHighestCSVMin_, btagsf, btagsf_bcUp, btagsf_bcDown, btagsf_lUp, btagsf_lDown) ; 

  } //// if !isData

  if ( isRegionB ) {
    h1_["RegB_HT"] -> Fill(htak4.getHT(), evtwt) ; 
    h1_["RegB_HT_wts"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHT*btagsf) ; 
    h1_["RegB_HT_btagsf_bcUp"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHT*btagsf_bcUp) ; 
    h1_["RegB_HT_btagsf_bcDown"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHT*btagsf_bcDown) ; 
    h1_["RegB_HT_btagsf_lUp"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHT*btagsf_lUp) ; 
    h1_["RegB_HT_btagsf_lDown"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHT*btagsf_lDown) ; 
    h1_["RegB_HT_toptagsfUp"] -> Fill(htak4.getHT(), evtwt*toptagsfUp*evtwtHT*btagsf) ; 
    h1_["RegB_HT_toptagsfDown"] -> Fill(htak4.getHT(), evtwt*toptagsfDown*evtwtHT*btagsf) ; 
    h1_["RegB_HT_htwtUp"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHTUp*btagsf) ; 
    h1_["RegB_HT_htwtDown"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHTDown*btagsf) ; 

    h1_["RegB_mtprime"] -> Fill(MtprimeDummy, evtwt) ; 
    h1_["RegB_mtprime_wts"] -> Fill(MtprimeDummy, evtwt*toptagsf*evtwtHT*btagsf) ; 
    h1_["RegB_mtprime_btagsf_bcUp"] -> Fill(MtprimeDummy, evtwt*toptagsf*evtwtHT*btagsf_bcUp) ; 
    h1_["RegB_mtprime_btagsf_bcDown"] -> Fill(MtprimeDummy, evtwt*toptagsf*evtwtHT*btagsf_bcDown) ; 
    h1_["RegB_mtprime_btagsf_lUp"] -> Fill(MtprimeDummy, evtwt*toptagsf*evtwtHT*btagsf_lUp) ; 
    h1_["RegB_mtprime_btagsf_lDown"] -> Fill(MtprimeDummy, evtwt*toptagsf*evtwtHT*btagsf_lDown) ; 
    h1_["RegB_mtprime_toptagsfUp"] -> Fill(MtprimeDummy, evtwt*toptagsfUp*evtwtHT*btagsf) ; 
    h1_["RegB_mtprime_toptagsfDown"] -> Fill(MtprimeDummy, evtwt*toptagsfDown*evtwtHT*btagsf) ; 
    h1_["RegB_mtprime_htwtUp"] -> Fill(MtprimeDummy, evtwt*toptagsf*evtwtHTUp*btagsf) ; 
    h1_["RegB_mtprime_htwtDown"] -> Fill(MtprimeDummy, evtwt*toptagsf*evtwtHTDown*btagsf) ; 

    h1_["RegB_mtprime_corr"] -> Fill(MtprimeDummy_corr, evtwt) ; 
    h1_["RegB_mtprime_corr_wts"] -> Fill(MtprimeDummy_corr, evtwt*toptagsf*evtwtHT*btagsf) ; 
    h1_["RegB_mtprime_corr_btagsf_bcUp"] -> Fill(MtprimeDummy_corr, evtwt*toptagsf*evtwtHT*btagsf_bcUp) ; 
    h1_["RegB_mtprime_corr_btagsf_bcDown"] -> Fill(MtprimeDummy_corr, evtwt*toptagsf*evtwtHT*btagsf_bcDown) ; 
    h1_["RegB_mtprime_corr_btagsf_lUp"] -> Fill(MtprimeDummy_corr, evtwt*toptagsf*evtwtHT*btagsf_lUp) ; 
    h1_["RegB_mtprime_corr_btagsf_lDown"] -> Fill(MtprimeDummy_corr, evtwt*toptagsf*evtwtHT*btagsf_lDown) ; 
    h1_["RegB_mtprime_corr_toptagsfUp"] -> Fill(MtprimeDummy_corr, evtwt*toptagsfUp*evtwtHT*btagsf) ; 
    h1_["RegB_mtprime_corr_toptagsfDown"] -> Fill(MtprimeDummy_corr, evtwt*toptagsfDown*evtwtHT*btagsf) ; 
    h1_["RegB_mtprime_corr_htwtUp"] -> Fill(MtprimeDummy_corr, evtwt*toptagsf*evtwtHTUp*btagsf) ; 
    h1_["RegB_mtprime_corr_htwtDown"] -> Fill(MtprimeDummy_corr, evtwt*toptagsf*evtwtHTDown*btagsf) ; 
  }

  if ( isRegionD ) {
    h1_["RegD_HT"] -> Fill(htak4.getHT(), evtwt) ; 
    h1_["RegD_HT_wts"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHT*btagsf) ; 
    h1_["RegD_HT_btagsf_bcUp"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHT*btagsf_bcUp) ; 
    h1_["RegD_HT_btagsf_bcDown"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHT*btagsf_bcDown) ; 
    h1_["RegD_HT_btagsf_lUp"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHT*btagsf_lUp) ; 
    h1_["RegD_HT_btagsf_lDown"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHT*btagsf_lDown) ; 
    h1_["RegD_HT_toptagsfUp"] -> Fill(htak4.getHT(), evtwt*toptagsfUp*evtwtHT*btagsf) ; 
    h1_["RegD_HT_toptagsfDown"] -> Fill(htak4.getHT(), evtwt*toptagsfDown*evtwtHT*btagsf) ; 
    h1_["RegD_HT_htwtUp"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHTUp*btagsf) ; 
    h1_["RegD_HT_htwtDown"] -> Fill(htak4.getHT(), evtwt*toptagsf*evtwtHTDown*btagsf) ; 

    h1_["RegD_mtprime"] -> Fill(Mtprime, evtwt) ; 
    h1_["RegD_mtprime_wts"] -> Fill(Mtprime, evtwt*toptagsf*evtwtHT*btagsf) ; 
    h1_["RegD_mtprime_btagsf_bcUp"] -> Fill(Mtprime, evtwt*toptagsf*evtwtHT*btagsf_bcUp) ; 
    h1_["RegD_mtprime_btagsf_bcDown"] -> Fill(Mtprime, evtwt*toptagsf*evtwtHT*btagsf_bcDown) ; 
    h1_["RegD_mtprime_btagsf_lUp"] -> Fill(Mtprime, evtwt*toptagsf*evtwtHT*btagsf_lUp) ; 
    h1_["RegD_mtprime_btagsf_lDown"] -> Fill(Mtprime, evtwt*toptagsf*evtwtHT*btagsf_lDown) ; 
    h1_["RegD_mtprime_toptagsfUp"] -> Fill(Mtprime, evtwt*toptagsfUp*evtwtHT*btagsf) ; 
    h1_["RegD_mtprime_toptagsfDown"] -> Fill(Mtprime, evtwt*toptagsfDown*evtwtHT*btagsf) ; 
    h1_["RegD_mtprime_htwtUp"] -> Fill(Mtprime, evtwt*toptagsf*evtwtHTUp*btagsf) ; 
    h1_["RegD_mtprime_htwtDown"] -> Fill(Mtprime, evtwt*toptagsf*evtwtHTDown*btagsf) ; 

    h1_["RegD_mtprime_corr"] -> Fill(Mtprime_corr, evtwt) ; 
    h1_["RegD_mtprime_corr_wts"] -> Fill(Mtprime_corr, evtwt*toptagsf*evtwtHT*btagsf) ; 
    h1_["RegD_mtprime_corr_btagsf_bcUp"] -> Fill(Mtprime_corr, evtwt*toptagsf*evtwtHT*btagsf_bcUp) ; 
    h1_["RegD_mtprime_corr_btagsf_bcDown"] -> Fill(Mtprime_corr, evtwt*toptagsf*evtwtHT*btagsf_bcDown) ; 
    h1_["RegD_mtprime_corr_btagsf_lUp"] -> Fill(Mtprime_corr, evtwt*toptagsf*evtwtHT*btagsf_lUp) ; 
    h1_["RegD_mtprime_corr_btagsf_lDown"] -> Fill(Mtprime_corr, evtwt*toptagsf*evtwtHT*btagsf_lDown) ; 
    h1_["RegD_mtprime_corr_toptagsfUp"] -> Fill(Mtprime_corr, evtwt*toptagsfUp*evtwtHT*btagsf) ; 
    h1_["RegD_mtprime_corr_toptagsfDown"] -> Fill(Mtprime_corr, evtwt*toptagsfDown*evtwtHT*btagsf) ; 
    h1_["RegD_mtprime_corr_htwtUp"] -> Fill(Mtprime_corr, evtwt*toptagsf*evtwtHTUp*btagsf) ; 
    h1_["RegD_mtprime_corr_htwtDown"] -> Fill(Mtprime_corr, evtwt*toptagsf*evtwtHTDown*btagsf) ; 
  }

  std::vector<std::pair<int, double>> lhe_ids_wts;
  for (auto idx : index(*h_lhewtids.product()) ) {
    int id = (*h_lhewtids.product()).at(idx.first) ; 
    double wt = (*h_lhewts.product()).at(idx.first) ; 
    lhe_ids_wts.push_back(std::make_pair(id, wt)) ; 
  }

  selectedevt_.runno_ = int(runno);
  selectedevt_.lumisec_ = int(lumisec);
  selectedevt_.evtno_ = int(evtno);
  selectedevt_.EvtWeight_ = double(*h_evtwtGen.product());
  selectedevt_.EvtWtPV_ = double(*h_evtwtPV.product()) ; 
  selectedevt_.EvtWtPVLow_ = double(*h_evtwtPVLow.product()) ; 
  selectedevt_.EvtWtPVHigh_ = double(*h_evtwtPVHigh.product()) ; 
  selectedevt_.EvtWtHT_ = evtwtHT;
  selectedevt_.EvtWtHTUp_ = evtwtHTUp;
  selectedevt_.EvtWtHTDown_ = evtwtHTDown;
  selectedevt_.npv_ = npv ; 
  selectedevt_.npuTrue_ = int(*h_npuTrue.product()) ; 
  selectedevt_.toptagsf_ = toptagsf ; 
  selectedevt_.toptagsf_Up_ = toptagsfUp ; 
  selectedevt_.toptagsf_Down_ = toptagsfDown ; 
  selectedevt_.btagsf_ = btagsf ; 
  selectedevt_.btagsf_bcUp_ = btagsf_bcUp ; 
  selectedevt_.btagsf_bcDown_ = btagsf_bcDown ; 
  selectedevt_.btagsf_lUp_ = btagsf_lUp ; 
  selectedevt_.btagsf_lDown_ = btagsf_lDown ; 
  selectedevt_.mtprime_ = Mtprime ;
  selectedevt_.mtprimeDummy_ = MtprimeDummy ;
  selectedevt_.mtprimeDummyA_ = MtprimeDummyA ;
  selectedevt_.mtprimeDummyC_ = MtprimeDummyC ;
  selectedevt_.ht_ = htak4.getHT();
  selectedevt_.nAK4_ = int(nAK4);
  selectedevt_.nAK8_ = int(nAK8);
  selectedevt_.isRegionA_ = isRegionA ; 
  selectedevt_.isRegionB_ = isRegionB ; 
  selectedevt_.isRegionC_ = isRegionC ; 
  selectedevt_.isRegionD_ = isRegionD ; 
  selectedevt_.isRegionNotABCD_ = isRegionNotABCD;
  selectedevt_.lhewts_ = lhe_ids_wts;
  selectedevt_.htHat_ = htHat ; 

  jets_.idxAK4             .clear() ; jets_.idxAK4             .reserve(goodAK4Jets.size()) ;   
  jets_.ptAK4              .clear() ; jets_.ptAK4              .reserve(goodAK4Jets.size()) ; 
  jets_.etaAK4             .clear() ; jets_.etaAK4             .reserve(goodAK4Jets.size()) ;   
  jets_.phiAK4             .clear() ; jets_.phiAK4             .reserve(goodAK4Jets.size()) ;   
  jets_.MAK4               .clear() ; jets_.MAK4               .reserve(goodAK4Jets.size()) ;   
  jets_.csvAK4             .clear() ; jets_.csvAK4             .reserve(goodAK4Jets.size()) ;   
  jets_.partonFlavourAK4   .clear() ; jets_.partonFlavourAK4   .reserve(goodAK4Jets.size()) ;   
  jets_.hadronFlavourAK4   .clear() ; jets_.hadronFlavourAK4   .reserve(goodAK4Jets.size()) ;   

  for (vlq::Jet jet : goodAK4Jets) {
    jets_.idxAK4             .push_back(jet.getIndex()) ; 
    jets_.ptAK4              .push_back(jet.getPt()) ; 
    jets_.etaAK4             .push_back(jet.getEta()) ; 
    jets_.phiAK4             .push_back(jet.getPhi()) ; 
    jets_.MAK4               .push_back(jet.getMass()) ; 
    jets_.csvAK4             .push_back(jet.getCSV());
    jets_.partonFlavourAK4   .push_back(jet.getPartonFlavour());
    jets_.hadronFlavourAK4   .push_back(jet.getHadronFlavour());
  }

  jets_.idxAK8             .clear() ; jets_.idxAK8             .reserve(goodAK8Jets.size()) ;  
  jets_.ptAK8              .clear() ; jets_.ptAK8              .reserve(goodAK8Jets.size()) ;  
  jets_.etaAK8             .clear() ; jets_.etaAK8             .reserve(goodAK8Jets.size()) ;  
  jets_.phiAK8             .clear() ; jets_.phiAK8             .reserve(goodAK8Jets.size()) ;  
  jets_.MAK8               .clear() ; jets_.MAK8               .reserve(goodAK8Jets.size()) ;  
  jets_.SoftDropMassAK8    .clear() ; jets_.SoftDropMassAK8    .reserve(goodAK8Jets.size()) ;  
  jets_.PrunedMassAK8      .clear() ; jets_.PrunedMassAK8      .reserve(goodAK8Jets.size()) ;  
  jets_.tau1AK8            .clear() ; jets_.tau1AK8            .reserve(goodAK8Jets.size()) ;  
  jets_.tau2AK8            .clear() ; jets_.tau2AK8            .reserve(goodAK8Jets.size()) ;  
  jets_.tau3AK8            .clear() ; jets_.tau3AK8            .reserve(goodAK8Jets.size()) ;  
  jets_.csvAK8             .clear() ; jets_.csvAK8             .reserve(goodAK8Jets.size()) ;  
  jets_.partonFlavourAK8   .clear() ; jets_.partonFlavourAK8   .reserve(goodAK8Jets.size()) ;  
  jets_.hadronFlavourAK8   .clear() ; jets_.hadronFlavourAK8   .reserve(goodAK8Jets.size()) ;  
  jets_.doubleBAK8         .clear() ; jets_.doubleBAK8         .reserve(goodAK8Jets.size()) ;  
  jets_.sj0CSVAK8          .clear() ; jets_.sj0CSVAK8          .reserve(goodAK8Jets.size()) ;  
  jets_.sj1CSVAK8          .clear() ; jets_.sj1CSVAK8          .reserve(goodAK8Jets.size()) ;  
  jets_.hadronFlavourSJ0AK8.clear() ; jets_.hadronFlavourSJ0AK8.reserve(goodAK8Jets.size()) ;  
  jets_.hadronFlavourSJ1AK8.clear() ; jets_.hadronFlavourSJ1AK8.reserve(goodAK8Jets.size()) ;  

  for (vlq::Jet jet : goodAK8Jets) {
    jets_.idxAK8             .push_back(jet.getIndex()) ; 
    jets_.ptAK8              .push_back(jet.getPt()) ; 
    jets_.etaAK8             .push_back(jet.getEta()) ; 
    jets_.phiAK8             .push_back(jet.getPhi()) ; 
    jets_.MAK8               .push_back(jet.getMass()) ; 
    jets_.SoftDropMassAK8    .push_back(jet.getSoftDropMass()) ;
    jets_.PrunedMassAK8      .push_back(jet.getPrunedMass()) ;
    jets_.tau1AK8            .push_back(jet.getTau1()) ;
    jets_.tau2AK8            .push_back(jet.getTau2()) ;
    jets_.tau3AK8            .push_back(jet.getTau3()) ;
    jets_.csvAK8             .push_back(jet.getCSV());
    jets_.partonFlavourAK8   .push_back(jet.getPartonFlavour());
    jets_.hadronFlavourAK8   .push_back(jet.getHadronFlavour());
    jets_.doubleBAK8         .push_back(jet.getDoubleBAK8()) ;
    jets_.sj0CSVAK8          .push_back(jet.getCSVSubjet0()) ;
    jets_.sj1CSVAK8          .push_back(jet.getCSVSubjet1()) ;
    jets_.hadronFlavourSJ0AK8.push_back(jet.getHadronFlavourSubjet0()) ;
    jets_.hadronFlavourSJ1AK8.push_back(jet.getHadronFlavourSubjet1()) ;
  }

  jets_.idxHTagged             .clear() ; jets_.idxHTagged             .reserve(goodHTaggedJets.size()) ;   
  jets_.ptHTagged              .clear() ; jets_.ptHTagged              .reserve(goodHTaggedJets.size()) ;   
  jets_.etaHTagged             .clear() ; jets_.etaHTagged             .reserve(goodHTaggedJets.size()) ;   
  jets_.phiHTagged             .clear() ; jets_.phiHTagged             .reserve(goodHTaggedJets.size()) ;   
  jets_.MHTagged               .clear() ; jets_.MHTagged               .reserve(goodHTaggedJets.size()) ;   
  jets_.SoftDropMassHTagged    .clear() ; jets_.SoftDropMassHTagged    .reserve(goodHTaggedJets.size()) ;   
  jets_.PrunedMassHTagged      .clear() ; jets_.PrunedMassHTagged      .reserve(goodHTaggedJets.size()) ;   
  jets_.tau1HTagged            .clear() ; jets_.tau1HTagged            .reserve(goodHTaggedJets.size()) ;   
  jets_.tau2HTagged            .clear() ; jets_.tau2HTagged            .reserve(goodHTaggedJets.size()) ;   
  jets_.tau3HTagged            .clear() ; jets_.tau3HTagged            .reserve(goodHTaggedJets.size()) ;   
  jets_.csvHTagged             .clear() ; jets_.csvHTagged             .reserve(goodHTaggedJets.size()) ;   
  jets_.partonFlavourHTagged   .clear() ; jets_.partonFlavourHTagged   .reserve(goodHTaggedJets.size()) ;   
  jets_.hadronFlavourHTagged   .clear() ; jets_.hadronFlavourHTagged   .reserve(goodHTaggedJets.size()) ;   
  jets_.doubleBHTagged         .clear() ; jets_.doubleBHTagged         .reserve(goodHTaggedJets.size()) ;   
  jets_.sj0CSVHTagged          .clear() ; jets_.sj0CSVHTagged          .reserve(goodHTaggedJets.size()) ;   
  jets_.sj1CSVHTagged          .clear() ; jets_.sj1CSVHTagged          .reserve(goodHTaggedJets.size()) ;   
  jets_.hadronFlavourSJ0HTagged.clear() ; jets_.hadronFlavourSJ0HTagged.reserve(goodHTaggedJets.size()) ;   
  jets_.hadronFlavourSJ1HTagged.clear() ; jets_.hadronFlavourSJ1HTagged.reserve(goodHTaggedJets.size()) ;   

  for (vlq::Jet jet : goodHTaggedJets) {
    jets_.idxHTagged             .push_back(jet.getIndex()) ; 
    jets_.ptHTagged              .push_back(jet.getPt()) ; 
    jets_.etaHTagged             .push_back(jet.getEta()) ; 
    jets_.phiHTagged             .push_back(jet.getPhi()) ; 
    jets_.MHTagged               .push_back(jet.getMass()) ; 
    jets_.SoftDropMassHTagged    .push_back(jet.getSoftDropMass()) ;
    jets_.PrunedMassHTagged      .push_back(jet.getPrunedMass()) ;
    jets_.tau1HTagged            .push_back(jet.getTau1()) ;
    jets_.tau2HTagged            .push_back(jet.getTau2()) ;
    jets_.tau3HTagged            .push_back(jet.getTau3()) ;
    jets_.csvHTagged             .push_back(jet.getCSV());
    jets_.partonFlavourHTagged   .push_back(jet.getPartonFlavour());
    jets_.hadronFlavourHTagged   .push_back(jet.getHadronFlavour());
    jets_.doubleBHTagged         .push_back(jet.getDoubleBAK8()) ;
    jets_.sj0CSVHTagged          .push_back(jet.getCSVSubjet0()) ;
    jets_.sj1CSVHTagged          .push_back(jet.getCSVSubjet1()) ;
    jets_.hadronFlavourSJ0HTagged.push_back(jet.getHadronFlavourSubjet0()) ;
    jets_.hadronFlavourSJ1HTagged.push_back(jet.getHadronFlavourSubjet1()) ;
  }

  jets_.idxAntiHTagged             .clear() ; jets_.idxAntiHTagged             .reserve(antiHTaggedJets.size()) ;   
  jets_.ptAntiHTagged              .clear() ; jets_.ptAntiHTagged              .reserve(antiHTaggedJets.size()) ;   
  jets_.etaAntiHTagged             .clear() ; jets_.etaAntiHTagged             .reserve(antiHTaggedJets.size()) ;   
  jets_.phiAntiHTagged             .clear() ; jets_.phiAntiHTagged             .reserve(antiHTaggedJets.size()) ;   
  jets_.MAntiHTagged               .clear() ; jets_.MAntiHTagged               .reserve(antiHTaggedJets.size()) ;   
  jets_.SoftDropMassAntiHTagged    .clear() ; jets_.SoftDropMassAntiHTagged    .reserve(antiHTaggedJets.size()) ;   
  jets_.PrunedMassAntiHTagged      .clear() ; jets_.PrunedMassAntiHTagged      .reserve(antiHTaggedJets.size()) ;   
  jets_.tau1AntiHTagged            .clear() ; jets_.tau1AntiHTagged            .reserve(antiHTaggedJets.size()) ;   
  jets_.tau2AntiHTagged            .clear() ; jets_.tau2AntiHTagged            .reserve(antiHTaggedJets.size()) ;   
  jets_.tau3AntiHTagged            .clear() ; jets_.tau3AntiHTagged            .reserve(antiHTaggedJets.size()) ;   
  jets_.csvAntiHTagged             .clear() ; jets_.csvAntiHTagged             .reserve(antiHTaggedJets.size()) ;   
  jets_.partonFlavourAntiHTagged   .clear() ; jets_.partonFlavourAntiHTagged   .reserve(antiHTaggedJets.size()) ;   
  jets_.hadronFlavourAntiHTagged   .clear() ; jets_.hadronFlavourAntiHTagged   .reserve(antiHTaggedJets.size()) ;   
  jets_.doubleBAntiHTagged         .clear() ; jets_.doubleBAntiHTagged         .reserve(antiHTaggedJets.size()) ;   
  jets_.sj0CSVAntiHTagged          .clear() ; jets_.sj0CSVAntiHTagged          .reserve(antiHTaggedJets.size()) ;   
  jets_.sj1CSVAntiHTagged          .clear() ; jets_.sj1CSVAntiHTagged          .reserve(antiHTaggedJets.size()) ;   
  jets_.hadronFlavourSJ0AntiHTagged.clear() ; jets_.hadronFlavourSJ0AntiHTagged.reserve(antiHTaggedJets.size()) ;   
  jets_.hadronFlavourSJ1AntiHTagged.clear() ; jets_.hadronFlavourSJ1AntiHTagged.reserve(antiHTaggedJets.size()) ;   

  for (vlq::Jet jet : antiHTaggedJets) {
    jets_.idxAntiHTagged             .push_back(jet.getIndex()) ; 
    jets_.ptAntiHTagged              .push_back(jet.getPt()) ; 
    jets_.etaAntiHTagged             .push_back(jet.getEta()) ; 
    jets_.phiAntiHTagged             .push_back(jet.getPhi()) ; 
    jets_.MAntiHTagged               .push_back(jet.getMass()) ; 
    jets_.SoftDropMassAntiHTagged    .push_back(jet.getSoftDropMass()) ;
    jets_.PrunedMassAntiHTagged      .push_back(jet.getPrunedMass()) ;
    jets_.tau1AntiHTagged            .push_back(jet.getTau1()) ;
    jets_.tau2AntiHTagged            .push_back(jet.getTau2()) ;
    jets_.tau3AntiHTagged            .push_back(jet.getTau3()) ;
    jets_.csvAntiHTagged             .push_back(jet.getCSV());
    jets_.partonFlavourAntiHTagged   .push_back(jet.getPartonFlavour());
    jets_.hadronFlavourAntiHTagged   .push_back(jet.getHadronFlavour());
    jets_.doubleBAntiHTagged         .push_back(jet.getDoubleBAK8()) ;
    jets_.sj0CSVAntiHTagged          .push_back(jet.getCSVSubjet0()) ;
    jets_.sj1CSVAntiHTagged          .push_back(jet.getCSVSubjet1()) ;
    jets_.hadronFlavourSJ0AntiHTagged.push_back(jet.getHadronFlavourSubjet0()) ;
    jets_.hadronFlavourSJ1AntiHTagged.push_back(jet.getHadronFlavourSubjet1()) ;
  }

  jets_.idxTopTagged             .clear() ; jets_.idxTopTagged             .reserve(goodTopTaggedJets.size()) ;   
  jets_.ptTopTagged              .clear() ; jets_.ptTopTagged              .reserve(goodTopTaggedJets.size()) ;   
  jets_.etaTopTagged             .clear() ; jets_.etaTopTagged             .reserve(goodTopTaggedJets.size()) ;   
  jets_.phiTopTagged             .clear() ; jets_.phiTopTagged             .reserve(goodTopTaggedJets.size()) ;   
  jets_.MTopTagged               .clear() ; jets_.MTopTagged               .reserve(goodTopTaggedJets.size()) ;   
  jets_.SoftDropMassTopTagged    .clear() ; jets_.SoftDropMassTopTagged    .reserve(goodTopTaggedJets.size()) ;   
  jets_.PrunedMassTopTagged      .clear() ; jets_.PrunedMassTopTagged      .reserve(goodTopTaggedJets.size()) ;   
  jets_.tau1TopTagged            .clear() ; jets_.tau1TopTagged            .reserve(goodTopTaggedJets.size()) ;   
  jets_.tau2TopTagged            .clear() ; jets_.tau2TopTagged            .reserve(goodTopTaggedJets.size()) ;   
  jets_.tau3TopTagged            .clear() ; jets_.tau3TopTagged            .reserve(goodTopTaggedJets.size()) ;   
  jets_.csvTopTagged             .clear() ; jets_.csvTopTagged             .reserve(goodTopTaggedJets.size()) ;   
  jets_.partonFlavourTopTagged   .clear() ; jets_.partonFlavourTopTagged   .reserve(goodTopTaggedJets.size()) ;   
  jets_.hadronFlavourTopTagged   .clear() ; jets_.hadronFlavourTopTagged   .reserve(goodTopTaggedJets.size()) ;   
  jets_.doubleBTopTagged         .clear() ; jets_.doubleBTopTagged         .reserve(goodTopTaggedJets.size()) ;   
  jets_.sj0CSVTopTagged          .clear() ; jets_.sj0CSVTopTagged          .reserve(goodTopTaggedJets.size()) ;   
  jets_.sj1CSVTopTagged          .clear() ; jets_.sj1CSVTopTagged          .reserve(goodTopTaggedJets.size()) ;   
  jets_.hadronFlavourSJ0TopTagged.clear() ; jets_.hadronFlavourSJ0TopTagged.reserve(goodTopTaggedJets.size()) ;   
  jets_.hadronFlavourSJ1TopTagged.clear() ; jets_.hadronFlavourSJ1TopTagged.reserve(goodTopTaggedJets.size()) ;   

  for (vlq::Jet jet : goodTopTaggedJets) {
    jets_.idxTopTagged             .push_back(jet.getIndex()) ; 
    jets_.ptTopTagged              .push_back(jet.getPt()) ; 
    jets_.etaTopTagged             .push_back(jet.getEta()) ; 
    jets_.phiTopTagged             .push_back(jet.getPhi()) ; 
    jets_.MTopTagged               .push_back(jet.getMass()) ; 
    jets_.SoftDropMassTopTagged    .push_back(jet.getSoftDropMass()) ;
    jets_.PrunedMassTopTagged      .push_back(jet.getPrunedMass()) ;
    jets_.tau1TopTagged            .push_back(jet.getTau1()) ;
    jets_.tau2TopTagged            .push_back(jet.getTau2()) ;
    jets_.tau3TopTagged            .push_back(jet.getTau3()) ;
    jets_.csvTopTagged             .push_back(jet.getCSV());
    jets_.partonFlavourTopTagged   .push_back(jet.getPartonFlavour());
    jets_.hadronFlavourTopTagged   .push_back(jet.getHadronFlavour());
    jets_.doubleBTopTagged         .push_back(jet.getDoubleBAK8()) ;
    jets_.sj0CSVTopTagged          .push_back(jet.getCSVSubjet0()) ;
    jets_.sj1CSVTopTagged          .push_back(jet.getCSVSubjet1()) ;
    jets_.hadronFlavourSJ0TopTagged.push_back(jet.getHadronFlavourSubjet0()) ;
    jets_.hadronFlavourSJ1TopTagged.push_back(jet.getHadronFlavourSubjet1()) ;
  }

  jets_.idxAntiTopTagged             .clear() ; jets_.idxAntiTopTagged             .reserve(antiTopTaggedJets.size()) ;   
  jets_.ptAntiTopTagged              .clear() ; jets_.ptAntiTopTagged              .reserve(antiTopTaggedJets.size()) ;   
  jets_.etaAntiTopTagged             .clear() ; jets_.etaAntiTopTagged             .reserve(antiTopTaggedJets.size()) ;   
  jets_.phiAntiTopTagged             .clear() ; jets_.phiAntiTopTagged             .reserve(antiTopTaggedJets.size()) ;   
  jets_.MAntiTopTagged               .clear() ; jets_.MAntiTopTagged               .reserve(antiTopTaggedJets.size()) ;   
  jets_.SoftDropMassAntiTopTagged    .clear() ; jets_.SoftDropMassAntiTopTagged    .reserve(antiTopTaggedJets.size()) ;   
  jets_.PrunedMassAntiTopTagged      .clear() ; jets_.PrunedMassAntiTopTagged      .reserve(antiTopTaggedJets.size()) ;   
  jets_.tau1AntiTopTagged            .clear() ; jets_.tau1AntiTopTagged            .reserve(antiTopTaggedJets.size()) ;   
  jets_.tau2AntiTopTagged            .clear() ; jets_.tau2AntiTopTagged            .reserve(antiTopTaggedJets.size()) ;   
  jets_.tau3AntiTopTagged            .clear() ; jets_.tau3AntiTopTagged            .reserve(antiTopTaggedJets.size()) ;   
  jets_.csvAntiTopTagged             .clear() ; jets_.csvAntiTopTagged             .reserve(antiTopTaggedJets.size()) ;   
  jets_.partonFlavourAntiTopTagged   .clear() ; jets_.partonFlavourAntiTopTagged   .reserve(antiTopTaggedJets.size()) ;   
  jets_.hadronFlavourAntiTopTagged   .clear() ; jets_.hadronFlavourAntiTopTagged   .reserve(antiTopTaggedJets.size()) ;   
  jets_.doubleBAntiTopTagged         .clear() ; jets_.doubleBAntiTopTagged         .reserve(antiTopTaggedJets.size()) ;   
  jets_.sj0CSVAntiTopTagged          .clear() ; jets_.sj0CSVAntiTopTagged          .reserve(antiTopTaggedJets.size()) ;   
  jets_.sj1CSVAntiTopTagged          .clear() ; jets_.sj1CSVAntiTopTagged          .reserve(antiTopTaggedJets.size()) ;   
  jets_.hadronFlavourSJ0AntiTopTagged.clear() ; jets_.hadronFlavourSJ0AntiTopTagged.reserve(antiTopTaggedJets.size()) ;   
  jets_.hadronFlavourSJ1AntiTopTagged.clear() ; jets_.hadronFlavourSJ1AntiTopTagged.reserve(antiTopTaggedJets.size()) ;   

  for (vlq::Jet jet : antiTopTaggedJets) {
    jets_.idxAntiTopTagged             .push_back(jet.getIndex()) ; 
    jets_.ptAntiTopTagged              .push_back(jet.getPt()) ; 
    jets_.etaAntiTopTagged             .push_back(jet.getEta()) ; 
    jets_.phiAntiTopTagged             .push_back(jet.getPhi()) ; 
    jets_.MAntiTopTagged               .push_back(jet.getMass()) ; 
    jets_.SoftDropMassAntiTopTagged    .push_back(jet.getSoftDropMass()) ;
    jets_.PrunedMassAntiTopTagged      .push_back(jet.getPrunedMass()) ;
    jets_.tau1AntiTopTagged            .push_back(jet.getTau1()) ;
    jets_.tau2AntiTopTagged            .push_back(jet.getTau2()) ;
    jets_.tau3AntiTopTagged            .push_back(jet.getTau3()) ;
    jets_.csvAntiTopTagged             .push_back(jet.getCSV());
    jets_.partonFlavourAntiTopTagged   .push_back(jet.getPartonFlavour());
    jets_.hadronFlavourAntiTopTagged   .push_back(jet.getHadronFlavour());
    jets_.doubleBAntiTopTagged         .push_back(jet.getDoubleBAK8()) ;
    jets_.sj0CSVAntiTopTagged          .push_back(jet.getCSVSubjet0()) ;
    jets_.sj1CSVAntiTopTagged          .push_back(jet.getCSVSubjet1()) ;
    jets_.hadronFlavourSJ0AntiTopTagged.push_back(jet.getHadronFlavourSubjet0()) ;
    jets_.hadronFlavourSJ1AntiTopTagged.push_back(jet.getHadronFlavourSubjet1()) ;
  }

  tree_->Fill();

  //// Lepton veto 
  vlq::ElectronCollection goodElectrons; 
  electronmaker(evt, goodElectrons) ;

  vlq::MuonCollection goodMuons; 
  muonmaker(evt, goodMuons) ; 

  h1_["nel"]->Fill(goodElectrons.size(),evtwt*toptagsf*evtwtHTDown*btagsf) ;
  h1_["nmu"]->Fill(goodMuons.size(),evtwt*toptagsf*evtwtHTDown*btagsf) ;

  int nel_passreliso(0), nmu_passreliso(0);

  for (vlq::Electron el : goodElectrons ) {
    h1_["elIsoDR"] ->Fill(el.getIso03(),evtwt*toptagsf*evtwtHTDown*btagsf) ; 
    //// Get 2D isolation
    std::pair<double, double> el_drmin_ptrel(Utilities::drmin_pTrel<vlq::Electron, vlq::Jet>(el, goodAK4Jets)) ;
    h2_["elIso2D"] -> Fill(el_drmin_ptrel.first, el_drmin_ptrel.second,evtwt*toptagsf*evtwtHTDown*btagsf) ;
    //// Get rel. iso
    if (el.getIso03() < 0.0571) ++nel_passreliso; 
  }

  for (vlq::Muon mu : goodMuons ) {
    h1_["muIsoDR"] ->Fill(mu.getIso04(),evtwt*toptagsf*evtwtHTDown*btagsf) ; 
    //// Get 2D isolation
    std::pair<double, double> mu_drmin_ptrel(Utilities::drmin_pTrel<vlq::Muon, vlq::Jet>(mu, goodAK4Jets)) ;
    h2_["muIso2D"] -> Fill(mu_drmin_ptrel.first, mu_drmin_ptrel.second,evtwt*toptagsf*evtwtHTDown*btagsf) ;
    //// Get rel. iso
    if ( mu.getIso04() < 0.3) ++nmu_passreliso ;
  }

  h1_["nelAfterRelIso"] -> Fill(nel_passreliso,evtwt*toptagsf*evtwtHTDown*btagsf) ;
  h1_["nmuAfterRelIso"] -> Fill(nmu_passreliso,evtwt*toptagsf*evtwtHTDown*btagsf) ;

  //// Apply 2D isolation
  CandidateCleaner cleanleptons(0.4,40); //// The second argument is for lepton 2D iso, setting to -1 disables it

  cleanleptons(goodElectrons,goodAK4Jets) ;
  cleanleptons(goodMuons,goodAK4Jets) ;

  cleanleptons(goodElectrons,goodAK8Jets) ;
  cleanleptons(goodMuons,goodAK8Jets) ;

  h1_["nelAfter2DIso"]->Fill(goodElectrons.size(),evtwt*toptagsf*evtwtHTDown*btagsf) ;
  h1_["nmuAfter2DIso"]->Fill(goodMuons.size(),evtwt*toptagsf*evtwtHTDown*btagsf) ;

  for (vlq::Electron el : goodElectrons ) {
    h1_["elIsoDRAfter2DIso"] ->Fill(el.getIso03(),evtwt*toptagsf*evtwtHTDown*btagsf) ; 
  }

  for (vlq::Muon mu : goodMuons ) {
    h1_["muIsoDRAfter2DIso"] ->Fill(mu.getIso04(),evtwt*toptagsf*evtwtHTDown*btagsf) ; 
  }

  return true;

}

// ------------ method called once each job just before starting event loop  ------------
void VLQAna::beginJob() {

  h1_["npv_noreweight"] = fs->make<TH1D>("npv_noreweight", ";N(PV);;", 51, -0.5, 50.5) ; 
  h1_["npv"] = fs->make<TH1D>("npv", ";N(PV);;", 51, -0.5, 50.5) ; 

  if ( !doPreselOnly_ ) {
    tree_ = fs->make<TTree>("tree", "TtHT") ; 
    selectedevt_.RegisterTree(tree_,"SelectedEvent") ; 
    jets_.RegisterTree(tree_,"JetInfo") ; 
  }

  if (storePreselEvts_ || doPreselOnly_) {

    h1_["Presel_HT"] = fs->make<TH1D>("Presel_HT", "H_{T};H_{T} [GeV];;",50,500,3000) ; 
    h1_["Presel_Mjj"] = fs->make<TH1D>("Presel_Mjj", "M_{j1,j2} [GeV];;",50,500,3000) ; 

    h1_["Presel_ptAK4_0"] = fs->make<TH1D>("Presel_ptAK4_0", "p_{T} AK4;p_{T} (1st AK4 jet) [GeV];;",50,0,2000) ; 
    h1_["Presel_ptAK4_1"] = fs->make<TH1D>("Presel_ptAK4_1", "p_{T} AK4;p_{T} (2nd AK4 jet) [GeV];;",50,0,2000) ; 
    h1_["Presel_ptAK4_2"] = fs->make<TH1D>("Presel_ptAK4_2", "p_{T} AK4;p_{T} (3rd AK4 jet) [GeV];;",50,0,2000) ; 
    h1_["Presel_ptAK4_3"] = fs->make<TH1D>("Presel_ptAK4_3", "p_{T} AK4;p_{T} (4th AK4 jet) [GeV];;",50,0,2000) ; 
    h1_["Presel_ptAK8_0"] = fs->make<TH1D>("Presel_ptAK8_0", "p_{T} AK8;p_{T} (1st AK8 jet) [GeV];;",50,0,2000) ; 

    h1_["Presel_etaAK4_0"] = fs->make<TH1D>("Presel_etaAK4_0", "#eta AK4;#eta (1st AK4 jet) [GeV];;",60,-6.,6.) ; 
    h1_["Presel_etaAK4_1"] = fs->make<TH1D>("Presel_etaAK4_1", "#eta AK4;#eta (2nd AK4 jet) [GeV];;",60,-6.,6.) ; 
    h1_["Presel_etaAK4_2"] = fs->make<TH1D>("Presel_etaAK4_2", "#eta AK4;#eta (3rd AK4 jet) [GeV];;",60,-6.,6.) ; 
    h1_["Presel_etaAK4_3"] = fs->make<TH1D>("Presel_etaAK4_3", "#eta AK4;#eta (4th AK4 jet) [GeV];;",60,-6.,6.) ; 
    h1_["Presel_etaAK8_0"] = fs->make<TH1D>("Presel_etaAK8_0", "#eta AK8;#eta (1st AK8 jet) [GeV];;",60,-6.,6.) ; 

    h1_["Presel_NAK4"] = fs->make<TH1D>("Presel_NAK4", ";N(AK4 jets);;", 21, -0.5, 20.5) ; 
    h1_["Presel_NAK8"] = fs->make<TH1D>("Presel_NAK8", ";N(AK4 jets);;", 11, -0.5, 10.5) ; 

    if ( doPreselOnly_ ) {
      h1_["cutflow"] = fs->make<TH1D>("cutflow", "cut flow", 5, 0.5, 5.5) ;  
      h1_["cutflow"] -> GetXaxis() -> SetBinLabel(1,  "All") ; 
      h1_["cutflow"] -> GetXaxis() -> SetBinLabel(2,  "Trig.+PV") ; 
      h1_["cutflow"] -> GetXaxis() -> SetBinLabel(3,  "N(AK4)>=4") ; 
      h1_["cutflow"] -> GetXaxis() -> SetBinLabel(4,  "H_{T}>1000GeV") ; 
      h1_["cutflow"] -> GetXaxis() -> SetBinLabel(5,  "N(AK8)>=1") ; 
      return ; 
    }

  }

  h1_["cutflow"] = fs->make<TH1D>("cutflow", "cut flow", 12, 0.5, 12.5) ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(1,  "All") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(2,  "Trig.+PV") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(3,  "N(AK4)>=4") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(4,  "H_{T}>H_{T}^{min}") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(5,  "N(AK8)>=1") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(6,  "p_{T}(AK8_{0})>p_{T}^{min}(AK8_{0})") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(7,  "N(H)>0") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(8,  "N(top)>0") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(9,  "RegionA") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(10, "RegionB") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(11, "RegionC") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(12, "RegionD") ; 

  h1_["RegD_HT"] = fs->make<TH1D>("RegD_HT", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegD_HT_wts"] = fs->make<TH1D>("RegD_HT_wts", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegD_HT_btagsf_bcUp"] = fs->make<TH1D>("RegD_HT_btagsf_bcUp", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegD_HT_btagsf_bcDown"] = fs->make<TH1D>("RegD_HT_btagsf_bcDown", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegD_HT_btagsf_lUp"] = fs->make<TH1D>("RegD_HT_btagsf_lUp", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegD_HT_btagsf_lDown"] = fs->make<TH1D>("RegD_HT_btagsf_lDown", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegD_HT_toptagsfUp"] = fs->make<TH1D>("RegD_HT_toptagsfUp", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegD_HT_toptagsfDown"] = fs->make<TH1D>("RegD_HT_toptagsfDown", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegD_HT_htwtUp"] = fs->make<TH1D>("RegD_HT_htwtUp", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegD_HT_htwtDown"] = fs->make<TH1D>("RegD_HT_htwtDown", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 

  h1_["RegB_HT"] = fs->make<TH1D>("RegB_HT", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegB_HT_wts"] = fs->make<TH1D>("RegB_HT_wts", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegB_HT_btagsf_bcUp"] = fs->make<TH1D>("RegB_HT_btagsf_bcUp", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegB_HT_btagsf_bcDown"] = fs->make<TH1D>("RegB_HT_btagsf_bcDown", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegB_HT_btagsf_lUp"] = fs->make<TH1D>("RegB_HT_btagsf_lUp", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegB_HT_btagsf_lDown"] = fs->make<TH1D>("RegB_HT_btagsf_lDown", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegB_HT_toptagsfUp"] = fs->make<TH1D>("RegB_HT_toptagsfUp", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegB_HT_toptagsfDown"] = fs->make<TH1D>("RegB_HT_toptagsfDown", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegB_HT_htwtUp"] = fs->make<TH1D>("RegB_HT_htwtUp", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 
  h1_["RegB_HT_htwtDown"] = fs->make<TH1D>("RegB_HT_htwtDown", "H_{T};H_{T} [GeV];;",50,1000,3000) ; 

  //////////

  h1_["RegD_mtprime"] = fs->make<TH1D>("RegD_mtprime", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_wts"] = fs->make<TH1D>("RegD_mtprime_wts", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_btagsf_bcUp"] = fs->make<TH1D>("RegD_mtprime_btagsf_bcUp", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_btagsf_bcDown"] = fs->make<TH1D>("RegD_mtprime_btagsf_bcDown", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_btagsf_lUp"] = fs->make<TH1D>("RegD_mtprime_btagsf_lUp", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_btagsf_lDown"] = fs->make<TH1D>("RegD_mtprime_btagsf_lDown", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_toptagsfUp"] = fs->make<TH1D>("RegD_mtprime_toptagsfUp", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_toptagsfDown"] = fs->make<TH1D>("RegD_mtprime_toptagsfDown", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_htwtUp"] = fs->make<TH1D>("RegD_mtprime_htwtUp", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_htwtDown"] = fs->make<TH1D>("RegD_mtprime_htwtDown", "M(T);M(T) [GeV];;",80,500,2500) ; 

  h1_["RegD_mtprime_corr"] = fs->make<TH1D>("RegD_mtprime_corr", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_corr_wts"] = fs->make<TH1D>("RegD_mtprime_corr_wts", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_corr_btagsf_bcUp"] = fs->make<TH1D>("RegD_mtprime_corr_btagsf_bcUp", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_corr_btagsf_bcDown"] = fs->make<TH1D>("RegD_mtprime_corr_btagsf_bcDown", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_corr_btagsf_lUp"] = fs->make<TH1D>("RegD_mtprime_corr_btagsf_lUp", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_corr_btagsf_lDown"] = fs->make<TH1D>("RegD_mtprime_corr_btagsf_lDown", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_corr_toptagsfUp"] = fs->make<TH1D>("RegD_mtprime_corr_toptagsfUp", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_corr_toptagsfDown"] = fs->make<TH1D>("RegD_mtprime_corr_toptagsfDown", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_corr_htwtUp"] = fs->make<TH1D>("RegD_mtprime_corr_htwtUp", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegD_mtprime_corr_htwtDown"] = fs->make<TH1D>("RegD_mtprime_corr_htwtDown", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 

  h1_["RegB_mtprime"] = fs->make<TH1D>("RegB_mtprime", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_wts"] = fs->make<TH1D>("RegB_mtprime_wts", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_btagsf_bcUp"] = fs->make<TH1D>("RegB_mtprime_btagsf_bcUp", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_btagsf_bcDown"] = fs->make<TH1D>("RegB_mtprime_btagsf_bcDown", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_btagsf_lUp"] = fs->make<TH1D>("RegB_mtprime_btagsf_lUp", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_btagsf_lDown"] = fs->make<TH1D>("RegB_mtprime_btagsf_lDown", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_toptagsfUp"] = fs->make<TH1D>("RegB_mtprime_toptagsfUp", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_toptagsfDown"] = fs->make<TH1D>("RegB_mtprime_toptagsfDown", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_htwtUp"] = fs->make<TH1D>("RegB_mtprime_htwtUp", "M(T);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_htwtDown"] = fs->make<TH1D>("RegB_mtprime_htwtDown", "M(T);M(T) [GeV];;",80,500,2500) ; 

  h1_["RegB_mtprime_corr"] = fs->make<TH1D>("RegB_mtprime_corr", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_corr_wts"] = fs->make<TH1D>("RegB_mtprime_corr_wts", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_corr_btagsf_bcUp"] = fs->make<TH1D>("RegB_mtprime_corr_btagsf_bcUp", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_corr_btagsf_bcDown"] = fs->make<TH1D>("RegB_mtprime_corr_btagsf_bcDown", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_corr_btagsf_lUp"] = fs->make<TH1D>("RegB_mtprime_corr_btagsf_lUp", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_corr_btagsf_lDown"] = fs->make<TH1D>("RegB_mtprime_corr_btagsf_lDown", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_corr_toptagsfUp"] = fs->make<TH1D>("RegB_mtprime_corr_toptagsfUp", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_corr_toptagsfDown"] = fs->make<TH1D>("RegB_mtprime_corr_toptagsfDown", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_corr_htwtUp"] = fs->make<TH1D>("RegB_mtprime_corr_htwtUp", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 
  h1_["RegB_mtprime_corr_htwtDown"] = fs->make<TH1D>("RegB_mtprime_corr_htwtDown", "M(T) - M(H-jet) - M(top-jet) + M(H) + M(top);M(T) [GeV];;",80,500,2500) ; 

  h1_["nel"] = fs->make<TH1D>("nel", "nel;N(electrons);Events;;",5,-0.5,4.5) ;
  h1_["nmu"] = fs->make<TH1D>("nmu", "nmu;N(muons);Events;;",5,-0.5,4.5) ;

  h1_["elIsoDR"] = fs->make<TH1D>("elIsoDR", "elIsoDR;Electron rel. iso in #DeltaR < 0.4, EA corr.;Electrons;;",200,0,10) ;
  h1_["muIsoDR"] = fs->make<TH1D>("muIsoDR", "muIsoDR;Muon rel. iso. in #DeltaR < 0.3, #delta#beta corr.;Muons;;",200,0.,10) ;

  h2_["elIso2D"] = fs->make<TH2D>("elIso2D", "elIso2D;#DeltaR(e,jet);p_{T}^{rel} [GeV];Electrons;",20,0,4,200,0,100) ;
  h2_["muIso2D"] = fs->make<TH2D>("muIso2D", "muIso2D;#DeltaR(#mu,jet);p_{T}^{rel} [GeV];Muons;",20,0,4,200,0.,100) ;

  h1_["nelAfterRelIso"] = fs->make<TH1D>("nelAfterRelIso", "nelAfterRelIso;N(electrons);Events;;",5,-0.5,4.5) ;
  h1_["nmuAfterRelIso"] = fs->make<TH1D>("nmuAfterRelIso", "nmuAfterRelIso;N(muons);Events;;",5,-0.5,4.5) ;

  h1_["nelAfter2DIso"] = fs->make<TH1D>("nelAfter2DIso", "nelAfter2DIso;N(electrons);Events;;",5,-0.5,4.5) ;
  h1_["nmuAfter2DIso"] = fs->make<TH1D>("nmuAfter2DIso", "nmuAfter2DIso;N(muons);Events;;",5,-0.5,4.5) ;

  h1_["elIsoDRAfter2DIso"] = fs->make<TH1D>("elIsoDRAfter2DIso", "elIsoDRAfter2DIso;Electron rel. iso in #DeltaR < 0.4, EA corr.;Events;;",200,0,10) ;
  h1_["muIsoDRAfter2DIso"] = fs->make<TH1D>("muIsoDRAfter2DIso", "muIsoDRAfter2DIso;Muon rel. iso. in #DeltaR < 0.3, #delta#beta corr.;Events;;",200,0.,10) ;

}

// ------------ method called once each job just after ending the event loop  ------------
void VLQAna::endJob() {
  return ; 
}

double VLQAna::getHTReweightingSF(double ht, double err) {
  double wt(1);
  std::unique_ptr<TF1> htReweightFun     = std::unique_ptr<TF1>( new TF1("htReweightFun", "1.25868 - (0.000210855*x)",1000, 10000) ) ; 
  std::unique_ptr<TF1> htReweightFunErr  = std::unique_ptr<TF1>( new TF1("htReweightFunErr", 
        "TMath::Sqrt( (0.00432715*0.00432715) + ((3.27325*x/1000000)*(3.27325*x/1000000)) - (2*x*0.00432715*3.27325*0.97547/1000000) )", 
        1000, 10000) ) ;
  wt = htReweightFun->Eval(ht) + (err*htReweightFunErr->Eval(ht)) ; 
  return wt;
}

DEFINE_FWK_MODULE(VLQAna);
