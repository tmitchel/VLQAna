#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Analysis/VLQAna/interface/JetMaker.h"
#include "Analysis/VLQAna/interface/HH4bTree.h"
#include "Analysis/VLQAna/interface/Utilities.h"
#include "Analysis/VLQAna/interface/BTagSFUtils.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TF1.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include <algorithm> 

class HH4b : public edm::EDFilter {
  public:
    explicit HH4b(const edm::ParameterSet&) ;
    ~HH4b() ; 
  private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    bool passHiggsTagging(vlq::Jet jet) ; 

    // ----------member data ---------------------------
    edm::EDGetTokenT<int>           t_runno      ;
    edm::EDGetTokenT<int>           t_lumisec    ;
    edm::EDGetTokenT<int>           t_evtno      ;
    edm::EDGetTokenT<bool>          t_isData     ;
    edm::EDGetTokenT<double>        t_evtwtGen   ;
    edm::EDGetTokenT<double>        t_evtwtPV    ;
    edm::EDGetTokenT<double>        t_evtwtPVLow ;
    edm::EDGetTokenT<double>        t_evtwtPVHigh;
    edm::EDGetTokenT<double>        t_htHat      ;
    edm::EDGetTokenT<int>           t_npuTrue    ;
    edm::EDGetTokenT<vector<int>>   t_lhewtids   ;
    edm::EDGetTokenT<vector<double>>t_lhewts     ;
    edm::EDGetTokenT<bool>          t_hltdecision;
    edm::EDGetTokenT<unsigned>      t_npv        ;
    JetMaker     jetAK8maker                     ; 
    JetMaker     jetHTaggedmaker                 ; 
    const std::string  btageffFile_              ;
    edm::Service<TFileService> fs                ; 
    std::map<std::string, TH1D*> h1_             ; 
    std::map<std::string, TH2D*> h2_             ; 
    hh4b::EventInfoBranches selectedevt; 
    hh4b::DiJetInfoBranches leadingdijets; 
    hh4b::JetInfoBranches ak8jets ; 
    hh4b::JetInfoBranches hjets ; 
    TTree* tree_ ; 
    const bool printEvtNo_;

    //std::unique_ptr<BTagSFUtils> btagsfutils_ ; 
    //std::unique_ptr<BTagSFUtils> doublebtagsfutils_ ; 

    std::ofstream outfile0b_ ; 
    std::ofstream outfile4b_ ; 

    static const double CSVv2L ; 
    static const double DoubleBL;
    static const double DoubleBM;
    static const double DoubleBT;
};

const double HH4b::CSVv2L = 0.460;
const double HH4b::DoubleBL = 0.30;
const double HH4b::DoubleBM = 0.60;
const double HH4b::DoubleBT = 0.80;

using namespace std; 

HH4b::HH4b(const edm::ParameterSet& iConfig) :
  t_runno               (consumes<int>           (iConfig.getParameter<edm::InputTag>("runno"      ))),   
  t_lumisec             (consumes<int>           (iConfig.getParameter<edm::InputTag>("lumisec"    ))),   
  t_evtno               (consumes<int>           (iConfig.getParameter<edm::InputTag>("evtno"      ))),   
  t_isData              (consumes<bool>          (iConfig.getParameter<edm::InputTag>("isData"     ))),   
  t_evtwtGen            (consumes<double>        (iConfig.getParameter<edm::InputTag>("evtwtGen"   ))),   
  t_evtwtPV             (consumes<double>        (iConfig.getParameter<edm::InputTag>("evtwtPV"    ))),   
  t_evtwtPVLow          (consumes<double>        (iConfig.getParameter<edm::InputTag>("evtwtPVLow" ))),   
  t_evtwtPVHigh         (consumes<double>        (iConfig.getParameter<edm::InputTag>("evtwtPVHigh"))),   
  t_htHat               (consumes<double>        (iConfig.getParameter<edm::InputTag>("htHat"      ))),   
  t_npuTrue             (consumes<int>           (iConfig.getParameter<edm::InputTag>("npuTrue"    ))),   
  t_lhewtids            (consumes<vector<int>>   (iConfig.getParameter<edm::InputTag>("lhewtids"   ))),   
  t_lhewts              (consumes<vector<double>>(iConfig.getParameter<edm::InputTag>("lhewts"     ))),   
  t_hltdecision         (consumes<bool>          (iConfig.getParameter<edm::InputTag>("hltdecision"))),   
  t_npv                 (consumes<unsigned>      (iConfig.getParameter<edm::InputTag>("npv"        ))),   
  jetAK8maker           (iConfig.getParameter<edm::ParameterSet> ("jetAK8selParams"),consumesCollector()), 
  jetHTaggedmaker       (iConfig.getParameter<edm::ParameterSet> ("jetHTaggedselParams"),consumesCollector()),
  btageffFile_          (iConfig.getParameter<std::string>       ("btageffFile")), 
  printEvtNo_           (iConfig.getParameter<bool>("printEvtNo"))  
  //btagsfutils_          (new BTagSFUtils()),
  //doublebtagsfutils_    (new BTagSFUtils())
{

}

HH4b::~HH4b() { }

bool HH4b::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<int>           h_runno         ; evt.getByToken(t_runno      , h_runno      ) ; 
  Handle<int>           h_lumisec       ; evt.getByToken(t_lumisec    , h_lumisec    ) ; 
  Handle<int>           h_evtno         ; evt.getByToken(t_evtno      , h_evtno      ) ; 
  Handle<bool>          h_isData        ; evt.getByToken(t_isData     , h_isData     ) ; 
  Handle<double>        h_evtwtGen      ; evt.getByToken(t_evtwtGen   , h_evtwtGen   ) ; 
  Handle<double>        h_evtwtPV       ; evt.getByToken(t_evtwtPV    , h_evtwtPV    ) ; 
  Handle<double>        h_evtwtPVLow    ; evt.getByToken(t_evtwtPVLow , h_evtwtPVLow ) ; 
  Handle<double>        h_evtwtPVHigh   ; evt.getByToken(t_evtwtPVHigh, h_evtwtPVHigh) ; 
  Handle<double>        h_htHat         ; evt.getByToken(t_htHat      , h_htHat      ) ; 
  Handle<int>           h_npuTrue       ; evt.getByToken(t_npuTrue    , h_npuTrue    ) ; 
  Handle<vector<int>>   h_lhewtids      ; evt.getByToken(t_lhewtids   , h_lhewtids   ) ; 
  Handle<vector<double>>h_lhewts        ; evt.getByToken(t_lhewts     , h_lhewts     ) ; 
  Handle<bool>          h_hltdecision   ; evt.getByToken(t_hltdecision, h_hltdecision) ; 
  Handle<unsigned>      h_npv           ; evt.getByToken(t_npv        , h_npv        ) ; 

  const int runno(*h_runno.product()) ; 
  const int lumisec(*h_lumisec.product()) ; 
  const int evtno(*h_evtno.product()) ; 
  const bool isData(*h_isData.product()) ; 
  const double evtwt((*h_evtwtGen.product()) * (*h_evtwtPV.product())) ; 
  const double htHat((*h_htHat.product())) ; 
  const bool hltdecision(*h_hltdecision.product()) ; 
  const unsigned npv(*h_npv.product()) ; 

  vlq::JetCollection  goodAK8Jets ;  
  jetAK8maker(evt, goodAK8Jets); 

  double aK8DEta(goodAK8Jets.size()>= 2 ? abs(goodAK8Jets.at(0).getEta() -  goodAK8Jets.at(1).getEta()) : 999999) ; 
  TLorentzVector p4null; p4null.SetPtEtaPhiE(0,0,0,0); 
  TLorentzVector p4_jet0 = goodAK8Jets.size()>=1 ? goodAK8Jets.at(0).getP4() : p4null ;
  TLorentzVector p4_jet1 = goodAK8Jets.size()>= 2 ? goodAK8Jets.at(1).getP4() : p4null ;
  double mjj(goodAK8Jets.size()>= 2 ? (p4_jet0 + p4_jet1).Mag() : -999999) ;
  double mj_0(goodAK8Jets.size()>= 2 ? goodAK8Jets.at(0).getSoftDropMass() : -999999) ;
  double mj_1(goodAK8Jets.size()>= 2 ? goodAK8Jets.at(1).getSoftDropMass() : -999999) ;
  double mjjred(goodAK8Jets.size()>= 2 ? mjj - mj_0 - mj_1 + 250. : -999999) ;
  double t21_0(goodAK8Jets.size()>= 2 ? 
      ( goodAK8Jets.at(0).getTau1() != 0 ? goodAK8Jets.at(0).getTau2() / goodAK8Jets.at(0).getTau1() : 999999) 
      : 999999 ) ; 
  double t21_1(goodAK8Jets.size()>= 2 ? 
      ( goodAK8Jets.at(1).getTau1() != 0 ? goodAK8Jets.at(1).getTau2() / goodAK8Jets.at(1).getTau1() : 999999) 
      : 999999 ) ; 
  int nsjBTagged(goodAK8Jets.size()>= 2 ? 
      goodAK8Jets.at(0).getNSubjetsBTaggedCSVL() + goodAK8Jets.at(1).getNSubjetsBTaggedCSVL() : -999999) ; 
  double doubleb0(goodAK8Jets.size()>= 1 ?
      goodAK8Jets.at(0).getDoubleBAK8() : -999999) ;
  double doubleb1(goodAK8Jets.size()>= 2 ?
      goodAK8Jets.at(1).getDoubleBAK8() : -999999) ;
  std::bitset<4>doublebcat0(1), doublebcat1(1);
  if (doubleb0 > DoubleBL) doublebcat0 = doublebcat0 << 1;
  if (doubleb0 > DoubleBM) doublebcat0 = doublebcat0 << 1;
  if (doubleb0 > DoubleBT) doublebcat0 = doublebcat0 << 1;
  if (doubleb1 > DoubleBL) doublebcat1 = doublebcat1 << 1;
  if (doubleb1 > DoubleBM) doublebcat1 = doublebcat1 << 1;
  if (doubleb1 > DoubleBT) doublebcat1 = doublebcat1 << 1;
  std::bitset<4>doublebcategory(doublebcat0 | doublebcat1);


  //// Event flags for cut flows
  const bool isHltNPV   (hltdecision && (npv > 0)) ; 
  const bool isAK8Pt    (goodAK8Jets.size()>= 2) ; 
  const bool isAK8Eta   (goodAK8Jets.size()>= 2 && ( abs(goodAK8Jets.at(0).getEta()) <= 2.4 && abs(goodAK8Jets.at(1).getEta()) <= 2.4) ); 
  const bool isAK8DEta  (aK8DEta < 1.3);
  const bool isAK8Mjj   (mjjred > 1100.);
  const bool isAK8Mj    ((mj_0 >= 105 && mj_0 <= 135) && (mj_1 >= 105 && mj_1 <= 135) ) ; 
  const bool isAK8t21   ((t21_0 < 0.55 && t21_1 < 0.55)); 
  const bool isEvt0b    (nsjBTagged == 0) ; 
  const bool isEvt1b    (nsjBTagged == 1) ; 
  const bool isEvt2b    (nsjBTagged == 2) ; 
  const bool isEvt3b    (nsjBTagged == 3) ; 
  const bool isEvt4b    (nsjBTagged == 4) ; 
  const bool isEvt3bHPHP(nsjBTagged == 3 && (t21_0 < 0.55 && t21_1 < 0.55));

  //// Cut flow
  h1_["cutflow"] -> Fill(1) ; 
  if ( isHltNPV  ) h1_["cutflow"] -> Fill(2) ; 
  if ( isHltNPV && isAK8Pt  ) h1_["cutflow"] -> Fill(3) ; 
  if ( isHltNPV && isAK8Pt && isAK8Eta  ) h1_["cutflow"] -> Fill(4) ; 
  if ( isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta ) h1_["cutflow"] -> Fill(5) ; 
  if ( isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mjj  ) h1_["cutflow"] -> Fill(6) ;
  if ( isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mjj && isAK8Mj  ) h1_["cutflow"] -> Fill(7) ;
  if ( isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mjj && isAK8Mj && isAK8t21  ) h1_["cutflow"] -> Fill(8) ;
  if ( isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mjj && isAK8Mj && isAK8t21 && isEvt0b  ) h1_["cutflow"] -> Fill( 9) ;
  if ( isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mjj && isAK8Mj && isAK8t21 && isEvt1b  ) h1_["cutflow"] -> Fill(10) ;
  if ( isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mjj && isAK8Mj && isAK8t21 && isEvt2b  ) h1_["cutflow"] -> Fill(11) ;
  if ( isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mjj && isAK8Mj && isAK8t21 && isEvt3b  ) h1_["cutflow"] -> Fill(12) ;
  if ( isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mjj && isAK8Mj && isAK8t21 && isEvt4b  ) h1_["cutflow"] -> Fill(13) ;
  if ( isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mjj && isAK8Mj && isAK8t21 && isEvt3bHPHP) h1_["cutflow"] -> Fill(14) ;

  if ( printEvtNo_ && isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mjj && isAK8Mj && isAK8t21 && isEvt0b  ) {
    outfile0b_ << runno << " " << lumisec << " " << evtno << " " 
      << goodAK8Jets.at(0).getCSVSubjet0() << " " << goodAK8Jets.at(0).getCSVSubjet1() << " " 
      << goodAK8Jets.at(1).getCSVSubjet0() << " " << goodAK8Jets.at(1).getCSVSubjet1() << " " 
      << std::endl ; 
  }

  if ( printEvtNo_ && isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mjj && isAK8Mj && isAK8t21 && isEvt4b  ) {
    outfile4b_ << runno << " " << lumisec << " " << evtno << " " 
      << goodAK8Jets.at(0).getCSVSubjet0() << " " << goodAK8Jets.at(0).getCSVSubjet1() << " " 
      << goodAK8Jets.at(1).getCSVSubjet0() << " " << goodAK8Jets.at(1).getCSVSubjet1() << " " 
      << std::endl ; 
  }

  std::vector<double>ak8jets_sjcsvs;
  std::vector<double>ak8jets_sjpts;
  std::vector<double>ak8jets_sjetas;
  std::vector<int>   ak8jets_sjflhads;

  std::vector<double>ak8jets_doublebs;
  std::vector<double>ak8jets_pts;
  std::vector<int>   ak8jets_flhads;

  ak8jets.njets = goodAK8Jets.size() ; 
  for (auto ijet : index(goodAK8Jets) ) {
    ak8jets.Index[ijet.first] = (ijet.second).getIndex() ;
    ak8jets.ParentIndex[ijet.first] = -1 ;
    ak8jets.Pt[ijet.first] = (ijet.second).getPt() ;
    ak8jets.Eta[ijet.first] = (ijet.second).getEta() ;
    ak8jets.Phi[ijet.first] = (ijet.second).getPhi() ;
    ak8jets.Energy[ijet.first] = (ijet.second).getEnergy() ;
    ak8jets.Mass[ijet.first] = (ijet.second).getMass() ;
    ak8jets.MassPruned[ijet.first] = (ijet.second).getPrunedMass() ;
    ak8jets.MassSoftDrop[ijet.first] = (ijet.second).getSoftDropMass() ;
    ak8jets.tau1[ijet.first] = (ijet.second).getTau1() ;
    ak8jets.tau2[ijet.first] = (ijet.second).getTau2() ;
    ak8jets.tau3[ijet.first] = (ijet.second).getTau3() ;
    ak8jets.hadFlavour[ijet.first] = (ijet.second).getHadronFlavour() ;
    ak8jets.CSVIVFv2[ijet.first] = (ijet.second).getCSV() ;
    ak8jets.groomedMassCorr[ijet.first] = (ijet.second).getGroomedMassCorr() ;
    ak8jets.nsubjets[ijet.first] = (ijet.second).getNSubjets() ;
    ak8jets.nsubjetsBTaggedCSVL[ijet.first] = (ijet.second).getNSubjetsBTaggedCSVL() ;
    leadingdijets.nsubjetsBTaggedCSVL_  += (ijet.second).getNSubjetsBTaggedCSVL() ;  
    ak8jets.hadFlavourSubjet0[ijet.first] = (ijet.second).getHadronFlavourSubjet0() ;
    ak8jets.hadFlavourSubjet1[ijet.first] = (ijet.second).getHadronFlavourSubjet1() ;
    ak8jets.ptSubjet0[ijet.first] = (ijet.second).getPtSubjet0() ;
    ak8jets.ptSubjet1[ijet.first] = (ijet.second).getPtSubjet1() ;
    ak8jets.etaSubjet0[ijet.first] = (ijet.second).getEtaSubjet0() ;
    ak8jets.etaSubjet1[ijet.first] = (ijet.second).getEtaSubjet1() ;
    ak8jets.phiSubjet0[ijet.first] = (ijet.second).getPhiSubjet0() ;
    ak8jets.phiSubjet1[ijet.first] = (ijet.second).getPhiSubjet1() ;
    ak8jets.energySubjet0[ijet.first] = (ijet.second).getEnergySubjet0() ;
    ak8jets.energySubjet1[ijet.first] = (ijet.second).getEnergySubjet1() ;
    ak8jets.doublesv[ijet.first] = (ijet.second).getDoubleBAK8();
    ak8jets.csvSubjet0[ijet.first] = (ijet.second).getCSVSubjet0() ;
    ak8jets.csvSubjet1[ijet.first] = (ijet.second).getCSVSubjet1() ;

    if ( !isData ) { 

      ak8jets_sjcsvs.push_back((ijet.second).getCSVSubjet0()) ; 
      ak8jets_sjcsvs.push_back((ijet.second).getCSVSubjet1()) ; 

      ak8jets_sjpts.push_back((ijet.second).getPtSubjet0()) ;  
      ak8jets_sjpts.push_back((ijet.second).getPtSubjet1()) ;  

      ak8jets_sjetas.push_back((ijet.second).getEtaSubjet0()) ;  
      ak8jets_sjetas.push_back((ijet.second).getEtaSubjet1()) ;  

      ak8jets_sjflhads.push_back(abs((ijet.second).getHadronFlavourSubjet0())) ; 
      ak8jets_sjflhads.push_back(abs((ijet.second).getHadronFlavourSubjet1())) ;

      ak8jets_doublebs.push_back((ijet.second).getDoubleBAK8());
      ak8jets_pts.push_back((ijet.second).getPt());
      ak8jets_flhads.push_back((ijet.second).getHadronFlavour());

    } //// if !isData do  b tag SFs

  } //// Loop over all AK8 jets

  //// Event selection
  if ( !isHltNPV || !isAK8Pt ) return false; 

  leadingdijets.doublebcategory_ = int(doublebcategory.to_ulong()) ; 

  //// For trigger turn on
  if ( isAK8Pt && isAK8Eta && isAK8DEta && (goodAK8Jets.at(0).getMass() > 50. || goodAK8Jets.at(1).getMass() > 50.)) {
    h1_["npv"] -> Fill(*h_npv.product()); 
    h1_["mjj"] -> Fill(mjj) ;
    h1_["mjjred"] -> Fill(mjjred) ;
    h2_["mjjred_deta"] -> Fill(mjjred, aK8DEta) ;
  }
  if ( isHltNPV && isAK8Pt && isAK8Eta && isAK8DEta && isAK8Mj  ) {
    h1_["npv_MJSel"] -> Fill(*h_npv.product()); 
    h1_["mjj_MJSel"] -> Fill(mjj) ;
    h1_["mjjred_MJSel"] -> Fill(mjjred) ;
    h2_["mjjred_deta_MJSel"] -> Fill(mjjred, aK8DEta) ;
  }

  //// Leading 2 jets pass Higgs tagging
  if ( !passHiggsTagging(goodAK8Jets.at(0)) || !passHiggsTagging(goodAK8Jets.at(1)) ) return false ;

  vlq::JetCollection goodHTaggedJets ;
  goodHTaggedJets.push_back(goodAK8Jets.at(0)); 
  goodHTaggedJets.push_back(goodAK8Jets.at(1)); 

  double btagsf(1) ;
  double btagsf_bcUp(1) ; 
  double btagsf_bcDown(1) ; 
  double btagsf_lUp(1) ; 
  double btagsf_lDown(1) ; 

  TLorentzVector p4_hjet0 = goodHTaggedJets.at(0).getP4() ; 
  TLorentzVector p4_hjet1 = goodHTaggedJets.at(1).getP4() ; 
  TLorentzVector p4_leading2hjets = p4_hjet0+p4_hjet1 ;

  leadingdijets.deta_leading2hjets_ = aK8DEta ; 
  leadingdijets.minv_leading2hjets_ = p4_leading2hjets.Mag();
  leadingdijets.minv_leading2hjets_subtr_ = p4_leading2hjets.Mag() 
    - (goodHTaggedJets.at(0).getMass() - 125.) 
    - (goodHTaggedJets.at(1).getMass() - 125.);
  leadingdijets.pt_leading2hjets_ = p4_leading2hjets.Pt();
  leadingdijets.eta_leading2hjets_ = p4_leading2hjets.Eta();
  leadingdijets.y_leading2hjets_ = p4_leading2hjets.Rapidity();
  leadingdijets.phi_leading2hjets_ = p4_leading2hjets.Phi();
  leadingdijets.energy_leading2hjets_ = p4_leading2hjets.E();

  std::vector<std::pair<int, double>> lhe_ids_wts;
  for (auto idx : index(*h_lhewtids.product()) ) {
    int id = (*h_lhewtids.product()).at(idx.first) ; 
    double wt = (*h_lhewts.product()).at(idx.first) ; 
    lhe_ids_wts.push_back(std::make_pair(id, wt)) ; 
  }

  selectedevt.runno_ = int(runno);
  selectedevt.lumisec_ = int(lumisec);
  selectedevt.evtno_ = int(evtno);
  selectedevt.evtwt_ = evtwt ; 
  selectedevt.evtwtPV_ = double(*h_evtwtPV.product()) ; 
  selectedevt.evtwtPVLow_ = double(*h_evtwtPVLow.product()) ; 
  selectedevt.evtwtPVHigh_ = double(*h_evtwtPVHigh.product()) ; 
  selectedevt.npv_ = npv ; 
  selectedevt.npuTrue_ = int(*h_npuTrue.product()) ; 
  selectedevt.lhewts_ = lhe_ids_wts;
  selectedevt.htHat_ = htHat ; 

  leadingdijets.nsubjetsBTaggedCSVL_ = 0 ;

  std::vector<double>hjets_sjcsvs;
  std::vector<double>hjets_sjpts;
  std::vector<double>hjets_sjetas;
  std::vector<int>hjets_sjflhads;

  std::vector<double>hjets_doublebs;
  std::vector<double>hjets_pts;
  std::vector<int>   hjets_flhads;

  hjets.njets = goodHTaggedJets.size() ; 
  for (auto ijet : index(goodHTaggedJets) ) {
    hjets.Index[ijet.first] = (ijet.second).getIndex() ;
    hjets.ParentIndex[ijet.first] = -1 ;
    hjets.Pt[ijet.first] = (ijet.second).getPt() ;
    hjets.Eta[ijet.first] = (ijet.second).getEta() ;
    hjets.Phi[ijet.first] = (ijet.second).getPhi() ;
    hjets.Energy[ijet.first] = (ijet.second).getEnergy() ;
    hjets.Mass[ijet.first] = (ijet.second).getMass() ;
    hjets.MassPruned[ijet.first] = (ijet.second).getPrunedMass() ;
    hjets.MassSoftDrop[ijet.first] = (ijet.second).getSoftDropMass() ;
    hjets.tau1[ijet.first] = (ijet.second).getTau1() ;
    hjets.tau2[ijet.first] = (ijet.second).getTau2() ;
    hjets.tau3[ijet.first] = (ijet.second).getTau3() ;
    hjets.hadFlavour[ijet.first] = (ijet.second).getHadronFlavour() ;
    hjets.CSVIVFv2[ijet.first] = (ijet.second).getCSV() ;
    hjets.groomedMassCorr[ijet.first] = (ijet.second).getGroomedMassCorr() ;
    hjets.nsubjets[ijet.first] = (ijet.second).getNSubjets() ;
    hjets.nsubjetsBTaggedCSVL[ijet.first] = (ijet.second).getNSubjetsBTaggedCSVL() ;
    leadingdijets.nsubjetsBTaggedCSVL_  += (ijet.second).getNSubjetsBTaggedCSVL() ;  
    hjets.hadFlavourSubjet0[ijet.first] = (ijet.second).getHadronFlavourSubjet0() ;
    hjets.hadFlavourSubjet1[ijet.first] = (ijet.second).getHadronFlavourSubjet1() ;
    hjets.ptSubjet0[ijet.first] = (ijet.second).getPtSubjet0() ;
    hjets.ptSubjet1[ijet.first] = (ijet.second).getPtSubjet1() ;
    hjets.etaSubjet0[ijet.first] = (ijet.second).getEtaSubjet0() ;
    hjets.etaSubjet1[ijet.first] = (ijet.second).getEtaSubjet1() ;
    hjets.phiSubjet0[ijet.first] = (ijet.second).getPhiSubjet0() ;
    hjets.phiSubjet1[ijet.first] = (ijet.second).getPhiSubjet1() ;
    hjets.energySubjet0[ijet.first] = (ijet.second).getEnergySubjet0() ;
    hjets.energySubjet1[ijet.first] = (ijet.second).getEnergySubjet1() ;
    hjets.doublesv[ijet.first] = (ijet.second).getDoubleBAK8();
    hjets.csvSubjet0[ijet.first] = (ijet.second).getCSVSubjet0() ;
    hjets.csvSubjet1[ijet.first] = (ijet.second).getCSVSubjet1() ;

    if ( !isData ) { 

      hjets_sjcsvs.push_back((ijet.second).getCSVSubjet0()) ; 
      hjets_sjcsvs.push_back((ijet.second).getCSVSubjet1()) ; 

      hjets_sjpts.push_back((ijet.second).getPtSubjet0()) ;  
      hjets_sjpts.push_back((ijet.second).getPtSubjet1()) ;  

      hjets_sjetas.push_back((ijet.second).getEtaSubjet0()) ;  
      hjets_sjetas.push_back((ijet.second).getEtaSubjet1()) ;  

      hjets_sjflhads.push_back(abs((ijet.second).getHadronFlavourSubjet0())) ; 
      hjets_sjflhads.push_back(abs((ijet.second).getHadronFlavourSubjet1())) ;

      hjets_doublebs.push_back((ijet.second).getDoubleBAK8());
      hjets_pts.push_back((ijet.second).getPt());
      hjets_flhads.push_back((ijet.second).getHadronFlavour());

    } //// if !isData do  b tag SFs

  } //// Loop over all Higgs jets

  //if ( !isData ) {
  //  btagsfutils_->getBTagSFs (hjets_sjcsvs, hjets_sjpts, hjets_sjetas, hjets_sjflhads, CSVv2L, 
  //      btagsf, btagsf_bcUp, btagsf_bcDown, btagsf_lUp, btagsf_lDown) ; 
  //  doublebtagsfutils_->getBTagSFs (hjets_sjcsvs, hjets_sjpts, hjets_sjetas, hjets_sjflhads, CSVv2L, 
  //      btagsf, btagsf_bcUp, btagsf_bcDown, btagsf_lUp, btagsf_lDown) ; 
  //}

  selectedevt.btagsf_ = btagsf;
  selectedevt.btagsf_bcUp_ = btagsf_bcUp ; 
  selectedevt.btagsf_bcDown_ = btagsf_bcDown ; 
  selectedevt.btagsf_lUp_ = btagsf_lUp ; 
  selectedevt.btagsf_lDown_ = btagsf_lDown ; 

  tree_->Fill();

  return true ; 
}

bool HH4b::passHiggsTagging(vlq::Jet jet) {
  bool passHiggsTagging(0);
  if (jet.getPt() > 300. 
      && abs(jet.getEta()) <= 2.4
      && abs(jet.getPrunedMass()) >= 90 
      //&& (jet.getTau1() == 0. || jet.getTau2()/jet.getTau1() < 0.75)
      && (jet.getTau1() == 0. || jet.getTau2()/jet.getTau1() < 1.00)
     ) passHiggsTagging = true ; 
  return passHiggsTagging;
}

// ------------ method called once each job just before starting event loop  ------------
void HH4b::beginJob() {

  h1_["cutflow"] = fs->make<TH1D>("cutflow", "cut flow", 14, 0.5, 14.5) ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(1,"All") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(2,"Trigger") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(3 ,"p_{T}") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(4,"|#eta|") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(5,"|#Delta#eta|") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(6,"M(jet_{0},jet(1))") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(7,"M(jets)") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(8,"#tau_{21}") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(9,"0b") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(10,"1b") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(11,"2b") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(12,"3b") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(13,"4b") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(14,"3b+HPHP") ; 

  h1_["npv"] = fs->make<TH1D>("npv", ";N(PV);;", 51, -0.5, 50.5) ; 
  h1_["mjj"] = fs->make<TH1D>("mjj", ";M(jj) [GeV];;", 300, 0., 3000.) ; 
  h1_["mjjred"] = fs->make<TH1D>("mjjred", ";M_{red}(jj) [GeV];;", 300, 0., 3000.) ; 
  h2_["mjjred_deta"] = fs->make<TH2D>("mjjred_deta", ";M_{red}(jj) [GeV]; #Delta#eta(jj);", 300, 0., 3000., 3, 0., 1.302) ; 

  h1_["npv_MJSel"] = fs->make<TH1D>("npv_MJSel", ";N(PV);;", 51, -0.5, 50.5) ; 
  h1_["mjj_MJSel"] = fs->make<TH1D>("mjj_MJSel", ";M(jj) [GeV];;", 300, 0., 3000.) ; 
  h1_["mjjred_MJSel"] = fs->make<TH1D>("mjjred_MJSel", ";M_{red}(jj) [GeV];;", 300, 0., 3000.) ; 
  h2_["mjjred_deta_MJSel"] = fs->make<TH2D>("mjjred_deta_MJSel", ";M_{red}(jj) [GeV]; #Delta#eta(jj);", 300, 0., 3000., 3, 0., 1.302) ; 

  tree_ = fs->make<TTree>("tree", "HH4b") ; 
  selectedevt.RegisterTree(tree_,"SelectedEvent") ; 
  leadingdijets.RegisterTree(tree_,"LeadingDiJets") ;
  ak8jets.RegisterTree(tree_,"AK8Jets") ; 
  hjets.RegisterTree(tree_,"HJets") ; 

  if ( printEvtNo_ ) {
    outfile0b_.open("EvtNo_0b.txt",ios::out) ; 
    outfile0b_ << " Cat0b: Runno " << " Lumisec " << " Evtno " 
      << " j0-sj0 CSVv2 " << " j0-sj1 CSVv2 " 
      << " j1-sj0 CSVv2 " << " j1-sj1 CSVv2 " 
      << std::endl ; 

    outfile4b_.open("EvtNo_4b.txt",ios::out) ; 
    outfile4b_ << " Cat4b: Runno " << " Lumisec " << " Evtno " 
      << " j0-sj0 CSVv2 " << " j0-sj1 CSVv2 " 
      << " j1-sj0 CSVv2 " << " j1-sj1 CSVv2 " 
      << std::endl ; 
  }

}

void HH4b::endJob() {

  outfile0b_.close() ;
  outfile4b_.close() ;
  return ; 
}

DEFINE_FWK_MODULE(HH4b);
