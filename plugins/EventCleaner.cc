#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h" 
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include "Analysis/VLQAna/interface/PickGenPart.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include <TFile.h>
#include <TH1D.h>

class EventCleaner : public edm::EDFilter {
  public:
    explicit EventCleaner(const edm::ParameterSet&) ;
    ~EventCleaner() ; 
  private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    virtual void beginRun(edm::Run const& run, edm::EventSetup const& es);
    virtual void endRun(edm::Run const& run, edm::EventSetup const& es);
    double GetLumiWeightsPVBased (const std::string file, const std::string hist, const unsigned npv) ; 

    typedef edm::Handle <std::vector<std::string>> hstring ; 
    typedef edm::Handle <std::vector<float>> hfloat ; 
    typedef edm::Handle <std::vector<int>> hint ; 

    edm::LumiReWeighting LumiWeights_;
    edm::LumiReWeighting LumiWeightsLow_;
    edm::LumiReWeighting LumiWeightsHigh_;

    edm::InputTag l_runno                                  ; 
    edm::InputTag l_lumisec                                ; 
    edm::InputTag l_evtno                                  ; 
    edm::InputTag l_trigName                               ; 
    edm::InputTag l_trigBit                                ; 
    edm::InputTag l_BadChargedCandidateFilter              ;
    edm::InputTag l_BadPFMuonFilter                        ;
    edm::InputTag l_metFiltersName                         ; 
    edm::InputTag l_metFiltersBit                          ; 
    edm::InputTag l_lheEvtProd                             ; 
    std::string   l_genEvtInfoProd                         ; 
    edm::InputTag l_vtxRho                                 ; 
    edm::InputTag l_vtxZ                                   ; 
    edm::InputTag l_vtxNdf                                 ; 
    edm::InputTag l_npv                                    ;
    edm::InputTag l_puNtrueInt                             ;
    const std::string hltORAND_                            ;
    std::vector<std::string> hltPaths_                     ; 
    std::vector<std::string> metFilters_                   ; 
    const bool isData_                                     ; 
    const bool doPUReweightingOfficial_                    ;
    const std::string file_PVWt_                           ; 
    const std::string file_PUDistData_                     ;
    const std::string file_PUDistDataLow_                  ;
    const std::string file_PUDistDataHigh_                 ;
    const std::string file_PUDistMC_                       ;
    const std::string hist_PVWt_                           ; 
    const std::string hist_PUDistData_                     ;
    const std::string hist_PUDistMC_                       ;
    const bool cleanEvents_                                ; 
    const bool storeLHEWts_                                ;
    edm::EDGetTokenT<GenEventInfoProduct> t_genEvtInfoProd ; 
    edm::EDGetTokenT<LHEEventProduct> t_lheEvtProd         ; 
    std::vector<std::string> lhe_weight_labels_            ;

    hstring h_trigName                                     ;
    hstring h_metFiltersName                               ;

    edm::ParameterSet TtZParams_                           ;
    edm::ParameterSet TtHParams_                           ;
    edm::ParameterSet TbWParams_                           ;
    edm::ParameterSet BbZParams_                           ;
    edm::ParameterSet BbHParams_                           ;
    edm::ParameterSet BtWParams_                           ;

    PickGenPart pickTtZ                                    ; 
    PickGenPart pickTtH                                    ; 
    PickGenPart pickTbW                                    ; 
    PickGenPart pickBbZ                                    ; 
    PickGenPart pickBbH                                    ; 
    PickGenPart pickBtW                                    ; 

};

EventCleaner::EventCleaner(const edm::ParameterSet& iConfig) :
  l_runno                    (iConfig.getParameter<edm::InputTag>            ("runnoLabel")),
  l_lumisec                  (iConfig.getParameter<edm::InputTag>            ("lumisecLabel")),
  l_evtno                    (iConfig.getParameter<edm::InputTag>            ("evtnoLabel")),
  l_trigName                 (iConfig.getParameter<edm::InputTag>            ("trigNameLabel")),
  l_trigBit                  (iConfig.getParameter<edm::InputTag>            ("trigBitLabel")),
  l_BadChargedCandidateFilter(iConfig.getParameter<edm::InputTag>            ("BadChargedCandidateFilter")),
  l_BadPFMuonFilter          (iConfig.getParameter<edm::InputTag>            ("BadPFMuonFilter")),             
  l_metFiltersName           (iConfig.getParameter<edm::InputTag>            ("metFiltersNameLabel")),
  l_metFiltersBit            (iConfig.getParameter<edm::InputTag>            ("metFiltersBitLabel")),
  l_lheEvtProd               (iConfig.getParameter<std::string>              ("lheProdName")),
  l_genEvtInfoProd           (iConfig.getParameter<std::string>              ("genEvtInfoProdName")),
  l_vtxRho                   (iConfig.getParameter<edm::InputTag>            ("vtxRhoLabel")),  
  l_vtxZ                     (iConfig.getParameter<edm::InputTag>            ("vtxZLabel")),  
  l_vtxNdf                   (iConfig.getParameter<edm::InputTag>            ("vtxNdfLabel")),  
  l_npv                      (iConfig.getParameter<edm::InputTag>            ("npvLabel")),  
  l_puNtrueInt               (iConfig.getParameter<edm::InputTag>            ("puNtrueIntLabel")),  
  hltORAND_                  (iConfig.getParameter<std::string>              ("hltORAND")), 
  hltPaths_                  (iConfig.getParameter<std::vector<std::string>> ("hltPaths")), 
  metFilters_                (iConfig.getParameter<std::vector<std::string>> ("metFilters")), 
  isData_                    (iConfig.getParameter<bool>                     ("isData")),
  doPUReweightingOfficial_   (iConfig.getParameter<bool>                     ("DoPUReweightingOfficial")),
  file_PVWt_                 (iConfig.getParameter<std::string>              ("File_PVWt")),
  file_PUDistData_           (iConfig.getParameter<std::string>              ("File_PUDistData")),
  file_PUDistDataLow_        (iConfig.getParameter<std::string>              ("File_PUDistDataLow")),
  file_PUDistDataHigh_       (iConfig.getParameter<std::string>              ("File_PUDistDataHigh")),
  file_PUDistMC_             (iConfig.getParameter<std::string>              ("File_PUDistMC")),
  hist_PVWt_                 (iConfig.getParameter<std::string>              ("Hist_PVWt")),
  hist_PUDistData_           (iConfig.getParameter<std::string>              ("Hist_PUDistData")),
  hist_PUDistMC_             (iConfig.getParameter<std::string>              ("Hist_PUDistMC")),
  cleanEvents_               (iConfig.getParameter<bool>                     ("cleanEvents")),
  storeLHEWts_               (iConfig.getParameter<bool>                     ("storeLHEWts")),
  TtZParams_                 (iConfig.getParameter<edm::ParameterSet>        ("TtZParams")),  
  TtHParams_                 (iConfig.getParameter<edm::ParameterSet>        ("TtHParams")),  
  TbWParams_                 (iConfig.getParameter<edm::ParameterSet>        ("TbWParams")),  
  BbZParams_                 (iConfig.getParameter<edm::ParameterSet>        ("BbZParams")),  
  BbHParams_                 (iConfig.getParameter<edm::ParameterSet>        ("BbHParams")),  
  BtWParams_                 (iConfig.getParameter<edm::ParameterSet>        ("BtWParams")),  
  pickTtZ                    (TtZParams_,consumesCollector()), 
  pickTtH                    (TtHParams_,consumesCollector()), 
  pickTbW                    (TbWParams_,consumesCollector()), 
  pickBbZ                    (BbZParams_,consumesCollector()), 
  pickBbH                    (BbHParams_,consumesCollector()), 
  pickBtW                    (BtWParams_,consumesCollector())  
{

  consumes<unsigned int>            (l_runno         ) ; 
  consumes<unsigned int>            (l_lumisec       ) ; 
  consumes<ULong64_t>               (l_evtno         ) ; 
  consumes<std::vector<std::string>, edm::InRun>(l_trigName) ; 
  consumes<std::vector<float>>      (l_trigBit       ) ; 
  consumes<bool>                    (l_BadChargedCandidateFilter) ; 
  consumes<bool>                    (l_BadPFMuonFilter) ; 
  consumes<std::vector<std::string>, edm::InRun>(l_metFiltersName) ; 
  consumes<std::vector<float>>      (l_metFiltersBit ) ; 
  consumes<std::vector<float>>      (l_vtxRho        ) ; 
  consumes<std::vector<float>>      (l_vtxZ          ) ; 
  consumes<std::vector<int>>        (l_vtxNdf        ) ; 
  consumes<int>                     (l_npv           ) ; 
  consumes<int>                     (l_puNtrueInt    ) ; 

  if ( !isData_ ) {
  }

  if (doPUReweightingOfficial_) {
    LumiWeights_     = edm::LumiReWeighting(file_PUDistMC_, file_PUDistData_    , hist_PUDistMC_, hist_PUDistData_) ;
    LumiWeightsLow_  = edm::LumiReWeighting(file_PUDistMC_, file_PUDistDataLow_ , hist_PUDistMC_, hist_PUDistData_) ;
    LumiWeightsHigh_ = edm::LumiReWeighting(file_PUDistMC_, file_PUDistDataHigh_, hist_PUDistMC_, hist_PUDistData_) ;
  }
  t_genEvtInfoProd = consumes<GenEventInfoProduct>(l_genEvtInfoProd) ; 
  if ( storeLHEWts_ ) {
    consumes<LHERunInfoProduct>({"externalLHEProducer"});
    t_lheEvtProd = consumes<LHEEventProduct>(l_lheEvtProd) ; 
  }
  produces<int>("evtno");
  produces<int>("lumisec");
  produces<int>("runno");
  produces<bool>("isData");
  produces<bool>("hltdecision");
  produces<std::string>("evttype"); 
  produces<double>("evtwtGen");
  produces<double>("evtwtPV");
  produces<double>("evtwtPVLow");
  produces<double>("evtwtPVHigh");
  produces<unsigned>("npv");
  produces<int>("npuTrue");
  produces<double>("htHat");
  produces<std::vector<int>>("lhewtids") ; 
  produces<std::vector<double>>("lhewts") ; 

}

EventCleaner::~EventCleaner() {}

bool EventCleaner::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;

  Handle<unsigned int> h_runno    ; evt.getByLabel (l_runno                  , h_runno                );
  Handle<unsigned int> h_lumisec  ; evt.getByLabel (l_lumisec                , h_lumisec              );
  Handle<ULong64_t>    h_evtno    ; evt.getByLabel (l_evtno                  , h_evtno                );
  Handle<bool>         h_BadChargedCandidateFilter; evt.getByLabel (l_BadChargedCandidateFilter, h_BadChargedCandidateFilter); 
  Handle<bool>         h_BadPFMuonFilter          ; evt.getByLabel (l_BadPFMuonFilter          , h_BadPFMuonFilter          ); 
  hfloat  h_trigBit               ; evt.getByLabel (l_trigBit                , h_trigBit              ); 
  hfloat  h_metFiltersBit         ; evt.getByLabel (l_metFiltersBit          , h_metFiltersBit        ); 
  hfloat  h_vtxRho                ; evt.getByLabel (l_vtxRho                 , h_vtxRho               );
  hfloat  h_vtxZ                  ; evt.getByLabel (l_vtxZ                   , h_vtxZ                 );
  hint    h_vtxNdf                ; evt.getByLabel (l_vtxNdf                 , h_vtxNdf               );
  Handle<int>   h_npv             ; evt.getByLabel (l_npv                    , h_npv                  );
  Handle<int>   h_puNtrueInt      ; evt.getByLabel (l_puNtrueInt             , h_puNtrueInt           );

  const int runno  ( (2*isData_ - 1) * (*h_runno) )   ; //// If MC, -ve sign for runno  
  const int lumisec( (2*isData_ - 1) * (*h_lumisec) ) ; //// If MC, -ve sign for lumisec
  const int evtno  ( (2*isData_ - 1) * (*h_evtno) )   ; //// If MC, -ve sign for evtno  

  const size_t ntrigs (hltPaths_.size());
  boost::dynamic_bitset<> hltdecisions(std::max(int(ntrigs),1));
  if ( ntrigs > 0 ) { 
    for ( size_t i = 0; i < ntrigs; ++i) {
      const string& myhltpath = hltPaths_.at(i);
      vector<string>::const_iterator it;
      for (it = h_trigName.product()->begin(); it != (h_trigName.product())->end(); ++it ) {
        if ( it->find(myhltpath) < std::string::npos) {
          hltdecisions[i] = (h_trigBit.product())->at( it - (h_trigName.product())->begin() ) ;
        }
      }
    }
  }
  else {
    hltdecisions[0] = 1;
  }
  if ( cleanEvents_ ) {
    if ( hltORAND_ == "OR" ) { 
      if ( !hltdecisions.any() ) {
        return false ; 
      }
    }
    else if ( hltORAND_ == "AND" ) {
      int hltdecision = 1;
      for (size_t i = 0; i < ntrigs; ++i) {
        hltdecision *= hltdecisions[i];
      }
      if (hltdecision == 0) return false;
    }
    else {
      edm::LogError("EventCleaner::hltdecison type set to") << hltORAND_ << ". Not acceptable!!! Set hltORAND to either OR or AND !!!!";
      return false;
    }
  }

  bool metfilterdecision(1) ; 
  if ( *h_BadChargedCandidateFilter.product() == false || *h_BadPFMuonFilter.product() == false) return false;
  if ( isData_ ) {
    for ( const string& metfilter : metFilters_ ) {
      vector<string>::const_iterator it ; 
      for (it = (h_metFiltersName.product())->begin(); it != (h_metFiltersName.product())->end(); ++it) {
        if ( it->find(metfilter) < std::string::npos) {
          metfilterdecision *= (h_metFiltersBit.product())->at( it - (h_metFiltersName.product())->begin() ) ; 
        }
      }
    }
    if ( !metfilterdecision ) return false ; 
  }

  const int npv(*h_npv) ; 

  double evtwtGen(1.0) ; 
  double htHat(0.0);
  std::vector<int> lhewtids;
  std::vector<double> lhewts;
  if ( !isData_ && storeLHEWts_ ) {

    Handle<LHEEventProduct> h_lheEvtProd;
    evt.getByToken(t_lheEvtProd, h_lheEvtProd); 
    double nominal_wt = h_lheEvtProd->weights()[0].wgt ; // h_lheEvtProd->hepeup().XWGTUP ; //h_lheEvtProd->originalXWGTUP();
    lhewtids.reserve(h_lheEvtProd->weights().size()) ;
    lhewts.reserve(h_lheEvtProd->weights().size()) ; 
    for (unsigned iwt = 0; iwt < h_lheEvtProd->weights().size(); ++iwt) {
      lhewtids.push_back(std::stoi(h_lheEvtProd->weights()[iwt].id)) ;
      lhewts.push_back(h_lheEvtProd->weights()[iwt].wgt / nominal_wt) ; 
    }

    std::vector<lhef::HEPEUP::FiveVector> lheParticles = h_lheEvtProd->hepeup().PUP;
    for(size_t idxPart = 0; idxPart < lheParticles.size();++idxPart){
      unsigned absPdgId = TMath::Abs(h_lheEvtProd->hepeup().IDUP[idxPart]);
      unsigned status = h_lheEvtProd->hepeup().ISTUP[idxPart];
      if(status==1 &&((absPdgId >=1 &&absPdgId<=6) || absPdgId == 21)){
        htHat += sqrt(pow(lheParticles[idxPart][0],2.) + pow(lheParticles[idxPart][1],2.));
      }
    }

  }
  else {
    lhewtids.push_back(-999999);
    lhewts.push_back(-999999) ;
  }

  int npuTrue(-1);
  double evtwtPV(1.0) ;
  double evtwtPVLow(1.0) ;
  double evtwtPVHigh(1.0) ;
  if ( !isData_ ) {
    npuTrue = *h_puNtrueInt ; 
    if ( doPUReweightingOfficial_ ) { 
      evtwtPV *= LumiWeights_.weight(npuTrue) ; 
      evtwtPVLow *= LumiWeightsLow_.weight(npuTrue) ; 
      evtwtPVHigh *= LumiWeightsHigh_.weight(npuTrue) ; 
    }
  }

  string evttype(isData_ ? "EvtType_Data" : "EvtType_MC");

  if ( !isData_ ) {
    vlq::GenParticleCollection vlqTtZ = pickTtZ(evt) ;  
    vlq::GenParticleCollection vlqTtH = pickTtH(evt) ;  
    vlq::GenParticleCollection vlqTbW = pickTbW(evt) ;  
    vlq::GenParticleCollection vlqBbZ = pickBbZ(evt) ;  
    vlq::GenParticleCollection vlqBbH = pickBbH(evt) ;  
    vlq::GenParticleCollection vlqBtW = pickBtW(evt) ;  

    if ( vlqTtZ.size() == 2 ) evttype = "EvtType_MC_tZtZ" ; 
    if ( vlqTtH.size() == 2 ) evttype = "EvtType_MC_tHtH" ; 
    if ( vlqTbW.size() == 2 ) evttype = "EvtType_MC_bWbW" ; 

    if ( vlqTtZ.size() == 1 && vlqTtH.size() == 1 ) evttype = "EvtType_MC_tZtH" ; 
    if ( vlqTtZ.size() == 1 && vlqTbW.size() == 1 ) evttype = "EvtType_MC_tZbW" ; 
    if ( vlqTtH.size() == 1 && vlqTbW.size() == 1 ) evttype = "EvtType_MC_tHbW" ; 

    if ( vlqBbZ.size() == 2 ) evttype = "EvtType_MC_bZbZ" ; 
    if ( vlqBbH.size() == 2 ) evttype = "EvtType_MC_bHbH" ; 
    if ( vlqBtW.size() == 2 ) evttype = "EvtType_MC_tWtW" ; 

    if ( vlqBbZ.size() == 1 && vlqBbH.size() == 1 ) evttype = "EvtType_MC_bZbH" ; 
    if ( vlqBbZ.size() == 1 && vlqBtW.size() == 1 ) evttype = "EvtType_MC_bZtW" ; 
    if ( vlqBbH.size() == 1 && vlqBtW.size() == 1 ) evttype = "EvtType_MC_bHtW" ; 
  }

  unique_ptr<int>ptr_evtno(new int(evtno)); 
  unique_ptr<int>ptr_lumisec(new int(lumisec)); 
  unique_ptr<int>ptr_runno(new int(runno)); 
  unique_ptr<bool>ptr_isData(new bool(isData_)); 
  unique_ptr<bool>ptr_hltdecision(new bool(hltdecisions.any())); 
  unique_ptr<string>ptr_evttype(new string(evttype)); 
  unique_ptr<double>ptr_evtwtGen(new double(evtwtGen)); 
  unique_ptr<double>ptr_evtwtPV(new double(evtwtPV)); 
  unique_ptr<double>ptr_evtwtPVLow(new double(evtwtPVLow)); 
  unique_ptr<double>ptr_evtwtPVHigh(new double(evtwtPVHigh)); 
  unique_ptr<unsigned>ptr_npv(new unsigned(npv)); 
  unique_ptr<int>ptr_npuTrue(new int(npuTrue)); 
  unique_ptr<double>ptr_htHat(new double(htHat)); 
  unique_ptr<std::vector<int>> ptr_lhewtids(new std::vector<int>(lhewtids));
  unique_ptr<std::vector<double>> ptr_lhewts(new std::vector<double>(lhewts));

  evt.put(std::move(ptr_evtno), "evtno");
  evt.put(std::move(ptr_lumisec), "lumisec");
  evt.put(std::move(ptr_runno), "runno");
  evt.put(std::move(ptr_isData), "isData");
  evt.put(std::move(ptr_hltdecision), "hltdecision");
  evt.put(std::move(ptr_evttype), "evttype");
  evt.put(std::move(ptr_evtwtGen), "evtwtGen");
  evt.put(std::move(ptr_evtwtPV), "evtwtPV");
  evt.put(std::move(ptr_evtwtPVLow), "evtwtPVLow");
  evt.put(std::move(ptr_evtwtPVHigh), "evtwtPVHigh");
  evt.put(std::move(ptr_npv), "npv");
  evt.put(std::move(ptr_npuTrue), "npuTrue");
  evt.put(std::move(ptr_htHat), "htHat");
  evt.put(std::move(ptr_lhewtids), "lhewtids") ;
  evt.put(std::move(ptr_lhewts), "lhewts") ;

  return true ; 

}

void EventCleaner::beginJob() {
}

void EventCleaner::endJob() {
}

void EventCleaner::beginRun(edm::Run const& run, edm::EventSetup const& es) {

  run.getByLabel (l_trigName               , h_trigName             );
  run.getByLabel (l_metFiltersName         , h_metFiltersName       );

}

void EventCleaner::endRun(edm::Run const& run, edm::EventSetup const& es) {
}

double EventCleaner::GetLumiWeightsPVBased (const std::string file, const std::string hist, const unsigned npv) { 
  double wtPU(1) ;
  TFile f(file.c_str()) ;
  TH1D* hwt = dynamic_cast<TH1D*>( f.Get( hist.c_str() ) ) ; 
  wtPU = npv > 0 && npv <= 60 ? hwt->GetBinContent(npv) : 1.; 
  delete hwt ; 
  f.Close() ; 
  return wtPU ;
}

DEFINE_FWK_MODULE(EventCleaner);
