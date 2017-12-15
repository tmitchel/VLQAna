#include <iostream>
#include <vector>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "AnalysisDataFormats/BoostedObjects/interface/ResolvedVjj.h"
#include "AnalysisDataFormats/BoostedObjects/interface/Jet.h"
#include "AnalysisDataFormats/BoostedObjects/interface/Candidate.h"

#include "Analysis/VLQAna/interface/JetMaker.h"
#include "Analysis/VLQAna/interface/JetID.h"
#include "Analysis/VLQAna/interface/PickGenPart.h"

#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>

using namespace std;

class MassReco : public edm::EDFilter {
public:
  explicit MassReco(const edm::ParameterSet&);
  ~MassReco();

private:
  virtual void beginJob() override;
  virtual bool filter(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  pair<double, double> vector_eval(vector<pair<double, double> >);
  double resolvedChi2(vector<TLorentzVector>, TLorentzVector, double, double);
  double boostedChi2(vector<TLorentzVector>, TLorentzVector, TLorentzVector, double, double, double);
  pair<double, double> doBoostedReco(vector<vlq::Jet>, TLorentzVector, double, TLorentzVector, double);
  pair<double, double> doResolvedReco(vector<vlq::Jet>, double, TLorentzVector);
  	
	edm::EDGetTokenT<vector<vlq::Jet> >  ak8_t;
  edm::EDGetTokenT<vector<vlq::Jet> >  ak4_t;
  edm::EDGetTokenT<vector<vlq::Jet> >  bjets_t;
  edm::EDGetTokenT<vector<vlq::Jet> >  zjets_t;
  edm::EDGetTokenT<vector<vlq::Jet> >  hjets_t;
  edm::EDGetTokenT<vector<vlq::Candidate> > zllcands_t;
  edm::EDGetTokenT<double>             st_t;
  edm::EDGetTokenT<double>             evtwt_t;
  edm::EDGetTokenT<vector<int>>        t_lhewtids;
  edm::EDGetTokenT<vector<double>>     t_lhewts;  
  edm::EDGetTokenT<double>             bTagwt_t;
  edm::EDGetTokenT<double>             sjbTagwt_t;


  const double ptMin_;
	const double STMaxControl_;
	const double STMin_;
  const string zdecaymode_;
  const string signalType_;
	const bool optimizeReco_;
	const bool controlReco_;
  const bool doSkim_;
  const bool isData_;
  const int pdfID_offset_;
  const int scale_offset_;
  const bool syst_;
  const bool vv_;
	PickGenPart genpart ;

  edm::Service<TFileService> fs;
  map<string, TH1D*> h1_;
  map<string, TH2D*> h2_;

};

MassReco::MassReco(const edm::ParameterSet& iConfig) :
	ak8_t         (consumes<vector<vlq::Jet> > (iConfig.getParameter<edm::InputTag>("ak8jets"))),
  ak4_t         (consumes<vector<vlq::Jet> > (iConfig.getParameter<edm::InputTag>("jets"))),
  bjets_t       (consumes<vector<vlq::Jet> > (iConfig.getParameter<edm::InputTag>("bjets"))),
  zjets_t       (consumes<vector<vlq::Jet> > (iConfig.getParameter<edm::InputTag>("zjets"))),
  hjets_t       (consumes<vector<vlq::Jet> > (iConfig.getParameter<edm::InputTag>("hjets"))),
  zllcands_t    (consumes<vector<vlq::Candidate> > (iConfig.getParameter<edm::InputTag>("zllcands"))),
  st_t          (consumes<double>           (iConfig.getParameter<edm::InputTag>("st"))),
  evtwt_t       (consumes<double>           (iConfig.getParameter<edm::InputTag>("Prewt"))),
  t_lhewtids              (consumes<vector<int> >     (iConfig.getParameter<edm::InputTag>("lhewtids"))),
  t_lhewts (consumes<vector<double> > (iConfig.getParameter<edm::InputTag>("lhewts"))),
  bTagwt_t      (consumes<double>           (iConfig.getParameter<edm::InputTag>("bTagwt"))),
  sjbTagwt_t    (consumes<double>           (iConfig.getParameter<edm::InputTag>("sjbTagwt"))),

 
  ptMin_     (iConfig.getParameter<double> ("ptMin")),
	STMaxControl_ (iConfig.getParameter<double> ("STMaxControl")),
	STMin_     (iConfig.getParameter<double> ("STMin")),
  zdecaymode_       (iConfig.getParameter<string> ("zdecaymode")),
  signalType_ (iConfig.getParameter<string> ("signalType")),
	optimizeReco_ (iConfig.getParameter<bool> ("optimizeReco")),
	controlReco_  (iConfig.getParameter<bool> ("controlReco")),
  doSkim_       (iConfig.getParameter<bool> ("doSkim")),
  isData_       (iConfig.getParameter<bool> ("isData")),
  pdfID_offset_ (iConfig.getParameter<int>  ("pdfID_offset")),
  scale_offset_ (iConfig.getParameter<int>  ("scale_offset")),
  syst_         (iConfig.getParameter<bool> ("syst")),
  vv_           (iConfig.getParameter<bool> ("vv")),

	genpart    (iConfig.getParameter<edm::ParameterSet>("genParams"),consumesCollector())

{}

MassReco::~MassReco() {}

bool MassReco::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
	
  edm::Handle<vector<vlq::Candidate> > zllcands_h; evt.getByToken(zllcands_t, zllcands_h);
	edm::Handle<vector<vlq::Jet> > ak8_h      ; evt.getByToken(ak8_t,  ak8_h)     ;
	edm::Handle<vector<vlq::Jet> > ak4_h      ; evt.getByToken(ak4_t,  ak4_h)     ;
  edm::Handle<vector<vlq::Jet> > bjets_h    ; evt.getByToken(bjets_t,  bjets_h) ;
  edm::Handle<vector<vlq::Jet> > zjets_h    ; evt.getByToken(zjets_t,  zjets_h) ;
  edm::Handle<vector<vlq::Jet> > hjets_h    ; evt.getByToken(hjets_t,  hjets_h) ;
  edm::Handle<double>            st_h       ; evt.getByToken(st_t,     st_h)    ;
  edm::Handle<double>            evtwt_h    ; evt.getByToken(evtwt_t, evtwt_h);
  Handle<vector<int> > h_lhewtids; evt.getByToken(t_lhewtids    ,h_lhewtids   ) ;
  Handle<vector<double> > h_lhewts;evt.getByToken(t_lhewts ,h_lhewts ) ;
  edm::Handle<double>            bTagwt_h   ; evt.getByToken(bTagwt_t, bTagwt_h);
  edm::Handle<double>            sjbTagwt_h ; evt.getByToken(sjbTagwt_t, sjbTagwt_h);

  double evtwt;
  if (evtwt_h.isValid() && bTagwt_h.isValid() && sjbTagwt_h.isValid())
   evtwt = *evtwt_h.product() * *bTagwt_h.product() * *sjbTagwt_h.product() ;
  else
    return false;
 
  vector<pair<int,double>> lhe_id_wts;
  if (!isData_){
    for (unsigned i=0; i<(*h_lhewtids.product()).size(); i++){
      int id = (*h_lhewtids.product()).at(i);
      double wt = (*h_lhewts.product()).at(i);
      lhe_id_wts.push_back(make_pair(id, wt));
    }
  }
  
  TLorentzVector zllcand = (*zllcands_h.product()).at(0).getP4();
  double ST = *st_h.product();

  vector<vlq::Jet> ak4s, bjets, ak8s, zjets, hjets;
  if (ak4_h.isValid()){
    ak4s = *ak4_h.product();
    bjets = *bjets_h.product();
  }
  else
    return false;
  if (ak8_h.isValid()) {
  	ak8s = *ak8_h.product();
    zjets = *zjets_h.product();
    hjets = *hjets_h.product();
  }

  h1_["pt_z"] -> Fill(zllcand.Pt(), evtwt);
  h1_["pt_ak41"] -> Fill(ak4s.at(0).getPt(), evtwt);
  if (hjets.size() > 0) {
    h1_["pt_hjet"] -> Fill(hjets.at(0).getPt(), evtwt);
    h1_["nhjet"] -> Fill(hjets.size(), evtwt);
  }
  if (zjets.size() > 0){
    h1_["pt_zjet"] -> Fill(zjets.at(0).getPt(), evtwt);
    h1_["nzjet"] -> Fill(zjets.size(), evtwt);
  }
  h1_["st"] -> Fill(ST, evtwt);
  h1_["nbjet"] -> Fill(bjets.size(), evtwt);
  h1_["nak4"] -> Fill(ak4s.size(), evtwt);


	double chiCut_ = 100000;

 	////////////////////////
	//Begin SR Mass Reco. //	///////////////////////////////////////////////////
	////////////////////////

	if (ak4s.at(0).getPt() > 100 && ak4s.at(1).getPt() > 50 && bjets.size() == 1 && ST > STMin_){

    h1_["recoCutflow"] -> Fill(1, evtwt);

		pair<double, double> resReco_bZ, boostReco_bZ;
  	pair<double, double> resReco_bH, boostReco_bH;
 
 		resReco_bZ.first = 9999;
  	resReco_bZ.second = -1;
  	resReco_bH.first = 9999;
  	resReco_bH.second = -1;
  	boostReco_bZ.first = 9999;
  	boostReco_bZ.second = -1;
  	boostReco_bH.first = 9999;
  	boostReco_bH.second = -1;

   	if (zjets.size() > 0)
    	boostReco_bZ = doBoostedReco(ak4s, zjets.at(0).getP4(), 91.2, zllcand, 150.);
   	if (hjets.size() > 0)
    	boostReco_bH = doBoostedReco(ak4s, hjets.at(0).getP4(), 125., zllcand, 150.);

    if (ak4s.size() > 3){
    	if (zjets.size() == 0)
      	resReco_bZ = doResolvedReco(ak4s, 91.2, zllcand);
    	if (hjets.size() == 0)
    		resReco_bH = doResolvedReco(ak4s, 125., zllcand);
  	}
    pair<double, double> comboBZ, comboBH;

  	if (resReco_bZ.second > 0 && resReco_bZ.first < chiCut_){
    	h1_["resReco_bZ_1b"]->Fill(resReco_bZ.second, evtwt);
      h1_["st_bZ_1b"] -> Fill(ST, evtwt);
      if (!syst_ && !vv_ && !isData_) {
        for (unsigned i = 0; i < 101; i++) {
          h1_[Form("resReco_bZ_1b_pdf%d", i+1)] -> Fill(resReco_bZ.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
          h1_[Form("st_bZ_1b_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);

        }
        for (unsigned i = 0; i < 9; i++) {
          h1_[Form("resReco_bZ_1b_scale%d", i+1)] -> Fill(resReco_bZ.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
          h1_[Form("st_bZ_1b_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second);
        }
      }
  	} 

  	if (resReco_bH.second > 0 && resReco_bH.first < chiCut_){
    	h1_["resReco_bH_1b"]->Fill(resReco_bH.second, evtwt);
      h1_["st_bH_1b"] -> Fill(ST, evtwt);
      if (!syst_ && !vv_ && !isData_) {
       for (unsigned i = 0; i < 101; i++) {
         h1_[Form("resReco_bH_1b_pdf%d", i+1)] -> Fill(resReco_bH.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
         h1_[Form("st_bH_1b_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
       }
       for (unsigned i = 0; i < 9; i++) {
         h1_[Form("resReco_bH_1b_scale%d", i+1)] -> Fill(resReco_bH.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
         h1_[Form("st_bH_1b_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
       }
      }
  	}

  	if (boostReco_bZ.second > 0 && boostReco_bZ.first < chiCut_){
      h1_["boostReco_bZ"] -> Fill(boostReco_bZ.second, evtwt);
      h1_["st_bZ_boost"] -> Fill(ST, evtwt);
      if (!syst_ && !vv_ && !isData_) {
        for (unsigned i = 0; i < 101; i++) {
          h1_[Form("boostReco_bZ_pdf%d", i+1)] -> Fill(boostReco_bZ.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
          h1_[Form("st_bZ_boost_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
        }
        for (unsigned i = 0; i < 9; i++) {
          h1_[Form("boostReco_bZ_scale%d", i+1)] -> Fill(boostReco_bZ.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
          h1_[Form("st_bZ_boost_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
        }
      }
  	}

  	if (boostReco_bH.second > 0 && boostReco_bH.first < chiCut_){
      h1_["boostReco_bH"] -> Fill(boostReco_bH.second, evtwt);
      h1_["st_bH_boost"] -> Fill(ST, evtwt);
      if (!syst_ && !vv_ && !isData_) {
        for (unsigned i = 0; i < 101; i++) {
          h1_[Form("boostReco_bH_pdf%d", i+1)] -> Fill(boostReco_bH.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
          h1_[Form("st_bH_boost_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
        }
        for (unsigned i = 0; i < 9; i++) {
          h1_[Form("boostReco_bH_scale%d", i+1)] -> Fill(boostReco_bH.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
          h1_[Form("st_bH_boost_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
        }
      }
  	}


//    // Find if bZbZ or bZbH channel has better reconstruction for each category individually
//    if (zjets.size() > 0 && hjets.size() == 0) {
//      h1_["boostReco"] -> Fill(boostReco_bZ.second, evtwt);
//      if (!syst_ && !vv_ && !isData_) {
//        for (unsigned i = 0; i < 101; i++)
//          h1_[Form("boostReco_pdf%d", i+1)] -> Fill(boostReco_bZ.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
//        for (unsigned i = 0; i < 9; i++)
//          h1_[Form("boostReco_scale%d", i+1)] -> Fill(boostReco_bZ.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
//      }
//    }
//
//    else if (zjets.size() == 0 && hjets.size() > 0) {
//      h1_["boostReco"] -> Fill(boostReco_bH.second, evtwt);
//      if (!syst_ && !vv_ && !isData_) {
//        for (unsigned i = 0; i < 101; i++)
//          h1_[Form("boostReco_pdf%d", i+1)] -> Fill(boostReco_bH.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
//        for (unsigned i = 0; i < 9; i++)
//          h1_[Form("boostReco_scale%d", i+1)] -> Fill(boostReco_bH.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
//      }
//    }
//
//    else if (zjets.size() > 0 && hjets.size() > 0) {
//      if (boostReco_bZ.first < boostReco_bH.first) {
//        h1_["boostReco"] -> Fill(boostReco_bZ.second, evtwt);
//        if (!syst_ && !vv_ && !isData_) {
//          for (unsigned i = 0; i < 101; i++)
//            h1_[Form("boostReco_pdf%d", i+1)] -> Fill(boostReco_bZ.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
//          for (unsigned i = 0; i < 9; i++)
//            h1_[Form("boostReco_scale%d", i+1)] -> Fill(boostReco_bZ.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
//        }
//      }
//      else {
//        h1_["boostReco"] -> Fill(boostReco_bH.second, evtwt);
//        if (!syst_ && !vv_ && !isData_) {
//          for (unsigned i = 0; i < 101; i++)
//            h1_[Form("boostReco_pdf%d", i+1)] -> Fill(boostReco_bH.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
//          for (unsigned i = 0; i < 9; i++)
//            h1_[Form("boostReco_scale%d", i+1)] -> Fill(boostReco_bH.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
//        }
//      }
//    }
//
//    if (zjets.size() == 0 && hjets.size() > 0) {
//      h1_["resReco"] -> Fill(resReco_bZ.second, evtwt);
//        if (!syst_ && !vv_ && !isData_) {
//          for (unsigned i = 0; i < 101; i++)
//            h1_[Form("resReco_pdf%d", i+1)] -> Fill(resReco_bZ.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
//          for (unsigned i = 0; i < 9; i++)
//            h1_[Form("resReco_scale%d", i+1)] -> Fill(resReco_bZ.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
//        }
//    }
//    else if (zjets.size() > 0 && hjets.size() == 0) {
//      h1_["resReco"] -> Fill(resReco_bH.second, evtwt);
//        if (!syst_ && !vv_ && !isData_) {
//          for (unsigned i = 0; i < 101; i++)
//            h1_[Form("resReco_pdf%d", i+1)] -> Fill(resReco_bH.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
//          for (unsigned i = 0; i < 9; i++)
//            h1_[Form("resReco_scale%d", i+1)] -> Fill(resReco_bH.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
//        }
//    }
//    else if (zjets.size() == 0 && hjets.size() == 0) {
//      if (resReco_bZ.first < resReco_bH.first) {
//        h1_["resReco"] -> Fill(resReco_bZ.second, evtwt);
//        if (!syst_ && !vv_ && !isData_) {
//          for (unsigned i = 0; i < 101; i++)
//            h1_[Form("resReco_pdf%d", i+1)] -> Fill(resReco_bZ.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
//          for (unsigned i = 0; i < 9; i++)
//            h1_[Form("resReco_scale%d", i+1)] -> Fill(resReco_bZ.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
//        }
//      }
//      else {
//        h1_["resReco"] -> Fill(resReco_bH.second, evtwt);
//        if (!syst_ && !vv_ && !isData_) {
//          for (unsigned i = 0; i < 101; i++)
//            h1_[Form("resReco_pdf%d", i+1)] -> Fill(resReco_bH.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
//          for (unsigned i = 0; i < 9; i++)
//            h1_[Form("resReco_scale%d", i+1)] -> Fill(resReco_bH.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
//        }
//      }
//    }






//    // Find if bZbZ or bZbH channel has better reconstruction for each category individually
//    if (zjets.size() > 0 && hjets.size() == 0) h1_["boostReco"] -> Fill(boostReco_bZ.second, evtwt);
//    else if (zjets.size() == 0 && hjets.size() > 0) h1_["boostReco"] -> Fill(boostReco_bH.second, evtwt);
//    else if (zjets.size() > 0 && hjets.size() > 0) {
//      if (boostReco_bZ.first < boostReco_bH.first) h1_["boostReco"] -> Fill(boostReco_bZ.second, evtwt);
//      else h1_["boostReco"] -> Fill(boostReco_bH.second, evtwt);
//    }
//
//    if (zjets.size() == 0 && hjets.size() > 0) h1_["resReco"] -> Fill(resReco_bZ.second, evtwt);
//    else if (zjets.size() > 0 && hjets.size() == 0) h1_["resReco"] -> Fill(resReco_bH.second, evtwt);
//    else if (zjets.size() == 0 && hjets.size() == 0) {
//      if (resReco_bZ.first < resReco_bH.first) h1_["resReco"] -> Fill(resReco_bZ.second, evtwt);
//      else h1_["resReco"] -> Fill(resReco_bH.second, evtwt);
//    }
//
//    // Find if boosted or resolved category does better for each channel individually
//    if (boostReco_bZ.first < 9999 && boostReco_bZ.first > 0) comboBZ = boostReco_bZ;
//    else if(resReco_bZ.first < 9999 && resReco_bZ.first > 0) comboBZ = resReco_bZ;
//
//    if (boostReco_bH.first < 9999 && boostReco_bH.first > 0) comboBH = boostReco_bH;
//    else if (resReco_bH.first < 9999 && resReco_bH.first > 0) comboBH = resReco_bH;
//
//
//    if (comboBZ.first > 0) h1_["comboReco_bZ"] -> Fill(comboBZ.second, evtwt);
//    if (comboBH.first > 0) h1_["comboReco_bH"] -> Fill(comboBH.second, evtwt);
//
//    // Find best combination of channel and category
//    if (comboBZ.first <= 0 && comboBH.first > 0)  h1_["comboReco"] -> Fill(comboBH.second, evtwt);
//    else if (comboBZ.first > 0 && comboBH.first <= 0) h1_["comboReco"] -> Fill(comboBZ.second, evtwt);
//    else if (comboBZ.first < comboBH.first)  h1_["comboReco"] -> Fill(comboBZ.second, evtwt);
//    else if (comboBZ.first > comboBH.first)  h1_["comboReco"] -> Fill(comboBH.second, evtwt);
    
  }


  else if (ak4s.at(0).getPt() > 100 && ak4s.at(1).getPt() > 50 && bjets.size() >= 2 && ST > STMin_){


		pair<double, double> resReco_bZ_2b, boostReco_bZ_2b;
  	pair<double, double> resReco_bH_2b, boostReco_bH_2b;
 
 		resReco_bZ_2b.first = 9999;
  	resReco_bZ_2b.second = -1;
  	resReco_bH_2b.first = 9999;
  	resReco_bH_2b.second = -1;
  	boostReco_bZ_2b.first = 9999;
  	boostReco_bZ_2b.second = -1;
  	boostReco_bH_2b.first = 9999;
  	boostReco_bH_2b.second = -1;

   	if (zjets.size() > 0)
    	boostReco_bZ_2b = doBoostedReco(ak4s, zjets.at(0).getP4(), 91.2, zllcand, 150.);
   	if (hjets.size() > 0)
    	boostReco_bH_2b = doBoostedReco(ak4s, hjets.at(0).getP4(), 125., zllcand, 150.);

    if (ak4s.size() > 3){
    	if (zjets.size() == 0)
      	resReco_bZ_2b = doResolvedReco(ak4s, 91.2, zllcand);
    	if (hjets.size() == 0)
    		resReco_bH_2b = doResolvedReco(ak4s, 125., zllcand);
  	}

    pair<double, double> comboBZ_2b, comboBH_2b;

  	if (resReco_bZ_2b.second > 0 && resReco_bZ_2b.first < chiCut_){
    	h1_["resReco_bZ_2b"]->Fill(resReco_bZ_2b.second, evtwt);
      h1_["st_bZ_2b"]->Fill(ST, evtwt);
      if (!syst_ && !vv_ && !isData_) {
      for (unsigned i = 0; i < 101; i++) {
        h1_[Form("resReco_bZ_2b_pdf%d", i+1)] -> Fill(resReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
        h1_[Form("st_bZ_2b_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
      }
      for (unsigned i = 0; i < 9; i++) {
        h1_[Form("resReco_bZ_2b_scale%d", i+1)] -> Fill(resReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
        h1_[Form("st_bZ_2b_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
      }
      }
  	} 

  	if (resReco_bH_2b.second > 0 && resReco_bH_2b.first < chiCut_){
    	h1_["resReco_bH_2b"]->Fill(resReco_bH_2b.second, evtwt);
    	h1_["st_bH_2b"]->Fill(ST, evtwt);
      if (!syst_ && !vv_ && !isData_) {
       for (unsigned i = 0; i < 101; i++) {
         h1_[Form("resReco_bH_2b_pdf%d", i+1)] -> Fill(resReco_bH_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
         h1_[Form("st_bH_2b_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
       }
       for (unsigned i = 0; i < 9; i++) {
         h1_[Form("resReco_bH_2b_scale%d", i+1)] -> Fill(resReco_bH_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
         h1_[Form("st_bH_2b_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
       }
      }
  	}

  	if (boostReco_bZ_2b.second > 0 && boostReco_bZ_2b.first < chiCut_){
      h1_["boostReco_bZ"] -> Fill(boostReco_bZ_2b.second, evtwt);
      h1_["st_bZ_boost"] -> Fill(ST, evtwt);
      if (!syst_ && !vv_ && !isData_) {
       for (unsigned i = 0; i < 101; i++) {
         h1_[Form("boostReco_bZ_pdf%d", i+1)] -> Fill(boostReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
         h1_[Form("st_bZ_boost_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
       }
      for (unsigned i = 0; i < 9; i++) {
         h1_[Form("boostReco_bZ_scale%d", i+1)] -> Fill(boostReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
         h1_[Form("st_bZ_boost_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
       }
      }
  	}

  	if (boostReco_bH_2b.second > 0 && boostReco_bH_2b.first < chiCut_){
      h1_["boostReco_bH"] -> Fill(boostReco_bH_2b.second, evtwt);
      h1_["st_bH_boost"] -> Fill(ST, evtwt);
      if (!syst_ && !vv_ && !isData_) {
        for (unsigned i = 0; i < 101; i++) {
          h1_[Form("boostReco_bH_pdf%d", i+1)] -> Fill(boostReco_bH_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
         h1_[Form("st_bH_boost_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
        }
        for (unsigned i = 0; i < 9; i++) {
          h1_[Form("boostReco_bH_scale%d", i+1)] -> Fill(boostReco_bH_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
         h1_[Form("st_bH_boost_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
        }
      }
  	}

 //   // Find if bZbZ or bZbH channel has better reconstruction for each category individually
 //   if (zjets.size() > 0 && hjets.size() == 0) {
 //     h1_["boostReco_2b"] -> Fill(boostReco_bZ_2b.second, evtwt);
 //     if (!syst_ && !vv_ && !isData_) {
 //       for (unsigned i = 0; i < 101; i++)
 //         h1_[Form("boostReco_2b_pdf%d", i+1)] -> Fill(boostReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
 //       for (unsigned i = 0; i < 9; i++)
 //         h1_[Form("boostReco_2b_scale%d", i+1)] -> Fill(boostReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
 //     }
 //   }

 //   else if (zjets.size() == 0 && hjets.size() > 0) {
 //     h1_["boostReco_2b"] -> Fill(boostReco_bH_2b.second, evtwt);
 //     if (!syst_ && !vv_ && !isData_) {
 //       for (unsigned i = 0; i < 101; i++)
 //         h1_[Form("boostReco_2b_pdf%d", i+1)] -> Fill(boostReco_bH_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
 //       for (unsigned i = 0; i < 9; i++)
 //         h1_[Form("boostReco_2b_scale%d", i+1)] -> Fill(boostReco_bH_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
 //     }
 //   }

 //   else if (zjets.size() > 0 && hjets.size() > 0) {
 //     if (boostReco_bZ_2b.first < boostReco_bH_2b.first) {
 //       h1_["boostReco_2b"] -> Fill(boostReco_bZ_2b.second, evtwt);
 //       if (!syst_ && !vv_ && !isData_) {
 //         for (unsigned i = 0; i < 101; i++)
 //           h1_[Form("boostReco_2b_pdf%d", i+1)] -> Fill(boostReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
 //         for (unsigned i = 0; i < 9; i++)
 //           h1_[Form("boostReco_2b_scale%d", i+1)] -> Fill(boostReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
 //       }
 //     }
 //     else {
 //       h1_["boostReco_2b"] -> Fill(boostReco_bH_2b.second, evtwt);
 //       if (!syst_ && !vv_ && !isData_) {
 //         for (unsigned i = 0; i < 101; i++)
 //           h1_[Form("boostReco_2b_pdf%d", i+1)] -> Fill(boostReco_bH_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
 //         for (unsigned i = 0; i < 9; i++)
 //           h1_[Form("boostReco_2b_scale%d", i+1)] -> Fill(boostReco_bH_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
 //       }
 //     }
 //   }

 //   if (zjets.size() == 0 && hjets.size() > 0) {
 //     h1_["resReco_2b"] -> Fill(resReco_bZ_2b.second, evtwt);
 //       if (!syst_ && !vv_ && !isData_) {
 //         for (unsigned i = 0; i < 101; i++)
 //           h1_[Form("resReco_2b_pdf%d", i+1)] -> Fill(resReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
 //         for (unsigned i = 0; i < 9; i++)
 //           h1_[Form("resReco_2b_scale%d", i+1)] -> Fill(resReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
 //       }
 //   }
 //   else if (zjets.size() > 0 && hjets.size() == 0) {
 //     h1_["resReco_2b"] -> Fill(resReco_bH_2b.second, evtwt);
 //       if (!syst_ && !vv_ && !isData_) {
 //         for (unsigned i = 0; i < 101; i++)
 //           h1_[Form("resReco_2b_pdf%d", i+1)] -> Fill(resReco_bH_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
 //         for (unsigned i = 0; i < 9; i++)
 //           h1_[Form("resReco_2b_scale%d", i+1)] -> Fill(resReco_bH_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
 //       }
 //   }
 //   else if (zjets.size() == 0 && hjets.size() == 0) {
 //     if (resReco_bZ_2b.first < resReco_bH_2b.first) {
 //       h1_["resReco_2b"] -> Fill(resReco_bZ_2b.second, evtwt);
 //       if (!syst_ && !vv_ && !isData_) {
 //         for (unsigned i = 0; i < 101; i++)
 //           h1_[Form("resReco_2b_pdf%d", i+1)] -> Fill(resReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
 //         for (unsigned i = 0; i < 9; i++)
 //           h1_[Form("resReco_2b_scale%d", i+1)] -> Fill(resReco_bZ_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
 //       }
 //     }
 //     else {
 //       h1_["resReco_2b"] -> Fill(resReco_bH_2b.second, evtwt);
 //       if (!syst_ && !vv_ && !isData_) {
 //         for (unsigned i = 0; i < 101; i++)
 //           h1_[Form("resReco_2b_pdf%d", i+1)] -> Fill(resReco_bH_2b.second, evtwt*lhe_id_wts.at(i+pdfID_offset_).second);
 //         for (unsigned i = 0; i < 9; i++)
 //           h1_[Form("resReco_2b_scale%d", i+1)] -> Fill(resReco_bH_2b.second, evtwt*lhe_id_wts.at(i+scale_offset_).second);
 //       }
 //     }
 //   }
	}

	///////////////////////
	//Begin Optimization //	/////////////////////////////////////////////
	///////////////////////

 	if (optimizeReco_){

		GenParticleCollection genPartsInfo;
		genPartsInfo = genpart(evt);
		TLorentzVector bGen, bbarGen, q1, q2, ZGen, HGen;
		TLorentzVector bJet, bbarJet, qJet, qbarJet, ZJet, HJet;
		TLorentzVector had_bjet, lep_bjet, had_bgen, lep_bgen;

		for (auto& gen : genPartsInfo){
			if (gen.getPdgID() == 5 && abs(gen.getMom0PdgID()) == 8000002 && gen.getPt() != 0)
				bGen = gen.getP4();
			if (gen.getPdgID() == -5 && abs(gen.getMom0PdgID()) == 8000002 && gen.getPt() != 0)
				bbarGen = gen.getP4();
			if (signalType_ == "EvtType_MC_bZbZ"){
				if (gen.getPdgID() >= 1 && gen.getPdgID() <= 5 && abs(gen.getMom0PdgID()) == 23 && gen.getPt() != 0)
					q1 = gen.getP4();
				if (gen.getPdgID() <= -1 && gen.getPdgID() >= -5 && abs(gen.getMom0PdgID()) == 23 && gen.getPt() != 0)
					q2 = gen.getP4();
			}
      else if (signalType_ == "EvtType_MC_bZbH"){ 
 				if (gen.getPdgID() >= 1 && gen.getPdgID() <= 5 && abs(gen.getMom0PdgID()) == 25 && gen.getPt() != 0)
					q1 = gen.getP4();
				if (gen.getPdgID() <= -1 && gen.getPdgID() >= -5 && abs(gen.getMom0PdgID()) == 25 && gen.getPt() != 0)
					q2 = gen.getP4();
			}
			if (abs(gen.getPdgID()) == 23 && abs(gen.getMom0PdgID()) == 8000002 && gen.getPt() != 0)
				ZGen = gen.getP4();
			if (abs(gen.getPdgID()) == 25 && abs(gen.getMom0PdgID()) == 8000002 && gen.getPt() != 0)
				HGen = gen.getP4();
			
		}
		
		for (auto& jet : ak4s){
			if (jet.getP4().DeltaR(bGen) < 0.3 && jet.getPt() != 0)
				bJet = jet.getP4();
			if (jet.getP4().DeltaR(bbarGen) < 0.3 && jet.getPt() != 0)
				bbarJet = jet.getP4();
			if (jet.getP4().DeltaR(q1) < 0.3 && jet.getPt() != 0)
				qJet = jet.getP4();
			if (jet.getP4().DeltaR(q2) < 0.3 && jet.getPt() != 0)
				qbarJet = jet.getP4();
		}
		
		for (auto& jet : zjets){
			if (jet.getP4().DeltaR(ZGen) < 0.3 && jet.getPt() != 0)
				ZJet = jet.getP4();
		}

		for (auto& jet : hjets){
			if (jet.getP4().DeltaR(HGen) < 0.3 && jet.getPt() != 0)
				HJet = jet.getP4();
		}

		double vlqMass_ = 1000;
		double bcheck = abs((bJet+qJet+qbarJet).M() - vlqMass_);
		double bbarcheck = abs((bbarJet+qJet+qbarJet).M() - vlqMass_);

    if (bcheck < bbarcheck){
      had_bjet = bJet;
      lep_bjet = bbarJet;
      had_bgen = bGen;
      lep_bgen = bbarGen;
    }
    else {
      had_bjet = bbarJet;
      lep_bjet = bJet;
      had_bgen = bbarGen;
      lep_bgen = bGen;
    }
		if (qJet != qbarJet){
			if (signalType_ == "EvtType_MC_bZbZ")
				h1_["Zresolution"]->Fill((qJet+qbarJet).M(), evtwt);
			if (signalType_ == "EvtType_MC_bZbH") 
				h1_["Hresolution"]->Fill((qJet+qbarJet).M(), evtwt);	
			h1_["resolveRes"]->Fill((qJet+qbarJet+had_bjet).M(), evtwt);

		}

		if (zjets.size() > 0){
			h1_["ztagRes"]->Fill(ZJet.M(), evtwt);
  		h1_["boostZRes"]->Fill((ZJet+had_bjet).M(), evtwt);
    }
    if (hjets.size() > 0){
      h1_["htagRes"] -> Fill(HJet.M(), evtwt);
  		h1_["boostHRes"]->Fill((HJet+had_bjet).M(), evtwt);
    }

  	h1_["Blepresolution"]->Fill((zllcand+lep_bjet).M(), evtwt);
	}


  return true;
}

void MassReco::beginJob(){

  h1_["recoCutflow"] = fs->make<TH1D>("recoCutlow", "Count", 2, -0.5, 1.5);

  h1_["resReco_bZ_1b"] = fs->make<TH1D>("resReco_bZ_1b", "Resolved Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["resReco_bH_1b"] = fs->make<TH1D>("resReco_bH_1b", "Resolved Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco_bZ_1b"] = fs->make<TH1D>("boostReco_bZ_1b", "Boosted Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco_bH_1b"] = fs->make<TH1D>("boostReco_bH_1b", "Boosted Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);

    h1_["pt_z"]   = fs->make<TH1D>("pt_z", "", 1000, 0., 3000);
    h1_["pt_ak41"]= fs->make<TH1D>("pt_ak41", "", 1000, 0., 3000);
    h1_["pt_hjet"]= fs->make<TH1D>("pt_hjet", "", 1000, 0., 3000);
    h1_["pt_zjet"]= fs->make<TH1D>("pt_zjet", "", 1000, 0., 3000);
    h1_["st"]     = fs->make<TH1D>("st", "", 1000, 0., 3000);
    h1_["nhjet"]  = fs->make<TH1D>("nhjet", "", 10, 0., 100);
    h1_["nzjet"]  = fs->make<TH1D>("nzjet", "", 10, 0., 10);
    h1_["nbjet"]  = fs->make<TH1D>("nbjet", "", 10, 0., 10);
    h1_["nak4"]   = fs->make<TH1D>("nak4", "", 10, 0., 10);

  h1_["st_bZ_boost"] = fs->make<TH1D>("st_bZ_boost", "ST", 100, 0., 4000.);
  h1_["st_bH_boost"] = fs->make<TH1D>("st_bH_boost", "ST", 100, 0., 4000.);
  h1_["st_bZ_1b"] = fs->make<TH1D>("st_bZ_1b", "ST", 100, 0., 4000.);
  h1_["st_bH_1b"] = fs->make<TH1D>("st_bH_1b", "ST", 100, 0., 4000.);
  h1_["st_bZ_2b"] = fs->make<TH1D>("st_bZ_2b", "ST", 100, 0., 4000.);
  h1_["st_bH_2b"] = fs->make<TH1D>("st_bH_2b", "ST", 100, 0., 4000.);

  if (!syst_ && !vv_ && !isData_) {
  for (unsigned i = 0; i < 101; i++) {
    string st_bZ_boost_name = Form("st_bZ_boost_pdf%d", i+1);
    string st_bH_boost_name = Form("st_bH_boost_pdf%d", i+1);
    string st_bZ_1b_name = Form("st_bZ_1b_pdf%d", i+1);
    string st_bZ_2b_name = Form("st_bZ_2b_pdf%d", i+1);
    string st_bH_1b_name = Form("st_bH_1b_pdf%d", i+1);
    string st_bH_2b_name = Form("st_bH_2b_pdf%d", i+1);

    string resReco_bZ_1b_name = Form("resReco_bZ_1b_pdf%d", i+1);
    string resReco_bH_1b_name = Form("resReco_bH_1b_pdf%d", i+1);
    string resReco_bZ_2b_name = Form("resReco_bZ_2b_pdf%d", i+1);
    string resReco_bH_2b_name = Form("resReco_bH_2b_pdf%d", i+1);
    string boostReco_bZ_1b_name = Form("boostReco_bZ_1b_pdf%d", i+1);
    string boostReco_bH_1b_name = Form("boostReco_bH_1b_pdf%d", i+1);
    string boostReco_bZ_2b_name = Form("boostReco_bZ_2b_pdf%d", i+1);
    string boostReco_bH_2b_name = Form("boostReco_bH_2b_pdf%d", i+1);

    h1_[st_bZ_boost_name.c_str()] = fs->make<TH1D>(st_bZ_boost_name.c_str(), "reco", 100, 0., 4000.);
    h1_[st_bH_boost_name.c_str()] = fs->make<TH1D>(st_bH_boost_name.c_str(), "reco", 100, 0., 4000.);
    h1_[st_bZ_1b_name.c_str()] = fs->make<TH1D>(st_bZ_1b_name.c_str(), "reco", 100, 0., 4000.);
    h1_[st_bZ_2b_name.c_str()] = fs->make<TH1D>(st_bZ_2b_name.c_str(), "reco", 100, 0., 4000.);
    h1_[st_bH_1b_name.c_str()] = fs->make<TH1D>(st_bH_1b_name.c_str(), "reco", 100, 0., 4000.);
    h1_[st_bH_2b_name.c_str()] = fs->make<TH1D>(st_bH_2b_name.c_str(), "reco", 100, 0., 4000.);

    h1_[resReco_bZ_1b_name.c_str()] = fs->make<TH1D>(resReco_bZ_1b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[resReco_bH_1b_name.c_str()] = fs->make<TH1D>(resReco_bH_1b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[resReco_bZ_2b_name.c_str()] = fs->make<TH1D>(resReco_bZ_2b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[resReco_bH_2b_name.c_str()] = fs->make<TH1D>(resReco_bH_2b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bZ_1b_name.c_str()] = fs->make<TH1D>(boostReco_bZ_1b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bH_1b_name.c_str()] = fs->make<TH1D>(boostReco_bH_1b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bZ_2b_name.c_str()] = fs->make<TH1D>(boostReco_bZ_2b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bH_2b_name.c_str()] = fs->make<TH1D>(boostReco_bH_2b_name.c_str(), "reco", 1000, 0., 3000.);

    string resReco_bZ_name = Form("resReco_bZ_pdf%d", i+1);
    string resReco_bH_name = Form("resReco_bH_pdf%d", i+1);
    string boostReco_bZ_name = Form("boostReco_bZ_pdf%d", i+1);
    string boostReco_bH_name = Form("boostReco_bH_pdf%d", i+1);
    h1_[resReco_bZ_name.c_str()] = fs->make<TH1D>(resReco_bZ_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[resReco_bH_name.c_str()] = fs->make<TH1D>(resReco_bH_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bZ_name.c_str()] = fs->make<TH1D>(boostReco_bZ_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bH_name.c_str()] = fs->make<TH1D>(boostReco_bH_name.c_str(), "reco", 1000, 0., 3000.);

    string resReco_name = Form("resReco_pdf%d", i+1);
    string boostReco_name = Form("boostReco_pdf%d", i+1);
    h1_[resReco_name.c_str()] = fs->make<TH1D>(resReco_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_name.c_str()] = fs->make<TH1D>(boostReco_name.c_str(), "reco", 1000, 0., 3000.);

    string resReco_2b_name = Form("resReco_2b_pdf%d", i+1);
    string boostReco_2b_name = Form("boostReco_2b_pdf%d", i+1);
    h1_[resReco_2b_name.c_str()] = fs->make<TH1D>(resReco_2b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_2b_name.c_str()] = fs->make<TH1D>(boostReco_2b_name.c_str(), "reco", 1000, 0., 3000.);
  }

  for (unsigned i = 0; i < 9; i++) {

    string st_bZ_boost_name = Form("st_bZ_boost_scale%d", i+1);
    string st_bH_boost_name = Form("st_bH_boost_scale%d", i+1);

    string st_bZ_1b_name = Form("st_bZ_1b_scale%d", i+1);
    string st_bZ_2b_name = Form("st_bZ_2b_scale%d", i+1);
    string st_bH_1b_name = Form("st_bH_1b_scale%d", i+1);
    string st_bH_2b_name = Form("st_bH_2b_scale%d", i+1);
   
    string resReco_bZ_1b_name = Form("resReco_bZ_1b_scale%d", i+1);
    string resReco_bH_1b_name = Form("resReco_bH_1b_scale%d", i+1);
    string resReco_bZ_2b_name = Form("resReco_bZ_2b_scale%d", i+1);
    string resReco_bH_2b_name = Form("resReco_bH_2b_scale%d", i+1);
    string boostReco_bZ_1b_name = Form("boostReco_bZ_1b_scale%d", i+1);
    string boostReco_bH_1b_name = Form("boostReco_bH_1b_scale%d", i+1);
    string boostReco_bZ_2b_name = Form("boostReco_bZ_2b_scale%d", i+1);
    string boostReco_bH_2b_name = Form("boostReco_bH_2b_scale%d", i+1);

    h1_[st_bZ_boost_name.c_str()] = fs->make<TH1D>(st_bZ_boost_name.c_str(), "reco", 100, 0., 4000.);
    h1_[st_bH_boost_name.c_str()] = fs->make<TH1D>(st_bH_boost_name.c_str(), "reco", 100, 0., 4000.);

    h1_[st_bZ_1b_name.c_str()] = fs->make<TH1D>(st_bZ_1b_name.c_str(), "reco", 100, 0., 4000.);
    h1_[st_bZ_2b_name.c_str()] = fs->make<TH1D>(st_bZ_2b_name.c_str(), "reco", 100, 0., 4000.);
    h1_[st_bH_1b_name.c_str()] = fs->make<TH1D>(st_bH_1b_name.c_str(), "reco", 100, 0., 4000.);
    h1_[st_bH_2b_name.c_str()] = fs->make<TH1D>(st_bH_2b_name.c_str(), "reco", 100, 0., 4000.);

    h1_[resReco_bZ_1b_name.c_str()] = fs->make<TH1D>(resReco_bZ_1b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[resReco_bH_1b_name.c_str()] = fs->make<TH1D>(resReco_bH_1b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[resReco_bZ_2b_name.c_str()] = fs->make<TH1D>(resReco_bZ_2b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[resReco_bH_2b_name.c_str()] = fs->make<TH1D>(resReco_bH_2b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bZ_1b_name.c_str()] = fs->make<TH1D>(boostReco_bZ_1b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bH_1b_name.c_str()] = fs->make<TH1D>(boostReco_bH_1b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bZ_2b_name.c_str()] = fs->make<TH1D>(boostReco_bZ_2b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bH_2b_name.c_str()] = fs->make<TH1D>(boostReco_bH_2b_name.c_str(), "reco", 1000, 0., 3000.);

    string resReco_bZ_name = Form("resReco_bZ_scale%d", i+1);
    string resReco_bH_name = Form("resReco_bH_scale%d", i+1);
    string boostReco_bZ_name = Form("boostReco_bZ_scale%d", i+1);
    string boostReco_bH_name = Form("boostReco_bH_scale%d", i+1);
    h1_[resReco_bZ_name.c_str()] = fs->make<TH1D>(resReco_bZ_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[resReco_bH_name.c_str()] = fs->make<TH1D>(resReco_bH_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bZ_name.c_str()] = fs->make<TH1D>(boostReco_bZ_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_bH_name.c_str()] = fs->make<TH1D>(boostReco_bH_name.c_str(), "reco", 1000, 0., 3000.);

    string resReco_name = Form("resReco_scale%d", i+1);
    string boostReco_name = Form("boostReco_scale%d", i+1);
    h1_[resReco_name.c_str()] = fs->make<TH1D>(resReco_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_name.c_str()] = fs->make<TH1D>(boostReco_name.c_str(), "reco", 1000, 0., 3000.);

    string resReco_2b_name = Form("resReco_2b_scale%d", i+1);
    string boostReco_2b_name = Form("boostReco_2b_scale%d", i+1);
    h1_[resReco_2b_name.c_str()] = fs->make<TH1D>(resReco_2b_name.c_str(), "reco", 1000, 0., 3000.);
    h1_[boostReco_2b_name.c_str()] = fs->make<TH1D>(boostReco_2b_name.c_str(), "reco", 1000, 0., 3000.);

  }
  }

  h1_["resReco_bZ"] = fs->make<TH1D>("resReco_bZ", "Resolved Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["resReco_bH"] = fs->make<TH1D>("resReco_bH", "Resolved Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco_bZ"] = fs->make<TH1D>("boostReco_bZ", "Boosted Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco_bH"] = fs->make<TH1D>("boostReco_bH", "Boosted Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco_bZ"] = fs->make<TH1D>("comboReco_bZ", "Combo Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco_bH"] = fs->make<TH1D>("comboReco_bH", "Combo Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco"] = fs->make<TH1D>("comboReco", "Combo Reconstruction;M_{#chi^{2}}(B);;", 1000, 0., 3000);

  h1_["resReco"] = fs->make<TH1D>("resReco", "Resolved Reconstruction;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco"] = fs->make<TH1D>("boostReco", "Boosted Reconstruction;M_{#chi^{2}}(B);;", 1000, 0., 3000);

  h1_["resReco_bZ_2b"] = fs->make<TH1D>("resReco_bZ_2b", "Resolved Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["resReco_bH_2b"] = fs->make<TH1D>("resReco_bH_2b", "Resolved Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco_bZ_2b"] = fs->make<TH1D>("boostReco_bZ_2b", "Boosted Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco_bH_2b"] = fs->make<TH1D>("boostReco_bH_2b", "Boosted Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco_bZ_2b"] = fs->make<TH1D>("comboReco_bZ_2b", "Combo Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco_bH_2b"] = fs->make<TH1D>("comboReco_bH_2b", "Combo Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco_2b"] = fs->make<TH1D>("comboReco_2b", "Combo Reconstruction;M_{#chi^{2}}(B);;", 1000, 0., 3000);

  h1_["resReco_2b"] = fs->make<TH1D>("resReco_2b", "Resolved Reconstruction;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco_2b"] = fs->make<TH1D>("boostReco_2b", "Boosted Reconstruction;M_{#chi^{2}}(B);;", 1000, 0., 3000);

	h1_["Zresolution"] = fs->make<TH1D>("Zresolution", "Reco. Z Mass;M [GeV];;", 100, 20, 160);
	h1_["Hresolution"] = fs->make<TH1D>("Hresolution", "Reco. H Mass;M [GeV];;", 100, 65, 185);
	h1_["resolveRes"] = fs->make<TH1D>("resolveRes", "Had. Reco. B Mass;M [GeV];;", 200, 0, 2000);
	h1_["Blepresolution"] = fs->make<TH1D>("Blepresolution", "Lep. Reco. B Mass;M [GeV];;", 200, 0, 2000);
	h1_["ztagRes"] = fs->make<TH1D>("ztagRes", "Tagged Z Mass;M [GeV];;", 100, 0, 160);
	h1_["htagRes"] = fs->make<TH1D>("htagRes", "Tagged H Mass;M [GeV];;", 100, 65, 185);
	h1_["boostZRes"] = fs->make<TH1D>("boostZRes", "Had. Reco. B Mass;M [GeV];;", 200, 0, 2000);
	h1_["boostHRes"] = fs->make<TH1D>("boostHRes", "Had. Reco. B Mass;M [GeV];;", 200, 0, 2000);

}
void MassReco::endJob(){
	return;
}

pair<double, double> MassReco::vector_eval(vector<pair<double, double> > vec){
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

double MassReco::resolvedChi2(vector<TLorentzVector> jets, TLorentzVector Leptons, double bosMass, double mass){

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

double MassReco::boostedChi2(vector<TLorentzVector> ak4Jets, TLorentzVector ak8Jet, TLorentzVector Leptons, double bosMass, double mass, double pT){
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
pair<double, double> MassReco::doBoostedReco(vector<vlq::Jet> ak4Jets, TLorentzVector fatJet, double bosMass, TLorentzVector Leptons, double pT){
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

pair<double, double> MassReco::doResolvedReco(vector<vlq::Jet> collection, double bosMass, TLorentzVector Leptons){
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


DEFINE_FWK_MODULE(MassReco);


  
 
