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
  edm::Handle<double>            bTagwt_h   ; evt.getByToken(bTagwt_t, bTagwt_h);
  edm::Handle<double>            sjbTagwt_h ; evt.getByToken(sjbTagwt_t, sjbTagwt_h);

  TLorentzVector zllcand = (*zllcands_h.product()).at(0).getP4();
	double evtwt;
  if (evtwt_h.isValid())
   evtwt = *evtwt_h.product() * *bTagwt_h.product() * *sjbTagwt_h.product() ;
  else
    return false;
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

  bool doZ = true;
  bool doH = true;

	double chiCut_ = 100000;

	for (auto& jet : ak8s){
		h1_["ak8subjetiness"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["ak8prunedMass"]->Fill(jet.getPrunedMass(), evtwt);
    h1_["ak8softDropMass"] -> Fill(jet.getSoftDropMass(), evtwt);
	}
	for (auto& jet : zjets){
		h1_["Zsubjetiness"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["ZprunedMass"]->Fill(jet.getPrunedMass(), evtwt);
    h1_["ZsoftDropMass"] -> Fill(jet.getSoftDropMass(), evtwt);
	}
	for (auto& jet : hjets){
		h1_["Hsubjetiness"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["HprunedMass"]->Fill(jet.getPrunedMass(), evtwt);
    h1_["HsoftDropMass"] -> Fill(jet.getSoftDropMass(), evtwt);
	}

 	////////////////////////
	//Begin SR Mass Reco. //	///////////////////////////////////////////////////
	////////////////////////

	if (ak4s.at(0).getPt() > 100 && ak4s.at(1).getPt() > 50 && bjets.size() == 1 && ST > STMin_){

    h1_["recoCutflow"] -> Fill(1, evtwt);

	for (auto& jet : ak8s){
		h1_["ak8subjetiness_SR"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["ak8prunedMass_SR"]->Fill(jet.getPrunedMass(), evtwt);
    h1_["ak8softDropMass_SR"] -> Fill(jet.getSoftDropMass(), evtwt);
	}
	for (auto& jet : zjets){
		h1_["Zsubjetiness_SR"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["ZprunedMass_SR"]->Fill(jet.getPrunedMass(), evtwt);
    h1_["ZsoftDropMass_SR"] -> Fill(jet.getSoftDropMass(), evtwt);
	}
	for (auto& jet : hjets){
		h1_["Hsubjetiness_SR"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["HprunedMass_SR"]->Fill(jet.getPrunedMass(), evtwt);
    h1_["HsoftDropMass_SR"] -> Fill(jet.getSoftDropMass(), evtwt);
	}

		pair<double, double> resReco_bZ, boostReco_bZ, mergeReco_bZ;
  	pair<double, double> resReco_bH, boostReco_bH, mergeReco_bH;
 
 		resReco_bZ.first = 9999;
  	resReco_bZ.second = -1;
  	resReco_bH.first = 9999;
  	resReco_bH.second = -1;
  	boostReco_bZ.first = 9999;
  	boostReco_bZ.second = -1;
  	boostReco_bH.first = 9999;
  	boostReco_bH.second = -1;
  	mergeReco_bZ.first = 9999;
  	mergeReco_bZ.second = -1;
  	mergeReco_bH.first = 9999;
  	mergeReco_bH.second = -1;

   	if (zjets.size() > 0 && doZ)
    	boostReco_bZ = doBoostedReco(ak4s, zjets.at(0).getP4(), 91.2, zllcand, 150.);
   	if (hjets.size() > 0 && doH)
    	boostReco_bH = doBoostedReco(ak4s, hjets.at(0).getP4(), 125., zllcand, 150.);

  	if (ak4s.size() > 3){
    	if (zjets.size() == 0 && doZ)
      	resReco_bZ = doResolvedReco(ak4s, 91.2, zllcand);
    	if (hjets.size() == 0 && doH)
    		resReco_bH = doResolvedReco(ak4s, 125., zllcand);
  	}

  	for (unsigned i=0; i<ak4s.size(); i++){
    	if (ak4s.at(i).getP4().M() > 65 && ak4s.at(i).getP4().M() < 105 && ak4s.at(i).getP4().Pt() > 200 && doZ && ak4s.size() < 4 && zjets.size() == 0)
      	mergeReco_bZ = doBoostedReco(ak4s, ak4s.at(i).getP4(), 91.2, zllcand, 150.);
    	if (ak4s.at(i).getP4().M() > 105 && ak4s.at(i).getP4().M() < 135 && ak4s.at(i).getP4().Pt() > 200 && doH && ak4s.size() < 4 && hjets.size() == 0)
    	mergeReco_bH = doBoostedReco(ak4s, ak4s.at(i).getP4(), 125., zllcand, 150.);
  	}

    pair<double, double> comboBZ, comboBH;

  	if (resReco_bZ.second > 0 && resReco_bZ.first < chiCut_){

    	h1_["resReco_bZ"]->Fill(resReco_bZ.second, evtwt);
      h2_["resbZ_2d"]->Fill(resReco_bZ.first, resReco_bZ.second);
    	h1_["resST_bZ"]->Fill(ST, evtwt);
  	} 
  	if (resReco_bH.second > 0 && resReco_bH.first < chiCut_){

    	h1_["resReco_bH"]->Fill(resReco_bH.second, evtwt);
      h2_["resbH_2d"]->Fill(resReco_bH.first, resReco_bH.second);
    	h1_["resST_bZ"]->Fill(ST, evtwt);
  	}
  	if (boostReco_bZ.second > 0 && boostReco_bZ.first < chiCut_){

    	h1_["boostReco_bZ"]->Fill(boostReco_bZ.second, evtwt);
    	h2_["boostbZ_2d"]->Fill(boostReco_bZ.first, boostReco_bZ.second);
      h1_["boostST_bZ"]->Fill(ST, evtwt);
  	}
  	if (boostReco_bH.second > 0 && boostReco_bH.first < chiCut_){

    	h1_["boostReco_bH"]->Fill(boostReco_bH.second, evtwt);
    	h2_["boostbH_2d"]->Fill(boostReco_bH.first, boostReco_bH.second);
      h1_["boostST_bH"]->Fill(ST, evtwt);
  	}
  	if (mergeReco_bZ.second > 0 && mergeReco_bZ.first < chiCut_){
 
    	h1_["mergeReco_bZ"]->Fill(mergeReco_bZ.second, evtwt);
    	h2_["mergebZ_2d"]->Fill(mergeReco_bZ.first, mergeReco_bZ.second);
      h1_["mergeST_bZ"]->Fill(ST, evtwt);
  	}
  	if (mergeReco_bH.second > 0 && mergeReco_bH.first < chiCut_){
 
    	h1_["mergeReco_bH"]->Fill(mergeReco_bH.second, evtwt);
    	h2_["mergebH_2d"]->Fill(mergeReco_bH.first, mergeReco_bH.second);
      h1_["mergeST_bH"]->Fill(ST, evtwt);
  	}

  	if (resReco_bZ.first < boostReco_bZ.first && resReco_bZ.first < mergeReco_bZ.first && resReco_bZ.first < chiCut_){
      comboBZ = resReco_bZ;
    	h1_["comboReco_bZ"]->Fill(resReco_bZ.second, evtwt);
      h2_["combobZ_2d"]->Fill(resReco_bZ.first, resReco_bZ.second);
    	h1_["comboST_bZ"]->Fill(ST, evtwt);
  	}
  	else if (boostReco_bZ.first < resReco_bZ.first && boostReco_bZ.first < mergeReco_bZ.first && boostReco_bZ.first < chiCut_){
      comboBZ = boostReco_bZ;
    	h1_["comboReco_bZ"]->Fill(boostReco_bZ.second, evtwt);
      h2_["combobZ_2d"]->Fill(boostReco_bZ.first, boostReco_bZ.second);
    	h1_["comboST_bZ"]->Fill(ST, evtwt);
  	}
  	else if (mergeReco_bZ.first < resReco_bZ.first && mergeReco_bZ.first < boostReco_bZ.first && mergeReco_bZ.first < chiCut_){
      comboBZ = mergeReco_bZ;
    	h1_["comboReco_bZ"]->Fill(mergeReco_bZ.second, evtwt);
      h2_["combobZ_2d"]->Fill(mergeReco_bZ.first, mergeReco_bZ.second);
    	h1_["comboST_bZ"]->Fill(ST, evtwt);
  	}

  	if (resReco_bH.first < boostReco_bH.first && resReco_bH.first < mergeReco_bH.first && resReco_bH.first < chiCut_){
      comboBH = resReco_bH; 
    	h1_["comboReco_bH"]->Fill(resReco_bH.second, evtwt);
      h2_["combobH_2d"]->Fill(resReco_bH.first, resReco_bH.second);
    	h1_["comboST_bH"]->Fill(ST, evtwt);
  	}
  	else if (boostReco_bH.first < resReco_bH.first && boostReco_bH.first < mergeReco_bH.first && boostReco_bH.first < chiCut_){
      comboBH = boostReco_bH;
    	h1_["comboReco_bH"]->Fill(boostReco_bH.second, evtwt);
      h2_["combobH_2d"]->Fill(boostReco_bH.first, boostReco_bH.second);
    	h1_["comboST_bH"]->Fill(ST, evtwt);
  	}
  	else if (mergeReco_bH.first < resReco_bH.first && mergeReco_bH.first < boostReco_bH.first && mergeReco_bH.first < chiCut_){
      comboBH = mergeReco_bH;
    	h1_["comboReco_bH"]->Fill(mergeReco_bH.second, evtwt);
      h2_["combobH_2d"]->Fill(mergeReco_bH.first, mergeReco_bH.second);
    	h1_["comboST_bH"]->Fill(ST, evtwt);
  	}

    if (comboBZ.first <= 0 && comboBH.first > 0){
      h1_["bhReco_final"] -> Fill(comboBH.second, evtwt);
      h1_["comboReco"] -> Fill(comboBH.second, evtwt);
    }
    else if (comboBH.first <= 0 && comboBZ.first > 0){
      h1_["bzReco_final"] -> Fill(comboBZ.second,evtwt);
      h1_["comboReco"] -> Fill(comboBZ.second, evtwt);
    }
    else if (comboBZ.first < comboBH.first){
      h1_["bzReco_final"] -> Fill(comboBZ.second,evtwt);
      h1_["comboReco"] -> Fill(comboBZ.second, evtwt);
    }
    else if (comboBH.first < comboBZ.first){
      h1_["bhReco_final"] -> Fill(comboBH.second,evtwt);
      h1_["comboReco"] -> Fill(comboBH.second, evtwt);
    }

 	 	h1_["ST"]->Fill(ST, evtwt);

	}
  else if (ak4s.at(0).getPt() > 100 && ak4s.at(1).getPt() > 50 && bjets.size() >= 2 && ST > STMin_){

		pair<double, double> resReco_bZ_2b, boostReco_bZ_2b, mergeReco_bZ_2b;
  	pair<double, double> resReco_bH_2b, boostReco_bH_2b, mergeReco_bH_2b;
 
 		resReco_bZ_2b.first = 9999;
  	resReco_bZ_2b.second = -1;
  	resReco_bH_2b.first = 9999;
  	resReco_bH_2b.second = -1;
  	boostReco_bZ_2b.first = 9999;
  	boostReco_bZ_2b.second = -1;
  	boostReco_bH_2b.first = 9999;
  	boostReco_bH_2b.second = -1;
  	mergeReco_bZ_2b.first = 9999;
  	mergeReco_bZ_2b.second = -1;
  	mergeReco_bH_2b.first = 9999;
  	mergeReco_bH_2b.second = -1;

   	if (zjets.size() > 0 && doZ)
    	boostReco_bZ_2b = doBoostedReco(ak4s, zjets.at(0).getP4(), 91.2, zllcand, 150.);
   	if (hjets.size() > 0 && doH)
    	boostReco_bH_2b = doBoostedReco(ak4s, hjets.at(0).getP4(), 125., zllcand, 150.);

  	if (ak4s.size() > 3){
    	if (zjets.size() == 0 && doZ)
      	resReco_bZ_2b = doResolvedReco(ak4s, 91.2, zllcand);
    	if (hjets.size() == 0 && doH)
    		resReco_bH_2b = doResolvedReco(ak4s, 125., zllcand);
  	}

  	for (unsigned i=0; i<ak4s.size(); i++){
    	if (ak4s.at(i).getP4().M() > 65 && ak4s.at(i).getP4().M() < 105 && ak4s.at(i).getP4().Pt() > 200 && doZ && ak4s.size() < 4 && zjets.size() == 0)
      	mergeReco_bZ_2b = doBoostedReco(ak4s, ak4s.at(i).getP4(), 91.2, zllcand, 150.);
    	if (ak4s.at(i).getP4().M() > 105 && ak4s.at(i).getP4().M() < 135 && ak4s.at(i).getP4().Pt() > 200 && doH && ak4s.size() < 4 && hjets.size() == 0)
      	mergeReco_bH_2b = doBoostedReco(ak4s, ak4s.at(i).getP4(), 125., zllcand, 150.);
  	}

    pair<double, double> comboBZ_2b, comboBH_2b;

  	if (resReco_bZ_2b.second > 0 && resReco_bZ_2b.first < chiCut_)
    	h1_["resReco_bZ_2b"]->Fill(resReco_bZ_2b.second, evtwt); 

  	if (resReco_bH_2b.second > 0 && resReco_bH_2b.first < chiCut_)
    	h1_["resReco_bH_2b"]->Fill(resReco_bH_2b.second, evtwt);

  	if (boostReco_bZ_2b.second > 0 && boostReco_bZ_2b.first < chiCut_)
    	h1_["boostReco_bZ_2b"]->Fill(boostReco_bZ_2b.second, evtwt);

  	if (boostReco_bH_2b.second > 0 && boostReco_bH_2b.first < chiCut_)
    	h1_["boostReco_bH_2b"]->Fill(boostReco_bH_2b.second, evtwt);

  	if (mergeReco_bZ_2b.second > 0 && mergeReco_bZ_2b.first < chiCut_)
    	h1_["mergeReco_bZ_2b"]->Fill(mergeReco_bZ_2b.second, evtwt);

  	if (mergeReco_bH_2b.second > 0 && mergeReco_bH_2b.first < chiCut_)
    	h1_["mergeReco_bH_2b"]->Fill(mergeReco_bH_2b.second, evtwt);

  	if (resReco_bZ_2b.first < boostReco_bZ_2b.first && resReco_bZ_2b.first < mergeReco_bZ_2b.first && resReco_bZ_2b.first < chiCut_){
      comboBZ_2b = resReco_bZ_2b;
    	h1_["comboReco_bZ_2b"]->Fill(resReco_bZ_2b.second, evtwt);
  	}
  	else if (boostReco_bZ_2b.first < resReco_bZ_2b.first && boostReco_bZ_2b.first < mergeReco_bZ_2b.first && boostReco_bZ_2b.first < chiCut_){
      comboBZ_2b = boostReco_bZ_2b;
    	h1_["comboReco_bZ_2b"]->Fill(boostReco_bZ_2b.second, evtwt);
  	}
  	else if (mergeReco_bZ_2b.first < resReco_bZ_2b.first && mergeReco_bZ_2b.first < boostReco_bZ_2b.first && mergeReco_bZ_2b.first < chiCut_){
      comboBZ_2b = mergeReco_bZ_2b;
    	h1_["comboReco_bZ_2b"]->Fill(mergeReco_bZ_2b.second, evtwt);
  	}

  	if (resReco_bH_2b.first < boostReco_bH_2b.first && resReco_bH_2b.first < mergeReco_bH_2b.first && resReco_bH_2b.first < chiCut_){
      comboBH_2b = resReco_bH_2b; 
    	h1_["comboReco_bH_2b"]->Fill(resReco_bH_2b.second, evtwt);
  	}
  	else if (boostReco_bH_2b.first < resReco_bH_2b.first && boostReco_bH_2b.first < mergeReco_bH_2b.first && boostReco_bH_2b.first < chiCut_){
      comboBH_2b = boostReco_bH_2b;
    	h1_["comboReco_bH_2b"]->Fill(boostReco_bH_2b.second, evtwt);
  	}
  	else if (mergeReco_bH_2b.first < resReco_bH_2b.first && mergeReco_bH_2b.first < boostReco_bH_2b.first && mergeReco_bH_2b.first < chiCut_){
      comboBH_2b = mergeReco_bH_2b;
    	h1_["comboReco_bH_2b"]->Fill(mergeReco_bH_2b.second, evtwt);
  	}

    if (comboBZ_2b.first <= 0 && comboBH_2b.first > 0){
      h1_["bhReco_final_2b"] -> Fill(comboBH_2b.second, evtwt);
      h1_["comboReco_2b"] -> Fill(comboBH_2b.second, evtwt);
    }
    else if (comboBH_2b.first <= 0 && comboBZ_2b.first > 0){
      h1_["bzReco_final_2b"] -> Fill(comboBZ_2b.second,evtwt);
      h1_["comboReco_2b"] -> Fill(comboBZ_2b.second, evtwt);
    }
    else if (comboBZ_2b.first < comboBH_2b.first){
      h1_["bzReco_final_2b"] -> Fill(comboBZ_2b.second,evtwt);
      h1_["comboReco_2b"] -> Fill(comboBZ_2b.second, evtwt);
    }
    else if (comboBH_2b.first < comboBZ_2b.first){
      h1_["bhReco_final_2b"] -> Fill(comboBH_2b.second,evtwt);
      h1_["comboReco_2b"] -> Fill(comboBH_2b.second, evtwt);
    }

	}

	////////////////////////
	//Begin CR Mass Reco. //	///////////////////////////////////////////////////
	////////////////////////

	
	if (bjets.size() > 0 && ST < STMaxControl_ &&  controlReco_){

	for (auto& jet : ak8s){
		h1_["ak8subjetiness_CR"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["ak8prunedMass_CR"]->Fill(jet.getPrunedMass(), evtwt);
    h1_["ak8softDropMass_CR"] -> Fill(jet.getSoftDropMass(), evtwt);
	}
	for (auto& jet : zjets){
		h1_["Zsubjetiness_CR"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["ZprunedMass_CR"]->Fill(jet.getPrunedMass(), evtwt);
    h1_["ZsoftDropMass_CR"] -> Fill(jet.getSoftDropMass(), evtwt);
	}
	for (auto& jet : hjets){
		h1_["Hsubjetiness_CR"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["HprunedMass_CR"]->Fill(jet.getPrunedMass(), evtwt);
    h1_["HsoftDropMass_CR"] -> Fill(jet.getSoftDropMass(), evtwt);
	}


		pair<double, double> resCon_bZ, boostCon_bZ, mergeCon_bZ;
		pair<double, double> resCon_bH, boostCon_bH, mergeCon_bH;
					
		resCon_bZ.first = 9999;
		resCon_bZ.second = -1;
		resCon_bH.first = 9999;
		resCon_bH.second = -1;
		boostCon_bZ.first = 9999;
		boostCon_bZ.second = -1;
		boostCon_bH.first = 9999;
		boostCon_bH.second = -1;
		mergeCon_bZ.first = 9999;
		mergeCon_bZ.second = -1;
		mergeCon_bH.first = 9999;
		mergeCon_bH.second = -1;

		if (zjets.size() > 0 && doZ)
			boostCon_bZ = doBoostedReco(ak4s, zjets.at(0).getP4(), 91.2, zllcand, 150.);
		if (hjets.size() > 0 && doH)
			boostCon_bH = doBoostedReco(ak4s, hjets.at(0).getP4(), 125., zllcand, 150.);

		if (ak4s.size() > 3){
			if (zjets.size() == 0 && doZ)
				resCon_bZ = doResolvedReco(ak4s, 91.2, zllcand);
			if (hjets.size() == 0 && doH)
			resCon_bH = doResolvedReco(ak4s, 125., zllcand);
		}


		for (unsigned i=0; i<ak4s.size(); i++){
			if (ak4s.at(i).getP4().M() > 65 && ak4s.at(i).getP4().M() < 105 && ak4s.at(i).getP4().Pt() > 200 && doZ && ak4s.size() < 4 && zjets.size() == 0)
				mergeCon_bZ = doBoostedReco(ak4s, ak4s.at(i).getP4(), 91.2, zllcand, 150.);
			if (ak4s.at(i).getP4().M() > 105 && ak4s.at(i).getP4().M() < 135 && ak4s.at(i).getP4().Pt() > 200 && doH && ak4s.size() < 4 && hjets.size() == 0)
			mergeCon_bH = doBoostedReco(ak4s, ak4s.at(i).getP4(), 125., zllcand, 150.);
		}

		if (resCon_bZ.second > 0 && resCon_bZ.first < chiCut_){
			h1_["resCon_bZ"]->Fill(resCon_bZ.second, evtwt);
			h1_["resSTCon_bZ"]->Fill(ST, evtwt);
		} 
		if (resCon_bH.second > 0 && resCon_bH.first < chiCut_){ 
			h1_["resCon_bH"]->Fill(resCon_bH.second, evtwt);
			h1_["resSTCon_bZ"]->Fill(ST, evtwt);
		}
		if (boostCon_bZ.second > 0 && boostCon_bZ.first < chiCut_){
			h1_["boostCon_bZ"]->Fill(boostCon_bZ.second, evtwt);
			h1_["boostSTCon_bZ"]->Fill(ST, evtwt);
		}
		if (boostCon_bH.second > 0 && boostCon_bH.first < chiCut_){
			h1_["boostCon_bH"]->Fill(boostCon_bH.second, evtwt);
			h1_["boostSTCon_bH"]->Fill(ST, evtwt);
		}
		if (mergeCon_bZ.second > 0 && mergeCon_bZ.first < chiCut_){ 
			h1_["mergeCon_bZ"]->Fill(mergeCon_bZ.second, evtwt);
			h1_["mergeSTCon_bZ"]->Fill(ST, evtwt);
		}
		if (mergeCon_bH.second > 0 && mergeCon_bH.first < chiCut_){ 
			h1_["mergeCon_bH"]->Fill(mergeCon_bH.second, evtwt);
			h1_["mergeSTCon_bH"]->Fill(ST, evtwt);
		}

		if (resCon_bZ.first < boostCon_bZ.first && resCon_bZ.first < mergeCon_bZ.first && resCon_bZ.first < chiCut_){
			h1_["comboCon_bZ"]->Fill(resCon_bZ.second, evtwt);
			h1_["comboSTCon_bZ"]->Fill(ST, evtwt);
		}
		else if (boostCon_bZ.first < resCon_bZ.first && boostCon_bZ.first < mergeCon_bZ.first && boostCon_bZ.first < chiCut_){
			h1_["comboCon_bZ"]->Fill(boostCon_bZ.second, evtwt);
			h1_["comboSTCon_bZ"]->Fill(ST, evtwt);
		}
		else if (mergeCon_bZ.first < resCon_bZ.first && mergeCon_bZ.first < boostCon_bZ.first && mergeCon_bZ.first < chiCut_){
			h1_["comboCon_bZ"]->Fill(mergeCon_bZ.second, evtwt);
			h1_["comboSTCon_bZ"]->Fill(ST, evtwt);
		}

		if (resCon_bH.first < boostCon_bH.first && resCon_bH.first < mergeCon_bH.first && resCon_bH.first < chiCut_){
			h1_["comboCon_bH"]->Fill(resCon_bH.second, evtwt);
			h1_["comboSTCon_bH"]->Fill(ST, evtwt);
		}
		else if (boostCon_bH.first < resCon_bH.first && boostCon_bH.first < mergeCon_bH.first && boostCon_bH.first < chiCut_){
			h1_["comboCon_bH"]->Fill(boostCon_bH.second, evtwt);
			h1_["comboSTCon_bH"]->Fill(ST, evtwt);
		}
		else if (mergeCon_bH.first < resCon_bH.first && mergeCon_bH.first < boostCon_bH.first && mergeCon_bH.first < chiCut_){
			h1_["comboCon_bH"]->Fill(mergeCon_bH.second, evtwt);
			h1_["comboSTCon_bH"]->Fill(ST, evtwt);
		}

		h1_["STCon"]->Fill(ST, evtwt);

	}

	///////////////////////
	//Begin Optimization //	/////////////////////////////////////////////
	///////////////////////

 	if (optimizeReco_){

		GenParticleCollection genPartsInfo;
		genPartsInfo = genpart(evt);
		TLorentzVector bGen, bbarGen, q1, q2, ZGen, HGen;
		TLorentzVector Zmerge, Hmerge;
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
			if (jet.getP4().DeltaR(ZGen) < 0.3 && jet.getPt() != 0)
				Zmerge = jet.getP4();
			if (jet.getP4().DeltaR(HGen) < 0.3 && jet.getPt() != 0)
				Hmerge = jet.getP4();
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
      h1_["resolve_genBres"] -> Fill((q1 + q2 + had_bgen).M(), evtwt);

		}

		if (zjets.size() > 0){
			h1_["ztagRes"]->Fill(ZJet.M(), evtwt);
  		h1_["boostZRes"]->Fill((ZJet+had_bjet).M(), evtwt);
    }
    if (hjets.size() > 0){
      h1_["htagRes"] -> Fill(HJet.M(), evtwt);
  		h1_["boostHRes"]->Fill((HJet+had_bjet).M(), evtwt);
    }

    if (Zmerge.Pt() > 0){
	  	h1_["Zmerge"]->Fill(Zmerge.M(), evtwt);
      h1_["mergeZRes"] -> Fill((Zmerge + had_bjet).M(), evtwt);
    }
    if (Hmerge.Pt() > 0){
	  	h1_["Hmerge"]->Fill(Hmerge.M(), evtwt);
      h1_["mergeHRes"] -> Fill((Hmerge + had_bjet).M(), evtwt);
    }

    if (ZGen.Pt() > 0){
  		h1_["ztagGen"]->Fill(ZGen.M(), evtwt);
      h1_["boost_genBZres"] -> Fill((ZGen + had_bgen).M(), evtwt);
    }
    if (HGen.Pt() > 0){
		  h1_["htagGen"]->Fill(HGen.M(), evtwt);
      h1_["boost_genBHres"] -> Fill((HGen + had_bgen).M(), evtwt);
    }

  	h1_["Blepresolution"]->Fill((zllcand+lep_bjet).M(), evtwt);
	}


  return true;
}

void MassReco::beginJob(){

  h1_["recoCutflow"] = fs->make<TH1D>("recoCutlow", "Count", 2, -0.5, 1.5);


  h1_["bzReco_final"] = fs->make<TH1D>("bzReco_final","Combo Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000.);
  h1_["bhReco_final"] = fs->make<TH1D>("bhReco_final","Combo Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000.);
  h1_["resReco_bZ"] = fs->make<TH1D>("resReco_bZ", "Resolved Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["resReco_bH"] = fs->make<TH1D>("resReco_bH", "Resolved Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco_bZ"] = fs->make<TH1D>("boostReco_bZ", "Boosted Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco_bH"] = fs->make<TH1D>("boostReco_bH", "Boosted Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["mergeReco_bZ"] = fs->make<TH1D>("mergeReco_bZ", "Merged Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["mergeReco_bH"] = fs->make<TH1D>("mergeReco_bH", "Merged Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco_bZ"] = fs->make<TH1D>("comboReco_bZ", "Combo Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco_bH"] = fs->make<TH1D>("comboReco_bH", "Combo Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco"] = fs->make<TH1D>("comboReco", "Combo Reconstruction;M_{#chi^{2}}(B);;", 1000, 0., 3000);

  h1_["bzReco_final_2b"] = fs->make<TH1D>("bzReco_final_2b","Combo Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000.);
  h1_["bhReco_final_2b"] = fs->make<TH1D>("bhReco_final_2b","Combo Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000.);
  h1_["resReco_bZ_2b"] = fs->make<TH1D>("resReco_bZ_2b", "Resolved Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["resReco_bH_2b"] = fs->make<TH1D>("resReco_bH_2b", "Resolved Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco_bZ_2b"] = fs->make<TH1D>("boostReco_bZ_2b", "Boosted Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostReco_bH_2b"] = fs->make<TH1D>("boostReco_bH_2b", "Boosted Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["mergeReco_bZ_2b"] = fs->make<TH1D>("mergeReco_bZ_2b", "Merged Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["mergeReco_bH_2b"] = fs->make<TH1D>("mergeReco_bH_2b", "Merged Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco_bZ_2b"] = fs->make<TH1D>("comboReco_bZ_2b", "Combo Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco_bH_2b"] = fs->make<TH1D>("comboReco_bH_2b", "Combo Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboReco_2b"] = fs->make<TH1D>("comboReco_2b", "Combo Reconstruction;M_{#chi^{2}}(B);;", 1000, 0., 3000);

  h1_["resST_bZ"] = fs->make<TH1D>("resST_bZ", "Resolved ST B->bZ;S_{T};;", 1000, 0., 3000);
  h1_["resST_bH"] = fs->make<TH1D>("resST_bH", "Resolved ST B->bH;S_{T};;", 1000, 0., 3000);
  h1_["boostST_bZ"] = fs->make<TH1D>("boostST_bZ", "Boosted ST B->bZ;S_{T};;", 1000, 0., 3000);
  h1_["boostST_bH"] = fs->make<TH1D>("boostST_bH", "Boosted ST B->bH;S_{T};;", 1000, 0., 3000);
  h1_["mergeST_bZ"] = fs->make<TH1D>("mergeST_bZ", "Merged ST B->bZ;S_{T};;", 1000, 0., 3000);
  h1_["mergeST_bH"] = fs->make<TH1D>("mergeST_bH", "Merged ST B->bH;S_{T};;", 1000, 0., 3000);
  h1_["comboST_bZ"] = fs->make<TH1D>("comboST_bZ", "Combo ST B->bZ;S_{T};;", 1000, 0., 3000);
  h1_["comboST_bH"] = fs->make<TH1D>("comboST_bH", "Combo ST B->bH;S_{T};;", 1000, 0., 3000);

  h2_["resbZ_2d"] = fs->make<TH2D>("resbZ_d2", ";M_{#chi^{2}}(B);Min(#chi^{2});", 100, 0., 50., 400, 0, 2000);
  h2_["mergebZ_2d"] = fs->make<TH2D>("mergebZ_d2", ";M_{#chi^{2}}(B);Min(#chi^{2});", 100, 0., 50., 400, 0, 2000);
  h2_["boostbZ_2d"] = fs->make<TH2D>("boostbZ_d2", ";M_{#chi^{2}}(B);Min(#chi^{2});", 100, 0., 50., 400, 0, 2000);
  h2_["combobZ_2d"] = fs->make<TH2D>("combobZ_d2", ";M_{#chi^{2}}(B);Min(#chi^{2});", 100, 0., 50., 400, 0, 2000);
  h2_["resbH_2d"] = fs->make<TH2D>("resbH_d2", ";M_{#chi^{2}}(B);Min(#chi^{2});", 100, 0., 50., 400, 0, 2000);
  h2_["mergebH_2d"] = fs->make<TH2D>("mergebH_d2", ";M_{#chi^{2}}(B);Min(#chi^{2});", 100, 0., 50., 400, 0, 2000);
  h2_["boostbH_2d"] = fs->make<TH2D>("boostbH_d2", ";M_{#chi^{2}}(B);Min(#chi^{2});", 100, 0., 50., 400, 0, 2000);
  h2_["combobH_2d"] = fs->make<TH2D>("combobH_d2", ";M_{#chi^{2}}(B);Min(#chi^{2});", 100, 0., 50., 400, 0, 2000);

  h1_["ST"] = fs->make<TH1D>("ST", "ST;S_{T};;", 1000, 0., 8000.);

	h1_["Zresolution"] = fs->make<TH1D>("Zresolution", "Reco. Z Mass;M [GeV];;", 100, 20, 160);
	h1_["Hresolution"] = fs->make<TH1D>("Hresolution", "Reco. H Mass;M [GeV];;", 100, 65, 185);
	h1_["resolveRes"] = fs->make<TH1D>("resolveRes", "Had. Reco. B Mass;M [GeV];;", 200, 0, 2000);
  h1_["resolve_genBres"] = fs->make<TH1D>("resolve_genBres", "Gen. Had. Reco. B Mass;M [GeV];;", 200,0,2000);
	h1_["Blepresolution"] = fs->make<TH1D>("Blepresolution", "Lep. Reco. B Mass;M [GeV];;", 200, 0, 2000);
	h1_["ztagGen"] = fs->make<TH1D>("ztagGen", "Gen. Ztagged Jet Mass;M [GeV];;", 100, 20, 160);
  h1_["htagGen"] = fs->make<TH1D>("htagGen", "Gen. Htagged Jet Mass; M[GeV];;", 100, 65, 185);
	h1_["ztagRes"] = fs->make<TH1D>("ztagRes", "Tagged Z Mass;M [GeV];;", 100, 0, 160);
	h1_["htagRes"] = fs->make<TH1D>("htagRes", "Tagged H Mass;M [GeV];;", 100, 65, 185);
	h1_["boostZRes"] = fs->make<TH1D>("boostZRes", "Had. Reco. B Mass;M [GeV];;", 200, 0, 2000);
	h1_["boostHRes"] = fs->make<TH1D>("boostHRes", "Had. Reco. B Mass;M [GeV];;", 200, 0, 2000);
	h1_["boost_genBZres"] = fs->make<TH1D>("boost_genBZres", "Gen. Had. Reco. B Mass;M [GeV];;", 200, 0, 2000);
	h1_["boost_genBHres"] = fs->make<TH1D>("boost_genBHres", "Gen. Had. Reco. B Mass;M [GeV];;", 200, 0, 2000);
	h1_["Zmerge"] = fs->make<TH1D>("Zmerge", "Merged Z jet;M [GeV];;", 100, 20, 160);
	h1_["Hmerge"] = fs->make<TH1D>("Hmerge", "Merged H jet;M [GeV];;", 100, 65, 185);
  h1_["mergeZRes"] = fs->make<TH1D>("mergeZRes", "Gen. Merged Z jet;M [GeV];;", 200,0,2000);
  h1_["mergeHRes"] = fs->make<TH1D>("mergeHRes", "Gen. Merged H jet;M [GeV];;", 200,0,2000);

	h1_["ak8subjetiness"] = fs->make<TH1D>("ak8subjetiness", "AK8 Subjetiness;#tau_{21};;", 20, 0., 1.);
	h1_["ak8prunedMass"] = fs->make<TH1D>("ak8prunedMass", "AK8 Pruned Mass;M;;", 100, 0., 250.);
  h1_["ak8softDropMass"] = fs->make<TH1D>("ak8softDropMass", "AK8 Soft Drop Mass;M;;", 100, 0., 250.);
	h1_["Zsubjetiness"] = fs->make<TH1D>("Zsubjetiness","Z subjetiness;M;;", 10, 0, 5);
	h1_["Hsubjetiness"] = fs->make<TH1D>("Hsubjetiness","H subjetiness;M;;", 10, 0, 5);
	h1_["ZprunedMass"] = fs->make<TH1D>("ZprunedMass","Z Pruned Mass;M;;", 50, 30, 150);
	h1_["HprunedMass"] = fs->make<TH1D>("HprunedMass","H Pruned Mass;M;;", 50, 60, 180);
  h1_["HsoftDropMass"] = fs->make<TH1D>("HsoftDropMass", "H Soft Drop Mass;M;;", 100, 0., 250.);
  h1_["ZsoftDropMass"] = fs->make<TH1D>("ZsoftDropMass", "Z Soft Drop Mass;M;;", 100, 0., 250.);

	h1_["ak8subjetiness_SR"] = fs->make<TH1D>("ak8subjetiness_SR", "AK8 Subjetiness;#tau_{21};;", 20, 0., 1.);
	h1_["ak8prunedMass_SR"] = fs->make<TH1D>("ak8prunedMass_SR", "AK8 Pruned Mass;M;;", 100, 0., 250.);
  h1_["ak8softDropMass_SR"] = fs->make<TH1D>("ak8softDropMass_SR", "AK8 Soft Drop Mass;M;;", 100, 0., 250.);
	h1_["Zsubjetiness_SR"] = fs->make<TH1D>("Zsubjetiness_SR","Z subjetiness;#tau_{21}M;;", 30, 0, 5);
	h1_["Hsubjetiness_SR"] = fs->make<TH1D>("Hsubjetiness_SR","H subjetiness;#tau_{21};;", 30, 0, 5);
	h1_["ZprunedMass_SR"] = fs->make<TH1D>("ZprunedMass_SR","Z Pruned Mass;M;;", 50, 30, 150);
	h1_["HprunedMass_SR"] = fs->make<TH1D>("HprunedMass_SR","H Pruned Mass;M;;", 50, 60, 180);
  h1_["HsoftDropMass_SR"] = fs->make<TH1D>("HsoftDropMass_SR", "H Soft Drop Mass;M;;", 100, 0., 250.);
  h1_["ZsoftDropMass_SR"] = fs->make<TH1D>("ZsoftDropMass_SR", "Z Soft Drop Mass;M;;", 100, 0., 250.);


	h1_["ak8subjetiness_CR"] = fs->make<TH1D>("ak8subjetiness_CR", "AK8 Subjetiness;#tau_{21};;", 20, 0., 1.);
	h1_["ak8prunedMass_CR"] = fs->make<TH1D>("ak8prunedMass_CR", "AK8 Pruned Mass;M;;", 100, 0., 250.);
  h1_["ak8softDropMass_CR"] = fs->make<TH1D>("ak8softDropMass_CR", "AK8 Soft Drop Mass;M;;", 100, 0., 250.);
	h1_["Zsubjetiness_CR"] = fs->make<TH1D>("Zsubjetiness_CR","Z subjetiness;#tau_{21};;", 30, 0, 5);
	h1_["Hsubjetiness_CR"] = fs->make<TH1D>("Hsubjetiness_CR","H subjetiness;#tau_{21};;", 30, 0, 5);
	h1_["ZprunedMass_CR"] = fs->make<TH1D>("ZprunedMass_CR","Z Pruned Mass;M;;", 50, 30, 150);
	h1_["HprunedMass_CR"] = fs->make<TH1D>("HprunedMass_CR","H Pruned Mass;M;;", 50, 60, 180);
  h1_["HsoftDropMass_CR"] = fs->make<TH1D>("HsoftDropMass_CR", "H Soft Drop Mass;M;;", 100, 0., 250.);
  h1_["ZsoftDropMass_CR"] = fs->make<TH1D>("ZsoftDropMass_CR", "Z Soft Drop Mass;M;;", 100, 0., 250.);

  h1_["resCon_bZ"] = fs->make<TH1D>("resCon_bZ", "Resolved Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["resCon_bH"] = fs->make<TH1D>("resCon_bH", "Resolved Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostCon_bZ"] = fs->make<TH1D>("boostCon_bZ", "Boosted Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["boostCon_bH"] = fs->make<TH1D>("boostCon_bH", "Boosted Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["mergeCon_bZ"] = fs->make<TH1D>("mergeCon_bZ", "Merged Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["mergeCon_bH"] = fs->make<TH1D>("mergeCon_bH", "Merged Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboCon_bZ"] = fs->make<TH1D>("comboCon_bZ", "Combo Reconstruction B->bZ;M_{#chi^{2}}(B);;", 1000, 0., 3000);
  h1_["comboCon_bH"] = fs->make<TH1D>("comboCon_bH", "Combo Reconstruction B->bH;M_{#chi^{2}}(B);;", 1000, 0., 3000);

  h1_["resSTCon_bZ"] = fs->make<TH1D>("resSTCon_bZ", "Resolved ST B->bZ;S_{T};;", 1000, 0., 3000);
  h1_["resSTCon_bH"] = fs->make<TH1D>("resSTCon_bH", "Resolved ST B->bH;S_{T};;", 1000, 0., 3000);
  h1_["boostSTCon_bZ"] = fs->make<TH1D>("boostSTCon_bZ", "Boosted ST B->bZ;S_{T};;", 1000, 0., 3000);
  h1_["boostSTCon_bH"] = fs->make<TH1D>("boostSTCon_bH", "Boosted ST B->bH;S_{T};;", 1000, 0., 3000);
  h1_["mergeSTCon_bZ"] = fs->make<TH1D>("mergeSTCon_bZ", "Merged ST B->bZ;S_{T};;", 1000, 0., 3000);
  h1_["mergeSTCon_bH"] = fs->make<TH1D>("mergeSTCon_bH", "Merged ST B->bH;S_{T};;", 1000, 0., 3000);
  h1_["comboSTCon_bZ"] = fs->make<TH1D>("comboSTCon_bZ", "Combo ST B->bZ;S_{T};;", 1000, 0., 3000);
  h1_["comboSTCon_bH"] = fs->make<TH1D>("comboSTCon_bH", "Combo ST B->bH;S_{T};;", 1000, 0., 3000);

  h1_["STCon"] = fs->make<TH1D>("STCon", "Sum Pt;S_{T};;", 1000, 0., 8000.);


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

  for (int mass = 0; mass <= 2000; mass+=5){

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


  
 
