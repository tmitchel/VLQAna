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
  edm::EDGetTokenT<double>             evtwt_t;
  edm::EDGetTokenT<double>             st_t;

  const double ptMin_;
	const double STMaxControl_;
	const double STMin_;
  const string zdecaymode_;
  const string signalType_;
	const bool optimizeReco_;
	const bool controlReco_;
	PickGenPart genpart ;

  edm::Service<TFileService> fs;
  map<string, TH1D*> h1_;

};

MassReco::MassReco(const edm::ParameterSet& iConfig) :
	ak8_t         (consumes<vector<vlq::Jet> > (iConfig.getParameter<edm::InputTag>("ak8jets"))),
  ak4_t         (consumes<vector<vlq::Jet> > (iConfig.getParameter<edm::InputTag>("jets"))),
  bjets_t       (consumes<vector<vlq::Jet> > (iConfig.getParameter<edm::InputTag>("bjets"))),
  zjets_t       (consumes<vector<vlq::Jet> > (iConfig.getParameter<edm::InputTag>("zjets"))),
  hjets_t       (consumes<vector<vlq::Jet> > (iConfig.getParameter<edm::InputTag>("hjets"))),
  zllcands_t    (consumes<vector<vlq::Candidate> > (iConfig.getParameter<edm::InputTag>("zllcands"))),
  evtwt_t       (consumes<double>           (iConfig.getParameter<edm::InputTag>("evtwt"))),
  st_t          (consumes<double>           (iConfig.getParameter<edm::InputTag>("st"))),
 
  ptMin_     (iConfig.getParameter<double> ("ptMin")),
	STMaxControl_ (iConfig.getParameter<double> ("STMaxControl")),
	STMin_     (iConfig.getParameter<double> ("STMin")),
  zdecaymode_       (iConfig.getParameter<string> ("zdecaymode")),
  signalType_ (iConfig.getParameter<string> ("signalType")),
	optimizeReco_ (iConfig.getParameter<bool> ("optimizeReco")),
	controlReco_  (iConfig.getParameter<bool> ("controlReco")),
	genpart    (iConfig.getParameter<edm::ParameterSet>("genParams"),consumesCollector())
{}

MassReco::~MassReco() {}

bool MassReco::filter(edm::Event& evt, const edm::EventSetup& iSetup) {

  //TLorentzVector leptons = (lepP4.at(0) + lepP4.at(1));
  //  cout << lepP4.at(0).E() << " " << lepP4.at(1).E() << endl;
	
	edm::Handle<vector<vlq::Candidate> > zllcands_h; evt.getByToken(zllcands_t, zllcands_h);
	edm::Handle<vector<vlq::Jet> > ak8_h      ; evt.getByToken(ak8_t,  ak8_h)     ;
	edm::Handle<vector<vlq::Jet> > ak4_h      ; evt.getByToken(ak4_t,  ak4_h)     ;
  edm::Handle<vector<vlq::Jet> > bjets_h    ; evt.getByToken(bjets_t,  bjets_h) ;
  edm::Handle<vector<vlq::Jet> > zjets_h    ; evt.getByToken(zjets_t,  zjets_h) ;
  edm::Handle<vector<vlq::Jet> > hjets_h    ; evt.getByToken(hjets_t,  hjets_h) ;
  edm::Handle<double>              evtwt_h  ; evt.getByToken(evtwt_t,  evtwt_h) ;
  edm::Handle<double>              st_h     ; evt.getByToken(st_t,     st_h)    ;
	TLorentzVector zllcand = (*zllcands_h.product()).at(0).getP4();
	double evtwt = *evtwt_h.product();
	double ST = *st_h.product();
	vector<vlq::Jet> ak8s = *ak8_h.product();
  vector<vlq::Jet> ak4s = *ak4_h.product();
  vector<vlq::Jet> bjets = *bjets_h.product();
  vector<vlq::Jet> zjets = *zjets_h.product();
  vector<vlq::Jet> hjets = *hjets_h.product();

	bool doZ = (signalType_ == "EvtType_MC_bZbZ" || signalType_ == "");
  bool doH = (signalType_ == "EvtType_MC_bZbH" || signalType_ == "");

	double chiCut_ = 20;

	for (auto& jet : ak8s){
		h1_["ak8subjetiness"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["ak8prunedMass"]->Fill(jet.getPrunedMass(), evtwt);
	}
	for (auto& jet : zjets){
		h1_["Zsubjetiness"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["ZprunedMass"]->Fill(jet.getPrunedMass(), evtwt);
	}
	for (auto& jet : hjets){
		h1_["Hsubjetiness"]->Fill(jet.getTau2()/jet.getTau1(), evtwt);
		h1_["HprunedMass"]->Fill(jet.getPrunedMass(), evtwt);
	}

 	////////////////////////
	//Begin SR Mass Reco. //	///////////////////////////////////////////////////
	////////////////////////

	if (ak4s.at(0).getPt() > 100 && ak4s.at(1).getPt() > 50 && bjets.size() > 0 && ST > STMin_){

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

  	if (resReco_bZ.second > 0 && resReco_bZ.first < chiCut_){
    	h1_["resReco_bZ"]->Fill(resReco_bZ.second, evtwt);
    	h1_["resST_bZ"]->Fill(ST, evtwt);
  	} 
  	if (resReco_bH.second > 0 && resReco_bH.first < chiCut_){ 
    	h1_["resReco_bH"]->Fill(resReco_bH.second, evtwt);
    	h1_["resST_bZ"]->Fill(ST, evtwt);
  	}
  	if (boostReco_bZ.second > 0 && boostReco_bZ.first < chiCut_){
    	h1_["boostReco_bZ"]->Fill(boostReco_bZ.second, evtwt);
    	h1_["boostST_bZ"]->Fill(ST, evtwt);
  	}
  	if (boostReco_bH.second > 0 && boostReco_bH.first < chiCut_){
    	h1_["boostReco_bH"]->Fill(boostReco_bH.second, evtwt);
    	h1_["boostST_bH"]->Fill(ST, evtwt);
  	}
  	if (mergeReco_bZ.second > 0 && mergeReco_bZ.first < chiCut_){ 
    	h1_["mergeReco_bZ"]->Fill(mergeReco_bZ.second, evtwt);
    	h1_["mergeST_bZ"]->Fill(ST, evtwt);
  	}
  	if (mergeReco_bH.second > 0 && mergeReco_bH.first < chiCut_){ 
    	h1_["mergeReco_bH"]->Fill(mergeReco_bH.second, evtwt);
    	h1_["mergeST_bH"]->Fill(ST, evtwt);
  	}

  	if (resReco_bZ.first < boostReco_bZ.first && resReco_bZ.first < mergeReco_bZ.first && resReco_bZ.first < chiCut_){
    	h1_["comboReco_bZ"]->Fill(resReco_bZ.second, evtwt);
    	h1_["comboST_bZ"]->Fill(ST, evtwt);
  	}
  	else if (boostReco_bZ.first < resReco_bZ.first && boostReco_bZ.first < mergeReco_bZ.first && boostReco_bZ.first < chiCut_){
    	h1_["comboReco_bZ"]->Fill(boostReco_bZ.second, evtwt);
    	h1_["comboST_bZ"]->Fill(ST, evtwt);
  	}
  	else if (mergeReco_bZ.first < resReco_bZ.first && mergeReco_bZ.first < boostReco_bZ.first && mergeReco_bZ.first < chiCut_){
    	h1_["comboReco_bZ"]->Fill(mergeReco_bZ.second, evtwt);
    	h1_["comboST_bZ"]->Fill(ST, evtwt);
  	}

  	if (resReco_bH.first < boostReco_bH.first && resReco_bH.first < mergeReco_bH.first && resReco_bH.first < chiCut_){
    	h1_["comboReco_bH"]->Fill(resReco_bH.second, evtwt);
    	h1_["comboST_bH"]->Fill(ST, evtwt);
  	}
  	else if (boostReco_bH.first < resReco_bH.first && boostReco_bH.first < mergeReco_bH.first && boostReco_bH.first < chiCut_){
    	h1_["comboReco_bH"]->Fill(boostReco_bH.second, evtwt);
    	h1_["comboST_bH"]->Fill(ST, evtwt);
  	}
  	else if (mergeReco_bH.first < resReco_bH.first && mergeReco_bH.first < boostReco_bH.first && mergeReco_bH.first < chiCut_){
    	h1_["comboReco_bH"]->Fill(mergeReco_bH.second, evtwt);
    	h1_["comboST_bH"]->Fill(ST, evtwt);
  	}

 	 	h1_["ST"]->Fill(ST, evtwt);

	}
	////////////////////////
	//Begin CR Mass Reco. //	///////////////////////////////////////////////////
	////////////////////////

	
	if (bjets.size() > 0 && ST < STMaxControl_ &&  controlReco_){

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
			if (signalType_ == "EvtType_MC_bZbH"){ 
 				if (gen.getPdgID() >= 1 && gen.getPdgID() <= 5 && abs(gen.getMom0PdgID()) == 25 && gen.getPt() != 0)
					q1 = gen.getP4();
				if (gen.getPdgID() <= -1 && gen.getPdgID() >= -5 && abs(gen.getMom0PdgID()) == 25 && gen.getPt() != 0)
					q2 = gen.getP4();
			}
			if (abs(gen.getPdgID()) == 23 && abs(gen.getMom0PdgID()) == 8000002)// && abs(gen.getDau0PdgID()) >=1 && abs(gen.getDau0PdgID()) <= 5 && gen.getPt() != 0)
				ZGen = gen.getP4();
			if (abs(gen.getPdgID()) == 25 && abs(gen.getMom0PdgID()) == 8000002)// && abs(gen.getDau0PdgID()) >=1 && abs(gen.getDau0PdgID()) <= 5 && gen.getPt() != 0)
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

		double vlqMass_ = 1200;
		double bcheck = abs((bJet+qJet+qbarJet).M() - vlqMass_);
		double bbarcheck = abs((bbarJet+qJet+qbarJet).M() - vlqMass_);

    if (bcheck < bbarcheck){
      had_bjet = bJet;
      lep_bjet = bbarJet;
    }
    else {
      had_bjet = bbarJet;
      lep_bjet = bJet;
    }
		if (qJet != qbarJet){
			if (signalType_ == "EvtType_MC_bZbZ")
				h1_["Zresolution"]->Fill((qJet+qbarJet).M(), evtwt);
			if (signalType_ == "EvtType_MC_bZbH") 
				h1_["Hresolution"]->Fill((qJet+qbarJet).M(), evtwt);	
			h1_["Bhadresolution"]->Fill((qJet+qbarJet+had_bjet).M(), evtwt);
			h1_["Blepresolution"]->Fill((zllcand+lep_bjet).M(), evtwt);

		}
		if (zjets.size() > 0)
			h1_["ztagPlain"]->Fill(ZJet.M(), evtwt);
		h1_["ztagRes"]->Fill(ZGen.M(), evtwt);
		h1_["htagRes"]->Fill(HJet.M(), evtwt);
		h1_["ztagRecores"]->Fill((ZJet+had_bjet).M(), evtwt);
		h1_["htagRecores"]->Fill((HJet+had_bjet).M(), evtwt);
		h1_["Zmerge"]->Fill(Zmerge.M(), evtwt);
		h1_["Hmerge"]->Fill(Hmerge.M(), evtwt);
	}


  return true;
}

void MassReco::beginJob(){
  h1_["resReco_bZ"] = fs->make<TH1D>("resReco_bZ", "Resolved Reconstruction B->bZ", 1000, 0., 3000);
  h1_["resReco_bH"] = fs->make<TH1D>("resReco_bH", "Resolved Reconstruction B->bH", 1000, 0., 3000);
  h1_["boostReco_bZ"] = fs->make<TH1D>("boostReco_bZ", "Boosted Reconstruction B->bZ", 1000, 0., 3000);
  h1_["boostReco_bH"] = fs->make<TH1D>("boostReco_bH", "Boosted Reconstruction B->bH", 1000, 0., 3000);
  h1_["mergeReco_bZ"] = fs->make<TH1D>("mergeReco_bZ", "Merged Reconstruction B->bZ", 1000, 0., 3000);
  h1_["mergeReco_bH"] = fs->make<TH1D>("mergeReco_bH", "Merged Reconstruction B->bH", 1000, 0., 3000);
  h1_["comboReco_bZ"] = fs->make<TH1D>("comboReco_bZ", "Combo Reconstruction B->bZ", 1000, 0., 3000);
  h1_["comboReco_bH"] = fs->make<TH1D>("comboReco_bH", "Combo Reconstruction B->bH", 1000, 0., 3000);

  h1_["resST_bZ"] = fs->make<TH1D>("resST_bZ", "Resolved ST B->bZ", 1000, 0., 3000);
  h1_["resST_bH"] = fs->make<TH1D>("resST_bH", "Resolved ST B->bH", 1000, 0., 3000);
  h1_["boostST_bZ"] = fs->make<TH1D>("boostST_bZ", "Boosted ST B->bZ", 1000, 0., 3000);
  h1_["boostST_bH"] = fs->make<TH1D>("boostST_bH", "Boosted ST B->bH", 1000, 0., 3000);
  h1_["mergeST_bZ"] = fs->make<TH1D>("mergeST_bZ", "Merged ST B->bZ", 1000, 0., 3000);
  h1_["mergeST_bH"] = fs->make<TH1D>("mergeST_bH", "Merged ST B->bH", 1000, 0., 3000);
  h1_["comboST_bZ"] = fs->make<TH1D>("comboST_bZ", "Combo ST B->bZ", 1000, 0., 3000);
  h1_["comboST_bH"] = fs->make<TH1D>("comboST_bH", "Combo ST B->bH", 1000, 0., 3000);

  h1_["ST"] = fs->make<TH1D>("ST", "Sum Pt", 1000, 0., 8000.);

	h1_["Zresolution"] = fs->make<TH1D>("Zresolution", "reco Z resolution", 100, 20, 160);
	h1_["Hresolution"] = fs->make<TH1D>("Hresolution", "reco H resolution", 100, 65, 185);
	h1_["Bhadresolution"] = fs->make<TH1D>("Bhadresolution", "had. reco. B resolution", 200, 0, 2000);
	h1_["Blepresolution"] = fs->make<TH1D>("Blepresolution", "lep. reco B resolution", 200, 0, 2000);
	h1_["ztagPlain"] = fs->make<TH1D>("ztagPlain", "ztagged jet mass", 100, 20, 160);
	h1_["ztagRes"] = fs->make<TH1D>("ztagRes", "tagged Z resolution", 100, 0, 160);
	h1_["htagRes"] = fs->make<TH1D>("htagRes", "tagged H resolution", 100, 65, 185);
	h1_["ztagRecores"] = fs->make<TH1D>("ztagRecores", "had reco B resolution w/ Ztag", 200, 0, 2000);
	h1_["htagRecores"] = fs->make<TH1D>("htagRecores", "had reco B resolution w/ Htag", 200, 0, 2000);
	h1_["Zmerge"] = fs->make<TH1D>("Zmerge", "merged Z jet", 100, 20, 160);
	h1_["Hmerge"] = fs->make<TH1D>("Hmerge", "merged H jet", 100, 65, 185);

	h1_["ak8subjetiness"] = fs->make<TH1D>("ak8subjetiness", "tau2/1", 20, 0., 1.);
	h1_["ak8prunedMass"] = fs->make<TH1D>("ak8prunedMass", "Pruned M", 100, 0., 250.);
	h1_["Zsubjetiness"] = fs->make<TH1D>("Zsubjetiness","Z subjetiness", 10, 0, 5);
	h1_["Hsubjetiness"] = fs->make<TH1D>("Hsubjetiness","H subjetiness", 10, 0, 5);
	h1_["ZprunedMass"] = fs->make<TH1D>("ZprunedMass","Z pruned Mass", 50, 30, 150);
	h1_["HprunedMass"] = fs->make<TH1D>("HprunedMass","H pruned Mass", 50, 60, 180);

  h1_["resCon_bZ"] = fs->make<TH1D>("resCon_bZ", "Resolved Reconstruction B->bZ", 1000, 0., 3000);
  h1_["resCon_bH"] = fs->make<TH1D>("resCon_bH", "Resolved Reconstruction B->bH", 1000, 0., 3000);
  h1_["boostCon_bZ"] = fs->make<TH1D>("boostCon_bZ", "Boosted Reconstruction B->bZ", 1000, 0., 3000);
  h1_["boostCon_bH"] = fs->make<TH1D>("boostCon_bH", "Boosted Reconstruction B->bH", 1000, 0., 3000);
  h1_["mergeCon_bZ"] = fs->make<TH1D>("mergeCon_bZ", "Merged Reconstruction B->bZ", 1000, 0., 3000);
  h1_["mergeCon_bH"] = fs->make<TH1D>("mergeCon_bH", "Merged Reconstruction B->bH", 1000, 0., 3000);
  h1_["comboCon_bZ"] = fs->make<TH1D>("comboCon_bZ", "Combo Reconstruction B->bZ", 1000, 0., 3000);
  h1_["comboCon_bH"] = fs->make<TH1D>("comboCon_bH", "Combo Reconstruction B->bH", 1000, 0., 3000);

  h1_["resSTCon_bZ"] = fs->make<TH1D>("resSTCon_bZ", "Resolved ST B->bZ", 1000, 0., 3000);
  h1_["resSTCon_bH"] = fs->make<TH1D>("resSTCon_bH", "Resolved ST B->bH", 1000, 0., 3000);
  h1_["boostSTCon_bZ"] = fs->make<TH1D>("boostSTCon_bZ", "Boosted ST B->bZ", 1000, 0., 3000);
  h1_["boostSTCon_bH"] = fs->make<TH1D>("boostSTCon_bH", "Boosted ST B->bH", 1000, 0., 3000);
  h1_["mergeSTCon_bZ"] = fs->make<TH1D>("mergeSTCon_bZ", "Merged ST B->bZ", 1000, 0., 3000);
  h1_["mergeSTCon_bH"] = fs->make<TH1D>("mergeSTCon_bH", "Merged ST B->bH", 1000, 0., 3000);
  h1_["comboSTCon_bZ"] = fs->make<TH1D>("comboSTCon_bZ", "Combo ST B->bZ", 1000, 0., 3000);
  h1_["comboSTCon_bH"] = fs->make<TH1D>("comboSTCon_bH", "Combo ST B->bH", 1000, 0., 3000);

  h1_["STCon"] = fs->make<TH1D>("STCon", "Sum Pt", 1000, 0., 8000.);


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


  
 
