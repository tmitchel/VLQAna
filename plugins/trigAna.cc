#include <iostream>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Analysis/VLQAna/interface/ElectronMaker.h"
#include "Analysis/VLQAna/interface/DileptonCandsProducer.h"
#include "Analysis/VLQAna/interface/CandidateCleaner.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>

class trigAna : public edm::EDFilter {
public:
    explicit trigAna(const edm::ParameterSet&);
    ~trigAna();

private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    edm::EDGetTokenT<bool> t_phltdec;
    edm::EDGetTokenT<int>  t_pevtno;   
    edm::EDGetTokenT<int>  t_prunno;   
    edm::EDGetTokenT<int>  t_plumisec;   

    edm::EDGetTokenT<bool> t_thltdec;
    edm::EDGetTokenT<int>  t_tevtno;   
    edm::EDGetTokenT<int>  t_trunno;   
    edm::EDGetTokenT<int>  t_tlumisec;   

    const bool isData;
    edm::ParameterSet DilepCandParams; 
    ElectronMaker electonmaker;

    edm::Service<TFileService> fs   ; 
    std::map<std::string, TH1D*> h1_; 
    std::map<std::string, TH2D*> h2_;

};

trigAna::trigAna(const edm::ParameterSet& iConfig) :
    t_phltdec   (consumes<bool> (iConfig.getParameter<edm::InputTag>("probe_hltdec"))),
    t_pevtno    (consumes<int>  (iConfig.getParameter<edm::InputTag>("probe_evtno"))),
    t_prunno    (consumes<int>  (iConfig.getParameter<edm::InputTag>("probe_runno"))),
    t_plumisec  (consumes<int>  (iConfig.getParameter<edm::InputTag>("probe_lumisec"))),

    t_thltdec   (consumes<bool> (iConfig.getParameter<edm::InputTag>("tag_hltdec"))),
    t_tevtno    (consumes<int>  (iConfig.getParameter<edm::InputTag>("tag_evtno"))),
    t_trunno    (consumes<int>  (iConfig.getParameter<edm::InputTag>("tag_runno"))),
    t_tlumisec  (consumes<int>  (iConfig.getParameter<edm::InputTag>("tag_lumisec"))),

    isData          (iConfig.getParameter<bool> ("isData")),
    DilepCandParams (iConfig.getParameter<edm::ParameterSet> ("DilepCandParams")),
    electonmaker    (iConfig.getParameter<edm::ParameterSet> ("elselParams"),consumesCollector())
{}

trigAna::~trigAna() {}

bool trigAna::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
    using namespace edm;

    Handle<bool> h_probe_hltdec   ; evt.getByToken(t_phltdec  , h_probe_hltdec  );
    Handle<int>  h_probe_evtno    ; evt.getByToken(t_pevtno   , h_probe_evtno   );
    Handle<int>  h_probe_runno    ; evt.getByToken(t_prunno   , h_probe_runno   );
    Handle<int>  h_probe_lumisec  ; evt.getByToken(t_plumisec , h_probe_lumisec );

    Handle<bool> h_tag_hltdec     ; evt.getByToken(t_thltdec  , h_tag_hltdec    );
    Handle<int>  h_tag_evtno      ; evt.getByToken(t_tevtno   , h_tag_evtno     );
    Handle<int>  h_tag_runno      ; evt.getByToken(t_trunno   , h_tag_runno     );
    Handle<int>  h_tag_lumisec    ; evt.getByToken(t_tlumisec , h_tag_lumisec   );

    double evtwt(1.);

    h1_["cutflow"] -> Fill(1, evtwt);

    // only interested in events passing the tag
    if (*h_tag_hltdec.product()) h1_["cutflow"] -> Fill(2, evtwt);
    else return false;

    vlq::ElectronCollection electrons;
    electonmaker(evt, electrons);

    vlq::CandidateCollection dilepton;
    DileptonCandsProducer dileptonsprod(DilepCandParams);
    dileptonsprod.operator()<vlq::ElectronCollection>(dilepton, electrons);

    // want a single dilepton pair
    if (dilepton.size() == 1) h1_["cutflow"] -> Fill(3, evtwt);
    else return false;

    // only look at events where a single electron passes the probe trig. 
    vector<double> highPt;
    for (auto& el : electrons) {
        if (el.getPt() > 115)
            highPt.push_back(el.getPt());
    }

    if (highPt.size() == 1) h1_["cutflow"] -> Fill(4, evtwt);
    else return false;

    h1_["pt_all"] -> Fill(electrons.at(0).getPt(), evtwt);
    h1_["eta_all"] -> Fill(electrons.at(0).getEta(), evtwt);
    h2_["all"] -> Fill(electrons.at(0).getPt(), electrons.at(0).getEta());
    h2_["fixed_all"] -> Fill(electrons.at(0).getPt(), electrons.at(0).getEta());


    if (*h_probe_hltdec.product()) {
        h1_["pt_pass_probe"] -> Fill(electrons.at(0).getPt(), evtwt);
        h1_["eta_pass_probe"] -> Fill(electrons.at(0).getEta(), evtwt);
        h2_["pass"] -> Fill(electrons.at(0).getPt(), electrons.at(0).getEta());
        h2_["fixed_pass"] -> Fill(electrons.at(0).getPt(), electrons.at(0).getEta());
    }
    else {
        h1_["pt_fail_probe"] -> Fill(electrons.at(0).getPt(), evtwt);
        h1_["eta_fail_probe"] -> Fill(electrons.at(0).getEta(), evtwt);
        h2_["fail"] -> Fill(electrons.at(0).getPt(), electrons.at(0).getEta());
        h2_["fixed_fail"] -> Fill(electrons.at(0).getPt(), electrons.at(0).getEta());
    }

    return true;
}

void trigAna::beginJob() {

    double ptbins[] = {115, 125, 140, 200, 600}; 
    double etabins[] = {-2.4, -2.1, -1.9, 0, 1.9, 2.1, 2.4};

    h1_["cutflow"] = fs->make<TH1D>("cutflow", "cutflow", 5, 0.5, 5.5);
    h1_["pt_all"] = fs->make<TH1D>("pt_all", "All Tag Events", 48, 120., 600);
    h1_["pt_pass_probe"] = fs->make<TH1D>("pt_pass_probe", "Probe Events", 48, 120., 600);
    h1_["pt_fail_probe"] = fs->make<TH1D>("pt_fail_probe", "Anti-Probe Events", 48, 120., 600);
    h1_["eta_all"] = fs->make<TH1D>("eta_all", "All Tag Events", 100, -4., 4.);
    h1_["eta_pass_probe"] = fs->make<TH1D>("eta_pass_probe", "Probe Events", 100, -4., 4);
    h1_["eta_fail_probe"] = fs->make<TH1D>("eta_fail_probe", "Anti-Probe Events", 100, -4., 4);

    h2_["fixed_all"] = fs->make<TH2D>("fixed_all", "All Events", 48, 120., 600., 100, -4., 4.);
    h2_["fixed_pass"] = fs->make<TH2D>("fixed_pass", "Probe Events", 48, 120., 600., 100, -4., 4.);
    h2_["fixed_fail"] = fs->make<TH2D>("fixed_fail", "Anti-Prove Events", 48, 120., 600., 100, -4., 4.);

    h2_["all"] = fs->make<TH2D>("2d_all", "All Events",  4, ptbins, 6, etabins);
    h2_["pass"] = fs->make<TH2D>("2d_pass", "Probe Events",  4, ptbins, 6, etabins);
    h2_["fail"] = fs->make<TH2D>("2d_fail", "Anti-Prove Events", 4, ptbins, 6, etabins);
}

void trigAna::endJob() {
    return;
}

DEFINE_FWK_MODULE(trigAna);

