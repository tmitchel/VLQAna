#ifndef ANALYSIS_VLQANA_TTHT_HH
#define ANALYSIS_VLQANA_TTHT_HH
#include <string>
#include <TTree.h>
#include <vector>
#include <utility>  

class TtHEventInfoBranches {
  public:

    int    runno_ ;
    int    lumisec_ ; 
    int    evtno_;    
    int    bcno_;
    int    time_;
    double htHat_; 
    double EvtWeight_ ; 
    double EvtWtPV_;
    double EvtWtPVLow_;
    double EvtWtPVHigh_;
    double EvtWtHT_;
    double EvtWtHTUp_;
    double EvtWtHTDown_;
    int    npv_;
    int    npuTrue_;
    double toptagsf_;
    double toptagsf_Up_;
    double toptagsf_Down_;
    double btagsf_;
    double btagsf_bcUp_;
    double btagsf_bcDown_;
    double btagsf_lUp_;
    double btagsf_lDown_;
    double mtprime_;
    double mtprimeDummy_;
    double mtprimeDummyA_;
    double mtprimeDummyC_;
    double ht_ ; 
    int    nAK4_;
    int    nAK8_;
    bool   isRegionA_;
    bool   isRegionB_;
    bool   isRegionC_;
    bool   isRegionD_;
    bool   isRegionNotABCD_;
    std::vector<std::pair<int, double> > lhewts_ ; 

    void RegisterTree(TTree* cutTree, std::string name="SelectedEvents") {
      cutTree->Branch((name+"_runno").c_str(), &runno_, (name+"_runno/I").c_str());
      cutTree->Branch((name+"_lumisec").c_str(), &lumisec_, (name+"_lumisec/I").c_str());
      cutTree->Branch((name+"_evtno").c_str(), &evtno_, (name+"_evtno/I").c_str());
      cutTree->Branch((name+"_bcno").c_str(), &bcno_, (name+"_bcno/I").c_str());
      cutTree->Branch((name+"_time").c_str(), &time_, (name+"_time/I").c_str());
      cutTree->Branch((name+"_htHat").c_str(), &htHat_, (name+"_htHat/D").c_str()) ; 
      cutTree->Branch((name+"_EvtWeight").c_str(), &EvtWeight_, "EvtWeight/D");
      cutTree->Branch((name+"_EvtWtPV").c_str(), &EvtWtPV_, "EvtWtPV/D");
      cutTree->Branch((name+"_EvtWtPVLow").c_str(), &EvtWtPVLow_, "EvtWtPVLow/D");
      cutTree->Branch((name+"_EvtWtPVHigh").c_str(), &EvtWtPVHigh_, "EvtWtPVHigh/D");
      cutTree->Branch((name+"_EvtWtHT").c_str(), &EvtWtHT_, "EvtWtHT/D");
      cutTree->Branch((name+"_EvtWtHTUp").c_str(), &EvtWtHTUp_, "EvtWtHTUp/D");
      cutTree->Branch((name+"_EvtWtHTDown").c_str(), &EvtWtHTDown_, "EvtWtHTDown/D");
      cutTree->Branch((name+"_npv").c_str(), &npv_, "npv/I");
      cutTree->Branch((name+"_npuTrue").c_str(), &npuTrue_, "npuTrue/I");
      cutTree->Branch((name+"_toptagsf").c_str(), &toptagsf_, "toptagsf/D");
      cutTree->Branch((name+"_toptagsf_Up").c_str(), &toptagsf_Up_, "toptagsf_Up/D");
      cutTree->Branch((name+"_toptagsf_Down").c_str(), &toptagsf_Down_, "toptagsf_Down/D");
      cutTree->Branch((name+"_btagsf").c_str(), &btagsf_, "btagsf/D");
      cutTree->Branch((name+"_btagsf_bcUp").c_str(), &btagsf_bcUp_, "btagsf_bcUp/D");
      cutTree->Branch((name+"_btagsf_bcDown").c_str(), &btagsf_bcDown_, "btagsf_bcDown/D");
      cutTree->Branch((name+"_btagsf_lUp").c_str(), &btagsf_lUp_, "btagsf_lUp/D");
      cutTree->Branch((name+"_btagsf_lDown").c_str(), &btagsf_lDown_, "btagsf_lDown/D");
      cutTree->Branch((name+"_mtprime").c_str(), &mtprime_, "mtprime/D");
      cutTree->Branch((name+"_mtprimeDummy").c_str(), &mtprimeDummy_, "mtprimeDummy/D");
      cutTree->Branch((name+"_mtprimeDummyA").c_str(), &mtprimeDummyA_, "mtprimeDummyA/D");
      cutTree->Branch((name+"_mtprimeDummyC").c_str(), &mtprimeDummyC_, "mtprimeDummyC/D");
      cutTree->Branch((name+"_ht").c_str(), &ht_, "ht/D");
      cutTree->Branch((name+"_nAK4").c_str(), &nAK4_, "nAK4/I");
      cutTree->Branch((name+"_nAK8").c_str(), &nAK8_, "nAK8/I");
      cutTree->Branch((name+"_isRegionA").c_str(), &isRegionA_, "isRegionA/O");
      cutTree->Branch((name+"_isRegionB").c_str(), &isRegionB_, "isRegionB/O");
      cutTree->Branch((name+"_isRegionC").c_str(), &isRegionC_, "isRegionC/O");
      cutTree->Branch((name+"_isRegionD").c_str(), &isRegionD_, "isRegionD/O");
      cutTree->Branch((name+"_isRegionNotABCD").c_str(), &isRegionNotABCD_, "isRegionNotABCD/O");
      cutTree->Branch((name+"_lhewts").c_str(), &lhewts_) ;  
    }

};

class TtHJetInfoBranches {
  public:

    std::vector<int>  idxAK4; 
    std::vector<double> ptAK4;
    std::vector<double> etaAK4;
    std::vector<double> phiAK4;
    std::vector<double> MAK4;
    std::vector<double> csvAK4;
    std::vector<double> partonFlavourAK4; 
    std::vector<double> hadronFlavourAK4; 

    std::vector<int> idxAK8; 
    std::vector<double> ptAK8;
    std::vector<double> etaAK8;
    std::vector<double> phiAK8;
    std::vector<double> MAK8;
    std::vector<double> SoftDropMassAK8;
    std::vector<double> PrunedMassAK8;
    std::vector<double> tau1AK8;
    std::vector<double> tau2AK8;
    std::vector<double> tau3AK8;
    std::vector<double> csvAK8;
    std::vector<double> partonFlavourAK8; 
    std::vector<double> hadronFlavourAK8; 

    std::vector<double> doubleBAK8;
    std::vector<double> sj0CSVAK8;
    std::vector<double> sj1CSVAK8;
    std::vector<double> hadronFlavourSJ0AK8; 
    std::vector<double> hadronFlavourSJ1AK8; 

    std::vector<int> idxHTagged; 
    std::vector<double> ptHTagged;
    std::vector<double> etaHTagged;
    std::vector<double> phiHTagged;
    std::vector<double> MHTagged;
    std::vector<double> SoftDropMassHTagged;
    std::vector<double> PrunedMassHTagged;
    std::vector<double> tau1HTagged;
    std::vector<double> tau2HTagged;
    std::vector<double> tau3HTagged;
    std::vector<double> csvHTagged;
    std::vector<double> partonFlavourHTagged; 
    std::vector<double> hadronFlavourHTagged; 

    std::vector<double> doubleBHTagged;
    std::vector<double> sj0CSVHTagged;
    std::vector<double> sj1CSVHTagged;
    std::vector<double> hadronFlavourSJ0HTagged; 
    std::vector<double> hadronFlavourSJ1HTagged; 

    std::vector<int> idxAntiHTagged; 
    std::vector<double> ptAntiHTagged;
    std::vector<double> etaAntiHTagged;
    std::vector<double> phiAntiHTagged;
    std::vector<double> MAntiHTagged;
    std::vector<double> SoftDropMassAntiHTagged;
    std::vector<double> PrunedMassAntiHTagged;
    std::vector<double> tau1AntiHTagged;
    std::vector<double> tau2AntiHTagged;
    std::vector<double> tau3AntiHTagged;
    std::vector<double> csvAntiHTagged;
    std::vector<double> partonFlavourAntiHTagged; 
    std::vector<double> hadronFlavourAntiHTagged; 

    std::vector<double> doubleBAntiHTagged;
    std::vector<double> sj0CSVAntiHTagged;
    std::vector<double> sj1CSVAntiHTagged;
    std::vector<double> hadronFlavourSJ0AntiHTagged; 
    std::vector<double> hadronFlavourSJ1AntiHTagged; 

    std::vector<int> idxTopTagged; 
    std::vector<double> ptTopTagged;
    std::vector<double> etaTopTagged;
    std::vector<double> phiTopTagged;
    std::vector<double> MTopTagged;
    std::vector<double> SoftDropMassTopTagged;
    std::vector<double> PrunedMassTopTagged;
    std::vector<double> tau1TopTagged;
    std::vector<double> tau2TopTagged;
    std::vector<double> tau3TopTagged;
    std::vector<double> csvTopTagged;
    std::vector<double> partonFlavourTopTagged; 
    std::vector<double> hadronFlavourTopTagged; 

    std::vector<double> doubleBTopTagged;
    std::vector<double> sj0CSVTopTagged;
    std::vector<double> sj1CSVTopTagged;
    std::vector<double> hadronFlavourSJ0TopTagged; 
    std::vector<double> hadronFlavourSJ1TopTagged; 

    std::vector<int> idxAntiTopTagged; 
    std::vector<double> ptAntiTopTagged;
    std::vector<double> etaAntiTopTagged;
    std::vector<double> phiAntiTopTagged;
    std::vector<double> MAntiTopTagged;
    std::vector<double> SoftDropMassAntiTopTagged;
    std::vector<double> PrunedMassAntiTopTagged;
    std::vector<double> tau1AntiTopTagged;
    std::vector<double> tau2AntiTopTagged;
    std::vector<double> tau3AntiTopTagged;
    std::vector<double> csvAntiTopTagged;
    std::vector<double> partonFlavourAntiTopTagged; 
    std::vector<double> hadronFlavourAntiTopTagged; 

    std::vector<double> doubleBAntiTopTagged;
    std::vector<double> sj0CSVAntiTopTagged;
    std::vector<double> sj1CSVAntiTopTagged;
    std::vector<double> hadronFlavourSJ0AntiTopTagged; 
    std::vector<double> hadronFlavourSJ1AntiTopTagged;
 
    void RegisterTree(TTree *cutTree, std::string name="JetInfo") {

      cutTree->Branch("idxAK4", &idxAK4); 
      cutTree->Branch("ptAK4", &ptAK4); 
      cutTree->Branch("etaAK4", &etaAK4);
      cutTree->Branch("phiAK4", &phiAK4);
      cutTree->Branch("MAK4", &MAK4);
      cutTree->Branch("csvAK4", &csvAK4);
      cutTree->Branch("partonFlavourAK4", &partonFlavourAK4);
      cutTree->Branch("hadronFlavourAK4", &hadronFlavourAK4);

      cutTree->Branch("idxAK8", &idxAK8); 
      cutTree->Branch("ptAK8", &ptAK8); 
      cutTree->Branch("etaAK8", &etaAK8);
      cutTree->Branch("phiAK8", &phiAK8);
      cutTree->Branch("MAK8", &MAK8);
      cutTree->Branch("SoftDropMassAK8", &SoftDropMassAK8);
      cutTree->Branch("PrunedMassAK8", &PrunedMassAK8);
      cutTree->Branch("tau1AK8", &tau1AK8);
      cutTree->Branch("tau2AK8", &tau2AK8);
      cutTree->Branch("tau3AK8", &tau3AK8);
      cutTree->Branch("csvAK8", &csvAK8);
      cutTree->Branch("partonFlavourAK8", &partonFlavourAK8);
      cutTree->Branch("hadronFlavourAK8", &hadronFlavourAK8);

      cutTree->Branch("doubleBAK8",&doubleBAK8);
      cutTree->Branch("sj0CSVAK8",&sj0CSVAK8);
      cutTree->Branch("sj1CSVAK8",&sj1CSVAK8);
      cutTree->Branch("hadronFlavourSJ0AK8",&hadronFlavourSJ0AK8);
      cutTree->Branch("hadronFlavourSJ1AK8",&hadronFlavourSJ1AK8);

      cutTree->Branch("idxHTagged", &idxHTagged); 
      cutTree->Branch("ptHTagged", &ptHTagged); 
      cutTree->Branch("etaHTagged", &etaHTagged);
      cutTree->Branch("phiHTagged", &phiHTagged);
      cutTree->Branch("MHTagged", &MHTagged);
      cutTree->Branch("SoftDropMassHTagged", &SoftDropMassHTagged);
      cutTree->Branch("PrunedMassHTagged", &PrunedMassHTagged);
      cutTree->Branch("tau1HTagged", &tau1HTagged);
      cutTree->Branch("tau2HTagged", &tau2HTagged);
      cutTree->Branch("tau3HTagged", &tau3HTagged);
      cutTree->Branch("csvHTagged", &csvHTagged);
      cutTree->Branch("partonFlavourHTagged", &partonFlavourHTagged);
      cutTree->Branch("hadronFlavourHTagged", &hadronFlavourHTagged);

      cutTree->Branch("doubleBHTagged",&doubleBHTagged);
      cutTree->Branch("sj0CSVHTagged",&sj0CSVHTagged);
      cutTree->Branch("sj1CSVHTagged",&sj1CSVHTagged);
      cutTree->Branch("hadronFlavourSJ0HTagged",&hadronFlavourSJ0HTagged);
      cutTree->Branch("hadronFlavourSJ1HTagged",&hadronFlavourSJ1HTagged);

      cutTree->Branch("idxAntiHTagged", &idxAntiHTagged); 
      cutTree->Branch("ptAntiHTagged", &ptAntiHTagged); 
      cutTree->Branch("etaAntiHTagged", &etaAntiHTagged);
      cutTree->Branch("phiAntiHTagged", &phiAntiHTagged);
      cutTree->Branch("MAntiHTagged", &MAntiHTagged);
      cutTree->Branch("SoftDropMassAntiHTagged", &SoftDropMassAntiHTagged);
      cutTree->Branch("PrunedMassAntiHTagged", &PrunedMassAntiHTagged);
      cutTree->Branch("tau1AntiHTagged", &tau1AntiHTagged);
      cutTree->Branch("tau2AntiHTagged", &tau2AntiHTagged);
      cutTree->Branch("tau3AntiHTagged", &tau3AntiHTagged);
      cutTree->Branch("csvAntiHTagged", &csvAntiHTagged);
      cutTree->Branch("partonFlavourAntiHTagged", &partonFlavourAntiHTagged);
      cutTree->Branch("hadronFlavourAntiHTagged", &hadronFlavourAntiHTagged);

      cutTree->Branch("doubleBAntiHTagged",&doubleBAntiHTagged);
      cutTree->Branch("sj0CSVAntiHTagged",&sj0CSVAntiHTagged);
      cutTree->Branch("sj1CSVAntiHTagged",&sj1CSVAntiHTagged);
      cutTree->Branch("hadronFlavourSJ0AntiHTagged",&hadronFlavourSJ0AntiHTagged);
      cutTree->Branch("hadronFlavourSJ1AntiHTagged",&hadronFlavourSJ1AntiHTagged);

      cutTree->Branch("idxTopTagged", &idxTopTagged); 
      cutTree->Branch("ptTopTagged", &ptTopTagged); 
      cutTree->Branch("etaTopTagged", &etaTopTagged);
      cutTree->Branch("phiTopTagged", &phiTopTagged);
      cutTree->Branch("MTopTagged", &MTopTagged);
      cutTree->Branch("SoftDropMassTopTagged", &SoftDropMassTopTagged);
      cutTree->Branch("PrunedMassTopTagged", &PrunedMassTopTagged);
      cutTree->Branch("tau1TopTagged", &tau1TopTagged);
      cutTree->Branch("tau2TopTagged", &tau2TopTagged);
      cutTree->Branch("tau3TopTagged", &tau3TopTagged);
      cutTree->Branch("csvTopTagged", &csvTopTagged);
      cutTree->Branch("partonFlavourTopTagged", &partonFlavourTopTagged);
      cutTree->Branch("hadronFlavourTopTagged", &hadronFlavourTopTagged);

      cutTree->Branch("doubleBTopTagged",&doubleBTopTagged);
      cutTree->Branch("sj0CSVTopTagged",&sj0CSVTopTagged);
      cutTree->Branch("sj1CSVTopTagged",&sj1CSVTopTagged);
      cutTree->Branch("hadronFlavourSJ0TopTagged",&hadronFlavourSJ0TopTagged);
      cutTree->Branch("hadronFlavourSJ1TopTagged",&hadronFlavourSJ1TopTagged);

      cutTree->Branch("idxAntiTopTagged", &idxAntiTopTagged); 
      cutTree->Branch("ptAntiTopTagged", &ptAntiTopTagged); 
      cutTree->Branch("etaAntiTopTagged", &etaAntiTopTagged);
      cutTree->Branch("phiAntiTopTagged", &phiAntiTopTagged);
      cutTree->Branch("MAntiTopTagged", &MAntiTopTagged);
      cutTree->Branch("SoftDropMassAntiTopTagged", &SoftDropMassAntiTopTagged);
      cutTree->Branch("PrunedMassAntiTopTagged", &PrunedMassAntiTopTagged);
      cutTree->Branch("tau1AntiTopTagged", &tau1AntiTopTagged);
      cutTree->Branch("tau2AntiTopTagged", &tau2AntiTopTagged);
      cutTree->Branch("tau3AntiTopTagged", &tau3AntiTopTagged);
      cutTree->Branch("csvAntiTopTagged", &csvAntiTopTagged);
      cutTree->Branch("partonFlavourAntiTopTagged", &partonFlavourAntiTopTagged);
      cutTree->Branch("hadronFlavourAntiTopTagged", &hadronFlavourAntiTopTagged);

      cutTree->Branch("doubleBAntiTopTagged",&doubleBAntiTopTagged);
      cutTree->Branch("sj0CSVAntiTopTagged",&sj0CSVAntiTopTagged);
      cutTree->Branch("sj1CSVAntiTopTagged",&sj1CSVAntiTopTagged);
      cutTree->Branch("hadronFlavourSJ0AntiTopTagged",&hadronFlavourSJ0AntiTopTagged);
      cutTree->Branch("hadronFlavourSJ1AntiTopTagged",&hadronFlavourSJ1AntiTopTagged);

    }

    //void ReadTree(TTree* tree, std::string name="JetInfo") {
    //}

};

#endif 
