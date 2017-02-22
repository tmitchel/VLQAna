#ifndef ANALYSIS_VLQANA_APPLYLEPTONTRIGSFS_HH
#define ANALYSIS_VLQANA_APPLYLEPTONTRIGSFS_HH

class ApplyLeptonTrigSFs {
  public:

    ApplyLeptonTrigSFs (edm::ParameterSet const& iConfig) : 
      sf_(1),
      zdecayMode_(iConfig.getParameter<std::string>("zdecayMode")),
      eltrigeffmap_(zdecayMode_ == "zelel" ? iConfig.getParameter<std::string>("eltrigeffmap") : "") 
    { }

    ~ApplyLeptonTrigSFs () {} 

    double operator () (double pt, double eta) {

      if (zdecayMode_ == "zelel") {

        if(pt < 120.) pt = 120.1;
        if(pt > 600.) pt = 599.9;
        std::unique_ptr<TFile>f_effmap = std::unique_ptr<TFile>(new TFile(eltrigeffmap_.c_str())) ;
        std::unique_ptr<TH2>h2_effmap = std::unique_ptr<TH2>(dynamic_cast<TH2*>(f_effmap->Get("eta_pt"))) ; 
        int binpt = h2_effmap->GetXaxis()->FindBin(pt);
        int bineta = h2_effmap->GetYaxis()->FindBin(eta);
        sf_ = h2_effmap->GetBinContent(binpt, bineta) ; 

      } //// zelel
      else if (zdecayMode_ == "zmumu") {

        eta = abs(eta);
        if(pt > 500.) pt = 500.0;
        if(pt > 26 && pt <= 30){
          if (eta <= 2.4 && eta > 2.1)       sf_ = 0.906216;
          else if (eta <= 2.1 && eta > 1.2)  sf_ = 0.982503;
          else if (eta <= 1.2 && eta > 0.9)  sf_ = 0.956337;
          else if (eta <= 0.9 && eta > 0.0)  sf_ = 0.980357;
        }

        else if(pt > 30 && pt <= 40){
          if (eta <= 2.4 && eta > 2.1)       sf_ = 0.944555;
          else if (eta <= 2.1 && eta > 1.2)  sf_ = 0.99548;
          else if (eta <= 1.2 && eta > 0.9)  sf_ = 0.965959;
          else if (eta <= 0.9 && eta > 0.0)  sf_ = 0.984096;
        }

        else if(pt > 40 && pt <= 50){
          if (eta <= 2.4 && eta > 2.1)       sf_ = 0.957532;
          else if (eta <= 2.1 && eta > 1.2)  sf_ = 0.99912;
          else if (eta <= 1.2 && eta > 0.9)  sf_ = 0.967858;
          else if (eta <= 0.9 && eta > 0.0)  sf_ = 0.984904;
        }

        else if(pt > 50 && pt <= 60){
          if (eta <= 2.4 && eta > 2.1)       sf_ = 0.960397;
          else if (eta <= 2.1 && eta > 1.2)  sf_ = 0.99912;
          else if (eta <= 1.2 && eta > 0.9)  sf_ = 0.968175;
          else if (eta <= 0.9 && eta > 0.0)  sf_ = 0.985122;
        }

        else if(pt > 60 && pt <= 120){
          if (eta <= 2.4 && eta > 2.1)       sf_ = 0.95241;
          else if (eta <= 2.1 && eta > 1.2)  sf_ = 0.998937;
          else if (eta <= 1.2 && eta > 0.9)  sf_ = 0.964898;
          else if (eta <= 0.9 && eta > 0.0)  sf_ = 0.983897;
        }

        else if(pt > 120 && pt <= 200){
          if (eta <= 2.4 && eta > 2.1)       sf_ = 0.977009;
          else if (eta <= 2.1 && eta > 1.2)  sf_ = 1.00552;
          else if (eta <= 1.2 && eta > 0.9)  sf_ = 0.946724;
          else if (eta <= 0.9 && eta > 0.0)  sf_ = 0.976031;
        }

        else if(pt > 200 && pt <= 500){
          if (eta <= 2.4 && eta > 2.1)       sf_ = 0.918103;
          else if (eta <= 2.1 && eta > 1.2)  sf_ = 0.982804;
          else if (eta <= 1.2 && eta > 0.9)  sf_ = 0.949605;
          else if (eta <= 0.9 && eta > 0.0)  sf_ = 0.983284;
        }

      } /// zmumu

      return sf_ ;
    } 

  private: 
    double sf_ ;
    const std::string zdecayMode_ ; 
    const std::string eltrigeffmap_ ; 
};
#endif 
