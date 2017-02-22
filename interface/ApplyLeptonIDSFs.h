#ifndef ANALYSIS_VLQANA_APPLYLEPTONIDSFS_HH
#define ANALYSIS_VLQANA_APPLYLEPTONIDSFS_HH

class ApplyLeptonIDSFs {
  public:
    enum LEPTONIDTYPES_t {LOOSE, TIGHT} ; 
    ApplyLeptonIDSFs (edm::ParameterSet const& pars) : 
      sf_(1),
      zdecayMode_(pars.getParameter<std::string>("zdecayMode")) 
  {
    std::string lepidtypestr = pars.getParameter<std::string>("lepidtype") ; 
    if ( lepidtypestr == "LOOSE" ) type_ = LOOSE ;  
    else if ( lepidtypestr == "TIGHT" ) type_ = TIGHT ; 
    else edm::LogError("ApplyLeptonSF") << " >>>> WrongElectronIdType: " << type_<< " Check lepton id type !!!" ; 
  }
    ~ApplyLeptonIDSFs () {} 
    double IDSF (double pt, double eta){
      if (type_ == TIGHT && zdecayMode_ == "zelel"){
        if(pt > 150.) pt = 150.; 
        if(pt > 10 && pt <= 20){
          if(eta <= -2.0          && eta > -2.5)    sf_ = 0.806567 ;
          else if(eta <= -1.566   && eta > -2.0)    sf_ = 0.828571 ;
          else if(eta <= -1.4442  && eta > -1.566)  sf_ = 1.03286;
          else if(eta <= -0.8     && eta > -1.4442) sf_ = 1.00775;
          else if(eta <= 0.0      && eta > -0.8)    sf_ = 0.940722;
          else if(eta <= 0.8      && eta > 0.0)     sf_ = 0.945876;
          else if(eta <= 1.4442   && eta > 0.8)     sf_ = 0.98977;
          else if(eta <= 1.566    && eta > 1.4442)  sf_ = 1.03398;
          else if(eta <= 2.0      && eta > 1.566)   sf_ = 0.827485;
          else if(eta < 2.5       && eta > 2.0)     sf_ = 0.797015;
        }
        else if(pt > 20 && pt <= 35){
          if(eta <= -2.0          && eta > -2.5)    sf_ = 0.882459  ;
          else if(eta <= -1.566   && eta > -2.0)    sf_ = 0.927487 ;
          else if(eta <= -1.4442  && eta > -1.566)  sf_ = 1.00759 ;
          else if(eta <= -0.8     && eta > -1.4442) sf_ = 0.972039 ;
          else if(eta <= 0.0      && eta > -0.8)    sf_ = 0.952998 ;
          else if(eta <= 0.8      && eta > 0.0)     sf_ = 0.981997 ;
          else if(eta <= 1.4442   && eta > 0.8)     sf_ = 0.975207 ;
          else if(eta <= 1.566    && eta > 1.4442)  sf_ = 0.974874;
          else if(eta <= 2.0      && eta > 1.566)   sf_ = 0.908938 ;
          else if(eta < 2.5       && eta > 2.0)     sf_ = 0.863322 ;
        }
        else if(pt > 35 && pt <= 50){
          if(eta <= -2.0          && eta > -2.5)    sf_ =  0.918841 ;
          else if(eta <= -1.566   && eta > -2.0)    sf_ =  0.966988;
          else if(eta <= -1.4442  && eta > -1.566)  sf_ =  0.988196 ;
          else if(eta <= -0.8     && eta > -1.4442) sf_ = 0.974667 ;
          else if(eta <= 0.0      && eta > -0.8)    sf_ = 0.953457 ;
          else if(eta <= 0.8      && eta > 0.0)     sf_ = 0.979812 ;
          else if(eta <= 1.4442   && eta > 0.8)     sf_ = 0.974633 ;
          else if(eta <= 1.566    && eta > 1.4442)  sf_ = 0.966387 ;
          else if(eta <= 2.0      && eta > 1.566)   sf_ = 0.957823 ;
          else if(eta < 2.5       && eta > 2.0)     sf_ = 0.907801 ;
        }
        else if(pt > 50 && pt <= 90){
          if(eta <= -2.0          && eta > -2.5)    sf_ = 0.940397 ;
          else if(eta <= -1.566   && eta > -2.0)    sf_ = 0.980964 ;
          else if(eta <= -1.4442  && eta > -1.566)  sf_ = 0.995392 ;
          else if(eta <= -0.8     && eta > -1.4442) sf_ = 0.971675 ;
          else if(eta <= 0.0      && eta > -0.8)    sf_ = 0.953086 ;
          else if(eta <= 0.8      && eta > 0.0)     sf_ = 0.977528 ;
          else if(eta <= 1.4442   && eta > 0.8)     sf_ = 0.978934 ;
          else if(eta <= 1.566    && eta > 1.4442)  sf_ = 0.979907 ;
          else if(eta <= 2.0      && eta > 1.566)   sf_ = 0.968553 ;
          else if(eta < 2.5       && eta > 2.0)     sf_ = 0.937622 ;
        }
        else if(pt > 90 && pt <=150){
          if(eta <= -2.0          && eta > -2.5)    sf_ = 1.05109 ;
          else if(eta <= -1.566   && eta > -2.0)    sf_ = 1.00592 ;
          else if(eta <= -1.4442  && eta > -1.566)  sf_ = 1.10411 ;
          else if(eta <= -0.8     && eta > -1.4442) sf_ = 0.989486 ;
          else if(eta <= 0.0      && eta > -0.8)    sf_ = 0.975496 ;
          else if(eta <= 0.8      && eta > 0.0)     sf_ = 1.01189 ;
          else if(eta <= 1.4442   && eta > 0.8)     sf_ = 1.01073 ;
          else if(eta <= 1.566    && eta > 1.4442)  sf_ = 1.00716 ;
          else if(eta <= 2.0      && eta > 1.566)   sf_ = 0.988263 ;
          else if(eta < 2.5       && eta > 2.0)     sf_ = 1.02136 ;
        }
      }//end TIGHT and zelel

      if (type_ == TIGHT && zdecayMode_ == "zmumu"){
        if(pt > 120.) pt = 120.;
        if (pt > 20 && pt <= 25){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 0.975305;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.986312;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 0.972437;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.986361;
        }
        else if (pt > 25 && pt <= 30){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 0.973942;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.986158;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 0.970745;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.983009;
        }
        else if (pt > 30 && pt <= 40){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 0.969471;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.987881;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 0.973607;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.984386;
        }
        else if (pt > 40 && pt <= 50){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 0.972574;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.990532;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 0.974868;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.9859;
        }
        else if (pt > 50 && pt <= 60){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 0.967713;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.985438;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 0.974566;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.981911;
        }
        else if (pt > 60 && pt <= 120){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 0.963148;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.988777;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 0.974918;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.992158;
        }
      }//TIGHT and zmumu
      return sf_ ; 
    }

    double IsoSF(double pt, double eta){
      if (type_ == TIGHT && zdecayMode_ == "zmumu"){
        if(pt > 120.) pt = 120.;
        if (pt > 20 && pt <= 25){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 0.983849;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.991505;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 0.995189;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.98476;
        }
        else if (pt > 25 && pt <= 30){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 0.991877;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.995906;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 1.00043;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.993309;
        }
        else if (pt > 30 && pt <= 40){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 0.996227;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.997965;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 0.999432;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.993728;
        }
        else if (pt > 40 && pt <= 50){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 0.998324;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.998013;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 0.997695;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.995297;
        }
        else if (pt > 50 && pt <= 60){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 0.998603;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.998345;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 0.999144;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.996805;
        }
        else if (pt > 60 && pt <= 120){
          if (abs(eta) <= 2.4 && abs(eta) > 2.1)       sf_ = 1.00152;
          else if (abs(eta) <= 2.1 && abs(eta) > 1.2)  sf_ = 0.999213;
          else if (abs(eta) <= 1.2 && abs(eta) > 0.9)  sf_ = 0.999122;
          else if (abs(eta) <= 0.9 && abs(eta) > 0.0)  sf_ = 0.99878;
        }
      }//TIGHT and zmumu                                                                                                                                                                                                       
      return sf_ ;
    }

  private: 
    double sf_ ;
    const std::string zdecayMode_ ; 
    LEPTONIDTYPES_t type_ ; 

};
#endif 
