#ifndef ANALYSIS_INTERFACE_APPLYJER_H
#define ANALYSIS_INTERFACE_APPLYJER_H

double ApplyJERMass(int jerShift) {
  double massjerscale(1.07) ; 
  if (jerShift == 2) massjerscale += 0.103 ;
  if (jerShift == -1) massjerscale -= 0.103 ;
  return massjerscale ; 
}

double ApplyJERp4 (double eta, int jerShift) {

  eta = abs(eta) ;
  if (eta >= 4.7) eta = 4.699 ;
  double jerscale(0) ; 

  if (eta >= 3.2 && eta < 5.0) {
    if (jerShift == 1) jerscale = 1.16  ;  
    else if (jerShift == 2) jerscale = 1.16 + 0.029 ;  
    else if (jerShift == -1) jerscale = 1.16 - 0.029 ;  
  }
  else if (eta >= 3.0 && eta < 3.2) {
    if (jerShift == 1) jerscale = 1.328  ;  
    else if (jerShift == 2) jerscale = 1.328 + 0.022;  
    else if (jerShift == -1) jerscale = 1.328 - 0.022;  
  }
  else if (eta >= 2.8 && eta < 3.0) {
    if (jerShift == 1) jerscale =  1.857 ;  
    else if (jerShift == 2) jerscale =  1.857 + 0.071;  
    else if (jerShift == -1) jerscale =  1.857 - 0.071;  
  }
  else if (eta >= 2.5 && eta < 2.8) {
    if (jerShift == 1) jerscale = 1.364  ;  
    else if (jerShift == 2) jerscale = 1.364 + 0.039;  
    else if (jerShift == -1) jerscale = 1.364 - 0.039;  
  }
  else if (eta >= 2.3 && eta < 2.5) {
    if (jerShift == 1) jerscale = 1.177 ;  
    else if (jerShift == 2) jerscale = 1.177 + 0.041;  
    else if (jerShift == -1) jerscale = 1.177 - 0.041;  
  }
  else if (eta >= 2.1 && eta < 2.3) {
    if (jerShift == 1) jerscale = 1.067 ;  
    else if (jerShift == 2) jerscale = 1.067 + 0.053;  
    else if (jerShift == -1) jerscale = 1.067 - 0.053 ;  
  }
  else if (eta >= 1.9 && eta < 2.1) {
    if (jerShift == 1) jerscale = 1.140  ;  
    else if (jerShift == 2) jerscale = 1.140 + 0.047;  
    else if (jerShift == -1) jerscale = 1.140 - 0.047;  
  }
  else if (eta >= 1.7 && eta < 1.9) {
    if (jerShift == 1) jerscale = 1.082  ;  
    else if (jerShift == 2) jerscale = 1.082 + 0.035;  
    else if (jerShift == -1) jerscale = 1.082 - 0.035 ;  
  }
  else if (eta >= 1.3 && eta < 1.7) {
    if (jerShift == 1) jerscale = 1.084  ;  
    else if (jerShift == 2) jerscale = 1.084 + 0.011;  
    else if (jerShift == -1) jerscale = 1.084 - 0.011 ;  
  }
  else if (eta >= 1.1 && eta < 1.3) {
    if (jerShift == 1) jerscale = 1.123  ;  
    else if (jerShift == 2) jerscale = 1.123 + 0.024;  
    else if (jerShift == -1) jerscale = 1.123 - 0.024;  
  }
  else if (eta >= 0.8 && eta < 1.1) {
    if (jerShift == 1) jerscale = 1.114  ;  
    else if (jerShift == 2) jerscale = 1.114 + 0.013;  
    else if (jerShift == -1) jerscale = 1.114 - 0.013;  
  }
  else if (eta >= 0.5 && eta < 0.8) {
    if (jerShift == 1) jerscale = 1.138  ;  
    else if (jerShift == 2) jerscale = 1.138 + 0.013;  
    else if (jerShift == -1) jerscale = 1.138 - 0.013;  
  }
  else if (eta >= 0.0 && eta < 0.5) {
    if (jerShift == 1) jerscale = 1.109  ;  
    else if (jerShift == 2) jerscale = 1.109 + 0.008;  
    else if (jerShift == -1) jerscale = 1.109 - 0.008;  
  }

  return jerscale ; 
  
//2016 JERs
//https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Spring16_25nsV10_MC/Spring16_25nsV10_MC_SF_AK4PFchs.txt
//Same SF for AF4PFchs and AK8PFchs 
  //{1 JetEta 0 None ScaleFactor}
  //-4.7  -3.2  3   1.160 1.131 1.189
  //-3.2  -3.0  3   1.328 1.306 1.350
  //-3.0  -2.8  3   1.857 1.786 1.928
  //-2.8  -2.5  3   1.364 1.325 1.403
  //-2.5  -2.3  3   1.177 1.136 1.218
  //-2.3  -2.1  3   1.067 1.014 1.120
  //-2.1  -1.9  3   1.140 1.093 1.187
  //-1.9  -1.7  3   1.082 1.047 1.117
  //-1.7  -1.3  3   1.084 1.073 1.095
  //-1.3  -1.1  3   1.123 1.099 1.147
  //-1.1  -0.8  3   1.114 1.101 1.127
  //-0.8  -0.5  3   1.138 1.125 1.151
  //-0.5   0.0  3   1.109 1.101 1.117
  //0.0   0.5  3   1.109 1.101 1.117
  //0.5   0.8  3   1.138 1.125 1.151
  //0.8   1.1  3   1.114 1.101 1.127
  //1.1   1.3  3   1.123 1.099 1.147
  //1.3   1.7  3   1.084 1.073 1.095
  //1.7   1.9  3   1.082 1.047 1.117
  //1.9   2.1  3   1.140 1.093 1.187
  //2.1   2.3  3   1.067 1.014 1.120
  //2.3   2.5  3   1.177 1.136 1.218
  //2.5   2.8  3   1.364 1.325 1.403
  //2.8   3.0  3   1.857 1.786 1.928
  //3.0   3.2  3   1.328 1.306 1.350
  //3.2   4.7  3   1.160 1.131 1.189
}

#endif 
