#ifndef ANALYSIS_VLQANA_ZCANDSPRODUCER_HH
#define ANALYSIS_VLQANA_ZCANDSPRODUCER_HH

#include <iostream>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisDataFormats/BoostedObjects/interface/Candidate.h"

using namespace std;

class ZCandsProducer {
  public:


  void operator()(int N, int K, vlq::JetCollection& goodJets, vlq::CandidateCollection& Z)
  {
    std::string bitmask(K, 1); // K leading 1's                                                                                                               
    bitmask.resize(N, 0); // N-K trailing 0'                                                                                                                  
    double sum = 0;
    TLorentzVector jetinfo;
    jetinfo.SetPtEtaPhiM(0,0,0,0);
    vector<double> A, B;
    vector<TLorentzVector> C,D; 
    // print integers and permute bitmask                                                                                                                     
    do {
      for (int i = 0; i < N; ++i) // [0..N-1] integers                                                                                                        
	{
	  
	  if (bitmask[i]) {//cout << " " << i; 
	    sum +=  goodJets.at(i).getMass();
	    jetinfo += goodJets.at(i).getP4();
	    // cout<< " jet Mass is " << goodJets.at(i).getMass()<<endl;
	    
	     //   cout<< " the Pt sum is = " << sum<<endl;
	    A.push_back(sum);
	    C.push_back(jetinfo);
	  }
	  
	  //  if(A[2]<120. && A[2]>240.){continue;} 
	}
     
      // cout<< "Total mass of 2 jets  is  = " << A[1]<<endl;
       // cout<< "Total mass of 2 jets  ****T LORENTZ ****  = " << C.at(1).Mag()<<endl;
       B.push_back(A[1]);
       D.push_back(C.at(1));      
//if(A[2]<120. && A[2]>240.){cout<< "Total mass of 3 jets  is ******  = " << A[2]<<endl; B.push_back(A[2]);}
      sum =0;
      jetinfo.SetPtEtaPhiM(0,0,0,0);
      A.clear();
      C.clear();
      // cout << endl;
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    for ( unsigned j=0; j<D.size(); j++){
      //cout << " j the element is " << D.at(j).M()<<endl;

      if (D.at(j).Pt()>100. && D.at(j).Mag()>=70 && D.at(j).Mag()<=120){

	//	cout << " j the element is *********" << D.at(j).Mag()<<endl;
       	TLorentzVector zp4 = D.at(j);
	vlq::Candidate z(zp4) ;
	//	cout << " Z mass  is *********" << z.getMass()<<endl;
       	Z.push_back(z) ; 
      }
    }
    
  }
  
   
    ~ZCandsProducer () {}  

    // private:
  
}; 
#endif 
