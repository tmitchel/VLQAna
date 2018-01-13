#ifndef ANALYSIS_VLQANA_TOPCANDSPRODUCER_HH
#define ANALYSIS_VLQANA_TOPCANDSPRODUCER_HH

#include <iostream>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisDataFormats/BoostedObjects/interface/Candidate.h"

using namespace std;

class TopCandsProducer {
  public:


  void operator()(int N, int K, vlq::JetCollection& goodJets, vlq::CandidateCollection& tops)
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
	    //    cout<< " jet Mass is " << goodJets.at(i).getMass()<<endl;
	    
	    //  cout<< " the Pt sum is = " << sum<<endl;
	    A.push_back(sum);
	    C.push_back(jetinfo);
	  }
	  
	  //  if(A[2]<120. && A[2]>240.){continue;} 
	}
     
      // cout<< "Total mass of 3 jets  is  = " << A[2]<<endl;
      // cout<< "Total mass of 3 jets  ****T LORENTZ ****  = " << C.at(2).Mag()<<endl;
      B.push_back(A[2]);
      D.push_back(C.at(2));      
//if(A[2]<120. && A[2]>240.){cout<< "Total mass of 3 jets  is ******  = " << A[2]<<endl; B.push_back(A[2]);}
      sum =0;
      jetinfo.SetPtEtaPhiM(0,0,0,0);
      A.clear();
      C.clear();
      // cout << endl;
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    for ( unsigned j=0; j<D.size(); j++){
      //cout << " j the element is " << D.at(j).M()<<endl;

      if (D.at(j).Pt()>100. && D.at(j).Mag()>=120 && D.at(j).Mag()<=240){

	//	cout << " j the element is *********" << D.at(j).Mag()<<endl;
       	TLorentzVector topp4 = D.at(j);
	vlq::Candidate top(topp4) ;
	//	cout << " top mass  is *********" << top.getMass()<<endl;
       	tops.push_back(top) ; 
      }
    }

  }

  
  void operator()(vlq::JetCollection& goodJets, vlq::CandidateCollection& W, vlq::CandidateCollection& B){
    for (unsigned i=0; i<goodJets.size(); i++){
      TLorentzVector jetinfo;
      jetinfo = goodJets.at(i).getP4();
      
      if (jetinfo.Pt()>100. && jetinfo.Mag()>=60 && jetinfo.Mag()<=140){
	
	//	cout << " W candidate jet mass is *********" << jetinfo.M()<<endl;                                                                         
	//	cout << " W candidate jet pt  is *********" << jetinfo.Pt()<<endl;
        
	vlq::Candidate ws(jetinfo);     
	W.push_back(ws) ;
	
      }
      else {
      	//cout << " other candidate jet mass is *********" << jetinfo.M()<<endl;                                                                           
	// cout << " other candidate jet pt  is *********" << jetinfo.Pt()<<endl;                                                                           
     	vlq::Candidate b(jetinfo);
     	B.push_back(b);
	
      }
      
    }
    
  }
  
    
    ~TopCandsProducer () {}  

    // private:
  
}; 
#endif 
