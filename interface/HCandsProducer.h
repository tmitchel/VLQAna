#ifndef ANALYSIS_VLQANA_HCANDSPRODUCER_HH
#define ANALYSIS_VLQANA_HCANDSPRODUCER_HH

#include <iostream>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisDataFormats/BoostedObjects/interface/Candidate.h"

using namespace std;

class HCandsProducer {
  public:


  void operator()(int N, int K, vlq::JetCollection& goodJets, vlq::CandidateCollection& H)
  {
    std::string bitmask(K, 1); // K leading 1's                                                                                                               
    bitmask.resize(N, 0); // N-K trailing 0'                                                                                                                  
    double sum = 0;
    TLorentzVector jetinfo;
    jetinfo.SetPtEtaPhiM(0,0,0,0);
    vector<double> A, B;
    vector<TLorentzVector> C,D; 
    vlq::JetCollection E;
    // print integers and permute bitmask                                                                                                                     
    do {
      for (int i = 0; i < N; ++i) // [0..N-1] integers                                                                                                        
	{
	  
	  if (bitmask[i]) {//cout << " " << i; 
	    //cout << i<<"th element CSV in D is  = "<< goodJets.at(i).getCSV()<<endl;
	    sum +=  goodJets.at(i).getMass();
	    jetinfo += goodJets.at(i).getP4();
	    // cout<< " jet Mass is " << goodJets.at(i).getMass()<<endl;
	    
	    // cout<< " the mass sum is = " << sum<<endl;
	    A.push_back(sum);
	    C.push_back(jetinfo);
	    E.push_back(goodJets.at(i));
	  }
	  
	  //  if(A[2]<120. && A[2]>240.){continue;} 
	}
     
      //cout<< "Total mass of 2 jets  is  = " << A[1]<<endl;
      //cout<< "Total mass of 2 jets  ****T LORENTZ ****  = " << C.at(1).Mag()<<endl;

      //making sure that from the sum at least one jet is batagged
      // cout << "1 st element CSV in D is  = "<< D.at(0).getCSV()<<endl; 
      //  E.push_back(goodJets.at(i));

      //cout << " jet collection size is " << E.size()<<endl;
      //for(unsigned i=0; i<E.size();i++){
      //cout << i <<"th jet CSV is =" << E.at(i).getCSV()<<endl;
      // }
       //making sure at least one jet is btagged
      if ((E.at(0).getCSV()>0.8484 && E.at(0).getPt()>50)|| (E.at(1).getCSV()>0.8484 && E.at(1).getPt()>50)){
      // if (E.at(0).getPt()>50 ||  E.at(1).getPt()>50){ 
      B.push_back(A[1]);
      D.push_back(C.at(1)); 
      //cout <<" 1st jet CSV is =" << E.at(0).getCSV()<<endl;
      //cout <<" 2nd jet CSV is =" << E.at(1).getCSV()<<endl;
      //cout<< "Total mass of 2 jets  is  %%%%%%%%%%%%%%%%%%%%%= " << A[1]<<endl;
      //cout<< "Total mass of 2 jets  ****T LORENTZ ****%%%%%%%%%%%%  = " << C.at(1).Mag()<<endl;
      }
        
//if(A[2]<120. && A[2]>240.){cout<< "Total mass of 3 jets  is ******  = " << A[2]<<endl; B.push_back(A[2]);}
      sum =0;
      jetinfo.SetPtEtaPhiM(0,0,0,0);
      A.clear();
      C.clear();
      E.clear();
      // cout << endl;
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    for ( unsigned j=0; j<D.size(); j++){
      // cout << " j the element is " << D.at(j).M()<<endl;

      if (D.at(j).Pt()>100. && D.at(j).Mag()>=80 && D.at(j).Mag()<=160){

	//	cout << " j the element is *********" << D.at(j).Mag()<<endl;
       	TLorentzVector Hp4 = D.at(j);
	vlq::Candidate h(Hp4) ;
	//	cout << " higgs mass  is *********" << h.getMass()<<endl;
       	H.push_back(h) ; 
      }
    }

  }

  
    
    ~HCandsProducer () {}  

    // private:
  
}; 
#endif 
