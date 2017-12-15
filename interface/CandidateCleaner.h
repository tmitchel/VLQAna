#ifndef ANALYSIS_VLQANA_CANDIDATECLEANER_HH
#define ANALYSIS_VLQANA_CANDIDATECLEANER_HH
#include <vector>
#include <algorithm>
#include "AnalysisDataFormats/BoostedObjects/interface/Candidate.h"

///
/// Class to clean one set of candidates T1 from another T2:
/// Configurables:
///   dr_ = Separation in DR(eta, phi) between the two cadidates.
///   ptrel_ = relative pt between the two candidates (used for lepton 2D isolation). 
///   If ptrel_ < 0 (negative), 2D isolation is not applied, and an ordinary DR cleaning is used.
///   If ptrel_ > 0, T1 collection are cleaned from T2 if dr(t1, t2) > dr_ *or* if pt(t1)/pt(t2) > ptrel_
///

class CandidateCleaner {
  public:
    CandidateCleaner(double dr, double ptrel) : dr_(dr), ptrel_(ptrel) {}
    ~CandidateCleaner() {}
    template <class T1, class T2>
      void operator () (std::vector<T1>& cleanedcands, std::vector<T2> othercands) {

        cleanedcands.erase( std::remove_if(cleanedcands.begin(), cleanedcands.end(),
              [othercands, this](T1 t) { 
                bool toremove(false);
                for (auto cand : othercands) {
                  double dr((t.getP4()).DeltaR(cand.getP4())) ;
                  double ptrel( ( ((t.getP4()).Vect()).Cross((cand.getP4()).Vect()) ).Mag()/(cand.getP4()).Mag() ) ; 
                  if ( dr < dr_ && ( ptrel_ > 0 ? ptrel < ptrel_ : true)  ) { toremove = true; return toremove; }
                }
                return toremove; 
              }
              )
            ,cleanedcands.end()  ) ; 
      }

  private:
    double dr_;
    double ptrel_;
}; 
#endif
