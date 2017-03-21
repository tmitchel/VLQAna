import FWCore.ParameterSet.Config as cms
from Analysis.VLQAna.JetMaker_cfi import *

hh4b = cms.EDFilter("HH4b",
    runno               = cms.InputTag("evtcleaner","runno"), 
    lumisec             = cms.InputTag("evtcleaner","lumisec"), 
    evtno               = cms.InputTag("evtcleaner","evtno"), 
    isData              = cms.InputTag("evtcleaner","isData"), 
    evtwtGen            = cms.InputTag("evtcleaner","evtwtGen"),
    evtwtPV             = cms.InputTag("evtcleaner","evtwtPV"),
    evtwtPVLow          = cms.InputTag("evtcleaner","evtwtPVLow"),
    evtwtPVHigh         = cms.InputTag("evtcleaner","evtwtPVHigh"),
    htHat               = cms.InputTag("evtcleaner","htHat"),
    lhewtids            = cms.InputTag("evtcleaner","lhewtids"),
    lhewts              = cms.InputTag("evtcleaner","lhewts"), 
    hltdecision         = cms.InputTag("evtcleaner","hltdecision"), 
    npv                 = cms.InputTag("evtcleaner","npv"),
    npuTrue             = cms.InputTag("evtcleaner","npuTrue"),
    jetAK8selParams  = defaultAK8CHSJetSelectionParameters.clone(
      jetPtMin                  = cms.double(300),
      JetIDParams  = defaultAK8CHSJetIDParameters.clone(
        quality = cms.string  ("TIGHTLEPVETO"),
        ), 
      jetAbsEtaMax = cms.double(100), 
      JetSubstrParams = defaultCHSJetSubstructureParameters.clone(
        jettau1Label            = cms.InputTag("jetsAK8CHS", "jetAK8CHStau1Puppi"), 
        jettau2Label            = cms.InputTag("jetsAK8CHS", "jetAK8CHStau2Puppi"),  
        jettau3Label            = cms.InputTag("jetsAK8CHS", "jetAK8CHStau3Puppi"),  
        jetSoftDropMassLabel    = cms.InputTag("jetsAK8CHS", "jetAK8CHSsoftDropMassPuppi"),
        ), 
        SDMassCorrWt              = cms.string("puppiCorr.root"),
      ),
    jetHTaggedselParams = defaultCHSHJetSelectionParameters.clone(
      JetIDParams  = defaultAK8CHSJetIDParameters.clone(
        quality = cms.string  ("TIGHTLEPVETO"),
        ), 
      jetSoftDropMassMin    = cms.double(100) ,
      jetSoftDropMassMax    = cms.double(135) ,
      JetSubstrParams = defaultCHSJetSubstructureParameters.clone(
        jettau1Label            = cms.InputTag("jetsAK8CHS", "jetAK8CHStau1Puppi"), 
        jettau2Label            = cms.InputTag("jetsAK8CHS", "jetAK8CHStau2Puppi"),  
        jettau3Label            = cms.InputTag("jetsAK8CHS", "jetAK8CHStau3Puppi"),  
        jetSoftDropMassLabel    = cms.InputTag("jetsAK8CHS", "jetAK8CHSsoftDropMassPuppi"),
        ), 
      SDMassCorrWt              = cms.string("puppiCorr.root"),
      ),
    btageffFile         = cms.string("btag-eff-subjet.root"), 
    printEvtNo          = cms.bool(False),
    )
