import FWCore.ParameterSet.Config as cms

from Analysis.VLQAna.PickGenPart_cfi import *
from Analysis.VLQAna.JetMaker_cfi import *
from Analysis.VLQAna.ElectronMaker_cfi import *
from Analysis.VLQAna.MuonMaker_cfi import *

ana = cms.EDFilter("VLQAna", 
    runno                      = cms.InputTag("evtcleaner","runno"), 
    lumisec                    = cms.InputTag("evtcleaner","lumisec"), 
    evtno                      = cms.InputTag("evtcleaner","evtno"), 
    isData                     = cms.InputTag("evtcleaner","isData"), 
    hltdecision                = cms.InputTag("evtcleaner","hltdecision"), 
    evttype                    = cms.InputTag("evtcleaner","evttype"),
    evtwtGen                   = cms.InputTag("evtcleaner","evtwtGen"),
    evtwtPV                    = cms.InputTag("evtcleaner","evtwtPV"),
    evtwtPVLow                 = cms.InputTag("evtcleaner","evtwtPVLow"),
    evtwtPVHigh                = cms.InputTag("evtcleaner","evtwtPVHigh"),
    npv                        = cms.InputTag("evtcleaner","npv"),
    npuTrue                    = cms.InputTag("evtcleaner","npuTrue"),
    htHat                      = cms.InputTag("evtcleaner","htHat"),
    lhewtids                   = cms.InputTag("evtcleaner","lhewtids"),
    lhewts                     = cms.InputTag("evtcleaner","lhewts"), 
    elselParams                = defaultElectronMakerParameters.clone(
      elPtMin = cms.double(50),
      applyIso = cms.bool(False), 
      ), 
    muselParams                = defaultMuonMakerParameters.clone(
      muidtype = cms.string("LOOSE"), 
      muPtMin = cms.double(47),
      muIsoMin = cms.double(0.00),
      muIsoMax = cms.double(1000), 
      ), 
    jetAK4selParams            = defaultAK4JetSelectionParameters,
    jetAK8selParams            = defaultAK8CHSJetSelectionParameters,
    jetHTaggedselParams        = defaultCHSHJetSelectionParameters,
    jetAntiHTaggedselParams    = defaultCHSHJetSelectionParameters.clone(
      subjetCSVMin = cms.double(-1000000) ,
      subjetCSVMax = defaultCHSHJetSelectionParameters.subjetCSVMin,
      ), 
    jetTopTaggedselParams      = defaultCHSTJetSelectionParameters,
    jetAntiTopTaggedselParams  = defaultCHSTJetSelectionParameters.clone(
      subjetHighestCSVMin = cms.double(-1000000),
      subjetHighestCSVMax = defaultCHSTJetSelectionParameters.subjetHighestCSVMin, 
      ),
    leadingJetPtMin            = cms.double  (400.), 
    leadingJetPrunedMassMin    = cms.double  (50.), 
    HTMin                      = cms.double  (0.), 
    doBTagSFUnc                = cms.bool(False), 
    storePreselEvts            = cms.bool(False), 
    doPreselOnly               = cms.bool(False), 
    )
