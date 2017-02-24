import FWCore.ParameterSet.Config as cms

from Analysis.VLQAna.PickGenPart_cfi import *

massReco = cms.EDFilter("MassReco",
    met      = cms.InputTag("metFull", "metFullPt"),
    elPt     = cms.InputTag("electrons", "elPt"),
    elEta    = cms.InputTag("electrons", "elEta"),
    elPhi    = cms.InputTag("electrons", "elPhi"),
    elE      = cms.InputTag("electrons", "elE"),
    muPt     = cms.InputTag("muons", "muPt"),
    muEta    = cms.InputTag("muons", "muEta"),
    muPhi    = cms.InputTag("muons", "muPhi"),
    muE      = cms.InputTag("muons", "muE"),
		ak8jets  = cms.InputTag("ana", "ak8jets"),
    jets     = cms.InputTag("ana", "jets"),
    bjets    = cms.InputTag("ana", "bjets"),
    zjets    = cms.InputTag("ana", "wjets"),
    hjets    = cms.InputTag("ana", "hjets"),
    zllcands = cms.InputTag("ana", "zllcands"),
    evtwt    = cms.InputTag("ana", "PreWeight"),
    st       = cms.InputTag("ana", "st"),
		STMin    = cms.double(1000.),
		STMaxControl = cms.double(700.),
    ptMin    = cms.double(0.),
    zdecaymode = cms.string('zmumu'),
    signalType = cms.string(''),
		optimizeReco = cms.bool(True),
		controlReco = cms.bool(False),
 		genParams = getSelectParams,	
 
   )
