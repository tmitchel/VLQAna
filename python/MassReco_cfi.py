import FWCore.ParameterSet.Config as cms

from Analysis.VLQAna.PickGenPart_cfi import *

massReco = cms.EDFilter("MassReco",
    lhewtids = cms.InputTag("evtcleaner", "lhewtids"),
    lhewts   = cms.InputTag("evtcleaner", "lhewts"),
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
    Prewt    = cms.InputTag("ana", "PreWeight"),
    bTagwt   = cms.InputTag("ana", "btagsf"),
    sjbTagwt  = cms.InputTag("ana", "sjbtagsf"),
    st       = cms.InputTag("ana", "st"),
		STMin    = cms.double(1000.),
		STMaxControl = cms.double(1000.),
    ptMin    = cms.double(0.),
    zdecaymode = cms.string('zmumu'),
    signalType = cms.string(''),
		optimizeReco = cms.bool(False),
		controlReco = cms.bool(False),
    doSkim   = cms.bool(False),
    isData   = cms.bool(False),
    pdfID_offset = cms.int32(1),
    scale_offset = cms.int32(1),
    

 		genParams = getSelectParams,	
 
   )
