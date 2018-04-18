import FWCore.ParameterSet.Config as cms

from Analysis.VLQAna.ZCandSelector_cfi import *
from Analysis.VLQAna.ElectronMaker_cfi import *
from Analysis.VLQAna.MuonMaker_cfi import *
from Analysis.VLQAna.JetMaker_cfi import *
from Analysis.VLQAna.PickGenPart_cfi import *


trigAna = cms.EDFilter("trigAna",
    probe_evtno         = cms.InputTag("probe", "evtno"),
    probe_runno         = cms.InputTag("probe", "runno"),
    probe_lumisec       = cms.InputTag("probe", "lumisec"),
    probe_hltdec        = cms.InputTag("probe", "hltdecision"),

    tag_evtno           = cms.InputTag("tag", "evtno"),
    tag_runno           = cms.InputTag("tag", "runno"),
    tag_lumisec         = cms.InputTag("tag", "lumisec"),
    tag_hltdec          = cms.InputTag("tag", "hltdecision"),

    reject_hltdec       = cms.InputTag("reject", "hltdecision"),

    DilepCandParams     = defaultZCandSelectionParameters.clone(
        massMin         = cms.double(60),
        massMax         = cms.double(120),
        ptMaxLeadingLep = cms.double(120),
        ptMin           = cms.double(0.),
        ),
    elselParams         = defaultElectronMakerParameters.clone(
        elPtMin         = cms.double(35),
        ),
    muselParams         = defaultMuonMakerParameters.clone(
        muPtMin         = cms.double(30),
        ),
    jetAK4selParams     = defaultAK4JetSelectionParameters.clone(
        jetAbsEtaMax    = cms.double(2.4),
        ),

    isData              = cms.bool(False),
    isPhoton            = cms.bool(False),
    htTrig              = cms.bool(False),
    muTrig              = cms.bool(False),
    genParams           = getSelectParams,

    )
