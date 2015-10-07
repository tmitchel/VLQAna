import FWCore.ParameterSet.Config as cms

from Analysis.VLQAna.ZCandSelector_cfi import *
from Analysis.VLQAna.ElectronSelector_cfi import *
from Analysis.VLQAna.MuonSelector_cfi import *
from Analysis.VLQAna.PickGenPart_cfi import *
from Analysis.VLQAna.JetSelector_cfi import *

ana = cms.EDFilter("OS2LAna", 
    trigNameLabel              = cms.InputTag("TriggerUserData", "triggerNameTree"), 
    trigBitLabel               = cms.InputTag("TriggerUserData", "triggerBitTree"), 
    #metFiltersNameLabel        = cms.InputTag("METUserData", "triggerNameTree"), 
    #metFiltersBitLabel         = cms.InputTag("METUserData", "triggerBitTree"), 
    #hbheNoiseFilterLabel       = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun1"), 
    genEvtInfoProdName         = cms.string('generator'),
    #vtxRhoLabel                = cms.InputTag("vertexInfo", "rho"),
    #vtxZLabel                  = cms.InputTag("vertexInfo", "z"),
    #vtxNdfLabel                = cms.InputTag("vertexInfo", "ndof"),
    npvLabel                   = cms.InputTag("eventUserData", "npv"),
    jetAK8PtLabel              = cms.InputTag("jetsAK8", "jetAK8Pt"),
    jetAK8EtaLabel             = cms.InputTag("jetsAK8", "jetAK8Eta"),
    jetAK8PhiLabel             = cms.InputTag("jetsAK8", "jetAK8Phi"),
    jetAK8MassLabel            = cms.InputTag("jetsAK8", "jetAK8Mass"),
    jetAK8FilteredMassLabel    = cms.InputTag("jetsAK8", "jetAK8filteredMass"),
    jetAK8TrimmedMassLabel     = cms.InputTag("jetsAK8", "jetAK8trimmedMass"),
    jetAK8PrunedMassLabel      = cms.InputTag("jetsAK8", "jetAK8prunedMass"),
    jetAK8SoftDropMassLabel    = cms.InputTag("jetsAK8", "jetAK8softDropMass"),
    jetAK8EnergyLabel          = cms.InputTag("jetsAK8", "jetAK8Energy"),
    jetAK8FlavourLabel         = cms.InputTag("jetsAK8", "jetAK8Flavour"),
    jetAK8CSVLabel             = cms.InputTag("jetsAK8", "jetAK8CSV"),
    jetAK8JECLabel             = cms.InputTag("jetsAK8", "jetAK8jecFactor0"),
    jetAK8AreaLabel            = cms.InputTag("jetsAK8", "jetAK8jetArea"),
    jetAK8Tau1Label            = cms.InputTag("jetsAK8", "jetAK8tau1"), 
    jetAK8Tau2Label            = cms.InputTag("jetsAK8", "jetAK8tau2"),  
    jetAK8Tau3Label            = cms.InputTag("jetsAK8", "jetAK8tau3"),  
    jetAK8nSubJetsLabel        = cms.InputTag("jetsAK8", "jetAK8nSubJets"),  
    jetAK8minmassLabel         = cms.InputTag("jetsAK8", "jetAK8minmass"),  
    jetAK8VSubjetIndex0Label   = cms.InputTag("jetsAK8", "jetAK8VSubjetIndex0"),  
    jetAK8VSubjetIndex1Label   = cms.InputTag("jetsAK8", "jetAK8VSubjetIndex1"),  
    jetAK8TopSubjetIndex0Label = cms.InputTag("jetsAK8", "jetAK8TopSubjetIndex0"),
    jetAK8TopSubjetIndex1Label = cms.InputTag("jetsAK8", "jetAK8TopSubjetIndex1"),
    jetAK8TopSubjetIndex2Label = cms.InputTag("jetsAK8", "jetAK8TopSubjetIndex2"),
    jetAK8TopSubjetIndex3Label = cms.InputTag("jetsAK8", "jetAK8TopSubjetIndex3"),
    subjetAK8BDiscLabel        = cms.InputTag("subjetsAK8", "subjetAK8CSV"),
    subjetAK8PtLabel           = cms.InputTag("subjetsAK8", "subjetAK8Pt"),
    subjetAK8EtaLabel          = cms.InputTag("subjetsAK8", "subjetAK8Eta"),
    subjetAK8PhiLabel          = cms.InputTag("subjetsAK8", "subjetAK8Phi"),
    subjetAK8MassLabel         = cms.InputTag("subjetsAK8", "subjetAK8Mass"),
    subjetCmsTopTagBDiscLabel  = cms.InputTag("subjetsCmsTopTag", "subjetsCmsTopTagCSV"),
    subjetCmsTopTagPtLabel     = cms.InputTag("subjetsCmsTopTag", "subjetsCmsTopTagPt"),
    subjetCmsTopTagEtaLabel    = cms.InputTag("subjetsCmsTopTag", "subjetsCmsTopTagEta"),
    subjetCmsTopTagPhiLabel    = cms.InputTag("subjetsCmsTopTag", "subjetsCmsTopTagPhi"),
    subjetCmsTopTagMassLabel   = cms.InputTag("subjetsCmsTopTag", "subjetsCmsTopTagMass"),
    jetAK4PtLabel              = cms.InputTag("jetsAK4", "jetAK4Pt"),
    jetAK4EtaLabel             = cms.InputTag("jetsAK4", "jetAK4Eta"),
    jetAK4PhiLabel             = cms.InputTag("jetsAK4", "jetAK4Phi"),
    jetAK4MassLabel            = cms.InputTag("jetsAK4", "jetAK4Mass"),
    jetAK4EnergyLabel          = cms.InputTag("jetsAK4", "jetAK4Energy"),
    jetAK4FlavourLabel         = cms.InputTag("jetsAK4", "jetAK4Flavour"),
    jetAK4CSVLabel             = cms.InputTag("jetsAK4", "jetAK4CSV"),
    jetAK4JECLabel             = cms.InputTag("jetsAK4", "jetAK4jecFactor0"),
    jetAK4nHadEnergyLabel      = cms.InputTag("jetsAK4", "jetAK4neutralHadronEnergy"),
    jetAK4nEMEnergyLabel       = cms.InputTag("jetsAK4", "jetAK4neutralEmEnergy"),
    jetAK4HFHadronEnergyLabel  = cms.InputTag("jetsAK4", "jetAK4HFHadronEnergy"),
    jetAK4cHadEnergyLabel      = cms.InputTag("jetsAK4", "jetAK4chargedHadronEnergy"),
    jetAK4cEMEnergyLabel       = cms.InputTag("jetsAK4", "jetAK4chargedEmEnergy"),
    jetAK4numDaughtersLabel    = cms.InputTag("jetsAK4", "jetAK4numberOfDaughters"),
    jetAK4cMultipLabel         = cms.InputTag("jetsAK4", "jetAK4chargedMultiplicity"),
    jetAK4YLabel               = cms.InputTag("jetsAK4", "jetAK4Y"),
    muPt                       = cms.InputTag("muons", "muPt"), 
    muEta                      = cms.InputTag("muons", "muEta"), 
    muY                        = cms.InputTag("muons", "muY"),
    muPhi                      = cms.InputTag("muons", "muPhi"), 
    muMass                     = cms.InputTag("muons", "muMass"), 
    muE                        = cms.InputTag("muons", "muE"), 
    muCharge                   = cms.InputTag("muons", "muCharge"), 
    elPt                       = cms.InputTag("electrons", "elPt"), 
    elEta                      = cms.InputTag("electrons", "elEta"), 
    elY                        = cms.InputTag("electrons", "elY"),
    elPhi                      = cms.InputTag("electrons", "elPhi"), 
    elMass                     = cms.InputTag("electrons", "elMass"), 
    elE                        = cms.InputTag("electrons", "elE"), 
    elCharge                   = cms.InputTag("electrons", "elCharge"), 
    met                        = cms.InputTag("met", "metPt"), 
    metPhi                     = cms.InputTag("met", "metPhi"), 
    #metNoHF                    = cms.InputTag("met", "metPt"), 
    #metNoHFPhi                 = cms.InputTag("met", "metPhi"), 
    HbbCandsLabel              = cms.InputTag("hbb", "HbbCandidates"),
    jetAK4AreaLabel            = cms.InputTag("jetsAK4", "jetAK4jetArea"),
    hltPaths                   = cms.vstring (
        ### Muon paths
        #"HLT_IsoMu20_eta2p1_v",
        #"HLT_Mu30_TkMu11_v", 
        #"HLT_Mu45_eta2p1_v", 
        #"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"
        ### Electron paths
        "HLT_Ele27_eta2p1_WPLoose_Gsf_v", 
        ), 
    metFilters                 = cms.vstring (
        "Flag_goodVertices",
        "Flag_CSCTightHaloFilter",
        "Flag_HBHENoiseFilter", 
        ), 
    DilepCandParams            = defaultZCandSelectionParameters.clone(), 
    ZCandParams                = defaultZCandSelectionParameters.clone(
        massMin = cms.double(60),
        massMax = cms.double(120),
        ), 
    BoostedZCandParams         = defaultZCandSelectionParameters.clone(
        massMin = cms.double(60),
        massMax = cms.double(120),
        ptMin = cms.double(150.),
        ), 
    GenHSelParams              = genPartParams.clone(), 
    AK4JetSelParams            = defaultAK4JetSelectionParameters.clone(), 
    BTaggedLooseAK4SelParams   = defaultBTaggedAK4JetSelectionParameters.clone(jetCSVDiscMin = cms.double(0.423),), 
    BTaggedMediumAK4SelParams   = defaultBTaggedAK4JetSelectionParameters.clone(), 
    mupselParams = defaultMuPSelectionParamaters.clone(), 
    mumselParams = defaultMuMSelectionParamaters.clone(), 
    elpselParams = defaultElPSelectionParamaters.clone(), 
    elmselParams = defaultElMSelectionParamaters.clone(), 
    AK8JetSelParams            = defaultAK8JetSelectionParameters.clone(
        jetPtMin   = cms.double(300),
        jetMassMin = cms.double(0) ,
        ), 
    TJetSelParams              = defaultTJetSelectionParameters.clone(), 
    HJetSelParams              = defaultHJetSelectionParameters.clone(), 
    WJetSelParams              = defaultWJetSelectionParameters.clone(), 
    HTMin                      = cms.double  (300.), 
    isData                     = cms.bool    (False), 
    )

