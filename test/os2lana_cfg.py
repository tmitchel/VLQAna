import os, sys
from readVars import *
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('isData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is data?"
    )
options.register('skim', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Skim events?"
    )
options.register('newJECPayloadNames', 'Summer16_23Sep2016V4_MC',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Select one of EvtType_MC_tZtZ, EvtType_MC_tZtH, EvtType_MC_tZbW, EvtType_MC_tHtH, EvtType_MC_tHbW, EvtType_MC_bWbW, EvtType_MC_bZbZ, EvtType_MC_bZbH, EvtType_MC_bZtW, EvtType_MC_bHbH, EvtType_MC_bHtW, EvtType_MC_tWtW" 
    )
options.register('maketree', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Skim events?"
    )
options.register('zdecaymode', 'zmumu',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Z->mumu or Z->elel? Choose: 'zmumu' or 'zelel'"
    )
options.register('lepID', 'TIGHT',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "lepton ID? Choose: 'TIGHT' or 'LOOSE'"
    )
options.register('outFileName', 'os2lana.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('doPUReweightingOfficial', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do pileup reweighting using official recipe"
    )
options.register('filterSignal', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Select only tZXX or bZXX modes"
    )
options.register('signalType', 'EvtType_MC_bZbZ',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Select one of EvtType_MC_tZtZ, EvtType_MC_tZtH, EvtType_MC_tZbW, EvtType_MC_tHtH, EvtType_MC_tHbW, EvtType_MC_bWbW, EvtType_MC_bZbZ, EvtType_MC_bZbH, EvtType_MC_bZtW, EvtType_MC_bHbH, EvtType_MC_bHtW, EvtType_MC_tWtW" 
    )
options.register('applyLeptonIDSFs', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply lepton SFs to the MC"
    )
options.register('applyBTagSFs', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply b-tagging SFs to the MC"
    )
options.register('btageffmap', "bpbpbZb800_bTagEff.root",#until new SFs arrive
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "ROOT file with Th2D histos of b tag effs for b,c, and light flavoured jets"
    )
options.register('applyLeptonTrigSFs', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply trigger SFs to the MC"
    )
options.register('applyDYNLOCorr', False, ### Set to true only for DY process ### Only EWK NLO k-factor is applied
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply DY EWK k-factor to DY MC"
    )
options.register('FileNames', 'FileNames_DY',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of list of input files"
    )
options.register('massReco', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Optimize mass reconstruction"
    )
options.register('controlReco', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'Do control region reco'
    )
options.register('syst', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do systematics"
    )
options.register('storeLHEWts', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Store LHE weights",
    )
options.register('optimizeReco', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'meh'
    )

options.setDefault('maxEvents', -1)
options.parseArguments()

dataPath = '../data/'

hltpaths = []
if options.isData:
  options.filterSignal = False 
  options.signalType = "" 
  options.applyLeptonIDSFs = False 
  options.applyLeptonTrigSFs= False
  options.applyBTagSFs   = False 
  options.applyDYNLOCorr = False 
  options.doPUReweightingOfficial = False

if options.zdecaymode == "zmumu":
  hltpaths = [
      #"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"
      "HLT_IsoMu24_v",
      "HLT_IsoTkMu24_v"
      ]
elif options.zdecaymode == "zelel":
  hltpaths = [
      #"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      #"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
      ]
else:
  sys.exit("!!!Error: Wrong Z decay mode option chosen. Choose either 'zmumu' or 'zelel'!!!") 

if options.skim: 
  hltpaths = []

if options.filterSignal == True: 
  print 'signal type = ', len(options.signalType), 'skim : ', options.skim
  if options.skim or options.maketree:
    if len(options.signalType) != 0: sys.exit("!!!Error: Please do not specify any signal type when skimming the signal MC!!!")
  elif len(options.signalType) == 0:
    sys.exit("!!!Error: Cannot keep signalType empty when filterSignal switched on!!!") 

print options

process = cms.Process("OS2LAna")

from inputFiles_cfi import * 

process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring(
    'root://eoscms.cern.ch//store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192852/0000/B2GEDMNtuple_1.root'
    #FileNames[options.FileNames]
 ) 
  )

if options.isData:
  import FWCore.PythonUtilities.LumiList as LumiList
  process.source.lumisToProcess = LumiList.LumiList(
      filename = "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
      ).getVLuminosityBlockRange()

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
      options.outFileName
      )
    )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options  = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Analysis.VLQAna.EventCleaner_cff")
process.evtcleaner.File_PUDistData      = cms.string(os.path.join(dataPath,'RunII2016Rereco_25ns_PUXsec69000nb.root'))
process.evtcleaner.File_PUDistDataLow   = cms.string(os.path.join(dataPath,'RunII2016Rereco_25ns_PUXsec65550nb.root'))
process.evtcleaner.File_PUDistDataHigh  = cms.string(os.path.join(dataPath,'RunII2016Rereco_25ns_PUXsec72450nb.root'))
process.evtcleaner.File_PUDistMC        = cms.string(os.path.join(dataPath,'PUDistMC_Summer2016_25ns_Moriond17MC_PoissonOOTPU.root'))
process.evtcleaner.isData = options.isData 
process.evtcleaner.hltPaths = cms.vstring (hltpaths)  
process.evtcleaner.DoPUReweightingOfficial = cms.bool(options.doPUReweightingOfficial)  
process.evtcleaner.storeLHEWts = options.storeLHEWts

from Analysis.VLQAna.OS2LAna_cfi import * 

### Z candidate and jet selections 
process.ana = ana.clone(
    filterSignal = cms.bool(options.filterSignal),
    signalType = cms.string(options.signalType),
    zdecayMode = cms.string(options.zdecaymode),
    applyLeptonIDSFs = cms.bool(options.applyLeptonIDSFs),
    applyLeptonTrigSFs = cms.bool(options.applyLeptonTrigSFs),
    applyBTagSFs = cms.bool(options.applyBTagSFs),
    btageffmap = cms.string(os.path.join(dataPath,options.btageffmap)), 
    applyDYNLOCorr = cms.bool(options.applyDYNLOCorr),
    skim = cms.bool(options.skim),
    isData = cms.bool(options.isData),
    maketree = cms.bool(options.maketree), 
    fnamebtagSF = cms.string(os.path.join(dataPath,'CSVv2_Moriond17_B_H.csv')),
    fnameSJbtagSF = cms.string(os.path.join(dataPath,'subjet_CSVv2_Moriond17_B_H.csv')),
    File_DYNLOCorr = cms.string(os.path.join(dataPath,'scalefactors_v4.root')),
    )
process.ana.genParams.debug = cms.bool(False)

process.ana.jetAK4selParams.jecUncPayloadName = cms.string(os.path.join(dataPath,options.newJECPayloadNames+"_Uncertainty_AK4PFchs.txt"))
if options.maketree: 
  process.ana.jetAK4BTaggedselParams.jetCSVDiscMin = cms.double(CSVv2L)

process.ana.jetAK4BTaggedselParams.jecUncPayloadName = cms.string(os.path.join(dataPath,options.newJECPayloadNames+"_Uncertainty_AK4PFchs.txt"))
process.ana.jetAK8selParams.jecUncPayloadName = cms.string(os.path.join(dataPath,options.newJECPayloadNames+"_Uncertainty_AK8PFchs.txt"))
process.ana.jetHTaggedselParams.jecUncPayloadName = cms.string(os.path.join(dataPath,options.newJECPayloadNames+"_Uncertainty_AK8PFchs.txt"))
process.ana.jetWTaggedselParams.jecUncPayloadName = cms.string(os.path.join(dataPath,options.newJECPayloadNames+"_Uncertainty_AK8PFchs.txt"))
process.ana.jetTopTaggedselParams.jecUncPayloadName = cms.string(os.path.join(dataPath,options.newJECPayloadNames+"_Uncertainty_AK8PFchs.txt"))
masscorr = ['L2Relative', 'L3Absolute']
for i,c in enumerate(masscorr): masscorr[i] = os.path.join(dataPath,options.newJECPayloadNames+"_"+c+"_AK8PFchs.txt")
print "jec payload ", masscorr
process.ana.jetAK8selParams.jecAK8GroomedPayloadNames = cms.vstring(masscorr)
process.ana.jetHTaggedselParams.jecAK8GroomedPayloadNames = cms.vstring(masscorr)
process.ana.jetWTaggedselParams.jecAK8GroomedPayloadNames = cms.vstring(masscorr)
process.ana.jetTopTaggedselParams.jecAK8GroomedPayloadNames = cms.vstring(masscorr)
if options.newJECPayloadNames != '':
  corrections = ['L1FastJet', 'L2Relative', 'L3Absolute']
  ak4chsCorr = []
  for c in corrections: ak4chsCorr.append(os.path.join(dataPath,options.newJECPayloadNames+"_"+c+"_AK4PFchs.txt"))
  ak8chsCorr = []
  for c in corrections: ak8chsCorr.append(os.path.join(dataPath,options.newJECPayloadNames+"_"+c+"_AK8PFchs.txt"))
  print "ak4chsCorr ", ak4chsCorr
  print "ak8chsCorr ", ak8chsCorr
  process.ana.jetAK4selParams.newJECPayloadNames = cms.vstring(ak4chsCorr)
  process.ana.jetAK4BTaggedselParams.newJECPayloadNames = cms.vstring(ak4chsCorr)
  process.ana.jetAK8selParams.newJECPayloadNames = cms.vstring(ak8chsCorr)
  process.ana.jetHTaggedselParams.newJECPayloadNames = cms.vstring(ak8chsCorr)
  process.ana.jetWTaggedselParams.newJECPayloadNames = cms.vstring(ak8chsCorr)
  process.ana.jetTopTaggedselParams.newJECPayloadNames = cms.vstring(ak8chsCorr)

process.ana.elselParams.elidtype = cms.string(options.lepID)
process.ana.muselParams.muidtype = cms.string(options.lepID)
process.ana.muselParams.muIsoMax = cms.double(0.15)
process.ana.lepIdSFsParams.lepidtype = cms.string(options.lepID)
process.ana.lepIdSFsParams.zdecayMode = cms.string(options.zdecaymode)
process.ana.lepTrigSFsParams.zdecayMode = cms.string(options.zdecaymode)
if options.zdecaymode == "zelel": process.ana.lepTrigSFsParams.eltrigeffmap = cms.string(os.path.join(dataPath,"ElectronTriggerSF.root")) 
if options.zdecaymode == "zelel": 
  process.ana.DilepCandParams.ptMaxLeadingLep = cms.double(120)
  process.ana.ZCandParams.ptMaxLeadingLep = cms.double(120)
elif options.zdecaymode == "zmumu": 
  process.ana.DilepCandParams.ptMaxLeadingLep = cms.double(45)
  process.ana.ZCandParams.ptMaxLeadingLep = cms.double(45)
process.ana.ZCandParams.ptMin = cms.double(100.)
process.ana.jetAK8selParams.jetPtMin = cms.double(200) 
process.ana.jetAK4BTaggedselParams.jetPtMin = cms.double(50) 
process.ana.NAK4Min = cms.uint32(3)
if options.maketree:
  process.ana.STMin = cms.double(0.)
  process.ana.HTMin = cms.double(500.)
else: 
  process.ana.STMin = cms.double(1000.)
  process.ana.HTMin = cms.double(200.)

process.ana.pdfID_offset = cms.int32(9)
process.ana.scale_offset = cms.int32(0)

process.load('Analysis.VLQAna.MassReco_cfi')
process.massReco.ptMin = cms.double(150.)
process.massReco.zdecaymode = cms.string(options.zdecaymode)
process.massReco.signalType = cms.string(options.signalType)
process.massReco.controlReco = cms.bool(options.controlReco)
process.massReco.optimizeReco = cms.bool(options.optimizeReco)
process.massReco.isData = cms.bool(options.isData)

if options.skim: 
  process.ana.jetAK4selParams.jetPtMin = cms.double(20) 
  process.ana.jetAK8selParams.jetPtMin = cms.double(170) 
  process.ana.jetWTaggedselParams.jetPtMin = cms.double(170) 
  process.ana.jetHTaggedselParams.jetPtMin = cms.double(170) 
  process.ana.STMin = cms.double(0.)
else: 
  if options.syst: 
    process.anabcUp = process.ana.clone(
        btagsf_bcUp = cms.bool(True),
        sbtagsf_bcUp = cms.bool(True),
        )
    process.anabcDown = process.ana.clone(
        btagsf_bcDown = cms.bool(True),
        sbtagsf_bcDown = cms.bool(True),
        )
    process.analightUp = process.ana.clone(
        btagsf_lUp = cms.bool(True),
        sbtagsf_lUp = cms.bool(True),
        )
    process.analightDown = process.ana.clone(
        btagsf_lDown = cms.bool(True),
        sbtagsf_lDown = cms.bool(True),
        )
    
    process.anaJecUp = process.ana.clone()
    process.anaJecUp.jetAK4selParams.jecShift = cms.double(1.)
    process.anaJecUp.jetAK4BTaggedselParams.jecShift = cms.double(1.)
    process.anaJecUp.jetAK8selParams.jecShift = cms.double(1.)
    process.anaJecUp.jetHTaggedselParams.jecShift = cms.double(1.)
    process.anaJecUp.jetWTaggedselParams.jecShift = cms.double(1.)
    process.anaJecUp.jetTopTaggedselParams.jecShift = cms.double(1.)

    process.anaJecDown = process.ana.clone()
    process.anaJecDown.jetAK4selParams.jecShift = cms.double(-1.)
    process.anaJecDown.jetAK4BTaggedselParams.jecShift = cms.double(-1.)
    process.anaJecDown.jetAK8selParams.jecShift = cms.double(-1.)
    process.anaJecDown.jetHTaggedselParams.jecShift = cms.double(-1.)
    process.anaJecDown.jetWTaggedselParams.jecShift = cms.double(-1.)
    process.anaJecDown.jetTopTaggedselParams.jecShift = cms.double(-1.)

    process.anaJerUp = process.ana.clone()
    process.anaJerUp.jetAK4selParams.jerShift = cms.int32(2)
    process.anaJerUp.jetAK4BTaggedselParams.jerShift = cms.int32(2)
    process.anaJerUp.jetAK8selParams.jerShift = cms.int32(2)
    process.anaJerUp.jetHTaggedselParams.jerShift = cms.int32(2)
    process.anaJerUp.jetWTaggedselParams.jerShift = cms.int32(2)
    process.anaJerUp.jetTopTaggedselParams.jerShift = cms.int32(2)

    process.anaJerDown = process.ana.clone()
    process.anaJerDown.jetAK4selParams.jerShift = cms.int32(-1)
    process.anaJerDown.jetAK4BTaggedselParams.jerShift = cms.int32(-1)
    process.anaJerDown.jetAK8selParams.jerShift = cms.int32(-1)
    process.anaJerDown.jetHTaggedselParams.jerShift = cms.int32(-1)
    process.anaJerDown.jetWTaggedselParams.jerShift = cms.int32(-1)
    process.anaJerDown.jetTopTaggedselParams.jerShift = cms.int32(-1)
    
    process.anaJmrUp = process.ana.clone()
    process.anaJmrUp.jetHTaggedselParams.jmrShift = cms.int32(2)
    process.anaJmrUp.jetWTaggedselParams.jmrShift = cms.int32(2)

    process.anaJmrDown = process.ana.clone()
    process.anaJmrDown.jetHTaggedselParams.jmrShift = cms.int32(-1)
    process.anaJmrDown.jetWTaggedselParams.jmrShift = cms.int32(-1)

    process.anaPileupUp = process.ana.clone(
        PileupUp = cms.bool(True),
        )
    process.anaPileupDown = process.ana.clone(
        PileupDown = cms.bool(True),
        )
    process.anaDYDown = process.ana.clone(
        DYDown = cms.bool(True),
        )

    process.anatauUp = process.ana.clone(
        tauShift = cms.int32(1),
        )

    process.anatauDown = process.ana.clone(
        tauShift = cms.int32(-1),
        )

if options.massReco and options.syst:
  process.recobcUp = process.massReco.clone(
      bTagwt = cms.InputTag("ana", "btagsfbcUp"),
      sjbTagwt = cms.InputTag("ana", "sjbtagsfbcUp")
      )
  process.recobcDown = process.massReco.clone(
      bTagwt = cms.InputTag("ana", "btagsfbcDown"),
      sjbTagwt = cms.InputTag("ana", "sjbtagsfbcDown"),
      )
  process.recolightUp = process.massReco.clone(
      bTagwt = cms.InputTag("ana", "btagsflUp"),
      sjbTagwt = cms.InputTag("ana", "sjbtagsflUp"),
      )
  process.recolightDown = process.massReco.clone(
      bTagwt = cms.InputTag("ana", "btagsflDown"),
      sjbTagwt = cms.InputTag("ana", "sjbtagsflDown"),
      )
  process.recoJecUp = process.massReco.clone(
      ak8jets = cms.InputTag("anaJecUp", "ak8jets"),
      jets  = cms.InputTag("anaJecUp", "jets"),
      bjets = cms.InputTag("anaJecUp", "bjets"),
      zjets = cms.InputTag("anaJecUp", "wjets"),
      hjets = cms.InputTag("anaJecUp", "hjets"),
      Prewt = cms.InputTag("anaJecUp", "PreWeight"),
      )
  process.recoJecDown = process.massReco.clone(
      ak8jets = cms.InputTag("anaJecDown","ak8jets"),
      jets  = cms.InputTag("anaJecDown","jets"),
      bjets = cms.InputTag("anaJecDown", "bjets"),
      zjets = cms.InputTag("anaJecDown", "wjets"),
      hjets = cms.InputTag("anaJecDown", "hjets"),
      Prewt = cms.InputTag("anaJecDown", "PreWeight"),
      )
  process.recoJerUp = process.massReco.clone(
      ak8jets = cms.InputTag("anaJerUp","ak8jets"),
      jets  = cms.InputTag("anaJerUp","jets"),
      bjets = cms.InputTag("anaJerUp", "bjets"),
      zjets = cms.InputTag("anaJerUp", "wjets"),
      hjets = cms.InputTag("anaJerUp", "hjets"),
      Prewt = cms.InputTag("anaJerUp", "PreWeight"),
      )
  process.recoJerDown = process.massReco.clone(
      ak8jets = cms.InputTag("anaJerDown","ak8jets"),
      jets  = cms.InputTag("anaJerDown","jets"),
      bjets = cms.InputTag("anaJerDown", "bjets"),
      zjets = cms.InputTag("anaJerDown", "wjets"),
      hjets = cms.InputTag("anaJerDown", "hjets"),
      Prewt = cms.InputTag("anaJerDown", "PreWeight"),
      )
  process.recoJmrUp = process.massReco.clone(
      ak8jets = cms.InputTag("anaJmrUp","ak8jets"),
      jets  = cms.InputTag("anaJmrUp","jets"),
      bjets = cms.InputTag("anaJmrUp", "bjets"),
      zjets = cms.InputTag("anaJmrUp", "wjets"),
      hjets = cms.InputTag("anaJmrUp", "hjets"),
      Prewt = cms.InputTag("anaJmrUp", "PreWeight"),
      )
  process.recoJmrDown = process.massReco.clone(
      ak8jets = cms.InputTag("anaJmrDown","ak8jets"),
      jets  = cms.InputTag("anaJmrDown","jets"),
      bjets = cms.InputTag("anaJmrDown", "bjets"),
      zjets = cms.InputTag("anaJmrDown", "wjets"),
      hjets = cms.InputTag("anaJmrDown", "hjets"),
      Prewt = cms.InputTag("anaJmrDown", "PreWeight"),
      )
  
  process.recoPileupUp    = process.massReco.clone(Prewt = cms.InputTag("anaPileupUp", "PreWeight"))
  process.recoPileupDown  = process.massReco.clone(Prewt = cms.InputTag("anaPileupDown", "PreWeight"))
  process.recotauUp       = process.massReco.clone(Prewt = cms.InputTag("anatauUp", "PreWeight"))
  process.recotauDown     = process.massReco.clone(Prewt = cms.InputTag("anatauDown", "PreWeight"))
  process.recoDYDown      = process.massReco.clone(Prewt = cms.InputTag("anaDYdown", "PreWeight"))


## Event counters
from Analysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone(isData=options.isData)
process.cleanedEvents = eventCounter.clone(isData=options.isData)
process.finalEvents = eventCounter.clone(isData=options.isData)

process.load("Analysis.VLQAna.VLQCandProducer_cff")

if options.syst and not options.skim and not options.massReco:
  process.p = cms.Path(
    process.allEvents
    *process.evtcleaner
    *process.cleanedEvents
    *cms.ignore(
      process.ana
      *process.finalEvents
      )
    *cms.ignore(process.anabcUp)
    *cms.ignore(process.anabcDown)
    *cms.ignore(process.analightUp)
    *cms.ignore(process.analightDown)
    *cms.ignore(process.anaJecUp)
    *cms.ignore(process.anaJecDown)
    *cms.ignore(process.anaJerUp)
    *cms.ignore(process.anaJerDown)
    *cms.ignore(process.anaPileupUp)
    *cms.ignore(process.anaPileupDown)
    *cms.ignore(process.anaDYDown)
    *cms.ignore(process.anatauUp)
    *cms.ignore(process.anatauDown) 
    *cms.ignore(process.anaJmrUp)
    *cms.ignore(process.anaJmrDown)

    )
elif options.massReco:
  if options.syst:
    process.p = cms.Path(
      process.allEvents
      *process.evtcleaner
      *process.cleanedEvents
      *cms.ignore(process.ana)
      
      *cms.ignore(process.anabcUp)
      *cms.ignore(process.anabcDown)
      *cms.ignore(process.analightUp)
      *cms.ignore(process.analightDown)
      *cms.ignore(process.anaJecUp)
      *cms.ignore(process.anaJecDown)
      *cms.ignore(process.anaJerUp)
      *cms.ignore(process.anaJerDown)
      *cms.ignore(process.anaPileupUp)
      *cms.ignore(process.anaPileupDown)
      *cms.ignore(process.anaDYDown)
      *cms.ignore(process.anatauUp)
      *cms.ignore(process.anatauDown)
      *cms.ignore(process.anaJmrUp)
      *cms.ignore(process.anaJmrDown)

      *cms.ignore(process.recobcUp)
      *cms.ignore(process.recobcDown)
      *cms.ignore(process.recolightUp)
      *cms.ignore(process.recolightDown)
      *cms.ignore(process.recoJecUp)
      *cms.ignore(process.recoJecDown)
      *cms.ignore(process.recoJerUp)
      *cms.ignore(process.recoJerDown)
      *cms.ignore(process.recoPileupUp)
      *cms.ignore(process.recoPileupDown)
      *cms.ignore(process.recotauUp)
      *cms.ignore(process.recotauDown) 
      *cms.ignore(process.recoJmrUp)
      *cms.ignore(process.recoJmrDown)

      *process.massReco
      *process.finalEvents

      )
  else:
    process.p = cms.Path(
      process.allEvents
      *process.evtcleaner
      *process.cleanedEvents
      *process.ana
      *process.massReco
      *process.finalEvents
      ) 
else:
  process.p = cms.Path(
    process.allEvents
    *process.evtcleaner
    *process.cleanedEvents
    *process.ana
    *process.finalEvents
    )

if options.skim: 
  outCommand = ['keep *', 'drop *_evtcleaner_*_*', 'drop *_TriggerResults_*_*']#remove unwanted new branches OS2LAna
  process.out = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p')
      ),
    fileName = cms.untracked.string('CondOutputSkim'),
    outputCommands = cms.untracked.vstring(outCommand )
    )

  process.outpath = cms.EndPath(process.out)


#process.schedule = cms.Schedule(process.p)

open('dump.py','w').write(process.dumpPython()) 

