import os, sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


options = VarParsing('analysis')
options.register('isData', ISDATA,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is data?"
    )
options.register('skim', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Skim events?"
    )
options.register('newJECPayloadNames', 'JECPAYLOAD',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Select one of EvtType_MC_tZtZ, EvtType_MC_tZtH, EvtType_MC_tZbW, EvtType_MC_tHtH, EvtType_MC_tHbW, EvtType_MC_bWbW, EvtType_MC_bZbZ, EvtType_MC_bZbH, EvtType_MC_bZtW, EvtType_MC_bHbH, EvtType_MC_bHtW, EvtType_MC_tWtW" 
    )
options.register('maketree', MAKETREE,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Skim events?"
    )
options.register('zdecaymode', 'ZDECAYMODE',
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
options.register('doPUReweightingOfficial', PUREWEIGHT,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do pileup reweighting using official recipe"
    )
options.register('filterSignal', FILTERSIGNAL,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Select only tZXX or bZXX modes"
    )
options.register('signalType', 'SIGNALTYPE',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Select one of EvtType_MC_tZtZ, EvtType_MC_tZtH, EvtType_MC_tZbW, EvtType_MC_tHtH, EvtType_MC_tHbW, EvtType_MC_bWbW, EvtType_MC_bZbZ, EvtType_MC_bZbH, EvtType_MC_bZtW, EvtType_MC_bHbH, EvtType_MC_bHtW, EvtType_MC_tWtW" 
    )
options.register('applyLeptonIDSFs', LEPTONIDSFS,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply lepton SFs to the MC"
    )
options.register('applyBTagSFs', BTAGSFS,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply b-tagging SFs to the MC"
    )
options.register('btageffmap', "BTAGEFFMAP",#until new SFs arrive
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "ROOT file with Th2D histos of b tag effs for b,c, and light flavoured jets"
    )
options.register('applyLeptonTrigSFs', LEPTONTRIGSFS,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply trigger SFs to the MC"
    )
options.register('applyDYNLOCorr', DYNLOCORR, ### Set to true only for DY process ### Only EWK NLO k-factor is applied
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply DY EWK k-factor to DY MC"
    )
options.register('FileNames', 'FileNames_DY',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of list of input files"
    )
options.register('massReco', MASSRECO,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Optimize mass reconstruction"
    )
options.register('controlReco', CONTROLRECO,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'Do control region reco'
    )
options.register('syst', SYST,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do systematics"
    )
options.register('storeLHEWts', LHEWTS,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Store LHE weights",
    )
options.register('optimizeReco', OPTIMIZERECO,
    VarParsing.multiplicity.singleton, 
    VarParsing.varType.bool,
    "do Gen optimization"
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

process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring(
    'path_to_sample'
    ) 
  )

#if options.isData:
#  import FWCore.PythonUtilities.LumiList as LumiList
#  process.source.lumisToProcess = LumiList.LumiList(
#      filename = "../data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
#      ).getVLuminosityBlockRange()

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
      'CondOutputEvt'
      #options.outFileName
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
    optimizeReco = cms.bool(options.optimizeReco),
    syst = cms.bool(False),
    vv = cms.bool(ISVV)
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
process.massReco.syst = cms.bool(False)
process.massReco.vv = cms.bool(ISVV)

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
        syst = cms.bool(True)
        )
    process.anabcDown = process.ana.clone(
        btagsf_bcDown = cms.bool(True),
        sbtagsf_bcDown = cms.bool(True),
        syst = cms.bool(True)
        )
    process.analightUp = process.ana.clone(
        btagsf_lUp = cms.bool(True),
        sbtagsf_lUp = cms.bool(True),
        syst = cms.bool(True)
        )
    process.analightDown = process.ana.clone(
        btagsf_lDown = cms.bool(True),
        sbtagsf_lDown = cms.bool(True),
        syst = cms.bool(True)
        )
    
    process.anaJecUp = process.ana.clone()
    process.anaJecUp.jetAK4selParams.jecShift = cms.double(1.)
    process.anaJecUp.jetAK4BTaggedselParams.jecShift = cms.double(1.)
    process.anaJecUp.jetAK8selParams.jecShift = cms.double(1.)
    process.anaJecUp.jetHTaggedselParams.jecShift = cms.double(1.)
    process.anaJecUp.jetWTaggedselParams.jecShift = cms.double(1.)
    process.anaJecUp.jetTopTaggedselParams.jecShift = cms.double(1.)
    process.anaJecUp.syst = cms.bool(True)

    process.anaJecDown = process.ana.clone()
    process.anaJecDown.jetAK4selParams.jecShift = cms.double(-1.)
    process.anaJecDown.jetAK4BTaggedselParams.jecShift = cms.double(-1.)
    process.anaJecDown.jetAK8selParams.jecShift = cms.double(-1.)
    process.anaJecDown.jetHTaggedselParams.jecShift = cms.double(-1.)
    process.anaJecDown.jetWTaggedselParams.jecShift = cms.double(-1.)
    process.anaJecDown.jetTopTaggedselParams.jecShift = cms.double(-1.)
    process.anaJecDown.syst = cms.bool(True)

    process.anaJerUp = process.ana.clone()
    process.anaJerUp.jetAK4selParams.jerShift = cms.int32(2)
    process.anaJerUp.jetAK4BTaggedselParams.jerShift = cms.int32(2)
    process.anaJerUp.jetAK8selParams.jerShift = cms.int32(2)
    process.anaJerUp.jetHTaggedselParams.jerShift = cms.int32(2)
    process.anaJerUp.jetWTaggedselParams.jerShift = cms.int32(2)
    process.anaJerUp.jetTopTaggedselParams.jerShift = cms.int32(2)
    process.anaJerUp.syst = cms.bool(True)

    process.anaJerDown = process.ana.clone()
    process.anaJerDown.jetAK4selParams.jerShift = cms.int32(-1)
    process.anaJerDown.jetAK4BTaggedselParams.jerShift = cms.int32(-1)
    process.anaJerDown.jetAK8selParams.jerShift = cms.int32(-1)
    process.anaJerDown.jetHTaggedselParams.jerShift = cms.int32(-1)
    process.anaJerDown.jetWTaggedselParams.jerShift = cms.int32(-1)
    process.anaJerDown.jetTopTaggedselParams.jerShift = cms.int32(-1)
    process.anaJerDown.syst = cms.bool(True)
    
    process.anaJmrUp = process.ana.clone()
    process.anaJmrUp.jetHTaggedselParams.jmrShift = cms.int32(2)
    process.anaJmrUp.jetWTaggedselParams.jmrShift = cms.int32(2)
    process.anaJmrUp.syst = cms.bool(True)

    process.anaJmrDown = process.ana.clone()
    process.anaJmrDown.jetHTaggedselParams.jmrShift = cms.int32(-1)
    process.anaJmrDown.jetWTaggedselParams.jmrShift = cms.int32(-1)
    process.anaJmrDown.syst = cms.bool(True)

    process.anaPileupUp = process.ana.clone(
        PileupUp = cms.bool(True),
        syst = cms.bool(True)
        )
    process.anaPileupDown = process.ana.clone(
        PileupDown = cms.bool(True),
        syst = cms.bool(True)
        )
    process.anaDYSF = process.ana.clone(
        applyDYNLOCorr = cms.bool(True),
        syst = cms.bool(True)
        )

    process.anatauUp = process.ana.clone(
        tauShift = cms.int32(1),
         syst = cms.bool(True)
        )

    process.anatauDown = process.ana.clone(
        tauShift = cms.int32(-1),
         syst = cms.bool(True)
        )

if options.massReco and options.syst:
  process.recobc__plus = process.massReco.clone(
      bTagwt = cms.InputTag("ana", "btagsfbcUp"),
      sjbTagwt = cms.InputTag("ana", "sjbtagsfbcUp"),
      syst = cms.bool(True)
      )
  process.recobc__minus = process.massReco.clone(
      bTagwt = cms.InputTag("ana", "btagsfbcDown"),
      sjbTagwt = cms.InputTag("ana", "sjbtagsfbcDown"),
      syst = cms.bool(True)
      )
  process.recolight__plus = process.massReco.clone(
      bTagwt = cms.InputTag("ana", "btagsflUp"),
      sjbTagwt = cms.InputTag("ana", "sjbtagsflUp"),
      syst = cms.bool(True)
      )
  process.recolight__minus = process.massReco.clone(
      bTagwt = cms.InputTag("ana", "btagsflDown"),
      sjbTagwt = cms.InputTag("ana", "sjbtagsflDown"),
      syst = cms.bool(True)
      )
  process.recoJec__plus = process.massReco.clone(
      ak8jets = cms.InputTag("anaJecUp", "ak8jets"),
      jets  = cms.InputTag("anaJecUp", "jets"),
      bjets = cms.InputTag("anaJecUp", "bjets"),
      zjets = cms.InputTag("anaJecUp", "wjets"),
      hjets = cms.InputTag("anaJecUp", "hjets"),
      Prewt = cms.InputTag("anaJecUp", "PreWeight"),
      syst = cms.bool(True)
      )
  process.recoJec__minus = process.massReco.clone(
      ak8jets = cms.InputTag("anaJecDown","ak8jets"),
      jets  = cms.InputTag("anaJecDown","jets"),
      bjets = cms.InputTag("anaJecDown", "bjets"),
      zjets = cms.InputTag("anaJecDown", "wjets"),
      hjets = cms.InputTag("anaJecDown", "hjets"),
      Prewt = cms.InputTag("anaJecDown", "PreWeight"),
      syst = cms.bool(True)
      )
  process.recoJer__plus = process.massReco.clone(
      ak8jets = cms.InputTag("anaJerUp","ak8jets"),
      jets  = cms.InputTag("anaJerUp","jets"),
      bjets = cms.InputTag("anaJerUp", "bjets"),
      zjets = cms.InputTag("anaJerUp", "wjets"),
      hjets = cms.InputTag("anaJerUp", "hjets"),
      Prewt = cms.InputTag("anaJerUp", "PreWeight"),
      syst = cms.bool(True)
      )
  process.recoJer__minus = process.massReco.clone(
      ak8jets = cms.InputTag("anaJerDown","ak8jets"),
      jets  = cms.InputTag("anaJerDown","jets"),
      bjets = cms.InputTag("anaJerDown", "bjets"),
      zjets = cms.InputTag("anaJerDown", "wjets"),
      hjets = cms.InputTag("anaJerDown", "hjets"),
      Prewt = cms.InputTag("anaJerDown", "PreWeight"),
      syst = cms.bool(True)
      )

  process.recoJmr__plus = process.massReco.clone(
      ak8jets = cms.InputTag("anaJmrUp","ak8jets"),
      jets  = cms.InputTag("anaJmrUp","jets"),
      bjets = cms.InputTag("anaJmrUp", "bjets"),
      zjets = cms.InputTag("anaJmrUp", "wjets"),
      hjets = cms.InputTag("anaJmrUp", "hjets"),
      Prewt = cms.InputTag("anaJmrUp", "PreWeight"),
      syst = cms.bool(True)
      )
  process.recoJmr__minus = process.massReco.clone(
      ak8jets = cms.InputTag("anaJmrDown","ak8jets"),
      jets  = cms.InputTag("anaJmrDown","jets"),
      bjets = cms.InputTag("anaJmrDown", "bjets"),
      zjets = cms.InputTag("anaJmrDown", "wjets"),
      hjets = cms.InputTag("anaJmrDown", "hjets"),
      Prewt = cms.InputTag("anaJmrDown", "PreWeight"),
      syst = cms.bool(True)
      )
  
  process.recoPileup__plus    = process.massReco.clone(Prewt = cms.InputTag("anaPileupUp", "PreWeight"))
  process.recoPileup__plus.syst = cms.bool(True)
  process.recoPileup__minus   = process.massReco.clone(Prewt = cms.InputTag("anaPileupDown", "PreWeight"))
  process.recoPileup__minus.syst = cms.bool(True)
  process.recotau__plus       = process.massReco.clone(Prewt = cms.InputTag("anatauUp", "PreWeight"))
  process.recotau__plus.syst = cms.bool(True)
  process.recotau__minus      = process.massReco.clone(Prewt = cms.InputTag("anatauDown", "PreWeight"))
  process.recotau__minus.syst = cms.bool(True)
  process.recoDYSF__minus       = process.massReco.clone(Prewt = cms.InputTag("anaDYSF", "PreWeight"))
  process.recoDYSF__minus.syst = cms.bool(True)


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
    *cms.ignore(process.anaDYSF)
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
      *cms.ignore(process.anaDYSF)
      *cms.ignore(process.anatauUp)
      *cms.ignore(process.anatauDown)
      *cms.ignore(process.anaJmrUp)
      *cms.ignore(process.anaJmrDown)

      *cms.ignore(process.recobc__plus)
      *cms.ignore(process.recobc__minus)
      *cms.ignore(process.recolight__plus)
      *cms.ignore(process.recolight__minus)
      *cms.ignore(process.recoJec__plus)
      *cms.ignore(process.recoJec__minus)
      *cms.ignore(process.recoJer__plus)
      *cms.ignore(process.recoJer__minus)
      *cms.ignore(process.recoPileup__plus)
      *cms.ignore(process.recoPileup__minus)
      *cms.ignore(process.recotau__plus)
      *cms.ignore(process.recotau__minus) 
      *cms.ignore(process.recoJmr__plus)
      *cms.ignore(process.recoJmr__minus)
      *cms.ignore(process.recoDYSF__minus)

      *process.massReco
      *process.finalEvents

      )
  else:
    process.p = cms.Path(
      process.allEvents
      *process.evtcleaner
      *process.cleanedEvents
      *cms.ignore(process.ana)
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

