import os, sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('isData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is data?"
    )
options.register('isPhoton', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "trigger for photon datastream"
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
options.register('signalType', '',
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
options.register('btageffmap', "ttjets_bTagEff.root",#until new SFs arrive
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
options.register('massReco', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Optimize mass reconstruction"
    )
options.register('syst', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do systematics"
    )
options.register('storeLHEWts', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    'Store LHE weights'
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
  if options.isPhoton:
    hltpaths = [
        'HLT_Photon175_v'
        ]
    rejectpaths = [
        'HLT_Ele115_CaloIdVT_GsfTrkIdT_v'
        ]
  else:
    hltpaths = [
        #"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        #"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
        "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
        "HLT_Photon175_v",
        ]
    rejectpaths = []
  
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

process = cms.Process("OS2LAnaOnSkim")

process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring(
#    FileNames[options.FileNames]
    'path_to_sample'
    ) 
  )

#if options.isData:
#  import FWCore.PythonUtilities.LumiList as LumiList
#  process.source.lumisToProcess = LumiList.LumiList(
#      filename = dataPath+"/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
#      ).getVLuminosityBlockRange()

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
      'CondOutputEvt'
      #      options.outFileName
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

if options.isPhoton:
  process.evtrejecter = process.evtcleaner.clone(
      hltPaths = cms.vstring(rejectpaths)
      )

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
    isPhoton = cms.bool(options.isPhoton),
    maketree = cms.bool(options.maketree), 
    fnamebtagSF = cms.string(os.path.join(dataPath,'CSVv2_Moriond17_B_H.csv')),
    fnameSJbtagSF = cms.string(os.path.join(dataPath,'subjet_CSVv2_Moriond17_B_H.csv')),
    File_DYNLOCorr = cms.string(os.path.join(dataPath,'scalefactors_v4.root')),
    syst = cms.bool(False),
    vv = cms.bool(False),
    DYSFSyst = cms.bool(True)
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
if options.zdecaymode == "zelel": process.ana.lepTrigSFsParams.eltrigeffmap = cms.string(os.path.join(dataPath,"ElectronTriggerSF_v2.root")) 
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
  process.ana.HTMin = cms.double(200.)
else: 
  process.ana.STMin = cms.double(1000.)
  process.ana.HTMin = cms.double(200.)

process.ana.pdfID_offset = cms.int32(9)
process.ana.scale_offset = cms.int32(0)

if options.skim: 
  process.ana.jetAK4selParams.jetPtMin = cms.double(20) 
  process.ana.jetAK8selParams.jetPtMin = cms.double(170) 
  process.ana.jetWTaggedselParams.jetPtMin = cms.double(170) 
  process.ana.jetHTaggedselParams.jetPtMin = cms.double(170) 
  process.ana.STMin = cms.double(0.)
else: 
  if options.syst: 
    process.anabcPlus = process.ana.clone(
        btagsf_bcUp = cms.bool(True),
        sbtagsf_bcUp = cms.bool(True),
        syst = cms.bool(True)
        )
    process.anabcMinus = process.ana.clone(
        btagsf_bcDown = cms.bool(True),
        sbtagsf_bcDown = cms.bool(True),
        syst = cms.bool(True)
        )
    process.analightPlus = process.ana.clone(
        btagsf_lUp = cms.bool(True),
        sbtagsf_lUp = cms.bool(True),
        syst = cms.bool(True)
        )
    process.analightMinus = process.ana.clone(
        btagsf_lDown = cms.bool(True),
        sbtagsf_lDown = cms.bool(True),
        syst = cms.bool(True)
        )
    
    process.anaJecPlus = process.ana.clone()
    process.anaJecPlus.jetAK4selParams.jecShift = cms.double(1.)
    process.anaJecPlus.jetAK4BTaggedselParams.jecShift = cms.double(1.)
    process.anaJecPlus.jetAK8selParams.jecShift = cms.double(1.)
    process.anaJecPlus.jetHTaggedselParams.jecShift = cms.double(1.)
    process.anaJecPlus.jetWTaggedselParams.jecShift = cms.double(1.)
    process.anaJecPlus.jetTopTaggedselParams.jecShift = cms.double(1.)
    process.anaJecPlus.syst = cms.bool(True)
    
    process.anaJecMinus = process.ana.clone()
    process.anaJecMinus.jetAK4selParams.jecShift = cms.double(-1.)
    process.anaJecMinus.jetAK4BTaggedselParams.jecShift = cms.double(-1.)
    process.anaJecMinus.jetAK8selParams.jecShift = cms.double(-1.)
    process.anaJecMinus.jetHTaggedselParams.jecShift = cms.double(-1.)
    process.anaJecMinus.jetWTaggedselParams.jecShift = cms.double(-1.)
    process.anaJecMinus.jetTopTaggedselParams.jecShift = cms.double(-1.)
    process.anaJecMinus.syst = cms.bool(True)

    process.anaJerPlus = process.ana.clone()
    process.anaJerPlus.jetAK4selParams.jerShift = cms.int32(2)
    process.anaJerPlus.jetAK4BTaggedselParams.jerShift = cms.int32(2)
    process.anaJerPlus.jetAK8selParams.jerShift = cms.int32(2)
    process.anaJerPlus.jetHTaggedselParams.jerShift = cms.int32(2)
    process.anaJerPlus.jetWTaggedselParams.jerShift = cms.int32(2)
    process.anaJerPlus.jetTopTaggedselParams.jerShift = cms.int32(2)
    process.anaJerPlus.syst = cms.bool(True)

    process.anaJerMinus = process.ana.clone()
    process.anaJerMinus.jetAK4selParams.jerShift = cms.int32(0)
    process.anaJerMinus.jetAK4BTaggedselParams.jerShift = cms.int32(0)
    process.anaJerMinus.jetAK8selParams.jerShift = cms.int32(0)
    process.anaJerMinus.jetHTaggedselParams.jerShift = cms.int32(0)
    process.anaJerMinus.jetWTaggedselParams.jerShift = cms.int32(0)
    process.anaJerMinus.jetTopTaggedselParams.jerShift = cms.int32(0)
    process.anaJerMinus.syst = cms.bool(True)
    
    process.anaJmrPlus = process.ana.clone()
    process.anaJmrPlus.jetHTaggedselParams.jmrShift = cms.int32(2)
    process.anaJmrPlus.jetWTaggedselParams.jmrShift = cms.int32(2)
    process.anaJmrPlus.syst = cms.bool(True)

    process.anaJmrMinus = process.ana.clone()
    process.anaJmrMinus.jetHTaggedselParams.jmrShift = cms.int32(-1)
    process.anaJmrMinus.jetWTaggedselParams.jmrShift = cms.int32(-1)
    process.anaJmrMinus.syst = cms.bool(True)
 
    process.anaPileupPlus = process.ana.clone(
        PileupUp = cms.bool(True),
        syst = cms.bool(True),
        )
    process.anaPileupMinus = process.ana.clone(
        PileupDown = cms.bool(True),
        syst = cms.bool(True),
        )
    process.anaDYSFMinus = process.ana.clone(
        DYSFSyst = cms.bool(False),
        syst = cms.bool(True),
        )
    process.anatau21Plus = process.ana.clone(
        tau21Shift = cms.int32(1),
        syst = cms.bool(True),
        )
    process.anatau21Minus = process.ana.clone(
        tau21Shift = cms.int32(-1),
        syst = cms.bool(True),
        )
    process.anatau32Plus = process.ana.clone(
        tau32Shift = cms.int32(1),
        syst = cms.bool(True),
        )
    process.anatau32Minus = process.ana.clone(
        tau32Shift = cms.int32(-1),
        syst = cms.bool(True),
        )
    process.anaElsfPlus = process.ana.clone(
        elSyst = cms.int32(1),
        syst = cms.bool(True)
        )
    process.anaElsfMinus = process.ana.clone(
        elSyst = cms.int32(-1),
        syst = cms.bool(True)
        )
    process.anaHiggsUp = process.ana.clone(
        higgsShift = cms.int32(1),
        syst = cms.bool(True),
        )
    process.anaHiggsDown = process.ana.clone(
        higgsShift = cms.int32(-1),
        syst = cms.bool(True),
        )
    process.anaTopUp = process.ana.clone(
        topShift = cms.int32(1),
        syst = cms.bool(True),
        )
    process.anaTopDown = process.ana.clone(
        topShift = cms.int32(-1),
        syst = cms.bool(True),
        )
    process.anaWZUp = process.ana.clone(
        WZShift = cms.int32(1),
        syst = cms.bool(True),
        )
    process.anaWZDown = process.ana.clone(
        WZShift = cms.int32(-1),
        syst = cms.bool(True),
        )

## Event counters
from Analysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone(isData=options.isData)
process.cleanedEvents = eventCounter.clone(isData=options.isData)
process.finalEvents = eventCounter.clone(isData=options.isData)

process.load("Analysis.VLQAna.VLQCandProducer_cff")

if options.syst and not options.skim and not options.isData:
  if options.isPhoton:
    process.p = cms.Path(
      process.allEvents
      *process.evtcleaner
      *process.evtrejecter
      *process.cleanedEvents
      *cms.ignore(
        process.ana
        *process.finalEvents
        )
      *cms.ignore(process.anabcPlus)
      *cms.ignore(process.anabcMinus)
      *cms.ignore(process.analightPlus)
      *cms.ignore(process.analightMinus)
      *cms.ignore(process.anaJecPlus)
      *cms.ignore(process.anaJecMinus)
      *cms.ignore(process.anaJerPlus)
      *cms.ignore(process.anaJerMinus)
      *cms.ignore(process.anaPileupPlus)
      *cms.ignore(process.anaPileupMinus)
      *cms.ignore(process.anatau21Plus)
      *cms.ignore(process.anatau21Minus)
      *cms.ignore(process.anatau32Plus)
      *cms.ignore(process.anatau32Minus)
      *cms.ignore(process.anaJmrPlus)
      *cms.ignore(process.anaJmrMinus)
      *cms.ignore(process.anaDYSFMinus)
      *cms.ignore(process.anaElsfPlus)
      *cms.ignore(process.anaElsfMinus)
      *cms.ignore(process.anaHiggsUp)
      *cms.ignore(process.anaHiggsDown)
      *cms.ignore(process.anaTopUp)
      *cms.ignore(process.anaTopDown)
      *cms.ignore(process.anaWZUp)
      *cms.ignore(process.anaWZDown)
      )
  else:
    process.p = cms.Path(
      process.allEvents
      *process.evtcleaner
      *process.cleanedEvents
      *cms.ignore(
        process.ana
        *process.finalEvents
        )
      *cms.ignore(process.anabcPlus)
      *cms.ignore(process.anabcMinus)
      *cms.ignore(process.analightPlus)
      *cms.ignore(process.analightMinus)
      *cms.ignore(process.anaJecPlus)
      *cms.ignore(process.anaJecMinus)
      *cms.ignore(process.anaJerPlus)
      *cms.ignore(process.anaJerMinus)
      *cms.ignore(process.anaPileupPlus)
      *cms.ignore(process.anaPileupMinus)
      *cms.ignore(process.anatau21Plus)
      *cms.ignore(process.anatau21Minus)
      *cms.ignore(process.anatau32Plus)
      *cms.ignore(process.anatau32Minus)
      *cms.ignore(process.anaJmrPlus)
      *cms.ignore(process.anaJmrMinus)
      *cms.ignore(process.anaDYSFMinus)
      *cms.ignore(process.anaElsfPlus)
      *cms.ignore(process.anaElsfMinus)
      *cms.ignore(process.anaHiggsUp)
      *cms.ignore(process.anaHiggsDown)
      *cms.ignore(process.anaTopUp)
      *cms.ignore(process.anaTopDown)
      *cms.ignore(process.anaWZUp)
      *cms.ignore(process.anaWZDown)
      )

elif options.massReco:
  process.load('Analysis.VLQAna.MassReco_cfi')
  process.massReco.ptMin = cms.double(150.)
  process.massReco.zdecaymode = cms.string(options.zdecaymode)
  process.massReco.signalType = cms.string(options.signalType)
  process.p = cms.Path(
    process.allEvents
    *process.evtcleaner
    *process.cleanedEvents
    *process.ana
    *process.massReco
    *process.finalEvents
    )
else:
  if options.isPhoton:
    process.p = cms.Path(
      process.allEvents
      *process.evtcleaner
      *process.evtrejecter
      *process.cleanedEvents
      *process.ana
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
