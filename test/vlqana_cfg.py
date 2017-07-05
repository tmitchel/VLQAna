import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('isData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is data?"
    )
options.register('outFileName', 'singleT',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('hltORAND', 'AND',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "OR or AND HLT paths?"
    )
options.register('cleanEvents', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Clean events using EventCleaner?"
    )
options.register('doPUReweightingOfficial', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do pileup reweighting using official recipe"
    )
options.register('doBTagSFUnc', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply b-tag SF uncertainties"
    )
options.register('jecShift', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "JEC shift"
    )
options.register('jerShift', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "JER shift"
    )
options.register('topTagtau32', 0.54,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Top-tagging tau32 cut"
    )
options.register('topTagBDisc', 0.79,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Top-tagging b-discriminator cut"
    )
options.register('HTMin', 900,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum HT"
    )
options.register('storePreselEvts', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Store pre-selected events after pre-selection", 
    )
options.register('doPreselOnly', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Only run pre-selections"
    )
options.register('storeLHEWts', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Store LHE wts?"
    )
options.register('FileNames', 'test',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of list of input files"
    )

options.setDefault('maxEvents', -1000)
options.parseArguments()

hltpaths = ["HLT_PFHT800_v"]

if options.isData:
    options.doBTagSFUnc = False 
    options.jerShift = 0 
    options.doPUReweightingOfficial=False 
    options.storeLHEWts=False

HTMin=1100
if options.storePreselEvts:
  HTMin = options.HTMin

process = cms.Process("VLQAna")

from inputFiles_cfi import * 
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      FileNames[options.FileNames]
      )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
      options.outFileName+".root"
      )
    )

process.load("Analysis.VLQAna.EventCleaner_cff") 
process.evtcleaner.hltORAND = cms.string (options.hltORAND)  
process.evtcleaner.hltPaths = cms.vstring (hltpaths)  
process.evtcleaner.cleanEvents = cms.bool(options.cleanEvents)
process.evtcleaner.isData = options.isData 
process.evtcleaner.DoPUReweightingOfficial = options.doPUReweightingOfficial
process.evtcleaner.storeLHEWts = options.storeLHEWts

from Analysis.VLQAna.VLQAna_cfi import *

if options.isData == False: ### Careful, to be reset when B2GAnaFW_v80X_v2.4 MC are used
  for par in ['jetAK8selParams', 'jetHTaggedselParams', 'jetAntiHTaggedselParams', 'jetTopTaggedselParams']:
    setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jettau1Label'         ,cms.InputTag("jetsAK8CHS", "jetAK8CHStau1"))
    setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jettau2Label'         ,cms.InputTag("jetsAK8CHS", "jetAK8CHStau2"))
    setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jettau3Label'         ,cms.InputTag("jetsAK8CHS", "jetAK8CHStau3"))
    setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jetPrunedMassLabel'   ,cms.InputTag("jetsAK8CHS", "jetAK8CHSprunedMass"))
    setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jetTrimmedMassLabel'  ,cms.InputTag("jetsAK8CHS", "jetAK8CHStrimmedMass"))
    setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jetFilteredMassLabel' ,cms.InputTag("jetsAK8CHS", "jetAK8CHSfilteredMass"))
    setattr(getattr(getattr(ana, par), 'JetSubstrParams'), 'jetSoftDropMassLabel' ,cms.InputTag("jetsAK8CHS", "jetAK8CHSsoftDropMass"))

process.ana = ana.clone()
process.ana.doBTagSFUnc = options.doBTagSFUnc
process.ana.jetAK4selParams.jecShift = options.jecShift 
process.ana.jetAK4selParams.jerShift = options.jerShift 
process.ana.jetAK8selParams.jecShift = options.jecShift 
process.ana.jetAK8selParams.jerShift = options.jerShift 
process.ana.jetHTaggedselParams.jecShift = options.jecShift 
process.ana.jetHTaggedselParams.jerShift = options.jerShift 
process.ana.jetTopTaggedselParams.jecShift = options.jecShift 
process.ana.jetTopTaggedselParams.jerShift = options.jerShift 
process.ana.jetAntiHTaggedselParams.jecShift = options.jecShift 
process.ana.jetAntiHTaggedselParams.jerShift = options.jerShift 
#process.ana.jetTopTaggedselParams.jettau3Bytau2Max = options.topTagtau32
#process.ana.jetTopTaggedselParams.subjetHighestCSVMin = options.topTagBDisc
process.ana.storePreselEvts = options.storePreselEvts
process.ana.doPreselOnly = options.doPreselOnly
process.ana.HTMin = HTMin

process.anaCHS = process.ana.clone()

from Analysis.VLQAna.JetMaker_cfi import *
process.anaPuppi = process.ana.clone(
    jetAK8selParams = defaultAK8PuppiJetSelectionParameters,
    jetHTaggedselParams = defaultPuppiHJetSelectionParameters,
    jetAntiHtaggedSelParams = defaultPuppiHJetSelectionParameters.clone(
                                  subjetCSVMin = cms.double(-1000000),
                                  subjetCSVMax = defaultPuppiHJetSelectionParameters.subjetCSVMin,
                                  ),
    jetTopTaggedselParams = defaultPuppiTJetSelectionParameters.clone(),
    
    )

process.anaDoubleB = process.ana.clone(
    jetHTaggedselParams = defaultCHSHJetSelectionParameters.clone(
                                  subjetCSVMin = cms.double(-1000000),
                                  subjetCSVMax = cms.double(1000000),
                                  jetDoubleBMin = cms.double(0.8),
                                  ),
    jetAntiHTaggedSelParams = defaultCHSHJetSelectionParameters.clone(
                                  subjetCSVMin = cms.double(-1000000),
                                  subjetCSVMax = cms.double(1000000),
                                  jetDoubleBMax = cms.double(0.8),
                                  ),
    )
 
from Analysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone(isData=options.isData)
process.cleanedEvents = eventCounter.clone(isData=options.isData)
process.finalEvents = eventCounter.clone(isData=options.isData)
process.finalEventsPuppi = eventCounter.clone(isData=options.isData)
process.finalEventsDoubleB = eventCounter.clone(isData=options.isData)

process.p = cms.Path(
    process.allEvents
    *process.evtcleaner
    *process.cleanedEvents
    *process.anaCHS 
    *process.finalEvents
    )

process.pPuppi = cms.Path(
    process.allEvents
    *process.evtcleaner
    *process.anaPuppi
    *process.finalEventsPuppi
    )

process.pDoubleB = cms.Path(
    process.allEvents
    *process.evtcleaner
    *process.anaDoubleB
    *process.finalEventsDoubleB
    )

process.out = cms.OutputModule("PoolOutputModule",
        SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('evtcleaner')),
            )

process.schedule = cms.Schedule(process.p)

#process.outpath = cms.EndPath(process.out)

open('dump.py','w').write(process.dumpPython())
