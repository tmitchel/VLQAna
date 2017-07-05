import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')

options.register('isData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is data?"
    )
options.register('outFileName', 'hh4bTrig',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('hltORAND', 'OR',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "OR or AND HLT paths?"
    )
options.register('cleanEvents', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Clean events using EventCleaner?"
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
options.register('doPUReweightingOfficial', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do pileup reweighting using official recipe"
    )
options.register('storeLHEWts', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Store LHE wts?"
    )
options.register('FileNames', 'FileNames',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of list of input files"
    )
options.register('printEvtNo', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print event numbers"
    )

options.setDefault('maxEvents', -1)
options.parseArguments()

hltpathsLoose = [
  #"HLT_Mu50_v",
  #"HLT_PFJet260_v",
  ]

hltpathsTest = [
            "HLT_PFJet260_v",
            #'HLT_PFHT800_v', 
            #'HLT_PFHT900_v', 
            #'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v', 
            #'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v', 
            #'HLT_AK8PFJet360_TrimMass30_v',
            #'HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v',
            #'AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v',
            ]

if options.isData: 
    options.jerShift = 0 
    options.doPUReweightingOfficial=False 
    options.storeLHEWts=False

process = cms.Process("Trig")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

from inputFiles_cfi import * 
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
      FileNames[options.FileNames]
    ) 
    )

process.TFileService = cms.Service("TFileService",
       fileName = cms.string(
         options.outFileName+".root"
         )
       )

process.load("Analysis.VLQAna.EventCleaner_cff") 
process.evtcleaner.File_PUDistData= cms.string('RunII2016Rereco_25ns_PUXsec69000nb.root')
process.evtcleaner.File_PUDistDataLow =  cms.string('RunII2016Rereco_25ns_PUXsec65550nb.root')
process.evtcleaner.File_PUDistDataHigh = cms.string('RunII2016Rereco_25ns_PUXsec72450nb.root')
process.evtcleaner.File_PUDistMC = cms.string('PUDistMC_Summer2016_25ns_Moriond17MC_PoissonOOTPU.root')
process.evtcleaner.hltORAND = cms.string (options.hltORAND)  
process.evtcleaner.hltPaths = cms.vstring (hltpathsLoose)  
process.evtcleaner.cleanEvents = cms.bool(options.cleanEvents)
process.evtcleaner.isData = options.isData 
process.evtcleaner.DoPUReweightingOfficial = options.doPUReweightingOfficial
process.evtcleaner.storeLHEWts = options.storeLHEWts

process.evtcleanerden = process.evtcleaner.clone(
    cleanEvents = cms.bool(True),
    )

process.evtcleanertest = process.evtcleaner.clone(
    hltPaths = cms.vstring (hltpathsTest),
    )

from Analysis.VLQAna.HH4b_cff import *
process.hh4b = hh4b.clone()

process.hh4btest = process.hh4b.clone(
    runno                      = cms.InputTag("evtcleanertest","runno"), 
    lumisec                    = cms.InputTag("evtcleanertest","lumisec"), 
    evtno                      = cms.InputTag("evtcleanertest","evtno"), 
    isData                     = cms.InputTag("evtcleanertest","isData"), 
    hltdecision                = cms.InputTag("evtcleanertest","hltdecision"), 
    evttype                    = cms.InputTag("evtcleanertest","evttype"),
    evtwtGen                   = cms.InputTag("evtcleanertest","evtwtGen"),
    evtwtPV                    = cms.InputTag("evtcleanertest","evtwtPV"),
    evtwtPVLow                 = cms.InputTag("evtcleanertest","evtwtPVLow"),
    evtwtPVHigh                = cms.InputTag("evtcleanertest","evtwtPVHigh"),
    npv                        = cms.InputTag("evtcleanertest","npv"),
    npuTrue                    = cms.InputTag("evtcleanertest","npuTrue"),
    htHat                      = cms.InputTag("evtcleanertest","htHat"),
    lhewtids                   = cms.InputTag("evtcleanertest","lhewtids"),
    lhewts                     = cms.InputTag("evtcleanertest","lhewts"), 
    )

from Analysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone(isData=options.isData)
process.cleanedEvents = eventCounter.clone(isData=options.isData)
process.finalEvents = eventCounter.clone(isData=options.isData)

process.p = cms.Path(
    process.allEvents
    *process.evtcleaner
    *cms.ignore(process.hh4b)
    *process.evtcleanerden
    *process.evtcleanertest
    *process.hh4btest
    * process.finalEvents
    )

open('dump.py','w').write(process.dumpPython())
