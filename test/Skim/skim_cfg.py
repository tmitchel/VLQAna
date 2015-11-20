import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('outFileName', 'SkimmedB2GEdmNtuples.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('isData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is data?"
    )
options.register('skimType', '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Skim type: CR_Zelel, CR_Zmumu, SR_Zelel, SR_Zmumu" 
    )
options.setDefault('maxEvents', 1000)
options.parseArguments()

###
###Skim variables:
if options.skimType == "":
    sys.exit("!!!!Error: Enter 'skimType' . Options are: ''CR_Zelel, 'CR_Zmumu', 'SR_Zelel', 'SR_Zmumu'.\n")
global zptMin, nzelel, nzmumu, nak4jets, nak8jets 
if options.skimType == 'CR_Zelel':
  zptMin= 0.
  nzelel=1
  nzmumu=0
  nak4jets=1
  nak8jets=0
elif options.skimType == 'CR_Zmumu':
  zptMin= 0.
  nzelel=0
  nzmumu=1
  nak4jets=1
  nak8jets=0
elif options.skimType == 'SR_Zelel':
  zptMin= 200.
  nzelel=1
  nzmumu=0
  nak4jets=1
  nak8jets=1
elif options.skimType == 'SR_Zmumu':
  zptMin= 200.
  nzelel=0
  nzmumu=1
  nak4jets=1
  nak8jets=1
else:
  sys.exit("!!!!Error: Wrong 'skimType' entered. Options are: ''CR_Zelel, 'CR_Zmumu', 'SR_Zelel', 'SR_Zmumu'.\n") 

###

process = cms.Process("Skim")

from inputFiles_cfi import *

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
      #FileNames_TprimeBToTH_M1200
      'file:/afs/cern.ch/user/d/devdatta/eos/cms/store/group/phys_b2g/vorobiev/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_Run2Spring15_25ns_v74x_V61/151018_123637/0000/B2GEDMNtuple_1.root'
      )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.TFileService = cms.Service("TFileService",
       fileName = cms.string(options.outFileName)
       )

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outFileName),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p')
      ),
    dropMetaData = cms.untracked.string('DROPPED'),
    outputCommands = cms.untracked.vstring("keep *",)
    )

from Analysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone(isData=options.isData)
process.finalEvents = eventCounter.clone(isData=options.isData)

from Analysis.VLQAna.Skim_cff import *
process.skim = skim.clone(
    ZCandParams = skim.ZCandParams.clone(ptMin =  cms.double(zptMin)), 
    elselParams = skim.elselParams.clone(useVID = cms.bool(options.isData)), 
    jetAK4selParams            = skim.jetAK4selParams.clone(
      jetPtMin = cms.double(40), 
      jetAbsEtaMax = cms.double(2.4),
      ),
    jetAK8selParams            = skim.jetAK8selParams.clone(
      jetPtMin = cms.double(170), 
      jetAbsEtaMax = cms.double(2.4),
      ),
    nzelel = nzelel,
    nzmumu = nzmumu,
    nak4jets = nak4jets,
    nak8jets = nak8jets, 
    )
#process.skim.elselParams.useVID = cms.bool(options.isData)

process.p = cms.Path(
    process.allEvents
    *process.skim
    *process.finalEvents
    ) 

#process.outpath = cms.EndPath(process.out)
open('dump.py','w').write(process.dumpPython())
