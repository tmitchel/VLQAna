import os, sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

#######################
## configure options ##
#######################
options = VarParsing('analysis')

options.register('isData', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is data?"
    )
options.register('FileNames', 'FileNames_DY',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of list of input files"
    )
options.register('htTrig', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "use high HT trigger"
    )
options.register('muTrig', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "use Mu trigger"
    )
options.register('isPhoton', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "is a photon"
    )

options.setDefault('maxEvents', -1)
options.parseArguments()

dataPath = '../data/'

################################################
## create the process and feed in input files ##
################################################
process = cms.Process('trigAna')
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        #FileNames[options.FileNames]
'root://cmseos.fnal.gov//store/user/dmendis/ntuples2016/reminiaod/batch2/Electrons/SingleEG1/SingleEG1_skim_984.root', 
'root://cmseos.fnal.gov//store/user/dmendis/ntuples2016/reminiaod/batch2/Electrons/SingleEG1/SingleEG1_skim_1.root', 
'root://cmseos.fnal.gov//store/user/dmendis/ntuples2016/reminiaod/batch2/Electrons/SingleEG1/SingleEG1_skim_235.root', 
'root://cmseos.fnal.gov//store/user/dmendis/ntuples2016/reminiaod/batch2/Electrons/SingleEG1/SingleEG1_skim_835.root', 
'root://cmseos.fnal.gov//store/user/dmendis/ntuples2016/reminiaod/batch2/Electrons/SingleEG1/SingleEG1_skim_395.root', 
'root://cmseos.fnal.gov//store/user/dmendis/ntuples2016/reminiaod/batch2/Electrons/SingleEG1/SingleEG1_skim_162.root', 
'root://cmseos.fnal.gov//store/user/dmendis/ntuples2016/reminiaod/batch2/Electrons/SingleEG1/SingleEG1_skim_2.root', 

        )
    )
process.TFileService = cms.Service("TFileService",
        fileName = cms.string(
                'output.root'
                      )
            )

############################
## setup reporting module ##
############################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#################################
## setup generic event cleaner ##
## not sure if PU is necessary ##
#################################
process.load("Analysis.VLQAna.EventCleaner_cff")
process.evtcleaner.File_PUDistData = cms.string(os.path.join(dataPath,'RunII2016Rereco_25ns_PUXsec69000nb.root'))
process.evtcleaner.File_PUDistDataLow   = cms.string(os.path.join(dataPath,'RunII2016Rereco_25ns_PUXsec65550nb.root'))
process.evtcleaner.File_PUDistDataHigh  = cms.string(os.path.join(dataPath,'RunII2016Rereco_25ns_PUXsec72450nb.root'))
process.evtcleaner.File_PUDistMC        = cms.string(os.path.join(dataPath,'PUDistMC_Summer2016_25ns_Moriond17MC_PoissonOOTPU.root'))
process.evtcleaner.isData = options.isData 

##################################
## setup tag and probe triggers ##
##################################
if options.htTrig and options.muTrig:
  sys.exit("NOPE. Can only choose one schema at a time")

if options.htTrig:
  tag_hltpath = ['HLT_PFHT900_v']
elif options.muTrig:
  tag_hltpath = [
      'HLT_IsoMu24_v*',
      'HLT_IsoTkMu24_v*'
      ]
else:
  tag_hltpath = ['HLT_Ele32_eta2p1_WPTight_Gsf_v']

probe_hltpath = [
  'HLT_Ele115_CaloIdVT_GsfTrkIdT_v',
  'HLT_Photon175_v'
  ]

reject_hltpath = [
  'HLT_Ele115_CaloIdVT_GsfTrkIdT_v',
    ]

process.tag = process.evtcleaner.clone(
    hltPaths = cms.vstring(tag_hltpath)
    )
process.probe = process.evtcleaner.clone(
    hltPaths = cms.vstring(probe_hltpath)
    )
process.reject = process.evtcleaner.clone(
    hltPaths = cms.vstring(reject_hltpath)
    )

from Analysis.VLQAna.trigger_cfi import *
process.trigAna = trigAna.clone(
    isData = cms.bool(options.isData),
    isPhoton = cms.bool(options.isPhoton),
    htTrig = cms.bool(options.htTrig),
    muTrig = cms.bool(options.muTrig),
    )

if options.muTrig:
  process.trigAna.elselParams.elPtMin = cms.double(120.)
  process.trigAna.muselParams.muPtMin = cms.double(30.)
else:
  process.trigAna.elselParams.elPtMin = cms.double(35.)

###################
## event counter ##
###################
#from Analysis.EventCounter.eventcounter_cfi import eventCounter
#process.allEvents = eventCounter.clone(isData=options.isData)

###########################
## setup processing path ##
###########################
process.p = cms.Path(
 #   process.allEvents
    process.tag
    *cms.ignore(process.probe)
    *cms.ignore(process.reject)
    *process.trigAna
    )

open('dump.py', 'w').write(process.dumpPython())
