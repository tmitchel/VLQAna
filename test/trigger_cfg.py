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
            'root://cms-xrd-global.cern.ch//store/user/oiorio/samples/May/17May/B2GAnaFW_80X_V3p1/SingleElectron/Run2016B/SingleElectron/Run2016B-03Feb2017_ver2-v2_B2GAnaFW_80X_V3p1/170517_123348/0000/B2GEDMNtuple_1.root',
'root://cms-xrd-global.cern.ch//store/user/oiorio/samples/May/17May/B2GAnaFW_80X_V3p1/SingleElectron/Run2016B/SingleElectron/Run2016B-03Feb2017_ver2-v2_B2GAnaFW_80X_V3p1/170517_123348/0000/B2GEDMNtuple_2.root',
'root://cms-xrd-global.cern.ch//store/user/oiorio/samples/May/17May/B2GAnaFW_80X_V3p1/SingleElectron/Run2016B/SingleElectron/Run2016B-03Feb2017_ver2-v2_B2GAnaFW_80X_V3p1/170517_123348/0000/B2GEDMNtuple_3.root',
'root://cms-xrd-global.cern.ch//store/user/oiorio/samples/May/17May/B2GAnaFW_80X_V3p1/SingleElectron/Run2016B/SingleElectron/Run2016B-03Feb2017_ver2-v2_B2GAnaFW_80X_V3p1/170517_123348/0000/B2GEDMNtuple_4.root',
'root://cms-xrd-global.cern.ch//store/user/oiorio/samples/May/17May/B2GAnaFW_80X_V3p1/SingleElectron/Run2016B/SingleElectron/Run2016B-03Feb2017_ver2-v2_B2GAnaFW_80X_V3p1/170517_123348/0000/B2GEDMNtuple_5.root'
 
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
probe_hltpath = ['HLT_Ele32_eta2p1_WPTight_Gsf_v']
tag_hltpath = ['HLT_Ele115_CaloIdVT_GsfTrkIdT_v']

process.tag = process.evtcleaner.clone(
    hltPaths = cms.vstring(probe_hltpath)
    )
process.probe = process.evtcleaner.clone(
    hltPaths = cms.vstring(tag_hltpath)
    )

from Analysis.VLQAna.trigger_cfi import *
process.trigAna = trigAna.clone(
    isData = cms.bool(options.isData)
    )

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
    *process.probe
    *process.trigAna
    )

open('dump.py', 'w').write(process.dumpPython())
