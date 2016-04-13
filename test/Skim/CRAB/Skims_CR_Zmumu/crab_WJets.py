from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'WJetsamcnlo_25ns_CR_Zmumu'
config.General.workArea = 'B2GEDMNtuplesSkim_CR_Zmumu_20Nov/'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'skim_cfg.py' 
config.JobType.pyCfgParams = ['isData=False', 'skimType=CR_Zmumu']

config.section_("Data")
config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/vorobiev-B2GAnaFW_Run2Spring15_25ns_v74x_V61-d5e3976e154b1b4043d031aaa9c2809b/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.ignoreLocality = False
config.Data.publication = True
config.Data.outLFNDirBase = '/store/group/lpcbprime/noreplica/skhalil/B2GEDMNtuplesSkim_CR_Zmumu_20Nov/'
# This string is used to construct the output dataset name

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'

config.section_('User')
