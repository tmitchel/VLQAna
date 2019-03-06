from WMCore.Configuration import Configuration
from glob import glob

config = Configuration()

config.section_("General")

config.General.requestName = DUMMY_NAME
config.General.workArea = DUMMY_WORKDIR
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'os2lana_cfg.py'
config.JobType.pyCfgParams = [DATA, MODE, 'filterSignal=False','applyBTagSFs=False', 'applyLeptonTrigSFs=False', 'applyDYNLOCorr=False','skim=True', 'maketree=False', 'doPUReweightingOfficial=False', 'applyLeptonIDSFs=False', 'massReco=False', 'syst=False', 'storeLHEWts=False', JECNAME, ISPHOTON]
#config.JobType.inputFiles = ['RunII2016Rereco_25ns_PUXsec69000.root','RunII2016Rereco_25ns_PUXsec65550.root','RunII2016Rereco_25ns_PUXsec72450.root','PUDistMC_Summer2016_25ns_Moriond17MC_PoissonOOTPU.root','CSVv2_Moriond17_B_H.csv','subjet_CSVv2_Moriond17_B_H.csv','scalefactors_v4.root','inputFiles_cfi.py','Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt','Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt','CSVv2_ichep.csv','Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt','Spring16_25nsV6_MC_Uncertainty_AK8PFchs.txt']
config.JobType.inputFiles = [ifile for ifile in glob('data/*')]

config.section_("Data")
config.Data.inputDataset = DUMMY_DATASET
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = DUMMY_NUMBER
config.Data.ignoreLocality = False
config.Data.publication = False
config.Data.outLFNDirBase = DUMMY_OUTPUT_PATH

config.section_("Site")
config.Site.storageSite = DUMMY_SITE
config.section_('User')

