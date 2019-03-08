import os
import sys
import math
import subprocess
from glob import glob
from shutil import copyfile
from optparse import OptionParser

parser = OptionParser("analysis")
parser.add_option("-e", "--makeEOS", dest = "makeEOS",
                  action = "store_true", default = False,
                  help = "make eos directories as well"
                  )
parser.add_option("-n", "--dirName", dest = "directoryName",
                  default = "newDir",
                  help = "name of new condor directory"
                  )
parser.add_option("-p", "--eosPath", dest = "eosPath",
                  default = "",
                  help = "path to new directory"
                  )
parser.add_option("-t", "--makeTree", dest = "makeTree",
                  action = "store_true", default = False,
                  help = "produce trees"
                  )
parser.add_option("-s", "--syst", dest = "syst",
                  action = "store_true", default = False,
                  help = "run systematics"
                  )
parser.add_option('-m', '--transfer', dest = 'transfer',
                  action = 'store', default = 'transfer_review',
                  help = 'eos storage for cfgs'
                  )

(options, args) = parser.parse_args()

sampleDict = {
              'ttbar_pythia' : [5, 984, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/161222_110143/"],
              'ttbar_pythia_ext' : [5, 1192, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_backup_v2p4/170304_173805/"],
              'ttZ' : [5, 216, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/ttZJets_13TeV_madgraphMLM/RunIISummer16MiniAODv2/ttZJets_13TeV_madgraphMLM/RunIISummer16MiniAODv2-80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/170126_232339/"],
              'tZq' : [5, 302, '/store/group/phys_b2g/B2GAnaFW_80X_V2p4/tZq_ll_4f_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_25ns_80X_v2p4/180508_112043/'],
              'ST_tWll' : [1, 1, '/store/group/phys_b2g/B2GAnaFW_80X_V2p4/ST_tWll_5f_LO_13TeV-MadGraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_25ns_80X_v2p4/180508_112943/'],
              'ttW' : [5, 203, '/store/group/phys_b2g/B2GAnaFW_80X_V2p4/ttWJets_13TeV_madgraphMLM/RunIISummer16MiniAODv2/ttWJets_13TeV_madgraphMLM/RunIISummer16MiniAODv2-80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/170126_232439/'],

              'dy_pt100to250' : [5, 34, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170121_174920/"],
              'dy_pt100to250_ext' : [5, 34, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_ext1_v2p4/170122_184903/"],
              'dy_pt100to250_ext5': [5, 1056, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_ext5_v2p4/170304_171916/"], 

              'dy_pt250to400' : [2, 11, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170121_175352/"],
              'dy_pt250to400_ext' : [2, 9, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_ext1_v2p4/170122_185048/"],
              'dy_pt250to400_ext5': [5, 379, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_ext5_v2p4/170304_172140/"], 

              'dy_pt400to650' : [2, 18, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170121_175516/"],
              'dy_pt400to650_ext' : [2, 8, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_ext1_v2p4/170122_185328/"],

              'dy_pt650toInf' : [2, 15, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170121_175719/"],
              'dy_pt650toInf_ext' : [2, 8, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_ext1_v2p4/170122_185543/"],

              'ZZ' : [2, 12, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/ZZ_TuneCUETP8M1_13TeV-pythia8/B2GAnaFW_Spring16MiniAODv2_Moriond17_v80x_v2p4/170122_181713/"],
              'ZZ_ext' : [5, 66, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/ZZ_TuneCUETP8M1_13TeV-pythia8/B2GAnaFW_Spring16MiniAODv2_Moriond17_v80x_ext1_v2p4/170122_181521/"],
              'WW' : [10, 116, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/WW_TuneCUETP8M1_13TeV-pythia8/B2GAnaFW_Spring16MiniAODv2_Moriond17_v80x_ext1_v2p4/170122_180941/"],
              'WZ' : [10, 128, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/WZ_TuneCUETP8M1_13TeV-pythia8/B2GAnaFW_Spring16MiniAODv2_Moriond17_v80x_ext1_v2p4/170122_181324/"],

              'bprime800_bZbZ' : [2, 26, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192852/"],

              'bprime900_bZbZ' : [2, 17, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192939/"],

              'bprime1000_bZbZ' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193050/"],

              'bprime1100_bZbZ' : [2, 28, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193140/"],

              'bprime1200_bZbZ' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193237/"],

              'bprime1300_bZbZ' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193347/"],

              'bprime1400_bZbZ' : [2, 27, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193455/"],

              'bprime1500_bZbZ' : [2, 18, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193653/"],

              'bprime1700_bZbZ' : [2, 23, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193854/"],

              'bprime1800_bZbZ' : [2, 20, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_194007/"],

              'bprime800_bZbH' : [2, 26, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192852/"],

              'bprime900_bZbH' : [2, 17, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192939/"],

              'bprime1000_bZbH' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193050/"],

              'bprime1100_bZbH' : [2, 28, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193140/"],

              'bprime1200_bZbH' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193237/"],

              'bprime1300_bZbH' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193347/"],

              'bprime1400_bZbH' : [2, 27, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193455/"],

              'bprime1500_bZbH' : [2, 18, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193653/"],

              'bprime1700_bZbH' : [2, 23, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193854/"],

              'bprime1800_bZbH' : [2, 20, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_194007/"],

              'bprime800_bHbH' : [2, 26, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192852/"],

              'bprime900_bHbH' : [2, 17, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192939/"],

              'bprime1000_bHbH' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193050/"],

              'bprime1100_bHbH' : [2, 28, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193140/"],

              'bprime1200_bHbH' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193237/"],

              'bprime1300_bHbH' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193347/"],

              'bprime1400_bHbH' : [2, 27, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193455/"],

              'bprime1500_bHbH' : [2, 18, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193653/"],

              'bprime1700_bHbH' : [2, 23, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193854/"],

              'bprime1800_bHbH' : [2, 20, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_194007/"],

              'bprime800_bHtW' : [2, 26, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192852/"],

              'bprime900_bHtW' : [2, 17, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192939/"],

              'bprime1000_bHtW' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193050/"],

              'bprime1100_bHtW' : [2, 28, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193140/"],

              'bprime1200_bHtW' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193237/"],

              'bprime1300_bHtW' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193347/"],

              'bprime1400_bHtW' : [2, 27, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193455/"],

              'bprime1500_bHtW' : [2, 18, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193653/"],

              'bprime1700_bHtW' : [2, 23, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193854/"],

              'bprime1800_bHtW' : [2, 20, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_194007/"],

              'bprime800_bZtW' : [2, 26, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192852/"],

              'bprime900_bZtW' : [2, 17, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192939/"],

              'bprime1000_bZtW' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193050/"],

              'bprime1100_bZtW' : [2, 28, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193140/"],

              'bprime1200_bZtW' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193237/"],

              'bprime1300_bZtW' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193347/"],

              'bprime1400_bZtW' : [2, 27, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193455/"],

              'bprime1500_bZtW' : [2, 18, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193653/"],

              'bprime1700_bZtW' : [2, 23, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193854/"],

              'bprime1800_bZtW' : [2, 20, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_194007/"],

              'bprime800_tWtW' : [2, 26, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192852/"],

              'bprime900_tWtW' : [2, 17, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_192939/"],

              'bprime1000_tWtW' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193050/"],

              'bprime1100_tWtW' : [2, 28, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193140/"],

              'bprime1200_tWtW' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193237/"],

              'bprime1300_tWtW' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193347/"],

              'bprime1400_tWtW' : [2, 27, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193455/"],

              'bprime1500_tWtW' : [2, 18, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193653/"],

              'bprime1700_tWtW' : [2, 23, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_193854/"],

              'bprime1800_tWtW' : [2, 20, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/BprimeBprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170122_194007/"],

#              'tprime800_tZtZ' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195633/"],
#
#              'tprime900_tZtZ' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195742/"],
#
#              'tprime1000_tZtZ' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195900/"],
#
#              'tprime1100_tZtZ' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_182558/"],
#
#              'tprime1200_tZtZ' : [2, 36, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200051/"],
#
#              'tprime1300_tZtZ' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200229/"],
#
#              'tprime1400_tZtZ' : [2, 33, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200324/"],
#
#              'tprime1500_tZtZ' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200414/"],
#
#              'tprime1600_tZtZ' : [2, 35, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200506/"],
#
#              'tprime1700_tZtZ' : [2, 29, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200554/"],
#
#              'tprime1800_tZtZ' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200646/"],
#
#              'tprime800_tZtH' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195633/"],
#
#              'tprime900_tZtH' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195742/"],
#
#              'tprime1000_tZtH' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195900/"],
#
#              'tprime1100_tZtH' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_182558/"],
#
#              'tprime1200_tZtH' : [2, 36, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200051/"],
#
#              'tprime1300_tZtH' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200229/"],
#
#              'tprime1400_tZtH' : [2, 33, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200324/"],
#
#              'tprime1500_tZtH' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200414/"],
#
#              'tprime1600_tZtH' : [2, 35, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200506/"],
#
#              'tprime1700_tZtH' : [2, 29, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200554/"],
#
#              'tprime1800_tZtH' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200646/"],
#
#              'tprime800_tZbW' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195633/"],
#
#              'tprime900_tZbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195742/"],
#
#              'tprime1000_tZbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195900/"],
#
#              'tprime1100_tZbW' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_182558/"],
#
#              'tprime1200_tZbW' : [2, 36, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200051/"],
#
#              'tprime1300_tZbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200229/"],
#
#              'tprime1400_tZbW' : [2, 33, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200324/"],
#
#              'tprime1500_tZbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200414/"],
#
#              'tprime1600_tZbW' : [2, 35, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200506/"],
#
#              'tprime1700_tZbW' : [2, 29, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200554/"],
#
#              'tprime1800_tZbW' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200646/"],
#
#              'tprime800_tHtH' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195633/"],
#
#              'tprime900_tHtH' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195742/"],
#
#              'tprime1000_tHtH' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195900/"],
#
#              'tprime1100_tHtH' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_182558/"],
#
#              'tprime1200_tHtH' : [2, 36, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200051/"],
#
#              'tprime1300_tHtH' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200229/"],
#
#              'tprime1400_tHtH' : [2, 33, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200324/"],
#
#              'tprime1500_tHtH' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200414/"],
#
#              'tprime1600_tHtH' : [2, 35, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200506/"],
#
#              'tprime1700_tHtH' : [2, 29, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200554/"],
#
#              'tprime1800_tHtH' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200646/"],
#
#              'tprime800_tHbW' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195633/"],
#
#              'tprime900_tHbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195742/"],
#
#              'tprime1000_tHbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195900/"],
#
#              'tprime1100_tHbW' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_182558/"],
#
#              'tprime1200_tHbW' : [2, 36, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200051/"],
#
#              'tprime1300_tHbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200229/"],
#
#              'tprime1400_tHbW' : [2, 33, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200324/"],
#
#              'tprime1500_tHbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200414/"],
#
#              'tprime1600_tHbW' : [2, 35, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200506/"],
#
#              'tprime1700_tHbW' : [2, 29, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200554/"],
#
#              'tprime1800_tHbW' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200646/"],
#
#              'tprime800_bWbW' : [2, 24, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195633/"],
#
#              'tprime900_bWbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195742/"],
#
#              'tprime1000_bWbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_195900/"],
#
#              'tprime1100_bWbW' : [2, 19, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_182558/"],
#
#              'tprime1200_bWbW' : [2, 36, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200051/"],
#
#              'tprime1300_bWbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200229/"],
#
#              'tprime1400_bWbW' : [2, 33, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200324/"],
#
#              'tprime1500_bWbW' : [2, 21, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200414/"],
#
#              'tprime1600_bWbW' : [2, 35, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200506/"],
#
#              'tprime1700_bWbW' : [2, 29, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200554/"],
#
#              'tprime1800_bWbW' : [2, 22, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170120_200646/"],
#
#              'SingleMuon_23Sep2016B' : [1, 4390, "/store/user/dmendis/ntuples2016/reminiaod/batch2/Muons/SingleMu1/"],
#              'SingleMuon_23Sep2016C' : [1, 1449, "/store/user/dmendis/ntuples2016/reminiaod/batch2/Muons/SingleMu2/"],
#              'SingleMuon_23Sep2016D' : [1, 2430, "/store/user/dmendis/ntuples2016/reminiaod/batch2/Muons/SingleMu3/"],
#              'SingleMuon_23Sep2016E' : [1, 2064, "/store/user/dmendis/ntuples2016/reminiaod/batch2/Muons/SingleMu4/"],
#              'SingleMuon_23Sep2016F' : [1, 1507, "/store/user/dmendis/ntuples2016/reminiaod/batch2/Muons/SingleMu5/"],
#              'SingleMuon_23Sep2016G' : [1, 3557, "/store/user/dmendis/ntuples2016/reminiaod/batch2/Muons/SingleMu6/"],
#              'SingleMuon_23Sep2016H' : [1, 3852, "/store/user/dmendis/ntuples2016/reminiaod/batch2/Muons/SingleMu7/"],
#              'SingleMuon_23Sep2016H_v2' : [1, 103, "/store/user/dmendis/ntuples2016/reminiaod/batch2/Muons/SingleMu7p1/"],
              }


pwd = os.getcwd()
test = pwd.split('test')[0]+'test/'


for sample in sorted(sampleDict.keys()):

  ## make necessary directories and get the number of jobs I am setting up for
  os.popen('mkdir -p '+options.directoryName+'/'+sample)
  os.popen('mkdir -p '+options.directoryName+'/output/'+sample)
  print sample
  if sampleDict[sample][1] % sampleDict[sample][0]  == 0:
    nentry = sampleDict[sample][1]/sampleDict[sample][0]
  else:
    nentry = math.floor(sampleDict[sample][1]/sampleDict[sample][0])+1

  nentry = int(nentry)

  ## make batch file
  batchInput = open('batchDummy.py')
  batchOutput = open(options.directoryName+'/batch_'+sample+'.jdl', 'w')
  for line in batchInput:
    line = line.replace('queue', str(nentry))
    line = line.replace('path', pwd+'/'+options.directoryName+'/output/'+sample)
    line = line.replace('exe', 'run_'+sample+'.sh')
    batchOutput.writelines(line)
  batchInput.close()
  batchOutput.close()

  ## make run file
  runInput = open('runDummy.py')
  runOutput = open(options.directoryName+'/run_'+sample+'.sh', 'w')
  for line in runInput: 
    line = line.replace('SAMPLE', sample)
    line = line.replace('PYTHON', pwd+'/'+options.directoryName)
    line = line.replace('OUTPUT', '/store/user/tmitchel/'+options.eosPath+'/'+options.directoryName+'/') 
    line = line.replace('MEH', options.directoryName)
#    line = line.replace('MEH', '')
    line = line.replace('TEST', test)
    runOutput.writelines(line)
  runInput.close()
  runOutput.close()

  optionList_mc = {'ISDATA' : False, 'ISPHOTON': False, 'JECPAYLOAD' : "Summer16_23Sep2016V4_MC", 'MAKETREE' : options.makeTree, 'ZDECAYMODE' : "zmumu", 'PUREWEIGHT' : True, 'FILTERSIGNAL' : False, 'SIGNALTYPE' : "", 'LEPTONIDSFS' : True, 'BTAGSFS' : True, 'BTAGEFFMAP' : "ttjets_bTagEff.root", 'LEPTONTRIGSFS' : True, 'DYNLOCORR' : False, 'MASSRECO' : False, 'CONTROLRECO' : False, 'SYST' : options.syst, 'LHEWTS' : True, 'OPTIMIZERECO' : False, 'PDFS' : 'bprime800_bZ', 'ISVV' : False, 'DODYSF' : False}

  optionList_data = {'ISDATA' : True, 'ISPHOTON': False, 'JECPAYLOAD' : "Summer16_23Sep2016V4_MC", 'MAKETREE' : options.makeTree, 'ZDECAYMODE' : "zmumu", 'PUREWEIGHT' : False, 'FILTERSIGNAL' : False, 'SIGNALTYPE' : "", 'LEPTONIDSFS' : False, 'BTAGSFS' : False, 'BTAGEFFMAP' : "ttjets_bTagEff.root", 'LEPTONTRIGSFS' : False, 'DYNLOCORR' : False, 'MASSRECO' : False, 'CONTROLRECO' : False, 'SYST' : options.syst, 'LHEWTS' : False, 'OPTIMIZERECO' : False, 'PDFS': 'none', 'ISVV' : False, 'DODYSF' : False}

  if 'Muon' in sample or 'Electron' in sample or 'Photon' in sample:
    if 'B' in sample or 'C' in sample or 'D' in sample:
      optionList_data['JECPAYLOAD'] = "Summer16_23Sep2016BCDV4_DATA"

    elif '6E' in sample or 'F' in sample:
      optionList_data['JECPAYLOAD'] = "Summer16_23Sep2016EFV4_DATA"

    elif 'G' in sample:
      optionList_data['JECPAYLOAD'] = "Summer16_23Sep2016GV4_DATA"

    elif '2016H' in sample:
      optionList_data['JECPAYLOAD'] = "Summer16_23Sep2016HV4_DATA"

  elif 'bprime' in sample or 'tprime' in sample:
    optionList_mc['FILTERSIGNAL'] = True
    optionList_mc['PDFS'] = "'"+str(sample)+"'"

    if optionList_mc['MAKETREE'] == False and 'bZbZ' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_bZbZ"
    elif optionList_mc['MAKETREE'] == False and 'bZbH' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_bZbH"
    elif optionList_mc['MAKETREE'] == False and 'bHbH' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_bHbH"
    elif optionList_mc['MAKETREE'] == False and 'bHtW' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_bHtW"
    elif optionList_mc['MAKETREE'] == False and 'bZtW' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_bZtW"
    elif optionList_mc['MAKETREE'] == False and 'tWtW' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_tWtW"

    elif optionList_mc['MAKETREE'] == False and 'tZtZ' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_tZtZ"
    elif optionList_mc['MAKETREE'] == False and 'tZtH' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_tZtH"
    elif optionList_mc['MAKETREE'] == False and 'tHtH' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_tHtH"
    elif optionList_mc['MAKETREE'] == False and 'tHbW' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_tHbW"
    elif optionList_mc['MAKETREE'] == False and 'tZbW' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_tZbW"
    elif optionList_mc['MAKETREE'] == False and 'bWbW' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_bWbW"


    if 'bprime800' in sample or 'bprime900' in sample or 'bprime1000' in sample:
      optionList_mc['BTAGEFFMAP'] = "bpbpbZb800_bTagEff.root"

    elif 'bprime1100' in sample or 'bprime1200' in sample or 'bprime1300' in sample or 'bprime1400' in sample:
      optionList_mc['BTAGEFFMAP'] = "bpbpbZb1200_bTagEff.root"
      
    elif 'bprime1500' in sample or 'bprime1700' in sample or 'bprime1800' in sample:
      optionList_mc['BTAGEFFMAP'] = "bpbpbZb1700_bTagEff.root"
    
    elif 'tprime800' in sample or 'tprime900' in sample or 'tprime1000' in sample:
      optionList_mc['BTAGEFFMAP'] = "bpbpbZb800_bTagEff.root"

    elif 'tprime1100' in sample or 'tprime1200' in sample or 'tprime1300' in sample or 'tprime1400' in sample:
      optionList_mc['BTAGEFFMAP'] = "bpbpbZb1200_bTagEff.root"
      
    elif 'tprime1500' in sample or 'tprime1600' in sample or 'tprime1700' in sample or 'tprime1800' in sample:
      optionList_mc['BTAGEFFMAP'] = "bpbpbZb1700_bTagEff.root"


  elif 'WW' in sample or 'WZ' in sample or 'ZZ' in sample:
    optionList_mc['BTAGEFFMAP'] = 'vv_bTagEff.root'
    optionList_mc['LHEWTS'] = False
    optionList_mc['PDFS'] = "'none'"
    optionList_mc['ISVV'] = True
    
  elif 'ttbar' in sample or 'ttZ' in sample:
    optionList_mc['BTAGEFFMAP'] = 'ttjets_bTagEff.root'
    optionList_mc['PDFS'] = "'ttbar'"
    
  elif 'dy' in sample:
    optionList_mc['BTAGEFFMAP'] = 'dy_bTagEff.root'
    optionList_mc['DYNLOCORR'] = True
    optionList_mc['PDFS'] = "'"+str(sample[:13])+"'"

  if 'Photon' in sample:
    optionList_data['ISPHOTON'] = True

  if 'Muon' in sample or 'Electron' in sample or 'Photon' in sample:
    optionDict = optionList_data
  else:
    optionDict = optionList_mc

  inputScript = open('dummy_os2lana.py')
  outputScript = open('dummy_temp.py', 'w')
  for line in inputScript:
    for key in optionDict:
      line = line.replace(key, str(optionDict[key]))
    outputScript.writelines(line)
  inputScript.close()
  outputScript.close()

  ## make python scripts to submit
  nfile = 1
  n = 1
  for entry in range(1, int(nentry)+2):
    inputfile = open('dummy_temp.py')
    outputfile = open(options.directoryName+'/'+sample+'/'+sample+'_'+str(entry)+'.py', 'w')
    director = 'root://eoscms.cern.ch/'
    if 'Muon' in sample or 'Electron' in sample or 'Photon' in sample:
      director = 'root://cmseos.fnal.gov/'
      #director = 'root://cms-xrd-global.cern.ch/'
      if 'Photon' in sample:
        replacer = director+sampleDict[sample][2]+'/os2lana_skim_'+str(n)+'.root'
        for i in range(1, sampleDict[sample][0]):
          n += 1
          replacer += "',\n'"+director+sampleDict[sample][2]+'/os2lana_skim_'+str(n)+'.root' 
      else:
        replacer = director+sampleDict[sample][2]+'/'+sampleDict[sample][2].split('/')[-2]+'_skim_'+str(n)+'.root'
        for i in range(1, sampleDict[sample][0]):
          n += 1
          replacer += "',\n'"+director+sampleDict[sample][2]+'/'+sampleDict[sample][2].split('/')[-2]+'_skim_'+str(n)+'.root' 

    else:
      if int(n/1000) < 10:
        replacer = director+sampleDict[sample][2]+'000'+str(int(n/1000))+'/B2GEDMNtuple_'+str(n)+'.root'
      else:
        replacer = director+sampleDict[sample][2]+'00'+str(int(n/1000))+'/B2GEDMNtuple_'+str(n)+'.root'
     
      for i in range(1, sampleDict[sample][0]):
        n += 1
        if n > sampleDict[sample][1]:
          break
        if int(n/1000) < 10:
          replacer+="',\n'"+director+sampleDict[sample][2]+'000'+str(int((n)/1000))+"/B2GEDMNtuple_"+str(n)+".root"
        else:
          replacer+="',\n'"+director+sampleDict[sample][2]+'00'+str(int((n)/1000))+"/B2GEDMNtuple_"+str(n)+".root"

#    if int(n/1000) < 10:
#      replacer = director+sampleDict[sample][2]+'000'+str(int(n/1000))+'/B2GEDMNtuple_'+str(n)+'.root'
#    else:
#      replacer = director+sampleDict[sample][2]+'00'+str(int(n/1000))+'/B2GEDMNtuple_'+str(n)+'.root'
#   
#    for i in range(1, sampleDict[sample][0]):
#      n += 1
#      if n > sampleDict[sample][1]:
#        break
#      if int(n/1000) < 10:
#        replacer+="',\n'"+director+sampleDict[sample][2]+'000'+str(int((n)/1000))+"/B2GEDMNtuple_"+str(n)+".root"
#      else:
#        replacer+="',\n'"+director+sampleDict[sample][2]+'00'+str(int((n)/1000))+"/B2GEDMNtuple_"+str(n)+".root"
    n+=1

    for line in inputfile:
      line = line.replace('path_to_sample', replacer)
      line = line.replace('CondOutputSkim', 'Skim_'+str(sample)+'_'+str(entry)+'.root')
      line = line.replace('CondOutputEvt', 'Evt_'+str(sample)+'_'+str(entry)+'.root')
      outputfile.writelines(line)

    nfile+=sampleDict[sample][0]
    inputfile.close()
    outputfile.close()
    if n > sampleDict[sample][1]:
      break

  ## make corresponding directory in eos if true
  if options.makeEOS:
    os.popen("eos root://cmseos.fnal.gov/ mkdir -p /store/user/tmitchel/"+options.eosPath+"/"+options.directoryName+"/"+sample)

  os.popen('tar -czf '+sample+'.tar.gz '+options.directoryName+'/'+sample)
  os.popen('xrdcp -f '+sample+'.tar.gz root://cmseos.fnal.gov//store/user/tmitchel/'+options.transfer+'/Zmumu/')
  os.popen('rm '+sample+'.tar.gz')
