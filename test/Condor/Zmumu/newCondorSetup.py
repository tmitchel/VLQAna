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

(options, args) = parser.parse_args()

# sampleList = ['ttbar_pythia', 'ttbar_pythia_ext', 'ZZ', 'ZZ_ext', 'WW', 'WZ', 'dy_pt100to250', 'dy_pt100to250_ext', 'dy_pt250to400', 'dy_pt250to400_ext', 'dy_pt400to650', 'dy_pt400to650_ext', 'dy_pt650toInf', 'dy_pt650toInf_ext', 'bprime800', 'bprime900', 'bprime1000', 'bprime1100', 'bprime1200', 'bprime1300', 'bprime1400', 'bprime1500', 'bprime1700', 'bprime1800', 'SingleMuon_PromptH_v2', 'SingleMuon_PromptH_v3', 'SingleMuon_23Sep2016B_v3', 'SingleMuon_23Sep2016C_v1', 'SingleMuon_23Sep2016D_v1', 'SingleMuon_23Sep2016E_v1', 'SingleMuon_23Sep2016F_v1', 'SingleMuon_23Sep2016G_v1']



sampleDict = {'ttbar_pythia' : [5, 984, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/161222_110143/"],

              'ttbar_pythia_ext' : [5, 1192, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_backup_v2p4/170304_173805/"],

              'dy_pt100to250' : [5, 34, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170121_174920/"],

              'dy_pt100to250_ext' : [5, 34, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_ext1_v2p4/170122_184903/"],

              'dy_pt250to400' : [2, 11, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p4/170121_175352/"],

              'dy_pt250to400_ext' : [2, 9, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_ext1_v2p4/170122_185048/"],

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


#              'SingleMuon_23Sep2016B_v3' : [5, 1277, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016B-23Sep2016-v3_B2GAnaFW_80X_v2p4/161221_152050/"],
#              'SingleMuon_23Sep2016C_v1' : [5, 713, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016C-23Sep2016-v1_B2GAnaFW_80X_v2p4/161221_155142/"],
#              'SingleMuon_23Sep2016D_v1' : [5, 1060, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016D-23Sep2016-v1_B2GAnaFW_80X_v2p4/161221_160134/"],
#              'SingleMuon_23Sep2016E_v1' : [5, 886, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016E-23Sep2016-v1_B2GAnaFW_80X_v2p4/161221_162521/"],
#              'SingleMuon_23Sep2016F_v1' : [5, 768, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016F-23Sep2016-v1_B2GAnaFW_80X_v2p4/161221_170915/"],
#              'SingleMuon_23Sep2016G_v1' : [5, 1430, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016G-23Sep2016-v1_B2GAnaFW_80X_v2p4/161221_173126/"],
#              'SingleMuon_PromptH_v2' : [5, 1928, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016H-PromptReco-v2_B2GAnaFW_80X_v2p4/161221_180745/"],
#              'SingleMuon_PromptH_v3' : [5, 52, "/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016H-PromptReco-v3_B2GAnaFW_80X_v2p4/161221_180928/"],
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
    line = line.replace('TEST', test)
    runOutput.writelines(line)
  runInput.close()
  runOutput.close()

  optionList_mc = {'ISDATA' : False, 'JECPAYLOAD' : "Summer16_23Sep2016V4_MC", 'MAKETREE' : options.makeTree, 'ZDECAYMODE' : "zmumu", 'PUREWEIGHT' : True, 'FILTERSIGNAL' : False, 'SIGNALTYPE' : "", 'LEPTONIDSFS' : True, 'BTAGSFS' : True, 'BTAGEFFMAP' : "ttjets_bTagEff.root", 'LEPTONTRIGSFS' : True, 'DYNLOCORR' : False, 'MASSRECO' : True, 'CONTROLRECO' : True, 'SYST' : options.syst, 'LHEWTS' : True, 'OPTIMIZERECO' : True, 'PDFS' : 'bprime800_bZ', 'ISVV' : False}

  optionList_data = {'ISDATA' : True, 'JECPAYLOAD' : "Summer16_23Sep2016V4_MC", 'MAKETREE' : options.makeTree, 'ZDECAYMODE' : "zmumu", 'PUREWEIGHT' : False, 'FILTERSIGNAL' : False, 'SIGNALTYPE' : "", 'LEPTONIDSFS' : False, 'BTAGSFS' : False, 'BTAGEFFMAP' : "ttjets_bTagEff.root", 'LEPTONTRIGSFS' : False, 'DYNLOCORR' : False, 'MASSRECO' : True, 'CONTROLRECO' : True, 'SYST' : options.syst, 'LHEWTS' : False, 'OPTIMIZERECO' : False, 'PDFS': 'none', 'ISVV' : False}

  if 'Muon' in sample:
    if 'B_v' in sample or 'C_v' in sample or 'D_v' in sample:
      optionList_data['JECPAYLOAD'] = "Summer16_23Sep2016BCDV4_DATA"

    elif 'E_v' in sample or 'F_v' in sample:
      optionList_data['JECPAYLOAD'] = "Summer16_23Sep2016EFV4_DATA"

    elif 'G_v' in sample:
      optionList_data['JECPAYLOAD'] = "Summer16_23Sep2016GV4_DATA"

    elif 'Prompt' in sample:
      optionList_data['JECPAYLOAD'] = "Summer16_23Sep2016HV4_DATA"

  elif 'bprime' in sample:
    optionList_mc['FILTERSIGNAL'] = True
    optionList_mc['PDFS'] = "'"+str(sample)+"'"

    if optionList_mc['MAKETREE'] == False and 'bZbZ' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_bZbZ"
    elif optionList_mc['MAKETREE'] == False and 'bZbH' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_bZbH"
    elif optionList_mc['MAKETREE'] == False and 'bHbH' in sample:
      optionList_mc['SIGNALTYPE'] = "EvtType_MC_bHbH"

    if 'bprime800' in sample or 'bprime900' in sample or 'bprime1000' in sample:
      optionList_mc['BTAGEFFMAP'] = "bpbpbZb800_bTagEff.root"

    elif 'bprime1100' in sample or 'bprime1200' in sample or 'bprime1300' in sample or 'bprime1400' in sample:
      optionList_mc['BTAGEFFMAP'] = "bpbpbZb1200_bTagEff.root"
      
    elif 'bprime1500' in sample or 'bprime1700' in sample or 'bprime1800' in sample:
      optionList_mc['BTAGEFFMAP'] = "bpbpbZb1700_bTagEff.root"
      
  elif 'WW' in sample or 'WZ' in sample or 'ZZ' in sample:
    optionList_mc['BTAGEFFMAP'] = 'vv_bTagEff.root'
    optionList_mc['LHEWTS'] = False
    optionList_mc['PDFS'] = "'none'"
    optionList_mc['ISVV'] = True
    
  elif 'ttbar' in sample:
    optionList_mc['BTAGEFFMAP'] = 'ttjets_bTagEff.root'
    optionList_mc['PDFS'] = "'ttbar'"
    
  elif 'dy' in sample:
    optionList_mc['BTAGEFFMAP'] = 'dy_bTagEff.root'
    optionList_mc['DYNLOCORR'] = False
    optionList_mc['PDFS'] = "'"+str(sample[:13])+"'"

  if 'Muon' in sample:
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
    # if 'Muon' in sample:
    #   if 'B_v' in sample or 'C_v' in sample or 'D_v' in sample:
    #     inputfile = open('dummy_MuonBCD.py')
    #   elif 'E_v' in sample or 'F_v' in sample:
    #     inputfile = open('dummy_MuonEF.py')
    #   elif 'G_v' in sample:
    #     inputfile = open('dummy_MuonG.py')
    #   elif 'Prompt' in sample:
    #     inputfile = open('dummy_MuonH.py')
    # elif 'bprime' in sample:
    #   if 'bprime800' in sample or 'bprime900' in sample or 'bprime1000' in sample:
    #     inputfile = open('dummy_lowSig.py')
    #   elif 'bprime1100' in sample or 'bprime1200' in sample or 'bprime1300' in sample or 'bprime1400' in sample:
    #     inputfile = open('dummy_medSig.py')
    #   elif 'bprime1500' in sample or 'bprime1700' in sample or 'bprime1800' in sample:
    #     inputfile = open('dummy_highSig.py')
    # elif 'WW' in sample or 'WZ' in sample or 'ZZ' in sample:
    #   inputfile = open('dummy_vv.py')
    # elif 'ttbar' in sample:
    #   inputfile = open('dummy_ttbar.py')
    # elif 'dy' in sample:
    #   inputfile = open('dummy_dy.py')

    inputfile = open('dummy_temp.py')
    outputfile = open(options.directoryName+'/'+sample+'/'+sample+'_'+str(entry)+'.py', 'w')
    if int(n/1000) < 10:
      replacer = 'root://eoscms.cern.ch/'+sampleDict[sample][2]+'000'+str(int(n/1000))+'/B2GEDMNtuple_'+str(n)+'.root'
    else:
      replacer = 'root://eoscms.cern.ch/'+sampleDict[sample][2]+'00'+str(int(n/1000))+'/B2GEDMNtuple_'+str(n)+'.root'

    for i in range(1, sampleDict[sample][0]):
      n += 1
      if n > sampleDict[sample][1]:
        break
      if int(n/1000) < 10:
        replacer+="',\n'root://eoscms.cern.ch/"+sampleDict[sample][2]+'000'+str(int((n)/1000))+"/B2GEDMNtuple_"+str(n)+".root"
      else:
        replacer+="',\n'root://eoscms.cern.ch/"+sampleDict[sample][2]+'00'+str(int((n)/1000))+"/B2GEDMNtuple_"+str(n)+".root"
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

