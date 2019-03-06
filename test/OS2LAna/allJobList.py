#! /bin/python

def list_Zelel_data():
    jobList = [
        ['/DoubleEG/skhi-RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0-642b48c1c707cb5cfa93fc168fde448b/USER', 'DoubleEG-Run2016', '20'],
        ]
    return jobList

def list_Zmumu_data():
    jobList = [
        ['/DoubleMuon/skhi-RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0-642b48c1c707cb5cfa93fc168fde448b/USER','DoubleMuon', '10'],
        ]
    return jobList

def list_photon_data():
  jobList = [
      ['/SinglePhoton/vorobiev-SinglePhoton_Run2016B-23Sep2016-v3_B2GAnaFW_80X_v2p4-6f92d8e9b2717da1daa65a0e07f84d5f/USER', 'Photon2016B', '10'],
      ['/SinglePhoton/vorobiev-SinglePhoton_Run2016D-23Sep2016-v1_B2GAnaFW_80X_v2p4-6f92d8e9b2717da1daa65a0e07f84d5f/USER', 'Photon2016D', '10'],
      ['/SinglePhoton/vorobiev-SinglePhoton_Run2016E-23Sep2016-v1_B2GAnaFW_80X_v2p4-961c7d882d8721e72fac616aaa90ecc1/USER', 'Photon2016E', '10'],
      ['/SinglePhoton/vorobiev-SinglePhoton_Run2016F-23Sep2016-v1_B2GAnaFW_80X_v2p4-e26b2444814d18badce3899570108664/USER', 'Photon2016F', '10'],
      ['/SinglePhoton/vorobiev-SinglePhoton_Run2016G-23Sep2016-v1_B2GAnaFW_80X_v2p4-b45603bfa955d854bdcb9af322c0b037/USER', 'Photon2016G', '10'],
      ['/SinglePhoton/vorobiev-SinglePhoton_Run2016H-PromptReco-v2_B2GAnaFW_80X_v2p4-376a23645e94877b22a7f32873431514/USER', 'Photon2016H_v2', '10'],
      ['/SinglePhoton/vorobiev-SinglePhoton_Run2016H-PromptReco-v3_B2GAnaFW_80X_v2p4-376a23645e94877b22a7f32873431514/USER', 'Photon2016H_v3', '10'],
      ]
  return jobList


def list_MC(): 
    jobList = [
        ['/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-a300e501c1a543433113bdc094d47173/USER', 'DY_Pt-100to250', '10'],
        ['/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-a300e501c1a543433113bdc094d47173/USER', 'DY_Pt-250to400', '10'],
        ['/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-a300e501c1a543433113bdc094d47173/USER', 'DY_Pt-400to650', '10'],
        ['/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-a300e501c1a543433113bdc094d47173/USER', 'DY_Pt-650ToInf', '10'],
        ['/TT_TuneCUETP8M1_13TeV-powheg-pythia8/grauco-B2GAnaFW_80X_V2p1-edbed0685401a5848e7d61871b3a63d8/USER', 'TTbar_TuneCUETP8M1_powheg', '20'],
        ['/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TTbar_TuneCUETP8M2T4_powheg', '10'],
        
        ['/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TprimeTprime_M-800', '2'],
        ['/TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TprimeTprime_M-900', '2'],
        ['/TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TprimeTprime_M-1000', '2'],
        ['/TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TprimeTprime_M-1100', '2'],
        ['/TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TprimeTprime_M-1200', '2'],
        ['/TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TprimeTprime_M-1300', '2'],
        ['/TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TprimeTprime_M-1400', '2'],
        ['/TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TprimeTprime_M-1500', '2'],
        ['/TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TprimeTprime_M-1600', '2'],
        ['/TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TprimeTprime_M-1700', '2'],
        ['/TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TprimeTprime_M-1800', '2'],
        ['/BprimeBprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'BprimeBprime_M-800', '2'],
        ['/BprimeBprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'BprimeBprime_M-900', '2'],
        ['/BprimeBprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'BprimeBprime_M-1000', '2'],
        ['/BprimeBprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'BprimeBprime_M-1100', '2'],
        ['/BprimeBprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'BprimeBprime_M-1200', '2'],
        ['/BprimeBprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'BprimeBprime_M-1300', '2'],
        ['/BprimeBprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'BprimeBprime_M-1400', '2'],
        ['/BprimeBprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'BprimeBprime_M-1500', '2'],
        ['/BprimeBprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'BprimeBprime_M-1600', '2'],
        ['/BprimeBprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'BprimeBprime_M-1700', '2'],
        ['/BprimeBprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'BprimeBprime_M-1800', '2'],
         
        ]   
    return jobList
