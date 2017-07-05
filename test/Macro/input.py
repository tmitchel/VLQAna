#! /usr/bin/env python

# =====================================================
#  INPUTS		
# =====================================================
path = '/uscms_data/d3/tmitchel/update80X/CMSSW_8_0_20/src/Analysis/VLQAna/test/Macro/rootFiles/Zelel/June6_dropMerge/btag/'
pathS = '/uscms_data/d2/skhalil/MyVLQAna2/CMSSW_7_4_15_patch1/src/Analysis/VLQAna/test/CRAB_On_Skim/Histo/'
pathR = '/uscms_data/d3/dmendis/Rachitha2/CMSSW_7_4_16_patch2/src/Analysis/VLQAna/test/CRAB_0n_Skim/Histo/'

ch = 'CR_Zelel:'
if 'elel' in ch:
    suf = '_el'
elif 'mumu' in ch:
    suf = '_mu'
#f_Data_Oct2015 = TFile(path+'')
#f_Data_PromptReco = TFile(path+'SingleEG.root')
f_Data_PromptReco = TFile(path+'data.root')

f_DY100to200 = TFile(path+'dy_pt100to250.root')
f_DY200to400 = TFile(path+'dy_pt250to400.root')
f_DY400to600 = TFile(path+'dy_pt400to650.root')
f_DY600to800 = TFile(path+'dy_pt650toInf.root')
#f_DYmcnlo    = TFile(path+'ele_dy_inclusive.root')

# # f_WJ100to200 = TFile(path+'ele_WJetsToLNu_HT-100To200_os2lana_v1_ele_.root')
# # f_WJ200to400 = TFile(path+'ele_WJetsToLNu_HT-200To400_os2lana_v1_ele_.root')
# # f_WJ400to600 = TFile(path+'ele_WJetsToLNu_HT-400To600_os2lana_v1_ele_.root')
# # f_WJ600to800 = TFile(path+'ele_WJetsToLNu_HT-600To800_os2lana_v1_ele_.root')
# # f_WJ800to1200 = TFile(path+'ele_WJetsToLNu_HT-800To1200_os2lana_v1_ele_.root')
# # f_WJ1200to2500 = TFile(path+'ele_WJetsToLNu_HT-1200To2500_os2lana_v1_ele_.root')
# # f_WJ2500toInf = TFile(path+'ele_WJetsToLNu_HT-2500ToInf_os2lana_v1_ele_.root')

# # f_ST_tW_top     = TFile(path+'ele_ST_tW_5f_top_powheg-pythia8_os2lana_v1_ele_.root')
# # f_ST_tW_antitop = TFile(path+'ele_ST_tW_5f_antitop_powheg-pythia8_os2lana_v1_ele_.root')
# # f_ST_t          = TFile(path+'ele_ST_t_4f_amcatnlo-pythia8_os2lana_v1_ele_.root')
# # f_ST_t_ex1      = TFile(path+'ele_ST_t_4f_amcatnlo-pythia8_ext1_os2lana_v1_ele_.root')
# # f_ST_s          = TFile(path+'ele_ST_s_4f_amcatnlo-pythia8_os2lana_v1_ele_.root')

f_ZZTo2L2Nu     = TFile(path+'ZZ.root')
f_WZTo2L2Q      = TFile(path+'WZ.root')
f_WWTo2L2Nu   = TFile(path+'WW.root')
#f_WZTo3LNu      = TFile(path+'ele_WZto3.root')
#f_ZZTo4L            = TFile(path+'ele_ZZto4.root')

#f_ttbar         = TFile(path+'ttbar_Tune1.root')
f_ttbar = TFile(path+'ttbar.root')
f_BpBp_bZbZ_800 = TFile(path+'bprime800_bZ.root')
f_BpBp_bZbZ_900 = TFile(path+'bprime900_bZ.root')
f_BpBp_bZbZ_1000 = TFile(path+'bprime1000_bZ.root')
f_BpBp_bZbZ_1100 = TFile(path+'bprime1100_bZ.root')
f_BpBp_bZbZ_1200 = TFile(path+'bprime1200_bZ.root')
f_BpBp_bZbZ_1300 = TFile(path+'bprime1300_bZ.root')
f_BpBp_bZbZ_1400 = TFile(path+'bprime1400_bZ.root')
f_BpBp_bZbZ_1500 = TFile(path+'bprime1500_bZ.root')
#f_BpBp_bZbZ_1600 = TFile(path+'bprime1600_bZ.root')
f_BpBp_bZbZ_1700 = TFile(path+'bprime1700_bZ.root')
f_BpBp_bZbZ_1800 = TFile(path+'bprime1800_bZ.root')
#
f_BpBp_bZbH_800 = TFile(path+'bprime800_bH.root')
f_BpBp_bZbH_900 = TFile(path+'bprime900_bH.root')
f_BpBp_bZbH_1000 = TFile(path+'bprime1000_bH.root')
f_BpBp_bZbH_1100 = TFile(path+'bprime1100_bH.root')
f_BpBp_bZbH_1200 = TFile(path+'bprime1200_bH.root')
f_BpBp_bZbH_1300 = TFile(path+'bprime1300_bH.root')
f_BpBp_bZbH_1400 = TFile(path+'bprime1400_bH.root')
f_BpBp_bZbH_1500 = TFile(path+'bprime1500_bH.root')
#f_BpBp_bZbH_1600 = TFile(path+'bprime1600_bH.root')
f_BpBp_bZbH_1700 = TFile(path+'bprime1700_bH.root')
f_BpBp_bZbH_1800 = TFile(path+'bprime1800_bH.root')


#===== cross sections (pb)==========

Top_xs            = 831.76  *gSF
DY100to200_xs     = 83.12#147.4   *gSF #*1.23
DY200to400_xs     = 3.047#40.99   *gSF #*1.23
DY400to600_xs     = 0.3921#5.678   *gSF #*1.23
DY600to800_xs     = 0.03636#1.363   *gSF #*1.23
DY800to1200_xs    = .6759   *gSF *1.23
DY1200to2500_xs   = 0.116   *gSF *1.23
DY2500toInf_xs    = .002592   *gSF *1.23

DY_xs             = 6025.2  *gSF
WJ100to200_xs     = 1345.0  *gSF *1.21 
WJ200to400_xs     = 359.7   *gSF *1.21
WJ400to600_xs     = 48.9    *gSF *1.21
WJ600to800_xs     = 12.05   *gSF *1.21
WJ800to1200_xs    = 5.501   *gSF *1.21
WJ1200to2500_xs   = 1.329   *gSF *1.21
WJ2500toInf_xs    = 0.03216 *gSF *1.21
ST_tW_top_xs      = 35.6    *gSF
ST_tW_antitop_xs  = 35.6    *gSF 
ST_t_xs           = 70.69   *gSF
ST_s_xs           = 3.36    *gSF 

ZZTo2L2Nu_xs      = 16.91#0.564   *gSF
WZTo2L2Q_xs       = 46.74    *gSF
WWTo2L2Nu_xs      = 118.7  *gSF
#WZTo3LNu_xs       = 4.42965 *gSF
#ZZTo4L_xs         = 1.212   *gSF
TpTp800_xs        = 0.196   *gSF  
TpTp900_xs        = 0.0903  *gSF
TpTp1000_xs       = 0.044   *gSF
TpTp1200_xs       = 0.0118  *gSF
TpTp1500_xs       = 0.00200 *gSF

BpBp800_xs        = 0.196   *gSF
BpBp900_xs        = 0.0903  *gSF
BpBp1000_xs       = 0.044   *gSF
BpBp1100_xs       = 0.0224  *gSF
BpBp1200_xs       = 0.0118  *gSF
BpBp1300_xs       = 0.00639 *gSF
BpBp1400_xs       = 0.00354 *gSF
BpBp1500_xs       = 0.00200 *gSF
BpBp1600_xs       = 0.001148*gSF
BpBp1700_xs       = 0.000666*gSF
BpBp1800_xs       = 0.000391*gSF

#BpBp800_xs = 1.
#BpBp900_xs = 1.
##BpBp1000_xs = 1.
#BpBp1100_xs = 1.
#BpBp1200_xs = 1.
#BpBp1300_xs = 1.
#BpBp1400_xs = 1.
#BpBp1500_xs = 1.
#BpBp1600_xs = 1.
#BpBp1700_xs = 1.
#BpBp1800_xs = 1.
TpTp800_xs= 1.
TpTp1000_xs = 1.
TpTp1200_xs = 1.

#===== Number of generated events ======

Top_num          = f_ttbar.Get("allEvents/hEventCount_nowt").GetBinContent(1)
DY100to200_num   = f_DY100to200.Get("allEvents/hEventCount_nowt").GetBinContent(1)
DY200to400_num   = f_DY200to400.Get("allEvents/hEventCount_nowt").GetBinContent(1)
DY400to600_num   = f_DY400to600.Get("allEvents/hEventCount_nowt").GetBinContent(1)
DY600to800_num   = f_DY600to800.Get("allEvents/hEventCount_nowt").GetBinContent(1)

ZZTo2L2Nu_num    =  f_ZZTo2L2Nu.Get("allEvents/hEventCount_nowt").GetBinContent(1)
WZTo2L2Q_num     =  f_WZTo2L2Q.Get("allEvents/hEventCount_nowt").GetBinContent(1)
WWTo2L2Nu_num    =  f_WWTo2L2Nu.Get("allEvents/hEventCount_nowt").GetBinContent(1)

BpBp800_bZ_num      =  f_BpBp_bZbZ_800.Get("ana/signalEvts").GetBinContent(1)
BpBp900_bZ_num      =  f_BpBp_bZbZ_900.Get("ana/signalEvts").GetBinContent(1)
BpBp1000_bZ_num     =  f_BpBp_bZbZ_1000.Get("ana/signalEvts").GetBinContent(1)
BpBp1100_bZ_num     =  f_BpBp_bZbZ_1100.Get("ana/signalEvts").GetBinContent(1)
BpBp1200_bZ_num     =  f_BpBp_bZbZ_1200.Get("ana/signalEvts").GetBinContent(1)
BpBp1300_bZ_num     =  f_BpBp_bZbZ_1300.Get("ana/signalEvts").GetBinContent(1)
BpBp1400_bZ_num     =  f_BpBp_bZbZ_1400.Get("ana/signalEvts").GetBinContent(1)
BpBp1500_bZ_num     =  f_BpBp_bZbZ_1500.Get("ana/signalEvts").GetBinContent(1)
#BpBp1600_bZ_num     =  f_BpBp_bZbZ_1600.Get("ana/signalEvts").GetBinContent(1)
BpBp1700_bZ_num     =  f_BpBp_bZbZ_1700.Get("ana/signalEvts").GetBinContent(1)
BpBp1800_bZ_num     =  f_BpBp_bZbZ_1800.Get("ana/signalEvts").GetBinContent(1)
BpBp800_bH_num      =  f_BpBp_bZbH_800.Get("ana/signalEvts").GetBinContent(1)
BpBp900_bH_num      =  f_BpBp_bZbH_900.Get("ana/signalEvts").GetBinContent(1)
BpBp1000_bH_num     =  f_BpBp_bZbH_1000.Get("ana/signalEvts").GetBinContent(1)
BpBp1100_bH_num     =  f_BpBp_bZbH_1100.Get("ana/signalEvts").GetBinContent(1)
BpBp1200_bH_num     =  f_BpBp_bZbH_1200.Get("ana/signalEvts").GetBinContent(1)
BpBp1300_bH_num     =  f_BpBp_bZbH_1300.Get("ana/signalEvts").GetBinContent(1)
BpBp1400_bH_num     =  f_BpBp_bZbH_1400.Get("ana/signalEvts").GetBinContent(1)
BpBp1500_bH_num     =  f_BpBp_bZbH_1500.Get("ana/signalEvts").GetBinContent(1)
#BpBp1600_bH_num     =  f_BpBp_bZbH_1600.Get("ana/signalEvts").GetBinContent(1)
BpBp1700_bH_num     =  f_BpBp_bZbH_1700.Get("ana/signalEvts").GetBinContent(1)
BpBp1800_bH_num     =  f_BpBp_bZbH_1800.Get("ana/signalEvts").GetBinContent(1)
#
#print BpBp800_bZ_num, BpBp800_bH_num


# Legend
leg = TLegend(0.74,0.90,0.96,0.50)
#leg = TLegend(0.5, .90, 0.74, .65)
leg.SetBorderSize(0)
leg.SetFillColor(10)
leg.SetLineColor(10)
leg.SetLineWidth(0)


# =====================================================
#  FUNCTIONS		
# =====================================================

def setTitle(hs,xTitle):
    y = hs.GetYaxis()
    x = hs.GetXaxis()
    y.SetTitle("Events / Bin")
    x.SetTitle(xTitle)
    y.SetLabelSize(0.05)
    y.SetTitleSize(0.07)
    y.SetTitleOffset(0.6)
    y.SetTitleFont(42)
    x.SetTitleSize(0.05)
    x.SetTitleFont(42)

def prepareRatio(h_ratio, h_ratiobkg, scale, xTitle):
    h_ratio.SetTitle("")
    h_ratio.GetYaxis().SetTitle("Data / Bkg")
    h_ratio.GetXaxis().SetTitle(xTitle)   
    h_ratio.SetMarkerStyle(8) 
    h_ratio.SetMaximum(2)#3
    h_ratio.SetMinimum(0)#-1
    h_ratio.GetYaxis().SetLabelSize(0.06*scale)#0.06
    h_ratio.GetYaxis().SetTitleOffset(1.00/scale*0.5)
    h_ratio.GetYaxis().SetTitleSize(0.08*scale)#0.08
    h_ratio.GetYaxis().SetTitleFont(42)
    h_ratio.GetXaxis().SetLabelSize(0.06*scale)#0.06
    h_ratio.GetXaxis().SetTitleOffset(0.45*scale)
    h_ratio.GetXaxis().SetTitleSize(0.09*scale)#0.09
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetXaxis().SetNdivisions(510)
    h_ratio.SetTickLength(0.06,"X")
    h_ratio.SetTickLength(0.05,"Y")

    ## The uncertainty band
    h_ratio_bkg.SetMarkerSize(0)
    h_ratio_bkg.SetFillColor(kGray+1)
    h_ratio_bkg.GetYaxis().SetLabelSize(0.6*scale)
    h_ratio_bkg.GetYaxis().SetTitleOffset(1.00/scale*0.6)
    h_ratio_bkg.GetYaxis().SetTitleSize(0.08*scale)
    h_ratio_bkg.GetYaxis().SetTitleFont(42)
    h_ratio_bkg.GetXaxis().SetLabelSize(0.08*scale)
    h_ratio_bkg.GetXaxis().SetTitleOffset(0.45*scale)
    h_ratio_bkg.GetXaxis().SetTitleSize(0.09*scale)
    h_ratio_bkg.GetYaxis().SetNdivisions(505)
    h_ratio_bkg.GetXaxis().SetNdivisions(510)
    h_ratio_bkg.SetTickLength(0.05,"X")
    h_ratio_bkg.SetTickLength(0.05,"y")
    h_ratio_bkg.SetTitle("")    
    h_ratio_bkg.SetMaximum(2)
    h_ratio_bkg.SetMinimum(0)
    
def overUnderFlow(hist):
    xbins = hist.GetNbinsX()
    hist.SetBinContent(xbins, hist.GetBinContent(xbins)+hist.GetBinContent(xbins+1))
    hist.SetBinContent(1, hist.GetBinContent(0)+hist.GetBinContent(1))
    hist.SetBinError(xbins, TMath.Sqrt(TMath.Power(hist.GetBinError(xbins),2)+TMath.Power(hist.GetBinError(xbins+1),2)))
    hist.SetBinError(1, TMath.Sqrt(TMath.Power(hist.GetBinError(0),2)+TMath.Power(hist.GetBinError(1),2)))
    hist.SetBinContent(xbins+1, 0.)
    hist.SetBinContent(0, 0.)
    hist.SetBinError(xbins+1, 0.)
    hist.SetBinError(0, 0.)
    
def setCosmetics(hist, legname, hname, color):
    hist.Rebin(rebinS)
    hist.SetLineColor(color)
    hist.SetName(hname)
    #hist.GetXaxis().SetRangeUser(0,600)
    if 'Data' in hname:
        leg.AddEntry(hist, legname, 'pl')
        hist.SetMarkerStyle(8)
    elif 'tZ' in hname or 'BB' in hname:
        hist.SetLineWidth(2)
        leg.AddEntry(hist, legname, 'l')
    else:
        hist.SetFillColor(color)
        leg.AddEntry(hist, legname, 'f')

        
def getHisto( label, leg, dir, var, Samples, color, verbose) :
    histos = []
    for iSample in Samples :
        ifile = iSample[0]
        xs = iSample[1]
        nevt = iSample[2]
        lumi = iSample[3]
        readname = dir+'/'+var
        hist  = ifile.Get( readname ).Clone()
        if verbose:
            print 'file: {0:<20}, histo:{1:<10}, integral before weight:{2:<3.3f}, nEntries:{3:<3.0f}, weight:{4:<2.3f}'.format(
                ifile.GetName(),    
                hist.GetName(),
                hist.Integral(), hist.GetEntries(), xs * lumi /nevt
                )
            print xs, lumi, nevt, xs*lumi/nevt
        hist.Sumw2()    
        # if ifile == 'f_Data_PromptReco':
        #     xs = 0
        #     hist.SetMarkerColorAlpha(color, 100)
        hist.Scale( xs * lumi /nevt)
        hist.Rebin(50)
        hist.SetAxisRange(0.,2000.,"X")
        histos.append( hist )
        
    histo = histos[0]
    setCosmetics(histo, leg, label+var, color) 
    
    for ihisto in range(1, len(histos) ):
        #print 'ihisto =', ihisto, 'integral', histos[ihisto].Integral(), ', entries', histos[ihisto].GetEntries()
        histo.Add( histos[ihisto] )
        #print 'after addition', histo.Integral()

    #if verbose:    
     #   print 'newName: {0:<5}, Entries:{1:5.2f},  newIntegral: {2:5.2f}'.format(label+var, histo.GetEntries(), histo.Integral() )   
    return histo

