#! /usr/bin/env python
import sys
import os
import subprocess
from array import array
from ROOT import TH1D,TH2D,TFile,TMath,TCanvas,THStack,TLegend,TPave,TLine,TLatex, TF1, TGraph, TMultiGraph
from ROOT import gROOT,gStyle,gPad,gStyle
from ROOT import Double,kBlue,kRed,kOrange,kMagenta,kYellow,kCyan,kGreen,kGray,kBlack,kTRUE

gROOT.Macro("~/rootlogon.C")
gStyle.SetOptStat(0)

# ===============
# options
# ===============
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--Lumi', metavar='D', type='float', action='store',
                  default=35900.,#12900., #2263., # 538.134,#2200.,
                  dest='Lumi',
                  help='Data Luminosity in pb-1')

parser.add_option('--globalSF', metavar='SF', type='float',
                  default = 1.0,
                  dest='globalSF',
                  help='Global trigger SF (%default default)')

parser.add_option('--plotDir', metavar='P', type='string', action='store',
                  default='New',
                  dest='plotDir',
                  help='output directory of plots')

parser.add_option('--skimType', metavar='S', type='string', action='store',
                  default='CR_Zmumu',
                  dest='skimType',
                  help='Skim type: CR_Zelel, CR_Zmumu, SR_Zelel, SR_Zmumu')

parser.add_option('--processDir', metavar='pD', type='string', action='store',
                  default='massReco',
                  dest='processDir',
                  help='directory to read histograms from')

parser.add_option('--var', metavar='T', type='string', action='store',
                  default='pt_zmumu_pre',#cutflow, st
                  dest='var',
                  help='variable to plot')

parser.add_option('--sys', metavar='T', type='string', action='store',
                  default='nominal',
                  dest='sys',
                  help='nominal, BTagSFup, BTagSFdown, ScaleUp, ScaleDown, MatchUp, MatchDown')

parser.add_option('--verbose',action='store_true',
                  default=False,
                  dest='verbose',
                  help='verbose switch')

parser.add_option('--rebin', metavar='T', type='int', action='store',
                  default='1',
                  dest='rebin',
                  help='rebin the histograms')

parser.add_option('--drawEff', action='store_true', 
                  default=False,
                  dest='drawEff',
                  help='draw eff and sig'
)
parser.add_option('--forLimits', action='store',
                  default=True,
                  dest='forLimits',
                  help='make template for limits'
)

(options,args) = parser.parse_args()
# ==========end: options =============
var = options.var
lumi = options.Lumi
gSF  = options.globalSF
rebinS = options.rebin
pDir = options.processDir
plotDir = options.plotDir
skimType = options.skimType
verbose = options.verbose
forLimits =  options.forLimits
pDir2 = 'anaH/sig'

if 'elel' in skimType: title = 'e^{#pm}e^{#mp}+jets'
elif 'mumu' in skimType: title = '#mu^{#pm}#mu^{#mp}+jets'
else: title = ''
# ==========add the input============

execfile("input.py")

# === prepare the input labels and legends ===========
dataLabel     = 'Data_'
topLabel      = 'Top_'
dyLabel       = 'DY_'
wjLabel       = 'WJets_'
sTLabel       = 'ST_'
vvLabel       = 'VV_'

dataLeg       = 'Data'
topLeg        = 't#bar{t}'
dyLeg         = 'Drell-Yan'
wjLeg         = 'W+Jets'
sTLeg         = 'Single top'
vvLeg         = 'Diboson'
# === create structure ============
data     = [
#            [f_Data_Oct2015, 1., 1., 1.], # this corresponds to .53 fb-1
            [f_Data_PromptReco, 0., 1., 1.],  #addign this should give 2.2 fb-1  
           ]

top      = [[f_ttbar,         Top_xs,            Top_num,            lumi]]

#dy       = [[f_DYmcnlo,       DY_xs,             DYmcnlo_num,        lumi]]

dy       = [
             [f_DY100to200,    DY100to200_xs,      DY100to200_num,    lumi],
             [f_DY200to400,    DY200to400_xs,      DY200to400_num,    lumi],
             [f_DY400to600,    DY400to600_xs,      DY400to600_num,    lumi],
             [f_DY600to800,    DY600to800_xs,      DY600to800_num,    lumi],
           ]

# wjets    = [
#             [f_WJ100to200,    WJ100to200_xs,      WJ100to200_num,    lumi],
#             [f_WJ200to400,    WJ200to400_xs,      WJ200to400_num,    lumi],
#             [f_WJ400to600,    WJ400to600_xs,      WJ400to600_num,    lumi],
#             [f_WJ600to800,    WJ600to800_xs,      WJ600to800_num,    lumi],     
#             [f_WJ800to1200,   WJ800to1200_xs,     WJ800to1200_num,   lumi],
#             [f_WJ1200to2500,  WJ1200to2500_xs,    WJ1200to2500_num,  lumi],
#             [f_WJ2500toInf,   WJ2500toInf_xs,     WJ2500toInf_num,   lumi], 
#             ]

# st       = [
#             [f_ST_tW_top,     ST_tW_top_xs,       ST_tW_top_num,     lumi],
#             [f_ST_tW_antitop, ST_tW_antitop_xs,   ST_tW_antitop_num, lumi],
#             [f_ST_t,          ST_t_xs,            ST_t_num,          lumi],
#             [f_ST_t_ex1,      ST_t_xs,            ST_t_ex1_num,      lumi],
#             [f_ST_s,          ST_s_xs,            ST_s_num,          lumi], 
#            ]

vv       = [
            [f_ZZTo2L2Nu,     ZZTo2L2Nu_xs,       ZZTo2L2Nu_num,    lumi],
            [f_WZTo2L2Q,      WZTo2L2Q_xs,        WZTo2L2Q_num,     lumi],
            [f_WWTo2L2Nu,     WWTo2L2Nu_xs,       WWTo2L2Nu_num,    lumi],
#            [f_WZTo3LNu,    WZTo3LNu_xs,     WZTo3LNu_num,    lumi],
##            [f_ZZTo4L,     ZZTo4L_xs,      ZZTo4L_num, lumi],
           ]

#tZtZ_800 = [[f_TpTp_tZtZ_800, TpTp800_xs,         TpTp800_num,       lumi]]
# #tZbW_800 = [[f_TpTp_tZbW_800, TpTp800_xs,         TpTp800_num,       lumi]]
# #tZtH_800 = [[f_TpTp_tZtH_800, TpTp800_xs,         TpTp800_num,       lumi]]
#tZtZ_1000 = [[f_TpTp_tZtZ_1000, TpTp1000_xs,         TpTp1000_num,       lumi]]
# #tZbW_1000 = [[f_TpTp_tZbW_1000, TpTp1000_xs,         TpTp1000_num,       lumi]]
# #tZtH_1000 = [[f_TpTp_tZtH_1000, TpTp1000_xs,         TpTp1000_num,       lumi]]
#tZtZ_1200 = [[f_TpTp_tZtZ_1200, TpTp1200_xs,         TpTp1200_num,       lumi]]
# #tZbW_1200 = [[f_TpTp_tZbW_1200, TpTp1200_xs,         TpTp1200_num,       lumi]]
# #tZtH_1200 = [[f_TpTp_tZtH_1200, TpTp1200_xs,         TpTp1200_num,       lumi]]

bZbZ_800 = [[f_BpBp_bZbZ_800, BpBp800_xs,         BpBp800_bZ_num,       lumi]]
bZbZ_900 = [[f_BpBp_bZbZ_900, BpBp900_xs,         BpBp900_bZ_num,       lumi]]
bZbZ_1000 = [[f_BpBp_bZbZ_1000, BpBp1000_xs,         BpBp1000_bZ_num,       lumi]]
bZbZ_1100 = [[f_BpBp_bZbZ_1100, BpBp1100_xs,         BpBp1100_bZ_num,       lumi]]
bZbZ_1200 = [[f_BpBp_bZbZ_1200, BpBp1200_xs,         BpBp1200_bZ_num,       lumi]]
bZbZ_1300 = [[f_BpBp_bZbZ_1300, BpBp1300_xs,         BpBp1300_bZ_num,       lumi]]
bZbZ_1400 = [[f_BpBp_bZbZ_1400, BpBp1400_xs,         BpBp1400_bZ_num,       lumi]]
bZbZ_1500 = [[f_BpBp_bZbZ_1500, BpBp1500_xs,         BpBp1500_bZ_num,       lumi]]
#bZbZ_1600 = [[f_BpBp_bZbZ_1600, BpBp1600_xs,         BpBp1600_bZ_num,       lumi]]
bZbZ_1700 = [[f_BpBp_bZbZ_1700, BpBp1700_xs,         BpBp1700_bZ_num,       lumi]]
bZbZ_1800 = [[f_BpBp_bZbZ_1800, BpBp1800_xs,         BpBp1800_bZ_num,       lumi]]

bZbH_800 = [[f_BpBp_bZbH_800, BpBp800_xs,         BpBp800_bH_num,       lumi]]
bZbH_900 = [[f_BpBp_bZbH_900, BpBp900_xs,         BpBp900_bH_num,       lumi]]
bZbH_1000 = [[f_BpBp_bZbH_1000, BpBp1000_xs,         BpBp1000_bH_num,       lumi]]
bZbH_1100 = [[f_BpBp_bZbH_1100, BpBp1100_xs,         BpBp1100_bH_num,       lumi]]
bZbH_1200 = [[f_BpBp_bZbH_1200, BpBp1200_xs,         BpBp1200_bH_num,       lumi]]
bZbH_1300 = [[f_BpBp_bZbH_1300, BpBp1300_xs,         BpBp1300_bH_num,       lumi]]
bZbH_1400 = [[f_BpBp_bZbH_1400, BpBp1400_xs,         BpBp1400_bH_num,       lumi]]
bZbH_1500 = [[f_BpBp_bZbH_1500, BpBp1500_xs,         BpBp1500_bH_num,       lumi]]
#bZbH_1600 = [[f_BpBp_bZbH_1600, BpBp1600_xs,         BpBp1600_bH_num,       lumi]]
bZbH_1700 = [[f_BpBp_bZbH_1700, BpBp1700_xs,         BpBp1700_bH_num,       lumi]]
bZbH_1800 = [[f_BpBp_bZbH_1800, BpBp1800_xs,         BpBp1800_bH_num,       lumi]]

h_data     = getHisto(dataLabel,       dataLeg,        pDir, var,  data,     kBlack,     verbose)
h_top      = getHisto(topLabel,        topLeg,         pDir, var,  top,      8,          verbose)
h_dy       = getHisto(dyLabel,         dyLeg,          pDir, var,  dy,       90,         verbose)
#h_wjets    = getHisto(wjLabel,         wjLeg,          pDir, var,  wjets,    kBlue,      verbose)
#h_st       = getHisto(sTLabel,         sTLeg,          pDir, var,  st,       kCyan,      verbose)
h_vv       = getHisto(vvLabel,         vvLeg,          pDir, var,  vv,       kBlue,       verbose)
#h_tZtZ_800 = getHisto('TT_tZtZ_M800_', 'TT_tZtZ_M800', pDir, var,  tZtZ_800, kRed,    verbose)
# # #h_tZbW_800 = getHisto('TT_tZbW_M800_', 'TT_tZbW_M800', pDir, var,  tZbW_800, kGreen+3,  verbose)
#h_tZtH_800 = getHisto('TT_tZtH_M800_', 'TT_tZtH_M800', pDir, var,  tZtH_800, kGreen+2, verbose)

#h_tZtZ_1000 = getHisto('TT_tZtZ_M1000_', 'TT_tZtZ_M1000', pDir, var, tZtZ_1000, kRed+2, verbose)
# # #tZbW_1000 = [[f_TpTp_tZbW_1000, TpTp1000_xs,         TpTp1000_num,       lumi]]
#h_tZtH_1000 =getHisto('TT_tZtH_M1000_', 'TT_tZtH_M1000', pDir, var, tZtH_1000, kYellow+4, verbose)

#h_tZtZ_1200 = getHisto('TT_tZtZ_M1200_', 'TT_tZtZ_M1200', pDir, var,  tZtZ_1200, kRed+2,    verbose)
# # #h_tZbW_1200 = getHisto('TT_tZbW_M1200_', 'TT_tZbW_M1200', pDir, var,  tZbW_1200, kBlue+3,  verbose)
#h_tZtH_1200 = getHisto('TT_tZtH_M1200_', 'TT_tZtH_M1200', pDir, var,  tZtH_1200, kBlue+2, verbose)

#h_bZbZ_800 = getHisto('BB_bZbZ_M800_', 'BB_bZbZ_M800', pDir, var,  bZbZ_800, kCyan,    verbose)
#h_bZbZ_900 = getHisto('BB_bZbZ_M900_', 'BB_bZbZ_M900', pDir, var,  bZbZ_900, kCyan,    verbose)
#h_bZbZ_1000 = getHisto('BB_bZbZ_M1000_', 'BB_bZbZ_M1000', pDir, var,  bZbZ_1000, kCyan+2,    verbose)
#h_bZbZ_1100 = getHisto('BB_bZbZ_M1100_', 'BB_bZbZ_M1100', pDir, var,  bZbZ_1100, kCyan,    verbose)
#h_bZbZ_1200 = getHisto('BB_bZbZ_M1200_', 'BB_bZbZ_M1200', pDir, var,  bZbZ_1200, kCyan+2,    verbose)
#h_bZbZ_1300 = getHisto('BB_bZbZ_M1300_', 'BB_bZbZ_M1300', pDir, var,  bZbZ_1300, kCyan,    verbose)
#h_bZbZ_1400 = getHisto('BB_bZbZ_M1400_', 'BB_bZbZ_M1400', pDir, var,  bZbZ_1400, kCyan,    verbose)
#h_bZbZ_1500 = getHisto('BB_bZbZ_M1500_', 'BB_bZbZ_M1500', pDir, var,  bZbZ_1500, kCyan,    verbose)
##h_bZbZ_1600 = getHisto('BB_bZbZ_M1600_', 'BB_bZbZ_M1600', pDir, var,  bZbZ_1600, kCyan,    verbose)
#h_bZbZ_1700 = getHisto('BB_bZbZ_M1700_', 'BB_bZbZ_M1700', pDir, var,  bZbZ_1700, kCyan,    verbose)
#h_bZbZ_1800 = getHisto('BB_bZbZ_M1800_', 'BB_bZbZ_M1800', pDir, var,  bZbZ_1800, kCyan,    verbose)
##
h_bZbH_800 = getHisto('BB_bZbH_M800_', 'BB_bZbH_M800', pDir, var,  bZbH_800, kRed,    verbose)
#h_bZbH_900 = getHisto('BB_bZbH_M900_', 'BB_bZbH_M900', pDir, var,  bZbH_900, kCyan,    verbose)
#h_bZbH_1000 = getHisto('BB_bZbH_M1000_', 'BB_bZbH_M1000', pDir, var,  bZbH_1000, kCyan+2,    verbose)
#h_bZbH_1100 = getHisto('BB_bZbH_M1100_', 'BB_bZbH_M1100', pDir, var,  bZbH_1100, kCyan,    verbose)
h_bZbH_1200 = getHisto('BB_bZbH_M1200_', 'BB_bZbH_M1200', pDir, var,  bZbH_1200, kRed+2,    verbose)
#h_bZbH_1300 = getHisto('BB_bZbH_M1300_', 'BB_bZbH_M1300', pDir, var,  bZbH_1300, kCyan,    verbose)
#h_bZbH_1400 = getHisto('BB_bZbH_M1400_', 'BB_bZbH_M1400', pDir, var,  bZbH_1400, kCyan,    verbose)
#h_bZbH_1500 = getHisto('BB_bZbH_M1500_', 'BB_bZbH_M1500', pDir, var,  bZbH_1500, kCyan,    verbose)
# #h_bZbH_1600 = getHisto('BB_bZbH_M1600_', 'BB_bZbH_M1600', pDir, var,  bZbH_1600, kCyan,    verbose)
#h_bZbH_1700 = getHisto('BB_bZbH_M1700_', 'BB_bZbH_M1700', pDir, var,  bZbH_1700, kCyan,    verbose)
#h_bZbH_1800 = getHisto('BB_bZbH_M1800_', 'BB_bZbH_M1800', pDir, var,  bZbH_1800, kCyan,    verbose)


#if forLimits:
#    h_tZtH_800 = getHisto('TT_tZtH_M800_', 'TT_tZtH_M800', pDir2, var,  tZtH_800, kGreen+2, verbose)
#    h_tZtH_1000 =getHisto('TT_tZtH_M1000_', 'TT_tZtH_M1000', pDir2, var, tZtH_1000, kYellow+4, verbose)
#    h_tZtH_1200 = getHisto('TT_tZtH_M1200_', 'TT_tZtH_M1200', pDir2, var,  tZtH_1200, kBlue+2, verbose)
#    h_bZbH_800 = getHisto('BB_bZbH_M800_', 'BB_bZbH_M800', pDir2, var,  bZbH_800, kRed, verbose)
#    h_bZbH_1000 = getHisto('BB_bZbH_M1000_', 'BB_bZbH_M1000', pDir2, var,  bZbH_1000, kRed+2,    verbose)
#    h_bZbH_1200 = getHisto('BB_bZbH_M1200_', 'BB_bZbH_M1200', pDir2, var,  bZbH_1200, kRed+4,    verbose)

templates = []
templates.append(h_dy)
templates.append(h_top)
templates.append(h_vv)
#templates.append(h_st)
#templates.append(h_wjets)
#templates.append(h_tZtZ_800)
# # #templates.append(h_tZbW_800)
#templates.append(h_tZtH_800)
#templates.append(h_tZtZ_1000)
 # # #templates.append(h_tZbW_1000)
#templates.append(h_tZtH_1000)
#templates.append(h_tZtZ_1200)
# # #templates.append(h_tZbW_1200)
# #templates.append(h_tZtH_1200)
#templates.append(h_bZbZ_800)
#templates.append(h_bZbZ_900)
#templates.append(h_bZbZ_1000)
#templates.append(h_bZbZ_1100)
#templates.append(h_bZbZ_1200)
#templates.append(h_bZbZ_1300)
#templates.append(h_bZbZ_1400)
#templates.append(h_bZbZ_1500)
##templates.append(h_bZbZ_1600)
#templates.append(h_bZbZ_1700)
#templates.append(h_bZbZ_1800)
####
templates.append(h_bZbH_800)
#templates.append(h_bZbH_900)
#templates.append(h_bZbH_1000)
#templates.append(h_bZbH_1100)
templates.append(h_bZbH_1200)
#templates.append(h_bZbH_1300)
#templates.append(h_bZbH_1400)
#templates.append(h_bZbH_1500)
#### #templates.append(h_bZbH_1600)
#templates.append(h_bZbH_1700)
#templates.append(h_bZbH_1800)

#get background uncertainty
h_bkg = h_top.Clone()
h_bkg.Reset()
h_bkg.SetName("total bkg")
h_bkg.Add(h_dy)
h_bkg.Add(h_top)
h_bkg.Add(h_vv)
#h_bkg.Add(h_wjets)
#h_bkg.Add(h_st)

#histo properties
nBins = h_bkg.GetNbinsX()
bMin = h_bkg.GetBinLowEdge(1)
bMax = h_bkg.GetBinLowEdge(nBins+1)
bin1 = h_bkg.GetXaxis().FindBin(bMin)
bin2 = h_bkg.GetXaxis().FindBin(bMax)

for ibin in range(0,nBins+1):    
    iTop     = h_top.GetBinContent(ibin)
    iDY      = h_dy.GetBinContent(ibin)
#    iWJ      = h_wjets.GetBinContent(ibin)
    iVV      = h_vv.GetBinContent(ibin)
    # stat error
    stat_err = (h_bkg.GetBinError(ibin))**2 
    # add approximate systematic uncertainty to each bin
    lumi_err = 0.026**2
    if pDir == 'ana/pre':
        btag_err = 0
    else:
        btag_err = 0.017**2

    ID_err   = 0.03**2
    JES_err  = 0.05**2
    #JER_err = .009*.009
    #Pileup_err = .07*.07
    
    dy_err   = (0.15*iDY)**2
    top_err  = (0.15*iTop)**2
#    st_err   = (0.3*iTop)**2
#    wjet_err = (0.1*iWJ)**2
    vv_err   = (0.2*iVV)**2

    new_err = stat_err + lumi_err + btag_err + ID_err + JES_err + dy_err + top_err + vv_err
    #new_err = stat_err + lumi_err + ID_err + JES_err + JER_err + Pileup_err + dy_err + top_err + vv_err# + wjet_err +st_err

    if h_bkg.GetBinError(ibin) != 0: h_bkg.SetBinError(ibin, TMath.Sqrt(new_err))

h_bkg.SetMarkerSize(0)
h_bkg.SetLineWidth(2)
h_bkg.SetFillColor(14)
h_bkg.SetLineColor(0)
h_bkg.SetFillStyle(3244)

#histogram to print the total background with stat uncertainty
h_tot = h_top.Clone()
h_tot.Reset()
h_tot.SetName("Total_"+h_tot.GetName().split('_',1)[1])
h_tot.Add(h_top)
h_tot.Add(h_dy)
h_tot.Add(h_vv)
#h_tot.Add(h_wjets)
#h_tot.Add(h_st)

print h_tot.GetName().split('_',1)[1]

if var == 'chi_mass_cnt':
    overUnderFlow(h_data)

#=========Drawing==============
integralError = Double(5)
#print the latex table:
print '\\begin{tabular}{|c|c| }'
print '\hline'
print 'Sample     & Events  \\\\ '
print '\hline'
count = 0

#f = TFile(plotDir+"/"+var+".root", "RECREATE")

if 'res' in var:
    suffix = 'Res'
elif 'boost' in var:
    suffix = 'Boost'
elif 'merge' in var:
    suffix = 'Merge'
elif 'combo' in var:
    suffix = 'Combo'

if var == 'resST':
    num = 6
elif var == 'boostST':
    num = 8
elif var == 'mergeST':
    num = 8
elif var == 'comboST':
    num = 8
elif var == 'resReco':
    num = 8
elif var == 'boostReco':
    num = 10
elif var == 'mergeReco':
    num = 10
elif var == 'comboReco':
    num = 10
elif var == 'resReco_bZ':
    num = 11
elif var == 'boostReco_bZ':
    num = 13
elif var == 'mergeReco_bZ':
    num = 13
elif var == 'comboReco_bZ':
    num = 13
elif var == 'resReco_bH':
    num = 11
elif var == 'boostReco_bH':
    num = 13
elif var == 'mergeReco_bH':
    num = 13
elif var == 'comboReco_bH':
    num = 13
else:
		num = 0

for ihist in templates :
    #if var != 'cutflow':
    #    overUnderFlow(ihist)
    if var == 'eventCount':
        print ihist.GetName()
        for i in range(1,4):
            print '\n\n\n\ region '+str(i)+' eventCount: '+str(ihist.GetBinContent(i))
    count = count+1
    if count == 3:
        print '\hline'
    if count == 13:
        print '\hline'
    ihist.IntegralAndError(bin1,bin2,integralError)
    if 'TT' in ihist.GetName() or 'BB' in ihist.GetName():
        print '{0:<5} & {1:<5.2f} $\pm$ {2:<5.2f} \\\\ '.format(ihist.GetName().split('_')[1]+ihist.GetName().split('_')[2], ihist.Integral(bin1,bin2), integralError)
        n=ihist.GetName()[:-num]
    else:      
        print '{0:<5} & {1:<5.2f} $\pm$ {2:<5.2f} \\\\ '.format(ihist.GetName().split('_')[0], ihist.Integral(bin1,bin2), integralError)
        n=ihist.GetName().split('_')[0]

    #ihist.SetName("diel_"+suffix+"__"+n)
    #ihist.Write()
h_tot.IntegralAndError(bin1, bin2, integralError)
print '\hline'
print '{0:<5} & {1:<5.2f} $\pm$ {2:<5.2f}\\\\ '.format('Tot Bkg', h_tot.Integral(bin1,bin2), integralError)
print '\hline'
print '{0:<5} & {1:<5.0f} \\\\ '.format(h_data.GetName().split('_')[0], h_data.Integral())
print '\end{tabular}'
#print 'bkg : ', h_bkg.Integral(ibin,bin2), 'tot : ', h_tot.Integral(ibin,bin2)

#h_data.SetName("diel_"+suffix+"__DATA")

#h_data.Write()
#f.Close()

hs = THStack("","")

for ihist in reversed(templates[0:3]):
    #ihist.GetXaxis().SetRangeUser(0, 600)
    hs.Add(ihist)
#    print 'histo added', ihist.GetName()

# Canvas
c1 = TCanvas('c1', 'c1', 800, 600)
c1.Divide(1,2)
scale = (1.0 - 0.3)/0.35

# prepare top pad for original plot
pad = c1.cd(1)
pad.SetPad(0, 0.3, 1, 1)
pad.SetTopMargin(0.1)
pad.SetBottomMargin(0.005)
t = pad.GetTopMargin()

# prepare the 2nd pad
pad = c1.cd(2)
pad.SetPad(0, 0.0, 1, 0.3)
pad.SetTopMargin(0.06)
pad.SetBottomMargin(0.4)
pad.SetTickx(1)
pad.SetTicky(1)
c1.cd(1)

t = c1.GetTopMargin()

maxh=0
for temp in templates:
    if temp.GetMaximum() > maxh:
        maxh = temp.GetMaximum()



if h_data.GetMaximum() > hs.GetMaximum():
   hs.SetMaximum(h_data.GetMaximum())
else:
   h_data.SetMaximum(hs.GetMaximum())

if maxh > hs.GetMaximum() and maxh > h_data.GetMaximum():
    hs.SetMaximum(maxh)
    h_data.SetMaximum(maxh)


hs.SetMinimum(0.1)
#gPad.SetLogy()

if var == 'cutflow':
    data_nbins = h_data.GetNbinsX()
#    print data_nbins
    for a in range(0, 1):
        h_data.SetBinContent(data_nbins-a, -1)

#h_data.GetXaxis().SetRangeUser(0, 400)
#hs.Draw()
#hs.GetXaxis().SetRangeUser(0, 400)
hs.Draw("Hist")
h_bkg.Draw("e2 same")
h_data.Draw("same")

for ihist in reversed(templates[3:]):
#    print 'overlaying, ', ihist.GetName() 
    #ihist.GetXaxis().SetRangeUser(0, 600)
    ihist.Draw("ehist same")

xTitle= h_top.GetXaxis().GetTitle()
yTitle= h_top.GetYaxis().GetTitle()

setTitle(hs, xTitle)
gPad.RedrawAxis()
ll = TLatex()
ll.SetNDC(kTRUE)
ll.SetTextSize(0.05)
ll.DrawLatex(0.78,0.92, "35.9 fb^{-1} (13 TeV)");#2.2

cms = TLatex()
cms.SetNDC(kTRUE)
cms.SetTextFont(61)
cms.SetTextSize(0.08)
cms.DrawLatex(0.12, 1-t+0.2*t,"CMS")

sel = TLatex()
sel.SetNDC(kTRUE)
sel.SetTextSize(0.065)

chan = TLatex()
chan.SetNDC(kTRUE)
chan.SetTextSize(0.065)
chan.DrawLatex(0.50, 0.76, title)

prel = TLatex()
prel.SetNDC(kTRUE)
prel.SetTextFont(52)
prel.SetTextSize(0.75*t*0.76)
prel.DrawLatex(0.22,0.92,"Preliminary")

leg.Draw()
gPad.RedrawAxis()


c1.cd(2)
# add the systematic band
h_ratio = h_data.Clone()
h_ratio_bkg = h_bkg.Clone()
h_ratio_bkg.SetDirectory(0)
h_ratio.SetDirectory(0)
h_ratio.Divide(h_data, h_tot)
h_ratio_bkg.Divide(h_bkg, h_tot)

for ibin in range(1, nBins+1):
    if h_bkg.GetBinContent(ibin) == 0: h_ratio_bkg.SetBinContent(ibin,1)

prepareRatio(h_ratio, h_ratio_bkg, scale, xTitle)

fit = TF1("f1", "pol1", 200, 4000)
line = TLine(bMin, 1, bMax, 1)
line.SetLineColor(kBlack)
h_ratio.Draw("")
h_ratio_bkg.Draw("e2same")
#h_ratio.Fit("f1", "R")
h_ratio.Draw("same")
#h_ratio_bkg.Draw("")
line.Draw()

gPad.RedrawAxis()

#create a directory if it doesn't exist
m_1 = 'mkdir '+plotDir
m_2 = 'mkdir '+plotDir+"/"+skimType
if not os.path.isdir(plotDir):
    subprocess.call( [m_1], shell=True )
if not os.path.isdir(plotDir+"/"+skimType):
    subprocess.call( [m_2], shell=True )    
    
if options.drawEff:
    execfile("drawEff.py")
else:
    c1.SaveAs(plotDir+"/"+skimType+"/"+var+"_.png")
    c1.SaveAs(plotDir+"/"+skimType+"/"+var+"_.pdf")
#raw_input("hold on")
