#! /usr/bin/env python
import sys
import os
import subprocess
from array import array
from ROOT import TH1D,TH2D,TFile,TMath,TCanvas,THStack,TLegend,TPave,TLine,TLatex
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
                  default= 35900.,#2200.,
                  dest='Lumi',
                  help='Data Luminosity in pb-1')

parser.add_option('--globalSF', metavar='SF', type='float',
                  default = 1.0,
                  dest='globalSF',
                  help='Global trigger SF (%default default)')

parser.add_option('--plotDir', metavar='P', type='string', action='store',
                  default='DecPlots',
                  dest='plotDir',
                  help='output directory of plots')

parser.add_option('--skimType', metavar='S', type='string', action='store',
                  default='CR_Zelel',
                  dest='skimType',
                  help='Skim type: CR_Zelel, CR_Zmumu, SR_Zelel, SR_Zmumu')

parser.add_option('--processDir', metavar='pD', type='string', action='store',
                  default='ana/cnt',
                  dest='processDir',
                  help='directory to read histograms from')

parser.add_option('--var', metavar='T', type='string', action='store',
                  default='st',#cutflow, st
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
dataDir = 'ana/'+pDir.split('/')[-1]

if 'elel' in skimType: title = 'e^{#pm}e^{#mp}+jets'
elif 'mumu' in skimType: title = '#mu^{#pm}#mu^{#mp}+jets'
else: title = ''
# ==========add the input ============

execfile("root_input.py")

# === prepare the input labels and legends ===========
dataLabel     = 'Data_'
topLabel      = 'Top_'
ttZLabel      = 'ttZ_'
dyLabel       = 'DY_'
wjLabel       = 'WJets_'
sTLabel       = 'ST_'
vvLabel       = 'VV_'

dataLeg       = 'Data'
topLeg        = 't#bar{t}'
ttZLeg        = 'ttZ'
dyLeg         = 'Drell-Yan'
wjLeg         = 'W+Jets'
sTLeg         = 'Single top'
vvLeg         = 'other backgrounds'
# === create structure ============
data     = [
            #[f_Data_Oct2015, 1., 1., 1.], # this corresponds to .53 fb-1
            [f_Data_PromptReco, 1., 1., 1.],  #addign this should give 2.2 fb-1  
           ]

top      = [[f_ttbar,         Top_xs,            Top_num,            lumi]]

ttZ      = [[f_ttZ,           ttZ_xs,            ttZ_num,            lumi]]

#dy       = [[f_DYmcnlo,       DY_xs,             DYmcnlo_num,        lumi]]

dy       = [
            [f_DY100to200,    DY100to200_xs,      DY100to200_num,    lumi],
            [f_DY200to400,    DY200to400_xs,      DY200to400_num,    lumi],
            [f_DY400to600,    DY400to600_xs,      DY400to600_num,    lumi],
            [f_DY600toInf,    DY600toInf_xs,      DY600toInf_num,    lumi],    
           ]

#wjets    = [
#            [f_WJ100to200,    WJ100to200_xs,      WJ100to200_num,    lumi],
#            [f_WJ200to400,    WJ200to400_xs,      WJ200to400_num,    lumi],
#            [f_WJ400to600,    WJ400to600_xs,      WJ400to600_num,    lumi],
#            [f_WJ600to800,    WJ600to800_xs,      WJ600to800_num,    lumi],     
#            [f_WJ800to1200,   WJ800to1200_xs,     WJ800to1200_num,   lumi],
#            [f_WJ1200to2500,  WJ1200to2500_xs,    WJ1200to2500_num,  lumi],
#            [f_WJ2500toInf,   WJ2500toInf_xs,     WJ2500toInf_num,   lumi], 
#            ]
#
#st       = [
#            [f_ST_tW_top,     ST_tW_top_xs,       ST_tW_top_num,     lumi],
#            [f_ST_tW_antitop, ST_tW_antitop_xs,   ST_tW_antitop_num, lumi],
#            [f_ST_t,          ST_t_xs,            ST_t_num,          lumi],
#            [f_ST_t_ex1,      ST_t_xs,            ST_t_ex1_num,      lumi],
#            [f_ST_s,          ST_s_xs,            ST_s_num,          lumi], 
#           ]

vv       = [
            [f_ZZTo2L2Nu,     ZZTo2L2Nu_xs,       ZZTo2L2Nu_num,    lumi],
            [f_WZTo2L2Q,      WZTo2L2Q_xs,        WZTo2L2Q_num,     lumi],
            [f_WWTo2L2Nu,     WWTo2L2Nu_xs,       WWTo2L2Nu_num,    lumi],
            [f_ttW,           ttW_xs,            ttW_num,            lumi],
            [f_tZq,           tZq_xs,            tZq_num,            lumi],
            [f_ST_twll,           twll_xs,            twll_num,            lumi],
           ]

bZbZ_800  = [[f_BpBp_bZbZ_800,  BpBp800_xs,  BpBp800_bZ_num,  lumi]]
bZbZ_1000 = [[f_BpBp_bZbZ_1000, BpBp1000_xs, BpBp1000_bZ_num, lumi]]
bZbZ_1200 = [[f_BpBp_bZbZ_1200, BpBp1200_xs, BpBp1200_bZ_num, lumi]]
bZbH_800  = [[f_BpBp_bZbH_800,  BpBp800_xs,  BpBp800_bZbH_num,  lumi]]
bZbH_1000 = [[f_BpBp_bZbH_1000, BpBp1000_xs, BpBp1000_bH_num, lumi]]
bZbH_1200 = [[f_BpBp_bZbH_1200, BpBp1200_xs, BpBp1200_bZbH_num, lumi]]


h_data     = getHisto(dataLabel,       dataLeg,        dataDir, var,  data,     kBlack,     verbose)
h_top      = getHisto(topLabel,        topLeg,         pDir, var,  top,      8,          verbose)
h_ttZ      = getHisto(ttZLabel,        ttZLeg,         pDir, var,  ttZ,      kMagenta+2,          verbose)
h_dy       = getHisto(dyLabel,         dyLeg,          pDir, var,  dy,       90,         verbose)
#h_wjets    = getHisto(wjLabel,         wjLeg,          pDir, var,  wjets,    kBlue,      verbose)
#h_st       = getHisto(sTLabel,         sTLeg,          pDir, var,  st,       kCyan,      verbose)
h_vv       = getHisto(vvLabel,         vvLeg,          pDir, var,  vv,       kBlue,       verbose)

h_bZbZ_800  = getHisto('BB_bZbZ_M800_', 'BB_bZbZ_M800', pDir, var, bZbZ_800, kCyan, verbose)
#h_bZbZ_1000  = getHisto('BB_bZbZ_M1000_', 'BB_bZbZ_M1000', pDir, var, bZbZ_1000, kCyan, verbose)
h_bZbZ_1200  = getHisto('BB_bZbZ_M1200_', 'BB_bZbZ_M1200', pDir, var, bZbZ_1200, kCyan+2, verbose)
#h_bZbH_800  = getHisto('BB_bZbH_M800_', 'BB_bZbH_M800', pDir, var, bZbH_800, kRed, verbose)
#h_bZbH_1000  = getHisto('BB_bZbH_M1000_', 'BB_bZbH_M1000', pDir, var, bZbH_1000, kCyan, verbose)
#h_bZbH_1200  = getHisto('BB_bZbH_M1200_', 'BB_bZbH_M1200', pDir, var, bZbH_1200, kRed+2, verbose)


templates = []
templates.append(h_dy)
templates.append(h_top)
templates.append(h_ttZ)
templates.append(h_vv)
#templates.append(h_st)
#templates.append(h_wjets)
templates.append(h_bZbZ_800)
#templates.append(h_bZbZ_1000)
templates.append(h_bZbZ_1200)
#templates.append(h_bZbH_800)
#templates.append(h_bZbH_1000)
#templates.append(h_bZbH_1200)

#f = TFile(plotDir+"/"+skimType+"/"+var+".root", "RECREATE")
#for ihist in templates:
#    ihist.Write()
#f.Close()

#get background uncertainty
h_bkg = h_top.Clone()
h_bkg.Reset()
h_bkg.SetName("total bkg")
h_bkg.Add(h_dy)
h_bkg.Add(h_top)
h_bkg.Add(h_ttZ)
h_bkg.Add(h_vv)

#histo properties
nBins = h_bkg.GetNbinsX()
bMin = h_bkg.GetBinLowEdge(1)
bMax = h_bkg.GetBinLowEdge(nBins+1)
bin1 = h_bkg.GetXaxis().FindBin(bMin)
bin2 = h_bkg.GetXaxis().FindBin(bMax)

for ibin in range(0,nBins+1):    
    iTop     = h_top.GetBinContent(ibin)
    ittZ     = h_ttZ.GetBinContent(ibin)
    iDY      = h_dy.GetBinContent(ibin)
    iVV      = h_vv.GetBinContent(ibin)
    # stat error
    stat_err = (h_bkg.GetBinError(ibin))**2 
    # add approximate systematic uncertainty to each bin
    lumi_err = 0.026**2
    btag_err = 0.017**2
    ID_err   = 0.03**2
    JES_err  = 0.05*0.05
    dy_err   = (0.15*iDY)**2
    top_err  = (0.15*iTop)**2
    vv_err   = (0.15*iVV)**2
    ttZ_err = (0.15*ittZ)**2
    new_err = stat_err  
#    new_err = stat_err + lumi_err + btag_err + ID_err + JES_err + dy_err + top_err + vv_err
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

print h_tot.GetName().split('_',1)[1]

## =========Drawing==============
integralError = Double(5)
# print the latex table:
print '\\begin{tabular}{|c|c| }'
print '\hline'
print 'Sample     & Events  \\\\ '
print '\hline'
count = 0
for ihist in templates :
    #if var != 'cutFlow':overUnderFlow(ihist)
    count = count+1
    if count == 10:
        print '\hline'
    if count == 14:
        print '\hline'
    ihist.IntegralAndError(bin1,bin2,integralError)
    if 'TT' in ihist.GetName() or 'BB' in ihist.GetName():
        print '{0:<5} & {1:<5.2f} $\pm$ {2:<5.2f} \\\\ '.format(ihist.GetName().split('_')[1], ihist.Integral(bin1,bin2), integralError)
    else:      
        print '{0:<5} & {1:<5.2f} $\pm$ {2:<5.2f} \\\\ '.format(ihist.GetName().split('_')[0], ihist.Integral(bin1,bin2), integralError)

totalIntegralError = Double(100)
h_tot.IntegralAndError(bin1, bin2, totalIntegralError)
print '\hline'
print '{0:<5} & {1:<5.0f} $\pm$ {2:<5.2f}\\\\ '.format('Tot Bkg', h_tot.Integral(bin1,bin2), totalIntegralError)
print '\hline'
print '{0:<5} & {1:<5.0f} \\\\ '.format(h_data.GetName().split('_')[0], h_data.Integral())
print '\end{tabular}'
#print 'bkg : ', h_bkg.Integral(ibin,bin2), 'tot : ', h_tot.Integral(ibin,bin2)

hs = THStack("","")

for ihist in reversed(templates[0:4]):
    #ihist.GetXaxis().SetRangeUser(0, 300)
    hs.Add(ihist)
    print 'histo added', ihist.GetName()

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

maxHist = 0.

for hist in reversed(templates[4:]):
  if h_dy.GetMaximum() > hist.GetMaximum():
    maxHist = hs.GetMaximum()
  else:
    maxHist = hist.GetMaximum()
    break

#if h_vv.Integral() == 0:
#  h_vv.SetBinContent(0,0.1)

hs.SetMaximum(maxHist*1.5)
hs.SetMinimum(0.1)
gPad.SetLogy()

#h_data.GetXaxis().SetRangeUser(0, 1)
#h_bkg.GetXaxis().SetRangeUser(0, 1)
hs.Draw("Hist")
#hs.GetXaxis().SetRangeUser(0, 1)
#h_data.SetBinContent(9,0)
#hs.Draw("Hist")
h_bkg.Draw("e2 same")
h_data.Draw("same")

for ihist in reversed(templates[3:]):
    print 'overlaying, ', ihist.GetName() 
    ihist.GetXaxis().SetRangeUser(0, 1)
    ihist.SetFillColor(0)
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

line = TLine(bMin, 1, bMax, 1)
line.SetLineColor(kBlack)
h_ratio.Draw("")
h_ratio_bkg.Draw("e2same")
h_ratio.Draw("same")
line.Draw()

gPad.RedrawAxis()

#create a directory if it doesn't exist
m_1 = 'mkdir '+plotDir
m_2 = 'mkdir '+plotDir+"/"+skimType
if not os.path.isdir(plotDir):
    subprocess.call( [m_1], shell=True )
if not os.path.isdir(plotDir+"/"+skimType):
    subprocess.call( [m_2], shell=True )    
    
c1.SaveAs(plotDir+"/"+skimType+"/"+var+"_.pdf")
#c1.SaveAs(plotDir+"/"+skimType+"/"+var+"_.gif")
#raw_input("hold on")
