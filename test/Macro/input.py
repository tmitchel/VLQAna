#! /usr/bin/env python

# =====================================================
#  INPUTS		
# =====================================================
path = '/uscms_data/d3/tmitchel/update80X/CMSSW_8_0_20/src/Analysis/VLQAna/test/Macro/rootFiles/Zmumu/Oct10_bugFix/'
ch = 'CR_Zmumu'

f_Data_PromptReco = TFile(path+'ttbar.root')

f_DY100to200 = TFile(path+'dy_pt100to250.root')
f_DY200to400 = TFile(path+'dy_pt250to400.root')
f_DY400to600 = TFile(path+'dy_pt400to650.root')
f_DY600toInf = TFile(path+'dy_pt650toInf.root')

f_ZZTo2L2Nu     = TFile(path+'ZZ.root')
f_WZTo2L2Q      = TFile(path+'WZ.root')
f_WWTo2L2Nu     = TFile(path+'WW.root')

f_ttbar         = TFile(path+'ttbar.root')

f_BpBp_bZbZ_800  = TFile(path+'bprime800_bZ.root')
#f_BpBp_bZbZ_1000  = TFile(path+'bprime1000_bZ.root')
f_BpBp_bZbZ_1200 = TFile(path+'bprime1200_bZ.root')
#f_BpBp_bZbH_800  = TFile(path+'bprime800_bH.root')
#f_BpBp_bZbH_1000 = TFile(path+'bprime1000_bH.root')
#f_BpBp_bZbH_1200 = TFile(path+'bprime1200_bH.root')


#===== cross sections (pb)==========

Top_xs            = 831.76  *gSF * 0.9
DY100to200_xs     = 83.12   *gSF 
DY200to400_xs     = 3.047   *gSF 
DY400to600_xs     = 0.3921  *gSF 
DY600toInf_xs     = 0.03636 *gSF 
ZZTo2L2Nu_xs      = 16.91   *gSF
WZTo2L2Q_xs       = 46.74   *gSF
WWTo2L2Nu_xs      = 118.7   *gSF

BpBp800_xs        = 0.196 * gSF
BpBp1000_xs       = 0.044 * gSF
BpBp1200_xs       = 0.0118 * gSF

#===== Number of generated events ======

Top_num          =  f_ttbar.Get("allEvents/hEventCount_nowt").GetBinContent(1)
DY100to200_num   =  f_DY100to200.Get("allEvents/hEventCount_nowt").GetBinContent(1)
DY200to400_num   =  f_DY200to400.Get("allEvents/hEventCount_nowt").GetBinContent(1)
DY400to600_num   =  f_DY400to600.Get("allEvents/hEventCount_nowt").GetBinContent(1)
DY600toInf_num   =  f_DY600toInf.Get("allEvents/hEventCount_nowt").GetBinContent(1)
ZZTo2L2Nu_num    =  f_ZZTo2L2Nu.Get("allEvents/hEventCount_nowt").GetBinContent(1)
WZTo2L2Q_num     =  f_WZTo2L2Q.Get("allEvents/hEventCount_nowt").GetBinContent(1)
WWTo2L2Nu_num    =  f_WWTo2L2Nu.Get("allEvents/hEventCount_nowt").GetBinContent(1)
BpBp800_bZ_num      =  f_BpBp_bZbZ_800.Get('ana/signalEvts').GetBinContent(1)
#BpBp1000_bZ_num     =  f_BpBp_bZbZ_1000.Get('ana/signalEvts').GetBinContent(1)
BpBp1200_bZ_num     =  f_BpBp_bZbZ_1200.Get('ana/signalEvts').GetBinContent(1)
#BpBp800_bH_num      =  f_BpBp_bZbH_800.Get('ana/signalEvts').GetBinContent(1)
#BpBp1000_bH_num     =  f_BpBp_bZbH_1000.Get('ana/signalEvts').GetBinContent(1)
#BpBp1200_bH_num     =  f_BpBp_bZbH_1200.Get('ana/signalEvts').GetBinContent(1)

# Legend
leg = TLegend(0.76,0.88,0.94,0.50)
#leg = TLegend(0.73, 0.92, 0.98, 0.55)
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
    h_ratio.SetMaximum(3)
    h_ratio.SetMinimum(-1)
    h_ratio.GetYaxis().SetLabelSize(0.06*scale)
    h_ratio.GetYaxis().SetTitleOffset(1.00/scale*0.5)
    h_ratio.GetYaxis().SetTitleSize(0.08*scale)
    h_ratio.GetYaxis().SetTitleFont(42)
    h_ratio.GetXaxis().SetLabelSize(0.06*scale)
    h_ratio.GetXaxis().SetTitleOffset(.45*scale)
    h_ratio.GetXaxis().SetTitleSize(0.09*scale)
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
    h_ratio_bkg.GetXaxis().SetLabelSize(0.06*scale)
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
    if 'Data' in hname:
        leg.AddEntry(hist, legname, 'pl')
        hist.SetMarkerStyle(8)
    elif 'bZ' in hname or 'bH' in hname:          
        hist.SetLineWidth(2)
        leg.AddEntry(hist, legname, 'l')
    else:
        print 'hi'
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
#        hist2 = TH1D('cutflow', 'cut flow', 9, 0.5, 9.5)
#        hist2.GetXaxis().SetBinLabel(1, 'Trig.')
#        hist2.GetXaxis().SetBinLabel(2, 'l^{+}l^{-}')
#        hist2.GetXaxis().SetBinLabel(3, '75 < M(l^{+}l^{-} < 105')
#        hist2.GetXaxis().SetBinLabel(4, 'H_{T} #geq 200')
#        hist2.GetXaxis().SetBinLabel(5, 'N (AK4) #geq 3')
#        hist2.GetXaxis().SetBinLabel(6, 'leading jet pt > 100')
#        hist2.GetXaxis().SetBinLabel(7, '2nd jet pt > 50')
#        hist2.GetXaxis().SetBinLabel(8, 'N(bjet) #geq 1')
#        hist2.GetXaxis().SetBinLabel(9, 'S_{T} #geq 1000')
#        for i in range(1, 10):
#          hist2.SetBinContent(i, hist.GetBinContent(i+1))
#        #hist2.GetXaxis().SetLabelSize(0.8)
#        hist2.GetXaxis().SetLabelOffset(.05) 

        if verbose:
            print 'file: {0:<20}, histo:{1:<10}, integral before weight:{2:<3.3f}, nEntries:{3:<3.0f}, weight:{4:<2.3f}'.format(
                ifile.GetName(),    
                hist.GetName(),
                hist.Integral(), hist.GetEntries(), xs * lumi /nevt
                )
        hist.Sumw2()    
        hist.Scale( xs * lumi /nevt)
        hist.Rebin(40)
        histos.append( hist )
        
    histo = histos[0]
    setCosmetics(histo, leg, label+var, color) 
    for ihisto in range(1, len(histos) ):
        #print 'ihisto =', ihisto, 'integral', histos[ihisto].Integral(), ', entries', histos[ihisto].GetEntries()
        histo.Add( histos[ihisto] )
        #print 'after addition', histo.Integral()
    if verbose:    
        print 'newName: {0:<5}, Entries:{1:5.2f},  newIntegral: {2:5.2f}'.format(label+var, histo.GetEntries(), histo.Integral() )   
    return histo

