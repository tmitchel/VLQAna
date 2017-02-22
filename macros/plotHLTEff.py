#!/usr/bin/env python

usage='''
Usage: 
python plotHLTEff.py singleT_Trig.root HT
python plotHLTEff.py singleT_Trig.root Mjj

The first argument is the filename, the second argument is the histogram name used to plot the efficiency
'''

import sys, ROOT
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetNdivisions(405,"x");
ROOT.gStyle.SetEndErrorSize(0.)
ROOT.gStyle.SetErrorX(0.001)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
  
def DivideGraphs(g1,g2):
  if not g1.IsA().InheritsFrom("TGraph") or not g2.IsA().InheritsFrom("TGraph"):
    print "!!!Error: Inputs require two TGraph objects!!!"
    return
  if g1.GetN() != g2.GetN():
    print "!!!Error: The two TGraphs have different numbers of points!!!"
    return 

  n = g1.GetN()
  mjj = []
  ratio = []
  for i in range(0,n):
    x = ROOT.Double(0)
    y1 = ROOT.Double(0)
    y2 = ROOT.Double(0)
    g1.GetPoint(i,x,y1)
    g2.GetPoint(i,x,y2)
    r = y1/y2
    mjj.append(x)
    ratio.append(r)
  gratio = ROOT.TGraph(n, np.asarray(mjj), np.asarray(ratio))
  gratio.SetName("gratio")
  return gratio

def compareHLTEff():

  var = sys.argv[1]
  run = sys.argv[2]

  fjetht = ROOT.TFile("JetHT_BtoH_trigeff_mjj.root","READ")
  fsinglemu = ROOT.TFile("SingleMu_BtoH_trigeff_mjj.root","READ")

  effjetht = fjetht.Get("divide_mjj_by_mjj")
  effjetht.SetName("JetHT")
  effjetht.SetTitle("JetHT: RunsB-H")

  effsinglemu = fsinglemu.Get("divide_mjj_by_mjj")
  effsinglemu.SetName("SingleMuon")
  effsinglemu.SetTitle("SingleMuon: RunsB-H")

  effjetht.SetLineColor(1)
  effjetht.SetMarkerColor(1)
  effjetht.SetMarkerStyle(20)

  effsinglemu.SetLineColor(2)
  effsinglemu.SetMarkerColor(2)
  effsinglemu.SetMarkerStyle(21)

  gratio = DivideGraphs(effsinglemu,effjetht)
  gratio.SetTitle("SingleMu/ JetHT")
  gratio.SetMarkerColor(1)
  gratio.SetMarkerStyle(20)

  c0 = ROOT.TCanvas("JetHT_SingleMuon_BtoH_"+var+"_trigEff","",800,600)
  c0.cd()
  pad0 = ROOT.TPad("pad0","",0.0,0.3,1.0,1.0)
  pad1 = ROOT.TPad("pad1","",0.0,0.0,1.0,0.3)

  pad0.Draw()
  pad0.cd()
  pad0.SetBottomMargin(0.05)

  effjetht.Draw("aep")
  effsinglemu.Draw("ep")

  h = effjetht.GetHistogram()
  h.GetXaxis().SetRangeUser(700,3000)
  h.GetXaxis().SetTitleSize(0)
  h.GetXaxis().SetLabelSize(0)

  leg = ROOT.TLegend(0.5,0.5,0.88,0.78,"","brNDC")
  leg.SetBorderSize(0)
  leg.SetFillColor(0)
  leg.SetTextSize(0.050)
  leg.SetHeader(run)
  leg.AddEntry(effsinglemu,effsinglemu.GetTitle(),"lpe")
  leg.AddEntry(effjetht   ,effjetht.GetTitle(),"lpe")
  leg.Draw()

  pad0.Update()
  c0.cd()

  pad1.Draw()
  pad1.cd()
  pad1.SetTopMargin(0.05)
  pad1.SetBottomMargin(0.3)
  gratio.Draw("ape")
  gratio.SetMinimum(0.6)
  gratio.SetMaximum(1.6)
  gratio.GetHistogram().GetXaxis().SetRangeUser(700,3000)
  gratio.GetHistogram().GetXaxis().SetTitle(h.GetXaxis().GetTitle())
  gratio.GetHistogram().GetYaxis().SetTitle(gratio.GetTitle())
  gratio.GetXaxis().SetTitleSize(0.12)
  gratio.GetXaxis().SetLabelSize(0.12)
  gratio.GetYaxis().SetTitleSize(0.12)
  gratio.GetYaxis().SetLabelSize(0.12)
  gratio.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset()*2.5/7.)
  gratio.GetYaxis().SetNdivisions(505)

  pad1.Update()

  c0.cd()
  c0.SaveAs(c0.GetName()+".pdf")
  c0.SaveAs(c0.GetName()+".png")

def plotHLTEff () :

  if len(sys.argv) < 3:
    sys.exit(usage)
  
  sys.argv.append('-b')
  
  f = ROOT.TFile(sys.argv[1], "README")
  
  hden = f.Get("/hh4b/"+sys.argv[2])
  hnum = f.Get("/hh4btest/"+sys.argv[2])
  
  nbins = int(sys.argv[3])
  
  run = sys.argv[4]
  
  hden.Rebin(nbins)
  hnum.Rebin(nbins)
  
  c0 = ROOT.TCanvas("JetHT_"+run+"_"+sys.argv[2]+"_all_pass","",800,600)
  c0.cd()
  c0.SetLogy()
  hden.SetLineColor(2)
  hnum.SetLineColor(3)
  hden.GetXaxis().SetRangeUser(500,3000)
  hden.GetYaxis().SetTitle("Events/"+str(int(hden.GetBinWidth(1)))+" GeV")
  hden.Draw("hist")
  hnum.Draw("histsame")
  c0.Update()
  leg = ROOT.TLegend(0.5,0.6,0.88,0.88,"","brNDC")
  leg.SetBorderSize(0)
  leg.SetFillColor(0)
  leg.SetTextSize(0.020)
  leg.SetHeader(run)
  leg.AddEntry(hden,"PFJet260 (prescaled)","l")
  leg.AddEntry(hnum,"PFJet260 (prescaled) && OR combination","l")
  leg.Draw()
  c0.Update()
  txt0 = ROOT.TPaveText(0.4,0.55,0.88,0.65,"brNDC")
  txt0.SetBorderSize(0)
  txt0.SetFillColor(0)
  txt0.SetTextFont(42);
  txt0.SetTextSize(0.02);
  txt0.AddText("Jet p_{T} > 200 GeV, |#eta| < 2.4")
  txt0.AddText("#Delta#eta| < 1.3 jet mass > 50 GeV")
  txt0.Draw()
  c0.Update()
  c0.RedrawAxis()
  c0.SaveAs(c0.GetName()+".pdf")
  c0.SaveAs(c0.GetName()+".root")
  c0.SaveAs(c0.GetName()+".C")
  
  c = ROOT.TCanvas("JetHT_"+run+"_"+"trigeff_"+sys.argv[2],"",800,600)
  c.cd()
  gr = ROOT.TGraphAsymmErrors(hnum, hden, 'cp')
  gr.SetLineColor(1)
  gr.SetMarkerColor(1)
  gr.SetMarkerStyle(20)
  gr.Draw("AEP")
  gr.GetHistogram().GetXaxis().SetRangeUser(500,3000)
  gr.GetHistogram().GetXaxis().SetTitle(hden.GetXaxis().GetTitle())
  gr.GetHistogram().GetYaxis().SetTitle("Efficiency")
  c.Update()
  leg = ROOT.TLegend(0.3,0.5,0.88,0.70,"","brNDC")
  leg.SetBorderSize(0)
  leg.SetFillColor(0)
  leg.SetTextSize(0.020)
  leg.SetHeader(run)
  leg.AddEntry(gr,"#splitline{PFHT800 OR AK8PFHT650_TrimR0p1PT0p03Mass50}{#splitline{ OR AK8PFJet360 OR PFHT650_WideJetMJJ900DEtaJJ1p5}{#splitline{AK8PFHT700_TrimR0p1PT0p03Mass50}{ OR AK8DiPFJet280_200_TrimMass30_BTagCSV_p20}}}","l")
  leg.Draw()
  c.Update()
  txt0 = ROOT.TPaveText(0.4,0.35,0.88,0.45,"brNDC")
  txt0.SetBorderSize(0)
  txt0.SetFillColor(0)
  txt0.SetTextFont(42);
  txt0.SetTextSize(0.02);
  txt0.AddText("Jet p_{T} > 200 GeV, |#eta| < 2.4")
  txt0.AddText("#Delta#eta| < 1.3 jet mass > 50 GeV")
  txt0.Draw()
  c.Update()
  c.RedrawAxis()
  c.SaveAs(c.GetName()+".pdf")
  
  fout = ROOT.TFile(c.GetName()+".root","RECREATE")
  fout.cd()
  gr.Write()
  fout.Close()

compareHLTEff()
