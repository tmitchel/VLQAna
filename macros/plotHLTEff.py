#!/usr/bin/env python

import os, sys, ROOT
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
  
def DivideGraphs(g1,g2): ### Returns g1/g2
  if not g1.IsA().InheritsFrom("TGraph") or not g2.IsA().InheritsFrom("TGraph"):
    print "!!!Error: Inputs require two TGraph objects!!!"
    return
  ##if g1.GetN() != g2.GetN():
  ##  print "!!!Error: The two TGraphs have different numbers of points!!! %i and %i" % (g1.GetN(), g2.GetN())
  ##  return 

  n1 = g1.GetN()
  n2 = g2.GetN()
  mjj = []
  ratio = []

  for i2 in range(0,n2):
    x2 = ROOT.Double(0)
    y2 = ROOT.Double(0)
    g2.GetPoint(i2,x2,y2)

    for i1 in range(0,n1):
      x1 = ROOT.Double(0)
      y1 = ROOT.Double(0)
      g1.GetPoint(i1,x1,y1)

      if x1 == x2 and y2 != 0:
        r = y1/y2
        mjj.append(x2)
        ratio.append(r)

  gratio = ROOT.TGraph(len(mjj), np.asarray(mjj), np.asarray(ratio))
  gratio.SetName("gratio")
  return gratio

def compareHLTEff(num, den, samplenum, sampleden, tagnum, tagden, etabin):

  var = sys.argv[1]
  run = sys.argv[2]

  if etabin == "py0": etalabel = "0.0 < |#eta| < 0.434"
  elif etabin == "py1": etalabel = "0.434 < |#eta| < 0.868"
  elif etabin == "py2": etalabel = "0.868 < |#eta| < 1.3"
  else: etalabel = "0.0 < |#eta| < 1.3"

  fjetht = ROOT.TFile(num,"READ")
  fqcd = ROOT.TFile(den,"READ")

  if etabin != "": effjetht = fjetht.Get("divide_"+samplenum+"__mjjred_deta_MJSel_num_"+etabin+"_by_"+samplenum+"__mjjred_deta_MJSel_den_"+etabin)
  else: effjetht = fjetht.Get("divide_"+samplenum+"__mjjred_MJSel_by_"+samplenum+"__mjjred_MJSel")
  effjetht.SetName(tagnum)
  effjetht.SetTitle(tagnum)

  if etabin != "": effqcd = fqcd.Get("divide_"+sampleden+"__mjjred_deta_MJSel_num_"+etabin+"_by_"+sampleden+"__mjjred_deta_MJSel_den_"+etabin)
  else: effqcd = fqcd.Get("divide_"+sampleden+"__mjjred_MJSel_by_"+sampleden+"__mjjred_MJSel")
  effqcd.SetName(tagden)
  effqcd.SetTitle(tagden)

  effjetht.SetLineColor(1)
  effjetht.SetMarkerColor(1)
  effjetht.SetMarkerStyle(20)

  effqcd.SetLineColor(2)
  effqcd.SetMarkerColor(2)
  effqcd.SetMarkerStyle(21)

  gratio = DivideGraphs(effjetht, effqcd)
  gratio.SetTitle(tagnum+"/"+tagden)
  gratio.SetMarkerColor(1)
  gratio.SetMarkerStyle(20)

  c0 = ROOT.TCanvas(tagnum+"_"+tagden+"_"+var+"_"+etabin+"_trigEff","",800,600)
  c0.cd()
  pad0 = ROOT.TPad("pad0","",0.0,0.3,1.0,1.0)
  pad1 = ROOT.TPad("pad1","",0.0,0.0,1.0,0.3)

  pad0.Draw()
  pad0.cd()
  pad0.SetBottomMargin(0.05)

  effjetht.Draw("aep")
  effqcd.Draw("ep")

  h = effjetht.GetHistogram()
  h.GetXaxis().SetRangeUser(700,3000)
  h.GetXaxis().SetTitleSize(0)
  h.GetXaxis().SetLabelSize(0)
  h.GetYaxis().SetTitle("Trigger efficiency")

  leg = ROOT.TLegend(0.5,0.5,0.88,0.78,"","brNDC")
  leg.SetBorderSize(0)
  leg.SetFillColor(0)
  leg.SetTextSize(0.050)
  leg.SetHeader(etalabel)
  leg.AddEntry(effqcd,effqcd.GetTitle(),"lpe")
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
  gratio.GetHistogram().GetXaxis().SetTitle("M_{red}(jj) [GeV]")
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

  fout = ROOT.TFile(c0.GetName()+".root","RECREATE")
  fout.cd()
  gratio.Write()
  fout.Close()

def plotHLTEff (den, num, sample, tag) :

  fden = ROOT.TFile(den, "READ")

  fnum = ROOT.TFile(num, "READ")

  hden = fden.Get(sample+"__mjjred_deta_MJSel")
  hnum = fnum.Get(sample+"__mjjred_deta_MJSel")

  hxden = fden.Get(sample+"__mjjred_MJSel")
  hxnum = fnum.Get(sample+"__mjjred_MJSel")

  hxden.Rebin(4)
  hxnum.Rebin(4)

  c0 = ROOT.TCanvas(tag+"_mjjred_MJSel_all_pass","",800,600)
  c0.cd()
  c0.SetLogy()
  hxden.SetLineColor(2)
  hxnum.SetLineColor(3)
  hxden.GetXaxis().SetRangeUser(500,3000)
  hxden.GetYaxis().SetTitle("Events/"+str(int(hden.GetYaxis().GetBinWidth(1)))+" GeV")
  hxden.Draw("hist")
  hxnum.Draw("histsame")
  c0.Update()
  leg = ROOT.TLegend(0.5,0.65,0.88,0.88,"","brNDC")
  leg.SetBorderSize(0)
  leg.SetFillColor(0)
  leg.SetTextSize(0.020)
  leg.SetHeader(tag)
  leg.AddEntry(hxden,"HLT_PFJet260","l")
  leg.AddEntry(hxnum,"HLT_PFJet260 OR combination","l")
  #leg.AddEntry(hden,"MC sel.","l")
  #leg.AddEntry(hnum,"MC sel. and HLT_PFJet260 ","l")
  leg.Draw()
  c0.Update()
  txt0 = ROOT.TPaveText(0.45,0.50,0.88,0.60,"brNDC")
  txt0.SetBorderSize(0)
  txt0.SetFillColor(0)
  txt0.SetTextFont(42);
  txt0.SetTextSize(0.02);
  txt0.AddText("Jet p_{T} > 300 GeV, |#eta| < 2.4")
  txt0.AddText("#Delta#eta| < 1.3 105 < jet SD mass < 135 GeV")
  txt0.Draw()
  c0.Update()
  c0.RedrawAxis()
  c0.SaveAs(c0.GetName()+".pdf")
  c0.SaveAs(c0.GetName()+".root")
  c0.SaveAs(c0.GetName()+".C")

  nx = hden.GetNbinsX()
  ny = hden.GetNbinsY()
  c1 = ROOT.TCanvas(tag+"_mjjred_MJSel_all_trigeff","",800,600)
  fout = ROOT.TFile(c1.GetName()+".root","RECREATE")
  c1.cd()
  h = hden.ProjectionX(hden.GetName()+"_py")
  h.Reset()
  h.GetXaxis().SetRangeUser(750,3000)
  h.GetXaxis().SetTitle(hden.GetXaxis().GetTitle())
  h.GetYaxis().SetTitle("Efficiency")
  h.SetLineColor(0)
  h.SetFillColor(0)
  h.SetMarkerColor(0)
  h.SetMarkerSize(0)
  h.Draw()
  gx = ROOT.TGraphAsymmErrors(hxnum, hxden, 'cp')
  leg = ROOT.TLegend(0.23,0.5,0.88,0.70,"","brNDC")
  gx.SetLineColor(6)
  gx.SetMarkerColor(6)
  gx.SetMarkerStyle(34)
  gx.Draw("EP")
  leg.SetBorderSize(0)
  leg.SetFillColor(10)
  leg.SetTextSize(0.020)
  leg.SetHeader(tag)
  leg.AddEntry(h,"#splitline{PFHT800 OR AK8PFHT650_TrimR0p1PT0p03Mass50}{#splitline{ OR AK8PFJet360 OR PFHT650_WideJetMJJ900DEtaJJ1p5}{#splitline{AK8PFHT700_TrimR0p1PT0p03Mass50}{ OR AK8DiPFJet280_200_TrimMass30_BTagCSV_p20}}}","l")
  leg.AddEntry(h,"","")
  leg.AddEntry(h,"","")
  leg.AddEntry(h,"","")
  leg.AddEntry(gx, "All #eta", "l")
  leg.Draw()
  c1.Update()
  txt0 = ROOT.TPaveText(0.4,0.35,0.88,0.45,"brNDC")
  txt0.SetBorderSize(0)
  txt0.SetFillColor(10)
  txt0.SetTextFont(42);
  txt0.SetTextSize(0.02);
  txt0.AddText("Jet p_{T} > 300 GeV, |#eta| < 2.4")
  txt0.AddText("#Delta#eta| < 1.3 105 < jet SD mass < 135 GeV")
  txt0.Draw()
  c1.Update()
  hnumx=[None]*ny
  hdenx=[None]*ny
  gr=[None]*ny
  print "nbinsy = ", ny
  for j in range(0,ny):
    print " bin = %i low %f high %f" % (j+1, hden.GetYaxis().GetBinLowEdge(j+1), hden.GetYaxis().GetBinUpEdge(j+1))
    hnumx[j] = hnum.ProjectionX(hnum.GetName()+"_num_py%i"%j, j, j+1)
    hdenx[j] = hden.ProjectionX(hden.GetName()+"_den_py%i"%j, j, j+1)
    hnumx[j].Rebin(4)
    hdenx[j].Rebin(4)
    gr[j] = ROOT.TGraphAsymmErrors(hnumx[j], hdenx[j], 'cp')
    #gr[j].SetName("greff_%i"%j)
    #gr = ROOT.TGraphAsymmErrors(hnum, hden, 'cp')
    gr[j].SetLineColor(j+1)
    gr[j].SetMarkerColor(j+1)
    gr[j].SetMarkerStyle(j+1)
    gr[j].Draw("EP")
    c1.Update()
    c1.RedrawAxis()
    c1.Update()
    leg.AddEntry(gr[j], "%f < |#eta| < %f" % (hden.GetYaxis().GetBinLowEdge(j+1), hden.GetYaxis().GetBinUpEdge(j+1)),"l")
    fout.cd()
    hnumx[j].Write()
    hdenx[j].Write()
    gr[j].Write()
  fout.cd()
  gx.Write()
  c1.SaveAs(c1.GetName()+".pdf")
  fout.Close()

plotHLTEff("/afs/cern.ch/user/d/devdatta/afswork/CMSREL/Analysis/CMSSW_8_0_20/src/Analysis/VLQAna/test/HH4bTrig_QCDHT_HLTPFJet260_13Mar2017/Final_histograms_hh4b.root", "/afs/cern.ch/user/d/devdatta/afswork/CMSREL/Analysis/CMSSW_8_0_20/src/Analysis/VLQAna/test/HH4bTrig_QCDHT_HLTPFJet260_13Mar2017/Final_histograms_hh4btest.root", "QCD", "QCD_HLTPFJet260Trig")

plotHLTEff("/afs/cern.ch/user/d/devdatta/afswork/CMSREL/Analysis/CMSSW_8_0_20/src/Analysis/VLQAna/test/HH4bTrig_JetHT_QCDHT_JetTrig_13Mar2017/Final_histograms_hh4b.root", "/afs/cern.ch/user/d/devdatta/afswork/CMSREL/Analysis/CMSSW_8_0_20/src/Analysis/VLQAna/test/HH4bTrig_JetHT_QCDHT_JetTrig_13Mar2017/Final_histograms_hh4btest.root", "QCD", "QCD_JetTrig")

plotHLTEff("/afs/cern.ch/user/d/devdatta/afswork/CMSREL/Analysis/CMSSW_8_0_20/src/Analysis/VLQAna/test/HH4bTrig_JetHT_QCDHT_JetTrig_13Mar2017/Final_histograms_hh4b.root", "/afs/cern.ch/user/d/devdatta/afswork/CMSREL/Analysis/CMSSW_8_0_20/src/Analysis/VLQAna/test/HH4bTrig_JetHT_QCDHT_JetTrig_13Mar2017/Final_histograms_hh4btest.root", "DATA", "JetHT_JetTrig")

compareHLTEff("JetHT_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_JetTrig_mjjred_MJSel_all_trigeff.root", "DATA", "QCD", "JetHT", "QCD","")
compareHLTEff("QCD_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_HLTPFJet260Trig_mjjred_MJSel_all_trigeff.root", "QCD", "QCD", "QCDJetTrigs", "QCDHLTPFJet260","")
compareHLTEff("JetHT_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_HLTPFJet260Trig_mjjred_MJSel_all_trigeff.root", "DATA", "QCD", "JetHT", "QCDHLTPFJet260","")

compareHLTEff("JetHT_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_JetTrig_mjjred_MJSel_all_trigeff.root", "DATA", "QCD", "JetHT", "QCD", "py0")
compareHLTEff("QCD_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_HLTPFJet260Trig_mjjred_MJSel_all_trigeff.root", "QCD", "QCD", "QCDJetTrigs", "QCDHLTPFJet260", "py0")
compareHLTEff("JetHT_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_HLTPFJet260Trig_mjjred_MJSel_all_trigeff.root", "DATA", "QCD", "JetHT", "QCDHLTPFJet260", "py0")

compareHLTEff("JetHT_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_JetTrig_mjjred_MJSel_all_trigeff.root", "DATA", "QCD", "JetHT", "QCD", "py1")
compareHLTEff("QCD_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_HLTPFJet260Trig_mjjred_MJSel_all_trigeff.root", "QCD", "QCD", "QCDJetTrigs", "QCDHLTPFJet260", "py1")
compareHLTEff("JetHT_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_HLTPFJet260Trig_mjjred_MJSel_all_trigeff.root", "DATA", "QCD", "JetHT", "QCDHLTPFJet260", "py1")

compareHLTEff("JetHT_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_JetTrig_mjjred_MJSel_all_trigeff.root", "DATA", "QCD", "JetHT", "QCD", "py2")
compareHLTEff("QCD_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_HLTPFJet260Trig_mjjred_MJSel_all_trigeff.root", "QCD", "QCD", "QCDJetTrigs", "QCDHLTPFJet260", "py2")
compareHLTEff("JetHT_JetTrig_mjjred_MJSel_all_trigeff.root", "QCD_HLTPFJet260Trig_mjjred_MJSel_all_trigeff.root", "DATA", "QCD", "JetHT", "QCDHLTPFJet260", "py2")
