#!/usr/bin/env python

### PDF4LHC recommendations for LHC Run II: http://arxiv.org/abs/1510.03865

import os, sys, ROOT
import math

def getPDFUncert(fin,cat): 

  idxpdfs=range(10,10+101)
  nwts=len(idxpdfs)

  allevt   = 0
  allevtup = [] 
  allevtlo = [] 

  passevt   = 0
  passevtup = []
  passevtlo = []

  f = ROOT.TFile.Open(fin, "READ")
  mychain = f.Get("hh4b/tree")
  entries=mychain.GetEntriesFast()

  for i in idxpdfs:

    evtwt=0
    evtwtpdfUp=0
    evtwtpdfLo=0
    passevtwt=0
    passevtwtpdfUp=0
    passevtwtpdfLo=0

    evts=50000
    ievt = 0
    for event in mychain:
      ievt = ievt+1
      if ievt > evts: break
    
      lhewts = event.SelectedEvent_lhewts
      lhewtcentral = float(lhewts[idxpdfs[0]].second)

      evtwt+=lhewtcentral

      lhewt = lhewts[i-1]
      if i%2 == 0:
        evtwtpdfUp += float(lhewt.second) 
      elif i%2 == 1:
        evtwtpdfLo += float(lhewt.second)
    
      if ( event.LeadingDiJets_minv_leading2hjets > 1100 and \
          event.LeadingDiJets_deta_leading2hjets < 1.3 and \
          event.AK8Jets_Pt[0] > 300 and event.AK8Jets_Pt[1] > 300 and \
          abs(event.AK8Jets_Eta[0]) < 2.4 and abs(event.AK8Jets_Eta[1]) < 2.4 and \
          (105. < event.AK8Jets_MassSoftDrop[0] < 135.) and (105. < event.AK8Jets_MassSoftDrop[1] < 135.) and \
          (event.AK8Jets_tau2[0]/event.AK8Jets_tau1[0] < 0.55) and (event.AK8Jets_tau2[1]/event.AK8Jets_tau1[1] < 0.55) ) :

        passTT = event.AK8Jets_doublesv[0] > 0.8 and event.AK8Jets_doublesv[1] > 0.8
        passLL = event.AK8Jets_doublesv[0] > 0.3 and event.AK8Jets_doublesv[1] > 0.3

        if (cat=="LL" and passLL and not passTT) or (cat=="TT" and passTT) : 
    
          passevtwt+=lhewtcentral
          if i%2 == 0:
            passevtwtpdfUp += float(lhewt.second) 
          elif i%2 == 1:
            passevtwtpdfLo += float(lhewt.second)

    allevt=evtwt
    allevtup.append(evtwtpdfUp)
    allevtlo.append(evtwtpdfLo)

    passevt=passevtwt
    passevtup.append(passevtwtpdfUp)
    passevtlo.append(passevtwtpdfLo)

  aenominal = passevt/allevt

  allevtup_sum=0
  for k in allevtup:
    allevtup_sum += pow(k - allevt,2.)
  allevtup_sum = math.sqrt(allevtup_sum/100)

  passevtup_sum=0
  for k in passevtup:
    passevtup_sum += pow(k - passevt,2.)
  passevtup_sum = math.sqrt(passevtup_sum/100)

  aepdfup = passevtup_sum/allevtup_sum

  allevtlo_sum=0
  for k in allevtlo:
    allevtlo_sum += pow(k - allevt,2.)
  allevtlo_sum = math.sqrt(allevtlo_sum/100)

  passevtlo_sum=0
  for k in passevtlo:
    passevtlo_sum += pow(k - passevt,2.)
  passevtlo_sum = math.sqrt(passevtlo_sum/100)

  aepdflo = passevtlo_sum/allevtlo_sum

  print " aenominal %f aepdfup %f aepdflo %f" % (aenominal, aepdfup, aepdflo)

  return aepdfup/aenominal, aepdflo/aenominal

from glob import glob

fs = glob("/afs/cern.ch/user/d/devdatta/afswork/CMSREL/Analysis/CMSSW_8_0_20/src/Analysis/VLQAna/test/HH4bTrees_15Mar2017/*/output/hh4b_0.root")

for f in fs:
  mass = f.split("/")[-3]
  print mass
  out = open("pdfuncert%s.log" % mass,"w")
  upLL, downLL = getPDFUncert(f,"LL")
  upTT, downTT = getPDFUncert(f,"TT")
  out.write("%s CMS_PDF_Scales lnN %1.3f/%1.3f  -  %1.3f/%1.3f  -\n" % (mass, upLL, downLL, upTT, downTT))
  out.close()
