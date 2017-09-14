#!/usr/bin/env python

import os, sys
from ROOT import TH2D, TFile, gROOT
from array import array
from FWCore.ParameterSet.VarParsing import VarParsing

gROOT.SetBatch(1)

options = VarParsing('analysis')
options.register('path', '/uscms_data/d3/tmitchel/80X/CMSSW_8_0_20/src/Analysis/VLQAna/macros/',
	VarParsing.multiplicity.singleton,
	VarParsing.varType.string,
	'path to files'
)
options.register('FileName', 'bprime1200_bZ_el',
	VarParsing.multiplicity.singleton,
	VarParsing.varType.string,	
	'file name'
)
options.parseArguments()

inputFile = TFile(options.path+options.FileName+'.root', 'READ')
outputFile = TFile(options.path+options.FileName+'_bTagEff.root', 'RECREATE')

bins = {
	'b':	[[0., 40., 60., 80., 100., 150., 200., 300., 400., 500., 600., 700., 800., 1000.],[0., 0.6, 1.2, 2.4]],
	'c':	[[0., 40., 60., 80., 100., 150., 200., 300., 1000.],[0., 0.6, 1.2, 2.4]],
	'l':	[[0., 40., 60., 80., 100., 150., 200., 1000.],[0., 0.6, 1.2, 2.4]]
}

for flav in ['b', 'c', 'l']:

	nIn = inputFile.Get('ana/pt_eta_'+flav+'_btagged')
	dIn = inputFile.Get('ana/pt_eta_'+flav+'_all')	
	
	xShift = dIn.GetXaxis().GetBinWidth(1)/2.
	yShift = dIn.GetYaxis().GetBinWidth(1)/2.
	
	xBins = array('d', bins[flav][0])
	yBins = array('d', bins[flav][1])

	dOut = TH2D('denom_'+flav, "Denominator for "+flav+' in '+options.FileName, (len(xBins)-1), xBins, (len(yBins)-1), yBins)
	nOut = TH2D('numer_'+flav, "Numerator for "+flav+' in '+options.FileName, (len(xBins)-1), xBins, (len(yBins)-1), yBins)
	effOut = TH2D('eff_'+flav, "Efficiency for "+flav+' in '+options.FileName, (len(xBins)-1), xBins, (len(yBins)-1), yBins)

	for i in range(1, dOut.GetXaxis().GetNbins()+1):
		for j in range(1, dOut.GetYaxis().GetNbins()+1):

			binXMin = dIn.GetXaxis().FindBin(dOut.GetXaxis().GetBinLowEdge(i)+xShift)
			binXMax = dIn.GetXaxis().FindBin(dOut.GetXaxis().GetBinUpEdge(i)-xShift)
			binYMinPos = dIn.GetYaxis().FindBin(dOut.GetYaxis().GetBinLowEdge(j)+yShift)
			binYMaxPos = dIn.GetYaxis().FindBin(dOut.GetYaxis().GetBinUpEdge(j)-yShift)
			binYMaxNeg = dIn.GetYaxis().FindBin(-dOut.GetYaxis().GetBinLowEdge(j)-yShift)
			binYMinNeg = dIn.GetYaxis().FindBin(-dOut.GetYaxis().GetBinUpEdge(j)+yShift)

			denom = dIn.Integral(binXMin,binXMax,binYMinPos,binYMaxPos)
			denom += dIn.Integral(binXMin,binXMax,binYMinNeg,binYMaxNeg)
			numer = nIn.Integral(binXMin,binXMax,binYMinPos,binYMaxPos)
			numer += nIn.Integral(binXMin,binXMax,binYMinNeg,binYMaxNeg)
	
			if (i==dOut.GetXaxis().GetNbins()):
				denom += dIn.Integral(binXMax+1,dIn.GetXaxis().GetNbins()+1,binYMinPos,binYMaxPos)
				denom += dIn.Integral(binXMax+1,dIn.GetXaxis().GetNbins()+1,binYMinNeg,binYMaxNeg)
				numer += nIn.Integral(binXMax+1,nIn.GetXaxis().GetNbins()+1,binYMinPos,binYMaxPos)
				numer += nIn.Integral(binXMax+1,nIn.GetXaxis().GetNbins()+1,binYMinNeg,binYMaxNeg)

			dOut.SetBinContent(i,j,denom)
			nOut.SetBinContent(i,j,numer)
			if denom != 0.:
				effOut.SetBinContent(i,j,numer/denom)

	for i in range(1, dOut.GetXaxis().GetNbins()+1):
		for j in range(1, dOut.GetYaxis().GetNbins()+1):
			efficiency = effOut.GetBinContent(i,j)
			if efficiency == 0 or efficiency == 1:
				print 'AHHHHHHHHHHH!!!! Bin(%i, %i) for %s jets has a b-tag efficiency of %.3f'%(i,j,flav,efficiency)
	
	for i in range(1, dOut.GetXaxis().GetNbins()+1):
		effOut.SetBinContent(i, dOut.GetYaxis().GetNbins()+1, effOut.GetBinContent(i, dOut.GetYaxis().GetNbins()))

	for i in range(1, dOut.GetYaxis().GetNbins()+1):
		effOut.SetBinContent(dOut.GetXaxis().GetNbins()+1,j, effOut.GetBinContent(dOut.GetXaxis().GetNbins(), j))

	outputFile.cd()

	dOut.Write()
	nOut.Write()
	effOut.Write()

outputFile.Close()




















































