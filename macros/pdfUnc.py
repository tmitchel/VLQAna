#!/usr/bin/env python

import ROOT
import math
from glob import glob
from optparse import OptionParser

parser = OptionParser('analysis')
parser.add_option('-d', '--dir', dest = 'dir',
                  help = 'directory to read from'
                  )

(options, args) = parser.parse_args()

def getPDF(fin):
  print fin
  ifile = ROOT.TFile(fin)
  pdfs = range(10, 111)

  upPdfs_pre = []
  upPdfs_post = []
  dnPdfs_pre = []
  dnPdfs_post = []

  for ipdf in pdfs:
    pdf = str(ipdf)
    
    if ipdf % 2 == 0:
      upPdfs_pre.append(ifile.Get('pdfProcess'+pdf+'/pre_pdf').Integral())
      upPdfs_post.append(ifile.Get('pdfProcess'+pdf+'/post_pdf').Integral())
    elif ipdf % 2 == 1:
      dnPdfs_pre.append(ifile.Get('pdfProcess'+pdf+'/pre_pdf').Integral())
      dnPdfs_post.append(ifile.Get('pdfProcess'+pdf+'/post_pdf').Integral())

  nom_pre = ifile.Get('pdfProcess10/pre_pdf').Integral()
  nom_post = ifile.Get('pdfProcess10/post_pdf').Integral()

  nom_ratio = nom_post / nom_pre
  preUp_sum, postUp_sum, preDn_sum, postDn_sum = 0, 0, 0, 0

  for preUp in upPdfs_pre:
    preUp_sum += pow(preUp - nom_pre, 2)
  preUp_sum = math.sqrt(preUp_sum/100)

  for postUp in upPdfs_post:
    postUp_sum += pow(postUp - nom_post, 2)
  postUp_sum = math.sqrt(postUp_sum/100)

  for preDn in dnPdfs_pre:
    preDn_sum += pow(preDn - nom_pre, 2)
  preDn_sum = math.sqrt(preDn_sum/100)

  for postDn in dnPdfs_post:
    postDn_sum += pow(postDn - nom_post, 2)
  postDn_sum = math.sqrt(postDn_sum/100)

  up_ratio = postUp_sum / preUp_sum
  dn_ratio = postDn_sum / preDn_sum

  return up_ratio/nom_ratio, dn_ratio/nom_ratio

fins = glob('*.root')
out = open('pdfuncerts_v3.log', 'w')
for fin in fins:
  up, down = getPDF(fin)
  #out.write('mass %s: \n \t scales: up - %1.3f, down - %1.3f\n\n' % (fin.split('/')[-1], up, down))
  out.write("%s %1.3f %1.3f  -\n" % (fin.split('/')[-1], up, down))
