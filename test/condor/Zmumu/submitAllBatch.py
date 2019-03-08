#!/usr/bin/env python

from subprocess import call
from glob import glob

fileList = [ifile for ifile in glob('*') if '.jdl' in ifile]

for ifile in fileList:
  call("condor_submit "+ifile, shell=True)
