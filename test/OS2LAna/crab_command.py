
#!/bin/python

from string import *
import subprocess
import os
from allJobList import *

from optparse import OptionParser
parser = OptionParser()

parser.add_option('--command', metavar='C', type='string', action='store',
                  default='submit',
                  dest='command',
                  help='crab command to execute')

parser.add_option('--channel', metavar='S', type='string', action='store',
                  default='Zmumu',
                  dest='channel',
                  help='skim type')

#parser.add_option('--signal', metavar='S', type='int', action='store',
#                  default = '0',
#                  dest='signal',
#                  help='check signal MC or not')

parser.add_option('--isData', metavar='I', type='int', action='store',
                  default=0,
                  dest='isData',
                  help='run on data or MC')

parser.add_option('--inputCfg', metavar='C', type='string', action='store',
                  default='crab_dummy_config.py',
                  dest='inputCfg',
                  help='input config tag to be used')

parser.add_option('--outLabel', metavar='L', type='string', action='store',
                  default='reskim',
                  dest='outLabel',
                  help='output tag to be used')

parser.add_option('--user', metavar='C', type='string', action='store',
                  default='tmitchel',
                  dest='user',
                  help='user subdir path')

(options, args) = parser.parse_args()
print options
command = options.command
user = options.user

jecdict = {
    'B': "Summer16_23Sep2016BCDV4_DATA",
    'C': "Summer16_23Sep2016BCDV4_DATA",
    'D': "Summer16_23Sep2016BCDV4_DATA",
    'E': "Summer16_23Sep2016EFV4_DATA",
    'F': "Summer16_23Sep2016EFV4_DATA",
    'G': "Summer16_23Sep2016GV4_DATA",
    'H': "Summer16_23Sep2016HV4_DATA",
    }

#set the parameter options in os2lana.py
if options.isData == 1:
    isData  = 'isData=True'
else:
    isData = 'isData=False'

isphoton = 'isPhoton=False'
   
if 'Zmumu' in options.channel:
    mode = 'zdecaymode=zmumu'
    print 'Submitting jobs for muon channel -------->'
    if options.isData:
        jobList = list_Zmumu_data()
    else: jobList = list_MC()
elif 'Zelel' in options.channel:
    mode = 'zdecaymode=zelel'
    print 'Submitting jobs for electron channel -------->'
    if options.isData:
        jobList = list_Zelel_data()
    else: jobList = list_MC()
elif 'photon' in options.channel:
  mode = 'zdecaymode=zelel'
  isphoton = '"isPhoton=True"'
  print 'Submitting jobs for photon channel ------->'
  if options.isData:
    jobList = list_photon_data()

for job in jobList:
    print job  
    print '------prepare to run on :  ' + job[0] + ' -----------'
    jecname = job[1].split('2016')[-1].split('_')[0]
    f = open(options.inputCfg, 'r')
    instring = f.read()
    baseList = job[0].split('/')
     
    outname = job[1] + '_' + options.outLabel
     #print outname
    a0 = instring.replace( 'DUMMY_WORKDIR', "'"+options.outLabel+"'")
    a1 = a0.replace( 'DUMMY_DATASET', "'"+job[0]+"'" )
    a2 = a1.replace( 'DUMMY_NUMBER',  job[2])
    a3 = a2.replace( 'DUMMY_NAME', "'"+outname+"'" )
    a4 = a3.replace( 'DUMMY_SITE',"'"+'T3_US_FNALLPC'+"'")
    a5 = a4.replace( 'DUMMY_OUTPUT_PATH', "'"+'/store/user/'+user+'/skim2018/'+options.channel+"/'")
    a6 = a5.replace( 'DATA', "'"+isData+"'")
    a7 = a6.replace( 'MODE', "'"+mode+"'")
   
    if 'Tprime' in baseList[1] or 'Bprime' in baseList[1]:
        a8 = a7.replace( 'FILTERSIGNAL', "'filterSignal=True'")     
    else:
        a8 = a7.replace( ', FILTERSIGNAL', '')      
    a9 = a8.replace( 'JECNAME', "'newJECPayloadNames="+jecdict[jecname]+"'")
    a10 = a9.replace( 'ISPHOTON', isphoton )
    #if 'DY' in  baseList[1]: a9 = a8.replace( 'DYNLO', "'applyDYNLOCorr=True'")
    #else: a9 = a8.replace( 'DYNLO', "'applyDYNLOCorr=False'")
    print '------ Config : ------- '
    print a10
    
    # write the crab config files .........
    c = 'mkdir '+options.channel
    if not os.path.isdir(options.channel):
        subprocess.call( [c], shell=True )    
    crabName = 'crab_' + outname + '.py'
    fout = open( options.channel+'/'+crabName, 'w')
    fout.write( a10 )
    fout.close()
    
    print '------ CRAB starting up! ------'
    # now run/check the crab jobs:
    if command == 'submit':
        exe = 'crab submit -c ' + options.channel+'/'+crabName
    elif command == 'report' or command == 'status':
        exe = 'crab '+options.command+' -d '+options.channel+'/'+crabName
    print exe
    print '--------------->'
    
    subprocess.call( [exe], shell=True )
