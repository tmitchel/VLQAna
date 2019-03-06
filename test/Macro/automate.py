#!/bin/python

import subprocess

options = [
['pt_elel-3j'],
['pt_el_lead-3j'],
['pt_el_2nd-3j'],
['eta_el_lead-3j'],
['eta_el_2nd-3j'],
['m_elel-3j'],
['dr_elel-3j'],
['npv_noweight'],
['npv'],
['nak4'],
['ht'],
['st'],
['met'],
['metPhi'],
['ptak4jet1'],
['etaak4jet1'],
['cvsak4jet1'],
['ptak4jet2'],
['etaak4jet2'],
['cvsak4jet2'],
['ptak4jet3'],
['etaak4jet3'],
['cvsak4jet3'],
['phi_jet1MET'],
['mass_zelel'],
['dr_elel'],
['pt_zelel'],
['pt_el1'],
['eta_el1'],
['pt_el2'],
['eta_el2'],
['nbjets'],
['ptbjetleading'],
['etabjetleading'],
['nak8'],
['nwjet'],
['nhjet'],
['ntjet'],
['ptak8leading'],
['etaak8leading'],
['prunedmak8leading'],
['softdropmak8leading'],
['ptak82nd'],
['etaak82nd'],
['prunedmak82nd'],
['softdropmak82nd'],
['subjetinessak8leading'],
['subjetinessak82nd'],
#['st_bZ_boost'],
#['st_bH_boost'],
#['st_bZ_1b'],
#['st_bH_1b'],
#['st_bZ_2b'],
#['st_bH_2b'],
#['boostReco_bZ'],
#['boostReco_bH'],
#['resReco_bZ_1b'],
#['resReco_bH_1b'],
#['resReco_bZ_2b'],
#['resReco_bH_2b'],
#['boostReco'],
#['resReco_1b'],
#['resReco_2b'],
['nak4-3j'],
['nak8-3j'],
['nhjets-3j'],
['nzjets-3j'],
['nbjets-3j'],
['pt_bjet-3j'],
['pt_ak4_lead-3j'],
['pt_ak4_2nd-3j'],
['pt_ak4_3rd-3j'],
['pt_ak4_4th-3j'],
['eta_ak4_lead-3'],
['eta_ak4_2nd-3j'],
['eta_ak4_3rd-3j'],
['eta_ak4_4th-3j'],
['st-3j'],
['ht-3j'],
['met-3j'],
['npv-3j'],

    
    ]

for option in options:
  option[0] += ''

command = 'python plot.py --var={0:s} --plotDir=refactor_v2'

for option in options :
    s = command.format(
        option[0]
        )

    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo %s"%s,""]                                                                      , shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( [s, ""]                                                                               , shell=True)
