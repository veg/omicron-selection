"""
Combine analysis results for Omicron selection

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2021-01-17)


"""

import argparse
import csv
import random
import os
import json
import sys
import datetime, math
from collections import Counter
from operator import itemgetter


arguments = argparse.ArgumentParser(description='Combine analysis results for Omicron selection')
arguments.add_argument('-i', '--input',  help = 'Directory to scan', required = True, type = str)
arguments.add_argument('-o', '--output',  help = 'Save JSON here', required = True, type = str)
settings = arguments.parse_args()

msa_name = lambda date: date.strftime ("%Y%m%d")
file_name = lambda b,f,e: os.path.join (b,f + e)

required_extensions = {
    "variants"  : ".S.json",
    "fel"       : ".S.uniq.fas.FEL.json",
    "meme"      : ".S.uniq.fas.MEME.json",
    "busted"    : ".S.uniq.fas.BUSTED.json",
    "slac"      : ".S.uniq.fas.SLAC.json",
    "bgm"       : ".S.uniq.fas.BGM.json",
    "Cluster 1" : ".cluster1.json",
    "Cluster 2" : ".cluster2.json",
    "Cluster 3" : ".cluster3.json"
}

max_date = datetime.datetime (2000,1,1)
combined_data                 = {}
max_p                         = 0.05

nucs = set (['A','C','G','T'])
ever_selected = set ()
site_info     = []
dates         = []
substitutions = None
clusters   = {}

def convert_d (d, key):
    result = []
    for i in range (len (d)):
        result.append ([d[str(i)][key],d[str(i)]['proportion']])
    return result

for root, dirs, files in os.walk(settings.input):
    dir_name = os.path.basename(root) 
    try:
        dir_date  = datetime.datetime.strptime (dir_name, "%Y-%m-%d")
        base_name = msa_name(dir_date)
        required_names = set ([base_name + re for re in required_extensions.values()])
        matched_names = [k for k in files if k in required_names]
        if len (matched_names) < len (required_names):
            continue
        
        is_max = False
        if dir_date > max_date:
            max_date = dir_date
            is_max = True
            
        directory_record = {}
        N_rich = set ()
        problematic = set ()
        
        with open (file_name (root, base_name, required_extensions["variants"]), "r") as fh:
            variant_json = json.load (fh)
            directory_record["N"] = variant_json["N"]
            directory_record["H"] = variant_json["H"]
            directory_record["HA"] = variant_json ["HA"]
            
            ## describe haplotype data; identify variants that are 
            # gappy or have too many Ns
             
            allN = 0
            for site, variants in enumerate (variant_json ["haplotypes"]):
                
                Ns   = 0
                gaps = 0
                fs   = 0
                amb  = 0
                
                for codon, ccount in variants.items():
                    if codon == 'NNN': 
                        Ns += ccount
                    elif codon == '---':
                        gaps += ccount
                    else:
                        codon_letters = set ([k for k in codon])
                        d = codon_letters.difference (nucs)
                        if len (d) > 0:
                            if '-' in d:
                                fs += ccount
                            else:
                                amb += ccount
                                
                if Ns * 10 > variant_json["H"] or amb * 20 > variant_json["H"]:
                    N_rich.add (site + 1)
                if gaps * 20 > variant_json["H"] and gaps * 3 < variant_json["H"]:
                    problematic.add (site + 1)
                if gaps * 2 > variant_json["H"]:
                    problematic.add (site + 1)
                if fs * 100 > variant_json["H"]:
                    problematic.add (site + 1)
                    
                allN += Ns
            
            directory_record ["NNN"] = allN / variant_json["H"] 
                    
        
            
        with open (file_name (root, base_name, required_extensions["slac"]), "r") as fh:
            slac_info = json.load (fh)
            L    = len (slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"])
            S    = sum ([k[2] for k in slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]]) / L
            S2   = sum ([k[2]*k[2] for k in slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]])
            NS   = sum ([k[3] for k in slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]]) / L
            NS2   = sum ([k[3]*k[3] for k in slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]])
            
            cutoff_S = S + math.sqrt (S2/L - S*S)*2
            cutoff_NS = NS + math.sqrt (NS2/L - NS*NS)*2
            
            outliers = set ([i+1 for i, k in enumerate (slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]) if k[2]>=cutoff_S and k[3] >= cutoff_NS])
            for site, slac_data in enumerate (slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]):
                if len (site_info) <= site:
                    site_info.append ({'fel' : [], 'meme' : [], 'slac' : []})
                
                site_info[site]['slac'].append (
                    [0,0]
                )
                
            for branch, slac_data in slac_info["branch attributes"]["0"].items():
                if branch[0:4] == "Node":
                    for site,k in enumerate (slac_data["nonsynonymous substitution count"][0]):
                        site_info[site]['slac'][-1][1] += k
                        site_info[site]['slac'][-1][0] += slac_data["synonymous substitution count"][0][site]
            
                
            problematic.update (outliers)
        
        directory_record ['issues'] = {'problematic' : sorted (list(problematic)), 'N rich' : sorted (list (N_rich))}
        
        with open (file_name (root, base_name, required_extensions["meme"]), "r") as fh:
            meme_info = json.load (fh)
            directory_record ['omega'] = {
                    "leaves" : meme_info["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *background*"][0][0],
                    "internal" : meme_info["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *test*"][0][0]
            }
            for site, meme_data in enumerate (meme_info["MLE"]["content"]["0"]):
                    
                site_info[site]['meme'].append (
                    [meme_data[6],meme_data[7],meme_data[3], meme_data[4]]
                    #pv, branches, omega+, p+
                )
                
                if meme_data[6] <= max_p:
                    ever_selected.add (site)

        with open (file_name (root, base_name, required_extensions["fel"]), "r") as fh:
            fel_info = json.load (fh)
            for site, fel_data in enumerate (fel_info["MLE"]["content"]["0"]):
                     
                site_info[site]['fel'].append (
                    [fel_data[4],fel_data[0],fel_data[1],fel_data[5],fel_data[6],fel_data[7],fel_data[8]]
                )
                #pv, alpha, beta, branch length, omega, LB, UB
                
                if fel_data[4] <= max_p:
                    ever_selected.add (site)


        
        if is_max:
            try:
                with open (file_name (root, base_name, ".subs.json"), "r") as fh:
                    substitutions = json.load (fh)
            except Exception as e:
                print (e, file = sys.stderr)
                
            
            for e in ["Cluster 1","Cluster 2","Cluster 3"]:
                with open (file_name (root, base_name, required_extensions[e])) as fh:
                    clusters[e] = json.load (fh)
        
            with open (file_name (root, base_name, required_extensions["bgm"]), "r") as fh:
                bgm_info = json.load (fh)
                pairs    = []
                if "MLE" in bgm_info:
                    for pair_info in bgm_info["MLE"]["content"]:
                        if pair_info [4] >= 0.8:
                            pairs.append (pair_info)
  
            with open (file_name (root, base_name, required_extensions["busted"]), "r") as fh:
                busted_info = json.load (fh)
                #print (busted_info["fits"]["Unconstrained model"]["Rate Distributions"])
                busted = {
                    'p' : busted_info ["test results"]["p-value"],
                    'leaves' : convert_d (busted_info["fits"]["Unconstrained model"]["Rate Distributions"]["Background"], "omega"),
                    'internal' : convert_d (busted_info["fits"]["Unconstrained model"]["Rate Distributions"]["Test"], "omega"),
                    'srv' : convert_d (busted_info["fits"]["Unconstrained model"]["Rate Distributions"]["Synonymous site-to-site rates"], "rate"),
                }
            counts = variant_json ["counts"]
            haplos = variant_json ["haplotypes"]
            haplos_all = variant_json ["all-haplotypes"]
            seq_dates = variant_json ["dates"]

              
            
        dates.append (base_name)
        combined_data [base_name] = directory_record
         
        
    except Exception as e:
        print (dir_name, e, file = sys.stderr)
        pass
             
sorted_dates = sorted ([[k,i] for i,k in enumerate (dates)], key = itemgetter (0))

combined_data ['analysis dates'] = [k[0] for k in sorted_dates]
combined_data ['busted'] = busted

site_selection_info = {}
            
for k in sorted (list (ever_selected)):
   site_selection_info [k+1] = {
        'fel' :  [site_info[k]['fel'][d[1]] for d in sorted_dates],
        'meme' : [site_info[k]['meme'][d[1]] for d in sorted_dates],
        'slac' : [site_info[k]['slac'][d[1]] for d in sorted_dates],
   }
                
combined_data ['p'] = max_p
combined_data ['sites'] = site_selection_info
combined_data ['pairs'] = pairs
combined_data ['counts'] = counts
combined_data ['haplos'] = haplos
combined_data ['haplos_all'] = haplos_all
combined_data ['dates'] = seq_dates
combined_data ['subs'] = substitutions
combined_data ['clusters'] = clusters
  
with open (settings.output, 'w') as fh:
    json.dump (combined_data, fh)
