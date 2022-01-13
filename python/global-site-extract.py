import sys, csv, json
from   collections import Counter
from   operator import itemgetter, attrgetter
import datetime
import scipy.stats
from math import log

from Bio import SeqIO

TT = {"TTT":"F","TCT":"S","TAT":"Y","TGT":"C","TTC":"F","TCC":"S","TAC":"Y","TGC":"C","TTA":"L","TCA":"S","TAA":"*","TGA":"*","TTG":"L","TCG":"S","TAG":"*","TGG":"W","CTT":"L","CCT":"P","CAT":"H","CGT":"R","CTC":"L","CCC":"P","CAC":"H","CGC":"R","CTA":"L","CCA":"P","CAA":"Q","CGA":"R","CTG":"L","CCG":"P","CAG":"Q","CGG":"R","ATT":"I","ACT":"T","AAT":"N","AGT":"S","ATC":"I","ACC":"T","AAC":"N","AGC":"S","ATA":"I","ACA":"T","AAA":"K","AGA":"R","ATG":"M","ACG":"T","AAG":"K","AGG":"R","GTT":"V","GCT":"A","GAT":"D","GGT":"G","GTC":"V","GCC":"A","GAC":"D","GGC":"G","GTA":"V","GCA":"A","GAA":"E","GGA":"G","GTG":"V","GCG":"A","GAG":"E","GGG":"G","---":"-","NNN":"?"}

sites = [int (k) for k in sys.argv[2].split (",")]
patterns_by_date = {}

counts_by_site = [Counter () for k in sites]
N = 0

pattern        =  [c for c in sys.argv[3]] if len (sys.argv) > 3 and len (sys.argv[3]) == 2  else None

with open(sys.argv[1]) as handle:
    counter = Counter()
    counter_aa = Counter()

    pattern_match = [0,0,0,0]

    for record in SeqIO.parse(handle, "fasta"):
        s = str (record.seq)
        try:
            D =  datetime.datetime.strptime (record.id.split ("|")[-1], "%Y-%m-%d").strftime ("%Y%m%d")
        except:
            D = None

        codon_pattern = ""
        aa_pattern    = ""
        p_idx         = 0
        
        for idx, i in enumerate (sites):
            codon = s[3*(i-1):3*i]
            codon_pattern += codon
            
            if codon in TT:
                aa = TT[codon]
            else:
                aa = '?'
                
            if pattern and aa == pattern[idx]:
                p_idx += (idx+1) 
                
                    
            aa_pattern += aa
            counts_by_site[idx][aa] += 1
            
        if aa_pattern not in patterns_by_date:
            patterns_by_date [aa_pattern] = Counter()
            
        pattern_match [p_idx] += 1
            
        patterns_by_date [aa_pattern][D] += 1
            
        counter[codon_pattern] += 1
        counter_aa[aa_pattern] += 1
        N += 1


if pattern:
    P1 = counts_by_site[0][pattern[0]]
    P2 = counts_by_site[1][pattern[1]]
    P12 = counter_aa [sys.argv[3]]
    
    if P1 > 0 and P2 > 0 and P12 > 0:
        print ("Pair log odds:", log (P12/N,2) - log (P1/N,2) - log (P2/N,2), file = sys.stderr)
        oddsratio, pvalue = scipy.stats.fisher_exact([[pattern_match[0], pattern_match[1]], [pattern_match[2], pattern_match[3]]])
        print ("OR = ", oddsratio, " P=", pvalue)
    else:
        print ("Log odds are not defined (or infinite)", file = sys.stderr)
    


json.dump (patterns_by_date, sys.stderr)

counter = sorted (counter.items(), key = itemgetter (1), reverse = True)            
            
for k, c in counter:
    print (k, c)
        
        
counter_aa = sorted (counter_aa.items(), key = itemgetter (1), reverse = True)            

for k, c in counter_aa:
    print (k, c)
              
