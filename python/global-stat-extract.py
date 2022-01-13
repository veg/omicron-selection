import sys, csv, json
from   collections import Counter
from   operator import itemgetter, attrgetter
import datetime

from Bio import SeqIO

dates = Counter ()
countries = Counter ()

with open(sys.argv[1]) as handle:

    for record in SeqIO.parse(handle, "fasta"):
        try:
            D =  datetime.datetime.strptime (record.id.split ("|")[-1], "%Y-%m-%d").strftime ("%Y%m%d")
        except:
            D = None
            
        dates [D] += 1
        countries [record.id.split ("/")[1]] += 1
        
              
 
dates = sorted (dates.items(), key = itemgetter (1), reverse = True)            
            
for k, c in dates:
    print (k, c)

countries = sorted (countries.items(), key = itemgetter (1), reverse = True)            
            
for k, c in countries:
    print (k, c)
              
