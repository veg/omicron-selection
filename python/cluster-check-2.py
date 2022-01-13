import sys, json
from Bio import SeqIO

sequence_ids = set ()
by_id = {}

def getN (id):
    return int(id.split ('||')[1])

M = 0

with open(sys.argv[2]) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        try:
            sid = record.id.upper()
            sequence_ids.add (sid)
            M += getN (sid)
        except Exception as e:
            pass
               

cluster_ids = set ()

with open (sys.argv[1], "r") as fh:
    cluster_json = json.load (fh)
    N = 0
    for k in cluster_json:
        for m in k['members']:
            cluster_ids.add (m.upper())
            by_id [m.split('|')[1].upper()] = m
            N += getN (m)
            
N2 = 0
M2 = 0            
            
for k in (sequence_ids):
    try:
        m = getN (k)
        n = getN (by_id [k.split('|')[1].upper()])
        #print (m, '=>', n)
        M2 += m
    except:
        print ("#### ERROR MISSING ID %s" % k)

print (len (sequence_ids), "sequences, ", N, "total in clusters", M, M2, "total in sequence", len (cluster_ids))