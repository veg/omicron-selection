import sys, json

blacklist = set ()

if len(sys.argv) > 2:
    with open (sys.argv[2], "r") as fh:
        cluster_json = json.load (fh)
        largest = max ([len(k['members']) for k in cluster_json])
        for k in cluster_json:
            if k["size"] != largest:
                for n in k["members"]:
                    blacklist.add (n)

   

with open (sys.argv[1], "r") as fh:
    cluster_json = json.load (fh)
    for k in cluster_json:
        N = 0
        for m in k['members']:
            s = m.split ("||")
            if (len (s) > 1):
                N += int (s[1])
            else:
                N += 1
        seq = k['centroid'].split ("\n")
        if seq[0].split ('>')[1] in blacklist:
            print ("%s" % "\n".join (seq), file = sys.stderr)
        else:
            id = seq[0].split ("||")[0] + "||" + str (N)
            #print (N, file = sys.stderr)
            print (id)
            print (seq[1].strip())
        #NN += N
    
