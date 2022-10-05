#Executable commands
P3=python3.9
BEALIGN=bealign 
BAM2MSA=bam2msa
HYPHY=hyphy
HYPHYMPI="mpirun -np 6 HYPHYMPI"
RXML=raxml-ng
TN93=tn93-cluster
FASTA_DIFF=fasta_diff


T="${3:-0.00075}"
#clustering threshold

T2="${4:-0.002}"
#outliner threshold

T3="${5:-3}"
#minimum number of sequences per cluster to include in final tree

REF="${6:-NONE}"
#add sequences from this file as references (outgroups)

FILE=$1
PREV=$2

#echo "Removing spaces from sequence names "
#awk '{ if ($0 ~ "^>") {sub(" ", "_"); print ;} else print;}' $FILE > ${FILE}.rn
#mv ${FILE}.rn ${FILE}

echo "Trimming down to the S neighborhood and removing spaces from sequence names"
#$P3 python/filter-sites.py $FILE  20000,26000 > ${FILE}.S.raw

if [ -z "$PREV" ] || [ $PREV == "NONE" ]
then
    echo "No previous run data"
    echo "Running bealign on the entire alignment"
    #$BEALIGN -r CoV2-S ${FILE}.S.raw ${FILE}.S.bam
    #$BAM2MSA ${FILE}.S.bam ${FILE}.S.msa

else
    echo "Extracting new/changed sequences"
    $FASTA_DIFF  -p remove -t id_sequence -m ${FILE}.S.raw ${PREV}.S.raw > ${FILE}.raw.diff
    echo "Running bealign on the new sequences"
    $BEALIGN -r CoV2-S ${FILE}.raw.diff  ${FILE}.S.diff.bam
    $BAM2MSA ${FILE}.S.diff.bam ${FILE}.S.diff.msa    
    echo "Merging MSA files"
    $FASTA_DIFF -p replace -t id_sequence -m ${FILE}.S.diff.msa    ${PREV}.S.msa > ${FILE}.S.msa  
fi


echo "Compressing to unique haplotypes"
#$P3 python/exact-copies.py  ${FILE}.S.msa > ${FILE}.u.clusters.json
#$P3 python/cluster-processor.py ${FILE}.u.clusters.json > ${FILE}.S.u.fas

#$TN93 -f -t 0.0 ${FILE}.S.u.fas > ${FILE}.t0.clusters.json
#$P3 python/cluster-processor.py ${FILE}.t0.clusters.json > ${FILE}.S.all.fas

#echo "Compressing to haplotypes at $T distance and checking for outliers"
#$TN93 -f -t $T ${FILE}.S.all.fas > ${FILE}.t1.clusters.json
#$P3 python/cluster-processor.py ${FILE}.t1.clusters.json > ${FILE}.S.uniq.fas
#$TN93 -f -t $T2 ${FILE}.S.uniq.fas > ${FILE}.t2.clusters.json
#$P3 python/cluster-processor.py ${FILE}.t1.clusters.json ${FILE}.t2.clusters.json > ${FILE}.S.uniq-all.fas 2> ${FILE}.S.blacklist.txt


echo "Rebuilding consensus and striking orphans"

$P3 python/cluster-processor-consensus.py $T3 ${FILE}.S.msa ${FILE}.S.uniq-all.fas ${FILE}.t1.clusters.json ${FILE}.t0.clusters.json ${FILE}.u.clusters.json > ${FILE}.S.uniq.fas 

if [ -z "$REF" ] || [ $REF == "NONE" ]
then
    echo "No reference sequences to add"
else
    echo "Appending reference sequences"
    cat $REF >> ${FILE}.S.uniq.fas 
fi

echo "Running raxml"
$RXML --redo --threads 5 --msa ${FILE}.S.uniq.fas --tree pars{5} --model GTR+G+I

echo "Labeling trees"
if [ -z "$REF" ] || [ $REF == "NONE" ]
then
    echo "No reference sequences to add"
    cp ${FILE}.S.uniq.fas.raxml.bestTree ${FILE}.S.labeled.nwk
    cp ${FILE}.S.uniq.fas.raxml.bestTree ${FILE}.S.labeled-internal.nwk
    ALL="All"
    INT="Internal"
else
    $HYPHY hbl/tree-labeler.bf --tree ${FILE}.S.uniq.fas.raxml.bestTree --output ${FILE}.S.labeled.nwk --output-internal ${FILE}.S.labeled-internal.nwk --reroot $REF
    ALL="CLADE"
    INT="CLADE"
fi
 
echo "Running SLAC"
$HYPHY slac --alignment ${FILE}.S.uniq.fas --kill-zero-lengths Constrain --tree ${FILE}.S.labeled.nwk --branches $ALL --samples 0


$HYPHY hbl/slac-mapper.bf ${FILE}.S.uniq.fas.SLAC.json ${FILE}.subs.json

echo "Running BGM"
$HYPHY bgm --alignment ${FILE}.S.uniq.fas --tree ${FILE}.S.labeled-internal.nwk --branches $INT --min-subs 2 --steps 1000000 --samples 1000 --burn-in 100000

echo "Running BUSTED[S]"
$HYPHY busted  --alignment ${FILE}.S.uniq.fas --tree ${FILE}.S.labeled-internal.nwk --branches $INT --rates 3 --starting-points 5

echo "Running FEL"


$HYPHYMPI fel --alignment ${FILE}.S.uniq.fas  --tree ${FILE}.S.labeled-internal.nwk --branches $INT  --ci Yes

echo "Running MEME"
$HYPHYMPI meme --alignment ${FILE}.S.uniq.fas  --tree ${FILE}.S.labeled-internal.nwk --branches $INT 

echo "Extracting variant counts"
$P3 python/countN.py ${FILE}.S.msa ${FILE}.S.uniq.fas ${FILE}.S.all.fas  > ${FILE}.S.json 2>${FILE}.S.genomes.json

echo "Extracting cluster traces"

$P3 python/global-site-extract.py ${FILE}.S.msa  339,371,373,375 2> ${FILE}.cluster1.json
$P3 python/global-site-extract.py ${FILE}.S.msa  493,496,498,505 2> ${FILE}.cluster2.json
$P3 python/global-site-extract.py ${FILE}.S.msa  764,856,954,969,981 2> ${FILE}.cluster3.json
$P3 python/global-site-extract.py ${FILE}.S.msa  417,452 2> ${FILE}.cluster4.json