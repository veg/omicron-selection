#Scripts for the selection analysis of the BA.1/BA.2 S gene sequences

Python, shell, and HyPhy scripts provided in this repository are used to preprocess, reduce, and analyze [GISAID](http://gisaid.org) data from Omicron (but can also do others) SARS-CoV-2 S-gene sequences and preparing data for visualization using [an ObservableHQ notebook](https://observablehq.com/@spond/ba1-selection)

This scripts are **fragile** -- they depend on the specific format of GISAID exports (e.g. metainformation in the sequence names) and the availability of tools on the system


#### Dependencies

1. `HyPhy ≥ 2.5.34` : [github.com/veg/hyphy](http://github.com/veg/hyphy)
2. `BioExt ≥ 0.20.4` for `bealign` and `bam2msa`: [github.com/veg/BioExt](http://github.com/veg/BioExt)
3. `TN93 ≥ 1.07` for `tn93-cluster` and `fasta_diff`: [https://github.com/veg/tn93](https://github.com/veg/tn93) 
4. `raxml-ng ≥ 1.03`: [github.com/amkozlov/raxml-ng](https://github.com/amkozlov/raxml-ng)

#### Workflow

* Obtain data from GISAID (or GenBank); whole genome SARS-CoV-2 FASTA sequences. Certain annotation features assume that the header of the sequence will include sampling date (split on `|` take third element) and sampling locations (split on `/` take second element).

`>hCoV-19/South_Africa/NCV289/2021|EPI_ISL_7852877|2021-12-02`

* Create a data directory (e.g. `BA1`), and within this directory (names are important!), create `YYYY-MM-DD` subdirectory (e.g. `2022-01-15`) to indicate the date when the analysis was run. The raw FASTA file needs to go into this directory names `YYYYMMDD` (no extension, e.g. `BA1/2022-01-15/20220115`).

* Run the analysis script (`run.sh`) on the file 

> Please adjust the executable command block at the top of the script to reflect where various executables live on your system

```
P3=python3.9
BEALIGN=bealign
BAM2MSA=bam2msa
HYPHY=hyphy
HYPHYMPI="mpirun -np 8 HYPHYMPI"
RXML=raxml-ng
TN93=tn93-cluster
FASTA_DIFF=fasta_diff
```

```
sh run.sh BA1/2022-01-15/20220115
```

If you have a previous analysis, you can also run the script like so 

```
sh run.sh BA1/2022-01-15/20220115 BA/2022-01-10/20220110
```

where `BA/2022-01-10/20220110` points to a past analysis. This option reduced some of the work but trying to reuse some of the previous results.

There are several other POSITIONAL arguments that can be passed tot `run.sh`

```
sh run.sh FILE PREV T T2 T3 REF
```

`T` controls how to create clusters of similar sequences. All sequences that are within `T` subs/site distance of each other become one cluster. By default `T=0.00075`

`T2` controls how to identify outlier/divergence sequences that are EXCLUDED from the analysis. All sequences that T2 or more units away from the largest (consensus) cluster are flagged as suspect. By default `T=0.002`

`T3` controls how to filter "small" clusters. Only sequences belonging to cluster with `T3` or more members will go into the tree (default = 5)

`REF` points to a file that contains outgroup sequences (e.g. BA.2 and BA.3)

Example calls 


```
sh run.sh BA.2/2022-01-13/20220113 data-BA.2/2022-01-05/20220105 0.0005 0.002 1 references/BA2.msa
```

```
sh run.sh BA.1/2022-01-13/20220113 data-BA.2/2022-01-12/20220112 0.00075 0.002 5 references/BA1.msa
```

* Create a JSON summary of all analysis in a single file

```
python3.9 python/generate-data-for-observable.py -i BA1 -o summary/BA1.json
```

* Navitage to [an ObservableHQ notebook](https://observablehq.com/@spond/omicron-selection-file) and load the resulting JSON file.
 