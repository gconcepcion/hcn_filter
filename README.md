# hcn_filter
High Copy Number filter for PacBio raw data

Use this script to filter out high copy number regions from PacBio raw data fasta files.


Dependencies:

Python/2.7.9
daligner/dazzlerDB

Usage:

$ python hcn_filter.py  /path/to/pacbio_subreads.fasta --genome_size 120 --ncpus 200 --debug 


* This script was developed for an SGE; though should be easily modifiable for most Job management systems.

