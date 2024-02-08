#!/bin/bash
#
# make a small test kraken2 database
#
# get from NCBI: NC_045512.2.fna
# conda install seqfu
#
seqfu cat --append "|kraken:taxid|2697049" NC_045512.2.fna.gz > NC_045512.2_taxid.fasta

exit

# create a new database and add fasta
time docker run --rm -u 1000:1000 -v $PWD:/mnt -w /mnt staphb/kraken2 kraken2-build --add-to-library NC_045512.2_taxid.fasta --db corona --threads 4

# build taxonomy
# this takes a while
time docker run --rm -u 1000:1000 -v $PWD:/mnt -w /mnt staphb/kraken2 kraken2-build --download-taxonomy --db corona

# build database
time docker run --rm -u 1000:1000 -v $PWD:/mnt -w /mnt staphb/kraken2 kraken2-build --build --db corona

# clean up
time docker run --rm -u 1000:1000 -v $PWD:/mnt -w /mnt staphb/kraken2 kraken2-build --clean --db corona

exit 0
