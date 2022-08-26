#!/bin/bash

mkdir tmp
mkdir marker_genes_from_genes
mkdir marker_genes_from_genome

python ./../extract_gtdb_mg.py --threads 4 --genome GUT_GENOME103755.fna --overwrite_tmp --tmp tmp/ -o marker_genes_from_genome/

python ./../extract_gtdb_mg.py --threads 4 --genome GUT_GENOME103755.fna --overwrite_tmp --tmp tmp/ -o marker_genes_from_genome/ --list_only


python ./../extract_gtdb_mg.py --threads 4 --aa GUT_GENES103755.faa --overwrite_tmp --tmp tmp/ -o marker_genes_from_genes/

python ./../extract_gtdb_mg.py --threads 4 --aa GUT_GENES103755.faa --overwrite_tmp --tmp tmp/ --nt GUT_GENES103755.fna -o marker_genes_from_genes/

python ./../extract_gtdb_mg.py --threads 4 --aa GUT_GENES103755.faa --overwrite_tmp --tmp tmp/ --nt GUT_GENES103755.fna -o marker_genes_from_genes/ --list_only

