#!/bin/bash

mkdir tmp
mkdir marker_genes_from_genes
mkdir marker_genes_from_genome

python ./../extract_gtdb_mg.py --genome GUT_GENOME103755.fna --tmp tmp/ -o marker_genes_from_genome/

python ./../extract_gtdb_mg.py --aa GUT_GENES103755.faa --tmp tmp/ -o marker_genes_from_genes/

python ./../extract_gtdb_mg.py --aa GUT_GENES103755.faa --nt GUT_GENES103755.fna -o marker_genes_from_genes/

