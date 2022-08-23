# Extract GTDB mg

Test file in test/test.sh

Input: genome(s) **OR** genes (amino acid) + (optional = nucleotide)

Output: marker genes according to GTDB

**Note:** To get nucleotide marker genes from existing gene predictions you need to provide options --aa and --nt



Requirements:

Python version >3

hmmsearch

prodigal



##### From genomes

This step first runs prodigal to get gene predictions

```bash
# From genomes
python ./../extract_gtdb_mg.py --genome GUT_GENOME103755.fna --tmp tmp/ -o marker_genes_from_genome/
```

---

This step requires that you have already run e.g. prodigal to get gene predictions.

##### From genes (aa)

```bash
# from genes
python ./../extract_gtdb_mg.py --aa GUT_GENES103755.faa --tmp tmp/ -o marker_genes_from_genes/
```



##### From genes (nt and aa) 

```bash
# from genomes
python ./../extract_gtdb_mg.py --aa GUT_GENES103755.faa --nt GUT_GENES103755.fna -o marker_genes_from_genes/
```

