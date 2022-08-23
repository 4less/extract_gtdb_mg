# Extract GTDB mg

This tool extracts GTDB marker genes (both nucleotide and aa) from a set of genomes or genes (protein sequence required from genes). GTDB-tk requires a downloaded database which adds verbosity and also does not output the nucleotide sequence of marker genes (only protein).

Code for this is taken from GTDB-tk and was changed/arranged accordingly.

#### Input

genome(s) **OR** genes (amino acid) + (optional = nucleotide)

#### Output

Output: marker genes according to GTDB



**Note:** To get nucleotide marker genes from existing gene predictions you need to provide options --aa and --nt



### Requirements:

Python version >3

hmmsearch

prodigal

## Options

```
usage: extract_gtdb_mg.py [-h] [-a AA] [-n NT] [-g GENOME] -o OUTPUT_DIR -p TMP [-k] [-x]

optional arguments:
  -h, --help            show this help message and exit
  -a AA, --aa AA        One or more files with genes (Protein) (Mutually exclusive to (--genomes))
  -n NT, --nt NT        One or more files with genes (DNA) (Mutually exclusive to (--genomes) and
                        requires --aa)
  -g GENOME, --genome GENOME
                        Genome file(s) (Mutually exclusive to (--aa/--nt))
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Marker_gene output folder
  -p TMP, --tmp TMP     Temporary folder
  -k, --keep_tmp        Keep temporary files
  -x, --overwrite_tmp   Ignore if tmp dir already exists
```





## Examples

Test file in test/test.sh

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

