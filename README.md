# bio_small_scripts
Miscellaneous scripts

## annotate_gff3.py
Add functional annotation to the 9th column of a gff3 file
```
$ ./annotate_gff3.py -h
usage: annotate_gff3.py [-h] -g GFF -b BLAST -f FUNC [-o OUT] [-p PREFIX]

add functional annotation to the 9th column of a gff3 file

optional arguments:
  -h, --help            show this help message and exit
  -g GFF, --gff GFF     Annotation in gff3 format (required)
  -b BLAST, --blast BLAST
                        tab-separated blast output (required) with "-outfmt "[6|7] qaccver saccver pident length
                        mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs""
  -f FUNC, --func FUNC  tab-separated function table (required)
  -o OUT, --out OUT     Annotation in gff3 format (default: stdout)
  -p PREFIX, --prefix PREFIX
                        locus ID prefix (default: gene)
```

## blast2dotplot.py 
```
$ ./blast2dotplot.py
usage: blast2dotplot.py [-h] --input FILE

Draw a dot plot based on a pairwise BLASTN/TBLASTX result

optional arguments:
  -h, --help            show this help message and exit
  --input FILE, --in FILE, -i FILE
                        input BLASTN/TBLASTX result file in XML format (-outfmt 5)
```
![NC_000913 3_NC_004337_BLASTN](https://user-images.githubusercontent.com/58936715/180001866-5f67723e-114d-42c1-aa75-8660e9aedba6.png)

## depth_alignment_breakpoint.py

## ddbj_to_gff3.py
```
$ ./ddbj_to_gff3.py -h
usage: ddbj_to_gff3.py [-h] -i INPUT [-o OUTPUT]

Convert DDBJ/GenBank flatfile into gff3

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input GenBank/DDBJ flatfile (required)
  -o OUTPUT, --output OUTPUT
                        output gff3 file (default:out.gff3)
```
## gb2nrfaa.py
Extract the longest isoforms of protein-coding genes from a NCBI RefSeq euaryotic genome assembly
```
$ ./gb2nrfaa.py -h
usage: gb2nrfaa.py [-h] -i INPUT [-o OUTPUT]

Extract the longest isoforms of protein-coding genes from a NCBI RefSeq euaryotic genome assembly

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        GenBank flat file format of the genomic sequence(s) (required)
  -o OUTPUT, --output OUTPUT
                        output fasta-formatted file (required)
```

## orfind.py
ORF prediction allowing CDS overlaps and GFF3 output
```
$ ./orfind.py -h
usage: orfind.py [-h] -i INPUT [-o OUT_GFF] [-a OUT_FAA] [-f OUT_FNA] [-g TRANS_TABLE] [-m MIN_AA_LEN]

Predict ORFs

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        sequence file in FASTA format (required)
  -o OUT_GFF, --out_gff OUT_GFF
                        output annotation in gff3 format (default: stdout)
  -a OUT_FAA, --out_faa OUT_FAA
                        output protein sequences in FASTA format (optional)
  -f OUT_FNA, --out_fna OUT_FNA
                        output CDS sequences in FASTA format (optional)
  -g TRANS_TABLE, --trans_table TRANS_TABLE
                        translation table (default: 1)
  -m MIN_AA_LEN, --min_aa_len MIN_AA_LEN
                        minimum protein length (default: 50)
```
## plot_circular_genome.py
Generate genome diagram in SVG
```
$ ./plot_circular_genome.py -h
usage: plot_circular_genome.py [-h] -i INPUT [-n NT] [-w WINDOW] [-s STEP]

Generate genome diagram in SVG.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Genbank/DDBJ flatfile (required)
  -n NT, --nt NT        dinucleotide (default: GC).
  -w WINDOW, --window WINDOW
                        window size (default: 1000)
  -s STEP, --step STEP  step size (default: 100)
```

## plot_skew.py
Generate dinucleotide skew plot(s) of FASTA format DNA sequences in SVG format. Plots are saved separately for each entry in a multifasta file
```
$ ./plot_skew.py
usage: plot_skew.py [-h] -i INPUT [-n NT] [-w WINDOW] [-s STEP]

Generate dinucleotide skew plot(s) of FASTA format DNA sequences in SVG format. Plots are saved separately for each entry in a multifasta file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Fasta (required)
  -n NT, --nt NT        dinucleotide (default: GC).
  -w WINDOW, --window WINDOW
                        window size (default: 1000)
  -s STEP, --step STEP  step size (default: 100)
```
Example: [Escherichia coli str. K-12 substr. MG1655, complete genome NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3?report=fasta)
```
$ ./plot_skew.py NC_000913.3.fasta # (output: NC_000913.3.svg)
```
![NC_000913 3](https://user-images.githubusercontent.com/58936715/180006455-c88b7461-6796-4517-aec3-b75281620441.png)

## vcf_to_sv_density_plot.py
