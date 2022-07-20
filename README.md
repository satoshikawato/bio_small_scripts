# bio_small_scripts
Miscellaneous scripts

## annotate_gff3.py
Add functional annotation to the 9th column of a gff3 file

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

## gb2nrfaa.py
Extract the longest isoforms of protein-coding genes from a NCBI RefSeq euaryotic genome assembly

## orfind.py

## plot_circular_genome.py
Generate genome diagram in SVG

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
$ ./plot_skew.py ./plot_skew.py NC_000913.3.fasta # (output: NC_000913.3.svg)
```
![NC_000913 3](https://user-images.githubusercontent.com/58936715/180006455-c88b7461-6796-4517-aec3-b75281620441.png)

## vcf_to_sv_density_plot.py
