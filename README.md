# bio_small_scripts
Small scripts to alleviate daily hassles in the dry lab.
See [here](https://qiita.com/satoshi_kawato) (in Japanese) for the details and motivation behind each piece.

2023-12-21: `plot_linear_genome.py` and `plot_circular_genome.py` will be merged into a single package, which will be made public soon.
## blast2bed.py
convert BLASTN/BLASTX/TBLASTX output into BED format
## Requirements
- Python3

```$ ./blast2bed.py
usage: blast2bed.py [-h] -i INPUT [-o OUTPUT] [-s SCORE] [-e EVALUE]

convert BLASTN/BLASTX/TBLASTX output into BED format

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        tab-separated blast output (required) with "-outfmt "[6|7]"
  -o OUTPUT, --output OUTPUT
                        output BED format file (default: stdout)
  -s SCORE, --score SCORE
                        score (default: 0)
  -e EVALUE, --evalue EVALUE
                        E-value threshold (default: 1e-30)

```

## gb2faa.py
Extract protein sequences from a genbank file downloaded from NCBI
### Requirements
- [Biopython](https://biopython.org/)
### Usage
```
$ ./gb2faa.py
usage: gb2faa.py [-h] -i INPUT [-o OUTPUT]

Extract protein sequences from a genbank file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        GenBank flat file format of the genomic sequence(s) (required)
  -o OUTPUT, --output OUTPUT
                        output fasta file (default: out.faa)
```

## annotate_gff3.py
Add functional annotation to the 9th column of a gff3 file
### Requirements
### Usage
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
Draw a dot plot based on a pairwise BLASTN/TBLASTX result
### Requirements
- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [SVGwrite](https://svgwrite.readthedocs.io/en/latest/)
### Usage
```
$ ./blast2dotplot.py
usage: blast2dotplot.py [-h] --input FILE

Draw a dot plot based on a pairwise BLASTN/TBLASTX result

optional arguments:
  -h, --help            show this help message and exit
  --input FILE, --in FILE, -i FILE
                        input BLASTN/TBLASTX result file in XML format (-outfmt 5)
```
![Query_1_Subject_1_BLASTN](https://user-images.githubusercontent.com/58936715/200174764-310d8112-2e50-42b6-9e31-b030a48648c2.svg)

## scaffold2contigs.py
### Requirements
- [Biopython](https://biopython.org/)
### Usage
```
$ ./scaffold2contigs.py
usage: scaffold2contigs.py [-h] --input INPUT [--output OUTPUT] [-d DIGIT]

Split scaffolds into contigs

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT, --in INPUT
                        Input FASTA file
  --output OUTPUT, -o OUTPUT, --out OUTPUT
                        output FASTA file (default: stdout)
  -d DIGIT, --digit DIGIT
                        number of digits for zero-padding (default:3)
```
```
$ less in.fa
>scaffold1
ACTGTGCATNNNNNNACGCTGCANnnNNCTGCAnnnCTGCAnnNNNNCTGCA
>scaffold2
ACGACGACGCGATAGAGnnnnnnAGACGAGAGNNNnnACGACGACG
```
```
$ ./scaffold2contigs.py -i in.fa
>scaffold1_001
ACTGTGCAT
>scaffold1_002
ACGCTGCA
>scaffold1_003
CTGCA
>scaffold1_004
CTGCA
>scaffold1_005
CTGCA
>scaffold2_001
ACGACGACGCGATAGAG
>scaffold2_002
AGACGAGAG
>scaffold2_003
ACGACGACG
```
## depth_alignment_breakpoint.py

## ddbj_to_gff3.py
### Requirements
- [Biopython](https://biopython.org/)
### Usage
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
### Requirements
- [Biopython](https://biopython.org/)
### Usage
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

##get_flanking_reads.py
Extract flanking reads from a BAM file
## Requirements
- [pysam](https://pysam.readthedocs.io/en/latest/)
```
$ ./get_flanking_reads.py
usage: get_flanking_reads.py [-h] --input FILE [--output FILE] -r REF [-w WINDOW] [-m MIN] [-f FLANK] [-p PRIME] [-q]

Extract flanking reads from a BAM file

optional arguments:
  -h, --help            show this help message and exit
  --input FILE, -i FILE, --in FILE
                        Input BAM file
  --output FILE, -o FILE, --out FILE, --output FILE
                        output txt file
  -r REF, --ref REF     reference entry name
  -w WINDOW, --window WINDOW
                        window (defalt:100)
  -m MIN, --min MIN     minimum outut read lenth threshold (defalt:200)
  -f FLANK, --flank FLANK
                        flanking bases (defalt:100)
  -p PRIME, --prime PRIME
                        5'/3'-end
  -q, --fastq           output fastq
```
## orfind.py
ORF prediction allowing CDS overlaps and GFF3 output
### Requirements
- [Biopython](https://biopython.org/)
### Usage
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
## plot_linear_genome.py
Generate plot in SVG
### Requirements
- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [SVGwrite](https://svgwrite.readthedocs.io/en/latest/)
### Usage
```
$ ./plot_linear_genome.py
usage: plot_linear_genome.py [-h] -i [INPUT ...] [-b [BLAST ...]] [-t TABLE] [-o OUTPUT] [-n NT] [-w WINDOW] [-s STEP]
                             [--separate_strands] [--show_gc] [--align_center] [--evalue EVALUE] [--bitscore BITSCORE]
                             [--identity IDENTITY]

Generate plot in SVG

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT ...], --input [INPUT ...]
                        genbank (required)
  -b [BLAST ...], --blast [BLAST ...]
                        input BLAST result file in tab-separated format (-outfmt 6 or 7) (optional)
  -t TABLE, --table TABLE
                        color table (optional)
  -o OUTPUT, --output OUTPUT
                        output prefix (default: diagram)
  -n NT, --nt NT        dinucleotide (default: GC).
  -w WINDOW, --window WINDOW
                        window size (default: 1000)
  -s STEP, --step STEP  step size (default: 100)
  --separate_strands    separate forward and reverse strands (default: False). Features of undefined strands are shown
                        on the forward strand.
  --show_gc             plot GC content below genome (default: False).
  --align_center        Align genomes to the center (default: False).
  --evalue EVALUE       evalue threshold (default=1e-2)
  --bitscore BITSCORE   bitscore threshold (default=50)
  --identity IDENTITY   identity threshold (default=0)
```
Example:[lymphocystis disease virus](https://en.wikipedia.org/wiki/Lymphocystis)

- [Lymphocystis disease virus 1, complete genome (NC_001824.1)](https://www.ncbi.nlm.nih.gov/nuccore/NC_001824.1)
- [Lymphocystis disease virus Sa isolate SA9, complete genome (NC_033423.1)](https://www.ncbi.nlm.nih.gov/nuccore/NC_033423.1)
- [Lymphocystis disease virus 2 LCDV-JP_Oita_2018 DNA, complete genome (LC534415.1)](https://www.ncbi.nlm.nih.gov/nuccore/LC534415.1)

```sh
$ tblastx -query NC_001824.fasta -subject NC_033423.fasta -outfmt 7 -out NC_001824_NC_033423.tblastx.out
$ tblastx -query NC_033423.fasta -subject LC534415.fasta -outfmt 7 -out NC_033423_LC534415.tblastx.out
```
```sh
./plot_linear_genome.py -i NC_001824.gb NC_033423.gb LC534415.gb -b NC_001824_NC_033423.tblastx.out NC_033423_LC534415.tblastx.out --separate_strands --align_center --evalue 1e-10 --bitscore 100 # output: out.svg
```
![out](https://user-images.githubusercontent.com/58936715/210162432-bb72385c-3942-4d83-bbd1-1b8ac9fca348.svg)


## plot_circular_genome.py
Generate genome diagram in SVG
### Requirements
- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [SVGwrite](https://svgwrite.readthedocs.io/en/latest/)
### Usage
```
$ ./plot_circular_genome.py -h
usage: plot_circular_genome.py [-h] -i INPUT [-t TABLE] [-n NT] [-w WINDOW] [-s STEP]

Generate genome diagram in SVG.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Genbank/DDBJ flatfile (required)
  -t TABLE, --table TABLE
                        color table (optional)
  -n NT, --nt NT        dinucleotide (default: GC).
  -w WINDOW, --window WINDOW
                        window size (default: 1000)
  -s STEP, --step STEP  step size (default: 100)
```
![NC_003225](https://user-images.githubusercontent.com/58936715/192824730-4fddc66c-8853-4ee7-83e6-58f945916411.svg)
![NZ_LR214945](https://user-images.githubusercontent.com/58936715/192824717-099c20e0-f976-431e-9881-f333ee10e597.svg)
## plot_skew.py
Generate dinucleotide skew plot(s) of FASTA format DNA sequences in SVG format. Plots are saved separately for each entry in a multifasta file
### Requirements
- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [Matplotlib](https://matplotlib.org/)
### Usage
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
![NC_000913 3](https://user-images.githubusercontent.com/58936715/192764584-f21297c2-8595-40c4-baae-20d6aab5575e.svg)
## vcf_to_sv_density_plot.py


## msa_to_txt.py
### Requirements
- [Biopython](https://biopython.org/)
### Usage
```sh
$ ./msa_to_txt.py
usage: aln_to_txt_wrap.py [-h] --input FILE [--output FILE] [-r REF] [-s START] [-e END] [-g GAP] [-w WRAP] [--gap_inclusive]

Convert FASTA-format multiple sequence alignment into a txt file. Assumes Courier New

optional arguments:
  -h, --help            show this help message and exit
  --input FILE, -i FILE, --in FILE
                        Input FASTA file
  --output FILE, -o FILE, --out FILE, --output FILE
                        output txt file
  -r REF, --ref REF     reference entry name
  -s START, --start START
                        start position
  -e END, --end END     end position
  -g GAP, --gap GAP     gap character (default: "-")
  -w WRAP, --wrap WRAP  line width (default: 100)
  --gap_inclusive       Gap inclusive (default: False).
```
```sh
$ less 16S.aligned.fasta

>NR_024570.1 Escherichia coli strain U 5/41 16S ribosomal RNA, partial sequence
---------AGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGG-AAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAG-CAC-AAAGAGGGGGACCTTAGGGC--------CTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCAACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCNGCGTGTATGAAGAAGGCCTTC-GGGTTGTAAAGTACTTTCAGCGGGGAGGAAG-GGAGTAAAGTTAATACCTTTGCTCATTGACGTTACC-CGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCA-GGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTT-GAGGCGTGGCTTCCGGANNTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAA-TGAATTGACGGGGGCC-GCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTT-CAGAGATGAGAATGTGCCT-----TCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGC-GGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAG-AATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACTTCGG-GAGGGCG----------------------------------------------------------------------------------
>NR_044682.2 Haemophilus influenzae strain 680 16S ribosomal RNA, partial sequence
A-ATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAACGGTAGCAGGAGAAAGCTTGCTTTCTTGCTGACGAGTGGCGGACGGGTGAGTAATGCTTGGG-AATCTGGCTTATGGAGGGGGATAACGACGGGAAACTGTCGCTAATACCGCGTATTATCGGAAG-ATG-AAAGTGCGGGACTGAGAGGC--------CGCATGCCATAGGATGAGCCCAAGTGGGATTAGGTAGTTGGTGGGGTAAATGCCTACCAAGCCTGCGATCTCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCGCNATGGGGGGAACCCTGACGCAGCCATGCCGCGTGAATGAAGAAGGCCTTC-GGGTTGTAAAGTTCTTTCGGTATTGAGGAAG-GTTGATGTGTTAATAGCACATCAAATTGACGTTAAA-TACAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGNGTGCGAGCGTTAATCGGAATAACTGGGCGTAAAGGGCACGCAGGCGGTTATTTAAGTGAGGTGTGAAAGCCCCGGGCTTAACCTGGGNATTGCATTTCAGACTGGGTAACTAGAGTACTTTAGGGAGGGGTAGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAATACCGAAGGCGAAGGCAGCCCCTTGGGAATGTACTGACGCTCA-TGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGCTGTCGATTTGGGGGTTGGGGTTT---AACTCTGGCACCCGTAGCTAACGTGATAAATCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAA-TGAATTGACGGGGGCCNGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTACTCTTGACATCCTAAGAAGAGCT-CAGAGATGAGCTTGTGCCT-----TCGGGAACTTAGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGC-GACTTGGTCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTNGGGATGACGTCAAGTCATCATGGCCCTTACGAGTAGGGCTACACACGTGCTACAATGGCGTATACAGAGGGAAGCGAAGCTGCGAGGTGGAGCGAATCTCATAAAGTACGTCTAAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGCGAATCAG-AATGTCGCGGTGAATACGTTCCCGGGCNTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGTACCAGAAGTAGATAGCTTAACCTTTT-GGAGGGCGTTTACCACGGTATGATTCATGACTGGGG-----------------------------------------------------
>NR_112116.2 Bacillus subtilis strain IAM 12118 16S ribosomal RNA, complete sequence
TTATCGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGG--ACAGATGGGAGCTTGCTCCCTGAT--GTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTGTAAGACTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATGGTTGTTTGAA-CCGCATGGTTCAAACATAAAAGGTGGCTTCGGCTACCACTTACAGATGGACCCGCGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCAACGATGCGTAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTTTTC-GGATCGTAAAGCTCTGTTGTTAGGGAAGAACAAGTACCGTTCGAATAGGGCGGTACCTTGACGGTACC-TAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGA-GGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTCAAA-GGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTCTGACAATCC-TAGAGATAGGACGTCCCCT-----TCGGGGGCAGAGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTGGATCTTAGTTGCCAGC--ATTCAGTTGGGCACTCTAAGGTGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACAGAACAAAGGGCAGCGAAACCGCGAGGTTAAGCCAATCCCACAAATCTGTTCTCAGTTCGGATCGCAGTCTGCAACTCGACTGCGTGAAGCTGGAATCGCTAGTAATCGCGGATCAG-CATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCTTTTAGGAGCCAGCCGCCGAAGGTGGGACAGATGATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTT
>NR_044761.1 Helicobacter pylori strain ATCC 43504 16S ribosomal RNA, partial sequence
TTTATGGAGAGTTTGATCCTGGCTCAGAGTGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAACGAT-GAAGCTTCTAGCTTGCTAGAGTGCTGATTAGTGGCGCACGGGTGAGTAACGCATAGGTCATGTGCCTCTTAGTTTGGGATAGCCATTGGAAACGATGATTAATACCAGATACTCCCTACGG-GGG---------------AAAGAT--------TTATCGCTAAGAGATCAGCCTATGTCCTATCAGCTTGTTGGTAAGGTAATGGCTTACCAAGGCTATGACGGGTATCCGGCCTGAGAGGGTGAACGGACACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTAGGGAATATTGCTCAATGGGGGAAACCCTGAAGCAGCAACGCCGCGTGGAGGATGAAGGTTTTA-GGATTGTAAACTCCTTTTGTTAGAGAAGATA--------------------------ATGACGGTATC-TAACGAATAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTACTCGGAATCACTGGGCGTAAAGAGCGCGTAGGCGGGATAGTCAGTCAGGTGTGAAATCCTATGGCTTAACCATAGAACTGCATTTGAAACTACTATTCTAGAGTGTGGGAGAGGTAGGTGGAATTCTTGGTGTAGGGGTAAAATCCGTAGAGATCAAGAGGAATACTCATTGCGAAGGCGACCTGCTGGAACATTACTGACGCTGATTGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGGATGCTAGTTGTTGGAGGGCTTAGTCTCTCCAGTAATGCAGCTAACGCATTAAGCATCCCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAA-GGAATAGACGGGGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGATACACGAAGAACCTTACCTAGGCTTGACATTGAGAGAATCCGC-TAGAAATAGTGGAGTGTCTAGCTTGCTAGACCTTGAAAACAGGTGCTGCACGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCCTTTCTTAGTTGCTAACAGGTTATGCTGAGAACTCTAAGGATACTGCCTCCG-TAAGGAGGAGGAAGGTGGGGACGACGTCAAGTCATCATGGCCCTTACGCCTAGGGCTACACACGTGCTACAATGGGGTGCACAAAGAGAAGCAATACTGTGAAGTGGAGCCAATCTT-CAAAACACCTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTGCATGAAGCTGGAATCGCTAGTAATCGCAAATCAGCCATGTTGCGGTGAATACGTTCCCGGGTCTTGTACTCACCGCCCGTCACACCATGGGAGTTGTGTTTGCCTTAAGTCAGGATGCTAAATT-------GGCTACTGCCCACGGCACACACAGCGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGTGAACCTGCGGCTGGATCACCTCCTT-
>NR_025900.1 Thermus aquaticus strain YT-1 16S ribosomal RNA, partial sequence
---------------------GCTCAGGGTGAACGCTGGCGGCGTGCCTAAGACATGCAAGTCGTGCGGG-CCGTGGGGTATCTCAC---------GGTCAGCGGCGGACGGGTGAGTAACGCGTGGGTGACCTACCCGGAAGAGGGGGACAACATGGGGAAACCCAGGCTAATCCCCCATGTGGACACATC-CTGTGGGGTGTGTTTAAAGGGTTT--------TGCCCGCTTCCGGATGGGCCCGCGTCCCATCAGCTAGTTGGTGGGGTAAGAGCCCACCAAGGCGACGACGGGTAGCCGGTCTGAGAGGACGGCCGGCCACAGGGGCACTGAGACACGGGCCCCACTCCTACGGGAGGCAGCAGTTAGGAATCTTCCGCAATGGGCGCAAGCCTGACGGAGCGACGCCGCTTGGAGGAGGAAGCCCTTC-GGGGTGTAAACTCCTGAACCCGGGACGAAAC--------CCCCGATGAGG----GGACTGACGGTACC--GGGGTAATAGCGCCGGCCAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGCGCGAGCGTTACCCGGATTTACTGGGCGTAAAGGGCGTGTAGGCGGCTTGGGGCGTCCCATGTGAAAGGCCACGGCTCAACCGTGGAGGAGCGTGGGATACGCTCAGGCTAGACGGTGGGAGAGGGTGGTGGAATTCCCGGAGTAGCGGTGAAATGCGCAGATACCGGGAGGAACGCCGATGGCGAAGGCAGCCACCTGGTCCACTCGTGACGCTGA-GGCGCGAAAGCGTGGGGAGCAAACCGGATTAGATACCCGGGTAGTCCACGCCCTAAACGATGCGCGCTAGGTCTCTGGG-------TTATCTGGGGGCCGAAGCTAACGCGTTAAGCGCGCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTCAAA-GGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGCCTTGACATGCTAGGGAACCTGGGTGAAAGCCTGGGGTGCCCCGCG-AGGGGAGCCCTAGCACAGGTGCTGCATGGCCGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTGCCGTTAGTTGCCAGCGGGTGAAGCCGGGCACTCTAACGGGACTGCCTGCG-AAAGCAGGAGGAAGGCGGGGACGACGTCTGGTCATCATGGCCCTTACGGCCTGGGCGACACACGTGCTACAATGCCCACTACAGAGCGAGGCGACCTGGCAACAGGGAGCGAATCGCAAAAAGGTGGGCGTAGTTCGGATTGGGGTCTGCAACCCGACCCCATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACGCCATGGGAGCGGGTTCTACCCGAAGTCGCCGGG--AGCCT----TAGGGCAGGCGCCGAGGGTAGGGCCCGTGACTGGGGCGAAGTCGTAACAAGGTAGCTGTACCG--------------------------
>NR_041751.1 Mycoplasma pneumoniae FH strain ATCC 15531 16S ribosomal RNA, partial sequence
-----------------------------TTAACGCTGGCGGCATGCCTAATACATGCAAGTCGATCGAA-AGTAGTAATACT---------------TTAGAGGCGAACGGGTGAGTAACACGTATCCAATCTACCTTATAATGGGGGATAACTAGTTGAAAGACTAGCTAATACCGCATAAGAACTTTGGTTCGCATGAATCAAAGTTGAAAGGACCTGCAAGGGTTCGTTATTTGATGAGGGTGCGCCATATCAGCTAGTTGGTGGGGTAACGGCCTACCAAGGCAATGACGTGTAGCTATGCTGAGAAGTAGAATAGCCACAATGGGACTGAGACACGGCCCATACTCCTACGGGAGGCAGCAGTAGGGAATTTTTCACAATGAGCGAAAGCTTGATGGAGCAATGCCGCGTGAACGATGAAGGTCTTTAAGATTGTAAAGTTCTTTTATTTGGGAAGAAT-GACTTTAGCAGGTAATGGCTAGAGTTTGACTGTACCATTTTGAATAAGTGACGACTAACTATGTGCCAGCAGTCGCGGTAATACATAGGTCGCAAGCGTTATCCGGATTTATTGGGCGTAAAGCAAGCGCAGGCGGATTGAAAAGTCTGGTGTTAAAGGCAGCTGCTTAACAGTTGTA-TGCATTGGAAACTATTAATCTAGAGTGTGGTAGGGAGTTTTGGAATTTCATGTGGAGCGGTGAAATGCGTAGATATATGAAGGAACACCAGTGGCGAAGGCGAAAACTTAGGCCATTACTGACGCTTA-GGCTTGAAAGTGTGGGGAGCAAATAGGATTAGATACCCTAGTAGTCCACACCGTAAACGATAGATACTAGCTGTCGGGGCG----ATCCCCTCGGTAGTGAAGTTAACACATTAAGTATCTCGCCTGGGTAGTACATTCGCAAGAATGAAACTCAAACGGAATTGACGGGGACCCGCACAAGTGGTGGAGCATGTTGCTTAATTCGACGGTACACGAAAAACCTTACCTAGACTTGACATCCTTGGCAAAGTTATGGAAACATAATGGAGGTT----------AACCGAGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCGTTAGTTAC----------------ATTGTCTAGCGAGACTGCTAATG-CAAATTGGAGGAAGGAAGGGATGACGTCAAATCATCATGCCCCTTATGTCTAGGGCTGCAAACGTGCTACAATGGCCAATACAAACAGTCGCCAGCTTGTAAAAGTGAGCAAATCTG-TAAAGTTGGTCTCAGTTCGGATTGAGGGCTGCAATTCGTCCTCATGAAGTCGGAATCACTAGTAATCGCGAATCAGCTATGTCGCGGTGAATACGTTCTCGGGTCTTGTACACACCGCCCGTCAAACTATGAAAGCTGGTAATATTTAAAAACGTGTTGCTAACCATTA-GGAAGCGCATGTCAAGGATAGCACCGGTGATTGGAGTTAAGTCGTAACAAGGTACCCCTACGAGAACGTGGGGGTGGATCACCTCCTTT
```
```
$ msa_to_txt.py -i 16S.aligned.fasta -o 16S.aligned.txt
$ less 16S.aligned.txt

                                          ......  ..************  **.*** ************. **            * ..  .                  
NR_024570.1        1 ---------AGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTTGCTGACG   91
NR_044682.2        1 A-ATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAACGGTAGCAGGAGAAAGCTTGCTTTCTTGCTGACG   99
NR_112116.2        1 TTATCGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGG--ACAGATGGGAGCTTGCTCCCTGAT--GTT   96
NR_044761.1        1 TTTATGGAGAGTTTGATCCTGGCTCAGAGTGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAACGAT-GAAGCTTCTAGCTTGCTAGAGTGCTGATT   99
NR_025900.1        1 ---------------------GCTCAGGGTGAACGCTGGCGGCGTGCCTAAGACATGCAAGTCGTGCGGG-CCGTGGGGTATCTCAC---------GGTC   69
NR_041751.1        1 -----------------------------TTAACGCTGGCGGCATGCCTAATACATGCAAGTCGATCGAA-AGTAGTAATACT---------------TT   55

                     ** **** ************  . * ..  * .* .*.    .   ****.*.*    .****.    ..****.**  .*     .        .         
NR_024570.1       92 AGTGGCGGACGGGTGAGTAATGTCTGGG-AAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAG-CAC-AAA  188
NR_044682.2      100 AGTGGCGGACGGGTGAGTAATGCTTGGG-AATCTGGCTTATGGAGGGGGATAACGACGGGAAACTGTCGCTAATACCGCGTATTATCGGAAG-ATG-AAA  196
NR_112116.2       97 AGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTGTAAGACTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATGGTTGTTTGAA-CCGCATG  195
NR_044761.1      100 AGTGGCGCACGGGTGAGTAACGCATAGGTCATGTGCCTCTTAGTTTGGGATAGCCATTGGAAACGATGATTAATACCAGATACTCCCTACGG-GGG----  194
NR_025900.1       70 AGCGGCGGACGGGTGAGTAACGCGTGGGTGACCTACCCGGAAGAGGGGGACAACATGGGGAAACCCAGGCTAATCCCCCATGTGGACACATC-CTGTGGG  168
NR_041751.1       56 AGAGGCGAACGGGTGAGTAACACGTATCCAATCTACCTTATAATGGGGGATAACTAGTTGAAAGACTAGCTAATACCGCATAAGAACTTTGGTTCGCATG  155

                                   .               ..     ***. ...   .    ** **.*.**.***. *****  **  ***.**.* . **    **.     
NR_024570.1      189 GAGGGGGACCTTAGGGC--------CTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAG  280
NR_044682.2      197 GTGCGGGACTGAGAGGC--------CGCATGCCATAGGATGAGCCCAAGTGGGATTAGGTAGTTGGTGGGGTAAATGCCTACCAAGCCTGCGATCTCTAG  288
NR_112116.2      196 GTTCAAACATAAAAGGTGGCTTCGGCTACCACTTACAGATGGACCCGCGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCAACGATGCGTAG  295
NR_044761.1      195 -----------AAAGAT--------TTATCGCTAAGAGATCAGCCTATGTCCTATCAGCTTGTTGGTAAGGTAATGGCTTACCAAGGCTATGACGGGTAT  275
NR_025900.1      169 GTGTGTTTAAAGGGTTT--------TGCCCGCTTCCGGATGGGCCCGCGTCCCATCAGCTAGTTGGTGGGGTAAGAGCCCACCAAGGCGACGACGGGTAG  260
NR_041751.1      156 AATCAAAGTTGAAAGGACCTGCAAGGGTTCGTTATTTGATGAGGGTGCGCCATATCAGCTAGTTGGTGGGGTAACGGCCTACCAAGGCAATGACGTGTAG  255

                     * .  ******.*  *. . *..*** .** ************ **. ********************* .***** ** * *.***.. * ** ..***     
NR_024570.1      281 CTGGTCTGAGAGGATGACCAGCAACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGA  380
NR_044682.2      289 CTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCGCNATGGGGGGAACCCTGA  388
NR_112116.2      296 CCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGA  395
NR_044761.1      276 CCGGCCTGAGAGGGTGAACGGACACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTAGGGAATATTGCTCAATGGGGGAAACCCTGA  375
NR_025900.1      261 CCGGTCTGAGAGGACGGCCGGCCACAGGGGCACTGAGACACGGGCCCCACTCCTACGGGAGGCAGCAGTTAGGAATCTTCCGCAATGGGCGCAAGCCTGA  360
NR_041751.1      256 CTATGCTGAGAAGTAGAATAGCCACAATGGGACTGAGACACGGCCCATACTCCTACGGGAGGCAGCAGTAGGGAATTTTTCACAATGAGCGAAAGCTTGA  355

                      * *** * **.**.** . ** ****.  **  .* ..***** . ** .     . .. .*.                            ****. **     
NR_024570.1      381 TGCAGCCATGCNGCGTGTATGAAGAAGGCCTTC-GGGTTGTAAAGTACTTTCAGCGGGGAGGAAG-GGAGTAAAGTTAATACCTTTGCTCATTGACGTTA  478
NR_044682.2      389 CGCAGCCATGCCGCGTGAATGAAGAAGGCCTTC-GGGTTGTAAAGTTCTTTCGGTATTGAGGAAG-GTTGATGTGTTAATAGCACATCAAATTGACGTTA  486
NR_112116.2      396 CGGAGCAACGCCGCGTGAGTGATGAAGGTTTTC-GGATCGTAAAGCTCTGTTGTTAGGGAAGAACAAGTACCGTTCGAATAGGGCGGTACCTTGACGGTA  494
NR_044761.1      376 AGCAGCAACGCCGCGTGGAGGATGAAGGTTTTA-GGATTGTAAACTCCTTTTGTTAGAGAAGATA--------------------------ATGACGGTA  448
NR_025900.1      361 CGGAGCGACGCCGCTTGGAGGAGGAAGCCCTTC-GGGGTGTAAACTCCTGAACCCGGGACGAAAC--------CCCCGATGAGG----GGACTGACGGTA  447
NR_041751.1      356 TGGAGCAATGCCGCGTGAACGATGAAGGTCTTTAAGATTGTAAAGTTCTTTTATTTGGGAAGAAT-GACTTTAGCAGGTAATGGCTAGAGTTTGACTGTA  454

                      .     ... .**.  **.*.**** .**********.***********. **.  ** ******.  **** * * *********** .. .* ****     
NR_024570.1      479 CC-CGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGC  577
NR_044682.2      487 AA-TACAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGNGTGCGAGCGTTAATCGGAATAACTGGGCGTAAAGGGCACGCAGGC  585
NR_112116.2      495 CC-TAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTCGCAGGC  593
NR_044761.1      449 TC-TAACGAATAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTACTCGGAATCACTGGGCGTAAAGAGCGCGTAGGC  547
NR_025900.1      448 CC--GGGGTAATAGCGCCGGCCAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGCGCGAGCGTTACCCGGATTTACTGGGCGTAAAGGGCGTGTAGGC  545
NR_041751.1      455 CCATTTTGAATAAGTGACGACTAACTATGTGCCAGCAGTCGCGGTAATACATAGGTCGCAAGCGTTATCCGGATTTATTGGGCGTAAAGCAAGCGCAGGC  554

                     ** ..    .**. . ***.***  *   .*** ***.   *    .*.*  .* **.      ** **..   . ** *    .* *****..  *.*.     
NR_024570.1      578 GGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGT  677
NR_044682.2      586 GGTTATTTAAGTGAGGTGTGAAAGCCCCGGGCTTAACCTGGGNATTGCATTTCAGACTGGGTAACTAGAGTACTTTAGGGAGGGGTAGAATTCCACGTGT  685
NR_112116.2      594 GGTTTCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGT  693
NR_044761.1      548 GGGATAGTCAGTCAGGTGTGAAATCCTATGGCTTAACCATAGAACTGCATTTGAAACTACTATTCTAGAGTGTGGGAGAGGTAGGTGGAATTCTTGGTGT  647
NR_025900.1      546 GGCTTGGGGCGTCCCATGTGAAAGGCCACGGCTCAACCGTGGAGGAGCGTGGGATACGCTCAGGCTAGACGGTGGGAGAGGGTGGTGGAATTCCCGGAGT  645
NR_041751.1      555 GGATTGAAAAGTCTGGTGTTAAAGGCAGCTGCTTAACAGTTGTA-TGCATTGGAAACTATTAATCTAGAGTGTGGTAGGGAGTTTTGGAATTTCATGTGG  653

                     **.***.****.**.*** *.  ..***** .*.  ..********  .    *.*      ..******* *  * ..*****.*********.**..*     
NR_024570.1      678 AGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCA-GGTGCGAAAGCGTGGGGAGCAAACAG  776
NR_044682.2      686 AGCGGTGAAATGCGTAGAGATGTGGAGGAATACCGAAGGCGAAGGCAGCCCCTTGGGAATGTACTGACGCTCA-TGTGCGAAAGCGTGGGGAGCAAACAG  784
NR_112116.2      694 AGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGA-GGAGCGAAAGCGTGGGGAGCGAACAG  792
NR_044761.1      648 AGGGGTAAAATCCGTAGAGATCAAGAGGAATACTCATTGCGAAGGCGACCTGCTGGAACATTACTGACGCTGATTGCGCGAAAGCGTGGGGAGCAAACAG  747
NR_025900.1      646 AGCGGTGAAATGCGCAGATACCGGGAGGAACGCCGATGGCGAAGGCAGCCACCTGGTCCACTCGTGACGCTGA-GGCGCGAAAGCGTGGGGAGCAAACCG  744
NR_041751.1      654 AGCGGTGAAATGCGTAGATATATGAAGGAACACCAGTGGCGAAGGCGAAAACTTAGGCCATTACTGACGCTTA-GGCTTGAAAGTGTGGGGAGCAAATAG  752

                     ************..*********.*. ******.*.    .* .  .   * *            .       * *. ****.. .***.    .*****     
NR_024570.1      777 GATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTT-GAGGCGTGGCTTCCGGANNTAACGCGTTAAGTCGACCGCCT  875
NR_044682.2      785 GATTAGATACCCTGGTAGTCCACGCTGTAAACGCTGTCGATTTGGGGGTTGGGGTTT---AACTCTGGCACCCGTAGCTAACGTGATAAATCGACCGCCT  881
NR_112116.2      793 GATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCT  892
NR_044761.1      748 GATTAGATACCCTGGTAGTCCACGCCCTAAACGATGGATGCTAGTTGTTGGAGGGCTTAGTCTCTCCAGTAATGCAGCTAACGCATTAAGCATCCCGCCT  847
NR_025900.1      745 GATTAGATACCCGGGTAGTCCACGCCCTAAACGATGCGCGCTAGGTCTCTGGG-------TTATCTGGGGGCCGAAGCTAACGCGTTAAGCGCGCCGCCT  837
NR_041751.1      753 GATTAGATACCCTAGTAGTCCACACCGTAAACGATAGATACTAGCTGTCGGGGCG----ATCCCCTCGGTAGTGAAGTTAACACATTAAGTATCTCGCCT  848

                     ***.*****.. ******  * *********  ****.******* ** *******.************.*.********* *  ** ****.*******     
NR_024570.1      876 GGGGAGTACGGCCGCAAGGTTAAAACTCAAA-TGAATTGACGGGGGCC-GCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTA  973
NR_044682.2      882 GGGGAGTACGGCCGCAAGGTTAAAACTCAAA-TGAATTGACGGGGGCCNGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTA  980
NR_112116.2      893 GGGGAGTACGGTCGCAAGACTGAAACTCAAA-GGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTA  991
NR_044761.1      848 GGGGAGTACGGTCGCAAGATTAAAACTCAAA-GGAATAGACGGGGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGATACACGAAGAACCTTA  946
NR_025900.1      838 GGGGAGTACGGCCGCAAGGCTGAAACTCAAA-GGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTA  936
NR_041751.1      849 GGGTAGTACATTCGCAAGAATGAAACTCAAACGGAATTGACGGGGACCCGCACAAGTGGTGGAGCATGTTGCTTAATTCGACGGTACACGAAAAACCTTA  948

                     **  . ******** .   * .        ** *       . . ..         .. *      ******* ****.** .*************** *     
NR_024570.1      974 CCTGGTCTTGACATCCACGGAAGTTTT-CAGAGATGAGAATGTGCCT-----TCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTG 1067
NR_044682.2      981 CCTACTCTTGACATCCTAAGAAGAGCT-CAGAGATGAGCTTGTGCCT-----TCGGGAACTTAGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTG 1074
NR_112116.2      992 CCAGGTCTTGACATCCTCTGACAATCC-TAGAGATAGGACGTCCCCT-----TCGGGGGCAGAGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCG 1085
NR_044761.1      947 CCTAGGCTTGACATTGAGAGAATCCGC-TAGAAATAGTGGAGTGTCTAGCTTGCTAGACCTTGAAAACAGGTGCTGCACGGCTGTCGTCAGCTCGTGTCG 1045
NR_025900.1      937 CCAGGCCTTGACATGCTAGGGAACCTGGGTGAAAGCCTGGGGTGCCCCGCG-AGGGGAGCCCTAGCACAGGTGCTGCATGGCCGTCGTCAGCTCGTGTCG 1035
NR_041751.1      949 CCTAGACTTGACATCCTTGGCAAAGTTATGGAAACATAATGGAGGTT----------AACCGAGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCG 1038

                     *** *********************************      ** ***.* . .       .  . . ..** *. * .*****.   *  **   ***     
NR_024570.1     1068 TGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGC-GGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGA 1166
NR_044682.2     1075 TGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGC-GACTTGGTCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGA 1173
NR_112116.2     1086 TGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTGGATCTTAGTTGCCAGC--ATTCAGTTGGGCACTCTAAGGTGACTGCCGGTGACAAACCGGA 1183
NR_044761.1     1046 TGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCCTTTCTTAGTTGCTAACAGGTTATGCTGAGAACTCTAAGGATACTGCCTCCG-TAAGGAGGA 1144
NR_025900.1     1036 TGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTGCCGTTAGTTGCCAGCGGGTGAAGCCGGGCACTCTAACGGGACTGCCTGCG-AAAGCAGGA 1134
NR_041751.1     1039 TGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCGTTAGTTAC----------------ATTGTCTAGCGAGACTGCTAATG-CAAATTGGA 1121

                     ******  **** ******.. ******** ****** * .  ****..**.*************.     *** *. *  ** *    *  *     **     
NR_024570.1     1167 GGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAG 1266
NR_044682.2     1174 GGAAGGTNGGGATGACGTCAAGTCATCATGGCCCTTACGAGTAGGGCTACACACGTGCTACAATGGCGTATACAGAGGGAAGCGAAGCTGCGAGGTGGAG 1273
NR_112116.2     1184 GGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACAGAACAAAGGGCAGCGAAACCGCGAGGTTAAG 1283
NR_044761.1     1145 GGAAGGTGGGGACGACGTCAAGTCATCATGGCCCTTACGCCTAGGGCTACACACGTGCTACAATGGGGTGCACAAAGAGAAGCAATACTGTGAAGTGGAG 1244
NR_025900.1     1135 GGAAGGCGGGGACGACGTCTGGTCATCATGGCCCTTACGGCCTGGGCGACACACGTGCTACAATGCCCACTACAGAGCGAGGCGACCTGGCAACAGGGAG 1234
NR_041751.1     1122 GGAAGGAAGGGATGACGTCAAATCATCATGCCCCTTATGTCTAGGGCTGCAAACGTGCTACAATGGCCAATACAAACAGTCGCCAGCTTGTAAAAGTGAG 1221

                     * .*.*    ***     .*  *** *****.*  * ******..** *  *.*****  ******.**********.. *****  *** ..*******     
NR_024570.1     1267 CGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAG-AATGCCACGGTGAA 1365
NR_044682.2     1274 CGAATCTCATAAAGTACGTCTAAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGCGAATCAG-AATGTCGCGGTGAA 1372
NR_112116.2     1284 CCAATCCCACAAATCTGTTCTCAGTTCGGATCGCAGTCTGCAACTCGACTGCGTGAAGCTGGAATCGCTAGTAATCGCGGATCAG-CATGCCGCGGTGAA 1382
NR_044761.1     1245 CCAATCTT-CAAAACACCTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTGCATGAAGCTGGAATCGCTAGTAATCGCAAATCAGCCATGTTGCGGTGAA 1343
NR_025900.1     1235 CGAATCGCAAAAAGGTGGGCGTAGTTCGGATTGGGGTCTGCAACCCGACCCCATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCCATGCCGCGGTGAA 1334
NR_041751.1     1222 CAAATCTG-TAAAGTTGGTCTCAGTTCGGATTGAGGGCTGCAATTCGTCCTCATGAAGTCGGAATCACTAGTAATCGCGAATCAGCTATGTCGCGGTGAA 1320

                     *******.**** .******.************..*.*.* .**  ...   .    **.. .      .*.. .       ...      . . .         
NR_024570.1     1366 TACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACTTCGG-GAGGGCG-------------- 1450
NR_044682.2     1373 TACGTTCCCGGGCNTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGTACCAGAAGTAGATAGCTTAACCTTTT-GGAGGGCGTTTACCACGGTAT 1471
NR_112116.2     1383 TACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCTTTTAGGAGCCAGCCGCCGAAGGTGG 1482
NR_044761.1     1344 TACGTTCCCGGGTCTTGTACTCACCGCCCGTCACACCATGGGAGTTGTGTTTGCCTTAAGTCAGGATGCTAAATT-------GGCTACTGCCCACGGCAC 1436
NR_025900.1     1335 TACGTTCCCGGGCCTTGTACACACCGCCCGTCACGCCATGGGAGCGGGTTCTACCCGAAGTCGCCGGG--AGCCT----TAGGGCAGGCGCCGAGGGTAG 1428
NR_041751.1     1321 TACGTTCTCGGGTCTTGTACACACCGCCCGTCAAACTATGAAAGCTGGTAATATTTAAAAACGTGTTGCTAACCATTA-GGAAGCGCATGTCAAGGATAG 1419

                            .. ... .                                                          
NR_024570.1     1451 -------------------------------------------------------------------- 1451
NR_044682.2     1472 GATTCATGACTGGGG----------------------------------------------------- 1486
NR_112116.2     1483 GACAGATGATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTT 1550
NR_044761.1     1437 ACACAGCGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGTGAACCTGCGGCTGGATCACCTCCTT- 1503
NR_025900.1     1429 GGCCCGTGACTGGGGCGAAGTCGTAACAAGGTAGCTGTACCG-------------------------- 1470
NR_041751.1     1420 CACCGGTGATTGGAGTTAAGTCGTAACAAGGTACCCCTACGAGAACGTGGGGGTGGATCACCTCCTTT 1487
```

