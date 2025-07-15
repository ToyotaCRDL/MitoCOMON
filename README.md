## MitoCOMON
MitoCOMON is a pipeline for constructing complete sequence(s) of mitochondrial DNA from several amplicons sequenced by long reads. It consists of two modules: MitoCOMON design is for designing primer sets to amplify mitochondrial DNA as several amplicons; MitoCOMON assembly is for assemlying the long reads of amplicons to build complete mitochondrial DNA sequence(s).

## Table of contents
- [Installation](#installation)
- [Quick start](#quick-start)
- [MitoCOMON design](#mitocomon-design)
  - [File preparation](#file-preparation)
  - [Prepare NCBI taxonomy ID](#prepare-ncbi-taxonomy-id)
  - [Run](#run)
- [MitoCOMON assembly](#mitocomon-assembly)
  - [Basic (Mammal sample)](#basic-mammal-sample)
  - [Basic (Bird sample)](#basic-bird-sample)
  - [Custom amplicon](#custom-amplicon)
  - [Test run](#test-run)
- [Options](#options)
  - [MitoCOMON design](#mitocomon-design-1)
  - [MitoCOMON assembly](#mitocomon-assembly-1)
- [License](#license)
- [Citation](#citation)


## Installation
```
(1) Clone this repository.
$ git clone https://github2.cae.tytlabs.co.jp/e1828/MitoCOMON.git

(2) Create an environment for running MitoCOMON with installing the dependent packages (You can use mamba instead of conda to finish faster).
$ conda create -n MitoCOMON3 python biopython cutadapt chopper minimap2 miniasm medaka mitos trnascan-se taxonkit seqkit mafft primer3

(3) Create another environment for running PrimerProspector for MitoCOMON design as it requires Python2.
$ conda create -n primerprospector python=2.7 primerprospector
```

For using MitoCOMON design, the main script of Primerprospector have to be modified. Modify the original script analyze_primers.py at Line# 586 (Add int()).
```
# analyze_primers.py is typically located at "path/to/anaconda_or_miniconda/envs/primerprospector/lib/python2.7/site-packages/primerprospector/analyze_primers.py"

return int(max(counts_all_max))
```

## Quick start
- MitoCOMON design (A pipeline designing primer sets from mitochondrial sequences.)
```
mitocomon_design.py -i mitochondria_database_seq.fasta -p file_prefix -t taxonomy_id#(e.g. 40674 for Mammalia) -d ./output_directory -e primerprospector
```

- MitoCOMON assembly (An assembling pipeline of long-read amplicons.)
```
mitocomon_assembly.py -i Long_read_file.fastq -p Prefix_of_output -m 300 -o path/to/output -a [mammalia|bird|custom] -t thread#
```

## MitoCOMON design
MitoCOMON design is a module to design a primer sets for amplifying whole mtDNA as several overlapping amplicons. It aligns mtDNA sequences of species in the taxonomy group of interest, detect and obtain the sequences of conserved regions, and filters the sequences by their thermodynamic parameter and specificity.

### File preparation
Mitochondria database sequences can be downloaded from RefSeq FTP server. Although the following shows the example of downloading mitochondrion.1.1.genomic.fna.gz, you can download the latest version if available.
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.1.genomic.fna.gz
gunzip mitochondrion.1.1.genomic.fna.gz
```

For using taxonkit to extract the sequences of the species in the target taxonomic group, download the file required for conversion between RefSeq ID and NCBI Taxonomy ID.
```
wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
```

### Prepare NCBI taxonomy ID
Then, prepare the NCBI's taxonomy ID of the taxonomic group of your interest. For example, taxonomy ID of mammalia is 40674, aves is 8782, and teleostei is 32443. The ID can be searched at [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy).

### Run
Using the prepared mtDNA sequence file and taxonomy ID, run the following command.
```
mitocomon_design.py -i /path/to/mitochondrion.1.1.genomic.fna -p file_prefix -t taxonomy_id# -d /path/to/nucl_gb.accession2taxid.gz -o ./output_directory -e primerprospector
```

The list of primer candidates is at ./output_directory/03_primer_candidates/file_prefix_bestregions_seq_primer3_ppout_candidates.txt. For construction of a primer set, candidate primer pairs are highly recommended to be tested by wet experiments.


## MitoCOMON assembly
MitoCOMON assembly is a module to assemble the long reads that were obtained by sequencing the amplicons amplified using the designed primers. After filtering the reads with their length and quality, it clusters and obtain the consensus sequence of each amplicon, and then assembles the sequences utilizing the overlapping sequences.

### Basic (Mammal sample)
```
mitocomon_assembly.py -i Long_read_file.fastq -p Prefix_of_output -m 300 -o path/to/output -a mammalia -t thread#
```

### Basic (Bird sample)
```
mitocomon_assembly.py -i Long_read_file.fastq -p Prefix_of_output -m 300 -o path/to/output -a bird -t thread#
```

### Custom amplicon
If you used primers designed by MitoCOMON design or other methods, provide a tab delimited custom amplicon file which lists the name of each amplicon, primer sequences to amplify the amplicon, and expected length of the amplicon.

Example (amplicon_file/mammal.txt)
```
BCM1	AAAGCAARGCACTGAAAATG	GGTTCGAWTCCTTCCTTTCTT	6800
BCM2	AAAGAGTTACTTTGATAGAGTAAATNA	TARTYTAATGAGTCGAAATCAYTT	6100
BCM3	AAGAAAGGAAGGAWTCGAACC	ATTACTTTTATTTGGAGTTGCACC	4300
BCM4	AARTGATTTCGACTCATTARAYTA	GTTTGCTGAAGATGGCGGTATATAGRC	6500
```

For running MitoCOMON assembly, provide the path custom amplicon file.
```
mitocomon_assembly.py -i Long_read_file.fastq -p Prefix_of_output -m 300 -o path/to/output -a custom -c custom_amplicon_file -t thread#
```

### Test run
Small fastq data of cattle mtDNA is available for testing the assembly. Expected to find fasta file of complete cattle mtDNA sequence (16340 bp) in the output folder.
```
mitocomon_assembly.py -i path/to/MitoCOMON/testdata/cattle_testdata.fastq -p mitocomon_test -m 300 -o mitocomon_test -a mammalia -t 36
```


## Options
### MitoCOMON design
```
usage: mitocomon_design.py [-h] -i INPUT [-o OUTPUTFOLDER] -p PREFIX -t TAXID
                           -d DATADIR [-z TAXKITOPTION] [-s ENTROPYTHRESHOLD]
                           [-b BASEPOSREFSEQID] -e PPENV

Designing primer sets for mitocommon.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Fasta file of Mitochondria sequences. Recommended to
                        download from RefSeq release. [Required]
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Save the results in the specified outputfolder.
                        Default = current working directory
  -p PREFIX, --prefix PREFIX
                        Prefix of the output file. [Required]
  -t TAXID, --taxid TAXID
                        Taxon ID of the organism of your interest (e.g.
                        Mammalia, Aves, etc.). You can search it at NCBI
                        Taxonomy or use taxonkit name2taxid. [Required]
  -d DATADIR, --datadir DATADIR
                        Path of nucl_gb.accession2taxid. Download from https:/
                        /ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb
                        .accession2taxid.gz [Required]
  -z TAXKITOPTION, --taxkitoption TAXKITOPTION
                        Additional options for running taxonkit. For example,
                        if you keep taxdump.tar.gz in the specific folder, you
                        can add "--data-dir ./PATH/TO/TAXDUMP.tar.gz".
                        Optional.
  -s ENTROPYTHRESHOLD, --entropythreshold ENTROPYTHRESHOLD
                        Threshold of shannon entropy score for determining
                        conserved bases. Default is 1.80. Optional.
  -b BASEPOSREFSEQID, --baseposRefSeqid BASEPOSREFSEQID
                        RefSeq ID of the entry you want to use for showing the
                        corresponding position of primer sequences. If not
                        set, the first RefSeq ID in
                        01_sequence_splitting/tmp/target_refseqid.txt will be
                        used. Optional.
  -e PPENV, --ppenv PPENV
                        The name of Python 2 environment where
                        PrimerProspector is installed. [Required]
```

### MitoCOMON assembly
```
usage: mitocomon_assembly.py [-h] -i INPUT [-o OUTPUTFOLDER] -p PREFIX
                             [-a AMPLICON] [-c CUSTOM_AMPLICON] [-t THREADS]
                             [-m MAXR] [--medaka_model MEDAKA_MODEL]
                             [-ldc LENGTH_DIFF_CONSENSUS]

mitocomon: A pipeline for construction of complete mitochondrial sequences
from long PCR fragment sequences

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Read file in fastq format. [Required]
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Save the results in the specified outputfolder.
                        Default = current working directory
  -p PREFIX, --prefix PREFIX
                        Prefix of the output file. [Required]
  -a AMPLICON, --amplicon AMPLICON
                        Amplicons to analyze. When using for amplicons with
                        custom primers, designate custom and use -c option.
                        Choose from mammal, bird, custom. Default = mammal
  -c CUSTOM_AMPLICON, --custom_amplicon CUSTOM_AMPLICON
                        File of amplicons with their names and primer pair
                        sequences.
  -t THREADS, --threads THREADS
                        Number of threads used for each procedure. Default =
                        1.
  -m MAXR, --maxr MAXR  Max number of reads used in the amplicon_sorter step.
                        Larger number recommended when some fragments showed
                        low concentration, but would take longer calculation
                        time. Default = 1000.
  --medaka_model MEDAKA_MODEL
                        The model to use in Medaka. Use if the automatic model
                        detection of Medaka does not work. Check the model
                        list in the help message of medaka_consensus.
  -ldc LENGTH_DIFF_CONSENSUS, --length_diff_consensus LENGTH_DIFF_CONSENSUS
                        Length difference consensus parameter for
                        amplicon_sorter. Decrease if you want to avoid nesting
                        of shorter amplicon sequence. Default = 40.
```

## License
See LICENSE.md.

## Citation
Under review.