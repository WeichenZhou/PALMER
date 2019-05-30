# PALMER

Pre-mAsking Long reads for Mobile Element inseRtion

* PALMER detects non-reference MEI events (LINE, Alu and SVA) and other insertions, by using the indexed reference-aligned BAM files from long-read technology as inputs. It uses the track from Repeatmasker (https://www.girinst.org/) to mask the portions of reads that aligned to these repeats, defines the significant characteristics of MEIs (TSD motifs, 5' inverted sequence, 3' transduction sequence, polyA-tail), and reports sequences for each insertion event.
* The ideal structure of an MEI event should be 5’-TSD-(5'inverted)-MEI-polyA-(TransD-polyA)-TSD-3’. 


Required resources:
```
  samtools/1.3.1  https://github.com/samtools/samtools
  ncbi-blast++/2.4.0  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

## Getting started

Download and Install
```
git clone https://github.com/mills-lab/PALMER.git
cd PALMER
make
```

Parameters
```
Usage:

--input
         aligned long-read sequencing BAM file with directory path

--workdir
         the user's working directory

--ref_ver (options: hg19, GRCh37, GRCh38 or other)
         reference genome used for the aligned file ('other' option for the cusmized genome out of hg19, GRCh37 or GRCh38)

--ref_fa
         indexed fasta file of reference genome fasta file with directory path used for the aligned bam file (wrong reference will cause error infromation)

--type (options: LINE, ALU, SVA, or CUSTOMIZED (if you want to setup your costomized sequence))
         type of MEIs or other kind of insertion to detect

--chr (default: ALL (for whole genome); options: chromosome1, chromosome2, ...chromosomeY)
         chromosome name for PALMER to run (if running for whole genome, don't need to assign). !!The chromosome names should be consistent with the ones in reference genome version!! e.g. for GRCh37, to run PALMER on chromosome1, the option should be '1', while for GRCh38 it should be 'chr1'

--custom_seq (default:no input)
         .fasta file with directory path to customize your insertion finding

--custom_index (default:no input; if you have both '--ref_ver other' and '--type LINE/ALU/SVA', you must give PALMER a index file (format: "CHR'	'START'	'END'	'MEI_NAME'
'" for each MEI to be masked in each line) for masking module; if you have --custom_seq parameter without --custom_index, PALMER will work without masking step)
         index file with directory path to mask the genome for your insertion finding

--TSD_finding (Fixed:TRUE for all MEIs ,or default: FALSE for CUSTOMIZED insertion)
         whether to run TSD motif finding module for your insertion calling

--output (default: output)
         prefix of output file
```

Example
```
./PALMER --input $DirPath/NA12878.washu.alignment_hs37d5.1.bam --workdir $DirPath/chr1.line.0406/ --ref_ver GRCh37 --output test --type LINE --chr chr3 --ref_fa $DirPath/hs37d5.fa
```
```
A callset of non-reference L1Hs in HG002, HG003, and HG004 [a Personal Genome Project trio derived from the Genome in a Bottle (GIAB) Consortium] using PALMER is available under:
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/PacBio_PALMER_11242017/
```

## Output 
We have two outputs: 'output_calls.txt' & 'output_TSD_reads.txt'.

'output_calls.txt' is the summary for all non-ref MEI calls.

'output_TSD_reads.txt' contains all details you want for the high confident (HC) supporting reads (SRs).

* By using raw sub-reads from a ~50x coverage PacBio genome, we recommend a cutoff for HC calls as ≥1 HC-SR and ≥5 SRs.

## Logs

**Ver1.4** Feb.27th.2019

* Highly improved calling for SVA.
* Now PALMER supports other reference based bam files besides GRCh37, GRCh38 and hg19.
* Time consuming: to run PALMER on chr1/GRCh37, calling would cost ~24 hours (LINE-1/GRCh37), ~28 hours (Alu) or ~4 hours (SVA), for 8gb running memory minimun.
* A fatal bug fixed.
* Optimized scripts and outputs.
* Minor bugs fixed.

**Ver1.3.3** Feb.3rd.2019 ^^^(*￣(oo)￣)^^^ Happy Lunar New Year! Year of the Pig!! ^^^(*￣(oo)￣)^^^ 

* A steady and sensitive version for detection all MEIs (LINE-1, Alu and SVA) in human genome.
* Time consuming: to run PALMER on chr1, calling would cost ~150 hours (LINE-1), ~20 hours (Alu) or ~4 hours (SVA), for 8gb running memory minimun. Right now, PALMER does not support multi-thread processing.
* Now PALMER can output whole structure of MEI sequence, including inserted main sequence as well as different characteristics (TSD, TD, polyA tail) that have been supported by previous version already.
* A fatal bug related to PacBio read name from fastq data fixed. 
* Minor bugs fixed.

**Ver1.3.0** Dec.22th.2018

* Frozen version for L1 false negative paper.
* Output 26mer sequence at 5' junction
* Minor bugs fixed

**Ver1.3** Dec.5th.2018

* Highly optimized performance of LINE-1 calling using raw sub-reads
	> Imported '5' inverted sequence detection' module (two priming mechanism induced)
	
	> Optimized 'CNV-related false positive exclusion' module by using raw sub-reads (deletion-, duplication-, insertion-, inversion-related false positives)
	
	> Optimized 'TSD finding' module
	
	> Optimized speed of calling MEIs (I/O related)
* Optimized output files
	> Add '5' inverted sequence' output
	
	> Add 'Length of poly-A tail' output
	
	> Add 'Number of high confident supporting reads' output
* Minor bugs fixed

**Ver1.2** Sep.5th.2018

* Better performance for Alu calling
* Import 'CNV-related false positive exclusion' module
* Import genotyping module (not online yet)
* Import 'Customized sequence finding and genome masking' moudule
* Several minor bugs fixed
* Optimized output files
* Optimized codes and annotations
* LICENSE added

**Ver1.1** Apr.24th.2018

* Alu and SAV detection module online.

**Ver1.0** Feb.14th.2018
