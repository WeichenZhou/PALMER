# PALMER

Pre-mAsking Long reads for Mobile Element inseRtion

* PALMER detects non-reference MEI events (LINE, Alu, SVA, and HERVK) and other insertions by using the indexed reference-aligned BAM/CRAM files from long-read technology as inputs. It masks the aligned portions of reads, defines the significant characteristics of MEIs (TSD motifs, 5' inverted sequence, 3' transduction sequence, polyA-tail), and reports sequences for each insertion event.
* The ideal structure of an MEI event would be 5’-TSD-(5'inverted)-MEI-polyA-(TransD-polyA)-TSD-3’.
* PALMER is able to detect other categories (e.g. numts) of non-reference insertion sequence under the customized setup by the user.

Required resources:
```
 samtools/1.3.1  https://github.com/samtools/samtools
 ncbi-blast++/2.10.0  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ (Lower version will introduce fatal bugs.)
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
USAGE:

Required

--input
         aligned long-read sequencing BAM file with directory path

--workdir
         the user's working directory. Please follow the format /your/woking/directory/ !!don't forget the last '/'!!

--ref_ver (options: hg19, GRCh37, GRCh38 or other)
         reference genome used for the aligned file ('other' option for the cusmized genome out of hg19, GRCh37 or GRCh38)

--ref_fa
         indexed fasta file of reference genome fasta file with directory path used for the aligned bam/cram file (wrong reference will cause error information)

--type (options: LINE, ALU, SVA, HERVK, or CUSTOMIZED (if you want to setup your costomized sequence))
         type of MEIs or other kinds of insertions to detect

--mode (options: raw, or asm)      
         type of input sequencing to be processed (raw: raw nanopore/PacBio-sub reads; asm: assembled contigs)

--chr (default: ALL (for whole genome, not recommended); options: chromosome1, chromosome2, ...chromosomeY)
         chromosome name for PALMER to run. !!The chromosome names should be consistent with the ones in reference genome version!! e.g. for GRCh37, to run PALMER on chromosome1, the option should be '1', while for GRCh38 it should be 'chr1'

Optional

--start (default: Null)
         start position in the genome for PALMER to run (default is null). !!It should go with --end if assigned

--end (default: Null)      
         end position in the genome for PALMER to run (default is null). !!It should go with --start if assigned
            
--custom_seq (default: Null)
         .fasta file with directory path to customize your insertion finding. e.g. NUMTs, MEIs in other species.

--TSD_finding (Fixed: TRUE for all MEIs ,or default: FALSE for CUSTOMIZED insertion)
         whether to run TSD motif finding module for your insertion calling

--len_custom_seq (MUST set up when activating TSD_finding for CUSTOMIZED insertion, otherwise CLOSED)
         interger value for the length of your customized sequence WITHOUT polyA tact

--L_len (default: 25bp)
         the minimum length of putative LINE-1 aligned to L1.3 sequences

--output (default: output)
         the prefix of the output file
```

Examples
```
1) Running PALMER on example PacBio subreads bam file under the 'example' folder to call LINE-1 insertions on GRCh38 genome
./PALMER --input $PALMER_Path/example/sample.bam --workdir $DirPath/ --ref_ver GRCh38 --output sample --type LINE --mode raw --chr chr19 --ref_fa $your.reference.file.path/GRCh38.fa

Results (sample_calls.txt & sample_TSD_reads.txt)  from example bam file can also be found under the 'example' folder.
```
```
2) Running PALMER on your aligned sequences on GRCh37 reference genome to call LINE-1 insertions in chromosome3 at position from 200,000 to 400,000
./PALMER --input $DirPath/your.bam.file --workdir $DirPath/ --ref_ver GRCh37 --output your.output.prefix --type LINE --mode raw --chr 3 --start 200000 --end 400000 --ref_fa $your.reference.file.path/hs37d5.fa
```
```
3) Running PALMER on your aligned assembled contigs in cram based on GRCh38 reference genome to call SVA insertions in chromosome3
./PALMER --input $DirPath/your.cram.file --workdir $DirPath/ --ref_ver GRCh38 --output your.output.prefix --type SVA --mode asm --chr chr3 --ref_fa $your.reference.file.path/GRCh38.fa
```
```
4) Running PALMER on your aligned bam to call Alu insertions in chromosome2a of Champanzee genome
./PALMER --input $DirPath/your.bam.file --workdir $DirPath/ --ref_ver other --output your.output.prefix --type ALU --mode raw --chr chr2a(chr.name.based.on.your.reference.fa) --ref_fa $your.reference.file.path/your.reference.fa 
```
```
5) Running PALMER on your aligned bam to call NumtS in chromosome5 of Champanzee genome
./PALMER --input $DirPath/your.bam.file --workdir $DirPath/ --ref_ver other --output your.output.prefix --chr chr5 --mode raw --ref_fa $your.reference.file.path/your.reference.fa --type CUSTOMIZED --custom_seq $your.custom_seq.file.path/Clint.mt 
```
```
6) Running PALMER on your aligned bam to call LINE-1 insertions in chromosomeX of mice genome
./PALMER --input $DirPath/your.bam.file --workdir $DirPath/ --output your.output.prefix --chr chrX --ref_ver other --mode raw --ref_fa $your.reference.file.path/your.reference.fa --type CUSTOMIZED --custom_seq $your.custom_seq.file.path/L1MdA_consensus.fa --TSD_finding TRUE --len_custom_seq (int)
```
```
7)
A callset of non-reference L1Hs in HG002, HG003, and HG004 [a Personal Genome Project trio derived from the Genome in a Bottle (GIAB) Consortium] using PALMER is available under:
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/PacBio_PALMER_11242017/
```

## Output 
We have two outputs: 'output_calls.txt' & 'output_TSD_reads.txt'.

'output_calls.txt' is the summary for all non-ref MEI calls.

'output_TSD_reads.txt' contains all details you want for the high confident (HC) supporting reads (SRs).

* By using raw sub-reads from a ~50x coverage PacBio genome, we recommend a cutoff for HC calls as ≥1 HC-SR and ≥5 (10% of the average coverage) SRs.


## Citation

For general use or LINE-1s:
* Weichen Zhou, Sarah B Emery, Diane A Flasch, Yifan Wang, Kenneth Y Kwan, Jeffrey M Kidd, John V Moran, Ryan E Mills,
[Identification and characterization of occult human-specific LINE-1 insertions using long-read sequencing technology](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz1173/5680708), 
Nucleic Acids Research, 2019, gkz1173, `https://doi.org/10.1093/nar/gkz1173`

For all MEIs:
* Torrin L. McDonald*,  Weichen Zhou*,  Christopher Castro,  Camille Mumm,  Jessica A. Switzenberg,  Ryan E. Mills,  Alan P. Boyle,
[Cas9 targeted enrichment of mobile elements using nanopore sequencing](https://www.nature.com/articles/s41467-021-23918-y), 
Nature Communications, 2021, `https://doi.org/10.1038/s41467-021-23918-y`

For PALMER2.0:
In preparation!

## Contact

* arthurz@med.umich.edu

## Logs
**Ver2.0.1** Nov.10th.2024! PALMER2.0.1 is online now!

* A frozen version for SMaHT and HGSVC3.

**Ver2.0.0** May.20th.2022! PALMER2.0.0 is online now!! 520 (｡・ω・｡)ﾉ♥♥♥♥♥♥♥♥♥♥!! 

* Capability of calling insertions in assembled contigs!!
* A couple of major bugs fixed!!
* Improved running time!!
* Minor bugs fixed.

**Ver1.7.2** Nov.28th.2020! Happy Thanksgiving!!

* Improved HIFI reads calling!!
* A couple of major bugs fixed!!
* Improved running time!!

**Ver1.7** Nov.11th.2020! Happy Singles Day & happy shopping!!

* Enabled HERV-K calling!!
* Enabled specific region calling!!
* Enabled cram file calling!!
* Minor bugs fixed.

**Ver1.6.2.Enhanced** Sep.27th.2020 by Jixing Guan

* Optimized PALMER and make samtools as build-in lib

**Ver1.6.2** May.19th.2020

* Fixed a bug that would crash the software when the read names are not unique in the raw fastq regarding the PacBio subreads.

**Ver1.6.1** May.19th.2020

* Fixed a bug when calling customized insertion sequences without TSD finding module.

**Ver1.6** May.11th.2020

* Frozen version for "Refining polymorphic retrotransposon insertions in human genomes". Good Luck!

**Ver1.5.1** May.7th.2020

* Highly optimized the performance of calling customized insertion sequences (non-humman genomes, non-MEIs)!!
* Sample bam file added, example updated, results from sample bam updated!!
* A fatal bug fixed calling L1NE-1 since Ver1.4.1.
* A fatal bug fixed related to the environment of computing clusters and the version of BLASTn. Now require the version ncbi-blast++/2.10.0.
* Minor bugs fixed.

**Ver1.5** May.4th.2020 "MAY THE FORCE BE WITH YOU!"

* Added one more option for the length of the customized insertion sequence.
* Optimized the performance for customized insertion sequence finding!
* Minor bugs fixed.

**Ver1.4.1** Nov.14th.2019

* Added one more option for an adjustable length of putative LINE-1 aligned to L1.3 sequences.

**Ver1.4** Feb.27th.2019

* Highly improved calling for SVA.
* Now PALMER supports other reference-based bam files besides GRCh37, GRCh38, and hg19.
* Time consumption: to run PALMER on chr1/GRCh37, calling would cost ~24 hours (LINE-1/GRCh37), ~28 hours (Alu) or ~4 hours (SVA), for 8gb running memory minimum.
* A fatal bug fixed.
* Optimized scripts and outputs.
* Minor bugs fixed.

**Ver1.3.3** Feb.3rd.2019 ^^^(*￣(oo)￣)^^^ Happy Lunar New Year! Year of the Pig!! ^^^(*￣(oo)￣)^^^ 

* A steady and sensitive version for detection all MEIs (LINE-1, Alu, and SVA) in human genome.
* Time consumption: to run PALMER on chr1, calling would cost ~150 hours (LINE-1), ~20 hours (Alu) or ~4 hours (SVA), for 8gb running memory minimum. Right now, PALMER does not support multi-thread processing.
* Now PALMER can output whole structure of MEI sequence, including inserted main sequence as well as different characteristics (TSD, TD, polyA tail) that have been supported by the previous version already.
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
* Import 'Customized sequence finding and genome masking' module
* Several minor bugs fixed
* Optimized output files
* Optimized codes and annotations
* LICENSE added

**Ver1.1** Apr.24th.2018

* Alu and SAV detection module online.

**Ver1.0** Feb.14th.2018
