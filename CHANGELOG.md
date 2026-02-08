# Changelog

**Ver2.3.2** Feb.7nd.2026! PALMER2.3.2

* Fixed a rare bug for reading the same cigar from different reads/contigs.
* Fix the --help information.

**Ver2.3.1** Jan.23rd.2026! PALMER2.3.1

* Fixed the bug for reading the CRAM files.
* Fixed the bug for TSD_finding=FALSE.

**Ver2.3** Dec.6th.2025! PALMER2.3 念头通达

* Reorganized output files and added columns in the output files for potential supporting reads from 5' end, 3' end, and the go-through.
* Optimized ALU calling in terms of running time and identity accuracy. 
* Optimized BLASTn calling.
* Implemented samtools API. 
* Implemented MSA for consensus INS_SEQ with hc supporting read weight=3.0
* Fixed a bug that overestimates the number of supporting reads, particularly in asm mode.
* Decreased the number of intermediate files.
* Running time improved.
* Minor format bugs fixed.
* Example updated.
* Scripts cleaned.
* README updated.

**Ver2.2** Nov.25th.2025! PALMER2.2 Happy Thanksgiving!

* Add a module for multi-threads. 
* Add a module to output to VCF files.
* Disable the function of blastn for reporting usage statistics to NCBI.
* Add a module for removing (or keeping) intermediate files to avoid crashing the file system.
* Add a genotyping module based on a Generalized Gaussian Mixture model. But still under development.
* Replace htslib with samtools
* Minor format bugs fixed.
* Example updated.
* Scripts cleaned.
* README updated.

**Ver2.1.1** Aug.20th.2025! PALMER2.1 

* Frozen version for SMaHT "Benchmarking and integrative detection of low frequency somatic mobile element insertions in human tissue". Good luck!
* A major bug in calling when the sequencing depth is very deep. It stole the supporting reads.
* Add an option to set MAPQ
* Minor bugs fixed

**Ver2.0.1** Nov.10th.2024! PALMER2.0.1 is online now!

* A frozen version for HGSVC3.

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

* Optimized PALMER and made samtools a built-in lib

**Ver1.6.2** May.19th.2020

* Fixed a bug that would crash the software when the read names are not unique in the raw fastq regarding the PacBio subreads.

**Ver1.6.1** May.19th.2020

* Fixed a bug when calling customized insertion sequences without the TSD finding module.

**Ver1.6** May.11th.2020

* Frozen version for "Refining polymorphic retrotransposon insertions in human genomes". Good Luck!

**Ver1.5.1** May.7th.2020

* Highly optimized the performance of calling customized insertion sequences (non-human genomes, non-MEIs)!!
* Sample bam file added, example updated, results from sample bam updated!!
* A fatal bug was fixed calling L1NE-1 since Ver1.4.1.
* A fatal bug was fixed related to the environment of computing clusters and the version of BLASTn. Now require the version ncbi-blast++/2.10.0.
* Minor bugs fixed.

**Ver1.5** May.4th.2020 "MAY THE FORCE BE WITH YOU!"

* Added one more option for the length of the customized insertion sequence.
* Optimized the performance for customized insertion sequence finding!
* Minor bugs fixed.

**Ver1.4.1** Nov.14th.2019

* Added one more option for an adjustable length of putative LINE-1 aligned to L1.3 sequences.

**Ver1.4** Feb.27th.2019

* Highly improved calling for SVA.
* Now PALMER supports other reference-based BAM files besides GRCh37, GRCh38, and hg19.
* Time consumption: to run PALMER on chr1/GRCh37, calling would cost ~24 hours (LINE-1/GRCh37), ~28 hours (Alu), or ~4 hours (SVA), for 8 GB running memory minimum.
* A fatal bug has been fixed.
* Optimized scripts and outputs.
* Minor bugs fixed.

**Ver1.3.3** Feb.3rd.2019 ^^^(*￣(oo)￣)^^^ Happy Lunar New Year! Year of the Pig!! ^^^(*￣(oo)￣)^^^ 

* A steady and sensitive version for detecting all MEIs (LINE-1, Alu, and SVA) in the human genome.
* Time consumption: to run PALMER on chr1, calling would cost ~150 hours (LINE-1), ~20 hours (Alu), or ~4 hours (SVA), for 8 GB running memory minimum. Right now, PALMER does not support multi-thread processing.
* Now PALMER can output the whole structure of the MEI sequence, including the inserted main sequence as well as different characteristics (TSD, TD, polyA tail) that have been supported by the previous version already.
* A fatal bug related to PacBio read name from fastq data has been fixed. 
* Minor bugs fixed.

**Ver1.3.0** Dec.22th.2018

* Frozen version for L1 false negative paper.
* Output 26mer sequence at 5' junction
* Minor bugs fixed

**Ver1.3** Dec.5th.2018

* Highly optimized performance of LINE-1 calling using raw sub-reads
	> Imported '5' inverted sequence detection' module (two-priming mechanism induced)
	
	> Optimized 'CNV-related false positive exclusion' module by using raw sub-reads (deletion-, duplication-, insertion-, inversion-related false positives)
	
	> Optimized 'TSD finding' module
	
	> Optimized speed of calling MEIs (I/O related)
* Optimized output files
	> Add '5' inverted sequence' output
	
	> Add 'Length of poly-A tail' output
	
	> Add 'Number of high confidence supporting reads' output
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


