## This is a stable fork from https://github.com/WeichenZhou/PALMER. If you have any issues, please visit that repository and subsequent updates and corrections will be pulled here.


# PALMER

Pre-mAsking Long reads for Mobile Element inseRtion

* PALMER detects non-reference MEI events (LINE, Alu and SVA) and other insertions, by using the indexed reference-aligned BAM files from long-read technology as inputs. It uses the track from [Repeatmasker](https://www.girinst.org/) to mask the portions of reads that aligned to these repeats, defines the significant characteristics of MEIs (TSD motifs, 5' inverted sequence, 3' transduction sequence, polyA-tail), and reports sequences for each insertion event.
* The ideal structure of an MEI event should be 5’-TSD-(5'inverted)-MEI-polyA-(TransD-polyA)-TSD-3’. 


Required resources:
```
 samtools/1.3.1  https://github.com/samtools/samtools
 ncbi-blast++/2.4.0  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
 git-lfs
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
         the user's working directory. Please follow the format /your/woking/directory/ !!don't forget the last '/'!!

--ref_ver (options: hg19, GRCh37, GRCh38 or other)
         reference genome used for the aligned file ('other' option for the cusmized genome out of hg19, GRCh37 or GRCh38)

--ref_fa
         indexed fasta file of reference genome fasta file with directory path used for the aligned bam file (wrong reference will cause error information)

--type (options: LINE, ALU, SVA, or CUSTOMIZED (if you want to setup your costomized sequence))
         type of MEIs or other kind of insertions to detect

--chr (default: ALL (for whole genome, not recommended); options: chromosome1, chromosome2, ...chromosomeY)
         chromosome name for PALMER to run. !!The chromosome names should be consistent with the ones in reference genome version!! e.g. for GRCh37, to run PALMER on chromosome1, the option should be '1', while for GRCh38 it should be 'chr1'

--custom_seq (default:no input)
         .fasta file with directory path to customize your insertion finding. e.g. NUMTs, MEIs in other species.

--custom_index (default:no input; if you have both '--ref_ver other' and '--type LINE/ALU/SVA', you must give PALMER a index file (format: "CHR'	'START'	'END'	'MEI_NAME'
'" for each MEI to be masked in each line) for masking module; if you have --custom_seq parameter without --custom_index, PALMER will work without the masking step)
         index file with directory path to mask the genome for your insertion finding

--TSD_finding (Fixed:TRUE for all MEIs ,or default: FALSE for CUSTOMIZED insertion)
         whether to run TSD motif finding module for your insertion calling

--len_custom_seq (MUST set up when activate TSD_finding for CUSTOMIZED insertion, otherwise CLOSED)
         integer value for the length of your customized sequence WITHOUT polyA tact

--L_len (default: 25bp)
         the minimum length of putative LINE-1 aligned to L1.3 sequences

--output (default: output)
         prefix of output file
```

Examples
```
1) Running PALMER on your aligned bam based on GRCh37 reference genome to call LINE-1 insertions in chromosome3
./PALMER --input $DirPath/your.bam.file --workdir $DirPath/ --ref_ver GRCh37 --output your.output.prefix --type LINE --chr 3 --ref_fa $your.reference.file.path/hs37d5.fa
```
```
2) Running PALMER on your aligned bam based on GRCh38 reference genome to call SVA insertions in chromosome3
./PALMER --input $DirPath/your.bam.file --workdir $DirPath/ --ref_ver GRCh38 --output your.output.prefix --type SVA --chr chr3 --ref_fa $your.reference.file.path/GRCh38.fa
```
```
3) Running PALMER on your aligned bam to call Alu insertions in chromosome2a of Champanzee genome
./PALMER --input $DirPath/your.bam.file --workdir $DirPath/ --ref_ver other --output your.output.prefix --type ALU --chr chr2a(chr.name.based.on.your.reference.fa) --ref_fa $your.reference.file.path/your.reference.fa --custom_index reference.Alu.coordinates.in.your.reference.Chimpanzee.genome 
```
```
4) Running PALMER on your aligned bam to call NumtS in chromosome5 of Champanzee genome
./PALMER --input $DirPath/your.bam.file --workdir $DirPath/ --ref_ver other --output your.output.prefix --chr chr5 --ref_fa $your.reference.file.path/your.reference.fa --type CUSTOMIZED --custom_seq $your.custom_seq.file.path/Clint.mt --custom_index $your.custom_index.file.path/Chimp_ref_NumtS.bed
```
```
5) Running PALMER on your aligned bam to call LINE-1 insertions in chromosomeX of mice genome
./PALMER --input $DirPath/your.bam.file --workdir $DirPath/ --output your.output.prefix --chr chrX --ref_ver other --ref_fa $your.reference.file.path/your.reference.fa --type CUSTOMIZED --custom_seq $your.custom_seq.file.path/L1MdA_consensus.fa --custom_index $your.custom_index.file.path/mm10_ucsc_repeatmasker_LINE.bed --TSD_finding TRUE --len_custom_seq (int)
```
```
6)
A callset of non-reference L1Hs in HG002, HG003, and HG004 [a Personal Genome Project trio derived from the Genome in a Bottle (GIAB) Consortium] using PALMER is available under:
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/PacBio_PALMER_11242017/
```

## Output 
We have two outputs: 'output_calls.txt' & 'output_TSD_reads.txt'.

'output_calls.txt' is the summary for all non-ref MEI calls.

'output_TSD_reads.txt' contains all details you want for the high confident (HC) supporting reads (SRs).

* By using raw sub-reads from a ~50x coverage PacBio genome, we recommend a cutoff for HC calls as ≥1 HC-SR and ≥5 SRs.


## Citation

* Weichen Zhou, Sarah B Emery, Diane A Flasch, Yifan Wang, Kenneth Y Kwan, Jeffrey M Kidd, John V Moran, Ryan E Mills,
[Identification and characterization of occult human-specific LINE-1 insertions using long-read sequencing technology](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz1173/5680708), 
Nucleic Acids Research, 2019, gkz1173, `https://doi.org/10.1093/nar/gkz1173`


## Logs

**Ver1.5** May.4th.2020 "MAY THE FORCE BE WITH YOU!"

* Added one more option for the length of customized insertion sequence.
* Optimized the performance for customized insertion sequence finding!
* Minor bugs fixed.

**Ver1.4.1** Nov.14th.2019

* Added one more option for an adjustable length of putative LINE-1 aligned to L1.3 sequences.

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
