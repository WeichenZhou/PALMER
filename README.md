# PALMER

Pre-mAsking Long reads for Mobile Element inseRtion

PALMER is used to detect non-reference MEI events within the masked sequence data. It uses the indexed reference-aligned BAM files from long-read technology as inputs. 

* It delineates large genome data into small bins (100kb bins), which will be processed separately, to increase the efficiency of program. 
* The track from Repeatmasker (https://www.girinst.org/) is used to mask the portions of reads that aligned to these repeats. Or you can customize your own repeats track to mask the genome.
* After obtaining the pre-masked reads, PALMER searches against the insertion sequence library (e.g. L1Hs sequence, GenBank Accession: L19088) by using Blastn. Then it uses the cigar information in the reads and joint the truncated segments from Blastn into one in a single read. These reads with insertion sequence are considered as supporting reads. 
* PALMER identifies the candidate TSD motif in 50bp 5’ upstream and 3kb 3’ downstream of insertion sequence for each read. 
* It then runs a module for filtering the candidate TSD motif and identifying transduction/polyA sequence. PALMER will define the supporting read with/without valid TSD motif and with/without valid transduction and polyA sequence. The ideal structure of an event having transduction sequence should be 5’-TSD-L1Hs-polyA-TransD-polyA-TSD-3’. 
* Afterwards, PALMER will cluster all supporting reads in one loci (or supporting one event) of the genome, as well as cluster the TSD motif among all supporting reads for one event and choose the most confident TSD motif with/without transduction/polyA sequence. PALMER will obtain two numbers for one event, the number of supporting reads and the number of supporting reads with predicted TSD motif. 
* Finally, PALMER will combine all events in each bin and output all candidate non-reference MEIs.


Required resources:
```
  samtools/1.3.1  https://github.com/samtools/samtools
  ncbi-blast++/2.4.0  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

## Getting started

Download and Install
```
git clone git@github.com:mills-lab/PALMER.git
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

--ref_ver (options: hg19, GRCh37 or GRCh38)
         reference genome used for the aligned file (only human genome rightnow)

--ref_fa
         indexed fasta file of reference genome fasta file with directory path used for the aligned file

--type (options: LINE, ALU or SVA)
         type of MEIs to detect

--chr (default: whole genome; options: chr1, chr2, ...chrY)
         chr name for PALMER to run (if running for whole genome, don't need to assign)

--custom_seq (default:no input)
         fasta file with directory path to customize your insertion finding

--custom_index (default:no input, if you have --custom_seq parameter without --custom_index, PALMER will work without masking step)
         index file with directory path to mask the genome for your insertion finding

--output (default: output)
         prefix of output file
```

Example
```
./PALMER --input $DirPath/NA12878.washu.alignment_hs37d5.1.bam --workdir $DirPath/chr1.line.0406/ --ref_ver GRCh37 --output test --type LINE --chr chr3 --ref_fa $DirPath/hs37d5.fa
```

## Output format 
We have two outputs: 

output.calls.txt & output.TSD_reads.txt.

## Logs
**Ver1.2** Sep.5th.2018

* Better performance for Alu calling
* Import 'CNV-related false positive exclusion' module
* Import genotyping module (not online yet)
* Import 'customized sequence finding and genome masking' moudle
* Several minor bugs fixed
* Optimized output files
* Optimized codes and annotations

**Ver1.1** Apr.24th.2018
* Alu and SAV detection module online.

**Ver1.0** Feb.14th.2018
