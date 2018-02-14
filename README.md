# PALMER

Pre-mAsking Long reads for Mobile Element inseRtion

PALMER is used to detect non-reference MEI events within the masked sequence data. It uses the reference-aligned BAM files from long-read technology as inputs. 

* It delineates large genome data into small bins (100kb bins), which will be piped into multi-threads and processed separately, to increase the efficiency of program. 
* The known repeats (L1NE-1, Alus or SVAs in reference) are used to mask the portions of reads that aligned to these repeats. 
* After obtaining the pre-masked reads, PALMER searches against the insertion sequence library (for L1Hs sequence, GenBank Accession: L19088)  by using Blastn. Then it uses the cigar information in the reads and joint the truncated segments from Blastn into one in a single read. These reads with insertion sequence are considered as supporting reads. 
* PALMER identifies the candidate TSD motif in 50bp 5’ upstream and 3kb 3’ downstream of insertion sequence for each read. 
* It then runs a module for filtering the candidate TSD motif and identifying transduction/polyA sequence. PALMER will define the supporting read with/without valid TSD motif and with/without valid transduction and polyA sequence. The ideal structure of an event having transduction sequence should be 5’-TSD-L1Hs-polyA-TransD-polyA-TSD-3’. 
* Afterwards, PALMER will cluster all supporting reads in one loci (or supporting one event) of the genome, as well as cluster the TSD motif among all supporting reads for one event and choose the most confident TSD motif with/without transduction/polyA sequence. PALMER will obtain two numbers for one event, the number of supporting reads and the number of supporting reads with predicted TSD motif. 
* Finally, PALMER will combine all events in each bin and output all candidate non-reference MEIs.


Required Resources:
```
  samtools/1.3.1  https://github.com/samtools/samtools
  ncbi-blast++/2.4.0  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

## Getting Started

Download and Install
```
git clone git@github.com:mills-lab/PALMER.git
cd PALMER
make
```

Parameters
```
  --input, input aligned long-read sequencing file
  --workdir, the user's working directory
  --ref, reference genome of the aligned file (options: GRCh37 or GRCh38)
  --type, type of MEIs to detect (options: LINE, ALU or SVA)
  --chr, chr name for PALMER to run (default: whole genome; options: chr1, chr2, ...chrY)
  --output, name of output file (default: output.txt)
```
