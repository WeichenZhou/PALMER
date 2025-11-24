////copyright by ArthurZhou @ UMich&Fudan&HUST
#include "common.hpp"

int samtools(string working_dir, string input_bam, string chr, string start, string end, string fasta, string MAPQ){

    string sys;
    sys="samtools view -q "+MAPQ+" -F 0x100 -F 0x200  -F 0x400 -T "+fasta+" "+input_bam+" "+chr+":"+start+"-"+end+" |sed -e 's/[ ][ ]*/_/g'  > "+working_dir+"region.sam";
    
    char *syst = new char[sys.length()+1];
    strcpy(syst, sys.c_str());
    system(syst);
    
    return 0;
}
