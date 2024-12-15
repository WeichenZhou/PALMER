////copyright by ArthurZhou @ UMich&Fudan&HUST
#include <stdlib.h>
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <cstdlib>

using namespace std;

int samtools(string working_dir, string input_bam, string chr, string start, string end, string fasta){
    
    //cout<<"Samtools Step is now running."<<endl;
    //##########hard code warming########## -q -F ##########
    //std::ios::sync_with_stdio(false);
    //std::cin.tie(0);
    
    string sys;
    sys="samtools view -q 10 -F 0x100 -F 0x200  -F 0x400 -T "+fasta+" "+input_bam+" "+chr+":"+start+"-"+end+" |sed -e 's/[ ][ ]*/_/g'  > "+working_dir+"region.sam";
    
    char *syst = new char[sys.length()+1];
    strcpy(syst, sys.c_str());
    system(syst);
    
    //string sys_replace;
    //sys_replace="sed -e 's/[ ][ ]*/_/g' "+working_dir+"region.pre.sam"+" > "+working_dir+"region.sam";
    
    ///char *syst_replace = new char[sys_replace.length()+1];
    //strcpy(syst_replace, sys_replace.c_str());
    //system(syst_replace);
    
    return 0;
}
