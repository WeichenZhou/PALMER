////copyright by ArthurZhou @ UMich&FUDAN&HUST
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

int RefSelector(string ref_file, string working_dir, string chr, string start, string end, int ref_n, string seq_in){
    
    string sys;
    if(ref_n==37){
        sys="samtools faidx "+ref_file+" "+chr+":"+start+"-"+end+" >"+working_dir+seq_in+".ref.fa";
    }
    else if(ref_n==38){
        sys="samtools faidx "+ref_file+" chr"+chr+":"+start+"-"+end+" >"+working_dir+seq_in+".ref.fa";
    }
    
    char *syst= new char[sys.length()+1];
    strcpy(syst, sys.c_str());
    system(syst);
    return 0;
}
