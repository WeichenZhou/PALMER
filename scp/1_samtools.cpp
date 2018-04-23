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

int samtools(string working_dir1, string working_dir, string input_bam, string chr, string start, string end){
    
    //cout<<"Samtools Step is now running."<<endl;
    
    string sys;
    sys="samtools view -q 10 -F 0x100 -F 0x200 -F 0x800 -F 0x400 "+input_bam+" "+chr+":"+start+"-"+end+" > "+working_dir+"region.sam";
    
    char *syst = new char[sys.length()+1];
    strcpy(syst, sys.c_str());
    system(syst);
    return 0;
}
