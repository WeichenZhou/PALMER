//copyright by ArthurZhou @ UMich&Fudan&HUST
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

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "1_samtools.cpp"
//#include "2_rm_selector.cpp"
#include "3_read_masker.cpp"
#include "4_blastn.cpp"
#include "5_blastn_caller.cpp"
#include "6_TSD_seq.cpp"
#include "7_FP_ex.cpp"
#include "8_calling.cpp"
using namespace std;

int tube(string working_dir, string input_bam, string chr, int start, int end, string type, int ref_n, string direc, string ref_fa, int tsd, int L_len, int cus_seq_len, string mode){
    
    //std::ios::sync_with_stdio(false);
    //std::cin.tie(0);
    
//building working directory
    int start1, end1;
    start1= start;
    end1 = end;
    stringstream ss1, ss2;
    ss1 << start1;
    string s_start =ss1.str();
    ss2 << end1;
    string s_end =ss2.str();
    
    string WD_tube=working_dir+chr+"_"+s_start+"_"+s_end+"/";
    
    string sys_mkdir="mkdir "+WD_tube;
    
    cout<<"Working in the direcotry "<<sys_mkdir<<"."<<endl;
    
    char *syst_mkdir =new char[sys_mkdir.length()+1];
    strcpy(syst_mkdir, sys_mkdir.c_str());
    system(syst_mkdir);
 
//1. samtools view
    
    samtools(WD_tube, input_bam, chr, s_start, s_end ,ref_fa);
    cout<<"Samtools Step for region "+chr+"_"+s_start+"_"+s_end+" now completed."<<endl;
    
//Repeat region output
    /*
    ofstream file1;
    
    string sys_RMloc = WD_tube+"RM.loc";
    char *syst_RMloc = new char[sys_RMloc.length()+1];
    strcpy(syst_RMloc, sys_RMloc.c_str());
    /*
    string chr_fix;
    chr_fix=chr;
    if(ref_n==37){
        chr_fix="chr"+chr;
    }
    
    file1.open(syst_RMloc);
    file1<<chr<<'\t'<<(start-100000)<<'\t'<<(end+100000)<<endl;
    */

//2. Repeat region selector
    
    //RMSelector(working_dir, WD_tube, sys_region);
    //cout<<"2. RMSelector Step for region "+chr+"_"+s_start+"_"+s_end+" is now completed."<<endl;
    
//3. Read masker
    
    ReadMasker(WD_tube, mode);
    //cout<<"3. Read Masker Step for region "+chr+"_"+s_start+"_"+s_end+" is now completed."<<endl;
    cout<<"Pre-masking step for "+chr+"_"+s_start+"_"+s_end+" completed."<<endl;
 
//4. Blastn
    
    blastn(WD_tube, type, direc);
    //cout<<"4. Blastn Step for region "+chr+"_"+s_start+"_"+s_end+" is now completed."<<endl;
    cout<<"Blastn Step for region "+chr+"_"+s_start+"_"+s_end+" completed."<<endl;
    
//5. Blastn caller
    
    BlastnCaller(WD_tube, chr, type, L_len, cus_seq_len);
    //cout<<"5. Blastn Caller Step for region "+chr+"_"+s_start+"_"+s_end+" is now completed."<<endl;
 
    cout<<"Single read calling step for "+chr+"_"+s_start+"_"+s_end+" completed."<<endl;

//6. TSD module
    
    tsd_module(WD_tube, type, tsd);
    if(tsd==1){
        cout<<"TSD_module step for "+chr+"_"+s_start+"_"+s_end+" completed."<<endl;
    }
    else if(tsd==0){
        cout<<"Skipped TSD_module step for "+chr+"_"+s_start+"_"+s_end+"."<<endl;
    }
    //getchar();

//7. Flase postive exclusion module
    
    fp_ex(WD_tube, ref_fa, chr, type, tsd);
    if(tsd==1){
        cout<<"False positive exclusion step for "+chr+"_"+s_start+"_"+s_end+" completed."<<endl;
    }
    else if(tsd==0){
        cout<<"Skipped false positive exclusion step for "+chr+"_"+s_start+"_"+s_end+"."<<endl;
    }
    
    //getchar();
    
//8. Calling module
    
    calling(WD_tube, type, tsd);
    cout<<"Calling step for "+chr+"_"+s_start+"_"+s_end+" completed."<<endl;
    
    return 0;
}
