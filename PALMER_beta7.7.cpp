 //copyright by ArthurZhou @ UMich&FUDAN&HUST
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
#include <numeric>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "scp/tube.cpp"

//--#include "index/LINEs.regions.GRCh37.h"
//#include "index/Alu.regions.GRCh37.h"
//#include "index/SVA.regions.GRCh37.h"
//---#include "index/region.split.index.GRCh37.h"
//#include "index/LINEs.regions.GRCh38.h"
//#include "index/Alu.regions.GRCh38.h"
//#include "index/SVA.regions.GRCh38.h"
//#include "index/region.split.index.GRCh38.h"
//--#include "RM/L1.3.fasta.h"
//--#include "RM/AluY.fasta.h"
//--#include "RM/SVA_F.fasta.h"

using namespace std;

int main(int argc, char *argv[]){

//parameters_start
    
    string T, WD, inputF, output, SP, ref, CHR;
    int NUM_threads=30;
    ifstream file1;
    //ofstream file2;
    int flag_wd=0;
    int flag_inputf=0;
    int flag_T=0;
    output="output.txt";
    SP="Human";
    //T="LINE";
    CHR="ALL";
    int ref_n=0;
    
    string dir;
    dir=argv[0];
    
    for(int i=1;i!=argc;i++){
        if(strncmp(argv[i],"--workdir",9)==0){
            WD=argv[i+1];
            flag_wd=1;
        }
        if(strncmp(argv[i],"--chr",5)==0){
            CHR=argv[i+1];
        }
        if(strncmp(argv[i],"--ref",5)==0){
            ref=argv[i+1];
            if(ref=="GRCh37") ref_n=37;
            else if(ref=="GRCh38") ref_n=38;
            else {
                cout<<"PLEASE INPUT A CORRECT REFERENCE :("<<endl;
            }
        }
        if(strncmp(argv[i],"--type",6)==0){
            T=argv[i+1];
            flag_T=1;
        }
        if(strncmp(argv[i],"--input",7)==0){
            flag_inputf=1;
            file1.open(argv[i+1],ios::in | ios::binary);
            if (!file1.is_open())
            {
                cout <<"CANNOT OPEN INPUT FILE :("<< endl;
                continue;
            }
            
            inputF=argv[i+1];
        }
        if(strncmp(argv[i],"--output",7)==0){
            output=argv[i+1];
        }
        if(strncmp(argv[i],"--species",9)==0){
            SP=argv[i+1];
            cout <<"WE DO NOT NEED SPECIES PARAMETER RIGHT NOW :)"<< endl;
        }
        if(strncmp(argv[i],"--thread",8)==0){
            //NUM_threads=argv[i+1];
            cout <<"WE DO NOT SUPPORT MULTITHREADS RIGHT NOW :)"<< endl;
        }
        
    }
    
    
    if(flag_T==0){
        cout<<"PLEASE ASSIGN A MEI TYPE! LINE/ALU/SVA"<<endl;
        exit(1);
    }
    if(flag_wd==0){
        cout<<"PLEASE SET UP A WORKING DIRECOTY!"<<endl;
        exit(1);
    }
    if(flag_inputf==0){
        cout<<"PLEASE INPUT A FILE!"<<endl;
        exit(1);
    }
    if(ref_n==0){
        cout<<"PLEASE ASSIGN A CORRECT REFERENCE! GRCh37/GRCh38"<<endl;
        exit(1);
    }
    
    string parameter[6];
    parameter[0]=T;
    parameter[1]=WD;
    parameter[2]=inputF;
    parameter[3]=output;
    parameter[4]=CHR;
    parameter[5]=ref;
    
    cout<<"Variant type is "<<parameter[0]<<endl;
    cout<<"Working directory is "<<parameter[1]<<endl;
    cout<<"Input file is "<<parameter[2]<<endl;
    cout<<"Output file is "<<parameter[1]<<parameter[3]<<endl;
    cout<<"Running on "<<parameter[4]<<endl;
    cout<<"ref is "<<parameter[5]<<endl;
    
//Buildup & index or HG version chose
    string sys_dir="dirname "+dir;
    char *syst_dir=new char[sys_dir.length()+1];
    strcpy(syst_dir, sys_dir.c_str());
    
    vector<string> dir_conv;
    dir_conv.clear();
    FILE *pp =popen(syst_dir,"r");
    char tmp[1024];
    while (fgets(tmp, sizeof(tmp), pp) != NULL) {
        if (tmp[strlen(tmp) - 1] == '\n') {
            tmp[strlen(tmp) - 1] = '\0';
        }
        dir_conv.push_back(tmp);
    }
    pclose(pp);
    string direc;
    direc=accumulate(dir_conv.begin(),dir_conv.end(),direc);
    
    string buildup=direc+"/index/";
    /*
    string sys_buildup="mkdir "+buildup;
    char *syst_buildup=new char[sys_buildup.length()+1];
    strcpy(syst_buildup, sys_buildup.c_str());
    
    system(syst_buildup);
    */
    //string sys_line_region=buildup+"LINEs.regions";
    /*
    char *syst_line_region =new char[sys_line_region.length()+1];
    strcpy(syst_line_region, sys_line_region.c_str());
    
    ofstream file11;
    file11.open(syst_line_region);
    */
    
    string sys_region_index,sys_line_region;
    
    if(ref_n==37){
        sys_region_index=buildup+"region.split.index.GRCh37";
        if(T=="LINE"){
            sys_line_region=buildup+"LINEs.regions.GRCh37";
        }
        if(T=="ALU"){
            sys_line_region=buildup+"Alu.regions.GRCh37";
        }
        if(T=="SVA"){
            sys_line_region=buildup+"SVA.regions.GRCh37";
        }
    }
    else if(ref_n==38){
        sys_region_index=buildup+"region.split.index.GRCh38";
        if(T=="LINE"){
            sys_line_region=buildup+"LINEs.regions.GRCh38";
        }
        if(T=="ALU"){
            sys_line_region=buildup+"Alu.regions.GRCh38";
        }
        if(T=="SVA"){
            sys_line_region=buildup+"SVA.regions.GRCh38";
        }
    }
    
    cout<<sys_region_index<<" "<<sys_line_region<<endl;
 //original

    //string sys_region_index=buildup+"region.split.index";
    char *syst_region_index =new char[sys_region_index.length()+1];
    strcpy(syst_region_index, sys_region_index.c_str());
    //ofstream file12;
    //file12.open(syst_region_index);
    /*
    if(ref_n==37){
        file12<<region_split_index_GRCh37;
    }
    else if(ref_n==38){
    //    file12<<region_split_index_GRCh38;
    }*/
/* 
    string sys_region_index=buildup+"region.split.index.test";
    char *syst_region_index =new char[sys_region_index.length()+1];
    strcpy(syst_region_index, sys_region_index.c_str());
    ofstream file12;
    file12.open(syst_region_index);
    file12<<region_split_index_test;
 */
    /*
    string sys_l=buildup+"L1.3.fasta";
    char *syst_l =new char[sys_l.length()+1];
    strcpy(syst_l, sys_l.c_str());
    ofstream file13;
    file13.open(syst_l);
    file13<<L1_3_fasta;
    */
    ifstream file2;
    file2.open(syst_region_index);
    
    if (!file2.is_open())
    {
        cout <<"CANNOT OPEN INDEX FILE"<< endl;
        //cout<<"YEA"<<endl;
        exit(1);
    }
    file2.close();
    file2.clear();
    file2.open(syst_region_index);
    /*
    string sys_makeref = "makeblastdb -in "+WD+"buildup/L1.3.fasta  -dbtype nucl -parse_seqids";
    char *syst_makeref = new char[sys_makeref.length()+1];
    strcpy(syst_makeref, sys_makeref.c_str());
    system(syst_makeref);
     */
//parameters_end
    
//multiple threads
    
    string input_index;
    int line_index=0;
    for(int i=0;!file2.eof();i++){
        //file2>>input_index;
        file2>>input_index;
        if(CHR=="ALL"){
            line_index=i;
        }
        else if(CHR=="chr1"&&input_index=="1"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr2"&&input_index=="2"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr3"&&input_index=="3"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr4"&&input_index=="4"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr5"&&input_index=="5"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr6"&&input_index=="6"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr7"&&input_index=="7"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr8"&&input_index=="8"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr9"&&input_index=="9"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr10"&&input_index=="10"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr11"&&input_index=="11"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr12"&&input_index=="12"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr13"&&input_index=="13"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr14"&&input_index=="14"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr15"&&input_index=="15"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr16"&&input_index=="16"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr17"&&input_index=="17"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr18"&&input_index=="18"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr19"&&input_index=="19"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr20"&&input_index=="20"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr21"&&input_index=="21"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chr22"&&input_index=="2"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chrX"&&input_index=="X"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chrY"&&input_index=="Y"&&ref_n==37){
            line_index=i;
        }
        else if(CHR=="chrM"&&input_index=="M"&&ref_n==37){
            line_index=i;
        }
        else if(CHR==input_index&&ref_n==38){
            line_index=i;
        }
        file2>>input_index;
        file2>>input_index;
        
    }
    
    cout<<"THERE ARE "<<line_index<<" REGIONS TO COUNT."<<endl;
    cout<<"Pre-masking step & single read calling step is initiated."<<endl;
    
    file2.close();
    file2.clear();
    file2.open(syst_region_index);
    
    int NUM_circle;
    NUM_circle=(line_index/(NUM_threads+1))+1;

    pid_t p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30;
    
    string input;
    string chr;
    int start, end;
    
    for(int i=0;i!=line_index;){
        for(int j=0;j!=NUM_circle;j++){
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p1=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p2=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p3=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p4=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p5=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p6=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p7=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p8=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p9=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p10=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p11=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p12=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p13=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p14=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p15=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p16=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p17=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p18=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p19=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p20=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p21=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p22=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p23=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p24=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p25=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p26=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p27=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p28=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p29=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            if(i!=line_index){
                i++;
                file2>>chr;
                file2>>start;
                file2>>end;
                int chr_index=0;
                if(CHR=="ALL"){
                    chr_index=1;
                }
                else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chr22"&& chr=="2"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR=="chrM"&& chr=="M"&&ref_n==37){
                    chr_index=1;
                }
                else if(CHR== chr&&ref_n==38){
                    chr_index=1;
                }
                if(chr_index==1){
                    if((p30=fork())==0){
                        tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n);
                        //cout<<"And this is process "<<getpid()<<endl;
                        return 0;
                    }
                }
            }else break;
            waitpid(p1,NULL,0);
            waitpid(p2,NULL,0);
            waitpid(p3,NULL,0);
            waitpid(p4,NULL,0);
            waitpid(p5,NULL,0);
            waitpid(p6,NULL,0);
            waitpid(p7,NULL,0);
            waitpid(p8,NULL,0);
            waitpid(p9,NULL,0);
            waitpid(p10,NULL,0);
            waitpid(p11,NULL,0);
            waitpid(p12,NULL,0);
            waitpid(p13,NULL,0);
            waitpid(p14,NULL,0);
            waitpid(p15,NULL,0);
            waitpid(p16,NULL,0);
            waitpid(p17,NULL,0);
            waitpid(p18,NULL,0);
            waitpid(p19,NULL,0);
            waitpid(p20,NULL,0);
            waitpid(p21,NULL,0);
            waitpid(p22,NULL,0);
            waitpid(p23,NULL,0);
            waitpid(p24,NULL,0);
            waitpid(p25,NULL,0);
            waitpid(p26,NULL,0);
            waitpid(p27,NULL,0);
            waitpid(p28,NULL,0);
            waitpid(p29,NULL,0);
            waitpid(p30,NULL,0);

        }
            waitpid(p1,NULL,0);
            waitpid(p2,NULL,0);
            waitpid(p3,NULL,0);
            waitpid(p4,NULL,0);
            waitpid(p5,NULL,0);
            waitpid(p6,NULL,0);
            waitpid(p7,NULL,0);
            waitpid(p8,NULL,0);
            waitpid(p9,NULL,0);
            waitpid(p10,NULL,0);
            waitpid(p11,NULL,0);
            waitpid(p12,NULL,0);
            waitpid(p13,NULL,0);
            waitpid(p14,NULL,0);
            waitpid(p15,NULL,0);
            waitpid(p16,NULL,0);
            waitpid(p17,NULL,0);
            waitpid(p18,NULL,0);
            waitpid(p19,NULL,0);
            waitpid(p20,NULL,0);
            waitpid(p21,NULL,0);
            waitpid(p22,NULL,0);
            waitpid(p23,NULL,0);
            waitpid(p24,NULL,0);
            waitpid(p25,NULL,0);
            waitpid(p26,NULL,0);
            waitpid(p27,NULL,0);
            waitpid(p28,NULL,0);
            waitpid(p29,NULL,0);
            waitpid(p30,NULL,0);
 
        
    }
    
//merge and calling

    cout<<"Merging step is initiated."<<endl;
    //mkdir
    
    string WD_chr1=WD+"chr1/";
    string WD_chr2=WD+"chr2/";
    string WD_chr3=WD+"chr3/";
    string WD_chr4=WD+"chr4/";
    string WD_chr5=WD+"chr5/";
    string WD_chr6=WD+"chr6/";
    string WD_chr7=WD+"chr7/";
    string WD_chr8=WD+"chr8/";
    string WD_chr9=WD+"chr9/";
    string WD_chr10=WD+"chr10/";
    string WD_chr11=WD+"chr11/";
    string WD_chr12=WD+"chr12/";
    string WD_chr13=WD+"chr13/";
    string WD_chr14=WD+"chr14/";
    string WD_chr15=WD+"chr15/";
    string WD_chr16=WD+"chr16/";
    string WD_chr17=WD+"chr17/";
    string WD_chr18=WD+"chr18/";
    string WD_chr19=WD+"chr19/";
    string WD_chr20=WD+"chr20/";
    string WD_chr21=WD+"chr21/";
    string WD_chr22=WD+"chr22/";
    
    string sys_WD_chr1="mkdir "+WD+"chr1/";
    string sys_WD_chr2="mkdir "+WD+"chr2/";
    string sys_WD_chr3="mkdir "+WD+"chr3/";
    string sys_WD_chr4="mkdir "+WD+"chr4/";
    string sys_WD_chr5="mkdir "+WD+"chr5/";
    string sys_WD_chr6="mkdir "+WD+"chr6/";
    string sys_WD_chr7="mkdir "+WD+"chr7/";
    string sys_WD_chr8="mkdir "+WD+"chr8/";
    string sys_WD_chr9="mkdir "+WD+"chr9/";
    string sys_WD_chr10="mkdir "+WD+"chr10/";
    string sys_WD_chr11="mkdir "+WD+"chr11/";
    string sys_WD_chr12="mkdir "+WD+"chr12/";
    string sys_WD_chr13="mkdir "+WD+"chr13/";
    string sys_WD_chr14="mkdir "+WD+"chr14/";
    string sys_WD_chr15="mkdir "+WD+"chr15/";
    string sys_WD_chr16="mkdir "+WD+"chr16/";
    string sys_WD_chr17="mkdir "+WD+"chr17/";
    string sys_WD_chr18="mkdir "+WD+"chr18/";
    string sys_WD_chr19="mkdir "+WD+"chr19/";
    string sys_WD_chr20="mkdir "+WD+"chr20/";
    string sys_WD_chr21="mkdir "+WD+"chr21/";
    string sys_WD_chr22="mkdir "+WD+"chr22/";
    
    char *syst_WD_chr1 =new char[sys_WD_chr1.length()+1];
    strcpy(syst_WD_chr1, sys_WD_chr1.c_str());
    system(syst_WD_chr1);
    
    char *syst_WD_chr2 =new char[sys_WD_chr2.length()+1];
    strcpy(syst_WD_chr2, sys_WD_chr2.c_str());
    system(syst_WD_chr2);
    
    char *syst_WD_chr3 =new char[sys_WD_chr3.length()+1];
    strcpy(syst_WD_chr3, sys_WD_chr3.c_str());
    system(syst_WD_chr3);
    
    char *syst_WD_chr4 =new char[sys_WD_chr4.length()+1];
    strcpy(syst_WD_chr4, sys_WD_chr4.c_str());
    system(syst_WD_chr4);
    
    char *syst_WD_chr5 =new char[sys_WD_chr5.length()+1];
    strcpy(syst_WD_chr5, sys_WD_chr5.c_str());
    system(syst_WD_chr5);
    
    char *syst_WD_chr6 =new char[sys_WD_chr6.length()+1];
    strcpy(syst_WD_chr6, sys_WD_chr6.c_str());
    system(syst_WD_chr6);
    
    char *syst_WD_chr7 =new char[sys_WD_chr7.length()+1];
    strcpy(syst_WD_chr7, sys_WD_chr7.c_str());
    system(syst_WD_chr7);
    
    char *syst_WD_chr8 =new char[sys_WD_chr8.length()+1];
    strcpy(syst_WD_chr8, sys_WD_chr8.c_str());
    system(syst_WD_chr8);
    
    char *syst_WD_chr9 =new char[sys_WD_chr9.length()+1];
    strcpy(syst_WD_chr9, sys_WD_chr9.c_str());
    system(syst_WD_chr9);
    
    char *syst_WD_chr10 =new char[sys_WD_chr10.length()+1];
    strcpy(syst_WD_chr10, sys_WD_chr10.c_str());
    system(syst_WD_chr10);
    
    char *syst_WD_chr11 =new char[sys_WD_chr11.length()+1];
    strcpy(syst_WD_chr11, sys_WD_chr11.c_str());
    system(syst_WD_chr11);
    
    char *syst_WD_chr12 =new char[sys_WD_chr12.length()+1];
    strcpy(syst_WD_chr12, sys_WD_chr12.c_str());
    system(syst_WD_chr12);
    
    char *syst_WD_chr13 =new char[sys_WD_chr13.length()+1];
    strcpy(syst_WD_chr13, sys_WD_chr13.c_str());
    system(syst_WD_chr13);
    
    char *syst_WD_chr14 =new char[sys_WD_chr1.length()+1];
    strcpy(syst_WD_chr14, sys_WD_chr14.c_str());
    system(syst_WD_chr14);
    
    char *syst_WD_chr15 =new char[sys_WD_chr15.length()+1];
    strcpy(syst_WD_chr15, sys_WD_chr15.c_str());
    system(syst_WD_chr15);
    
    char *syst_WD_chr16 =new char[sys_WD_chr16.length()+1];
    strcpy(syst_WD_chr16, sys_WD_chr16.c_str());
    system(syst_WD_chr16);
    
    char *syst_WD_chr17 =new char[sys_WD_chr17.length()+1];
    strcpy(syst_WD_chr17, sys_WD_chr17.c_str());
    system(syst_WD_chr17);
    
    char *syst_WD_chr18 =new char[sys_WD_chr18.length()+1];
    strcpy(syst_WD_chr18, sys_WD_chr18.c_str());
    system(syst_WD_chr18);
    
    char *syst_WD_chr19 =new char[sys_WD_chr19.length()+1];
    strcpy(syst_WD_chr19, sys_WD_chr19.c_str());
    system(syst_WD_chr19);
    
    char *syst_WD_chr20 =new char[sys_WD_chr20.length()+1];
    strcpy(syst_WD_chr20, sys_WD_chr20.c_str());
    system(syst_WD_chr20);
    
    char *syst_WD_chr21 =new char[sys_WD_chr21.length()+1];
    strcpy(syst_WD_chr21, sys_WD_chr21.c_str());
    system(syst_WD_chr21);
    
    char *syst_WD_chr22 =new char[sys_WD_chr22.length()+1];
    strcpy(syst_WD_chr22, sys_WD_chr22.c_str());
    system(syst_WD_chr22);
    
    file2.close();
    file2.clear();
    file2.open(syst_region_index);
    
    //merge
    for(int i=0;i!=line_index;i++){
        file2>>chr;
        file2>>start;
        file2>>end;
        
        stringstream ss1, ss2;
        ss1 << start;
        string s_start =ss1.str();
        ss2 << end;
        string s_end =ss2.str();

        string sys_merge="cat "+WD+chr+"_"+s_start+"_"+s_end+"/calls.txt >> "+WD+"chr"+chr+"/calls.txt";
        //cout<<sys_merge<<endl;
        char *syst_merge = new char[sys_merge.length()+1];
        strcpy(syst_merge, sys_merge.c_str());
        system(syst_merge);
        
        string sys_merge_blastn="cat "+WD+chr+"_"+s_start+"_"+s_end+"/TSD_output.txt >> "+WD+"chr"+chr+"/TSD_output.txt";
        char *syst_merge_blastn = new char[sys_merge_blastn.length()+1];
        strcpy(syst_merge_blastn, sys_merge_blastn.c_str());
        system(syst_merge_blastn);
        
    }
    cout<<"Merging step completed."<<endl;
 /*
    //calling & collaps
    cout<<"Calling step is initiated."<<endl;
    
        if((p1=fork())==0){
            cout<<" Calling results in chr1."<<endl;
            calling(WD_chr1);
            //collaps(WD_chr1);
            return 0;
        }
        if((p2=fork())==0){
            cout<<" Calling results in chr2."<<endl;
            calling(WD_chr2);
            //collaps(WD_chr2);
            return 0;
        }
        if((p3=fork())==0){
            cout<<" Calling results in chr3."<<endl;
            calling(WD_chr3);
            //collaps(WD_chr3);
            return 0;
        }
        if((p4=fork())==0){
            cout<<" Calling results in chr4."<<endl;
            calling(WD_chr4);
            //collaps(WD_chr4);
            return 0;
        }
        if((p5=fork())==0){
            cout<<" Calling results in chr5."<<endl;
            calling(WD_chr5);
            //collaps(WD_chr5);
            return 0;
        }
        if((p6=fork())==0){
            cout<<" Calling results in chr6."<<endl;
            calling(WD_chr6);
            //collaps(WD_chr6);
            return 0;
        }
        if((p7=fork())==0){
            cout<<" Calling results in chr7."<<endl;
            calling(WD_chr7);
            //collaps(WD_chr7);
            return 0;
        }
        if((p8=fork())==0){
            cout<<" Calling results in chr8."<<endl;
            calling(WD_chr8);
            //collaps(WD_chr8);
            return 0;
        }
        if((p9=fork())==0){
            cout<<" Calling results in chr9."<<endl;
            calling(WD_chr9);
            //collaps(WD_chr9);
            return 0;
        }
        if((p10=fork())==0){
            cout<<" Calling results in chr10."<<endl;
            calling(WD_chr10);
            //collaps(WD_chr10);
            return 0;
        }
        if((p11=fork())==0){
            cout<<" Calling results in chr11."<<endl;
            calling(WD_chr11);
            //collaps(WD_chr11);
            return 0;
        }
        if((p12=fork())==0){
            cout<<" Calling results in chr12."<<endl;
            calling(WD_chr12);
            //collaps(WD_chr12);
            return 0;
        }
        if((p13=fork())==0){
            cout<<" Calling results in chr13."<<endl;
            calling(WD_chr13);
            //collaps(WD_chr13);
            return 0;
        }
        if((p14=fork())==0){
            cout<<" Calling results in chr14."<<endl;
            calling(WD_chr14);
            //collaps(WD_chr14);
            return 0;
        }
        if((p15=fork())==0){
            cout<<" Calling results in chr15."<<endl;
            calling(WD_chr15);
            //collaps(WD_chr15);
            return 0;
        }
        if((p16=fork())==0){
            cout<<" Calling results in chr16."<<endl;
            calling(WD_chr16);
            //collaps(WD_chr16);
            return 0;
        }
        if((p17=fork())==0){
            cout<<" Calling results in chr17."<<endl;
            calling(WD_chr17);
            //collaps(WD_chr17);
            return 0;
        }
        if((p18=fork())==0){
            cout<<" Calling results in chr18."<<endl;
            calling(WD_chr18);
            //collaps(WD_chr18);
            return 0;
        }
        if((p19=fork())==0){
            cout<<" Calling results in chr19."<<endl;
            calling(WD_chr19);
            //collaps(WD_chr19);
            return 0;
        }
        if((p20=fork())==0){
            cout<<" Calling results in chr20."<<endl;
            calling(WD_chr20);
            //collaps(WD_chr20);
            return 0;
        }
        if((p21=fork())==0){
            cout<<" Calling results in chr21."<<endl;
            calling(WD_chr21);
            //collaps(WD_chr21);
            return 0;
        }
        if((p22=fork())==0){
            cout<<" Calling results in chr22."<<endl;
            calling(WD_chr22);
            //collaps(WD_chr22);
            return 0;
        }
        waitpid(p1,NULL,0);
        waitpid(p2,NULL,0);
        waitpid(p3,NULL,0);
        waitpid(p4,NULL,0);
        waitpid(p5,NULL,0);
        waitpid(p6,NULL,0);
        waitpid(p7,NULL,0);
        waitpid(p8,NULL,0);
        waitpid(p9,NULL,0);
        waitpid(p10,NULL,0);
        waitpid(p11,NULL,0);
        waitpid(p12,NULL,0);
        waitpid(p13,NULL,0);
        waitpid(p14,NULL,0);
        waitpid(p15,NULL,0);
        waitpid(p16,NULL,0);
        waitpid(p17,NULL,0);
        waitpid(p18,NULL,0);
        waitpid(p19,NULL,0);
        waitpid(p20,NULL,0);
        waitpid(p21,NULL,0);
        waitpid(p22,NULL,0);
 */   
    cout<<"Almost finished there."<<endl;
    /*
    string sys_output1 = "cat "+WD_chr1+"collaps.txt >> "+WD+output;
    char *syst_output1 = new char[sys_output1.length()+1];
    strcpy(syst_output1, sys_output1.c_str());
    system(syst_output1);
    
    string sys_output2 = "cat "+WD_chr2+"collaps.txt >> "+WD+output;
    char *syst_output2 = new char[sys_output2.length()+1];
    strcpy(syst_output2, sys_output2.c_str());
    system(syst_output2);
    
    string sys_output3 = "cat "+WD_chr3+"collaps.txt >> "+WD+output;
    char *syst_output3 = new char[sys_output3.length()+1];
    strcpy(syst_output3, sys_output3.c_str());
    system(syst_output3);
    
    string sys_output4 = "cat "+WD_chr4+"collaps.txt >> "+WD+output;
    char *syst_output4 = new char[sys_output4.length()+1];
    strcpy(syst_output4, sys_output4.c_str());
    system(syst_output4);
    
    string sys_output5 = "cat "+WD_chr5+"collaps.txt >> "+WD+output;
    char *syst_output5 = new char[sys_output5.length()+1];
    strcpy(syst_output5, sys_output5.c_str());
    system(syst_output5);
    
    string sys_output6 = "cat "+WD_chr6+"collaps.txt >> "+WD+output;
    char *syst_output6 = new char[sys_output6.length()+1];
    strcpy(syst_output6, sys_output6.c_str());
    system(syst_output6);
    
    string sys_output7 = "cat "+WD_chr7+"collaps.txt >> "+WD+output;
    char *syst_output7 = new char[sys_output7.length()+1];
    strcpy(syst_output7, sys_output7.c_str());
    system(syst_output7);
    
    string sys_output8 = "cat "+WD_chr8+"collaps.txt >> "+WD+output;
    char *syst_output8 = new char[sys_output8.length()+1];
    strcpy(syst_output8, sys_output8.c_str());
    system(syst_output8);
    
    string sys_output9 = "cat "+WD_chr9+"collaps.txt >> "+WD+output;
    char *syst_output9 = new char[sys_output9.length()+1];
    strcpy(syst_output9, sys_output9.c_str());
    system(syst_output9);
    
    string sys_output10 = "cat "+WD_chr10+"collaps.txt >> "+WD+output;
    char *syst_output10 = new char[sys_output10.length()+1];
    strcpy(syst_output10, sys_output10.c_str());
    system(syst_output10);
    
    string sys_output11 = "cat "+WD_chr11+"collaps.txt >> "+WD+output;
    char *syst_output11 = new char[sys_output11.length()+1];
    strcpy(syst_output11, sys_output11.c_str());
    system(syst_output11);
    
    string sys_output12 = "cat "+WD_chr12+"collaps.txt >> "+WD+output;
    char *syst_output12 = new char[sys_output12.length()+1];
    strcpy(syst_output12, sys_output12.c_str());
    system(syst_output12);
    
    string sys_output13 = "cat "+WD_chr13+"collaps.txt >> "+WD+output;
    char *syst_output13 = new char[sys_output13.length()+1];
    strcpy(syst_output13, sys_output13.c_str());
    system(syst_output13);
    
    string sys_output14 = "cat "+WD_chr14+"collaps.txt >> "+WD+output;
    char *syst_output14 = new char[sys_output14.length()+1];
    strcpy(syst_output14, sys_output14.c_str());
    system(syst_output14);
    
    string sys_output15 = "cat "+WD_chr15+"collaps.txt >> "+WD+output;
    char *syst_output15 = new char[sys_output15.length()+1];
    strcpy(syst_output15, sys_output15.c_str());
    system(syst_output15);
    
    string sys_output16 = "cat "+WD_chr16+"collaps.txt >> "+WD+output;
    char *syst_output16 = new char[sys_output16.length()+1];
    strcpy(syst_output16, sys_output16.c_str());
    system(syst_output16);
    
    string sys_output17 = "cat "+WD_chr17+"collaps.txt >> "+WD+output;
    char *syst_output17 = new char[sys_output17.length()+1];
    strcpy(syst_output17, sys_output17.c_str());
    system(syst_output17);
    
    string sys_output18 = "cat "+WD_chr18+"collaps.txt >> "+WD+output;
    char *syst_output18 = new char[sys_output18.length()+1];
    strcpy(syst_output18, sys_output18.c_str());
    system(syst_output18);
    
    string sys_output19 = "cat "+WD_chr19+"collaps.txt >> "+WD+output;
    char *syst_output19 = new char[sys_output19.length()+1];
    strcpy(syst_output19, sys_output19.c_str());
    system(syst_output19);
    
    string sys_output20 = "cat "+WD_chr20+"collaps.txt >> "+WD+output;
    char *syst_output20 = new char[sys_output20.length()+1];
    strcpy(syst_output20, sys_output20.c_str());
    system(syst_output20);
    
    string sys_output21 = "cat "+WD_chr21+"collaps.txt >> "+WD+output;
    char *syst_output21 = new char[sys_output21.length()+1];
    strcpy(syst_output21, sys_output21.c_str());
    system(syst_output21);
    
    string sys_output22 = "cat "+WD_chr22+"collaps.txt >> "+WD+output;
    char *syst_output22 = new char[sys_output22.length()+1];
    strcpy(syst_output22, sys_output22.c_str());
    system(syst_output22);
    */
    cout<<"Final calls finished."<<endl;
    cout<<"Results are in "+WD+output<<endl;
    
}
