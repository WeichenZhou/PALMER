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

using namespace std;

int main(int argc, char *argv[]){

//parameters_start
    
    string T, WD, inputF, output, SP, ref, CHR, ref_fa, cus, cusin, tsd_pass;
    //int NUM_threads=30;
    int NUM_threads=10;
    ifstream file1;
    ifstream file11;
    //ifstream file12;
    ifstream file13;
    //ofstream file2;
    int flag_wd=0;
    int flag_inputf=0;
    int flag_T=0;
    int flag_reffa=0;
    int flag_cus=0;
    int flag_cusin=0;
    int flag_tsd=1;
    int help=0;
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
        if(strncmp(argv[i],"--ref_ver",9)==0){
            ref=argv[i+1];
            if(ref=="GRCh37") ref_n=37;
            else if(ref=="GRCh38") ref_n=38;
            else if(ref=="hg19") ref_n=19;
            else if(ref=="other") ref_n=-1;
            else {
                cout<<"PLEASE INPUT A CORRECT CATEGORY OF REFERENCE :("<<endl;
            }
        }
        if(strncmp(argv[i],"--type",6)==0){
            T=argv[i+1];
            flag_T=1;
            if(T=="CUSTOMIZED"){
                flag_T=2;
            }
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
        if(strncmp(argv[i],"--ref_fa",8)==0){
            flag_reffa=1;
            file11.open(argv[i+1],ios::in | ios::binary);
            if (!file11.is_open())
            {
                cout <<"CANNOT OPEN REFERENCE FILE :("  << endl;
                cout <<"PLEASE ASSIGN A REFERENCE FILE."<< endl;
                //exit(1);
                flag_reffa=0;
            }
            ref_fa=argv[i+1];
        }
        
        if(strncmp(argv[i],"--custom_seq",12)==0){
            flag_cus=1;
            /*
            file12.open(argv[i+1],ios::in | ios::binary);
            if (!file12.is_open())
            {
                cout <<"CANNOT OPEN CUSTOMIZED FILE :("<< endl;
                cout<<"PLEASE ASSIGN A CUSTOMIZED FASTA FILE."<<endl;
                //exit(1);
                flag_cus=0;
            }*/
            cus=argv[i+1];
            
        }
        if(strncmp(argv[i],"--custom_index",14)==0){
            flag_cusin=1;
            if(flag_cus==1){
                 file13.open(argv[i+1],ios::in | ios::binary);
                 if (!file13.is_open())
                 {
                 cout <<"You are not assigning any index file for masking module."<< endl;
                 cout <<"Custmized finding will initiate without masking module."<<endl;
                 //exit(1);
                 flag_cusin=0;
                 }
            }
            cusin=argv[i+1];
        }
        if(strncmp(argv[i],"--TSD_finding",13)==0){
            
            tsd_pass=argv[i+1];
        }
        if(strncmp(argv[i],"--help",6)==0){
            
            help=1;
        }
    }
    
    if(T=="ALU"||T=="SVA"||T=="LINE"){
        if(tsd_pass=="FALSE"){
            cout<<"RETROTRANSPOSOTION EVENTS REQUIRE TSD MOTIF FINDING. GO ON WITH TRUE ANYWAY."<<endl;
        }
        //else if(tsd_pass=="TRUE")
    }
    else if (T=="CUSTOMIZED"){
        flag_tsd=0;
        if(tsd_pass=="TRUE"){
            flag_tsd=1;
        }
        if(flag_tsd==1){
            cout<<"Customized insertion finding will go on with TSD finding module."<<endl;
        }
    }
    
    
    if((flag_T==0&&flag_cus==0)||(flag_T==1&&flag_cus==1)||flag_wd==0||flag_inputf==0||ref_n==0||flag_reffa==0||help==1||(flag_T==1&&ref_n==-1&&flag_cusin==0)||(flag_T==2&&flag_cus==0)){
        if(flag_T==0&&flag_cus==0){
            cout<<"***ERROR*** PLEASE ASSIGN A MEI TYPE! LINE/ALU/SVA"<<endl;}
        if(flag_T==1&&flag_cus==1){
            cout<<"***ERROR*** PLEASE ASSIGN A MEI TYPE WITHOUT YOUR CUSTOMIZED SEQUENCE"<<endl;}
        if(flag_T==2&&flag_cus==0){
            cout<<"***ERROR*** PLEASE ASSIGN 'CUSTOMIZED' TYPE YOUR CUSTOMIZED SEQUENCE"<<endl;}
        if(flag_wd==0){
            cout<<"***ERROR*** PLEASE SET UP A WORKING DIRECOTY!"<<endl;}
        if(flag_inputf==0){
            cout<<"***ERROR*** PLEASE INPUT A FILE!"<<endl;}
        if(ref_n==0||flag_reffa==0){
            cout<<"***ERROR*** PLEASE ASSIGN A CORRECT REFERENCE version/fasta!"<<endl;}
        if(flag_T==1&&ref_n==-1&&flag_cusin==0){
            cout<<"***ERROR*** PLEASE ASSIGN A INDEX FILE FOR RUNNING MASKING MODULE ON YOUR MEI CALLING ON OTHER REFERENCE!"<<endl;}
        cout<<endl;
        cout<<"***WELCOME***"<<endl;
        cout<<"***PALMER:Pre-mAsking Long reads for Mobile Element inseRtion***"<<endl;
        cout<<"Version: 1.4"<<endl;
        cout<<"Presented by Weichen Zhou @ Mills Lab. Feb.28th.2019"<<endl;
        cout<<endl;
        cout<<"Usage:"<<endl;
        cout<<endl;
        cout<<"--input"<<endl;
        cout<<"         aligned long-read sequencing BAM file with directory path"<<endl;
        cout<<endl;
        cout<<"--workdir"<<endl;
        cout<<"         the user's working directory"<<endl;
        cout<<endl;
        cout<<"--ref_ver (options: hg19, GRCh37, GRCh38 or other)"<<endl;
        cout<<"         reference genome used for the aligned file ('other' option for the cusmized genome out of hg19, GRCh37 or GRCh38)"<<endl;
        cout<<endl;
        cout<<"--ref_fa"<<endl;
        cout<<"         indexed fasta file of reference genome fasta file with directory path used for the aligned bam file (wrong reference will cause error infromation)"<<endl;
        cout<<endl;
        cout<<"--type (options: LINE, ALU, SVA, or CUSTOMIZED (if you want to setup your costomized sequence))"<<endl;
        cout<<"         type of MEIs or other kind of insertion to detect"<<endl;
        cout<<endl;
        cout<<"--chr (default: ALL (for whole genome); options: chromosome1, chromosome2, ...chromosomeY)"<<endl;
        cout<<"         chromosome name for PALMER to run (if running for whole genome, don't need to assign). !!The chromosome names should be consistent with the ones in reference genome version!! e.g. for GRCh37, to run PALMER on chromosome1, the option should be '1', while for GRCh38 it should be 'chr1'"<<endl;
        cout<<endl;
        
        cout<<"--custom_seq (default:no input)"<<endl;
        cout<<"         .fasta file with directory path to customize your insertion finding"<<endl;
        cout<<endl;
        cout<<"--custom_index (default:no input; if you have both '--ref_ver other' and '--type LINE/ALU/SVA', you must give PALMER a index file (format: \"CHR'\t'START'\t'END'\t'MEI_NAME'\n'\" for each MEI to be masked in each line) for masking module; if you have --custom_seq parameter without --custom_index, PALMER will work without masking step)"<<endl;
        cout<<"         index file with directory path to mask the genome for your insertion finding"<<endl;
        cout<<endl;
        cout<<"--TSD_finding (Fixed:TRUE for all MEIs ,or default: FALSE for CUSTOMIZED insertion)"<<endl;
        cout<<"         whether to run TSD motif finding module for your insertion calling"<<endl;
        cout<<endl;
        
        cout<<"--output (default: output)"<<endl;
        cout<<"         prefix of output file"<<endl;
        cout<<endl;
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
    
    string sys_region_index,sys_line_region, sys_region_index_chr, sys_region_index_length;
    
    if(flag_cusin==1){
        sys_line_region=cusin;
    }
    else {
        sys_line_region="NULL";
    }
    
    
    if(T=="CUSTOMIZED"){
        T=cus;
    }
    
    
    if(ref_n==37){
        
        sys_region_index=buildup+"region.split.index.GRCh37";
        
        sys_region_index_chr=buildup+"region.split.index.chr.GRCh37";
        sys_region_index_length=buildup+"region.split.index.length.GRCh37";
        
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
        
        sys_region_index_chr=buildup+"region.split.index.chr.GRCh38";
        sys_region_index_length=buildup+"region.split.index.length.GRCh38";
        
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
    else if(ref_n==19){
        sys_region_index=buildup+"region.split.index.hg19";
        
        sys_region_index_chr=buildup+"region.split.index.chr.hg19";
        sys_region_index_length=buildup+"region.split.index.length.hg19";
        
        if(T=="LINE"){
            sys_line_region=buildup+"LINEs.regions.hg19";
        }
        if(T=="ALU"){
            sys_line_region=buildup+"Alu.regions.hg19";
        }
        if(T=="SVA"){
            sys_line_region=buildup+"SVA.regions.hg19";
        }
    }
    
    
//reference_index buildup
    else if(ref_n==-1){
        //sys_line_region=cusin;
        
        
        buildup=WD+"index/";
        
        string sys_build;
        sys_build="mkdir "+buildup;
        
        char *syst_build = new char[sys_build.length()+1];
        strcpy(syst_build, sys_build.c_str());
        system(syst_build);
        
        string sys1;
        sys1="samtools view "+inputF+" -H |grep \"@SQ\" | awk -F \":|\t\" '{ print $3}' >"+buildup+"chr.list";
        
        string sys2;
        sys2="samtools view "+inputF+" -H |grep \"@SQ\" | awk -F \":|\t\" '{ print $5}' >"+buildup+"length.list";
        
        char *syst1 = new char[sys1.length()+1];
        strcpy(syst1, sys1.c_str());
        system(syst1);
        
        char *syst2 = new char[sys2.length()+1];
        strcpy(syst2, sys2.c_str());
        system(syst2);
        
        sys_region_index_chr=buildup+"chr.list";
        sys_region_index_length=buildup+"length.list";
        sys_region_index=buildup+"region.split.index";
    }
    //cout<<sys_region_index<<" "<<sys_line_region<<endl;
 //original
    char *syst_region_index_chr =new char[sys_region_index_chr.length()+1];
    strcpy(syst_region_index_chr, sys_region_index_chr.c_str());
    char *syst_region_index_length =new char[sys_region_index_length.length()+1];
    strcpy(syst_region_index_length, sys_region_index_length.c_str());
    
    ifstream file91;
    ifstream file92;
    file91.open(syst_region_index_chr);
    file92.open(syst_region_index_length);
    
    int line_chr;
    int line_len;
    string input_inde;
    for(int i=0;!file91.eof();i++){
        getline(file91,input_inde);
        line_chr=i;
    }
    for(int i=0;!file92.eof();i++){
        getline(file92,input_inde);
        line_len=i;
    }
    
    //cout<<"line_chr="<<line_chr<<endl;
    //cout<<"line_length="<<line_len<<endl;
    
    file92.close();
    file91.close();
    file92.clear();
    file91.clear();
    
    file91.open(syst_region_index_chr);
    file92.open(syst_region_index_length);
    
    string chr_inde[line_chr];
    int len_inde[line_len];
    
    for(int i=0;i!=line_chr;i++){
        file91>>chr_inde[i];
    }
    for(int i=0;i!=line_len;i++){
        file92>>len_inde[i];
    }
    
    file92.close();
    file91.close();
    file92.clear();
    file91.clear();
    ifstream file2;
    
    //char *syst_region_index ;
    char *syst_region_index =new char[sys_region_index.length()+1];
    strcpy(syst_region_index, sys_region_index.c_str());
    
    
    if(ref_n==-1){
        
        ofstream file93;
        file93.open(syst_region_index,ios::trunc);
        
        int bin=1000000;
        for(int i=0;i!=line_chr;i++){
            int j=i;
            int sec;
                //last;
            sec=int(len_inde[j]/bin);
                //last=len_inde[j]%bin;
            for(int k=0;k!=sec;k++){
                file93<<chr_inde[i]<<'\t'<<(1+k*bin)<<'\t'<<(k*bin+bin)<<endl;
            }
            file93<<chr_inde[i]<<'\t'<<(1+sec*bin)<<'\t'<<len_inde[j]<<endl;
            
        }
        //ifstream file2;
        file93.close();
        file93.clear();
    }
    
    else {
    //string sys_region_index=buildup+"region.split.index";
        
        file2.open(syst_region_index);
        
        if (!file2.is_open())
        {
            cout <<"CANNOT OPEN REF INDEX FILE"<< endl;
            //cout<<"YEA"<<endl;
            exit(1);
        }
        file2.close();
        file2.clear();
        //file2.open(syst_region_index);
    }
    
    
//parameters_end
    
//multiple threads
    file2.open(syst_region_index);
    string input_index;
    int line_index=0;
    for(int i=1;!file2.eof();){
        //file2>>input_index;
        file2>>input_index;
        if(CHR=="ALL"){
            line_index=i;
            i++;
        }
        else if(CHR==input_index){
            line_index=i;
            i++;
        }
        file2>>input_index;
        file2>>input_index;
    }
    
    //line_index=line_index+1;
    
    cout<<"THERE ARE "<<line_index<<" REGIONS TO COUNT."<<endl;
    cout<<"Pre-masking step & single read calling step is initiated."<<endl;
    
    file2.close();
    file2.clear();
    file2.open(syst_region_index);
    
    int NUM_circle;
    NUM_circle=(line_index/(NUM_threads+1))+1;

    //pid_t p1, p2, p3, p4, p5, p6, p7, p8, p9, p10;
    //, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30;
    
    string input;
    string chr;
    int start, end;

//no multithread
    for(int i=0;i!=line_index;){
        //cout<<"right call"<<endl;
        //cout<<chr<<endl;
        
        file2>>chr;
        file2>>start;
        file2>>end;
        int chr_index=0;
        if(CHR=="ALL"){
            chr_index=1;
        }
        else if(CHR==chr){
            chr_index=1;
            
        }
        if(chr_index==1){
            i++;
            
            //getchar();
            //****
            tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n, direc, ref_fa, flag_tsd);
        }
        /*ver1.3
        if(chr_index==1){
            tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n, direc, ref_file);
        }*/
    }
    
   
//merge and calling

    cout<<"Merging step is initiated."<<endl;
    //mkdir
    /*
    for(int i=0;i!=line_chr;i++){
        string WD_chr=WD+chr_inde[i]+"/";
        string WD_chr="mkdir "+WD+chr_inde[i]+"/";
        
        if(CHR==chr_inde[i]||CHR=="ALL"){
            char *syst_WD_chr =new char[sys_WD_chr.length()+1];
            strcpy(syst_WD_chr, sys_WD_chr.c_str());
            system(syst_WD_chr);
        }
    }
    */
    file2.close();
    file2.clear();
    file2.open(syst_region_index);
    
    //merge
    
    string sys_final_title = WD+output+"_calls.txt";
    char *syst_final_title = new char[sys_final_title.length()+1];
    strcpy(syst_final_title, sys_final_title.c_str());
    ofstream file3;
    file3.open(syst_final_title,ios::trunc);
    
    file3<<"cluster_id"<<'\t'<<"chr start1"<<'\t'<<"start2"<<'\t'<<"end1"<<'\t'<<"end2"<<'\t'<<"start1_inVariant"<<'\t'<<"start2_inVariant"<<'\t'<<"end1_inVariant"<<'\t'<<"end2_inVariant"<<'\t'<<"Confident_supporting_reads"<<'\t'<<"Potential_supporting_reads"<<'\t'<<"Ptential_segmental_supporting_reads"<<'\t'<<"orientation"<<'\t'<<"polyA-tail_size"<<'\t'<<"5'_TSD_size"<<'\t'<<"3'_TSD_size"<<'\t'<<"Predicted_transD_size"<<'\t'<<"Has_5'_inverted_sequence?"<<'\t'<<"5'_inverted_seq_end"<<'\t'<<"5'_seq_start"<<endl;
    
    string sys_final_tsd_title = WD+output+"_TSD_reads.txt";
    char *syst_final_tsd_title = new char[sys_final_tsd_title.length()+1];
    strcpy(syst_final_tsd_title, sys_final_tsd_title.c_str());
    ofstream file31;
    file31.open(syst_final_tsd_title,ios::trunc);
    
    file31<<"cluster_id"<<'\t'<<"read_name.info"<<'\t'<<"5'_TSD"<<'\t'<<"3'_TSD"<<'\t'<<"Predicted_transD"<<'\t'<<"Unique_26mer_at_5'junction"<<'\t'<<"Whole_insertion_seq"<<endl;
    
    for(int i=0;i!=line_index;){
        file2>>chr;
        file2>>start;
        file2>>end;
        
        stringstream ss1, ss2;
        ss1 << start;
        string s_start =ss1.str();
        ss2 << end;
        string s_end =ss2.str();
        
        if(CHR==chr){
            string sys_final="cat "+WD+chr+"_"+s_start+"_"+s_end+"/calls.txt >> "+sys_final_title;
            //cout<<sys_final<<endl;
            char *syst_final = new char[sys_final.length()+1];
            strcpy(syst_final, sys_final.c_str());
            system(syst_final);
            
            string sys_final_tsd="cat "+WD+chr+"_"+s_start+"_"+s_end+"/TSD_output.txt >> "+sys_final_tsd_title;
            //cout<<sys_final<<endl;
            char *syst_final_tsd = new char[sys_final_tsd.length()+1];
            strcpy(syst_final_tsd, sys_final_tsd.c_str());
            system(syst_final_tsd);
            i++;
        }
    }
    
    cout<<"Merging step completed."<<endl;
    cout<<"Final calls finished."<<endl;
    cout<<"Results are in "+WD+output<<endl;
    
}
