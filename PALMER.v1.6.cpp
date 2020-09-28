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

#include "extension/samview.h"

using namespace std;

int main(int argc, char *argv[]){

    Samview samview = Samview();
//parameters_start
    //std::ios::sync_with_stdio(false);
    //std::cin.tie(0);
    
    string T, WD, inputF, output, SP, ref, CHR, ref_fa, cus, cusin, tsd_pass, LL_len, s_cus_seq_len;
    
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
    int L_len=25;
    int seq_len=-1;
    int flag_cus_seq_len=0;
    output="output.txt";
    SP="Human";
    //T="LINE";
    CHR="ALL";
    int ref_n=0;
    
    string dir;
    dir=argv[0];
    
    for(int i=1;i!=argc;++i){
        if(strncmp(argv[i],"--workdir",9)==0){
            WD=argv[i+1];
            flag_wd=1;
        }
        if(strncmp(argv[i],"--chr",5)==0){
            CHR=argv[i+1];
        }
        
        if(strncmp(argv[i],"--L_len",7)==0){
            LL_len=argv[i+1];
            int L_len = 6023- std::stoi(LL_len);
            //L_len=L_len;
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
        if(strncmp(argv[i],"--len_custom_seq",16)==0){
            s_cus_seq_len=argv[i+1];
            
            seq_len = std::stoi(s_cus_seq_len);
            //strstream sss;
            //sss << string_cus_seq_len;
            //sss >> cus_seq_len;
            //cus_seq_len=stoi(string_cus_seq_len)
            flag_cus_seq_len=1;
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
    
    
    if((flag_T==0&&flag_cus==0)||(flag_T==1&&flag_cus==1)||flag_wd==0||flag_inputf==0||ref_n==0||flag_reffa==0||help==1||(flag_T==1&&ref_n==-1&&flag_cusin==0)||(flag_T==2&&flag_cus==0)||(flag_T==2&&flag_cus==1&&flag_cusin==0)||(flag_T==2&&flag_cus==1&&flag_cusin==1&&flag_tsd==1&flag_cus_seq_len==0)){
        if(flag_T==0&&flag_cus==0){
            cout<<"***ERROR*** PLEASE ASSIGN A MEI TYPE! LINE/ALU/SVA"<<endl;}
        if(flag_T==1&&flag_cus==1){
            cout<<"***ERROR*** PLEASE ASSIGN A MEI TYPE WITHOUT YOUR CUSTOMIZED SEQUENCE"<<endl;}
        if(flag_T==2&&flag_cus==0){
            cout<<"***ERROR*** PLEASE ASSIGN 'CUSTOMIZED' TYPE YOUR CUSTOMIZED SEQUENCE"<<endl;}
        if(flag_wd==0){
            cout<<"***ERROR*** PLEASE SET UP A WORKING DIRECTORY!"<<endl;}
        if(flag_inputf==0){
            cout<<"***ERROR*** PLEASE INPUT A FILE!"<<endl;}
        if(ref_n==0||flag_reffa==0){
            cout<<"***ERROR*** PLEASE ASSIGN A CORRECT REFERENCE version/fasta!"<<endl;}
        if(flag_T==1&&ref_n==-1&&flag_cusin==0){
            cout<<"***ERROR*** PLEASE ASSIGN A INDEX FILE FOR RUNNING MASKING MODULE ON YOUR MEI CALLING ON OTHER REFERENCE!"<<endl;}
        if(flag_T==2&&flag_cus==1&&flag_cusin==0){
        cout<<"***ERROR*** PLEASE ASSIGN 'CUSTOMIZED_INDEX' WHILE YOU CHOOSE YOUR CUSTOMIZED TYPE"<<endl;}
        if(flag_T==2&&flag_cus==1&&flag_cusin==1&&flag_tsd==1&flag_cus_seq_len==0){
        cout<<"***ERROR*** PLEASE ASSIGN 'CUSTOMIZD_SEQUENCE_LENGTH' WHILE YOU ACTIVATE TSD_FINDING FOR YOUR CUSTOMIZED TYPE"<<endl;}
        
        cout<<endl;
        cout<<"***WELCOME***"<<endl;
        cout<<"***PALMER:Pre-mAsking Long reads for Mobile Element inseRtion***"<<endl;
        cout<<"Version: 1.6.2"<<endl;
        cout<<"Presented by Weichen Zhou @ Mills Lab. Feb.28th.2019"<<endl;
        cout<<endl;
        cout<<"Usage:"<<endl;
        cout<<endl;
        cout<<"--input"<<endl;
        cout<<"         aligned long-read sequencing BAM file with directory path"<<endl;
        cout<<endl;
        cout<<"--workdir"<<endl;
        cout<<"         the user's working directory. Please follow the format /your/woking/directory/ !!don't forget the last '/'!!"<<endl;
        cout<<endl;
        cout<<"--ref_ver (options: hg19, GRCh37, GRCh38 or other)"<<endl;
        cout<<"         reference genome used for the aligned file ('other' option for the cusmized genome out of hg19, GRCh37 or GRCh38)"<<endl;
        cout<<endl;
        cout<<"--ref_fa"<<endl;
        cout<<"         indexed fasta file of reference genome fasta file with directory path used for the aligned bam file (wrong reference will cause error information)"<<endl;
        cout<<endl;
        cout<<"--type (options: LINE, ALU, SVA, or CUSTOMIZED (if you want to setup your costomized sequence))"<<endl;
        cout<<"         type of MEIs or other kinds of insertions to detect"<<endl;
        cout<<endl;
        cout<<"--chr (default: ALL (for whole genome); options: chromosome1, chromosome2, ...chromosomeY)"<<endl;
        cout<<"         chromosome name for PALMER to run (if running for whole genome, don't need to assign). !!The chromosome names should be consistent with the ones in reference genome version!! e.g. for GRCh37, to run PALMER on chromosome1, the option should be '1', while for GRCh38 it should be 'chr1'"<<endl;
        cout<<endl;
        
        cout<<"--custom_seq (default:no input)"<<endl;
        cout<<"         .fasta file with directory path to customize your insertion finding. e.g. NUMTs, MEIs in other species."<<endl;
        cout<<endl;
        cout<<"--custom_index (default:no input; if you have both '--ref_ver other' and '--type LINE/ALU/SVA', you must give PALMER a index file (format: \"CHR'\t'START'\t'END'\t'MEI_NAME'\n'\" for each MEI to be masked in each line) for masking module; if you have --custom_seq parameter without --custom_index, PALMER will work without masking step)"<<endl;
        cout<<"         index file with directory path to mask the genome for your insertion finding"<<endl;
        cout<<endl;
        cout<<"--TSD_finding (Fixed: TRUE for all MEIs ,or default: FALSE for CUSTOMIZED insertion)"<<endl;
        cout<<"         whether to run TSD motif finding module for your insertion calling"<<endl;
        cout<<endl;
        
        cout<<"--len_custom_seq (MUST set up when activate TSD_finding for CUSTOMIZED insertion, otherwise CLOSED)"<<endl;
        cout<<"         interger value for the length of your customized sequence WITHOUT polyA tact"<<endl;
        cout<<endl;
        
        cout<<"--L_len (default: 25bp)"<<endl;
        cout<<"         the minimum length of putative LINE-1 aligned to L1.3 sequences"<<endl;
        cout<<endl;
        
        cout<<"--output (default: output)"<<endl;
        cout<<"         the prefix of the output file"<<endl;
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
    
    vector<string> dir_conv;
    dir_conv.clear();
    FILE *pp =popen(sys_dir.c_str(),"r");
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
    
    //cout<<T<<endl;
    //cout<<seq_len<<endl;
    //cout<<flag_cus_seq_len<<endl;
    //getchar();
    
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
        
        if(samview.SamViewHeaderOnly(inputF.c_str())!=0){
            cout << "CANNOT READ HEADER INFO FROM FILE<" << inputF << ">" << endl;
            return 1;
        };
        
        sys_region_index=buildup+"region.split.index";
    }
    ifstream file2;
    
    if(ref_n==-1){
        
        ofstream file93;
        file93.open(sys_region_index.c_str(),ios::trunc);
        
        int bin=1000000;
        for(int i=0;i!=samview.headerChr.size();++i){
            int j=i;
            int sec;
                //last;
            sec=int(samview.headerLength[j]/bin);
                //last=len_inde[j]%bin;
            for(int k=0;k!=sec;k++){
                file93<<samview.headerChr[i]<<'\t'<<(1+k*bin)<<'\t'<<(k*bin+bin)<<endl;
            }
            file93<<samview.headerChr[i]<<'\t'<<(1+sec*bin)<<'\t'<<samview.headerLength[j]<<endl;
            
        }
        //ifstream file2;
        file93.close();
        file93.clear();
    }
    
    else {
    //string sys_region_index=buildup+"region.split.index";
        file2.open(sys_region_index.c_str());
        
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
    file2.open(sys_region_index.c_str());
    string input_index;
    int line_index=0;
    for(int i=1;!file2.eof();){
        //file2>>input_index;
        file2>>input_index;
        if(CHR=="ALL"){
            line_index=i;
            ++i;
        }
        else if(CHR==input_index){
            line_index=i;
            ++i;
        }
        file2>>input_index;
        file2>>input_index;
    }
    
    //line_index=line_index+1;
    
    cout<<"THERE ARE "<<line_index<<" REGIONS TO COUNT."<<endl;
    cout<<"Pre-masking step & single read calling step is initiated."<<endl;
    
    file2.close();
    file2.clear();
    file2.open(sys_region_index.c_str());
    
    int NUM_circle;
    NUM_circle=(line_index/(NUM_threads+1))+1;

    string input;
    string chr;
    int start, end;

    //no multithread
    for(int i=0;i!=line_index;){
        file2>>chr;
        file2>>start;
        file2>>end;

        if ( CHR=="ALL" || CHR==chr ) {
            ++i;
            tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n, direc, ref_fa, flag_tsd, L_len, seq_len, &samview);
        }
    }
   
    //merge and calling
    cout<<"Merging step is initiated."<<endl;
    //mkdir
    /*
    for(int i=0;i!=line_chr;++i){
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
    file2.open(sys_region_index.c_str());
    
    //merge
    
    string sys_final_title = WD+output+"_calls.txt";
    ofstream file3;
    file3.open(sys_final_title.c_str(),ios::trunc);
    
    file3<<"cluster_id"<<'\t'<<"chr start1"<<'\t'<<"start2"<<'\t'<<"end1"<<'\t'<<"end2"<<'\t'<<"start1_inVariant"<<'\t'<<"start2_inVariant"<<'\t'<<"end1_inVariant"<<'\t'<<"end2_inVariant"<<'\t'<<"Confident_supporting_reads"<<'\t'<<"Potential_supporting_reads"<<'\t'<<"Ptential_segmental_supporting_reads"<<'\t'<<"orientation"<<'\t'<<"polyA-tail_size"<<'\t'<<"5'_TSD_size"<<'\t'<<"3'_TSD_size"<<'\t'<<"Predicted_transD_size"<<'\t'<<"Has_5'_inverted_sequence?"<<'\t'<<"5'_inverted_seq_end"<<'\t'<<"5'_seq_start"<<endl;
    
    string sys_final_tsd_title = WD+output+"_TSD_reads.txt";
    ofstream file31;
    file31.open(sys_final_tsd_title.c_str(),ios::trunc);
    
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
            system(sys_final.c_str());
            
            string sys_final_tsd="cat "+WD+chr+"_"+s_start+"_"+s_end+"/TSD_output.txt >> "+sys_final_tsd_title;
            system(sys_final_tsd.c_str());
            ++i;
        }
    }
    
    cout<<"Merging step completed."<<endl;
    
    //colapse for the redundant calls
    
    //kmer
    
    cout<<"Final calls finished."<<endl;
    cout<<"Results are in "+WD+output<<endl;
}
