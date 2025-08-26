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
#include <numeric>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
using namespace std;

int fp_ex(string WD_dir, string fasta, string chr, string t, int tsd_index){
    
    //std::ios::sync_with_stdio(false);
    //std::cin.tie(0);
    int BIN_buff=10;
    int BIN_5=50;
    int J_BIN=50;
    int BIN_3=3000;
    int J_BIN_mer=13;
    
    if (t=="ALU"){
        BIN_3=150;
        //BIN_5=50;
    }
    else if (t=="SVA"){
        BIN_3=2500;
        BIN_5=2000;
    }
    else if (t=="HERVK"){
        BIN_3=50;
        BIN_5=50;
        tsd_index=0;
    }
    
    ifstream file2;
    //ifstream file3;
    
    string sys_junc = WD_dir+"read_result_junction.txt";
    char *syst_junc = new char[sys_junc.length()+1];
    strcpy(syst_junc, sys_junc.c_str());
    file2.open(syst_junc);
    
    
    if (!file2.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'read_result_junction.txt'"<< endl;
        //exit(1);
        return 0;
    }
    
    int line;
    string input;
    for(int i=0;!file2.eof();++i){
        getline(file2,input);
        line=i;
    }
    
    file2.close();
    file2.clear();
    file2.open(syst_junc);
    
    string **info;
    info=new string *[line];
    for(int i=0;i!=line;++i) info[i]=new string[5];
    
    int **loc;
    loc=new int*[line];
    for(int i=0;i!=line;++i) loc[i]=new int[7];
    
    int **loc_TP;
    loc_TP=new int*[line];
    for(int i=0;i!=line;++i) loc_TP[i]=new int[7];
    
    for(int i=0;i!=line;++i){
        file2>>info[i][0];  //probe name
        file2>>loc[i][0];   //L_loc
        file2>>loc[i][1];   //L_loc
        file2>>loc[i][2];   //R_loc
        file2>>loc[i][3];   //R_loc
        file2>>loc[i][4];   //G_loc
        file2>>loc[i][5];   //G_loc
        file2>>info[i][1];  //chr
        file2>>info[i][2];  //orientation
        file2>>loc[i][6];  //read length
        file2>>info[i][3]; //junc 5' seq
        file2>>info[i][4]; //junc 3' seq
        file2>>loc_TP[i][0];
        file2>>loc_TP[i][1];
        file2>>loc_TP[i][2];
        file2>>loc_TP[i][3];
        file2>>loc_TP[i][4];
        file2>>loc_TP[i][5];
        file2>>loc_TP[i][6];
    }
    
    ifstream file1;
    
    string sys_input = WD_dir+"TSD_blastn_pre.txt";
    char *syst_input = new char[sys_input.length()+1];
    strcpy(syst_input, sys_input.c_str());
    file1.open(syst_input);
    
    ofstream file11;
    string sys_output = WD_dir+"TSD_blastn.txt";
    char *syst_output = new char[sys_output.length()+1];
    strcpy(syst_output, sys_output.c_str());
    file11.open(syst_output);
    
    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'TSD_blastn_pre.txt'"<< endl;
        //exit(1);
        return 0;
    }
    
    
    int line_tsd;
    for(int i=0;!file1.eof();++i){
        getline(file1,input);
        line_tsd=i;
    }
    file1.close();
    file1.clear();
    file1.open(syst_input);
    
    string *info_tsd;
    info_tsd= new string[line_tsd];
    
    int **loc_tsd;
    loc_tsd=new int*[line_tsd];
    for(int i=0;i!=line_tsd;++i) loc_tsd[i]=new int[7];
    
    string **loc_tsd_fp;
    loc_tsd_fp=new string*[line_tsd];
    for(int i=0;i!=line_tsd;++i) loc_tsd_fp[i]=new string[2];
    
    string *kmer_tsd;
    kmer_tsd= new string[line_tsd];
    
    for(int i=0;i!=line_tsd;++i){
        file1>>info_tsd[i];
        file1>>loc_tsd[i][0];
        file1>>loc_tsd[i][1];
        file1>>loc_tsd[i][2];
        file1>>loc_tsd[i][3];
        file1>>loc_tsd[i][4];
        file1>>loc_tsd[i][5];
        loc_tsd[i][6]=-1;
        //loc_tsd[i][7]=-1;
        //loc_tsd[i][8]=-1;
        loc_tsd_fp[i][0]="";
        loc_tsd_fp[i][1]="";
        kmer_tsd[i]="";
    }

    
    //cout<<"ready to process this"<<endl;
//FP_ex module
    if(tsd_index==1){
        for(int i=0;i!=line;++i){
            
            string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
            stringstream ss_0;
            ss_0.clear();
            ss_0<<loc[i][0];
            loc_0=ss_0.str();
            stringstream ss_1;
            ss_1.clear();
            ss_1<<loc[i][1];
            loc_1=ss_1.str();
            stringstream ss_2;
            ss_2.clear();
            ss_2<<loc[i][2];
            loc_2=ss_2.str();
            stringstream ss_3;
            ss_3.clear();
            ss_3<<loc[i][3];
            loc_3=ss_3.str();
            stringstream ss_4;
            ss_4.clear();
            ss_4<<loc[i][4];
            loc_4=ss_4.str();
            stringstream ss_5;
            ss_5.clear();
            ss_5<<loc[i][5];
            loc_5=ss_5.str();
            
            string loc_TP_0, loc_TP_1, loc_TP_2, loc_TP_3, loc_TP_4, loc_TP_5, loc_TP_6;
            stringstream ss_TP_0;
            ss_TP_0.clear();
            ss_TP_0<<loc_TP[i][0];
            loc_TP_0=ss_TP_0.str();
            stringstream ss_TP_1;
            ss_TP_1.clear();
            ss_TP_1<<loc_TP[i][1];
            loc_TP_1=ss_TP_1.str();
            stringstream ss_TP_2;
            ss_TP_2.clear();
            ss_TP_2<<loc_TP[i][2];
            loc_TP_2=ss_TP_2.str();
            stringstream ss_TP_3;
            ss_TP_3.clear();
            ss_TP_3<<loc_TP[i][3];
            loc_TP_3=ss_TP_3.str();
            stringstream ss_TP_4;
            ss_TP_4.clear();
            ss_TP_4<<loc_TP[i][4];
            loc_TP_4=ss_TP_4.str();
            stringstream ss_TP_5;
            ss_TP_5.clear();
            ss_TP_5<<loc_TP[i][5];
            loc_TP_5=ss_TP_5.str();
            stringstream ss_TP_6;
            ss_TP_6.clear();
            ss_TP_6<<loc_TP[i][6];
            loc_TP_6=ss_TP_6.str();
            
            string seq_index;
            
            seq_index=info[i][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
     
            
    //ref pull out for 5mer
            /*
            int start, end;
            start=loc[i][4]-loc[i][6];
            if(start<=0){
                start=1;
            }
            end=loc[i][5]+loc[i][6];
            
            if(loc[i][6]<6000){
                start=loc[i][4]-6000;
                end=loc[i][5]+6000;
            }
            
            string s_start, s_end;
            stringstream ss_s;
            ss_s.clear();
            ss_s<<start;
            s_start=ss_s.str();
            stringstream ss_e;
            ss_e.clear();
            ss_e<<end;
            s_end=ss_e.str();
            */
            //No need ref.fasta for now
            /*
            ofstream file20;
            string ref_file;
            ref_file=WD_dir+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str()+".ref.fasta";
            //cout<<ref_file<<endl;
            char *syst_ref_file = new char[ref_file.length()+1];
            strcpy(syst_ref_file, ref_file.c_str());
            file20.open(syst_ref_file);
            
            string sys_ref;
            sys_ref="samtools faidx "+fasta+" "+chr+":"+s_start+"-"+s_end+" > "+ref_file;
            //cout<<sys_ref<<endl;
            
            char *syst_ref = new char[sys_ref.length()+1];
            strcpy(syst_ref, sys_ref.c_str());
            system(syst_ref);
            //getchar();
            */
//Four -fold or seven -fold
    //ref pull out for fake junc
            int start_junc, end_junc;
            start_junc=loc[i][4]-J_BIN*7;
            if(start_junc<=0){
                start_junc=1;
            }
            end_junc=loc[i][5]+J_BIN*7;
            
            string s_start_junc, s_end_junc;
            stringstream ss_s_junc;
            ss_s_junc.clear();
            ss_s_junc<<start_junc;
            s_start_junc=ss_s_junc.str();
            stringstream ss_e_junc;
            ss_e_junc.clear();
            ss_e_junc<<end_junc;
            s_end_junc=ss_e_junc.str();
            
            ofstream file21;
            string ref_junc_file;
            ref_junc_file=WD_dir+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str()+".junc.ref.fasta";
            //cout<<ref_file<<endl;
            char *syst_ref_junc_file = new char[ref_junc_file.length()+1];
            strcpy(syst_ref_junc_file, ref_junc_file.c_str());
            file21.open(syst_ref_junc_file);
            
            string sys_ref_junc;
            sys_ref_junc="samtools faidx "+fasta+" "+chr+":"+s_start_junc+"-"+s_end_junc+" > "+ref_junc_file;
            //cout<<sys_ref<<endl;
            
            char *syst_ref_junc = new char[sys_ref_junc.length()+1];
            strcpy(syst_ref_junc, sys_ref_junc.c_str());
            system(syst_ref_junc);

    //FP exclude module
            //3' 26mer identify based on TSD
            //FP construct junction module

            for(int w=0;w!=line_tsd;++w){
                if(info_tsd[w]==seq_index){
                    
                    //string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                    string loc_tsd_0, loc_tsd_1, loc_tsd_2, loc_tsd_3, loc_tsd_4, loc_tsd_5;
                    stringstream ss_tsd_0;
                    ss_tsd_0.clear();
                    ss_tsd_0<<loc_tsd[w][0];
                    loc_tsd_0=ss_tsd_0.str();
                    stringstream ss_tsd_1;
                    ss_tsd_1.clear();
                    ss_tsd_1<<loc_tsd[w][1];
                    loc_tsd_1=ss_tsd_1.str();
                    stringstream ss_tsd_2;
                    ss_tsd_2.clear();
                    ss_tsd_2<<loc_tsd[w][2];
                    loc_tsd_2=ss_tsd_2.str();
                    stringstream ss_tsd_3;
                    ss_tsd_3.clear();
                    ss_tsd_3<<loc_tsd[w][3];
                    loc_tsd_3=ss_tsd_3.str();
                    
                    char *seq3= new char[info[i][4].length()+1];
                    strcpy(seq3, info[i][4].c_str());
                    
                    char *seq5= new char[info[i][3].length()+1];
                    strcpy(seq5, info[i][3].c_str());
                    
     //FP construct junction module
                    
                    //ofstream file22;
                    //string sys_3_kmer = WD_dir+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str()+"."+loc_tsd_0.c_str()+"."+loc_tsd_1.c_str()+"."+loc_tsd_2.c_str()+"."+loc_tsd_3.c_str()+".fake_junction_1.fasta";
                    //char *syst_3_kmer = new char[sys_3_kmer.length()+1];
                    //strcpy(syst_3_kmer, sys_3_kmer.c_str());
                    //file22.open(syst_3_kmer);
                    //file22<<">"<<seq_index<<"."<<loc_tsd_0.c_str()<<"."<<loc_tsd_1.c_str()<<"."<<loc_tsd_2.c_str()<<"."<<loc_tsd_3.c_str()<<".fake_junction_1.fasta"<<endl;
                    
                    string name_tag_1=seq_index+"."+loc_tsd_0.c_str()+"."+loc_tsd_1.c_str()+"."+loc_tsd_2.c_str()+"."+loc_tsd_3.c_str()+".fake_junction_1.fasta";
                    
                    //ofstream file222;
                    //string sys_3_kmer_2 = WD_dir+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str()+"."+loc_tsd_0.c_str()+"."+loc_tsd_1.c_str()+"."+loc_tsd_2.c_str()+"."+loc_tsd_3.c_str()+".fake_junction_2.fasta";
                    //char *syst_3_kmer_2 = new char[sys_3_kmer_2.length()+1];
                    //strcpy(syst_3_kmer_2, sys_3_kmer_2.c_str());
                    //file222.open(syst_3_kmer_2);
                    //file222<<">"<<seq_index<<"."<<loc_tsd_0.c_str()<<"."<<loc_tsd_1.c_str()<<"."<<loc_tsd_2.c_str()<<"."<<loc_tsd_3.c_str()<<".fake_junction_2.fasta"<<endl;
                    
                    string name_tag_2=seq_index+"."+loc_tsd_0.c_str()+"."+loc_tsd_1.c_str()+"."+loc_tsd_2.c_str()+"."+loc_tsd_3.c_str()+".fake_junction_2.fasta";
                    //cout<<name_tag_2<<endl;
                    /*
                    if(info[i][2]=="+"){
                        for(int n=loc_tsd[w][1]-1-J_BIN+1+J_BIN;n!=(loc_tsd[w][1]+(J_BIN*3)/2)&&n<=info[i][3].length();++n){
                            file22<<seq5[n];
                        }
                       */
                    int fix_5;
                    //fix_5=1+loc_tsd[w][4]-BIN_5+J_BIN;
                    //if(BIN_5>=(info[i][3].length()-J_BIN_mer)) fix_5=0;
                    //else if (BIN_5<(info[i][3].length()-J_BIN_mer)) fix_5=info[i][3].length()-J_BIN_mer-BIN_5-1;
                    fix_5=info[i][3].length()-J_BIN_mer-BIN_5-1;
                    
                    int fix_3;
                    //if(BIN_3>=info[i][4].length()) fix_3=0;
                    //else if(BIN_3<info[i][4].length()) fix_3=info[i][4].length()-BIN_3-1;
                    fix_3=info[i][4].length()-BIN_3-1;
                    
                    string fasta5_str="";
                    string fasta3_str="";
                    
                    if(info[i][2]=="+"){
                        stringstream fa_5,fa_3;
                        string fa_s5,fa_s3;
                        fa_3.clear();
                        fa_5.clear();
                        int patch_n=loc_tsd[w][1]-1-(J_BIN*3)/2+fix_5;
                        if(patch_n<0) {
                            //cout<<"patch_n="<<patch_n<<endl;
                            patch_n=0;
                        }
                        for(int n=patch_n;n!=(loc_tsd[w][1]+fix_5)&&n<=info[i][3].length();++n){
                            if(seq5[n]!='A'&&seq5[n]!='T'&&seq5[n]!='G'&&seq5[n]!='C'){
                                fasta5_str=fasta5_str+"N";
                            }
                            else {
                                fa_5.clear();
                                //file22<<seq5[n];
                                fa_5<<seq5[n];
                                fa_5>>fa_s5;
                                fasta5_str=fasta5_str+fa_s5;
                            }
                            
                        }
                        for(int n=loc_tsd[w][3];n!=(loc_tsd[w][3]+(J_BIN*3)/2)&&n<=info[i][4].length();++n){
                            if(seq3[n]!='A'&&seq3[n]!='T'&&seq3[n]!='G'&&seq3[n]!='C'){
                                fasta3_str=fasta3_str+"N";
                            }
                            else {
                                fa_3.clear();
                                //file222<<seq3[n];
                                fa_3<<seq3[n];
                                fa_3>>fa_s3;
                                fasta3_str=fasta3_str+fa_s3;
                            }
                        }
                    }
                    /*
                    else if(info[i][2]=="-"){
                        for(int n=loc_tsd[w][2]-1-J_BIN+J_BIN;n!=(loc_tsd[w][2]-1+(J_BIN*3)/2)&&n<=info[i][4].length();++n){
                            file222<<seq3[n];
                        }
                        */
                    
                    else if(info[i][2]=="-"){
                        stringstream fa_5,fa_3;
                        string fa_s5,fa_s3;
                        fa_3.clear();
                        fa_5.clear();
                        int patch_n_2=loc_tsd[w][2]-1-(J_BIN*3)/2+fix_3-1;
                        if(patch_n_2<0) {
                            //cout<<"patch_n_2="<<patch_n_2<<endl;
                            patch_n_2=0;
                        }
                        for(int n=patch_n_2;n!=(loc_tsd[w][2]+fix_3-1)&&n<=info[i][4].length();++n){
                            if(seq3[n]!='A'&&seq3[n]!='T'&&seq3[n]!='G'&&seq3[n]!='C'){
                                fasta3_str=fasta3_str+"N";
                            }
                            else {
                                fa_3.clear();
                                //file222<<seq3[n];
                                fa_3<<seq3[n];
                                fa_3>>fa_s3;
                                fasta3_str=fasta3_str+fa_s3;
                            }
                        }
                        for(int n=loc_tsd[w][0]-1+J_BIN_mer;n!=(loc_tsd[w][0]+(J_BIN*3)/2-1+J_BIN_mer)&&n<=info[i][3].length();++n){
                            if(seq5[n]!='A'&&seq5[n]!='T'&&seq5[n]!='G'&&seq5[n]!='C'){
                                fasta5_str=fasta5_str+"N";
                            }
                            else {
                                fa_5.clear();
                                //file22<<seq5[n];
                                fa_5<<seq5[n];
                                fa_5>>fa_s5;
                                fasta5_str=fasta5_str+fa_s5;
                            }
                        }
                    }
                    if(fasta3_str==""||fasta3_str==" "){
                        fasta3_str=fasta3_str+"N";
                    }
                    if(fasta5_str==""||fasta5_str==" "){
                        fasta5_str=fasta5_str+"N";
                    }
                    //file22<<endl;
                    //file222<<endl;
                    int fasta5_str_len=int(fasta5_str.length()*60/75);
                    int fasta3_str_len=int(fasta3_str.length()*60/75);
                    string fasta5_str_length=to_string(fasta5_str_len);
                    string fasta3_str_length=to_string(fasta3_str_len);
                    
                    string sys_3_blast = "bash -c \"blastn -evalue 0.05 -task blastn -query <(printf "+(fasta5_str)+") -subject "+ref_junc_file+" -dust no -outfmt \\\"7 std\\\" |grep -v \\\"#\\\" | awk '{if(\\\$3>=80&&\\\$4>="+fasta5_str_length+"&&(\\\$10-\\\$9)>0) print \\\"1\\\"}' |wc -l \"";
                    //> "+WD_dir+"read_result_junc_fake.txt";
                    //cout<<sys_3_blast<<endl;
                    
                    char *syst_3_blast = new char[sys_3_blast.length()+1];
                    strcpy(syst_3_blast, sys_3_blast.c_str());
                    vector<string> pp_3_int;
                    pp_3_int.clear();
                    FILE *pp_3 =popen(syst_3_blast,"r");
                    char *tmp=new char[1024];
                    while (fgets(tmp, sizeof(tmp), pp_3) != NULL) {
                        if (tmp[strlen(tmp) - 1] == '\n') {
                            tmp[strlen(tmp) - 1] = '\0';
                        }
                        //cout<<tmp<<endl;
                        pp_3_int.push_back(tmp);
                        
                    }
                    pclose(pp_3);
                    //loc_tsd_fp[w][0]=pp_3_int;
                    loc_tsd_fp[w][0]=accumulate(pp_3_int.begin(),pp_3_int.end(),loc_tsd_fp[w][0]);
                    
                    string sys_3_blast_2 = "bash -c \"blastn -evalue 0.05 -task blastn -query <(printf "+(fasta3_str)+") -subject "+ref_junc_file+" -dust no -outfmt \\\"7 std\\\" |grep -v \\\"#\\\" | awk '{if(\\\$3>=80&&\\\$4>="+fasta3_str_length+"&&(\\\$10-\\\$9)>0) print \\\"1\\\"}' |wc -l \"";
                    //> "+WD_dir+"read_result_junc_fake.txt";
                    //cout<<sys_3_blast_2<<endl;
                    
                    char *syst_3_blast_2 = new char[sys_3_blast_2.length()+1];
                    strcpy(syst_3_blast_2, sys_3_blast_2.c_str());
                    vector<string> pp_3_int_2;
                    pp_3_int_2.clear();
                    FILE *pp_3_2 =popen(syst_3_blast_2,"r");
                    char *tmp_2=new char[1024];
                    while (fgets(tmp_2, sizeof(tmp_2), pp_3_2) != NULL) {
                        if (tmp_2[strlen(tmp_2) - 1] == '\n') {
                            tmp_2[strlen(tmp_2) - 1] = '\0';
                        }
                        //cout<<tmp_2<<endl;
                        pp_3_int_2.push_back(tmp_2);
                    }
                    pclose(pp_3_2);
                    loc_tsd_fp[w][1]=accumulate(pp_3_int_2.begin(),pp_3_int_2.end(),loc_tsd_fp[w][1]);
                    //loc_tsd_fp[w][1]=pp_3_int_2;
                    
                    if(loc_tsd_fp[w][0]==""){
                        loc_tsd_fp[w][0]="-1";
                    }
                    if(loc_tsd_fp[w][1]==""){
                        loc_tsd_fp[w][1]="-1";
                    }
                    
                    delete [] tmp;
                    delete [] tmp_2;
                    /*
                    string sys_3_blast = "blastn -evalue 0.05 -task blastn -query "+sys_3_kmer+" -subject "+ref_junc_file+" -dust no -outfmt \"7 std\" |grep -v \"#\" | awk '{if($3>=80&&$4>=60&&($10-$9)>0) print \"1\"}' |wc -l > "+WD_dir+"read_result_junc_fake.txt";
                    char *syst_3_blast = new char[sys_3_blast.length()+1];
                    strcpy(syst_3_blast, sys_3_blast.c_str());
                    system(syst_3_blast);
                    
                    string sys_3_bl = WD_dir+"read_result_junc_fake.txt";
                    char *syst_3_bl = new char[sys_3_bl.length()+1];
                    strcpy(syst_3_bl, sys_3_bl.c_str());
                    ifstream file23;
                    file23.open(syst_3_bl);

                    file23>>loc_tsd[w][7];
                    
                    string sys_3_blast_2 = "blastn -evalue 0.05 -task blastn -query "+sys_3_kmer_2+" -subject "+ref_junc_file+" -dust no -outfmt \"7 std\" |grep -v \"#\" | awk '{if($3>=80&&$4>=60&&($10-$9)>0) print \"1\"}' |wc -l > "+WD_dir+"read_result_junc_fake_2.txt";
                    char *syst_3_blast_2 = new char[sys_3_blast_2.length()+1];
                    strcpy(syst_3_blast_2, sys_3_blast_2.c_str());
                    system(syst_3_blast_2);
                    
                    string sys_3_bl_2 = WD_dir+"read_result_junc_fake_2.txt";
                    char *syst_3_bl_2 = new char[sys_3_bl_2.length()+1];
                    strcpy(syst_3_bl_2, sys_3_bl_2.c_str());
                    ifstream file232;
                    file232.open(syst_3_bl_2);
                    
                    file232>>loc_tsd[w][8];
                    */
    //5' 26mer identify
                    /*
                    ofstream file42;
                    string sys_5_kmer = WD_dir+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str()+"."+loc_tsd_0.c_str()+"."+loc_tsd_1.c_str()+"."+loc_tsd_2.c_str()+"."+loc_tsd_3.c_str()+".5kmer.fasta";
                    char *syst_5_kmer = new char[sys_5_kmer.length()+1];
                    strcpy(syst_5_kmer, sys_5_kmer.c_str());
                    file42.open(syst_5_kmer);
                    
                    file42<<">"<<seq_index<<"."<<loc_tsd_0.c_str()<<"."<<loc_tsd_1.c_str()<<"."<<loc_tsd_2.c_str()<<"."<<loc_tsd_3.c_str()<<".5kmer.fasta"<<endl;
                    */
                    if(info[i][2]=="+"){
                        stringstream kmer_ss;
                        string kmer_s;
                        kmer_ss.clear();
                        for(int n=loc_tsd[w][1]-1-J_BIN_mer+1+J_BIN;n!=(loc_tsd[w][1]+J_BIN+J_BIN_mer)&&n<=info[i][3].length();++n){
                            //file42<<seq5[n];
                            if(seq5[n]!='A'&&seq5[n]!='T'&&seq5[n]!='G'&&seq5[n]!='C'){
                                kmer_tsd[w]=kmer_tsd[w]+"N";
                            }
                            else {
                                kmer_ss.clear();
                                kmer_ss<<seq5[n];
                                kmer_ss>>kmer_s;
                                kmer_tsd[w]=kmer_tsd[w]+kmer_s;
                            }
                        }
                        if(kmer_tsd[w]==""){
                            kmer_tsd[w]="N";
                        }
                        
                    }
                    else if(info[i][2]=="-"){
                        
                        stringstream kmer_ss;
                        string kmer_s;
                        kmer_ss.clear();
                        for(int n=loc_tsd[w][0]+(BIN_5+1-loc_tsd[w][4])-1+J_BIN_mer-J_BIN_mer;n!=(loc_tsd[w][0]+(BIN_5+1-loc_tsd[w][4])+J_BIN_mer+J_BIN_mer-1)&&n<=info[i][3].length();++n){
                            //file42<<seq5[n];
                            if(seq5[n]!='A'&&seq5[n]!='T'&&seq5[n]!='G'&&seq5[n]!='C'){
                                kmer_tsd[w]=kmer_tsd[w]+"N";
                            }
                            else {
                                kmer_ss.clear();
                                kmer_ss<<seq5[n];
                                kmer_ss>>kmer_s;
                                kmer_tsd[w]=kmer_tsd[w]+kmer_s;
                            }
                        }
                        if(kmer_tsd[w]==""){
                            kmer_tsd[w]="N";
                        }
                    }
                    //file42<<endl;
                    
                    //No need comparison blastn for now
                    /*
                    string sys_5_blast = "blastn -evalue 0.05 -task blastn -query "+sys_5_kmer+" -subject "+ref_file+" -dust no -outfmt \"7 std\" |grep -v \"#\" | awk '{if($3>=85&&$4>=22&&($10-$9)>0) print \"1\"}' |wc -l > "+WD_dir+"read_result_junc_5kmers.txt";
                    char *syst_5_blast = new char[sys_5_blast.length()+1];
                    strcpy(syst_5_blast, sys_5_blast.c_str());
                    system(syst_5_blast);
                    
                    string sys_5_bl = WD_dir+"read_result_junc_5kmers.txt";
                    char *syst_5_bl = new char[sys_5_bl.length()+1];
                    strcpy(syst_5_bl, sys_5_bl.c_str());
                    ifstream file43;
                    file43.open(syst_5_bl);
                    
                    file43>>loc_tsd[w][6];
                    */
                    loc_tsd[w][6]=0;

                    delete [] seq3;
                    delete [] seq5;
                    //file43.close();
                    //file43.clear();
                    //file22.close();
                    //file23.close();
                    //file22.clear();
                    //file23.clear();
                    //file222.close();
                    //file232.close();
                    //file222.clear();
                    //file232.clear();
                    //file42.close();
                    //file42.clear();
                    
                }
            }
            //file20.clear();
            //file20.close();
            file21.close();
            file21.clear();
        }
        
        for(int w=0;w!=line_tsd;++w){
            file11<<info_tsd[w]<<'\t'<<loc_tsd[w][0]<<'\t'<<loc_tsd[w][1]<<'\t'<<loc_tsd[w][2]<<'\t'<<loc_tsd[w][3]<<'\t'<<loc_tsd[w][4]<<'\t'<<loc_tsd[w][5]<<'\t'<<loc_tsd[w][6]<<'\t'<<loc_tsd_fp[w][0]<<'\t'<<loc_tsd_fp[w][1]<<'\t'<<kmer_tsd[w]<<endl;
        }
    }
    
    else if(tsd_index==0){
        for(int w=0;w!=line_tsd;++w){
            file11<<info_tsd[w]<<'\t'<<loc_tsd[w][0]<<'\t'<<loc_tsd[w][1]<<'\t'<<loc_tsd[w][2]<<'\t'<<loc_tsd[w][3]<<'\t'<<loc_tsd[w][4]<<'\t'<<loc_tsd[w][5]<<'\t'<<"0"<<'\t'<<"1"<<'\t'<<"1"<<'\t'<<"N"<<endl;
        }
    }
    
    for(int i=0;i!=line;++i){
        delete [] info[i];
        delete [] loc[i];
        delete [] loc_TP[i];
    }
    delete [] info;
    delete [] loc;
    
    for(int i=0;i!=line_tsd;++i){
        //delete [] info_tsd[i];
        delete [] loc_tsd[i];
        delete [] loc_tsd_fp[i];
    }
    delete [] info_tsd;
    delete [] loc_tsd;
    delete [] loc_TP;
    delete [] loc_tsd_fp;
    return 0;
    
}
