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
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
using namespace std;

int fp_ex(string WD_dir, string fasta, string chr){

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
    for(int i=0;!file2.eof();i++){
        getline(file2,input);
        line=i;
    }
    
    file2.close();
    file2.clear();
    file2.open(syst_junc);
    
    string **info;
    info=new string *[line];
    for(int i=0;i!=line;i++) info[i]=new string[5];
    
    int **loc;
    loc=new int*[line];
    for(int i=0;i!=line;i++) loc[i]=new int[7];
    
    for(int i=0;i!=line;i++){
        file2>>info[i][0];  //probe name
        file2>>loc[i][0];   //L_loc
        file2>>loc[i][1];   //L_loc
        file2>>loc[i][2];   //R_loc
        file2>>loc[i][3];   //R_loc
        file2>>loc[i][4];   //G_loc
        file2>>loc[i][5];   //G_loc
        file2>>info[i][1];  //chr
        file2>>info[i][2];  //orientation
        file2>>info[i][3]; //junc 5' seq
        file2>>info[i][4]; //junc 3' seq
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
    
    int line_tsd;
    for(int i=0;!file1.eof();i++){
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
    for(int i=0;i!=line_tsd;i++) loc_tsd[i]=new int[8];
    
    for(int i=0;i!=line_tsd;i++){
        file1>>info_tsd[i];
        file1>>loc_tsd[i][0];
        file1>>loc_tsd[i][1];
        file1>>loc_tsd[i][2];
        file1>>loc_tsd[i][3];
        file1>>loc_tsd[i][4];
        file1>>loc_tsd[i][5];
        loc_tsd[i][6]=-1;
        loc_tsd[i][7]=-1;
        //loc_tsd[i][8]=0;
    }
    
    ofstream file31;
    string sys_5_bl = WD_dir+"read_result_junc_5kmers.txt";
    char *syst_5_bl = new char[sys_5_bl.length()+1];
    strcpy(syst_5_bl, sys_5_bl.c_str());
    file31.open(syst_5_bl,ios::trunc);
    file31.clear();
    file31.close();
    //FP_ex module
    for(int i=0;i!=line;i++){
        
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
        
        string seq_index;
        
        seq_index=info[i][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2];
        
        //ref pull out
        int start, end;
        start=loc[i][4]-6000;
        if(start<=0){
            start=1;
        }
        end=loc[i][5]+6000;
        
        string s_start, s_end;
        stringstream ss_s;
        ss_s.clear();
        ss_s<<start;
        s_start=ss_s.str();
        stringstream ss_e;
        ss_e.clear();
        ss_e<<end;
        s_end=ss_e.str();
        
        ofstream file20;
        string ref_file;
        ref_file=WD_dir+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2]+".ref.fasta";
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
        
        //5' 26mer identify
        ofstream file21;
        string sys_5_kmer = WD_dir+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2]+".5kmer.fasta";
        char *syst_5_kmer = new char[sys_5_kmer.length()+1];
        strcpy(syst_5_kmer, sys_5_kmer.c_str());
        file21.open(syst_5_kmer);
        
        file21<<">"<<seq_index<<endl;
        file21<<info[i][3]<<endl;
        
        string sys_5_blast = "blastn -evalue 0.05 -task blastn -query "+sys_5_kmer+" -subject "+ref_file+"  -outfmt \"7 std\" |grep -v \"#\" | awk '{if($3>=90&&($10-$9)>0) print \"1\"}' |wc -l >> "+WD_dir+"read_result_junc_5kmers.txt";
        char *syst_5_blast = new char[sys_5_blast.length()+1];
        strcpy(syst_5_blast, sys_5_blast.c_str());
        system(syst_5_blast);
        //cout<<sys_5_blast<<endl;
        //cout<<seq_index<<endl;
        
        //3' 26mer identify based on TSD
        for(int w=0;w!=line_tsd;w++){
            if(info_tsd[w]==seq_index){
                
                string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                stringstream ss_0;
                ss_0.clear();
                ss_0<<loc_tsd[w][0];
                loc_0=ss_0.str();
                stringstream ss_1;
                ss_1.clear();
                ss_1<<loc_tsd[w][1];
                loc_1=ss_1.str();
                stringstream ss_2;
                ss_2.clear();
                ss_2<<loc_tsd[w][2];
                loc_2=ss_2.str();
                stringstream ss_3;
                ss_3.clear();
                ss_3<<loc_tsd[w][3];
                loc_3=ss_3.str();
                
                ofstream file22;
                string sys_3_kmer = WD_dir+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+".3kmer.fasta";
                char *syst_3_kmer = new char[sys_3_kmer.length()+1];
                strcpy(syst_3_kmer, sys_3_kmer.c_str());
                file22.open(syst_3_kmer);
                file22<<">"<<seq_index<<"."<<loc_0.c_str()<<"."<<loc_1.c_str()<<"."<<loc_2.c_str()<<"."<<loc_3.c_str()<<".3kmer.fasta"<<endl;
                
                char *seq3= new char[info[i][4].length()+1];
                strcpy(seq3, info[i][4].c_str());
                
                if(info[i][2]=="+"){
                    for(int n=loc_tsd[w][2];n!=(loc_tsd[w][2]+26)&&n<info[i][4].length();n++){
                        file22<<seq3[n];
                    }
                }
                else if(info[i][2]=="-"){
                    for(int n=loc_tsd[w][3]-13;n!=(loc_tsd[w][3]+13)&&n<info[i][4].length();n++){
                        file22<<seq3[n];
                    }
                }
                file22<<endl;
                
                string sys_3_blast = "blastn -evalue 0.05 -task blastn -query "+sys_3_kmer+" -subject "+ref_file+"  -outfmt \"7 std\" |grep -v \"#\" | awk '{if($3>=90&&($10-$9)>0) print \"1\"}' |wc -l > "+WD_dir+"read_result_junc_3kmers.txt";
                char *syst_3_blast = new char[sys_3_blast.length()+1];
                strcpy(syst_3_blast, sys_3_blast.c_str());
                system(syst_3_blast);
                
                //cout<<sys_3_blast<<endl;
                
                
                string sys_3_bl = WD_dir+"read_result_junc_3kmers.txt";
                char *syst_3_bl = new char[sys_3_bl.length()+1];
                strcpy(syst_3_bl, sys_3_bl.c_str());
                ifstream file23;
                file23.open(syst_3_bl);
                
                //input bl_3;
                file23>>loc_tsd[w][7];
                
                
                file22.close();
                file23.close();
                file22.clear();
                file23.clear();
            }
        }
        
        file20.close();
        file21.close();
        file20.clear();
        file21.clear();
    }
    
    ifstream file32;
    file32.open(syst_5_bl);
    for(int i=0;i!=line;i++){
        
        int int_5_bl;
        file32>>int_5_bl;
        
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
        
        string seq_index;
        
        seq_index=info[i][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2];
        for(int w=0;w!=line_tsd;w++){
            if(info_tsd[w]==seq_index){
                loc_tsd[w][6]=int_5_bl;
            }
        }
    }
    
    for(int w=0;w!=line_tsd;w++){
        file11<<info_tsd[w]<<'\t'<<loc_tsd[w][0]<<'\t'<<loc_tsd[w][1]<<'\t'<<loc_tsd[w][2]<<'\t'<<loc_tsd[w][3]<<'\t'<<loc_tsd[w][4]<<'\t'<<loc_tsd[w][5]<<'\t'<<loc_tsd[w][6]<<'\t'<<loc_tsd[w][7]<<endl;
    }
    
    for(int i=0;i!=line;i++){
        delete [] info[i];
        delete [] loc[i];
    }
    delete [] info;
    delete [] loc;
    
    for(int i=0;i!=line_tsd;i++){
        //delete [] info_tsd[i];
        delete [] loc_tsd[i];
    }
    delete [] info_tsd;
    delete [] loc_tsd;
    
    return 0;
    
}