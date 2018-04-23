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

int tsd_module(string WD_dir){
//for test
//int main(){
    
    //string WD_dir="/Users/zhouweichen/Documents/workspace/16.10.08.Pacbio_analysis/17.11.24.PALMER/v1.1.1/scp/IO_test/";
    
    ifstream file1;
    ifstream file2;
    
    string sys_readres = WD_dir + "read_result.txt";
    char * syst_readres = new char [sys_readres.length()+1];
    strcpy (syst_readres, sys_readres.c_str());
    
    file1.open(syst_readres);

    string sys_reg = WD_dir+"region.sam";
    char *syst_reg = new char[sys_reg.length()+1];
    strcpy(syst_reg, sys_reg.c_str());
    file2.open(syst_reg);
    
    int line_read;
    int line_region;
    
    string input;
    for(int i=0;!file1.eof();i++){
        getline(file1,input);
        line_read=i;
    }
    for(int i=0;!file2.eof();i++){
        getline(file2,input);
        line_region=i;
    }
    file1.close();
    file1.clear();
    file2.close();
    file2.clear();
    file1.open(syst_readres);
    file2.open(syst_reg);
    
    string **SEQ;
    SEQ=new string*[line_region];
    for(int i=0;i!=line_region;i++) SEQ[i]=new string[2];
    string *name;
    name=new string[line_read];
    
    string **info;
    info=new string*[line_read];
    for(int i=0;i!=line_read;i++) info[i]=new string[2];
    
    int **loc;
    loc=new int*[line_read];
    for(int i=0;i!=line_read;i++) loc[i]=new int[6];
    
    for(int i=0;i!=line_read;i++){
        file1>>name[i];
        file1>>loc[i][0];
        file1>>loc[i][1];
        file1>>loc[i][2];
        file1>>loc[i][3];
        file1>>loc[i][4];
        file1>>loc[i][5];
        file1>>info[i][0];
        file1>>info[i][1];
    }
    for(int i=0;i!=line_region;i++){
        file2>>SEQ[i][0];
        file2>>input;
        file2>>input;
        file2>>input;
        file2>>input;
        file2>>input;
        file2>>input;
        file2>>input;
        file2>>input;
        file2>>SEQ[i][1];
        getline(file2,input);
    }
    ofstream file3;
    ifstream file4;
    
    ofstream file5;
    ifstream file6;
    
    ofstream file8;
    ifstream file9;

    ofstream file7;
    ofstream file11;
    
    ofstream file12;
    ifstream file13;
    
    string sys_tsd = WD_dir+"read_result_TSD.txt";
    char *syst_tsd = new char[sys_tsd.length()+1];
    strcpy(syst_tsd, sys_tsd.c_str());
    file7.open(syst_tsd);
    
    //string sys_seq = WD_dir+"SEQ";
    //char *syst_seq = new char[sys_seq.length()+1];
    //strcpy(syst_seq, sys_seq.c_str());
    
    //string sys_readtag = WD_dir+"readtag";
    //char *syst_readtag = new char[sys_readtag.length()+1];
    //strcpy(syst_readtag, sys_readtag.c_str());
    
    string sys_tsd_info = WD_dir+"read_result_info.txt";
    char *syst_tsd_info = new char[sys_tsd_info.length()+1];
    strcpy(syst_tsd_info, sys_tsd_info.c_str());
    file11.open(syst_tsd_info);
    
    int J_BIN=13;
    int BIN=6;
    int BIN_3=3000;

    for(int i=0;i!=line_read;i++){
        int js1, js2, je1, je2;
        js1=loc[i][2]-J_BIN;
        if(js1<=0) js1=1;
        je1=loc[i][2]+J_BIN-1;
        js2=loc[i][3]-J_BIN+1;
        je2=loc[i][3]+J_BIN;
        
        int l1,l2;
        l1=je1-js1+1;
        l2=je2-js2+1;
        char *right_j;
        right_j=new char[l2];
        char *left_j;
        left_j=new char[l1];
        for(int j=0;j!=l1;j++) left_j[j]='N';
        for(int j=0;j!=l2;j++) right_j[j]='N';
        
        int s1,s2,e1,e2;
        s1=loc[i][2]-45;
        if(s1<=0) s1=1;
        e1=loc[i][2]+5;
        s2=loc[i][3]-5;
        e2=loc[i][3]+45;

        int s_3_1,s_3_2,e_3_1,e_3_2;
        s_3_1=loc[i][2]-BIN_3;
        if(s_3_1<=0) s_3_1=1;
        e_3_1=loc[i][2];
        s_3_2=loc[i][3];
        e_3_2=loc[i][3]+BIN_3;
        
        int l3,l4,l5,l6;
        l3=e1-s1+1;
        l4=e2-s2+1;
        l5=e_3_1-s_3_1+1;
        l6=e_3_2-s_3_2+1;
        char *left;
        char *left_3;
        char *right;
        char *right_3;
        left=new char[l3];
        right=new char[l4];
        left_3=new char[l5];
        right_3=new char[l6];
        for(int j=0;j!=l3;j++) left[j]='N';
        for(int j=0;j!=l4;j++) right[j]='N';
        for(int j=0;j!=l5;j++) left_3[j]='N';
        for(int j=0;j!=l6;j++) right_3[j]='N';
        
        for(int j=0;j!=line_region;j++){
            if(SEQ[j][0]==name[i]){
                //file3.open(syst_seq);
                //file3<<SEQ[j][1]<<"E"<<endl;;
                //file4.open(syst_seq);
                
                string seq_1;
                seq_1=SEQ[j][1]+"E";
                stringstream ss_seq;
                ss_seq.clear();
                ss_seq.str(seq_1);
                
                char seq;
                ss_seq>>seq;
                
                int n=1;
                int le1,le2,le3,le4,le5,le6;
                le1=le2=le3=le4=le5=le6=0;
                for(;seq!='E';){
                    if(n>=js1&&n<=je1) {left_j[le1]=seq;le1=le1+1;}
                    if(n>=js2&&n<=je2) {right_j[le2]=seq;le2=le2+1;}
                    
                    if(n>=s1&&n<=e1) {left[le3]=seq;le3=le3+1;}
                    if(n>=s2&&n<=e2) {right[le4]=seq;le4=le4+1;}
                    if(n>=s_3_1&&n<=e_3_1) {left_3[le5]=seq;le5=le5+1;}
                    if(n>=s_3_2&&n<=e_3_2) {right_3[le6]=seq;le6=le6+1;}
                    
                    //if(n>=loc[i][0]&&n<=loc[i][1]) {ins[le5]=seq;le5=le5+1;}
                    n++;
                    ss_seq>>seq;
                }
                //file4.close();
                //file4.clear();
                //file3.close();
                //file3.clear();
            }
        }
        
        file11<<name[i]<<'\t'<<loc[i][0]<<'\t'<<loc[i][1]<<'\t'<<loc[i][2]<<'\t'<<loc[i][3]<<'\t'<<loc[i][4]<<'\t'<<loc[i][5]<<'\t'<<info[i][0]<<'\t'<<info[i][1]<<'\t';
        if(info[i][1]=="+"){
            for(int w=0;w!=l1;w++){
                file11<<left_j[w];
            }
            file11<<endl;
        }
        if(info[i][1]=="-"){
            for(int w=0;w!=l2;w++){
                file11<<right_j[w];
            }
            file11<<endl;
        }
        
        file7<<name[i]<<'\t'<<loc[i][0]<<'\t'<<loc[i][1]<<'\t'<<loc[i][2]<<'\t'<<loc[i][3]<<'\t'<<loc[i][4]<<'\t'<<loc[i][5]<<'\t'<<info[i][0]<<'\t'<<info[i][1]<<'\t';
        
        
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
        
        seq_index="_";
        seq_index=seq_index+loc_0.c_str()+"_"+loc_1.c_str()+"_"+loc_2.c_str()+"_"+loc_3.c_str()+"_"+loc_4.c_str()+"_"+loc_5.c_str()+"_"+info[i][0]+"_"+info[i][1];
        
        
        //file12.open(syst_readtag);
        //file12<<"_"<<loc[i][0]<<"_"<<loc[i][1]<<"_"<<loc[i][2]<<"_"<<loc[i][3]<<"_"<<loc[i][4]<<"_"<<loc[i][5]<<"_"<<info[i][0]<<"_"<<info[i][1]<<endl;
        
        //string seq_index;
        //file13.open(syst_readtag);
        //getline(file13,seq_index);
        
        //file12.close();
        //file12.clear();
        //file13.close();
        //file13.clear();
        
        string sys_seq_5 = WD_dir+seq_index+".5.fasta";
        char *syst_seq_5 = new char[sys_seq_5.length()+1];
        strcpy(syst_seq_5, sys_seq_5.c_str());
        file5.open(syst_seq_5);
        
        string sys_seq_3 = WD_dir+seq_index+".3.fasta";
        char *syst_seq_3 = new char[sys_seq_3.length()+1];
        strcpy(syst_seq_3, sys_seq_3.c_str());
        file8.open(syst_seq_3);
        
        file5<<">"<<seq_index<<endl;
        file8<<">"<<seq_index<<endl;
    
        if(info[i][1]=="+"){
            
            for(int w=0;w!=l3;w++){
                file7<<left[w];
                file5<<left[w];
            }
            file7<<'\t';
            
            for(int w=0;w!=l6;w++){
                file8<<right_3[w];
                file7<<right_3[w];
            }
            file7<<endl;
            
            file8.close();
            file5.close();
            file8.clear();
            file5.clear();
            
        }
        if(info[i][1]=="-"){
            
            for(int w=0;w!=l4;w++){
                file7<<right[w];
                file5<<right[w];
            }
            file7<<'\t';
            
            for(int w=0;w!=l5;w++){
                file8<<left_3[w];
                file7<<left_3[w];
            }
            file7<<endl;
            
            file8.close();
            file5.close();
            file8.clear();
            file5.clear();
        
        }

        
//Possible TSD finding & filtering
        int len_5, len_3;
        stringstream ss,ss1;
        string s_len,s_len1;
        
        if(info[i][1]=="+"){
            ss<<51-e1+s1-1;
            ss>>s_len;
            ss1<<e1-s1+1;
            ss1>>s_len1;
            string sys_blastn = "blastn -task blastn -query "+sys_seq_5+" -subject "+sys_seq_3+" -word_size 6 -dust no -outfmt \"7 std\" |grep -v \"#\" | awk '{if($3>=80&&$4>=6&&($10-$9)>0) print \""+name[i]+seq_index+"\","+(s_len)+"+$7,"+(s_len)+"+$8,$9,$10,"+s_len1+",\"3001\"}' >> "+WD_dir+"TSD_blastn.txt";
            char *syst_blastn = new char[sys_blastn.length()+1];
            strcpy(syst_blastn, sys_blastn.c_str());
            system(syst_blastn);
        }
        else if(info[i][1]=="-"){
            ss<<3001-e_3_1+s_3_1-1;
            ss>>s_len;
            ss1<<e_3_1-s_3_1+1;
            ss1>>s_len1;
            string sys_blastn = "blastn -task blastn -query "+sys_seq_5+" -subject "+sys_seq_3+" -word_size 6 -dust no -outfmt \"7 std\" |grep -v \"#\" | awk '{if($3>=85&&$4>=6&&($10-$9)>0) print \""+name[i]+seq_index+"\",$7,$8,"+(s_len)+"+$9,"+(s_len)+"+$10,\"51\","+s_len1+"}' >> "+WD_dir+"TSD_blastn.txt";
            char *syst_blastn = new char[sys_blastn.length()+1];
            strcpy(syst_blastn, sys_blastn.c_str());
            system(syst_blastn);
        }
        
        
    }
    return 0;
}
