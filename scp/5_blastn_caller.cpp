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

int BlastnCaller(string WD_dir, string chr_fix, string t){
//test
//int main(){
    //string WD_dir="/Users/zhouweichen/Documents/workspace/16.10.08.Pacbio_analysis/17.11.24.PALMER/v1.1.1/scp/IO_test/";
    //string chr_fix="chr1";
    //string t="LINE";
    
    
    //cout<<"Blastn Caller Step is now running."<<endl;
    string sys_blastncaller;
    
    sys_blastncaller = "cat "+WD_dir+"blastn.txt |grep -v \"#\" > "+WD_dir+"blastn_refine.txt";
    
    
    char *syst_blastncaller = new char[sys_blastncaller.length()+1];
    strcpy(syst_blastncaller, sys_blastncaller.c_str());
    
    system(syst_blastncaller);

//Blastn Caller
    string sys_cigar = WD_dir+"cigar.2";
    char *syst_cigar =  new char[sys_cigar.length()+1];
    strcpy(syst_cigar, sys_cigar.c_str());
    
    ifstream file1;
    file1.open(syst_cigar);
    
    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'cigar.2'"<< endl;
        //exit(1);
        return 0;
    }
    
    string sys_blastnrefine = WD_dir+"blastn_refine.txt";
    char *syst_blastnrefine =  new char[sys_blastnrefine.length()+1];
    strcpy(syst_blastnrefine, sys_blastnrefine.c_str());
    
    ifstream file2;
    file2.open(syst_blastnrefine);
    
    if (!file2.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'blastn_refine.txt'"<< endl;
        //exit(1);
        return 0;
    }

    string sys_selecinfo = WD_dir+"selected.reads.info";
    char *syst_selecinfo =  new char[sys_selecinfo.length()+1];
    strcpy(syst_selecinfo, sys_selecinfo.c_str());
    
    ifstream file3;
    file3.open(syst_selecinfo);
    
    if (!file3.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'selected.reads.info'"<< endl;
        //exit(1);
        return 0;
    }

    //string sys_cigar3 = WD_dir+"cigar.3";
    //char *syst_cigar3 =  new char[sys_cigar3.length()+1];
    //strcpy(syst_cigar3, sys_cigar3.c_str());
    
    //ofstream file4;
    //file4.open(syst_cigar3);
    
    string sys_readresult = WD_dir+"read_result.txt";
    char *syst_readresult =  new char[sys_readresult.length()+1];
    strcpy(syst_readresult, sys_readresult.c_str());
    
    ofstream file5;
    file5.open(syst_readresult);
    
    ifstream file6;
    
    int line;
    string input;
    for(int i=0;!file3.eof();i++){
        getline(file3,input);
        line=i;
    }
    int blast;
    for(int i=0;!file2.eof();i++){
        getline(file2,input);
        blast=i;
    }
    
    file3.close();
    file3.clear();
    file3.open(syst_selecinfo);
    
    file2.close();
    file2.clear();
    file2.open(syst_blastnrefine);
    
    int **bla;
    bla=new int*[blast];
    for(int i=0;i!=blast;i++) bla[i]=new int[7];
    
    string *bla_name;
    bla_name=new string[blast];
    
    string *orient;
    orient = new string[blast];
    
    string **read;
    read=new string*[line];
    for(int i=0;i!=line;i++) read[i]=new string[3];
    int *read_loc;
    read_loc= new int[line];
    
    for(int i=0;i!=blast;i++){
        file2>>input;
        file2>>bla_name[i];     //read name
        file2>>input;
        file2>>bla[i][0];       //L1 s
        file2>>bla[i][1];       //L1 e
        file2>>bla[i][2];       //read s
        file2>>bla[i][3];       //read e
        orient[i]="+";          //orientation
        bla[i][4]=0;            //insert site s
        bla[i][5]=0;            //insert site e
        bla[i][6]=0;            //flag_bn
    }
    for(int i=0;i!=line;i++){
        file1>>read[i][2];      //cigar
        file3>>read[i][0];      //read name
        file3>>read[i][1];      //chr
        file3>>read_loc[i];     //pos
    }
    for(int i=0;i!=blast;i++){
        int buff;
        int insert_s=0;
        int insert_e=0;
        int L1_s=bla[i][0];
        int L1_e=bla[i][1];
        if(L1_s>L1_e){
            buff=L1_s;
            L1_s=L1_e;
            L1_e=buff;
        }
        int r_s=bla[i][2];
        int r_e=bla[i][3];
        if(r_s>r_e){
            buff=r_s;
            r_s=r_e;
            r_e=buff;
            orient[i]="-";
        }
        int M=0;
        int flag_r=0;
        
        for(int j=0;j!=line&&flag_r==0;j++){
            if(bla_name[i]==read[j][0]){
                flag_r=1;
                
                //file4<<read[j][2]<<endl;
                //file6.open(syst_cigar3);
                char cig;
                int number;
                int bit;        //tag in ref
                bit=read_loc[j]; //from the start pos
                insert_s=bit;
                insert_e=bit;
                int k=0;        //tag in read
                
                stringstream ss_cigar3;
                ss_cigar3.clear();
                ss_cigar3.str(read[j][2]);
                
                ss_cigar3>>number;
                ss_cigar3>>cig;
                
                for(;!(cig=='E'&&number==0);){
                    if(cig=='M'||cig=='X'||cig=='='){
                        k=k+number;
                        //M=M+number;
                        if(k<r_s){
                            insert_s=bit+number;
                            insert_e=bit+number;
                        }
                        else if(k>=r_s&&k<=r_e&&(k-number)<=r_s){
                            insert_s=bit+r_s-k+number;
                            insert_e=insert_s;
                            M=M+k-r_s;
                        }
                        else if(k>=r_s&&k<=r_e&&(k-number)>r_s){
                            //insert_s=bit;
                            insert_e=insert_s+number;
                            M=M+number;
                        }
                        else if(k>r_e&&(k-number)<=r_e&&(k-number)>=r_s){
                            insert_e=bit+r_e-k+number;
                            M=M+r_e-k+number;
                        }
                        else if(k>r_e&&(k-number)>=r_e){
                            break;
                        }
                        
                        bit=bit+number;
                    }
                    else if(cig=='S'||cig=='I'){
                        k=k+number;
                    }
                    else if(cig=='D'){
                        bit=bit+number;
                    }
                    ss_cigar3>>number;
                    ss_cigar3>>cig;
                }
                //file6.close();
                //file6.clear();
                //file4.close();
                //file4.clear();
                //file4.open(syst_cigar3);
            }
        }
        if(M>=int(0.5*(r_e-r_s))){ 
            bla[i][6]=1;
        }
        
        bla[i][0]=L1_s;
        bla[i][1]=L1_e;
        bla[i][2]=r_s;
        bla[i][3]=r_e;
        bla[i][4]=insert_s;
        bla[i][5]=insert_e;
        //file5<<bla_name[i]<<'\t'<<L1_s<<'\t'<<L1_e<<'\t'<<r_s<<'\t'<<r_e<<'\t'<<insert_s<<'\t'<<insert_e<<'\t'<<"chr"<<'\t'<<orient[i]<<chr<<endl;
    }
    
//read_connector
    
    //***********customized value***********
    int S=100;//Segmental hit gap in one read
    if (t=="ALU"){
        S=6;
    }
    
    for(int i=0;i!=blast;i++){
        if(bla[i][6]==0){
            int flag_bn=1;
            int insert_s=bla[i][4];
            int insert_e=bla[i][5];
            int L1_s=bla[i][0];
            int L1_e=bla[i][1];
            int r_s=bla[i][2];
            int r_e=bla[i][3];
            int L1_s_ex1=L1_s-S;
            int L1_s_ex2=L1_s+S;
            int L1_e_ex1=L1_e-S;
            int L1_e_ex2=L1_e+S;
            int r_s_ex1=r_s-S;
            int r_s_ex2=r_s+S;
            int r_e_ex1=r_e-S;
            int r_e_ex2=r_e+S;
            bla[i][6]=1;
            
            for(;flag_bn==1;){
                flag_bn=0;
                for(int j=0;j!=blast;j++){
                    if(bla_name[i]==bla_name[j]&&orient[i]==orient[j]&&bla[j][6]==0){
                        if(insert_e==bla[j][5]&&insert_s==bla[j][4]){
                            if(bla[j][0]<=L1_e_ex2&&bla[j][0]>=L1_e_ex1){
                                if(orient[i]=="+"&&bla[j][2]<=r_e_ex2&&bla[j][2]>=r_e_ex1){
                                    flag_bn=1;
                                    bla[j][6]=1;
                                    L1_e=bla[j][1];
                                    L1_e_ex1=L1_e-S;
                                    L1_e_ex2=L1_e+S;
                                    
                                    r_e=bla[j][3];
                                    r_e_ex1=r_e-S;
                                    r_e_ex2=r_e+S;
                                }
                                if(orient[i]=="-"&&bla[j][3]<=r_s_ex2&&bla[j][3]>=r_s_ex1){
                                    flag_bn=1;
                                    bla[j][6]=1;
                                    
                                    L1_e=bla[j][1];
                                    L1_e_ex1=L1_e-S;
                                    L1_e_ex2=L1_e+S;
                                    
                                    r_s=bla[j][2];
                                    r_s_ex1=r_s-S;
                                    r_s_ex2=r_s+S;

                                }
                            }
                            
                            else if(bla[j][1]<=L1_s_ex2&&bla[j][1]>=L1_s_ex1){
                                if(orient[i]=="+"&&bla[j][3]<=r_s_ex2&&bla[j][3]>=r_s_ex1){
                                    flag_bn=1;
                                    bla[j][6]=1;
                                    
                                    L1_s=bla[j][0];
                                    L1_s_ex1=L1_s-S;
                                    L1_s_ex2=L1_s+S;
                                    
                                    r_s=bla[j][2];
                                    r_s_ex1=r_s-S;
                                    r_s_ex2=r_s+S;
                                }
                                if(orient[i]=="-"&&bla[j][2]<=r_e_ex2&&bla[j][2]>=r_e_ex1){
                                    flag_bn=1;
                                    bla[j][6]=1;
                                    
                                    L1_s=bla[j][0];
                                    L1_s_ex1=L1_s-S;
                                    L1_s_ex2=L1_s+S;
                                    
                                    r_e=bla[j][3];
                                    r_e_ex1=r_e-S;
                                    r_e_ex2=r_e+S;
                                }
                            }
                            
                        }
                    }
                }
            }
            
            
            if(t=="LINE"&&L1_s<=5998&&L1_e>=6022){
                file5<<bla_name[i]<<'\t'<<L1_s<<'\t'<<L1_e<<'\t'<<r_s<<'\t'<<r_e<<'\t'<<insert_s<<'\t'<<insert_e<<'\t'<<chr_fix<<'\t'<<orient[i]<<endl;
            }
            else if(t!="LINE"){
                file5<<bla_name[i]<<'\t'<<L1_s<<'\t'<<L1_e<<'\t'<<r_s<<'\t'<<r_e<<'\t'<<insert_s<<'\t'<<insert_e<<'\t'<<chr_fix<<'\t'<<orient[i]<<endl;
            }
        }
    }
    
    for(int i=0;i!=blast;i++){
        delete [] bla[i];
        //delete [] sam_loc[i];
    }
    delete [] bla;
    
    for(int i=0;i!=line;i++){
        delete [] read[i];
        //delete [] sam_loc[i];
    }
    delete [] read;
    
    delete [] bla_name;
    delete [] orient;
    delete [] read_loc;
    
    return 0;
}
