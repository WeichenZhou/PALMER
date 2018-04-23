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

int ReadMasker(string WD_dir){
//int main(){
    //string WD_dir="/Users/zhouweichen/Documents/workspace/16.10.08.Pacbio_analysis/17.11.24.PALMER/v1.1.1/scp/IO_test/";
    //cout<<"Read Masker Step is now running."<<endl;
    
    ifstream file1;
    ifstream file2;
    ifstream file3;
    ofstream file5;
    ofstream file6;
    ofstream file7;
    
    string sys_regionsam = WD_dir+"region.sam";
    char *syst_regionsam = new char[sys_regionsam.length()+1];
    strcpy(syst_regionsam, sys_regionsam.c_str());
    
    file1.open(syst_regionsam);
    
    string sys_RMselec = WD_dir+"RM.selected";
    char *syst_RMselec = new char[sys_RMselec.length()+1];
    strcpy(syst_RMselec, sys_RMselec.c_str());
    
    file2.open(syst_RMselec);
    
    string sys_cigar2 = WD_dir+"cigar.2";
    char *syst_cigar2 = new char[sys_cigar2.length()+1];
    strcpy(syst_cigar2, sys_cigar2.c_str());
    
    //string sys_cigar1 = WD_dir+"cigar.1";
    //char *syst_cigar1 = new char[sys_cigar1.length()+1];
    //strcpy(syst_cigar1, sys_cigar1.c_str());
    
    file5.open(syst_cigar2);
    
    string sys_SEQ = WD_dir+"SEQ.masked";
    char *syst_SEQ = new char[sys_SEQ.length()+1];
    strcpy(syst_SEQ, sys_SEQ.c_str());
    
    file6.open(syst_SEQ);
    
    string sys_selecinfo = WD_dir+"selected.reads.info";
    char *syst_selecinfo = new char[sys_selecinfo.length()+1];
    strcpy(syst_selecinfo, sys_selecinfo.c_str());
    
    file7.open(syst_selecinfo);
    
    if (!file1.is_open()||!file2.is_open())
    {
        cout <<"CANNOT OPEN FILES"<< endl;
        exit(1);
    }
    
    string input;
    int input_int;
    
    int line=0;
    
    for(int i=0;!file1.eof();i++){
        getline(file1,input);
        line=i;
    }
    //cout<<line<<endl;
    
    int line_rm=0;
    for(int i=0;!file2.eof();i++){
        getline(file2,input);
        line_rm=i;
    }
    //cout<<line_rm<<endl;
    file1.close();
    file1.clear();
    file1.open(syst_regionsam);
    
    file2.close();
    file2.clear();
    file2.open(syst_RMselec);
    
    string **sam_info;
    sam_info=new string*[line];
    for(int i=0;i!=line;i++) sam_info[i]=new string[6];
    
    int **sam_loc;
    sam_loc=new int*[line];
    for(int i=0;i!=line;i++) sam_loc[i]=new int[5];
    
    for(int i=0;i!=line;i++){
        file1>>sam_info[i][0];        //QNAME
        file1>>sam_loc[i][0];         //FLAG
        file1>>sam_info[i][1];        //RNAME
        file1>>sam_loc[i][1];         //POS
        file1>>sam_loc[i][2];         //MAPQ
        file1>>sam_info[i][2];        //CIGAR
        file1>>sam_info[i][3];        //RNEXT
        file1>>sam_loc[i][3];         //PNEXT
        file1>>sam_loc[i][4];         //TLEN
        file1>>sam_info[i][4];        //SEQ
        getline(file1,sam_info[i][5]);//QUAL
        
    }
    
    int **rm_loc;
    rm_loc=new int*[line_rm];
    for(int i=0;i!=line_rm;i++) rm_loc[i]=new int[2];
    
    string **rm_info;
    rm_info=new string*[line_rm];
    for(int i=0;i!=line_rm;i++) rm_info[i]=new string[4];
    
    for(int i=0;i!=line_rm;i++){
        file2>>rm_info[i][0];           //chr
        file2>>rm_loc[i][0];            //start
        file2>>rm_loc[i][1];            //end
        file2>>rm_info[i][1];           //class
        file2>>rm_info[i][2];           //family
        file2>>rm_info[i][3];           //name
    }
//cout<<"region reading complete"<<endl;
    for(int i=0;i!=line;i++){
        //cout<<i<<endl;
        int number=0;
        int sc1, sc2, I, M, D, F, h1, h2, eq;
        sc1=sc2=I=M=D=F=h1=h2=eq=0;
        char cig;
        int read_s, read_e;
        read_s=sam_loc[i][1];
        
        
        //ofstream file4;
        //file4.open(syst_cigar1);
        //cout<<sam_info[i][2]<<"0E"<<endl;
        
        //file4<<sam_info[i][2]<<"0E"<<endl; //CIGAR file
        //file4.close();
        //file4.clear();
        //file3.open(syst_cigar1);
        
        string cigar1;
        cigar1=sam_info[i][2]+"0E";
        stringstream ss_cigar1;
        ss_cigar1.clear();
        ss_cigar1.str(cigar1);
        
//reads cigar&seq
        
        const int le=sam_info[i][4].length();
        char seq[le];
        const char *cseq=sam_info[i][4].c_str();
        strcpy(seq,cseq);
        
        ss_cigar1>>number;
        ss_cigar1>>cig;
        
        for(;!(cig=='E'&&number==0);){
            if(cig=='S'&&F==0) {sc1=sc1+number;F=1;}
            else if(cig=='H'&&F==0) h1=h1+number;
            else if (cig=='M') {M=M+number;F=1;}
            else if (cig=='='||cig=='X') {eq=eq+number;F=1;}
            else if (cig=='I') {I=I+number;F=1;}
            else if (cig=='D') {D=D+number;F=1;}
            else if (cig=='S'&&F==1) sc2=sc2+number;
            else if (cig=='H'&&F==1) h2=h2+number;
            ss_cigar1>>number;
            ss_cigar1>>cig;
        }
        
        //file3.close();
        //file3.clear();
        read_e=sam_loc[i][1]+M+D+eq;
        
        
        if(sc1>=50||sc2>=50||I>=50){
            
//mask
            int start=0;
            int end=0;
            
            for(int j=0;j!=line_rm;j++){
                if(!((read_s>rm_loc[j][1])||(read_e<rm_loc[j][0]))){
                    
                    stringstream ss_cigar1;
                    ss_cigar1.clear();
                    ss_cigar1.str(cigar1);
                    
                    //file3.open(syst_cigar1);
                    if(read_s>rm_loc[j][0]) start=read_s;
                    else start= rm_loc[j][0];
                    if(read_e>rm_loc[j][1]) end=rm_loc[j][1];
                    else end=read_e;
                    
                    ss_cigar1>>number;
                    ss_cigar1>>cig;
                    int bit=read_s; //bit = reference location
                    int k=0;        //k = read seq location
                    
                    
                    for(;!(cig=='E'&&number==0);){
                        if(cig=='M'||cig=='='||cig=='X') {
                            bit=bit+number;
                            
                            //cout<<start<<" "<<end<<endl;
                            //getchar();
                            //cout<<bit<<" "<<k<<" "<<le<<endl;
                            //cout<<(k+number-1)<<" "<<(k+bit-start-1)<<" "<<(k+end-bit+number-1)<<endl;
                            //cout<<"M"<<endl;
                            //cout<<cig<<" "<<number<<endl;
                            if(bit>=start&&bit<=end&&(bit-number)>=start){
                                //cout<<"M1"<<endl;
                                for(int x=0;x!=number;x++){
                                    seq[k+x]='N';
                                }
                            }
                            else if(bit>=start&&bit<=end&&(bit-number)<start){
                                //cout<<"M2"<<endl;
                                for(int x=start-bit+number;x!=number;x++){
                                    seq[k+x]='N';
                                }
                            }
                            else if((bit-number)>=start&&(bit-number)<=end&&bit>end){
                                //cout<<"M3"<<endl;
                                for(int x=0;x!=end-bit+number+1;x++){
                                    seq[k+x]='N';
                                }
                            }
                            k=k+number;
                        }
                        else if(cig=='D') {bit=bit+number;
                            //cout<<"D"<<endl;
                        }
                        else if(cig=='S'||cig=='I') {k=k+number;
                            //cout<<"S/I"<<endl;
                        }
                        //cout<<cig<<" "<<number<<endl;
                        //cout<<bit<<" "<<k<<" "<<le<<endl;
                        ss_cigar1>>number;
                        ss_cigar1>>cig;
                    }
                    //file3.close();
                    //file3.clear();
                }
            }
//output
            file5<<sam_info[i][2]<<"0E"<<endl;
            file6<<">"<<sam_info[i][0]<<'\t'<<sam_loc[i][1]<<endl;
            for(int x=0;x!=le;x++) file6<<seq[x];
            file6<<endl;
            file7<<sam_info[i][0]<<'\t'<<sam_info[i][1]<<'\t'<<sam_loc[i][1]<<endl;
        }
    }
    return 0;
}
