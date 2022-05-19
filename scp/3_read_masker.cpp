////copyright by ArthurZhou @ UMich&Fudan&HUST
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

int ReadMasker(string WD_dir, string mode){
    
    //std::ios::sync_with_stdio(false);
    //std::cin.tie(0);
    
    ifstream file1;
    //ifstream file2;
    ifstream file3;
    ofstream file5;
    ofstream file6;
    ofstream file7;
    
    string sys_regionsam = WD_dir+"region.sam";
    char *syst_regionsam = new char[sys_regionsam.length()+1];
    strcpy(syst_regionsam, sys_regionsam.c_str());
    
    file1.open(syst_regionsam);
    
    /*
    string sys_RMselec = WD_dir+"RM.selected";
    char *syst_RMselec = new char[sys_RMselec.length()+1];
    strcpy(syst_RMselec, sys_RMselec.c_str());
    
    file2.open(syst_RMselec);
    */
    
    string sys_cigar2 = WD_dir+"cigar.2";
    char *syst_cigar2 = new char[sys_cigar2.length()+1];
    strcpy(syst_cigar2, sys_cigar2.c_str());
    
    file5.open(syst_cigar2);
    
    string sys_SEQ = WD_dir+"SEQ.masked";
    char *syst_SEQ = new char[sys_SEQ.length()+1];
    strcpy(syst_SEQ, sys_SEQ.c_str());
    
    file6.open(syst_SEQ);
    
    string sys_selecinfo = WD_dir+"selected.reads.info";
    char *syst_selecinfo = new char[sys_selecinfo.length()+1];
    strcpy(syst_selecinfo, sys_selecinfo.c_str());
    
    file7.open(syst_selecinfo);
    
    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'region.sam'"<< endl;
        //exit(1);
        return 0;
    }
    
    string input;
    int input_int;
    
    int line=0;
    
    for(int i=0;!file1.eof();++i){
        getline(file1,input);
        line=i;
    }
    file1.close();
    file1.clear();
    file1.open(syst_regionsam);

    string **sam_info;
    sam_info=new string*[line];
    for(int i=0;i!=line;++i) sam_info[i]=new string[6];
    
    int **sam_loc;
    sam_loc=new int*[line];
    for(int i=0;i!=line;++i) sam_loc[i]=new int[5];
    
    for(int i=0;i!=line;++i){
        file1>>sam_info[i][0];        //QNAME
        //getline(file1,sam_info[i][0],'\t'); //QNAME
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
    
    /*
    int line_rm;
    int **rm_loc;
    string **rm_info;
    
    if (!file2.is_open()){
        line_rm=1;
        
        rm_loc=new int*[line_rm];
        for(int i=0;i!=line_rm;++i) rm_loc[i]=new int[2];
        
        rm_info=new string*[line_rm];
        for(int i=0;i!=line_rm;++i) rm_info[i]=new string[4];
        
        rm_info[0][0]="null";           //chr
        rm_loc[0][0]=0;            //start
        rm_loc[0][1]=0;            //end
        rm_info[0][1]="null";           //class
        rm_info[0][2]="null";           //family
        rm_info[0][3]="null";          //name
        
    }
    else {
        line_rm=0;
        for(int i=0;!file2.eof();++i){
            getline(file2,input);
            line_rm=i;
        }
        file2.close();
        file2.clear();
        file2.open(syst_RMselec);
        
        rm_loc=new int*[line_rm];
        for(int i=0;i!=line_rm;++i) rm_loc[i]=new int[2];
        
        rm_info=new string*[line_rm];
        for(int i=0;i!=line_rm;++i) rm_info[i]=new string[4];
        
        for(int i=0;i!=line_rm;++i){
            file2>>rm_info[i][0];           //chr
            file2>>rm_loc[i][0];            //start
            file2>>rm_loc[i][1];            //end
            file2>>rm_info[i][1];           //class
            file2>>rm_info[i][2];           //family
            file2>>rm_info[i][3];           //name
        }
    }
     */
    
//cout<<"region reading complete"<<endl;
    for(int i=0;i!=line;++i){
        //cout<<i<<endl;
        int number=0;
        int sc1, sc2, I, M, D, F, h1, h2, eq;
        sc1=sc2=I=M=D=F=h1=h2=eq=0;
        char cig;
        int read_s, read_e;
        read_s=sam_loc[i][1];
       
        string cigar1;
        cigar1=sam_info[i][2]+"0E";
        stringstream ss_cigar1;
        ss_cigar1.clear();
        ss_cigar1.str(cigar1);
        
//reads cigar&seq
        int le=sam_info[i][4].length();
        char *seq =new char[sam_info[i][4].length()+1];
        strcpy(seq,sam_info[i][4].c_str());
        
        /*
        int le=sam_info[i][4].length();
        char seq[le];
        char *cseq=sam_info[i][4].c_str();
        strcpy(seq,cseq);
        */
        ss_cigar1>>number;
        ss_cigar1>>cig;
        
        for(;!(cig=='E'&&number==0);){
            if(cig=='S'&&F==0) {sc1=sc1+number;F=1;}
            else if(cig=='H'&&F==0) h1=h1+number;
            else if (cig=='M') {M=M+number;F=1;}
            else if (cig=='='||cig=='X') {eq=eq+number;F=1;}
            else if (cig=='I') {I=I+number;F=1;}
            else if (cig=='D'||cig=='N') {D=D+number;F=1;}
            else if (cig=='S'&&F==1) sc2=sc2+number;
            else if (cig=='H'&&F==1) h2=h2+number;
            ss_cigar1>>number;
            ss_cigar1>>cig;
        }
        read_e=sam_loc[i][1]+M+D+eq-1;
        
        
        //if(sc1>=50||sc2>=50||I>=50){
            
//mask
            int start=0;
            int end=0;
            
            /*
            for(int j=0;j!=line_rm;++j){
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
                            
                            if(bit>=start&&bit<=end&&(bit-number)>=start){
                                for(int x=0;x!=number;++x){
                                    seq[k+x]='N';
                                }
                            }
                            else if(bit>=start&&bit<=end&&(bit-number)<start){
                                for(int x=start-bit+number;x!=number;++x){
                                    seq[k+x]='N';
                                }
                            }
                            else if((bit-number)>=start&&(bit-number)<=end&&bit>end){
                                for(int x=0;x!=end-bit+number+1;++x){
                                    seq[k+x]='N';
                                }
                            }
                            else if((bit-number)<start&&bit>end){
                                for(int x=start-bit+number;x!=end-bit+number+1;++x){
                                    seq[k+x]='N';
                                }
                            }
                            k=k+number;
                        }
                        else if(cig=='D'||cig=='N') {bit=bit+number;
                        }
                        else if(cig=='S'||cig=='I') {k=k+number;
                        }
                        ss_cigar1>>number;
                        ss_cigar1>>cig;
                    }
                }
            }
            */
            
            stringstream ss_cigar2;
            ss_cigar2.clear();
            ss_cigar2.str(cigar1);
            
            //file3.open(syst_cigar1);
            start=read_s;
            end=read_e;
            
            ss_cigar2>>number;
            ss_cigar2>>cig;
            
            int bit_s=0;
            int k_s=0;
            for(;!(cig=='E'&&number==0);){
                /*if(cig=='S') {
                   
                    for(int x=0;x!=number;++x){
                        seq[k_s+x]='N';
                    }
                    
                    k_s=k_s+number;
                }
                
                ss_cigar1>>number;
                ss_cigar1>>cig;
                */
                
                if(cig=='M'||cig=='=') {
                    
                    for(int x=0;x!=number;++x){
                        seq[k_s+x]='N';
                    }
                    
                    bit_s=bit_s+number;
                    k_s=k_s+number;
                }
                else if(cig=='X') {
                    bit_s=bit_s+number;
                    k_s=k_s+number;
                }
                else if(cig=='D'||cig=='N') {bit_s=bit_s+number;
                }
                
                else if(cig=='I') {k_s=k_s+number;
                }
                
                else if(cig=='S'&&mode=="asm") {
                   
                    for(int x=0;x!=number;++x){
                        seq[k_s+x]='N';
                    }
                    
                    k_s=k_s+number;
                }
                
                else if(cig=='S'){
                    k_s=k_s+number;
                }
                
                ss_cigar2>>number;
                ss_cigar2>>cig;
                
            }
            
            
            
            
//output
            file5<<sam_info[i][2]<<"0E"<<'\t'<<i<<endl;
            file6<<">"<<sam_info[i][0]<<"_"<<sam_loc[i][1]<<"_"<<i<<endl;
            for(int x=0;x!=le;++x) file6<<seq[x];
            file6<<endl;
            file7<<sam_info[i][0]<<'\t'<<sam_info[i][1]<<'\t'<<sam_loc[i][1]<<'\t'<<le<<'\t'<<i<<endl;
        //}
        delete [] seq;
    }
    
    for(int i=0;i!=line;++i){
        delete [] sam_info[i];
        delete [] sam_loc[i];
    }
    delete [] sam_info;
    delete [] sam_loc;
    
    /*
    for(int i=0;i!=line_rm;++i){
        delete [] rm_info[i];
        delete [] rm_loc[i];
    }
    delete [] rm_info;
    delete [] rm_loc;
    */
    return 0;
}
