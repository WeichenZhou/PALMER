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

using namespace std;

int calling(string WD_dir, string t, int tsd_index){
    
    //std::ios::sync_with_stdio(false);
    //std::cin.tie(0);
    
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
        BIN_5=4000;
    }
    else if (t=="HERVK"){
        BIN_3=50;
        BIN_5=50;
        tsd_index=0;
    }
    
    int le_5,le_3;
    le_5=BIN_5+1;
    le_3=BIN_3+1;
    
    ifstream file1;
    
    string sys_input = WD_dir+"read_result_TSD.txt";
    char *syst_input = new char[sys_input.length()+1];
    strcpy(syst_input, sys_input.c_str());
    
    file1.open(syst_input);
    
    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'read_result_TSD.txt'"<< endl;
        //exit(1);
        return 0;
    }
    
    ifstream file99;
    string sys_line = WD_dir+"read_result_ins_seq.txt";
    char *syst_line = new char[sys_line.length()+1];
    strcpy(syst_line, sys_line.c_str());
    file99.open(syst_line);
    
    if (!file99.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'read_result_ins_seq.txt'"<< endl;
        //exit(1);
        return 0;
    }
    
    ofstream file2;
    
    string sys_calls = WD_dir+"calls.txt";
    char *syst_calls = new char[sys_calls.length()+1];
    strcpy(syst_calls, sys_calls.c_str());
    file2.open(syst_calls);
    
    int line;
    string input;
    for(int i=0;!file1.eof();++i){
        getline(file1,input);
        line=i;
    }
    
    file1.close();
    file1.clear();
    file1.open(syst_input);
    
    int line_seq;
    //string input;
    for(int i=0;!file99.eof();++i){
        getline(file99,input);
        line_seq=i;
    }
    
    file99.close();
    file99.clear();
    file99.open(syst_line);
    
    string **info_line;
    info_line=new string *[line_seq];
    for(int i=0;i!=line_seq;++i) info_line[i]=new string[3];
    
    int **loc_TP_line;
    loc_TP_line=new int*[line_seq];
    for(int i=0;i!=line_seq;++i) loc_TP_line[i]=new int[7];
    
    for(int i=0;i!=line;++i){
        file99>>info_line[i][0];  //probe name
        file99>>input;
        file99>>input;
        file99>>input;
        file99>>input;
        file99>>input;
        file99>>input;
        file99>>input;
        file99>>info_line[i][1]; //orientation
        file99>>input;
        file99>>info_line[i][2]; //seq
        file99>>loc_TP_line[i][0];
        file99>>loc_TP_line[i][1];
        file99>>loc_TP_line[i][2];
        file99>>loc_TP_line[i][3];
        file99>>loc_TP_line[i][4];
        file99>>loc_TP_line[i][5];
        file99>>loc_TP_line[i][6];
    }
    
    
    string **info;
    info=new string *[line];
    for(int i=0;i!=line;++i) info[i]=new string[4];
    
    int **loc;
    loc=new int*[line];
    for(int i=0;i!=line;++i) loc[i]=new int[7];
    
    int **loc_TP;
    loc_TP=new int*[line];
    for(int i=0;i!=line;++i) loc_TP[i]=new int[7];
    
    string *orien;
    orien= new string[line];
    
    for(int i=0;i!=line;++i){
        file1>>info[i][0];  //probe name
        file1>>loc[i][0];   //L_loc
        file1>>loc[i][1];   //L_loc
        file1>>loc[i][2];   //R_loc
        file1>>loc[i][3];   //R_loc
        file1>>loc[i][4];   //G_loc
        file1>>loc[i][5];   //G_loc
        file1>>info[i][1];  //chr
        file1>>orien[i];    //orientation
        file1>>info[i][2]; //TSD 5' seq
        file1>>info[i][3]; //TSD 3' seq
        loc[i][6]=0;
        file1>>loc_TP[i][0];
        file1>>loc_TP[i][1];
        file1>>loc_TP[i][2];
        file1>>loc_TP[i][3];
        file1>>loc_TP[i][4];
        file1>>loc_TP[i][5];
        file1>>loc_TP[i][6];
    }
    
    //TSD_calling
    ifstream file3;
    string sys_input_TSD = WD_dir+"TSD_blastn.txt";
    char *syst_input_TSD = new char[sys_input_TSD.length()+1];
    strcpy(syst_input_TSD, sys_input_TSD.c_str());
    file3.open(syst_input_TSD);
    
    if (!file3.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'TSD_blastn.txt'"<< endl;
        //exit(1);
        return 0;
    }
    
    ofstream file4;
    string sys_output_TSD = WD_dir+"TSD_output.txt";
    char *syst_output_TSD = new char[sys_output_TSD.length()+1];
    strcpy(syst_output_TSD, sys_output_TSD.c_str());
    file4.open(syst_output_TSD);
    
    int line_tsd;
    for(int i=0;!file3.eof();++i){
        getline(file3,input);
        line_tsd=i;
    }
    file3.close();
    file3.clear();
    file3.open(syst_input_TSD);
    
    string *info_tsd;
    info_tsd= new string[line_tsd];
    
    int **loc_tsd;
    loc_tsd=new int*[line_tsd];
    for(int i=0;i!=line_tsd;++i) loc_tsd[i]=new int[12];
    
    string *kmer_tsd;
    kmer_tsd= new string[line_tsd];
    
    for(int i=0;i!=line_tsd;++i){
        file3>>info_tsd[i];
        file3>>loc_tsd[i][0];
        file3>>loc_tsd[i][1];
        file3>>loc_tsd[i][2];
        file3>>loc_tsd[i][3];
        file3>>loc_tsd[i][6];
        file3>>loc_tsd[i][7];
        loc_tsd[i][4]=0;
        loc_tsd[i][5]=0;
        loc_tsd[i][8]=0;
        file3>>loc_tsd[i][9]; //5' FP kmer
        file3>>loc_tsd[i][10];  //junction FP
        file3>>loc_tsd[i][11];  //junction FP _2
        file3>>kmer_tsd[i];
    }
    
    //calling
    
    //adjacent gap values
    //***********parameters detected***********
    int S=150;
    int L=50;
    
    if(t=="ALU"){
        S=15;
        L=5;
    }
    else if(t=="SVA"){
        S=75;
        L=25;
    }
    else if(t=="HERVK"){
        S=225;
        L=75;
    }
    
    for(int i=0;i!=line;++i){
        if(loc[i][6]!=-1){
            
            int start1=loc[i][4]-S;
            int start2=loc[i][4]+S;
            int end2=loc[i][5]+S;
            int end1=loc[i][5]-S;
            int L1_s1=loc[i][0]-L;
            int L1_s2=loc[i][0]+L;
            int L1_e1=loc[i][1]-L;
            int L1_e2=loc[i][1]+L;
            
            //int TP=loc_TP[i][0];
            int start_TP1=loc_TP[i][5]-S;
            int start_TP2=loc_TP[i][5]+S;
            int end_TP2=loc_TP[i][6]+S;
            int end_TP1=loc_TP[i][6]-S;
            int L1_s_TP1=loc_TP[i][1]-L;
            int L1_s_TP2=loc_TP[i][1]+L;
            int L1_e_TP1=loc_TP[i][2]-L;
            int L1_e_TP2=loc_TP[i][2]+L;
            
            int number=0;
            int number_true=0;
            int number_all=1;
            loc[i][6]=-1;
            int flag=1;
            
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
            
            string seq_index_a;
            
            seq_index_a=info[i][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+orien[i]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
            
            for(int w=0;w!=line_tsd;++w){
                if(info_tsd[w]==seq_index_a){
                    loc_tsd[w][8]=1;
                }
            }
            
            int flag_left=0;
            int flag_right=0;
            
            for(;flag==1;){
                flag=0;
                
                for(int j=0;j!=line;++j){

                    if(info[i][0]==info[j][0]&&loc[j][6]==0&&info[i][1]==info[j][1]&&loc[i][0]==loc[j][0]&&loc[i][1]==loc[j][1]&&loc[i][2]==loc[j][2]&&loc[i][3]==loc[j][3]&&loc[i][4]==loc[j][4]&&loc[i][5]==loc[j][5]&&orien[i]==orien[j]&&loc_TP[i][0]==loc_TP[j][0]&&loc_TP[i][1]==loc_TP[j][1]&&loc_TP[i][2]==loc_TP[j][2]&&loc_TP[i][3]==loc_TP[j][3]&&loc_TP[i][4]==loc_TP[j][4]&&loc_TP[i][5]==loc_TP[j][5]&&loc_TP[i][6]==loc_TP[j][6]){
                        loc[j][6]=-1;
                        
                    }
                    else if(info[i][0]!=info[j][0]&&loc[j][6]==0&&info[i][1]==info[j][1]&&orien[i]==orien[j]&&loc_TP[i][0]==loc_TP[j][0]){
                        if((loc[j][4]>=start1&&loc[j][4]<=start2)||(loc[j][5]<=end2&&loc[j][5]>=end1)){
                            
                            //left&right
                            if((loc[j][0]>=L1_s1&&loc[j][0]<=L1_s2)&&(loc[j][1]>=L1_e1&&loc[j][1]<=L1_e2)&&(loc_TP[j][1]>=L1_s_TP1&&loc_TP[j][1]<=L1_s_TP2)&&(loc_TP[j][2]>=L1_e_TP1&&loc_TP[j][2]<=L1_e_TP2)){
                                if(loc[j][4]-S<start1) start1=loc[j][4]-S;
                                if(loc[j][4]+S>start2) start2=loc[j][4]+S;
                                if(loc[j][5]-S<end1) end1=loc[j][5]-S;
                                if(loc[j][5]+S>end2) end2=loc[j][5]+S;
                                
                                if(loc[j][0]-L<L1_s1) L1_s1=loc[j][0]-L;
                                if(loc[j][0]+L>L1_s2) L1_s2=loc[j][0]+L;
                                if(loc[j][1]-L<L1_e1) L1_e1=loc[j][1]-L;
                                if(loc[j][1]+L>L1_e2) L1_e2=loc[j][1]+L;
                                
                                if(loc_TP[j][5]-S<start_TP1) start_TP1=loc_TP[j][5]-S;
                                if(loc_TP[j][5]+S>start_TP2) start_TP2=loc_TP[j][5]+S;
                                if(loc_TP[j][6]-S<end_TP1) end_TP1=loc_TP[j][6]-S;
                                if(loc_TP[j][6]+S>end_TP2) end_TP2=loc_TP[j][6]+S;
                                
                                if(loc_TP[j][1]-L<L1_s_TP1) L1_s_TP1=loc_TP[j][1]-L;
                                if(loc_TP[j][1]+L>L1_s_TP2) L1_s_TP2=loc_TP[j][1]+L;
                                if(loc_TP[j][2]-L<L1_e_TP1) L1_e_TP1=loc_TP[j][2]-L;
                                if(loc_TP[j][2]+L>L1_e_TP2) L1_e_TP2=loc_TP[j][2]+L;
                                
                                flag=1;
                                loc[j][6]=-1;
                                number_all=number_all+1;
                                
                                string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                                stringstream ss_0;
                                ss_0.clear();
                                ss_0<<loc[j][0];
                                loc_0=ss_0.str();
                                stringstream ss_1;
                                ss_1.clear();
                                ss_1<<loc[j][1];
                                loc_1=ss_1.str();
                                stringstream ss_2;
                                ss_2.clear();
                                ss_2<<loc[j][2];
                                loc_2=ss_2.str();
                                stringstream ss_3;
                                ss_3.clear();
                                ss_3<<loc[j][3];
                                loc_3=ss_3.str();
                                stringstream ss_4;
                                ss_4.clear();
                                ss_4<<loc[j][4];
                                loc_4=ss_4.str();
                                stringstream ss_5;
                                ss_5.clear();
                                ss_5<<loc[j][5];
                                loc_5=ss_5.str();
                                
                                string loc_TP_0, loc_TP_1, loc_TP_2, loc_TP_3, loc_TP_4, loc_TP_5, loc_TP_6;
                                stringstream ss_TP_0;
                                ss_TP_0.clear();
                                ss_TP_0<<loc_TP[j][0];
                                loc_TP_0=ss_TP_0.str();
                                stringstream ss_TP_1;
                                ss_TP_1.clear();
                                ss_TP_1<<loc_TP[j][1];
                                loc_TP_1=ss_TP_1.str();
                                stringstream ss_TP_2;
                                ss_TP_2.clear();
                                ss_TP_2<<loc_TP[j][2];
                                loc_TP_2=ss_TP_2.str();
                                stringstream ss_TP_3;
                                ss_TP_3.clear();
                                ss_TP_3<<loc_TP[j][3];
                                loc_TP_3=ss_TP_3.str();
                                stringstream ss_TP_4;
                                ss_TP_4.clear();
                                ss_TP_4<<loc_TP[j][4];
                                loc_TP_4=ss_TP_4.str();
                                stringstream ss_TP_5;
                                ss_TP_5.clear();
                                ss_TP_5<<loc_TP[j][5];
                                loc_TP_5=ss_TP_5.str();
                                stringstream ss_TP_6;
                                ss_TP_6.clear();
                                ss_TP_6<<loc_TP[j][6];
                                loc_TP_6=ss_TP_6.str();
                                
                                string seq_index_b;
                                
                                seq_index_b=info[j][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                                
                                for(int w=0;w!=line_tsd;++w){
                                    if(info_tsd[w]==seq_index_b){
                                        loc_tsd[w][8]=1;
                                    }
                                }
                                
                            }
                            //left
                            else if((loc[j][0]>=L1_s1&&loc[j][0]<=L1_s2)&&(loc[j][1]<L1_e1||loc[j][1]>L1_e2)){
                                if(loc[j][1]>(L1_e1+L)&&flag_right==0) {
                                    L1_e1=loc[j][1]-L;
                                    L1_e2=loc[j][1]+L;
                                    flag_left=1;
                                    if(loc[j][4]-S<start1) start1=loc[j][4]-S;
                                    if(loc[j][4]+S>start2) start2=loc[j][4]+S;
                                    if(loc[j][5]-S<end1) end1=loc[j][5]-S;
                                    if(loc[j][5]+S>end2) end2=loc[j][5]+S;
                                    
                                    if(loc[j][0]-L<L1_s1) L1_s1=loc[j][0]-L;
                                    if(loc[j][0]+L>L1_s2) L1_s2=loc[j][0]+L;
                                    flag=1;
                                    loc[j][6]=-1;
                                    number_all=number_all+1;
                                    
                                    string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                                    stringstream ss_0;
                                    ss_0.clear();
                                    ss_0<<loc[j][0];
                                    loc_0=ss_0.str();
                                    stringstream ss_1;
                                    ss_1.clear();
                                    ss_1<<loc[j][1];
                                    loc_1=ss_1.str();
                                    stringstream ss_2;
                                    ss_2.clear();
                                    ss_2<<loc[j][2];
                                    loc_2=ss_2.str();
                                    stringstream ss_3;
                                    ss_3.clear();
                                    ss_3<<loc[j][3];
                                    loc_3=ss_3.str();
                                    stringstream ss_4;
                                    ss_4.clear();
                                    ss_4<<loc[j][4];
                                    loc_4=ss_4.str();
                                    stringstream ss_5;
                                    ss_5.clear();
                                    ss_5<<loc[j][5];
                                    loc_5=ss_5.str();
                                    
                                    string loc_TP_0, loc_TP_1, loc_TP_2, loc_TP_3, loc_TP_4, loc_TP_5, loc_TP_6;
                                    stringstream ss_TP_0;
                                    ss_TP_0.clear();
                                    ss_TP_0<<loc_TP[j][0];
                                    loc_TP_0=ss_TP_0.str();
                                    stringstream ss_TP_1;
                                    ss_TP_1.clear();
                                    ss_TP_1<<loc_TP[j][1];
                                    loc_TP_1=ss_TP_1.str();
                                    stringstream ss_TP_2;
                                    ss_TP_2.clear();
                                    ss_TP_2<<loc_TP[j][2];
                                    loc_TP_2=ss_TP_2.str();
                                    stringstream ss_TP_3;
                                    ss_TP_3.clear();
                                    ss_TP_3<<loc_TP[j][3];
                                    loc_TP_3=ss_TP_3.str();
                                    stringstream ss_TP_4;
                                    ss_TP_4.clear();
                                    ss_TP_4<<loc_TP[j][4];
                                    loc_TP_4=ss_TP_4.str();
                                    stringstream ss_TP_5;
                                    ss_TP_5.clear();
                                    ss_TP_5<<loc_TP[j][5];
                                    loc_TP_5=ss_TP_5.str();
                                    stringstream ss_TP_6;
                                    ss_TP_6.clear();
                                    ss_TP_6<<loc_TP[j][6];
                                    loc_TP_6=ss_TP_6.str();
                                    
                                    string seq_index_b;
                                    
                                    seq_index_b=info[j][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                                    
                                    
                                    for(int w=0;w!=line_tsd;++w){
                                        if(info_tsd[w]==seq_index_b){
                                            loc_tsd[w][8]=1;
                                        }
                                    }
                                    
                                }
                                else if(loc[j][1]<=(L1_e1+L)){
                                    flag_left=1;
                                    if(loc[j][4]-S<start1) start1=loc[j][4]-S;
                                    if(loc[j][4]+S>start2) start2=loc[j][4]+S;
                                    if(loc[j][5]-S<end1) end1=loc[j][5]-S;
                                    if(loc[j][5]+S>end2) end2=loc[j][5]+S;
                                    
                                    if(loc[j][0]-L<L1_s1) L1_s1=loc[j][0]-L;
                                    if(loc[j][0]+L>L1_s2) L1_s2=loc[j][0]+L;
                                    flag=1;
                                    loc[j][6]=-1;
                                    number_all=number_all+1;
                                    
                                    string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                                    stringstream ss_0;
                                    ss_0.clear();
                                    ss_0<<loc[j][0];
                                    loc_0=ss_0.str();
                                    stringstream ss_1;
                                    ss_1.clear();
                                    ss_1<<loc[j][1];
                                    loc_1=ss_1.str();
                                    stringstream ss_2;
                                    ss_2.clear();
                                    ss_2<<loc[j][2];
                                    loc_2=ss_2.str();
                                    stringstream ss_3;
                                    ss_3.clear();
                                    ss_3<<loc[j][3];
                                    loc_3=ss_3.str();
                                    stringstream ss_4;
                                    ss_4.clear();
                                    ss_4<<loc[j][4];
                                    loc_4=ss_4.str();
                                    stringstream ss_5;
                                    ss_5.clear();
                                    ss_5<<loc[j][5];
                                    loc_5=ss_5.str();
                                    
                                    string loc_TP_0, loc_TP_1, loc_TP_2, loc_TP_3, loc_TP_4, loc_TP_5, loc_TP_6;
                                    stringstream ss_TP_0;
                                    ss_TP_0.clear();
                                    ss_TP_0<<loc_TP[j][0];
                                    loc_TP_0=ss_TP_0.str();
                                    stringstream ss_TP_1;
                                    ss_TP_1.clear();
                                    ss_TP_1<<loc_TP[j][1];
                                    loc_TP_1=ss_TP_1.str();
                                    stringstream ss_TP_2;
                                    ss_TP_2.clear();
                                    ss_TP_2<<loc_TP[j][2];
                                    loc_TP_2=ss_TP_2.str();
                                    stringstream ss_TP_3;
                                    ss_TP_3.clear();
                                    ss_TP_3<<loc_TP[j][3];
                                    loc_TP_3=ss_TP_3.str();
                                    stringstream ss_TP_4;
                                    ss_TP_4.clear();
                                    ss_TP_4<<loc_TP[j][4];
                                    loc_TP_4=ss_TP_4.str();
                                    stringstream ss_TP_5;
                                    ss_TP_5.clear();
                                    ss_TP_5<<loc_TP[j][5];
                                    loc_TP_5=ss_TP_5.str();
                                    stringstream ss_TP_6;
                                    ss_TP_6.clear();
                                    ss_TP_6<<loc_TP[j][6];
                                    loc_TP_6=ss_TP_6.str();
                                    
                                    string seq_index_b;
                                    
                                    seq_index_b=info[j][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                                    
                                    
                                    for(int w=0;w!=line_tsd;++w){
                                        if(info_tsd[w]==seq_index_b){
                                            loc_tsd[w][8]=1;
                                        }
                                    }
                                }
                            }
                            //right
                            else if((loc[j][0]<L1_s1||loc[j][0]>L1_s2)&&(loc[j][1]>=L1_e1&&loc[j][1]<=L1_e2)){
                                if(loc[j][0]<(L1_s1+L)&&flag_left==0) {
                                    L1_s1=loc[j][0]-L;
                                    L1_s2=loc[j][0]+L;
                                    flag_right=1;
                                    if(loc[j][4]-S<start1) start1=loc[j][4]-S;
                                    if(loc[j][4]+S>start2) start2=loc[j][4]+S;
                                    if(loc[j][5]-S<end1) end1=loc[j][5]-S;
                                    if(loc[j][5]+S>end2) end2=loc[j][5]+S;
                                    
                                    if(loc[j][1]-L<L1_e1) L1_e1=loc[j][1]-L;
                                    if(loc[j][1]+L>L1_e2) L1_e2=loc[j][1]+L;
                                    flag=1;
                                    loc[j][6]=-1;
                                    number_all=number_all+1;
                                    
                                    string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                                    stringstream ss_0;
                                    ss_0.clear();
                                    ss_0<<loc[j][0];
                                    loc_0=ss_0.str();
                                    stringstream ss_1;
                                    ss_1.clear();
                                    ss_1<<loc[j][1];
                                    loc_1=ss_1.str();
                                    stringstream ss_2;
                                    ss_2.clear();
                                    ss_2<<loc[j][2];
                                    loc_2=ss_2.str();
                                    stringstream ss_3;
                                    ss_3.clear();
                                    ss_3<<loc[j][3];
                                    loc_3=ss_3.str();
                                    stringstream ss_4;
                                    ss_4.clear();
                                    ss_4<<loc[j][4];
                                    loc_4=ss_4.str();
                                    stringstream ss_5;
                                    ss_5.clear();
                                    ss_5<<loc[j][5];
                                    loc_5=ss_5.str();
                                    
                                    string loc_TP_0, loc_TP_1, loc_TP_2, loc_TP_3, loc_TP_4, loc_TP_5, loc_TP_6;
                                    stringstream ss_TP_0;
                                    ss_TP_0.clear();
                                    ss_TP_0<<loc_TP[j][0];
                                    loc_TP_0=ss_TP_0.str();
                                    stringstream ss_TP_1;
                                    ss_TP_1.clear();
                                    ss_TP_1<<loc_TP[j][1];
                                    loc_TP_1=ss_TP_1.str();
                                    stringstream ss_TP_2;
                                    ss_TP_2.clear();
                                    ss_TP_2<<loc_TP[j][2];
                                    loc_TP_2=ss_TP_2.str();
                                    stringstream ss_TP_3;
                                    ss_TP_3.clear();
                                    ss_TP_3<<loc_TP[j][3];
                                    loc_TP_3=ss_TP_3.str();
                                    stringstream ss_TP_4;
                                    ss_TP_4.clear();
                                    ss_TP_4<<loc_TP[j][4];
                                    loc_TP_4=ss_TP_4.str();
                                    stringstream ss_TP_5;
                                    ss_TP_5.clear();
                                    ss_TP_5<<loc_TP[j][5];
                                    loc_TP_5=ss_TP_5.str();
                                    stringstream ss_TP_6;
                                    ss_TP_6.clear();
                                    ss_TP_6<<loc_TP[j][6];
                                    loc_TP_6=ss_TP_6.str();
                                    
                                    string seq_index_b;
                                    
                                    seq_index_b=info[j][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                                    
                                    
                                    for(int w=0;w!=line_tsd;++w){
                                        if(info_tsd[w]==seq_index_b){
                                            loc_tsd[w][8]=1;
                                        }
                                    }
                                    
                                }
                                else if(loc[j][0]>(L1_s1+L)){
                                    flag_right=1;
                                    if(loc[j][4]-S<start1) start1=loc[j][4]-S;
                                    if(loc[j][4]+S>start2) start2=loc[j][4]+S;
                                    if(loc[j][5]-S<end1) end1=loc[j][5]-S;
                                    if(loc[j][5]+S>end2) end2=loc[j][5]+S;
                                    
                                    if(loc[j][1]-L<L1_e1) L1_e1=loc[j][1]-L;
                                    if(loc[j][1]+L>L1_e2) L1_e2=loc[j][1]+L;
                                    flag=1;
                                    loc[j][6]=-1;
                                    number_all=number_all+1;
                                    
                                    string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                                    stringstream ss_0;
                                    ss_0.clear();
                                    ss_0<<loc[j][0];
                                    loc_0=ss_0.str();
                                    stringstream ss_1;
                                    ss_1.clear();
                                    ss_1<<loc[j][1];
                                    loc_1=ss_1.str();
                                    stringstream ss_2;
                                    ss_2.clear();
                                    ss_2<<loc[j][2];
                                    loc_2=ss_2.str();
                                    stringstream ss_3;
                                    ss_3.clear();
                                    ss_3<<loc[j][3];
                                    loc_3=ss_3.str();
                                    stringstream ss_4;
                                    ss_4.clear();
                                    ss_4<<loc[j][4];
                                    loc_4=ss_4.str();
                                    stringstream ss_5;
                                    ss_5.clear();
                                    ss_5<<loc[j][5];
                                    loc_5=ss_5.str();
                                    
                                    string loc_TP_0, loc_TP_1, loc_TP_2, loc_TP_3, loc_TP_4, loc_TP_5, loc_TP_6;
                                    stringstream ss_TP_0;
                                    ss_TP_0.clear();
                                    ss_TP_0<<loc_TP[j][0];
                                    loc_TP_0=ss_TP_0.str();
                                    stringstream ss_TP_1;
                                    ss_TP_1.clear();
                                    ss_TP_1<<loc_TP[j][1];
                                    loc_TP_1=ss_TP_1.str();
                                    stringstream ss_TP_2;
                                    ss_TP_2.clear();
                                    ss_TP_2<<loc_TP[j][2];
                                    loc_TP_2=ss_TP_2.str();
                                    stringstream ss_TP_3;
                                    ss_TP_3.clear();
                                    ss_TP_3<<loc_TP[j][3];
                                    loc_TP_3=ss_TP_3.str();
                                    stringstream ss_TP_4;
                                    ss_TP_4.clear();
                                    ss_TP_4<<loc_TP[j][4];
                                    loc_TP_4=ss_TP_4.str();
                                    stringstream ss_TP_5;
                                    ss_TP_5.clear();
                                    ss_TP_5<<loc_TP[j][5];
                                    loc_TP_5=ss_TP_5.str();
                                    stringstream ss_TP_6;
                                    ss_TP_6.clear();
                                    ss_TP_6<<loc_TP[j][6];
                                    loc_TP_6=ss_TP_6.str();
                                    
                                    string seq_index_b;
                                    
                                    seq_index_b=info[j][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                                    
                                    
                                    for(int w=0;w!=line_tsd;++w){
                                        if(info_tsd[w]==seq_index_b){
                                            loc_tsd[w][8]=1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            //calculate the number of 'left&right'
            for(int j=0;j!=line;++j){
            
                
          //for SVA
                if(t=="SVA"){
                    if((loc[j][1]>=L1_e1&&loc[j][1]<=L1_e2)&&(loc[j][4]>=start1&&loc[j][4]<=start2)&&(loc[j][5]>=end1&&loc[j][5]<=end2)&&loc[j][6]==-1&&info[i][1]==info[j][1]&&orien[i]==orien[j]&&loc_TP[i][0]==loc_TP[j][0]&&(loc_TP[j][1]>=L1_s_TP1&&loc_TP[j][1]<=L1_s_TP2)&&(loc_TP[j][2]>=L1_e_TP1&&loc_TP[j][2]<=L1_e_TP2)&&(loc_TP[j][5]>=start_TP1&&loc_TP[j][5]<=start_TP2)&&(loc_TP[j][6]>=end_TP1&&loc_TP[j][6]<=end_TP2)){
                        
                        string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                        stringstream ss_0;
                        ss_0.clear();
                        ss_0<<loc[j][0];
                        loc_0=ss_0.str();
                        stringstream ss_1;
                        ss_1.clear();
                        ss_1<<loc[j][1];
                        loc_1=ss_1.str();
                        stringstream ss_2;
                        ss_2.clear();
                        ss_2<<loc[j][2];
                        loc_2=ss_2.str();
                        stringstream ss_3;
                        ss_3.clear();
                        ss_3<<loc[j][3];
                        loc_3=ss_3.str();
                        stringstream ss_4;
                        ss_4.clear();
                        ss_4<<loc[j][4];
                        loc_4=ss_4.str();
                        stringstream ss_5;
                        ss_5.clear();
                        ss_5<<loc[j][5];
                        loc_5=ss_5.str();
                        
                        string loc_TP_0, loc_TP_1, loc_TP_2, loc_TP_3, loc_TP_4, loc_TP_5, loc_TP_6;
                        stringstream ss_TP_0;
                        ss_TP_0.clear();
                        ss_TP_0<<loc_TP[j][0];
                        loc_TP_0=ss_TP_0.str();
                        stringstream ss_TP_1;
                        ss_TP_1.clear();
                        ss_TP_1<<loc_TP[j][1];
                        loc_TP_1=ss_TP_1.str();
                        stringstream ss_TP_2;
                        ss_TP_2.clear();
                        ss_TP_2<<loc_TP[j][2];
                        loc_TP_2=ss_TP_2.str();
                        stringstream ss_TP_3;
                        ss_TP_3.clear();
                        ss_TP_3<<loc_TP[j][3];
                        loc_TP_3=ss_TP_3.str();
                        stringstream ss_TP_4;
                        ss_TP_4.clear();
                        ss_TP_4<<loc_TP[j][4];
                        loc_TP_4=ss_TP_4.str();
                        stringstream ss_TP_5;
                        ss_TP_5.clear();
                        ss_TP_5<<loc_TP[j][5];
                        loc_TP_5=ss_TP_5.str();
                        stringstream ss_TP_6;
                        ss_TP_6.clear();
                        ss_TP_6<<loc_TP[j][6];
                        loc_TP_6=ss_TP_6.str();
                        
                        string seq_index_c;
                        
                        seq_index_c=info[j][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                        
                        
                        int flag_number=0;
                        int flag_true_number=0;
                        for(int w=0;w!=line_tsd;++w){
                            if(info_tsd[w]==seq_index_c&&loc_tsd[w][8]==1){
                                loc_tsd[w][4]=1;
                                flag_number=1;
                                //*******
                                if(loc_tsd[w][10]>0&&loc_tsd[w][11]>0){
                                    //&&loc_tsd[w][10]==0){
                                    flag_true_number=1;
                                }
                            }
                        }
                        if(flag_number==1){
                            number=number+1;
                            //cout<<seq_index<<endl;
                        }
                        if(flag_true_number==1){
                            number_true=number_true+1;
                        }
                    }
                }
                
                
//for ALU and LINE and CUSTOMIZED
                else {
                    if((loc[j][0]>=L1_s1&&loc[j][0]<=L1_s2)&&(loc[j][1]>=L1_e1&&loc[j][1]<=L1_e2)&&(loc[j][4]>=start1&&loc[j][4]<=start2)&&(loc[j][5]>=end1&&loc[j][5]<=end2)&&loc[j][6]==-1&&info[i][1]==info[j][1]&&orien[i]==orien[j]&&loc_TP[i][0]==loc_TP[j][0]&&(loc_TP[j][1]>=L1_s_TP1&&loc_TP[j][1]<=L1_s_TP2)&&(loc_TP[j][2]>=L1_e_TP1&&loc_TP[j][2]<=L1_e_TP2)&&(loc_TP[j][5]>=start_TP1&&loc_TP[j][5]<=start_TP2)&&(loc_TP[j][6]>=end_TP1&&loc_TP[j][6]<=end_TP2)){
                        
                        //cout<<"get in one time"<<endl;
                        
                        string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                        stringstream ss_0;
                        ss_0.clear();
                        ss_0<<loc[j][0];
                        loc_0=ss_0.str();
                        stringstream ss_1;
                        ss_1.clear();
                        ss_1<<loc[j][1];
                        loc_1=ss_1.str();
                        stringstream ss_2;
                        ss_2.clear();
                        ss_2<<loc[j][2];
                        loc_2=ss_2.str();
                        stringstream ss_3;
                        ss_3.clear();
                        ss_3<<loc[j][3];
                        loc_3=ss_3.str();
                        stringstream ss_4;
                        ss_4.clear();
                        ss_4<<loc[j][4];
                        loc_4=ss_4.str();
                        stringstream ss_5;
                        ss_5.clear();
                        ss_5<<loc[j][5];
                        loc_5=ss_5.str();
                        
                        string loc_TP_0, loc_TP_1, loc_TP_2, loc_TP_3, loc_TP_4, loc_TP_5, loc_TP_6;
                        stringstream ss_TP_0;
                        ss_TP_0.clear();
                        ss_TP_0<<loc_TP[j][0];
                        loc_TP_0=ss_TP_0.str();
                        stringstream ss_TP_1;
                        ss_TP_1.clear();
                        ss_TP_1<<loc_TP[j][1];
                        loc_TP_1=ss_TP_1.str();
                        stringstream ss_TP_2;
                        ss_TP_2.clear();
                        ss_TP_2<<loc_TP[j][2];
                        loc_TP_2=ss_TP_2.str();
                        stringstream ss_TP_3;
                        ss_TP_3.clear();
                        ss_TP_3<<loc_TP[j][3];
                        loc_TP_3=ss_TP_3.str();
                        stringstream ss_TP_4;
                        ss_TP_4.clear();
                        ss_TP_4<<loc_TP[j][4];
                        loc_TP_4=ss_TP_4.str();
                        stringstream ss_TP_5;
                        ss_TP_5.clear();
                        ss_TP_5<<loc_TP[j][5];
                        loc_TP_5=ss_TP_5.str();
                        stringstream ss_TP_6;
                        ss_TP_6.clear();
                        ss_TP_6<<loc_TP[j][6];
                        loc_TP_6=ss_TP_6.str();
                        
                        string seq_index_c;
                        
                        seq_index_c=info[j][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                        
                        
                        int flag_number=0;
                        int flag_true_number=0;
                        //cout<<"Try to get in as "<<seq_index_c<<endl;
                        for(int w=0;w!=line_tsd;++w){
                            
                            //cout<<"Try to test as "<<info_tsd[w]<<" and "<<loc_tsd[w][8]<<endl;
                            
                            
                            if(info_tsd[w]==seq_index_c&&loc_tsd[w][8]==1){
                                loc_tsd[w][4]=1;
                                flag_number=1;
                                //cout<<"get in one time "<<seq_index_c<<endl;
                                //*******
                                if(loc_tsd[w][10]>0&&loc_tsd[w][11]>0){
                                    //&&loc_tsd[w][10]==0){
                                    flag_true_number=1;
                                }
                            }
                        }
                        if(flag_number==1){
                            number=number+1;
                            //cout<<seq_index<<endl;
                        }
                        if(flag_true_number==1){
                            number_true=number_true+1;
                        }
                    }
                }
            }
            
//note for above: we only need "number_all" here
            //cout<<"number is "<<number<<endl;
            //cout<<"number_ture is "<<number_true<<endl;
            //cout<<"number_all is "<<number_all<<endl;
            
            
            if(tsd_index==1){
            
//TSD_indentifier
            string **name_tsd;
            name_tsd =new string*[line_tsd];
            for(int j=0;j!=line_tsd;++j) name_tsd[j]= new string[number];
            
            int **ls;
            ls =new int*[line_tsd];
            for(int j=0;j!=line_tsd;++j) ls[j]= new int[number];
            int **le;
            le =new int*[line_tsd];
            for(int j=0;j!=line_tsd;++j) le[j]= new int[number];
            int **rs;
            rs =new int*[line_tsd];
            for(int j=0;j!=line_tsd;++j) rs[j]= new int[number];
            int **re;
            re =new int*[line_tsd];
            for(int j=0;j!=line_tsd;++j) re[j]= new int[number];
            int **len_5;
            len_5 =new int*[line_tsd];
            for(int j=0;j!=line_tsd;++j) len_5[j]= new int[number];
            int **len_3;
            len_3 =new int*[line_tsd];
            for(int j=0;j!=line_tsd;++j) len_3[j]= new int[number];
            
            string **kmerseq_tsd;
            kmerseq_tsd =new string*[line_tsd];
            for(int j=0;j!=line_tsd;++j) kmerseq_tsd[j]= new string[number];
            
            //loc_tsd[][5]
            for(int w=0;w!=line_tsd;++w){
//*********
                if(loc_tsd[w][4]==1&&loc_tsd[w][10]>0&&loc_tsd[w][11]>0){
                   //&&loc_tsd[w][10]==0){
                    loc_tsd[w][5]=1;
                    
                    name_tsd[w][loc_tsd[w][5]-1]=info_tsd[w];
                    ls[w][loc_tsd[w][5]-1]=loc_tsd[w][0];
                    le[w][loc_tsd[w][5]-1]=loc_tsd[w][1];
                    rs[w][loc_tsd[w][5]-1]=loc_tsd[w][2];
                    re[w][loc_tsd[w][5]-1]=loc_tsd[w][3];
                    len_5[w][loc_tsd[w][5]-1]=loc_tsd[w][6];
                    len_3[w][loc_tsd[w][5]-1]=loc_tsd[w][7];
                    kmerseq_tsd[w][loc_tsd[w][5]-1]=kmer_tsd[w];
                    
                    for(int w2=0;w2!=line_tsd;++w2){
                        
                        int flag_name=0;
                        for(int w_name=0;w_name!=loc_tsd[w][5];w_name++){
                            if(name_tsd[w][w_name]==info_tsd[w2]){
                                flag_name=1;
                                break;
                            }
                        }
//**********
                        if(loc_tsd[w2][4]==1&&loc_tsd[w2][10]>0&&flag_name==0&&loc_tsd[w2][11]>0){
                           //&&loc_tsd[w2][10]==0&&flag_name==0){
                            if((!((loc_tsd[w2][0]>loc_tsd[w][1])||(loc_tsd[w2][1]<loc_tsd[w][0])))&&(!((loc_tsd[w2][2]>loc_tsd[w][3])||(loc_tsd[w2][3]<loc_tsd[w][2])))){
                                int ls_tsd=loc_tsd[w][0];
                                int le_tsd=loc_tsd[w][1];
                                int rs_tsd=loc_tsd[w][2];
                                int re_tsd=loc_tsd[w][3];
                                
                                if(loc_tsd[w2][0]>ls_tsd) ls_tsd=loc_tsd[w2][0];
                                if(loc_tsd[w2][2]>rs_tsd) rs_tsd=loc_tsd[w2][2];
                                if(loc_tsd[w2][1]<le_tsd) le_tsd=loc_tsd[w2][1];
                                if(loc_tsd[w2][3]<re_tsd) re_tsd=loc_tsd[w2][3];
                                
                                //the minimum overlap should be larger than 50%
                                if(((le_tsd-ls_tsd+1)*2>=(loc_tsd[w][1]-loc_tsd[w][0]+1))&&((re_tsd-rs_tsd+1)*2>=(loc_tsd[w][3]-loc_tsd[w][2]+1))){
                                    loc_tsd[w][5]=1+loc_tsd[w][5];
                                    name_tsd[w][loc_tsd[w][5]-1]=info_tsd[w2];
                                    ls[w][loc_tsd[w][5]-1]=loc_tsd[w2][0];
                                    le[w][loc_tsd[w][5]-1]=loc_tsd[w2][1];
                                    rs[w][loc_tsd[w][5]-1]=loc_tsd[w2][2];
                                    re[w][loc_tsd[w][5]-1]=loc_tsd[w2][3];
                                    len_5[w][loc_tsd[w][5]-1]=loc_tsd[w2][6];
                                    len_3[w][loc_tsd[w][5]-1]=loc_tsd[w2][7];
                                    kmerseq_tsd[w][loc_tsd[w][5]-1]=kmer_tsd[w2];
                                }
                            }
                        }
                    }
                    //cout<<"name_tsd[w][loc_tsd[w][5]-1]="<<name_tsd[w][loc_tsd[w][5]-1]<<endl;
                    //cout<<"loc_tsd[w][5]="<<loc_tsd[w][5]<<endl;
                }
            }
            
            
           
            
            //transD ployA
            int flag_ployA=0;
            int flag_trans=1;
            
            for(;flag_ployA==0&&flag_trans==1;){
                
                //best subreads cluster candidate
                int p=0;
                for(int w=0;w!=line_tsd;++w){
//************
                    if(loc_tsd[w][4]==1&&loc_tsd[w][10]>0&&loc_tsd[w][11]>0){
                       //&&loc_tsd[w][10]==0){
                        if(loc_tsd[w][5]>p) p=loc_tsd[w][5];
                    }
                }
                
                //cout<<"p="<<p<<endl;
                if(p==0){
                    
                    file2<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<'\t'<<info[i][1]<<'\t'<<start1+S<<'\t'<<start2-S<<'\t'<<end1+S<<'\t'<<end2-S<<'\t'<<L1_s1+L<<'\t'<<L1_s2-L<<'\t'<<L1_e1+L<<'\t'<<L1_e2-L<<'\t'<<"0"<<'\t'<<number_all<<'\t'<<(number_all-number)<<'\t'<<orien[i]<<'\t'<<"-"<<'\t'<<"0"<<'\t'<<"0"<<'\t'<<"0"<<'\t'<<loc_TP[i][0]<<'\t'<<int((L1_s_TP1+L1_s_TP2)/2)<<'\t'<<int((L1_e_TP1+L1_e_TP2)/2)<<endl;
                    file4<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<'\t'<<seq_index_a<<'\t'<<"N"<<'\t'<<"N"<<'\t'<<"N"<<'\t'<<"N"<<'\t'<<info_line[i][2]<<endl;
                    
                    
                    for(int j=0;j!=line_tsd;++j) {delete [] name_tsd[j];}
                    delete [] name_tsd;
                    
                    for(int j=0;j!=line_tsd;++j) {delete [] ls[j];}
                    delete [] ls;
                    for(int j=0;j!=line_tsd;++j) {delete [] le[j];}
                    delete [] le;
                    for(int j=0;j!=line_tsd;++j) {delete [] re[j];}
                    delete [] re;
                    for(int j=0;j!=line_tsd;++j) {delete [] rs[j];}
                    delete [] rs;
                    for(int j=0;j!=line_tsd;++j) {delete [] len_5[j];}
                    delete [] len_5;
                    for(int j=0;j!=line_tsd;++j) {delete [] len_3[j];}
                    delete [] len_3;
                    
                    for(int j=0;j!=line_tsd;++j) {delete [] kmerseq_tsd[j];}
                    delete [] kmerseq_tsd;
                    
                    for(int w=0;w!=line_tsd;++w){
                        loc_tsd[w][5]=0;
                        loc_tsd[w][4]=0;
                        loc_tsd[w][8]=0;
                    }
                    
                    flag_trans=0;
                    continue;
                }
                
                //find the best subread fit
                int ls_buff, le_buff, rs_buff, re_buff;
                int number_buff;
                
                if(orien[i]=="+"){
                    ls_buff=0;
                    le_buff=0;
                    rs_buff=le_3;
                    re_buff=le_3;
                    number_buff=0;
                }
                if(orien[i]=="-"){
                    ls_buff=le_5;
                    le_buff=le_5;
                    rs_buff=0;
                    re_buff=0;
                    number_buff=0;
                }
                
                //find the most connected-SR and nearest one
                for(int w2=0;w2!=line_tsd;++w2){
//*******
                    if(loc_tsd[w2][4]==1&&loc_tsd[w2][5]==p&&loc_tsd[w2][10]>0&&loc_tsd[w2][11]>0){
                       //&&loc_tsd[w2][10]==0){
                        if(orien[i]=="+"){
                            if((rs_buff+re_buff)>(loc_tsd[w2][2]+loc_tsd[w2][3])){
                                rs_buff=loc_tsd[w2][2];
                                re_buff=loc_tsd[w2][3];
                                number_buff=w2;
                            }
                        }
                        else if(orien[i]=="-"){
                            if((rs_buff+re_buff)<(loc_tsd[w2][2]+loc_tsd[w2][3])){
                                rs_buff=loc_tsd[w2][2];
                                re_buff=loc_tsd[w2][3];
                                number_buff=w2;
                            }
                        }
                    }
                }
                
                
                //***********Hard code detected***********
                //Paramenter value: default length to consider a TransD
                if(orien[i]=="+"&&rs_buff<=50){
                    flag_trans=0;
                }
                else if(orien[i]=="-"&&re_buff>=(BIN_3-50)){
                    flag_trans=0;
                }
                
                if(flag_trans==0){
                    //subreads output
                    int l_tsd=0;
                    int r_tsd=0;
                    int trans_tsd=0;
                    int polyA_le=0;
                    //int flag_singleA=0;
                    
                    for(int w=0;w!=p;++w){
                        int ls_new, le_new, rs_new, re_new;
                        
                        if(orien[i]=="+"){
                            ls_new=ls[number_buff][w]-le_5+len_5[number_buff][w];
                            le_new=le[number_buff][w]-le_5+len_5[number_buff][w];
                            rs_new=rs[number_buff][w];
                            re_new=re[number_buff][w];
                        }
                        else if(orien[i]=="-"){
                            ls_new=ls[number_buff][w];
                            le_new=le[number_buff][w];
                            rs_new=rs[number_buff][w]-le_3+len_3[number_buff][w];
                            re_new=re[number_buff][w]-le_3+len_3[number_buff][w];
                        }
                        
                        //TSD seq
                        int read_seq=-1;
                        string name_tag;
                        for(int j=0;j!=line;++j){
                            
                            string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                            stringstream ss_0;
                            ss_0.clear();
                            ss_0<<loc[j][0];
                            loc_0=ss_0.str();
                            stringstream ss_1;
                            ss_1.clear();
                            ss_1<<loc[j][1];
                            loc_1=ss_1.str();
                            stringstream ss_2;
                            ss_2.clear();
                            ss_2<<loc[j][2];
                            loc_2=ss_2.str();
                            stringstream ss_3;
                            ss_3.clear();
                            ss_3<<loc[j][3];
                            loc_3=ss_3.str();
                            stringstream ss_4;
                            ss_4.clear();
                            ss_4<<loc[j][4];
                            loc_4=ss_4.str();
                            stringstream ss_5;
                            ss_5.clear();
                            ss_5<<loc[j][5];
                            loc_5=ss_5.str();
                            
                            string loc_TP_0, loc_TP_1, loc_TP_2, loc_TP_3, loc_TP_4, loc_TP_5, loc_TP_6;
                            stringstream ss_TP_0;
                            ss_TP_0.clear();
                            ss_TP_0<<loc_TP[j][0];
                            loc_TP_0=ss_TP_0.str();
                            stringstream ss_TP_1;
                            ss_TP_1.clear();
                            ss_TP_1<<loc_TP[j][1];
                            loc_TP_1=ss_TP_1.str();
                            stringstream ss_TP_2;
                            ss_TP_2.clear();
                            ss_TP_2<<loc_TP[j][2];
                            loc_TP_2=ss_TP_2.str();
                            stringstream ss_TP_3;
                            ss_TP_3.clear();
                            ss_TP_3<<loc_TP[j][3];
                            loc_TP_3=ss_TP_3.str();
                            stringstream ss_TP_4;
                            ss_TP_4.clear();
                            ss_TP_4<<loc_TP[j][4];
                            loc_TP_4=ss_TP_4.str();
                            stringstream ss_TP_5;
                            ss_TP_5.clear();
                            ss_TP_5<<loc_TP[j][5];
                            loc_TP_5=ss_TP_5.str();
                            stringstream ss_TP_6;
                            ss_TP_6.clear();
                            ss_TP_6<<loc_TP[j][6];
                            loc_TP_6=ss_TP_6.str();
                            
                            string seq_index,seq_index_2;
                            
                            seq_index=info[j][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                            
                            seq_index_2="";
                            seq_index_2=seq_index_2+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                            
                            if(name_tsd[number_buff][w]==seq_index){
                                read_seq=j;
                                name_tag=seq_index_2;
                            }
                        }
                        if(read_seq==-1){
                            loc_tsd[number_buff][4]=0;
                            continue;
                            
                        }
                        
                        char *seq5 = new char[info[read_seq][2].length()+1];
                        strcpy(seq5, info[read_seq][2].c_str());
                        
                        char *seq3 = new char[info[read_seq][3].length()+1];
                        strcpy(seq3, info[read_seq][3].c_str());
                        
                        char *seq_line = new char[info_line[read_seq][2].length()+1];
                        strcpy(seq_line, info_line[read_seq][2].c_str());
                        
                        file4<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<'\t'<<name_tsd[number_buff][w]<<'\t';
                        
                        if(orien[i]=="+"){
                            for(int n=ls_new-1;n!=le_new;++n){
                                file4<<seq5[n];
                            }
                            file4<<'\t';
                            
                            for(int n=rs_new-1;n!=re_new;++n){
                                file4<<seq3[n];
                            }
                            file4<<'\t';
                            for(int n=0;n!=rs_new-1;++n){
                                file4<<seq3[n];
                            }
                        }
                        
                        else if(orien[i]=="-"){
                            for(int n=ls_new-1;n!=le_new;++n){
                                file4<<seq5[n];
                            }
                            file4<<'\t';
                            
                            for(int n=rs_new-1;n!=re_new;++n){
                                file4<<seq3[n];
                            }
                            file4<<'\t';
                            for(int n=re_new;n!=info[read_seq][3].length();++n){
                                file4<<seq3[n];
                            }
                        }
//5' 26mer output
                        file4<<'\t'<<kmerseq_tsd[number_buff][w]<<'\t';
                        /*
                        string kmer_0, kmer_1, kmer_2, kmer_3;
                        stringstream ss_kmer_0;
                        ss_kmer_0.clear();
                        ss_kmer_0<<ls[number_buff][w];
                        kmer_0=ss_kmer_0.str();
                        
                        stringstream ss_kmer_1;
                        ss_kmer_1.clear();
                        ss_kmer_1<<le[number_buff][w];
                        kmer_1=ss_kmer_1.str();
                        
                        stringstream ss_kmer_2;
                        ss_kmer_2.clear();
                        ss_kmer_2<<rs[number_buff][w];
                        kmer_2=ss_kmer_2.str();
                        
                        stringstream ss_kmer_3;
                        ss_kmer_3.clear();
                        ss_kmer_3<<re[number_buff][w];
                        kmer_3=ss_kmer_3.str();
                        
                        string kmer_index;
                        kmer_index=name_tag+"."+kmer_0.c_str()+"."+kmer_1.c_str()+"."+kmer_2.c_str()+"."+kmer_3.c_str()+".5kmer.fasta";
                        
                        //cout<<kmer_index<<endl;
                        
                        ifstream file26;
                        
                        string sys_kmer_index = WD_dir+kmer_index;
                        char *syst_kmer_index = new char[sys_kmer_index.length()+1];
                        strcpy(syst_kmer_index, sys_kmer_index.c_str());
                        
                        file26.open(syst_kmer_index);
                        
                        string input_kmer;
                        file26>>input_kmer;
                        file26>>input_kmer;
                        
                        file4<<'\t'<<input_kmer<<'\t';
                        
                        file26.close();
                        file26.clear();
                        */
//whole insertin sequence output
                        
                        if(orien[i]=="+"){
                            for(int n=le_new;n!=info[read_seq][2].length();++n){
                                file4<<seq5[n];
                            }
//5bp
                            for(int n=5;n!=info_line[read_seq][2].length()-1-5;++n){
                                file4<<seq_line[n];
                            }
                            
                            for(int n=0;n!=rs_new-1;++n){
                                file4<<seq3[n];
                            }
                        }
                        
                        else if(orien[i]=="-"){
                            
                            for(int n=re_new;n!=info[read_seq][3].length();++n){
                                file4<<seq3[n];
                            }
//5bp
                            for(int n=1+5;n!=info_line[read_seq][2].length()-1-5;++n){
                                file4<<seq_line[n];
                            }
                            
                            for(int n=0;n!=ls_new-1;++n){
                                file4<<seq5[n];
                            }
                            
                        }
                        
                        
                        file4<<endl;
                        if(t=="LINE"){
                            polyA_le=polyA_le+loc[read_seq][1]-6022;
                        }
                        else if (t=="ALU"){
                            polyA_le=polyA_le+loc[read_seq][1]-282;
                        }
                        else if (t=="SVA"){
                            polyA_le=polyA_le+loc[read_seq][1]-1352;
                        }
                        else if (t=="HERVK"){
                            polyA_le=polyA_le+loc[read_seq][1]-9472;
                        }
                        //polyA_le=polyA_le+loc[read_seq][1]-6025;
                        l_tsd=l_tsd+le_new-ls_new;
                        r_tsd=r_tsd+re_new-rs_new;
                        if(orien[i]=="+"){
                            trans_tsd=trans_tsd+rs[number_buff][w];
                        }
                        if(orien[i]=="-"){
                            trans_tsd=trans_tsd+le_3-+re[number_buff][w];
                        }
                        
                        delete [] seq3;
                        delete [] seq5;
                        delete [] seq_line;
                    }
                    
                    //TSD & trans length
                    double l_len=0;
                    double trans_len=0;
                    double r_len=0;
                    double pa_len=0;
                    
                    if(p!=0){
                        l_len=1+l_tsd/p;
                        r_len=1+r_tsd/p;
                        pa_len=polyA_le/p;
                        trans_len=trans_tsd/(p);
                    }
                    
                    
                    //cluster output
                    if(pa_len<0){
                        file2<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<'\t'<<info[i][1]<<'\t'<<start1+S<<'\t'<<start2-S<<'\t'<<end1+S<<'\t'<<end2-S<<'\t'<<L1_s1+L<<'\t'<<L1_s2-L<<'\t'<<L1_e1+L<<'\t'<<L1_e2-L<<'\t'<<p<<'\t'<<number_all<<'\t'<<(number_all-number)<<'\t'<<orien[i]<<'\t'<<"-"<<'\t'<<l_len<<'\t'<<r_len<<'\t'<<trans_len<<'\t'<<loc_TP[i][0]<<'\t'<<int((L1_s_TP1+L1_s_TP2)/2)<<'\t'<<int((L1_e_TP1+L1_e_TP2)/2)<<endl;
                        
                    }
                    else {
                        file2<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<'\t'<<info[i][1]<<'\t'<<start1+S<<'\t'<<start2-S<<'\t'<<end1+S<<'\t'<<end2-S<<'\t'<<L1_s1+L<<'\t'<<L1_s2-L<<'\t'<<L1_e1+L<<'\t'<<L1_e2-L<<'\t'<<p<<'\t'<<number_all<<'\t'<<(number_all-number)<<'\t'<<orien[i]<<'\t'<<pa_len<<'\t'<<l_len<<'\t'<<r_len<<'\t'<<trans_len<<'\t'<<loc_TP[i][0]<<'\t'<<int((L1_s_TP1+L1_s_TP2)/2)<<'\t'<<int((L1_e_TP1+L1_e_TP2)/2)<<endl;
                    }
                    
                    //return
                    for(int w=0;w!=line_tsd;++w){
                        loc_tsd[w][5]=0;
                        loc_tsd[w][4]=0;
                        loc_tsd[w][8]=0;
                    }
                    
                    for(int j=0;j!=line_tsd;++j) {delete [] name_tsd[j];}
                    delete [] name_tsd;
                    
                    for(int j=0;j!=line_tsd;++j) {delete [] ls[j];}
                    delete [] ls;
                    for(int j=0;j!=line_tsd;++j) {delete [] le[j];}
                    delete [] le;
                    for(int j=0;j!=line_tsd;++j) {delete [] re[j];}
                    delete [] re;
                    for(int j=0;j!=line_tsd;++j) {delete [] rs[j];}
                    delete [] rs;
                    for(int j=0;j!=line_tsd;++j) {delete [] len_5[j];}
                    delete [] len_5;
                    for(int j=0;j!=line_tsd;++j) {delete [] len_3[j];}
                    delete [] len_3;
                    
                    for(int j=0;j!=line_tsd;++j) {delete [] kmerseq_tsd[j];}
                    delete [] kmerseq_tsd;
                }
                
                else {
                    //subreads output
                    int l_tsd=0;
                    int r_tsd=0;
                    int trans_tsd=0;
                    int polyA_le=0;
                    int flag_PA=0;
                    
                    // polyA
                    for(int w=0;w!=p;++w){
                        int ls_new, le_new, rs_new, re_new;
                        
                        if(orien[i]=="+"){
                            ls_new=ls[number_buff][w]-le_5+len_5[number_buff][w];
                            le_new=le[number_buff][w]-le_5+len_5[number_buff][w];
                            rs_new=rs[number_buff][w];
                            re_new=re[number_buff][w];
                        }
                        else if(orien[i]=="-"){
                            ls_new=ls[number_buff][w];
                            le_new=le[number_buff][w];
                            rs_new=rs[number_buff][w]-le_3+len_3[number_buff][w];
                            re_new=re[number_buff][w]-le_3+len_3[number_buff][w];
                        }
                        //TSD seq
                        int read_seq=-1;
                        for(int j=0;j!=line;++j){
                            
                            string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                            stringstream ss_0;
                            ss_0.clear();
                            ss_0<<loc[j][0];
                            loc_0=ss_0.str();
                            stringstream ss_1;
                            ss_1.clear();
                            ss_1<<loc[j][1];
                            loc_1=ss_1.str();
                            stringstream ss_2;
                            ss_2.clear();
                            ss_2<<loc[j][2];
                            loc_2=ss_2.str();
                            stringstream ss_3;
                            ss_3.clear();
                            ss_3<<loc[j][3];
                            loc_3=ss_3.str();
                            stringstream ss_4;
                            ss_4.clear();
                            ss_4<<loc[j][4];
                            loc_4=ss_4.str();
                            stringstream ss_5;
                            ss_5.clear();
                            ss_5<<loc[j][5];
                            loc_5=ss_5.str();
                            
                            string loc_TP_0, loc_TP_1, loc_TP_2, loc_TP_3, loc_TP_4, loc_TP_5, loc_TP_6;
                            stringstream ss_TP_0;
                            ss_TP_0.clear();
                            ss_TP_0<<loc_TP[j][0];
                            loc_TP_0=ss_TP_0.str();
                            stringstream ss_TP_1;
                            ss_TP_1.clear();
                            ss_TP_1<<loc_TP[j][1];
                            loc_TP_1=ss_TP_1.str();
                            stringstream ss_TP_2;
                            ss_TP_2.clear();
                            ss_TP_2<<loc_TP[j][2];
                            loc_TP_2=ss_TP_2.str();
                            stringstream ss_TP_3;
                            ss_TP_3.clear();
                            ss_TP_3<<loc_TP[j][3];
                            loc_TP_3=ss_TP_3.str();
                            stringstream ss_TP_4;
                            ss_TP_4.clear();
                            ss_TP_4<<loc_TP[j][4];
                            loc_TP_4=ss_TP_4.str();
                            stringstream ss_TP_5;
                            ss_TP_5.clear();
                            ss_TP_5<<loc_TP[j][5];
                            loc_TP_5=ss_TP_5.str();
                            stringstream ss_TP_6;
                            ss_TP_6.clear();
                            ss_TP_6<<loc_TP[j][6];
                            loc_TP_6=ss_TP_6.str();
                            
                            string seq_index;
                            
                            seq_index=info[j][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                            
                            if(name_tsd[number_buff][w]==seq_index){
                                read_seq=j;
                            }
                        }
                        if(read_seq==-1){
                            loc_tsd[number_buff][4]=0;
                            continue;
                            
                        }
                        
                        char *seq3 = new char[info[read_seq][3].length()+1];
                        strcpy(seq3, info[read_seq][3].c_str());

                        //char *seqA = new char[15];
                        
                        int singleA=0;
                        
                        if(orien[i]=="+"){
                            for(int n=rs_new-10;n!=rs_new;++n){
                                if(seq3[n]=='A'||seq3[n]=='a'){
                                    singleA=singleA+1;
                                }
                            }
                        }
                        
                        else if(orien[i]=="-"){
                            for(int n=re_new;n!=re_new+10;++n){
                                if(seq3[n]=='T'||seq3[n]=='t'){
                                    singleA=singleA+1;
                                }
                            }
                        }
                        
                        if(singleA>=7){
                            flag_PA=flag_PA+1;
                            
                        }
                        
                        delete [] seq3;
                    }
                    
                    loc_tsd[number_buff][4]=0;
                    if(p==1&&flag_PA==1){
                        flag_ployA=1;
                    }
                    else if(p==1&&flag_PA==0){
                        flag_ployA=0;
                        loc_tsd[number_buff][4]=0;
                        continue;
                    }
                    else if(p>1&&flag_PA>1){
                        flag_ployA=1;
                    }
                    else {
                        flag_ployA=0;
                        loc_tsd[number_buff][4]=0;
                        continue;
                    }
                    
                    for(int w=0;w!=p;++w){
                        int ls_new, le_new, rs_new, re_new;
                        
                        if(orien[i]=="+"){
                            ls_new=ls[number_buff][w]-le_5+len_5[number_buff][w];
                            le_new=le[number_buff][w]-le_5+len_5[number_buff][w];
                            rs_new=rs[number_buff][w];
                            re_new=re[number_buff][w];
                        }
                        else if(orien[i]=="-"){
                            ls_new=ls[number_buff][w];
                            le_new=le[number_buff][w];
                            rs_new=rs[number_buff][w]-le_3+len_3[number_buff][w];
                            re_new=re[number_buff][w]-le_3+len_3[number_buff][w];
                        }
                        //TSD seq
                        int read_seq=-1;
                        string name_tag;
                        for(int j=0;j!=line;++j){
                            
                            string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                            stringstream ss_0;
                            ss_0.clear();
                            ss_0<<loc[j][0];
                            loc_0=ss_0.str();
                            stringstream ss_1;
                            ss_1.clear();
                            ss_1<<loc[j][1];
                            loc_1=ss_1.str();
                            stringstream ss_2;
                            ss_2.clear();
                            ss_2<<loc[j][2];
                            loc_2=ss_2.str();
                            stringstream ss_3;
                            ss_3.clear();
                            ss_3<<loc[j][3];
                            loc_3=ss_3.str();
                            stringstream ss_4;
                            ss_4.clear();
                            ss_4<<loc[j][4];
                            loc_4=ss_4.str();
                            stringstream ss_5;
                            ss_5.clear();
                            ss_5<<loc[j][5];
                            loc_5=ss_5.str();
                            
                            string loc_TP_0, loc_TP_1, loc_TP_2, loc_TP_3, loc_TP_4, loc_TP_5, loc_TP_6;
                            stringstream ss_TP_0;
                            ss_TP_0.clear();
                            ss_TP_0<<loc_TP[j][0];
                            loc_TP_0=ss_TP_0.str();
                            stringstream ss_TP_1;
                            ss_TP_1.clear();
                            ss_TP_1<<loc_TP[j][1];
                            loc_TP_1=ss_TP_1.str();
                            stringstream ss_TP_2;
                            ss_TP_2.clear();
                            ss_TP_2<<loc_TP[j][2];
                            loc_TP_2=ss_TP_2.str();
                            stringstream ss_TP_3;
                            ss_TP_3.clear();
                            ss_TP_3<<loc_TP[j][3];
                            loc_TP_3=ss_TP_3.str();
                            stringstream ss_TP_4;
                            ss_TP_4.clear();
                            ss_TP_4<<loc_TP[j][4];
                            loc_TP_4=ss_TP_4.str();
                            stringstream ss_TP_5;
                            ss_TP_5.clear();
                            ss_TP_5<<loc_TP[j][5];
                            loc_TP_5=ss_TP_5.str();
                            stringstream ss_TP_6;
                            ss_TP_6.clear();
                            ss_TP_6<<loc_TP[j][6];
                            loc_TP_6=ss_TP_6.str();
                            
                            string seq_index,seq_index_2;
                            
                            seq_index=info[j][0]+"."+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                            
                            seq_index_2="";
                            seq_index_2=seq_index_2+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[j][1]+"."+orien[j]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str();
                            
                            if(name_tsd[number_buff][w]==seq_index){
                                read_seq=j;
                                name_tag=seq_index_2;
                            }
                        }
                        if(read_seq==-1){
                            loc_tsd[number_buff][4]=0;
                            continue;
                            
                        }
                        
                        char *seq5 = new char[info[read_seq][2].length()+1];
                        strcpy(seq5, info[read_seq][2].c_str());
                        
                        char *seq3 = new char[info[read_seq][3].length()+1];
                        strcpy(seq3, info[read_seq][3].c_str());
                        
                        char *seq_line = new char[info_line[read_seq][2].length()+1];
                        strcpy(seq_line, info_line[read_seq][2].c_str());
                        
                        file4<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<'\t'<<name_tsd[number_buff][w]<<'\t';
                        
                        if(orien[i]=="+"){
                            for(int n=ls_new-1;n!=le_new;++n){
                                file4<<seq5[n];
                            }
                            file4<<'\t';
                            
                            for(int n=rs_new-1;n!=re_new;++n){
                                file4<<seq3[n];
                            }
                            file4<<'\t';
                            for(int n=0;n!=rs_new-1;++n){
                                file4<<seq3[n];
                            }
                        }
                        
                        else if(orien[i]=="-"){
                            for(int n=ls_new-1;n!=le_new;++n){
                                file4<<seq5[n];
                            }
                            file4<<'\t';
                            
                            for(int n=rs_new-1;n!=re_new;++n){
                                file4<<seq3[n];
                            }
                            file4<<'\t';
                            for(int n=re_new;n!=info[read_seq][3].length();++n){
                                file4<<seq3[n];
                            }
                        }
                        
//5' 26mer output
                        file4<<'\t'<<kmerseq_tsd[number_buff][w]<<'\t';
                        
                        /*
                        string kmer_0, kmer_1, kmer_2, kmer_3;
                        stringstream ss_kmer_0;
                        ss_kmer_0.clear();
                        ss_kmer_0<<ls[number_buff][w];
                        kmer_0=ss_kmer_0.str();
                        
                        stringstream ss_kmer_1;
                        ss_kmer_1.clear();
                        ss_kmer_1<<le[number_buff][w];
                        kmer_1=ss_kmer_1.str();
                        
                        stringstream ss_kmer_2;
                        ss_kmer_2.clear();
                        ss_kmer_2<<rs[number_buff][w];
                        kmer_2=ss_kmer_2.str();
                        
                        stringstream ss_kmer_3;
                        ss_kmer_3.clear();
                        ss_kmer_3<<re[number_buff][w];
                        kmer_3=ss_kmer_3.str();
                        
                        string kmer_index;
                        kmer_index=name_tag+"."+kmer_0.c_str()+"."+kmer_1.c_str()+"."+kmer_2.c_str()+"."+kmer_3.c_str()+".5kmer.fasta";
                        
                        //cout<<kmer_index<<endl;
                        
                        ifstream file26;
                        
                        string sys_kmer_index = WD_dir+kmer_index;
                        char *syst_kmer_index = new char[sys_kmer_index.length()+1];
                        strcpy(syst_kmer_index, sys_kmer_index.c_str());
                        
                        file26.open(syst_kmer_index);
                        
                        string input_kmer;
                        file26>>input_kmer;
                        file26>>input_kmer;
                        
                        file4<<'\t'<<input_kmer<<'\t';
                        
                        file26.close();
                        file26.clear();
                        */
//whole insertin sequence output
                        
                        if(orien[i]=="+"){
                            for(int n=le_new;n!=info[read_seq][2].length();++n){
                                file4<<seq5[n];
                            }
                            
                            for(int n=0;n!=info_line[read_seq][2].length()-1;++n){
                                file4<<seq_line[n];
                            }
                            
                            for(int n=0;n!=rs_new-1;++n){
                                file4<<seq3[n];
                            }
                        }
                        
                        else if(orien[i]=="-"){
                            
                            for(int n=re_new;n!=info[read_seq][3].length();++n){
                                file4<<seq3[n];
                            }
                            
                            for(int n=1;n!=info_line[read_seq][2].length()-1;++n){
                                file4<<seq_line[n];
                            }
                            
                            for(int n=0;n!=ls_new-1;++n){
                                file4<<seq5[n];
                            }
                            
                        }
                        
                        
                        file4<<endl;
                        if(t=="LINE"){
                            polyA_le=polyA_le+loc[read_seq][1]-6022;
                        }
                        else if (t=="ALU"){
                            polyA_le=polyA_le+loc[read_seq][1]-282;
                        }
                        else if (t=="SVA"){
                            polyA_le=polyA_le+loc[read_seq][1]-1352;
                        }
                        else if (t=="HERVK"){
                            polyA_le=polyA_le+loc[read_seq][1]-9472;
                        }
                        l_tsd=l_tsd+le_new-ls_new;
                        r_tsd=r_tsd+re_new-rs_new;
                        if(orien[i]=="+"){
                            trans_tsd=trans_tsd+rs[number_buff][w];
                        }
                        if(orien[i]=="-"){
                            trans_tsd=trans_tsd+le_3-re[number_buff][w];
                        }
                        
                        delete [] seq3;
                        delete [] seq5;
                        delete [] seq_line;
                    }
                    
                    //TSD & trans length
                    double l_len=0;
                    double trans_len=0;
                    double r_len=0;
                    double pa_len=0;
                    
                    if(p!=0){
                        l_len=1+l_tsd/p;
                        r_len=1+r_tsd/p;
                        pa_len=polyA_le/p;
                        trans_len=trans_tsd/(p);
                    }
                    
                    
                    //cluster output
                    if(pa_len<0){
                        file2<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<'\t'<<info[i][1]<<'\t'<<start1+S<<'\t'<<start2-S<<'\t'<<end1+S<<'\t'<<end2-S<<'\t'<<L1_s1+L<<'\t'<<L1_s2-L<<'\t'<<L1_e1+L<<'\t'<<L1_e2-L<<'\t'<<p<<'\t'<<number_all<<'\t'<<(number_all-number)<<'\t'<<orien[i]<<'\t'<<"-"<<'\t'<<l_len<<'\t'<<r_len<<'\t'<<trans_len<<'\t'<<loc_TP[i][0]<<'\t'<<int((L1_s_TP1+L1_s_TP2)/2)<<'\t'<<int((L1_e_TP1+L1_e_TP2)/2)<<endl;
                        
                    }
                    else {
                        file2<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<'\t'<<info[i][1]<<'\t'<<start1+S<<'\t'<<start2-S<<'\t'<<end1+S<<'\t'<<end2-S<<'\t'<<L1_s1+L<<'\t'<<L1_s2-L<<'\t'<<L1_e1+L<<'\t'<<L1_e2-L<<'\t'<<p<<'\t'<<number_all<<'\t'<<(number_all-number)<<'\t'<<orien[i]<<'\t'<<pa_len<<'\t'<<l_len<<'\t'<<r_len<<'\t'<<trans_len<<'\t'<<loc_TP[i][0]<<'\t'<<int((L1_s_TP1+L1_s_TP2)/2)<<'\t'<<int((L1_e_TP1+L1_e_TP2)/2)<<endl;
                    }
                    
                    //return
                    for(int w=0;w!=line_tsd;++w){
                        loc_tsd[w][5]=0;
                        loc_tsd[w][4]=0;
                        loc_tsd[w][8]=0;
                    }
                    
                    for(int j=0;j!=line_tsd;++j) {delete [] name_tsd[j];}
                    delete [] name_tsd;
                    
                    for(int j=0;j!=line_tsd;++j) {delete [] ls[j];}
                    delete [] ls;
                    for(int j=0;j!=line_tsd;++j) {delete [] le[j];}
                    delete [] le;
                    for(int j=0;j!=line_tsd;++j) {delete [] re[j];}
                    delete [] re;
                    for(int j=0;j!=line_tsd;++j) {delete [] rs[j];}
                    delete [] rs;
                    for(int j=0;j!=line_tsd;++j) {delete [] len_5[j];}
                    delete [] len_5;
                    for(int j=0;j!=line_tsd;++j) {delete [] len_3[j];}
                    delete [] len_3;
                    
                    for(int j=0;j!=line_tsd;++j) {delete [] kmerseq_tsd[j];}
                    delete [] kmerseq_tsd;
                    
                }
            }
            }
            else if (tsd_index==0){
                file2<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<'\t'<<info[i][1]<<'\t'<<start1+S<<'\t'<<start2-S<<'\t'<<end1+S<<'\t'<<end2-S<<'\t'<<L1_s1+L<<'\t'<<L1_s2-L<<'\t'<<L1_e1+L<<'\t'<<L1_e2-L<<'\t'<<"0"<<'\t'<<number_all<<'\t'<<(number_all-number)<<'\t'<<orien[i]<<'\t'<<"-"<<'\t'<<"0"<<'\t'<<"0"<<'\t'<<"0"<<'\t'<<loc_TP[i][0]<<'\t'<<int((L1_s_TP1+L1_s_TP2)/2)<<'\t'<<int((L1_e_TP1+L1_e_TP2)/2)<<endl;
                file4<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<'\t'<<seq_index_a<<'\t'<<"N"<<'\t'<<"N"<<'\t'<<"N"<<'\t'<<"N"<<'\t'<<info_line[i][2]<<endl;
            }
                
        }
    }
    
    for(int i=0;i!=line;++i){
        delete [] info[i];
        delete [] loc[i];
        delete [] loc_TP[i];
        delete [] loc_TP_line[i];
        delete [] info_line[i];
    }
    delete [] info;
    delete [] loc;
    delete [] loc_TP_line;
    delete [] info_line;
    delete [] loc_TP;
    
    for(int i=0;i!=line_tsd;++i){
        delete [] loc_tsd[i];
        //delete [] loc[i];
    }
    delete [] loc_tsd;
    delete [] orien;
    delete [] info_tsd;
    delete [] kmer_tsd;
    
    return 0;
}

