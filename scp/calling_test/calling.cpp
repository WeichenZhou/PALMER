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

using namespace std;

//int calling(string WD_dir){
//test
int main(){
    string WD_dir="/home/arthurz/arthur_remflux_scratch/16.12.30.MEI_caller_forPB/v7.5_170407/scp/calling_test/";
    
    ifstream file1;
    
    string sys_input = WD_dir+"read_result_TSD.txt";
    char *syst_input = new char[sys_input.length()+1];
    strcpy(syst_input, sys_input.c_str());
    
    file1.open(syst_input);
    
    ofstream file2;
    
    string sys_calls = WD_dir+"calls.txt";
    char *syst_calls = new char[sys_calls.length()+1];
    strcpy(syst_calls, sys_calls.c_str());
    file2.open(syst_calls);
    
    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE"<< endl;
        exit(1);
    }
    
    int line;
    string input;
    for(int i=0;!file1.eof();i++){
        getline(file1,input);
        line=i;
    }
    
    file1.close();
    file1.clear();
    file1.open(syst_input);
    
    string **info;
    info=new string *[line];
    for(int i=0;i!=line;i++) info[i]=new string[4];
    
    int **loc;
    loc=new int*[line];
    for(int i=0;i!=line;i++) loc[i]=new int[7];
    
    string *orien;
    orien= new string[line];
    
    for(int i=0;i!=line;i++){
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
    }
    
    //TSD_calling
    ifstream file3;
    string sys_input_TSD = WD_dir+"TSD_blastn.txt";
    char *syst_input_TSD = new char[sys_input_TSD.length()+1];
    strcpy(syst_input_TSD, sys_input_TSD.c_str());
    file3.open(syst_input_TSD);
    
    if (!file3.is_open())
    {
        cout <<"CANNOT OPEN FILE"<< endl;
        exit(1);
    }
    
    ofstream file4;
    string sys_output_TSD = WD_dir+"TSD_output.txt";
    char *syst_output_TSD = new char[sys_output_TSD.length()+1];
    strcpy(syst_output_TSD, sys_output_TSD.c_str());
    file4.open(syst_output_TSD);
    
    int line_tsd;
    for(int i=0;!file3.eof();i++){
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
    for(int i=0;i!=line_tsd;i++) loc_tsd[i]=new int[9];
    
    for(int i=0;i!=line_tsd;i++){
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
    }
    
    string sys_readtag = WD_dir+"readtag";
    char *syst_readtag = new char[sys_readtag.length()+1];
    strcpy(syst_readtag, sys_readtag.c_str());
    
    ofstream file12;
    ifstream file13;
    
    
    //calling
    const int S=150;
    const int L=50;
    
    for(int i=0;i!=line;i++){
        if(loc[i][6]!=-1){
            //cout<<"OK"<<endl;
            //getchar();
            
            int start1=loc[i][4]-S;
            int start2=loc[i][4]+S;
            int end2=loc[i][5]+S;
            int end1=loc[i][5]-S;
            int L1_s1=loc[i][0]-L;
            int L1_s2=loc[i][0]+L;
            int L1_e1=loc[i][1]-L;
            int L1_e2=loc[i][1]+L;
            
            int number=0;
            int number_all=1;
            loc[i][6]=-1;
            int flag=1;
            
            file12.open(syst_readtag);
            file12<<info[i][0]<<"_"<<loc[i][0]<<"_"<<loc[i][1]<<"_"<<loc[i][2]<<"_"<<loc[i][3]<<"_"<<loc[i][4]<<"_"<<loc[i][5]<<"_"<<info[i][1]<<"_"<<orien[i]<<endl;
            
            string seq_index;
            file13.open(syst_readtag);
            getline(file13,seq_index);
            
            file12.close();
            file12.clear();
            file13.close();
            file13.clear();
            
            for(int w=0;w!=line_tsd;w++){
                if(info_tsd[w]==seq_index){
                    loc_tsd[w][8]=1;
                }
            }
            
            int flag_left=0;
            int flag_right=0;
            
            for(;flag==1;){
                flag=0;
                
                for(int j=0;j!=line;j++){

                    if(info[i][0]==info[j][0]&&loc[j][6]==0&&info[i][1]==info[j][1]&&loc[i][0]==loc[j][0]&&loc[i][1]==loc[j][1]&&loc[i][2]==loc[j][2]&&loc[i][3]==loc[j][3]&&loc[i][4]==loc[j][4]&&loc[i][5]==loc[j][5]&&orien[i]==orien[j]){
                        loc[j][6]=-1;
                        
                    }
                    else if(info[i][0]!=info[j][0]&&loc[j][6]==0&&info[i][1]==info[j][1]&&orien[i]==orien[j]){
                        if((loc[j][4]>=start1&&loc[j][4]<=start2)||(loc[j][5]<=end2&&loc[j][5]>=end1)){
                            
                            //left&right
                            if((loc[j][0]>=L1_s1&&loc[j][0]<=L1_s2)&&(loc[j][1]>=L1_e1&&loc[j][1]<=L1_e2)){
                                if(loc[j][4]-S<start1) start1=loc[j][4]-S;
                                if(loc[j][4]+S>start2) start2=loc[j][4]+S;
                                if(loc[j][5]-S<end1) end1=loc[j][5]-S;
                                if(loc[j][5]+S>end2) end2=loc[j][5]+S;
                                
                                if(loc[j][0]-L<L1_s1) L1_s1=loc[j][0]-L;
                                if(loc[j][0]+L>L1_s2) L1_s2=loc[j][0]+L;
                                if(loc[j][1]-L<L1_e1) L1_e1=loc[j][1]-L;
                                if(loc[j][1]+L>L1_e2) L1_e2=loc[j][1]+L;
                                
                                flag=1;
                                loc[j][6]=-1;
                                number_all=number_all+1;
                                
                                file12.open(syst_readtag);
                                file12<<info[j][0]<<"_"<<loc[j][0]<<"_"<<loc[j][1]<<"_"<<loc[j][2]<<"_"<<loc[j][3]<<"_"<<loc[j][4]<<"_"<<loc[j][5]<<"_"<<info[j][1]<<"_"<<orien[j]<<endl;
                                
                                string seq_index;
                                file13.open(syst_readtag);
                                getline(file13,seq_index);
                                
                                file12.close();
                                file12.clear();
                                file13.close();
                                file13.clear();
                                
                                for(int w=0;w!=line_tsd;w++){
                                    if(info_tsd[w]==seq_index){
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
                                    
                                    file12.open(syst_readtag);
                                    file12<<info[j][0]<<"_"<<loc[j][0]<<"_"<<loc[j][1]<<"_"<<loc[j][2]<<"_"<<loc[j][3]<<"_"<<loc[j][4]<<"_"<<loc[j][5]<<"_"<<info[j][1]<<"_"<<orien[j]<<endl;
                                    
                                    string seq_index;
                                    file13.open(syst_readtag);
                                    getline(file13,seq_index);
                                    
                                    file12.close();
                                    file12.clear();
                                    file13.close();
                                    file13.clear();
                                    
                                    for(int w=0;w!=line_tsd;w++){
                                        if(info_tsd[w]==seq_index){
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
                                    
                                    file12.open(syst_readtag);
                                    file12<<info[j][0]<<"_"<<loc[j][0]<<"_"<<loc[j][1]<<"_"<<loc[j][2]<<"_"<<loc[j][3]<<"_"<<loc[j][4]<<"_"<<loc[j][5]<<"_"<<info[j][1]<<"_"<<orien[j]<<endl;
                                    
                                    string seq_index;
                                    file13.open(syst_readtag);
                                    getline(file13,seq_index);
                                    
                                    file12.close();
                                    file12.clear();
                                    file13.close();
                                    file13.clear();
                                    
                                    for(int w=0;w!=line_tsd;w++){
                                        if(info_tsd[w]==seq_index){
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
                                    
                                    file12.open(syst_readtag);
                                    file12<<info[j][0]<<"_"<<loc[j][0]<<"_"<<loc[j][1]<<"_"<<loc[j][2]<<"_"<<loc[j][3]<<"_"<<loc[j][4]<<"_"<<loc[j][5]<<"_"<<info[j][1]<<"_"<<orien[j]<<endl;
                                    
                                    string seq_index;
                                    file13.open(syst_readtag);
                                    getline(file13,seq_index);
                                    
                                    file12.close();
                                    file12.clear();
                                    file13.close();
                                    file13.clear();
                                    
                                    for(int w=0;w!=line_tsd;w++){
                                        if(info_tsd[w]==seq_index){
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
                                    
                                    file12.open(syst_readtag);
                                    file12<<info[j][0]<<"_"<<loc[j][0]<<"_"<<loc[j][1]<<"_"<<loc[j][2]<<"_"<<loc[j][3]<<"_"<<loc[j][4]<<"_"<<loc[j][5]<<"_"<<info[j][1]<<"_"<<orien[j]<<endl;
                                    
                                    string seq_index;
                                    file13.open(syst_readtag);
                                    getline(file13,seq_index);
                                    
                                    file12.close();
                                    file12.clear();
                                    file13.close();
                                    file13.clear();
                                    
                                    for(int w=0;w!=line_tsd;w++){
                                        if(info_tsd[w]==seq_index){
                                            loc_tsd[w][8]=1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            for(int j=0;j!=line;j++){
                if((loc[j][0]>=L1_s1&&loc[j][0]<=L1_s2)&&(loc[j][1]>=L1_e1&&loc[j][1]<=L1_e2)&&(loc[j][4]>=start1&&loc[j][4]<=start2)&&(loc[j][5]>=end1&&loc[j][5]<=end2)&&loc[j][6]==-1&&info[i][1]==info[j][1]&&orien[i]==orien[j]){
                    file12.open(syst_readtag);
                    file12<<info[j][0]<<"_"<<loc[j][0]<<"_"<<loc[j][1]<<"_"<<loc[j][2]<<"_"<<loc[j][3]<<"_"<<loc[j][4]<<"_"<<loc[j][5]<<"_"<<info[j][1]<<"_"<<orien[j]<<endl;
                    
                    string seq_index;
                    file13.open(syst_readtag);
                    getline(file13,seq_index);
                    
                    file12.close();
                    file12.clear();
                    file13.close();
                    file13.clear();
                    
                    int flag_number=0;
                    for(int w=0;w!=line_tsd;w++){
                        if(info_tsd[w]==seq_index&&loc_tsd[w][8]==1){
                            loc_tsd[w][4]=1;
                            flag_number=1;
                            
                        }
                    }
                    if(flag_number==1){
                        number=number+1;
                    }
                }
            }
            //TSD_indentifier
            
            string **name_tsd;
            name_tsd =new string*[line_tsd];
            for(int j=0;j!=line_tsd;j++) name_tsd[j]= new string[number];
            
            int **ls;
            ls =new int*[line_tsd];
            for(int j=0;j!=line_tsd;j++) ls[j]= new int[number];
            int **le;
            le =new int*[line_tsd];
            for(int j=0;j!=line_tsd;j++) le[j]= new int[number];
            int **rs;
            rs =new int*[line_tsd];
            for(int j=0;j!=line_tsd;j++) rs[j]= new int[number];
            int **re;
            re =new int*[line_tsd];
            for(int j=0;j!=line_tsd;j++) re[j]= new int[number];
            int **len_5;
            len_5 =new int*[line_tsd];
            for(int j=0;j!=line_tsd;j++) len_5[j]= new int[number];
            int **len_3;
            len_3 =new int*[line_tsd];
            for(int j=0;j!=line_tsd;j++) len_3[j]= new int[number];
            
            //loc_tsd[][5]
            for(int w=0;w!=line_tsd;w++){
                if(loc_tsd[w][4]==1){
                    loc_tsd[w][5]=1;
                    
                    name_tsd[w][loc_tsd[w][5]-1]=info_tsd[w];
                    ls[w][loc_tsd[w][5]-1]=loc_tsd[w][0];
                    le[w][loc_tsd[w][5]-1]=loc_tsd[w][1];
                    rs[w][loc_tsd[w][5]-1]=loc_tsd[w][2];
                    re[w][loc_tsd[w][5]-1]=loc_tsd[w][3];
                    len_5[w][loc_tsd[w][5]-1]=loc_tsd[w][6];
                    len_3[w][loc_tsd[w][5]-1]=loc_tsd[w][7];
                    
                    for(int w2=0;w2!=line_tsd;w2++){
                        
                        int flag_name=0;
                        for(int w_name=0;w_name!=loc_tsd[w][5];w_name++){
                            if(name_tsd[w][w_name]==info_tsd[w2]){
                                flag_name=1;
                                break;
                            }
                        }
                        if(loc_tsd[w2][4]==1&&flag_name==0){
                            
                            if((!((loc_tsd[w2][0]>loc_tsd[w][1])||(loc_tsd[w2][1]<loc_tsd[w][0])))&&(!((loc_tsd[w2][2]>loc_tsd[w][3])||(loc_tsd[w2][3]<loc_tsd[w][2])))){
                                int ls_tsd=loc_tsd[w][0];
                                int le_tsd=loc_tsd[w][1];
                                int rs_tsd=loc_tsd[w][2];
                                int re_tsd=loc_tsd[w][3];
                                
                                if(loc_tsd[w2][0]>ls_tsd) ls_tsd=loc_tsd[w2][0];
                                if(loc_tsd[w2][2]>rs_tsd) rs_tsd=loc_tsd[w2][2];
                                if(loc_tsd[w2][1]<le_tsd) le_tsd=loc_tsd[w2][1];
                                if(loc_tsd[w2][3]<re_tsd) re_tsd=loc_tsd[w2][3];
                                
                                if(((le_tsd-ls_tsd)*2>(loc_tsd[w][1]-loc_tsd[w][0]))&&((re_tsd-rs_tsd)*2>(loc_tsd[w][3]-loc_tsd[w][2]))){
                                    loc_tsd[w][5]=1+loc_tsd[w][5];
                                    name_tsd[w][loc_tsd[w][5]-1]=info_tsd[w2];
                                    ls[w][loc_tsd[w][5]-1]=loc_tsd[w2][0];
                                    le[w][loc_tsd[w][5]-1]=loc_tsd[w2][1];
                                    rs[w][loc_tsd[w][5]-1]=loc_tsd[w2][2];
                                    re[w][loc_tsd[w][5]-1]=loc_tsd[w2][3];
                                    len_5[w][loc_tsd[w][5]-1]=loc_tsd[w2][6];
                                    len_3[w][loc_tsd[w][5]-1]=loc_tsd[w2][7];
                                    
                                }
                            }
                        }
                    }
                }
            }
            
            
            
            //transD ployA
            int flag_ployA=0;
            int flag_trans=1;
            
            for(;flag_ployA==0&&flag_trans==1;){
                
                //best subreads cluster candidate
                int p=0;
                for(int w=0;w!=line_tsd;w++){
                    if(loc_tsd[w][4]==1){
                        if(loc_tsd[w][5]>p) p=loc_tsd[w][5];
                    }
                }
                
                //cout<<"p="<<p<<endl;
                
                if(p==0){
                    for(int j=0;j!=line_tsd;j++) {delete [] name_tsd[j];}
                    delete [] name_tsd;
                    
                    for(int j=0;j!=line_tsd;j++) {delete [] ls[j];}
                    delete [] ls;
                    for(int j=0;j!=line_tsd;j++) {delete [] le[j];}
                    delete [] le;
                    for(int j=0;j!=line_tsd;j++) {delete [] re[j];}
                    delete [] re;
                    for(int j=0;j!=line_tsd;j++) {delete [] rs[j];}
                    delete [] rs;
                    for(int j=0;j!=line_tsd;j++) {delete [] len_5[j];}
                    delete [] len_5;
                    for(int j=0;j!=line_tsd;j++) {delete [] len_3[j];}
                    delete [] len_3;
                    for(int w=0;w!=line_tsd;w++){
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
                    rs_buff=3001;
                    re_buff=3001;
                    number_buff=0;
                }
                if(orien[i]=="-"){
                    ls_buff=51;
                    le_buff=51;
                    rs_buff=0;
                    re_buff=0;
                    number_buff=0;
                }
                
                //find the most connected-SR and nearest one
                for(int w2=0;w2!=line_tsd;w2++){
                    if(loc_tsd[w2][4]==1&&loc_tsd[w2][5]==p){
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
                
                if(orien[i]=="+"&&re_buff<=50){
                    flag_trans=0;
                }
                else if(orien[i]=="-"&&rs_buff>=2950){
                    flag_trans=0;
                }
                
                if(flag_trans==0){
                    //subreads output
                    int l_tsd=0;
                    int r_tsd=0;
                    int trans_tsd=0;
                    
                    //int flag_singleA=0;
                    
                    for(int w=0;w!=p;w++){
                        int ls_new, le_new, rs_new, re_new;
                        
                        if(orien[i]=="+"){
                            ls_new=ls[number_buff][w]-51+len_5[number_buff][w];
                            le_new=le[number_buff][w]-51+len_5[number_buff][w];
                            rs_new=rs[number_buff][w];
                            re_new=re[number_buff][w];
                        }
                        else if(orien[i]=="-"){
                            ls_new=ls[number_buff][w];
                            le_new=le[number_buff][w];
                            rs_new=rs[number_buff][w]-3001+len_3[number_buff][w];
                            re_new=re[number_buff][w]-3001+len_3[number_buff][w];
                        }
                        
                        //TSD seq
                        int read_seq=-1;
                        for(int j=0;j!=line;j++){
                            file12.open(syst_readtag);
                            file12<<info[j][0]<<"_"<<loc[j][0]<<"_"<<loc[j][1]<<"_"<<loc[j][2]<<"_"<<loc[j][3]<<"_"<<loc[j][4]<<"_"<<loc[j][5]<<"_"<<info[j][1]<<"_"<<orien[j]<<endl;
                            
                            string seq_index;
                            file13.open(syst_readtag);
                            getline(file13,seq_index);
                            
                            file12.close();
                            file12.clear();
                            file13.close();
                            file13.clear();
                            if(name_tsd[number_buff][w]==seq_index){
                                read_seq=j;
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
                        
                        file4<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<"cluster"<<i<<'\t'<<name_tsd[number_buff][w]<<'\t'<<ls_new<<'\t'<<le_new<<'\t'<<rs_new<<'\t'<<re_new<<'\t';
                        
                        if(orien[i]=="+"){
                            for(int n=ls_new-1;n!=le_new;n++){
                                file4<<seq5[n];
                            }
                            file4<<'\t';
                            
                            for(int n=rs_new-1;n!=re_new;n++){
                                file4<<seq3[n];
                            }
                            file4<<'\t';
                            for(int n=0;n!=rs_new;n++){
                                file4<<seq3[n];
                            }
                        }
                        
                        else if(orien[i]=="-"){
                            for(int n=ls_new-1;n!=le_new;n++){
                                file4<<seq5[n];
                            }
                            file4<<'\t';
                            
                            for(int n=rs_new-1;n!=re_new;n++){
                                file4<<seq3[n];
                            }
                            file4<<'\t';
                            for(int n=re_new;n!=info[read_seq][3].length();n++){
                                file4<<seq3[n];
                            }
                        }
                        
                        file4<<endl;
                        
                        l_tsd=l_tsd+le_new-ls_new;
                        r_tsd=r_tsd+re_new-rs_new;
                        if(orien[i]=="+"){
                            trans_tsd=trans_tsd+rs[number_buff][w];
                        }
                        if(orien[i]=="-"){
                            trans_tsd=trans_tsd+3001-+re[number_buff][w];
                        }
                        
                        delete [] seq3;
                        delete [] seq5;
                    }
                    
                    //TSD & trans length
                    double l_len=0;
                    double trans_len=0;
                    double r_len=0;
                    
                    if(p!=0){
                        l_len=1+l_tsd/p;
                        r_len=1+r_tsd/p;
                        
                        trans_len=trans_tsd/(p);
                    }
                    
                    
                    //cluster output
                    file2<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<"cluster"<<i<<'\t'<<info[i][1]<<'\t'<<start1+S<<'\t'<<start2-S<<'\t'<<end1+S<<'\t'<<end2-S<<'\t'<<L1_s1+L<<'\t'<<L1_s2-L<<'\t'<<L1_e1+L<<'\t'<<L1_e2-L<<'\t'<<number_all<<'\t'<<p<<'\t'<<orien[i]<<'\t'<<l_len<<'\t'<<r_len<<'\t'<<trans_len<<endl;
                    
                    //return
                    for(int w=0;w!=line_tsd;w++){
                        loc_tsd[w][5]=0;
                        loc_tsd[w][4]=0;
                        loc_tsd[w][8]=0;
                    }
                    
                    for(int j=0;j!=line_tsd;j++) {delete [] name_tsd[j];}
                    delete [] name_tsd;
                    
                    for(int j=0;j!=line_tsd;j++) {delete [] ls[j];}
                    delete [] ls;
                    for(int j=0;j!=line_tsd;j++) {delete [] le[j];}
                    delete [] le;
                    for(int j=0;j!=line_tsd;j++) {delete [] re[j];}
                    delete [] re;
                    for(int j=0;j!=line_tsd;j++) {delete [] rs[j];}
                    delete [] rs;
                    for(int j=0;j!=line_tsd;j++) {delete [] len_5[j];}
                    delete [] len_5;
                    for(int j=0;j!=line_tsd;j++) {delete [] len_3[j];}
                    delete [] len_3;
                }
                
                else {
                    //subreads output
                    int l_tsd=0;
                    int r_tsd=0;
                    int trans_tsd=0;
                    
                    int flag_PA=0;
                    
                    // polyA
                    for(int w=0;w!=p;w++){
                        int ls_new, le_new, rs_new, re_new;
                        
                        if(orien[i]=="+"){
                            ls_new=ls[number_buff][w]-51+len_5[number_buff][w];
                            le_new=le[number_buff][w]-51+len_5[number_buff][w];
                            rs_new=rs[number_buff][w];
                            re_new=re[number_buff][w];
                        }
                        else if(orien[i]=="-"){
                            ls_new=ls[number_buff][w];
                            le_new=le[number_buff][w];
                            rs_new=rs[number_buff][w]-3001+len_3[number_buff][w];
                            re_new=re[number_buff][w]-3001+len_3[number_buff][w];
                        }
                        //TSD seq
                        int read_seq=-1;
                        for(int j=0;j!=line;j++){
                            file12.open(syst_readtag);
                            file12<<info[j][0]<<"_"<<loc[j][0]<<"_"<<loc[j][1]<<"_"<<loc[j][2]<<"_"<<loc[j][3]<<"_"<<loc[j][4]<<"_"<<loc[j][5]<<"_"<<info[j][1]<<"_"<<orien[j]<<endl;
                            
                            string seq_index;
                            file13.open(syst_readtag);
                            getline(file13,seq_index);
                            
                            file12.close();
                            file12.clear();
                            file13.close();
                            file13.clear();
                            if(name_tsd[number_buff][w]==seq_index){
                                read_seq=j;
                            }
                        }
                        if(read_seq==-1){
                            loc_tsd[number_buff][4]=0;
                            continue;
                            
                        }
                        
                        char *seqA = new char[15];
                        
                        int singleA=0;
                        
                        if(orien[i]=="+"){
                            for(int n=rs_new-15;n!=rs_new;n++){
                                if(seqA[n]=='A'||seqA[n]=='a'){
                                    singleA=singleA+1;
                                }
                            }
                        }
                        
                        else if(orien[i]=="-"){
                            for(int n=re_new;n!=re_new+15;n++){
                                if(seqA[n]=='T'||seqA[n]=='t'){
                                    singleA=singleA+1;
                                }
                            }
                        }
                        
                        if(singleA>7){
                            flag_PA=flag_PA+1;
                            
                        }
                        
                        delete [] seqA;
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
                    
                    for(int w=0;w!=p;w++){
                        int ls_new, le_new, rs_new, re_new;
                        
                        if(orien[i]=="+"){
                            ls_new=ls[number_buff][w]-51+len_5[number_buff][w];
                            le_new=le[number_buff][w]-51+len_5[number_buff][w];
                            rs_new=rs[number_buff][w];
                            re_new=re[number_buff][w];
                        }
                        else if(orien[i]=="-"){
                            ls_new=ls[number_buff][w];
                            le_new=le[number_buff][w];
                            rs_new=rs[number_buff][w]-3001+len_3[number_buff][w];
                            re_new=re[number_buff][w]-3001+len_3[number_buff][w];
                        }
                        //TSD seq
                        int read_seq=-1;
                        for(int j=0;j!=line;j++){
                            file12.open(syst_readtag);
                            file12<<info[j][0]<<"_"<<loc[j][0]<<"_"<<loc[j][1]<<"_"<<loc[j][2]<<"_"<<loc[j][3]<<"_"<<loc[j][4]<<"_"<<loc[j][5]<<"_"<<info[j][1]<<"_"<<orien[j]<<endl;
                            
                            string seq_index;
                            file13.open(syst_readtag);
                            getline(file13,seq_index);
                            
                            file12.close();
                            file12.clear();
                            file13.close();
                            file13.clear();
                            if(name_tsd[number_buff][w]==seq_index){
                                read_seq=j;
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
                        
                        file4<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<"cluster"<<i<<'\t'<<name_tsd[number_buff][w]<<'\t'<<ls_new<<'\t'<<le_new<<'\t'<<rs_new<<'\t'<<re_new<<'\t';
                        
                        if(orien[i]=="+"){
                            for(int n=ls_new-1;n!=le_new;n++){
                                file4<<seq5[n];
                            }
                            file4<<'\t';
                            
                            for(int n=rs_new-1;n!=re_new;n++){
                                file4<<seq3[n];
                            }
                            file4<<'\t';
                            for(int n=0;n!=rs_new;n++){
                                file4<<seq3[n];
                            }
                        }
                        
                        else if(orien[i]=="-"){
                            for(int n=ls_new-1;n!=le_new;n++){
                                file4<<seq5[n];
                            }
                            file4<<'\t';
                            
                            for(int n=rs_new-1;n!=re_new;n++){
                                file4<<seq3[n];
                            }
                            file4<<'\t';
                            for(int n=re_new;n!=info[read_seq][3].length();n++){
                                file4<<seq3[n];
                            }
                        }
                        
                        file4<<endl;
                        
                        l_tsd=l_tsd+le_new-ls_new;
                        r_tsd=r_tsd+re_new-rs_new;
                        if(orien[i]=="+"){
                            trans_tsd=trans_tsd+rs[number_buff][w];
                        }
                        if(orien[i]=="-"){
                            trans_tsd=trans_tsd+3001-re[number_buff][w];
                        }
                        
                        delete [] seq3;
                        delete [] seq5;
                    }
                    
                    //TSD & trans length
                    double l_len=0;
                    double trans_len=0;
                    double r_len=0;
                    
                    if(p!=0){
                        l_len=1+l_tsd/p;
                        r_len=1+r_tsd/p;
                        
                        trans_len=trans_tsd/(p);
                    }
                    
                    
                    //cluster output
                    file2<<"cluster"<<i<<"_"<<info[i][1]<<"_"<<start1+S<<"_"<<start2-S<<"_"<<end1+S<<"_"<<end2-S<<"cluster"<<i<<'\t'<<info[i][1]<<'\t'<<start1+S<<'\t'<<start2-S<<'\t'<<end1+S<<'\t'<<end2-S<<'\t'<<L1_s1+L<<'\t'<<L1_s2-L<<'\t'<<L1_e1+L<<'\t'<<L1_e2-L<<'\t'<<number<<'\t'<<p<<'\t'<<orien[i]<<'\t'<<l_len<<'\t'<<r_len<<'\t'<<trans_len<<endl;
                    
                    //return
                    for(int w=0;w!=line_tsd;w++){
                        loc_tsd[w][5]=0;
                        loc_tsd[w][4]=0;
                        loc_tsd[w][8]=0;
                    }
                    
                    for(int j=0;j!=line_tsd;j++) {delete [] name_tsd[j];}
                    delete [] name_tsd;
                    
                    for(int j=0;j!=line_tsd;j++) {delete [] ls[j];}
                    delete [] ls;
                    for(int j=0;j!=line_tsd;j++) {delete [] le[j];}
                    delete [] le;
                    for(int j=0;j!=line_tsd;j++) {delete [] re[j];}
                    delete [] re;
                    for(int j=0;j!=line_tsd;j++) {delete [] rs[j];}
                    delete [] rs;
                    for(int j=0;j!=line_tsd;j++) {delete [] len_5[j];}
                    delete [] len_5;
                    for(int j=0;j!=line_tsd;j++) {delete [] len_3[j];}
                    delete [] len_3;
                    
                }                
            }
        }
    }
}

