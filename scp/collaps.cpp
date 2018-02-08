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

int collaps(string WD_dir){
    
    ifstream file1;
    ofstream file5;
    
    string sys_calls = WD_dir+"calls.txt";
    char *syst_calls = new char[sys_calls.length()+1];
    strcpy(syst_calls, sys_calls.c_str());
    
    file1.open(syst_calls);
    
    string sys_col = WD_dir+"collaps.txt";
    char *syst_col = new char[sys_col.length()+1];
    strcpy(syst_col, sys_col.c_str());
    
    file5.open(syst_col);
    
    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE"<< endl;
        exit(0);
    }
    
    string input;
    int line1,line2;
    
    for(int i=0;!file1.eof();i++){
        getline(file1,input);
        line1=i;
    }
    
    file1.close();
    file1.clear();
    file1.open(syst_calls);
    
    string *info;
    info= new string[line1];
    int **loc;
    loc=new int*[line1];
    for(int i=0;i!=line1;i++)loc[i]=new int[5];
    
    int **loc_ex;
    loc_ex=new int*[line1];
    for(int i=0;i!=line1;i++)loc_ex[i]=new int[4];
    
    string *orie;
    orie= new string[line1];
    
    string chr;
    int s,e, s_ex1, s_ex2, e_ex1, e_ex2;
    int r,tsd;
    int r_a,tsd_a;
    for(int i=0;i!=line1;i++){
        file1>>info[i];
        file1>>loc[i][0];
        file1>>input;
        file1>>input;
        file1>>loc[i][1];
        file1>>loc_ex[i][0];
        file1>>loc_ex[i][1];
        file1>>loc_ex[i][2];
        file1>>loc_ex[i][3];
        file1>>loc[i][2];
        file1>>loc[i][3];
        file1>>orie[i];
        loc[i][4]=0;
    }
    
    int L=50;
    
    for(int i=0;i!=line1;i++){
        if(loc[i][4]!=1){
            s=loc[i][0];
            e=loc[i][1];
            loc[i][4]=1;
            s_ex1=loc_ex[i][0]-L;
            s_ex2=loc_ex[i][1]+L;
            e_ex1=loc_ex[i][2]-L;
            e_ex2=loc_ex[i][3]+L;
            r=loc[i][2];
            tsd=loc[i][3];
            r_a=loc[i][2];
            tsd_a=loc[i][3];
            int flag=0;
            for(;flag==0;){
                flag=1;
                
                for(int j=0;j!=line1;j++){
                    if(loc[j][4]!=1){
                        if(info[i]==info[j]&&orie[i]==orie[j]){
                            if(!((s>loc[j][1])||(e<loc[j][0]))){
                                if(!((s_ex1>loc_ex[j][1])||(s_ex2<loc_ex[j][0]))&&!((e_ex1>loc_ex[j][3])||(e_ex2<loc_ex[j][2]))){
                                    
                                    if(loc[i][0]>loc[j][0]){
                                        s=loc[j][0];
                                    }
                                    if(loc[i][1]<loc[j][1]){
                                        e=loc[j][1];
                                    }
                                    flag=0;
                                    r_a=loc[j][2]+r_a;
                                    tsd_a=loc[j][3]+tsd_a;
                                    if(r<loc[j][2]) r=loc[j][2];
                                    if(tsd<loc[j][3]) tsd=loc[j][3];
                                    loc[j][4]=1;
                                    
                                    if(loc_ex[i][0]>loc_ex[j][0]){
                                        s_ex1=loc_ex[j][0]-L;
                                    }
                                    if(loc_ex[i][2]>loc_ex[j][2]){
                                        e_ex1=loc_ex[j][2]-L;
                                    }
                                    if(loc_ex[i][1]<loc_ex[j][1]){
                                        s_ex2=loc_ex[j][1]+L;
                                    }
                                    if(loc_ex[i][3]<loc_ex[j][3]){
                                        e_ex2=loc_ex[j][3]+L;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            //file5<<info[i]<<'\t'<<s<<'\t'<<e<<'\t'<<r<<'\t'<<tsd<<'\t'<<r_a<<'\t'<<tsd_a<<endl;
            file5<<info[i]<<'\t'<<s<<'\t'<<e<<'\t'<<s_ex1+L<<'\t'<<s_ex2-L<<'\t'<<e_ex1+L<<'\t'<<e_ex2-L<<'\t'<<r<<'\t'<<tsd<<'\t'<<r_a<<'\t'<<tsd_a<<'\t'<<orie[i]<<endl;
        }
    }
    
}
