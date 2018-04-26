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

int RMSelector(string WD, string WD_dir, string sys_region){
    
    //cout<<"RMSelector Step is now running."<<endl;
    
    //string sys_region=WD+"buildup/LINEs.regions";
    char *syst_region =new char[sys_region.length()+1];
    strcpy(syst_region, sys_region.c_str());
    
    ifstream file1;
    file1.open(syst_region);
    
    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'sys_region'."<< endl;
        //exit(1);
        return 0;
    }
    
    string input;
    int line;
    for(int i=0;!file1.eof();i++){
        getline(file1,input);
        line=i;
    }
        
    string sys_RMloc = WD_dir+"RM.loc";
    char *syst_RMloc = new char[sys_RMloc.length()+1];
    strcpy(syst_RMloc, sys_RMloc.c_str());
    
    ifstream file3;
    file3.open(syst_RMloc);
    
    string sys_RM_selescted = WD_dir+"RM.selected";
    char *syst_RM_selescted = new char[sys_RM_selescted.length()+1];
    strcpy(syst_RM_selescted, sys_RM_selescted.c_str());
    
    ofstream file2;
    file2.open(syst_RM_selescted);
    
    string chr;
    int start, end;
    file3>>chr;
    file3>>start;
    file3>>end;
    
    string **rm_info;
    rm_info= new string*[line];
    for(int i=0;i!=line;i++) rm_info[i]=new string[4];
    
    int **rm_loc;
    rm_loc=new int*[line];
    for(int i=0;i!=line;i++) rm_loc[i]=new int[2];
    file1.close();
    file1.clear();
    file1.open(syst_region);
    for(int i=0;i!=line;i++){
        //file1>>input;
        //file1>>input;
        //file1>>input;
        //file1>>input;
        //file1>>input;
        file1>>rm_info[i][0];
        if(rm_info[i][0]=="1"){
            rm_info[i][0]="chr1";
        }
        else if(rm_info[i][0]=="2"){
            rm_info[i][0]="chr2";
        }
        else if(rm_info[i][0]=="3"){
            rm_info[i][0]="chr3";
        }
        else if(rm_info[i][0]=="4"){
            rm_info[i][0]="chr4";
        }
        else if(rm_info[i][0]=="5"){
            rm_info[i][0]="chr5";
        }
        else if(rm_info[i][0]=="6"){
            rm_info[i][0]="chr6";
        }
        else if(rm_info[i][0]=="7"){
            rm_info[i][0]="chr7";
        }
        else if(rm_info[i][0]=="8"){
            rm_info[i][0]="chr8";
        }
        else if(rm_info[i][0]=="9"){
            rm_info[i][0]="chr9";
        }
        else if(rm_info[i][0]=="10"){
            rm_info[i][0]="chr10";
        }
        else if(rm_info[i][0]=="11"){
            rm_info[i][0]="chr11";
        }
        else if(rm_info[i][0]=="12"){
            rm_info[i][0]="chr12";
        }
        else if(rm_info[i][0]=="13"){
            rm_info[i][0]="chr13";
        }
        else if(rm_info[i][0]=="14"){
            rm_info[i][0]="chr14";
        }
        else if(rm_info[i][0]=="15"){
            rm_info[i][0]="chr15";
        }
        else if(rm_info[i][0]=="16"){
            rm_info[i][0]="chr16";
        }
        else if(rm_info[i][0]=="17"){
            rm_info[i][0]="chr17";
        }
        else if(rm_info[i][0]=="18"){
            rm_info[i][0]="chr18";
        }
        else if(rm_info[i][0]=="19"){
            rm_info[i][0]="chr19";
        }
        else if(rm_info[i][0]=="20"){
            rm_info[i][0]="chr20";
        }
        else if(rm_info[i][0]=="21"){
            rm_info[i][0]="chr21";
        }
        else if(rm_info[i][0]=="22"){
            rm_info[i][0]="chr22";
        }
        else if(rm_info[i][0]=="Y"){
            rm_info[i][0]="chrY";
        }
        else if(rm_info[i][0]=="X"){
            rm_info[i][0]="chrX";
        }
        file1>>rm_loc[i][0];
        file1>>rm_loc[i][1];
        //file1>>input;
        //file1>>input;
        file1>>rm_info[i][3];
        rm_info[i][1]="RM";
        rm_info[i][2]="RM";
        //getline(file1,input);
    }
    
    for(int i=0;i!=line;i++){
        if(rm_info[i][0]==chr){
            if(!((rm_loc[i][0]>end)||(rm_loc[i][1]<start))){
                file2<<rm_info[i][0]<<'\t'<<rm_loc[i][0]<<'\t'<<rm_loc[i][1]<<'\t'<<rm_info[i][1]<<'\t'<<rm_info[i][2]<<'\t'<<rm_info[i][3]<<endl;
            }
        }
    }
    
    for(int i=0;i!=line;i++){
        delete [] rm_info[i];
        delete [] rm_loc[i];
    }
    delete [] rm_info;
    delete [] rm_loc;
    
    return 0;
}
