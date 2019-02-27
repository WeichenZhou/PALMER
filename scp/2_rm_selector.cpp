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
    
    char *syst_region =new char[sys_region.length()+1];
    strcpy(syst_region, sys_region.c_str());
    
    if(sys_region=="NULL"){
        return 0;
    }
    else {
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
            
            file1>>rm_info[i][0];
            file1>>rm_loc[i][0];
            file1>>rm_loc[i][1];
            file1>>rm_info[i][3];
            rm_info[i][1]="RM";
            rm_info[i][2]="RM";
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
    
}
