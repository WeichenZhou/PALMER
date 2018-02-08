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

int blastn(string WD, string WD_dir){
    
    //cout<<"Blastn Step is now running."<<endl;
    
    string sys_blast = WD_dir+"SEQ.masked";
    char *syst_blast =  new char[sys_blast.length()+1];
    strcpy(syst_blast, sys_blast.c_str());
    
    ifstream file1;
    file1.open(syst_blast);
    
    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE"<< endl;
        exit(1);
    }
    
    string sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+WD+"buildup/L1.3.fasta -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" >  "+WD_dir+"blastn.txt";
    char *syst_blastn = new char[sys_blastn.length()+1];
    strcpy(syst_blastn, sys_blastn.c_str());
    
    system(syst_blastn);
}
