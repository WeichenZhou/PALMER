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

int blastn(string WD, string WD_dir, string t, string direc){
    
    //cout<<"Blastn Step is now running."<<endl;
    
    string sys_blast = WD_dir+"SEQ.masked";
    char *syst_blast =  new char[sys_blast.length()+1];
    strcpy(syst_blast, sys_blast.c_str());
    
    ifstream file1;
    file1.open(syst_blast);
    
    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'SEQ.masked'"<< endl;
        //exit(1);
        return 0;
    }
    string sys_blastn;
    if(t=="LINE"){
         sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+direc+"/lib/L1.3.fasta -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" >  "+WD_dir+"blastn.txt";
    }
    else if(t=="ALU"){
        sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+direc+"/lib/AluY.fasta -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" >  "+WD_dir+"blastn.txt";
    }
    else if(t=="SVA"){
        sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+direc+"/lib/SVA_F.fasta -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" >  "+WD_dir+"blastn.txt";
    }
    
    char *syst_blastn = new char[sys_blastn.length()+1];
    strcpy(syst_blastn, sys_blastn.c_str());
    
    system(syst_blastn);
    return 0;
}
