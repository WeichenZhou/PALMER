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

int blastn(string WD_dir, string t, string direc){
    
    //std::ios::sync_with_stdio(false);
    //std::cin.tie(0);
    
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
    
    /*
    string sys_blastnrefine = WD_dir+"blastn_refine.txt";
    char *syst_blastnrefine =  new char[sys_blastnrefine.length()+1];
    strcpy(syst_blastnrefine, sys_blastnrefine.c_str());
    
    ofstream file2;
    file2.open(syst_blastnrefine);
    
    if(t=="LINE"){
         sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+direc+"/lib/L1.3.fasta -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" |grep -v \"#\" >  "+WD_dir+"blastn_refine.txt";
    }
    else if(t=="ALU"){
        sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+direc+"/lib/AluY.fasta -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" |grep -v \"#\" >  "+WD_dir+"blastn_refine.txt";
    }
    else if(t=="SVA"){
        sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+direc+"/lib/SVA_F.fasta -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" |grep -v \"#\" >  "+WD_dir+"blastn_refine.txt";
    }
    else {
        sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+t+" -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" |grep -v \"#\" >  "+WD_dir+"blastn_refine.txt";
    }
     */
    
    if(t=="LINE"){
        sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+direc+"/lib/L1.3.fasta -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" >  "+WD_dir+"blastn.txt";
    }
    else if(t=="ALU"){
        sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+direc+"/lib/AluY.fasta -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" >  "+WD_dir+"blastn.txt";
    }
    else if(t=="SVA"){
        sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+direc+"/lib/SVA_F.fasta -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" >  "+WD_dir+"blastn.txt";
    }
    else if(t=="HERVK"){
        sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+direc+"/lib/HERVK.fasta -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" >  "+WD_dir+"blastn.txt";
    }
    else {
        sys_blastn = "blastn -evalue 0.001 -task blastn -gapopen 4 -query "+t+" -subject "+WD_dir+"SEQ.masked -outfmt \"7 qacc sacc evalue qstart qend sstart send\" >  "+WD_dir+"blastn.txt";
    }
    
    char *syst_blastn = new char[sys_blastn.length()+1];
    strcpy(syst_blastn, sys_blastn.c_str());
    //cout<<sys_blastn<<endl;
    system(syst_blastn);
    return 0;
}
