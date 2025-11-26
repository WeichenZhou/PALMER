//copyright by ArthurZhou @ UMich&Fudan&HUST
#include "common.hpp"

int blastn(string WD_dir, string t, string direc){
    
    string subject_path = WD_dir + "SEQ.masked";
    ifstream file1(subject_path);
    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'SEQ.masked'"<< endl;
        return 0;
    }

    string query;
    if(t=="LINE"){
        query = direc+"/lib/L1.3.fasta";
    }
    else if(t=="ALU"){
        query = direc+"/lib/AluY.fasta";
    }
    else if(t=="SVA"){
        query = direc+"/lib/SVA_F.fasta";
    }
    else if(t=="HERVK"){
        query = direc+"/lib/HERVK.fasta";
    }
    else {
        query = t;
    }

    vector<string> args = {"blastn", "-evalue", "0.001", "-task", "blastn", "-gapopen", "4",
                            "-max_target_seqs", "10000000", "-query", query, "-subject", subject_path,
                            "-outfmt", "7 qacc sacc evalue qstart qend sstart send"};

    vector<pair<string, string>> env{{"BLAST_USAGE_REPORT", "0"}};
    string output_path = WD_dir + "blastn.txt";
    run_process(args, output_path, env);
    return 0;
}
