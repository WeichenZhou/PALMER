//copyright by ArthurZhou @ UMich&Fudan&HUST
#include "common.hpp"
#include <fcntl.h>

namespace {

int RunBlastn(const vector<string> &args, const string &output_path) {
    pid_t pid = fork();
    if (pid < 0) {
        cerr << "Failed to fork for blastn call" << endl;
        return -1;
    }

    if (pid == 0) {
        setenv("BLAST_USAGE_REPORT", "0", 1);

        int fd = open(output_path.c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0644);
        if (fd < 0) {
            _exit(1);
        }
        dup2(fd, STDOUT_FILENO);
        close(fd);

        vector<char *> c_args;
        c_args.reserve(args.size() + 1);
        for (const auto &arg : args) {
            c_args.push_back(const_cast<char *>(arg.c_str()));
        }
        c_args.push_back(nullptr);

        execvp("blastn", c_args.data());
        _exit(127);
    }

    int status = 0;
    waitpid(pid, &status, 0);
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        cerr << "blastn exited with status " << status << endl;
        return -1;
    }

    return 0;
}

} // namespace

int blastn(string WD_dir, string t, string direc){
    
    string sys_blast = WD_dir+"SEQ.masked";
    ifstream file1;
    file1.open(sys_blast);

    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'SEQ.masked'"<< endl;
        //exit(1);
        return 0;
    }

    if (file1.peek() == ifstream::traits_type::eof()) {
        cout << "Skipping blastn: masked read file '" << sys_blast << "' is empty; no reads available for BLAST alignment." << endl;
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

    vector<string> args = {
        "blastn",
        "-evalue", "0.001",
        "-task", "blastn",
        "-gapopen", "4",
        "-max_target_seqs", "10000000",
        "-query", query,
        "-subject", WD_dir+"SEQ.masked",
        "-outfmt", "6 qacc sacc pident qstart qend sstart send"
    };

    if(t=="ALU"){
        args.emplace_back("-perc_identity");
        args.emplace_back("90");
    }

    string output_file = WD_dir+"blastn_refine.txt";

    if (RunBlastn(args, output_file) != 0) {
        return 0;
    }

    return 0;
}
