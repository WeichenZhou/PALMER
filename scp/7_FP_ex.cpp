//copyright by ArthurZhou @ UMich&Fudan&HUST
#include "common.hpp"
#include <htslib/faidx.h>
#include <array>
#include <filesystem>
#include <utility>
#include <vector>
#include <numeric>
#include <sstream>

namespace {

struct TempFastaFile {
    std::filesystem::path path;
    bool created{false};

    explicit TempFastaFile(const std::string &prefix) {
        path = std::filesystem::temp_directory_path() /
               std::filesystem::path(prefix + "-" + std::filesystem::unique_path("%%%%%%%%").string() + ".fasta");
    }

    bool write(const std::string &header, const std::string &sequence) {
        ofstream out(path);
        if (!out.is_open()) {
            return false;
        }

        out << '>' << header << '\n' << sequence << '\n';
        created = true;
        return static_cast<bool>(out);
    }

    ~TempFastaFile() {
        if (created) {
            std::error_code ec;
            std::filesystem::remove(path, ec);
        }
    }
};

string fetch_region_sequence(faidx_t *fai, const string &region) {
    int seq_len = 0;
    char *seq = fai_fetch(fai, region.c_str(), &seq_len);
    if (!seq) {
        return {};
    }

    string sequence(seq, seq_len);
    free(seq);
    return sequence;
}

int count_blast_hits(const std::string &query_path, const std::string &subject_path, int min_length) {
    const std::string command =
        "BLAST_USAGE_REPORT=0 blastn -evalue 0.05 -task blastn -query " + query_path +
        " -subject " + subject_path + " -dust no -outfmt 7 std";

    FILE *pp = popen(command.c_str(), "r");
    if (!pp) {
        return -1;
    }

    string line;
    int count = 0;
    char buffer[4096];
    while (fgets(buffer, sizeof(buffer), pp)) {
        line.assign(buffer);
        if (line.empty() || line[0] == '#') {
            continue;
        }

        istringstream iss(line);
        string qseqid, sseqid;
        double pident = 0.0;
        int align_len = 0;
        int mismatch = 0;
        int gapopen = 0;
        int qstart = 0;
        int qend = 0;
        int sstart = 0;
        int send = 0;
        double evalue = 0.0;
        double bitscore = 0.0;

        // qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        if (iss >> qseqid >> sseqid >> pident >> align_len >> mismatch >> gapopen >> qstart >> qend >> sstart >> send >> evalue >> bitscore) {
            if (pident >= 80.0 && align_len >= min_length && (send - sstart) > 0) {
                ++count;
            }
        }
    }

    pclose(pp);
    return count;
}

string to_dot_join(const vector<string> &parts) {
    if (parts.empty()) {
        return {};
    }
    string result = parts.front();
    for (size_t idx = 1; idx < parts.size(); ++idx) {
        result.push_back('.');
        result += parts[idx];
    }
    return result;
}

}
}

int fp_ex(string WD_dir, string fasta, string chr, string t, int tsd_index){
    
    //std::ios::sync_with_stdio(false);
    //std::cin.tie(0);
    int BIN_buff=10;
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
        BIN_5=2000;
    }
    else if (t=="HERVK"){
        BIN_3=50;
        BIN_5=50;
        tsd_index=0;
    }
    
    string sys_junc = WD_dir+"read_result_junction.txt";
    ifstream file2(sys_junc);

    if (!file2.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'read_result_junction.txt'"<< endl;
        //exit(1);
        return 0;
    }

    vector<array<string,5>> info;
    vector<array<int,7>> loc;
    vector<array<int,7>> loc_TP;
    info.reserve(1024);
    loc.reserve(1024);
    loc_TP.reserve(1024);

    while (file2)
    {
        array<string,5> info_entry;
        array<int,7> loc_entry{};
        array<int,7> loc_TP_entry{};
        if (!(file2>>info_entry[0]>>loc_entry[0]>>loc_entry[1]>>loc_entry[2]>>loc_entry[3]
              >>loc_entry[4]>>loc_entry[5]>>info_entry[1]>>info_entry[2]>>loc_entry[6]
              >>info_entry[3]>>info_entry[4]
              >>loc_TP_entry[0]>>loc_TP_entry[1]>>loc_TP_entry[2]>>loc_TP_entry[3]
              >>loc_TP_entry[4]>>loc_TP_entry[5]>>loc_TP_entry[6]))
        {
            break;
        }
        info.push_back(move(info_entry));
        loc.push_back(move(loc_entry));
        loc_TP.push_back(move(loc_TP_entry));
    }

    ifstream file1;

    string sys_input = WD_dir+"TSD_blastn_pre.txt";
    file1.open(sys_input);

    ofstream file11;
    string sys_output = WD_dir+"TSD_blastn.txt";
    file11.open(sys_output);

    if (!file1.is_open())
    {
        cout <<"CANNOT OPEN FILE, 'TSD_blastn_pre.txt'"<< endl;
        //exit(1);
        return 0;
    }

    vector<string> info_tsd;
    vector<array<int,7>> loc_tsd;
    vector<array<string,2>> loc_tsd_fp;
    vector<string> kmer_tsd;
    info_tsd.reserve(1024);
    loc_tsd.reserve(1024);
    loc_tsd_fp.reserve(1024);
    kmer_tsd.reserve(1024);

    while (file1)
    {
        string info_entry;
        array<int,7> loc_entry{};
        if (!(file1>>info_entry>>loc_entry[0]>>loc_entry[1]>>loc_entry[2]>>loc_entry[3]
              >>loc_entry[4]>>loc_entry[5]))
        {
            break;
        }
        loc_entry[6] = -1;
        info_tsd.push_back(move(info_entry));
        loc_tsd.push_back(loc_entry);
        loc_tsd_fp.push_back({"",""});
        kmer_tsd.push_back("");
    }

    const int line = static_cast<int>(info.size());
    const int line_tsd = static_cast<int>(info_tsd.size());

    faidx_t *fai = fai_load(fasta.c_str());
    if (!fai) {
        if (fai_build(fasta.c_str()) != 0) {
            cerr << "Failed to build FAI index for " << fasta << endl;
            return 0;
        }
        fai = fai_load(fasta.c_str());
        if (!fai) {
            cerr << "Failed to load FAI index for " << fasta << endl;
            return 0;
        }
    }

    
    //cout<<"ready to process this"<<endl;
//FP_ex module
    if(tsd_index==1){
        for(int i=0;i!=line;++i){
            
            vector<string> seq_parts = {
                info[i][0], to_string(loc[i][0]), to_string(loc[i][1]), to_string(loc[i][2]), to_string(loc[i][3]),
                to_string(loc[i][4]), to_string(loc[i][5]), info[i][1], info[i][2], to_string(loc_TP[i][0]), to_string(loc_TP[i][1]),
                to_string(loc_TP[i][2]), to_string(loc_TP[i][3]), to_string(loc_TP[i][4]), to_string(loc_TP[i][5]), to_string(loc_TP[i][6])
            };

            string seq_index = to_dot_join(seq_parts);
     
            
    //ref pull out for 5mer
            /*
            int start, end;
            start=loc[i][4]-loc[i][6];
            if(start<=0){
                start=1;
            }
            end=loc[i][5]+loc[i][6];
            
            if(loc[i][6]<6000){
                start=loc[i][4]-6000;
                end=loc[i][5]+6000;
            }
            
            string s_start, s_end;
            stringstream ss_s;
            ss_s.clear();
            ss_s<<start;
            s_start=ss_s.str();
            stringstream ss_e;
            ss_e.clear();
            ss_e<<end;
            s_end=ss_e.str();
            */
            //No need ref.fasta for now
            /*
            ofstream file20;
            string ref_file;
            ref_file=WD_dir+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str()+".ref.fasta";
            //cout<<ref_file<<endl;
            char *syst_ref_file = new char[ref_file.length()+1];
            strcpy(syst_ref_file, ref_file.c_str());
            file20.open(syst_ref_file);
            
            string sys_ref;
            sys_ref="samtools faidx "+fasta+" "+chr+":"+s_start+"-"+s_end+" > "+ref_file;
            //cout<<sys_ref<<endl;
            
            char *syst_ref = new char[sys_ref.length()+1];
            strcpy(syst_ref, sys_ref.c_str());
            system(syst_ref);
            //getchar();
            */
//Four -fold or seven -fold
    //ref pull out for fake junc
            int start_junc, end_junc;
            start_junc=loc[i][4]-J_BIN*7;
            if(start_junc<=0){
                start_junc=1;
            }
            end_junc=loc[i][5]+J_BIN*7;
            
            string s_start_junc, s_end_junc;
            stringstream ss_s_junc;
            ss_s_junc.clear();
            ss_s_junc<<start_junc;
            s_start_junc=ss_s_junc.str();
            stringstream ss_e_junc;
            ss_e_junc.clear();
            ss_e_junc<<end_junc;
            s_end_junc=ss_e_junc.str();

            const string region_string = chr + ":" + s_start_junc + "-" + s_end_junc;
            const string ref_sequence = fetch_region_sequence(fai, region_string);
            TempFastaFile ref_junc_file("palmer-ref");
            if (ref_sequence.empty() || !ref_junc_file.write(region_string, ref_sequence)) {
                cerr << "FAILED TO PREPARE FASTA REGION " << region_string << endl;
                continue;
            }

    //FP exclude module
            //3' 26mer identify based on TSD
            //FP construct junction module

            for(int w=0;w!=line_tsd;++w){
                if(info_tsd[w]==seq_index){
                    
                    //string loc_0, loc_1, loc_2, loc_3, loc_4, loc_5;
                    string loc_tsd_0, loc_tsd_1, loc_tsd_2, loc_tsd_3, loc_tsd_4, loc_tsd_5;
                    stringstream ss_tsd_0;
                    ss_tsd_0.clear();
                    ss_tsd_0<<loc_tsd[w][0];
                    loc_tsd_0=ss_tsd_0.str();
                    stringstream ss_tsd_1;
                    ss_tsd_1.clear();
                    ss_tsd_1<<loc_tsd[w][1];
                    loc_tsd_1=ss_tsd_1.str();
                    stringstream ss_tsd_2;
                    ss_tsd_2.clear();
                    ss_tsd_2<<loc_tsd[w][2];
                    loc_tsd_2=ss_tsd_2.str();
                    stringstream ss_tsd_3;
                    ss_tsd_3.clear();
                    ss_tsd_3<<loc_tsd[w][3];
                    loc_tsd_3=ss_tsd_3.str();
                    
                    const string &seq3 = info[i][4];
                    const string &seq5 = info[i][3];
                    
     //FP construct junction module
                    
                    //ofstream file22;
                    //string sys_3_kmer = WD_dir+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str()+"."+loc_tsd_0.c_str()+"."+loc_tsd_1.c_str()+"."+loc_tsd_2.c_str()+"."+loc_tsd_3.c_str()+".fake_junction_1.fasta";
                    //char *syst_3_kmer = new char[sys_3_kmer.length()+1];
                    //strcpy(syst_3_kmer, sys_3_kmer.c_str());
                    //file22.open(syst_3_kmer);
                    //file22<<">"<<seq_index<<"."<<loc_tsd_0.c_str()<<"."<<loc_tsd_1.c_str()<<"."<<loc_tsd_2.c_str()<<"."<<loc_tsd_3.c_str()<<".fake_junction_1.fasta"<<endl;
                    
                    string name_tag_1=seq_index+"."+loc_tsd_0.c_str()+"."+loc_tsd_1.c_str()+"."+loc_tsd_2.c_str()+"."+loc_tsd_3.c_str()+".fake_junction_1.fasta";
                    
                    //ofstream file222;
                    //string sys_3_kmer_2 = WD_dir+loc_0.c_str()+"."+loc_1.c_str()+"."+loc_2.c_str()+"."+loc_3.c_str()+"."+loc_4.c_str()+"."+loc_5.c_str()+"."+info[i][1]+"."+info[i][2]+"."+loc_TP_0.c_str()+"."+loc_TP_1.c_str()+"."+loc_TP_2.c_str()+"."+loc_TP_3.c_str()+"."+loc_TP_4.c_str()+"."+loc_TP_5.c_str()+"."+loc_TP_6.c_str()+"."+loc_tsd_0.c_str()+"."+loc_tsd_1.c_str()+"."+loc_tsd_2.c_str()+"."+loc_tsd_3.c_str()+".fake_junction_2.fasta";
                    //char *syst_3_kmer_2 = new char[sys_3_kmer_2.length()+1];
                    //strcpy(syst_3_kmer_2, sys_3_kmer_2.c_str());
                    //file222.open(syst_3_kmer_2);
                    //file222<<">"<<seq_index<<"."<<loc_tsd_0.c_str()<<"."<<loc_tsd_1.c_str()<<"."<<loc_tsd_2.c_str()<<"."<<loc_tsd_3.c_str()<<".fake_junction_2.fasta"<<endl;
                    
                    string name_tag_2=seq_index+"."+loc_tsd_0.c_str()+"."+loc_tsd_1.c_str()+"."+loc_tsd_2.c_str()+"."+loc_tsd_3.c_str()+".fake_junction_2.fasta";
                    //cout<<name_tag_2<<endl;
                    /*
                    if(info[i][2]=="+"){
                        for(int n=loc_tsd[w][1]-1-J_BIN+1+J_BIN;n!=(loc_tsd[w][1]+(J_BIN*3)/2)&&n<=info[i][3].length();++n){
                            file22<<seq5[n];
                        }
                       */
                    int fix_5;
                    //fix_5=1+loc_tsd[w][4]-BIN_5+J_BIN;
                    //if(BIN_5>=(info[i][3].length()-J_BIN_mer)) fix_5=0;
                    //else if (BIN_5<(info[i][3].length()-J_BIN_mer)) fix_5=info[i][3].length()-J_BIN_mer-BIN_5-1;
                    fix_5=info[i][3].length()-J_BIN_mer-BIN_5-1;
                    
                    int fix_3;
                    //if(BIN_3>=info[i][4].length()) fix_3=0;
                    //else if(BIN_3<info[i][4].length()) fix_3=info[i][4].length()-BIN_3-1;
                    fix_3=info[i][4].length()-BIN_3-1;
                    
                    string fasta5_str="";
                    string fasta3_str="";
                    
                    if(info[i][2]=="+"){
                        stringstream fa_5,fa_3;
                        string fa_s5,fa_s3;
                        fa_3.clear();
                        fa_5.clear();
                        int patch_n=loc_tsd[w][1]-1-(J_BIN*3)/2+fix_5;
                        if(patch_n<0) {
                            //cout<<"patch_n="<<patch_n<<endl;
                            patch_n=0;
                        }
                        for(int n=patch_n;n!=(loc_tsd[w][1]+fix_5)&&n<=info[i][3].length();++n){
                            if(seq5[n]!='A'&&seq5[n]!='T'&&seq5[n]!='G'&&seq5[n]!='C'){
                                fasta5_str=fasta5_str+"N";
                            }
                            else {
                                fa_5.clear();
                                //file22<<seq5[n];
                                fa_5<<seq5[n];
                                fa_5>>fa_s5;
                                fasta5_str=fasta5_str+fa_s5;
                            }
                            
                        }
                        for(int n=loc_tsd[w][3];n!=(loc_tsd[w][3]+(J_BIN*3)/2)&&n<=info[i][4].length();++n){
                            if(seq3[n]!='A'&&seq3[n]!='T'&&seq3[n]!='G'&&seq3[n]!='C'){
                                fasta3_str=fasta3_str+"N";
                            }
                            else {
                                fa_3.clear();
                                //file222<<seq3[n];
                                fa_3<<seq3[n];
                                fa_3>>fa_s3;
                                fasta3_str=fasta3_str+fa_s3;
                            }
                        }
                    }
                    /*
                    else if(info[i][2]=="-"){
                        for(int n=loc_tsd[w][2]-1-J_BIN+J_BIN;n!=(loc_tsd[w][2]-1+(J_BIN*3)/2)&&n<=info[i][4].length();++n){
                            file222<<seq3[n];
                        }
                        */
                    
                    else if(info[i][2]=="-"){
                        stringstream fa_5,fa_3;
                        string fa_s5,fa_s3;
                        fa_3.clear();
                        fa_5.clear();
                        int patch_n_2=loc_tsd[w][2]-1-(J_BIN*3)/2+fix_3-1;
                        if(patch_n_2<0) {
                            //cout<<"patch_n_2="<<patch_n_2<<endl;
                            patch_n_2=0;
                        }
                        for(int n=patch_n_2;n!=(loc_tsd[w][2]+fix_3-1)&&n<=info[i][4].length();++n){
                            if(seq3[n]!='A'&&seq3[n]!='T'&&seq3[n]!='G'&&seq3[n]!='C'){
                                fasta3_str=fasta3_str+"N";
                            }
                            else {
                                fa_3.clear();
                                //file222<<seq3[n];
                                fa_3<<seq3[n];
                                fa_3>>fa_s3;
                                fasta3_str=fasta3_str+fa_s3;
                            }
                        }
                        for(int n=loc_tsd[w][0]-1+J_BIN_mer;n!=(loc_tsd[w][0]+(J_BIN*3)/2-1+J_BIN_mer)&&n<=info[i][3].length();++n){
                            if(seq5[n]!='A'&&seq5[n]!='T'&&seq5[n]!='G'&&seq5[n]!='C'){
                                fasta5_str=fasta5_str+"N";
                            }
                            else {
                                fa_5.clear();
                                //file22<<seq5[n];
                                fa_5<<seq5[n];
                                fa_5>>fa_s5;
                                fasta5_str=fasta5_str+fa_s5;
                            }
                        }
                    }
                    if(fasta3_str==""||fasta3_str==" "){
                        fasta3_str=fasta3_str+"N";
                    }
                    if(fasta5_str==""||fasta5_str==" "){
                        fasta5_str=fasta5_str+"N";
                    }
                    //file22<<endl;
                    //file222<<endl;
                    int fasta5_str_len=int(fasta5_str.length()*60/75);
                    int fasta3_str_len=int(fasta3_str.length()*60/75);
                    string fasta5_str_length=to_string(fasta5_str_len);
                    string fasta3_str_length=to_string(fasta3_str_len);
                    
                    TempFastaFile query5_file("palmer-query5");
                    if (!query5_file.write(seq_index, fasta5_str)) {
                        cerr << "FAILED TO WRITE QUERY FASTA FOR " << seq_index << endl;
                        continue;
                    }

                    const int blast5_hits = count_blast_hits(query5_file.path.string(), ref_junc_file.path.string(), fasta5_str_len);
                    if (blast5_hits >= 0) {
                        loc_tsd_fp[w][0] = to_string(blast5_hits);
                    }

                    TempFastaFile query3_file("palmer-query3");
                    if (!query3_file.write(seq_index, fasta3_str)) {
                        cerr << "FAILED TO WRITE QUERY FASTA FOR " << seq_index << endl;
                        continue;
                    }

                    const int blast3_hits = count_blast_hits(query3_file.path.string(), ref_junc_file.path.string(), fasta3_str_len);
                    if (blast3_hits >= 0) {
                        loc_tsd_fp[w][1] = to_string(blast3_hits);
                    }

                    if(loc_tsd_fp[w][0]==""){
                        loc_tsd_fp[w][0]="-1";
                    }
                    if(loc_tsd_fp[w][1]==""){
                        loc_tsd_fp[w][1]="-1";
                    }

                    
                    if(info[i][2]=="+"){
                        stringstream kmer_ss;
                        string kmer_s;
                        kmer_ss.clear();
                        for(int n=loc_tsd[w][1]-1-J_BIN_mer+1+J_BIN;n!=(loc_tsd[w][1]+J_BIN+J_BIN_mer)&&n<=info[i][3].length();++n){
                            //file42<<seq5[n];
                            if(seq5[n]!='A'&&seq5[n]!='T'&&seq5[n]!='G'&&seq5[n]!='C'){
                                kmer_tsd[w]=kmer_tsd[w]+"N";
                            }
                            else {
                                kmer_ss.clear();
                                kmer_ss<<seq5[n];
                                kmer_ss>>kmer_s;
                                kmer_tsd[w]=kmer_tsd[w]+kmer_s;
                            }
                        }
                        if(kmer_tsd[w]==""){
                            kmer_tsd[w]="N";
                        }
                        
                    }
                    else if(info[i][2]=="-"){
                        
                        stringstream kmer_ss;
                        string kmer_s;
                        kmer_ss.clear();
                        for(int n=loc_tsd[w][0]+(BIN_5+1-loc_tsd[w][4])-1+J_BIN_mer-J_BIN_mer;n!=(loc_tsd[w][0]+(BIN_5+1-loc_tsd[w][4])+J_BIN_mer+J_BIN_mer-1)&&n<=info[i][3].length();++n){
                            //file42<<seq5[n];
                            if(seq5[n]!='A'&&seq5[n]!='T'&&seq5[n]!='G'&&seq5[n]!='C'){
                                kmer_tsd[w]=kmer_tsd[w]+"N";
                            }
                            else {
                                kmer_ss.clear();
                                kmer_ss<<seq5[n];
                                kmer_ss>>kmer_s;
                                kmer_tsd[w]=kmer_tsd[w]+kmer_s;
                            }
                        }
                        if(kmer_tsd[w]==""){
                            kmer_tsd[w]="N";
                        }
                    }

                    loc_tsd[w][6]=0;
                    
                }
            }
            //file20.clear();
            //file20.close();
        }
        
        for(int w=0;w!=line_tsd;++w){
            file11<<info_tsd[w]<<'\t'<<loc_tsd[w][0]<<'\t'<<loc_tsd[w][1]<<'\t'<<loc_tsd[w][2]<<'\t'<<loc_tsd[w][3]<<'\t'<<loc_tsd[w][4]<<'\t'<<loc_tsd[w][5]<<'\t'<<loc_tsd[w][6]<<'\t'<<loc_tsd_fp[w][0]<<'\t'<<loc_tsd_fp[w][1]<<'\t'<<kmer_tsd[w]<<endl;
        }
    }
    
    else if(tsd_index==0){
        for(int w=0;w!=line_tsd;++w){
            file11<<info_tsd[w]<<'\t'<<loc_tsd[w][0]<<'\t'<<loc_tsd[w][1]<<'\t'<<loc_tsd[w][2]<<'\t'<<loc_tsd[w][3]<<'\t'<<loc_tsd[w][4]<<'\t'<<loc_tsd[w][5]<<'\t'<<"0"<<'\t'<<"1"<<'\t'<<"1"<<'\t'<<"N"<<endl;
        }
    }
    
    // containers clean up automatically
    fai_destroy(fai);
    return 0;

}
