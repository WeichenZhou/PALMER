////copyright by ArthurZhou @ UMich&Fudan&HUST
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
#include <numeric>

#include <unordered_map>
#include <utility>
#include <ctime>
#include <cctype>
#include <climits>

#include <thread>
#include <atomic>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "scp/tube.cpp"

using namespace std;

//VCF
struct TsdInfo {
    string read_name;
    string tsd5;
    string tsd3;
    string predicted_transd;
    string junction_26mer;
    string insertion_seq;
    int n_count = 0;
};

struct CollapsedRecord {
    string chrom;
    long long pos = 0;
    long long end = 0;
    long long start_invar = 0;
    long long end_invar = 0;
    string ref_base;
    string alt_symbol;
    string info_prefix;
    vector<string> cluster_ids;
    vector<TsdInfo> tsd_candidates;
};

long long safe_average(const string &a, const string &b) {
    bool a_empty = a.empty();
    bool b_empty = b.empty();
    if (a_empty && b_empty) {
        return 0;
    }

    long long a_val = a_empty ? 0 : atoll(a.c_str());
    long long b_val = b_empty ? 0 : atoll(b.c_str());

    if (a_empty) return b_val;
    if (b_empty) return a_val;
    return (a_val + b_val) / 2;
}

string normalize_field(const string &value) {
    return value == "N" ? string("NA") : value;
}

string current_date() {
    time_t now = time(nullptr);
    tm *local_tm = localtime(&now);
    char buffer[16];
    strftime(buffer, sizeof(buffer), "%Y%m%d", local_tm);
    return string(buffer);
}

string reference_basename(const string &reference_path) {
    if (reference_path.empty()) {
        return string("NA");
    }

    string normalized = reference_path;
    const string file_prefix = "file://";
    if (normalized.rfind(file_prefix, 0) == 0) {
        normalized = normalized.substr(file_prefix.size());
    }

    size_t slash_pos = normalized.find_last_of("/\\");
    if (slash_pos == string::npos) {
        return normalized.empty() ? string("NA") : normalized;
    }

    string base = normalized.substr(slash_pos + 1);
    return base.empty() ? string("NA") : base;
}

vector<pair<string, long long>> load_contigs_from_fai(const string &reference_path) {
    vector<pair<string, long long>> contigs;

    if (reference_path.empty()) {
        return contigs;
    }

    string fai_path = reference_path + ".fai";
    ifstream fai_stream(fai_path);
    if (!fai_stream.is_open()) {
        return contigs;
    }

    string fai_line;
    while (getline(fai_stream, fai_line)) {
        if (fai_line.empty()) continue;

        stringstream fai_ss(fai_line);
        string contig_name;
        string contig_length;
        getline(fai_ss, contig_name, '\t');
        getline(fai_ss, contig_length, '\t');
        if (!contig_name.empty() && !contig_length.empty()) {
            contigs.push_back({contig_name, atoll(contig_length.c_str())});
        }
    }

    return contigs;
}

vector<pair<string, long long>> load_contigs_from_fasta(const string &reference_path) {
    vector<pair<string, long long>> contigs;

    if (reference_path.empty()) {
        return contigs;
    }

    ifstream fasta_stream(reference_path);
    if (!fasta_stream.is_open()) {
        return contigs;
    }

    string line;
    string current_contig;
    long long current_length = 0;

    auto flush_contig = [&](const string &name, long long length) {
        if (!name.empty()) {
            contigs.push_back({name, length});
        }
    };

    while (getline(fasta_stream, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            flush_contig(current_contig, current_length);
            current_contig.clear();
            current_length = 0;

            size_t space_pos = line.find_first_of(" \t");
            if (space_pos == string::npos) {
                current_contig = line.substr(1);
            } else {
                current_contig = line.substr(1, space_pos - 1);
            }
        } else {
            current_length += static_cast<long long>(line.size());
        }
    }

    flush_contig(current_contig, current_length);
    return contigs;
}

unordered_map<string, string> load_reference_sequences(const string &reference_path) {
    unordered_map<string, string> sequences;

    if (reference_path.empty()) {
        return sequences;
    }

    ifstream fasta_stream(reference_path);
    if (!fasta_stream.is_open()) {
        return sequences;
    }

    string line;
    string current_contig;
    string current_sequence;

    auto flush_sequence = [&](const string &name, const string &sequence) {
        if (!name.empty()) {
            sequences[name] = sequence;
        }
    };

    while (getline(fasta_stream, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            flush_sequence(current_contig, current_sequence);
            current_contig.clear();
            current_sequence.clear();

            size_t space_pos = line.find_first_of(" \t");
            if (space_pos == string::npos) {
                current_contig = line.substr(1);
            } else {
                current_contig = line.substr(1, space_pos - 1);
            }
        } else {
            current_sequence += line;
        }
    }

    flush_sequence(current_contig, current_sequence);
    return sequences;
}

char reference_base(const unordered_map<string, string> &reference_sequences, const string &chrom, long long position) {
    if (position <= 0) {
        return 'N';
    }

    auto it = reference_sequences.find(chrom);
    if (it == reference_sequences.end()) {
        return 'N';
    }

    const string &sequence = it->second;
    if (static_cast<size_t>(position) > sequence.size()) {
        return 'N';
    }

    char base = sequence[static_cast<size_t>(position - 1)];
    return static_cast<char>(toupper(base));
}

int missing_fields(const TsdInfo &info) {
    int missing = 0;
    missing += (info.read_name == "NA");
    missing += (info.tsd5 == "NA");
    missing += (info.tsd3 == "NA");
    missing += (info.predicted_transd == "NA");
    missing += (info.junction_26mer == "NA");
    missing += (info.insertion_seq == "NA");
    return missing;
}

void write_vcf_output(const string &calls_path, const string &tsd_path, const string &vcf_path, const string &ins_type, const string &reference_path, const string &command_line) {

    unordered_map<string, TsdInfo> tsd_records;

    static bool seeded_rng = false;
    if (!seeded_rng) {
        srand(static_cast<unsigned>(time(nullptr)));
        seeded_rng = true;
    }

    ifstream tsd_stream(tsd_path);
    if (tsd_stream.is_open()) {
        string tsd_line;
        getline(tsd_stream, tsd_line); // header
        while (getline(tsd_stream, tsd_line)) {
            if (tsd_line.empty()) continue;

            stringstream tsd_ss(tsd_line);
            vector<string> fields;
            string field;
            while (getline(tsd_ss, field, '\t')) {
                fields.push_back(field);
            }

            if (fields.size() < 7) continue;

            int n_count = 0;
            for (size_t idx = 2; idx <= 6; ++idx) {
                if (fields[idx] == "N") {
                    ++n_count;
                }
            }

            auto existing = tsd_records.find(fields[0]);
            bool should_replace = false;
            if (existing == tsd_records.end()) {
                should_replace = true;
            } else if (n_count < existing->second.n_count) {
                should_replace = true;
            } else if (n_count == existing->second.n_count) {
                should_replace = (rand() % 2) == 0;
            }

            if (should_replace) {
                TsdInfo info;
                info.read_name = normalize_field(fields[1]);
                info.tsd5 = normalize_field(fields[2]);
                info.tsd3 = normalize_field(fields[3]);
                info.predicted_transd = normalize_field(fields[4]);
                info.junction_26mer = normalize_field(fields[5]);
                info.insertion_seq = normalize_field(fields[6]);
                info.n_count = n_count;
                tsd_records[fields[0]] = info;
            }
        }
    }

    ifstream calls_stream(calls_path);
    if (!calls_stream.is_open()) {
        cerr << "CANNOT OPEN MERGED CALLS FILE: " << calls_path << endl;
        return;
    }

    ofstream vcf_stream(vcf_path);
    if (!vcf_stream.is_open()) {
        cerr << "CANNOT OPEN VCF OUTPUT FILE: " << vcf_path << endl;
        return;
    }

    string file_date = current_date();
    vector<pair<string, long long>> contigs = load_contigs_from_fai(reference_path);
    if (contigs.empty()) {
        contigs = load_contigs_from_fasta(reference_path);
    }
    unordered_map<string, string> reference_sequences = load_reference_sequences(reference_path);
    string reference_header = reference_basename(reference_path);
    string assembly_path = reference_header == "NA" ? string("NA") : string("file:") + reference_header;

    vcf_stream << "##fileformat=VCFv4.2" << endl;
    vcf_stream << "##fileDate=" << file_date << endl;
    vcf_stream << "##source=PALMER" << endl;
    vcf_stream << "##reference=" << reference_header << endl;
    for (const auto &contig : contigs) {
        vcf_stream << "##contig=<ID=" << contig.first
                   << ",assembly=" << assembly_path
                   << ",length=" << contig.second
                   << ",species=>" << endl;
    }
    vcf_stream << "##cmdline=" << (command_line.empty() ? string("NA") : command_line) << endl;
    vcf_stream << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of the insertion interval\">" << endl;
    vcf_stream << "##INFO=<ID=SUBTYPE,Number=1,Type=String,Description=\"Insertion subtype (matches --type)\">" << endl;
    vcf_stream << "##INFO=<ID=CS,Number=1,Type=Integer,Description=\"Confident supporting reads\">" << endl;
    vcf_stream << "##INFO=<ID=PS,Number=1,Type=Integer,Description=\"Potential supporting reads\">" << endl;
    vcf_stream << "##INFO=<ID=SEG,Number=1,Type=Integer,Description=\"Potential segmental supporting reads\">" << endl;
    vcf_stream << "##INFO=<ID=ORIENTATION,Number=1,Type=String,Description=\"Insertion orientation\">" << endl;
    vcf_stream << "##INFO=<ID=POLYA,Number=1,Type=Integer,Description=\"polyA tail size\">" << endl;
    vcf_stream << "##INFO=<ID=TSD5,Number=1,Type=Integer,Description=\"5' TSD size\">" << endl;
    vcf_stream << "##INFO=<ID=TSD3,Number=1,Type=Integer,Description=\"3' TSD size\">" << endl;
    vcf_stream << "##INFO=<ID=TRANSD,Number=1,Type=Integer,Description=\"Predicted 3' transduction size\">" << endl;
    vcf_stream << "##INFO=<ID=INV5,Number=1,Type=String,Description=\"Has 5' inverted sequence\">" << endl;
    vcf_stream << "##INFO=<ID=INV5_END,Number=1,Type=Integer,Description=\"End position of 5' inverted sequence\">" << endl;
    vcf_stream << "##INFO=<ID=INV5_START,Number=1,Type=Integer,Description=\"Start position of 5' inverted sequence\">" << endl;
    vcf_stream << "##INFO=<ID=START_INVAR,Number=1,Type=Integer,Description=\"Average start position within insertion sequence\">" << endl;
    vcf_stream << "##INFO=<ID=END_INVAR,Number=1,Type=Integer,Description=\"Average end position within insertion sequence\">" << endl;
    vcf_stream << "##INFO=<ID=TSD_READ,Number=1,Type=String,Description=\"Representative read supporting TSD\">" << endl;
    vcf_stream << "##INFO=<ID=TSD5_SEQ,Number=1,Type=String,Description=\"5' TSD sequence (NA if unavailable)\">" << endl;
    vcf_stream << "##INFO=<ID=TSD3_SEQ,Number=1,Type=String,Description=\"3' TSD sequence (NA if unavailable)\">" << endl;
    vcf_stream << "##INFO=<ID=TRANSD_READ,Number=1,Type=String,Description=\"Predicted transduction sequence (NA if unavailable)\">" << endl;
    vcf_stream << "##INFO=<ID=JUNC_26MER,Number=1,Type=String,Description=\"Unique 26-mer at 5' junction (NA if unavailable)\">" << endl;
    vcf_stream << "##INFO=<ID=INS_SEQ,Number=1,Type=String,Description=\"Representative insertion sequence (NA if unavailable)\">" << endl;
    vcf_stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << endl;

    string line;
    getline(calls_stream, line); // Skip header

    vector<CollapsedRecord> collapsed_records;
    unordered_map<string, size_t> key_to_index;

    while (getline(calls_stream, line)) {
        if (line.empty()) {
            continue;
        }

        stringstream ss(line);
        vector<string> fields;
        string field;
        while (getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 21) {
            continue;
        }

        string chrom = fields[1];
        long long pos_val = safe_average(fields[2], fields[3]);
        long long end_val = safe_average(fields[4], fields[5]);
        long long start_invariant = safe_average(fields[6], fields[7]);
        long long end_invariant = safe_average(fields[8], fields[9]);
        string end = to_string(end_val);
        string ref_base(1, reference_base(reference_sequences, chrom, pos_val));
        string alt_symbol = "<INS>";
        if (!ins_type.empty()) {
            alt_symbol = "<INS:" + ins_type + ">";
        }

        string info_prefix = "SVTYPE=INS";
        info_prefix += ";SUBTYPE=" + (ins_type.empty() ? string("NA") : ins_type);
        info_prefix += ";END=" + end;
        info_prefix += ";CS=" + fields[10];
        info_prefix += ";PS=" + fields[11];
        info_prefix += ";SEG=" + fields[12];
        info_prefix += ";ORIENTATION=" + fields[13];
        info_prefix += ";POLYA=" + fields[14];
        info_prefix += ";TSD5=" + fields[15];
        info_prefix += ";TSD3=" + fields[16];
        info_prefix += ";TRANSD=" + fields[17];
        info_prefix += ";INV5=" + fields[18];
        info_prefix += ";INV5_END=" + fields[19];
        info_prefix += ";INV5_START=" + fields[20];
        info_prefix += ";START_INVAR=" + to_string(start_invariant);
        info_prefix += ";END_INVAR=" + to_string(end_invariant);

        auto tsd_it = tsd_records.find(fields[0]);
        TsdInfo tsd_info;
        if (tsd_it != tsd_records.end()) {
            tsd_info = tsd_it->second;
        } else {
            tsd_info.read_name = "NA";
            tsd_info.tsd5 = "NA";
            tsd_info.tsd3 = "NA";
            tsd_info.predicted_transd = "NA";
            tsd_info.junction_26mer = "NA";
            tsd_info.insertion_seq = "NA";
        }

        string key;
        for (size_t i = 1; i < fields.size(); ++i) {
            if (i > 1) key += "|";
            key += fields[i];
        }

        auto idx_it = key_to_index.find(key);
        if (idx_it == key_to_index.end()) {
            CollapsedRecord record;
            record.chrom = chrom;
            record.pos = pos_val;
            record.end = end_val;
            record.start_invar = start_invariant;
            record.end_invar = end_invariant;
            record.ref_base = ref_base;
            record.alt_symbol = alt_symbol;
            record.info_prefix = info_prefix;
            record.cluster_ids.push_back(fields[0]);
            record.tsd_candidates.push_back(tsd_info);
            collapsed_records.push_back(record);
            key_to_index[key] = collapsed_records.size() - 1;
        } else {
            CollapsedRecord &record = collapsed_records[idx_it->second];
            record.cluster_ids.push_back(fields[0]);
            record.tsd_candidates.push_back(tsd_info);
        }
    }

    for (const auto &record : collapsed_records) {
        string id_field;
        for (size_t i = 0; i < record.cluster_ids.size(); ++i) {
            if (i > 0) id_field += ";";
            id_field += record.cluster_ids[i];
        }

        TsdInfo tsd_info;
        tsd_info.read_name = "NA";
        tsd_info.tsd5 = "NA";
        tsd_info.tsd3 = "NA";
        tsd_info.predicted_transd = "NA";
        tsd_info.junction_26mer = "NA";
        tsd_info.insertion_seq = "NA";

        int best_missing = INT_MAX;
        for (const auto &candidate : record.tsd_candidates) {
            int missing = missing_fields(candidate);
            if (missing < best_missing) {
                best_missing = missing;
                tsd_info = candidate;
            }
        }

        string info = record.info_prefix;
        info += ";TSD_READ=" + tsd_info.read_name;
        info += ";TSD5_SEQ=" + tsd_info.tsd5;
        info += ";TSD3_SEQ=" + tsd_info.tsd3;
        info += ";TRANSD_READ=" + tsd_info.predicted_transd;
        info += ";JUNC_26MER=" + tsd_info.junction_26mer;
        info += ";INS_SEQ=" + tsd_info.insertion_seq;

        vcf_stream << record.chrom << '\t'
                   << record.pos << '\t'
                   << id_field << '\t'
                   << record.ref_base << '\t'
                   << record.alt_symbol << '\t'
                   << '.' << '\t'
                   << "PASS" << '\t'
                   << info << endl;
    }

    cout << "VCF output written to " << vcf_path << endl;
}



//main
int main(int argc, char *argv[]){

//parameters_start
    //std::ios::sync_with_stdio(false);
    //std::cin.tie(0);
    
    string T, WD, inputF, output, SP, ref, CHR, region_START, region_END, ref_fa, cus, tsd_pass, LL_len, s_cus_seq_len, mode, mapq;
    
    //int NUM_threads=30;
    //int NUM_threads=10;
    unsigned int hardware_threads = std::thread::hardware_concurrency();
    int NUM_threads = (hardware_threads == 0) ? 10 : static_cast<int>(hardware_threads);
    
    ifstream file1;
    ifstream file11;
    //ifstream file12;
    ifstream file13;
    //ofstream file2;
    int flag_wd=0;
    int flag_inputf=0;
    int flag_T=0;
    int flag_reffa=0;
    int flag_cus=0;
    //int flag_cusin=0;
    int flag_tsd=1;
    int help=0;
    int L_len=25;
    int seq_len=-1;
    int flag_cus_seq_len=0;
    
    int intermediate=0;
    int flag_intermediate=0;
    int invalid_intermediate=0;
    
    int flag_s=0;
    int flag_e=0;
    int flag_m=0;
    output="output.txt";
    SP="Human";
    //T="LINE";
    CHR="ALL";
    mode="NA";
    mapq="10";

    int ref_n=0;
    //int flag_mapq=0;
    int mapq_int=10;

    string dir;
    dir=argv[0];
    
    for(int i=1;i!=argc;++i){
        if(strncmp(argv[i],"--workdir",9)==0){
            WD=argv[i+1];
            flag_wd=1;
        }
        if(strncmp(argv[i],"--chr",5)==0){
            CHR=argv[i+1];
        }
        
        if(strncmp(argv[i],"--start",7)==0){
            region_START=argv[i+1];
            flag_s=1;
        }
        
        if(strncmp(argv[i],"--end",5)==0){
            region_END=argv[i+1];
            flag_e=1;
        }
        
        if(strncmp(argv[i],"--L_len",7)==0){
            LL_len=argv[i+1];
            int L_len = 6023- std::stoi(LL_len);
            //L_len=L_len;
        }
    
    if(strncmp(argv[i],"--mapq",6)==0){
            mapq=argv[i+1];
            mapq_int = std::stoi(mapq);
        if (mapq_int<0 || mapq_int>100)
        {
                cout<<"PLEASE INPUT A CORRECT RANGE OF MAPQ :("<<endl;
            }
        }

        if(strncmp(argv[i],"--ref_ver",9)==0){
            ref=argv[i+1];
            if(ref=="GRCh37") ref_n=37;
            else if(ref=="GRCh38") ref_n=38;
            else if(ref=="hg19") ref_n=19;
            else if(ref=="other") ref_n=-1;
            else {
                cout<<"PLEASE INPUT A CORRECT CATEGORY OF REFERENCE :("<<endl;
            }
        }
        if(strncmp(argv[i],"--type",6)==0){
            T=argv[i+1];
            flag_T=1;
            if(T=="CUSTOMIZED"){
                flag_T=2;
            }
        }
        
        if(strncmp(argv[i],"--intermediate",14)==0){
            string intermediate_value = "0";
            flag_intermediate=1;

            if((i+1)<argc && argv[i+1][0]!='-'){
                intermediate_value = argv[i+1];
            }

            if(intermediate_value=="0"){
                intermediate=0;
            }
            else if(intermediate_value=="1"){
                intermediate=1;
            }
            else {
                invalid_intermediate=1;
            }
        }
        
        
        if(strncmp(argv[i],"--mode",6)==0){
            mode=argv[i+1];
            if(mode=="asm"){
                flag_m=2;
            }
            else if(mode=="raw"){
                flag_m=1;
            }
            else {
                cout<<"PLEASE INPUT A CORRECT DATA FORMAT MODE :("<<endl;
            }
        }
        
        if(strncmp(argv[i],"--input",7)==0){
            flag_inputf=1;
            file1.open(argv[i+1],ios::in | ios::binary);
            if (!file1.is_open())
            {
                cout <<"CANNOT OPEN INPUT FILE :("<< endl;
                continue;
            }
            
            inputF=argv[i+1];
        }
        if(strncmp(argv[i],"--output",7)==0){
            output=argv[i+1];
        }
        if(strncmp(argv[i],"--species",9)==0){
            SP=argv[i+1];
            cout <<"WE DO NOT NEED SPECIES PARAMETER RIGHT NOW :)"<< endl;
        }
        if(strncmp(argv[i],"--thread",8)==0){
            //NUM_threads=argv[i+1];
            //cout <<"WE DO NOT SUPPORT MULTITHREADS RIGHT NOW :)"<< endl;
            if(i+1<argc){
                int parsed_threads = atoi(argv[i+1]);
                if(parsed_threads<=0){
                    cout <<"PLEASE INPUT A POSITIVE NUMBER OF THREADS :("<< endl;
                }
                else {
                    NUM_threads = parsed_threads;
                }
            }
            else {
                cout <<"PLEASE ASSIGN A NUMBER AFTER --thread :("<< endl;
            }
        }
        if(strncmp(argv[i],"--ref_fa",8)==0){
            flag_reffa=1;
            file11.open(argv[i+1],ios::in | ios::binary);
            if (!file11.is_open())
            {
                cout <<"CANNOT OPEN REFERENCE FILE :("  << endl;
                cout <<"PLEASE ASSIGN A REFERENCE FILE."<< endl;
                //exit(1);
                flag_reffa=0;
            }
            ref_fa=argv[i+1];
        }
        
        if(strncmp(argv[i],"--custom_seq",12)==0){
            flag_cus=1;
            /*
            file12.open(argv[i+1],ios::in | ios::binary);
            if (!file12.is_open())
            {
                cout <<"CANNOT OPEN CUSTOMIZED FILE :("<< endl;
                cout<<"PLEASE ASSIGN A CUSTOMIZED FASTA FILE."<<endl;
                //exit(1);
                flag_cus=0;
            }*/
            cus=argv[i+1];
            
        }
        /*if(strncmp(argv[i],"--custom_index",14)==0){
            flag_cusin=1;
            if(flag_cus==1){
                 file13.open(argv[i+1],ios::in | ios::binary);
                 if (!file13.is_open())
                 {
                 cout <<"You are not assigning any index file for masking module."<< endl;
                 cout <<"Custmized finding will initiate without masking module."<<endl;
                 //exit(1);
                 flag_cusin=0;
                 }
            }
            cusin=argv[i+1];
        }*/
        if(strncmp(argv[i],"--TSD_finding",13)==0){
            
            tsd_pass=argv[i+1];
        }
        if(strncmp(argv[i],"--len_custom_seq",16)==0){
            s_cus_seq_len=argv[i+1];
            
            seq_len = std::stoi(s_cus_seq_len);
            //strstream sss;
            //sss << string_cus_seq_len;
            //sss >> cus_seq_len;
            //cus_seq_len=stoi(string_cus_seq_len)
            flag_cus_seq_len=1;
        }
        if(strncmp(argv[i],"--help",6)==0){
            
            help=1;
        }
    }
    
    if(T=="ALU"||T=="SVA"||T=="LINE"||T=="HERVK"){
        if(tsd_pass=="FALSE"){
            cout<<"MEI EVENTS REQUIRE TSD MOTIF FINDING. GO ON WITH TRUE ANYWAY."<<endl;
        }
        //else if(tsd_pass=="TRUE")
    }
    else if (T=="CUSTOMIZED"){
        flag_tsd=0;
        if(tsd_pass=="TRUE"){
            flag_tsd=1;
        }
        if(flag_tsd==1){
            cout<<"Customized insertion finding will go on with TSD finding module."<<endl;
        }
    }
    
    if(mode=="NA"){
        cout<<"A INPUT MODE FOR DATA FORMAT IS REQUIRED --mode."<<endl;
    }
    
    
    if((flag_T==0&&flag_cus==0)||(flag_T==1&&flag_cus==1)||flag_wd==0||flag_inputf==0||ref_n==0||flag_reffa==0||help==1||(flag_T==2&&flag_cus==0)||(flag_T==2&&flag_cus==1&&flag_tsd==1&flag_cus_seq_len==0)||(flag_m==0)||(mapq_int<0 || mapq_int>100)||invalid_intermediate==1){
        if(flag_T==0&&flag_cus==0){
            cout<<"***ERROR*** PLEASE ASSIGN A MEI TYPE! LINE/ALU/SVA/HERVK"<<endl;}
        if(flag_T==1&&flag_cus==1){
            cout<<"***ERROR*** PLEASE ASSIGN A MEI TYPE WITHOUT YOUR CUSTOMIZED SEQUENCE"<<endl;}
        if(flag_s!=flag_e){
            cout<<"***ERROR*** PLEASE ASSIGN START AND END AT THE SAME TIME"<<endl;}
        if(CHR=="ALL"&&(flag_s==1||flag_e==1)){
            cout<<"***ERROR*** CANNOT RUN WHOLE GENOME WITH ASSIGNING START OR END"<<endl;}
        if(flag_T==2&&flag_cus==0){
            cout<<"***ERROR*** PLEASE ASSIGN 'CUSTOMIZED' TYPE YOUR CUSTOMIZED SEQUENCE"<<endl;}
        if(flag_wd==0){
            cout<<"***ERROR*** PLEASE SET UP A WORKING DIRECTORY!"<<endl;}
        if(flag_inputf==0){
            cout<<"***ERROR*** PLEASE INPUT A FILE!"<<endl;}
        if(ref_n==0||flag_reffa==0){
            cout<<"***ERROR*** PLEASE ASSIGN A CORRECT REFERENCE version/fasta!"<<endl;}

        if(flag_T==2&&flag_cus==1&&flag_tsd==1&flag_cus_seq_len==0){
        cout<<"***ERROR*** PLEASE ASSIGN 'CUSTOMIZD_SEQUENCE_LENGTH' WHILE YOU ACTIVATE TSD_FINDING FOR YOUR CUSTOMIZED TYPE"<<endl;}
        if(flag_m==0){
            cout<<"***ERROR*** A INPUT MODE FOR DATA FORMAT IS REQUIRED"<<endl;}
        if(mapq_int<0 || mapq_int>100){
            cout<<"***ERROR*** A CORRECT RANGE OF MAPQ IS REQUIRED"<<endl;}
        if(invalid_intermediate==1){
            cout<<"***ERROR*** PLEASE SET --intermediate TO 0 (DELETE) OR 1 (KEEP)"<<endl;}
        cout<<endl;
        cout<<"***WELCOME***"<<endl;
        cout<<"***PALMER:Pre-mAsking Long reads for Mobile Element inseRtion***"<<endl;
        cout<<"Version: 2.0.0"<<endl;
        cout<<"Presented by Weichen Zhou @ Mills Lab. May.20th.2022"<<endl;
        cout<<endl;
        cout<<"USAGE:"<<endl;
        cout<<endl;
        cout<<"Required"<<endl;
        cout<<endl;
        cout<<"--input"<<endl;
        cout<<"         aligned long-read sequencing BAM file with directory path"<<endl;
        cout<<endl;
        cout<<"--workdir"<<endl;
        cout<<"         the user's working directory. Please follow the format /your/woking/directory/ !!don't forget the last '/'!!"<<endl;
        cout<<endl;
        cout<<"--ref_ver (options: hg19, GRCh37, GRCh38 or other)"<<endl;
        cout<<"         reference genome used for the aligned file ('other' option for the cusmized genome out of hg19, GRCh37 or GRCh38)"<<endl;
        cout<<endl;
        cout<<"--ref_fa"<<endl;
        cout<<"         indexed fasta file of reference genome fasta file with directory path used for the aligned bam file (wrong reference will cause error information)"<<endl;
        cout<<endl;
        cout<<"--type (options: LINE, ALU, SVA, HERVK, or CUSTOMIZED (if you want to setup your costomized sequence))"<<endl;
        cout<<"         type of MEIs or other kinds of insertions to detect"<<endl;
        cout<<endl;
        
        cout<<"--mode (options: raw, or asm)"<<endl;
        cout<<"         type of input sequencing to be processed (raw: raw nanopore/PacBio-sub reads; asm: assembled contigs)"<<endl;
        cout<<endl;
        
        cout<<"--chr (default: ALL (for whole genome); recommended options: chromosome1, chromosome2, ...chromosomeY)"<<endl;
        cout<<"         chromosome name for PALMER to run (if running for whole genome, don't need to assign). !!The chromosome names should be consistent with the ones in reference genome version!! e.g. for GRCh37, to run PALMER on chromosome1, the option should be '1', while for GRCh38 it should be 'chr1'"<<endl;
        cout<<endl;
        
        cout<<"Optional"<<endl;
        cout<<endl;
        
        cout<<"--start (default: Null)"<<endl;
        cout<<"         start position in the genome for PALMER to run (default is null). !!It should go with --end if assigned"<<endl;
        cout<<endl;
        
        cout<<"--end (default: Null)"<<endl;
        cout<<"         end position in the genome for PALMER to run (default is null). !!It should go with --start if assigned"<<endl;
        cout<<endl;
        
        cout<<"--custom_seq (default: Null)"<<endl;
        cout<<"         .fasta file with directory path to customize your insertion finding. e.g. NUMTs, MEIs in other species."<<endl;
        cout<<endl;
        //cout<<"--custom_index (default: Null; if you have both '--ref_ver other' and '--type LINE/ALU/SVA/HERVK', you must give PALMER a index file (format: \"CHR'\t'START'\t'END'\t'MEI_NAME'\n'\" for each MEI to be masked in each line) for masking module; if you have --custom_seq parameter without --custom_index, PALMER will work without masking step)"<<endl;
        //cout<<"         index file with directory path to mask the genome for your insertion finding"<<endl;
        cout<<endl;
        cout<<"--TSD_finding (Fixed: TRUE for all MEIs ,or default: FALSE for CUSTOMIZED insertion)"<<endl;
        cout<<"         whether to run TSD motif finding module for your insertion calling"<<endl;
        cout<<endl;
        

        cout<<"--mapq (default: MAPQ=10)"<<endl;
        cout<<"         the minimum MAPQ of the read for PALMER to process"<<endl;
        cout<<endl;
        
        cout<<"--thread (default: number of available hardware threads)"<<endl;
        cout<<"        number of concurrent region workers to launch during preprocessing and calling"<<endl;
        cout<<endl;
        
        cout<<"--len_custom_seq (MUST set up when activate TSD_finding for CUSTOMIZED insertion, otherwise CLOSED)"<<endl;
        cout<<"         interger value for the length of your customized sequence WITHOUT polyA tact"<<endl;
        cout<<endl;
        
        cout<<"--L_len (default: 25bp)"<<endl;
        cout<<"         the minimum length of putative LINE-1 aligned to L1.3 sequences"<<endl;
        cout<<endl;
        
        cout<<"--output (default: output)"<<endl;
        cout<<"         the prefix of the output file"<<endl;
        cout<<endl;
        
        cout<<"--intermediate (default: 0)"<<endl;
        cout<<"         set to 0 to delete intermediate subfolders after completion; set to 1 to retain them"<<endl;
        cout<<endl;
        
        exit(1);
    }
    

    cout<<"Variant type is "<<T<<endl;
    cout<<"Working directory is "<<WD<<endl;
    cout<<"Input file is "<<inputF<<endl;
    if(flag_m==2){
        cout<<"It is assigned as assembled contigs."<<endl;
    }
    else if(flag_m==1){
        cout<<"It is assigned as raw long-reads."<<endl;
    }
    
    cout<<"Output file is "<<WD<<output<<endl;
    if(flag_s==1&&flag_e==1){
        cout<<"Running on "<<CHR<<":"<<region_START<<"-"<<region_END<<endl;
    }
    else {
        cout<<"Running on "<<CHR<<endl;
    }
    cout<<"ref is "<<ref<<endl;
    
//Buildup & index or HG version chose
    string sys_dir="dirname "+dir;
    char *syst_dir=new char[sys_dir.length()+1];
    strcpy(syst_dir, sys_dir.c_str());
    
    vector<string> dir_conv;
    dir_conv.clear();
    FILE *pp =popen(syst_dir,"r");
    char tmp[1024];
    while (fgets(tmp, sizeof(tmp), pp) != NULL) {
        if (tmp[strlen(tmp) - 1] == '\n') {
            tmp[strlen(tmp) - 1] = '\0';
        }
        dir_conv.push_back(tmp);
    }
    pclose(pp);
    string direc;
    direc=accumulate(dir_conv.begin(),dir_conv.end(),direc);
    
    string buildup=direc+"/index/";
    
    string sys_region_index, sys_region_index_chr, sys_region_index_length;
    
   /* if(flag_cusin==1){
        sys_line_region=cusin;
    }
    else {
        sys_line_region="NULL";
    }*/
    
    
    if(T=="CUSTOMIZED"){
        T=cus;
    }
    
    //cout<<T<<endl;
    //cout<<seq_len<<endl;
    //cout<<flag_cus_seq_len<<endl;
    //getchar();
    
    if(ref_n==37){
        
        sys_region_index=buildup+"region.split.index.GRCh37";
        
        sys_region_index_chr=buildup+"region.split.index.chr.GRCh37";
        sys_region_index_length=buildup+"region.split.index.length.GRCh37";
        
        /*if(T=="LINE"){
            sys_line_region=buildup+"LINEs.regions.GRCh37";
        }
        if(T=="ALU"){
            sys_line_region=buildup+"Alu.regions.GRCh37";
        }
        if(T=="SVA"){
            sys_line_region=buildup+"SVA.regions.GRCh37";
        }
        if(T=="HERVK"){
            sys_line_region=buildup+"HERVK.regions.GRCh37";
        }*/
    }
    else if(ref_n==38){
        sys_region_index=buildup+"region.split.index.GRCh38";
        
        sys_region_index_chr=buildup+"region.split.index.chr.GRCh38";
        sys_region_index_length=buildup+"region.split.index.length.GRCh38";
        
        /*if(T=="LINE"){
            sys_line_region=buildup+"LINEs.regions.GRCh38";
        }
        if(T=="ALU"){
            sys_line_region=buildup+"Alu.regions.GRCh38";
        }
        if(T=="SVA"){
            sys_line_region=buildup+"SVA.regions.GRCh38";
        }
        if(T=="HERVK"){
            sys_line_region=buildup+"HERVK.regions.GRCh38";
        }*/
    }
    else if(ref_n==19){
        sys_region_index=buildup+"region.split.index.hg19";
        
        sys_region_index_chr=buildup+"region.split.index.chr.hg19";
        sys_region_index_length=buildup+"region.split.index.length.hg19";
        
        /*if(T=="LINE"){
            sys_line_region=buildup+"LINEs.regions.hg19";
        }
        if(T=="ALU"){
            sys_line_region=buildup+"Alu.regions.hg19";
        }
        if(T=="SVA"){
            sys_line_region=buildup+"SVA.regions.hg19";
        }
        if(T=="HERVK"){
            sys_line_region=buildup+"HERVK.regions.hg19";
        }*/
    }
    
    
//reference_index buildup
    else if(ref_n==-1){
        //sys_line_region=cusin;
        
        
        buildup=WD+"index/";
        
        string sys_build;
        sys_build="mkdir "+buildup;
        
        char *syst_build = new char[sys_build.length()+1];
        strcpy(syst_build, sys_build.c_str());
        system(syst_build);
        
        string sys1;
        sys1="samtools view "+inputF+" -H |grep \"@SQ\" | awk -F \":|\t\" '{ print $3}' >"+buildup+"chr.list";
        
        string sys2;
        sys2="samtools view "+inputF+" -H |grep \"@SQ\" | awk -F \":|\t\" '{ print $5}' >"+buildup+"length.list";
        
        char *syst1 = new char[sys1.length()+1];
        strcpy(syst1, sys1.c_str());
        system(syst1);
        
        char *syst2 = new char[sys2.length()+1];
        strcpy(syst2, sys2.c_str());
        system(syst2);
        
        sys_region_index_chr=buildup+"chr.list";
        sys_region_index_length=buildup+"length.list";
        sys_region_index=buildup+"region.split.index";
    }
    //cout<<sys_region_index<<" "<<sys_line_region<<endl;
 //original
    char *syst_region_index_chr =new char[sys_region_index_chr.length()+1];
    strcpy(syst_region_index_chr, sys_region_index_chr.c_str());
    char *syst_region_index_length =new char[sys_region_index_length.length()+1];
    strcpy(syst_region_index_length, sys_region_index_length.c_str());
    
    ifstream file91;
    ifstream file92;
    file91.open(syst_region_index_chr);
    file92.open(syst_region_index_length);
    
    int line_chr;
    int line_len;
    string input_inde;
    for(int i=0;!file91.eof();++i){
        getline(file91,input_inde);
        line_chr=i;
    }
    for(int i=0;!file92.eof();++i){
        getline(file92,input_inde);
        line_len=i;
    }
    
    //cout<<"line_chr="<<line_chr<<endl;
    //cout<<"line_length="<<line_len<<endl;
    
    file92.close();
    file91.close();
    file92.clear();
    file91.clear();
    
    file91.open(syst_region_index_chr);
    file92.open(syst_region_index_length);
    
    string chr_inde[line_chr];
    int len_inde[line_len];
    
    for(int i=0;i!=line_chr;++i){
        file91>>chr_inde[i];
    }
    for(int i=0;i!=line_len;++i){
        file92>>len_inde[i];
    }
    
    file92.close();
    file91.close();
    file92.clear();
    file91.clear();
    ifstream file2;
    
    //char *syst_region_index ;
    char *syst_region_index =new char[sys_region_index.length()+1];
    strcpy(syst_region_index, sys_region_index.c_str());
    
    
    if(ref_n==-1){
        
        ofstream file93;
        file93.open(syst_region_index,ios::trunc);
        
        int bin=1000000;
        for(int i=0;i!=line_chr;++i){
            int j=i;
            int sec;
                //last;
            sec=int(len_inde[j]/bin);
                //last=len_inde[j]%bin;
            for(int k=0;k!=sec;k++){
                file93<<chr_inde[i]<<'\t'<<(1+k*bin)<<'\t'<<(k*bin+bin)<<endl;
            }
            file93<<chr_inde[i]<<'\t'<<(1+sec*bin)<<'\t'<<len_inde[j]<<endl;
            
        }
        //ifstream file2;
        file93.close();
        file93.clear();
    }
    
    else {
    //string sys_region_index=buildup+"region.split.index";
        
        file2.open(syst_region_index);
        
        if (!file2.is_open())
        {
            cout <<"CANNOT OPEN REF INDEX FILE"<< endl;
            //cout<<"YEA"<<endl;
            exit(1);
        }
        file2.close();
        file2.clear();
        //file2.open(syst_region_index);
    }
    
    
//parameters_end
    
//multiple threads
    //int input_int1;
    //int input_int2;
    int region_s=0;
    int region_e=0;
    if(flag_s==1&&flag_e==1){
        region_s=stoi(region_START);
        region_e=stoi(region_END);
    }

    struct Region {
        string chr;
        int start;
        int end;
    };

    vector<Region> regions;
    regions.reserve(1024);
    
    file2.open(syst_region_index);

    string chr;
    int start, end;
    while(file2>>chr>>start>>end){
        bool include_region=false;
        if(CHR=="ALL"){
                    include_region=true;
                }
        else if(CHR==chr){
            if(flag_s==1&&flag_e==1){
                if(!((region_s>end)||(region_e<start))){
                    include_region=true;
                }
            }
            else if(flag_s==0&&flag_e==0){
                include_region=true;
            }
        }
        
        if(include_region){
            regions.push_back({chr,start,end});
        }
    }

        
    /*
    string input_index;
    int line_index=0;
    for(int i=1;!file2.eof();){
        //file2>>input_index;
        file2>>input_index;
        if(CHR=="ALL"){
            line_index=i;
            ++i;
            file2>>input_int1;
            file2>>input_int2;
        }
        else if(CHR==input_index&&flag_s==0&&flag_e==0){
            line_index=i;
            ++i;
            file2>>input_int1;
            file2>>input_int2;
        }
        else if(CHR==input_index&&flag_s==1&&flag_e==1){
            file2>>input_int1;
            file2>>input_int2;
            if(!((region_s>input_int2)||(region_e<input_int1))){
                line_index=i;
                ++i;
            }
        }
        
    }*/
    
    
    
    
    
    //line_index=line_index+1;
    
    //cout<<"THERE ARE "<<line_index<<" REGIONS TO COUNT."<<endl;
    //cout<<"Pre-masking step & single read calling step is initiated."<<endl;
    
    file2.close();
    file2.clear();
    //file2.open(syst_region_index);
    
    //int NUM_circle;
    //NUM_circle=(line_index/(NUM_threads+1))+1;

    //pid_t p1, p2, p3, p4, p5, p6, p7, p8, p9, p10;
    //, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30;
 
    //ADDITONAL
    size_t line_index = regions.size();
    cout<<"THERE ARE "<<line_index<<" REGIONS TO COUNT."<<endl;
    if(line_index==0){
        cout<<"NO REGION MATCHES THE GIVEN PARAMETERS."<<endl;
    }
    else{
        cout<<"Pre-masking step & single read calling step is initiated."<<endl;
    }
    
    if(line_index>0){
        unsigned int worker_count = static_cast<unsigned int>(NUM_threads);
        if(worker_count==0){
            worker_count=1;
        }
        if(worker_count>regions.size()){
            worker_count=static_cast<unsigned int>(regions.size());
        }
        
        cout<<"Launching processing with "<<worker_count<<" thread(s)."<<endl;
        
        std::atomic<size_t> next_region(0);
        auto worker=[&](){
            while(true){
                size_t idx = next_region.fetch_add(1);
                if(idx>=regions.size()){
                    break;
                }
                const Region &region = regions[idx];
                tube(WD, inputF, region.chr, region.start, region.end, T, ref_n, direc, ref_fa, flag_tsd, L_len, seq_len, mode, mapq_int, intermediate);
            }
        };
        
        vector<std::thread> threads;
        threads.reserve(worker_count);
        for(unsigned int i=0;i<worker_count;++i){
            threads.emplace_back(worker);
        }
        for(auto &thread:threads){
            thread.join();
        }
    }
    
    
    /*
    string input;
    string chr;
    int start, end;

//no multithread
    for(int i=0;i!=line_index;){
        //cout<<"right call"<<endl;
        //cout<<chr<<endl;
        
        file2>>chr;
        file2>>start;
        file2>>end;
        int chr_index=0;
        if(CHR=="ALL"){
            chr_index=1;
        }
        else if(CHR==chr&&flag_s==0&&flag_e==0){
            chr_index=1;
            
        }
        else if(CHR==chr&&flag_s==1&&flag_e==1){
            if(!((region_s>end)||(region_e<start))){
                chr_index=1;
            }
        }
        
        
        if(chr_index==1){
            ++i;
            
            //getchar();
            //****
            //cout<<flag_tsd<<endl;
            tube(WD, inputF, chr, start, end, T, ref_n, direc, ref_fa, flag_tsd, L_len, seq_len, mode, mapq_int);
        }
        /*ver1.3
        if(chr_index==1){
            tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n, direc, ref_file);
        }
    }*/
    
   
//merge and calling

    cout<<"Merging step is initiated."<<endl;
    //mkdir
    /*
    for(int i=0;i!=line_chr;++i){
        string WD_chr=WD+chr_inde[i]+"/";
        string WD_chr="mkdir "+WD+chr_inde[i]+"/";
        
        if(CHR==chr_inde[i]||CHR=="ALL"){
            char *syst_WD_chr =new char[sys_WD_chr.length()+1];
            strcpy(syst_WD_chr, sys_WD_chr.c_str());
            system(syst_WD_chr);
        }
    }
    */
    /*
    file2.close();
    file2.clear();
    file2.open(syst_region_index);
    */
    //merge
    
    string sys_final_title = WD+output+"_calls.txt";
    char *syst_final_title = new char[sys_final_title.length()+1];
    strcpy(syst_final_title, sys_final_title.c_str());
    ofstream file3;
    file3.open(syst_final_title,ios::trunc);
    
    file3<<"cluster_id"<<'\t'<<"chr"<<'\t'<<"start1"<<'\t'<<"start2"<<'\t'<<"end1"<<'\t'<<"end2"<<'\t'<<"start1_inVariant"<<'\t'<<"start2_inVariant"<<'\t'<<"end1_inVariant"<<'\t'<<"end2_inVariant"<<'\t'<<"Confident_supporting_reads"<<'\t'<<"Potential_supporting_reads"<<'\t'<<"Ptential_segmental_supporting_reads"<<'\t'<<"orientation"<<'\t'<<"polyA-tail_size"<<'\t'<<"5'_TSD_size"<<'\t'<<"3'_TSD_size"<<'\t'<<"Predicted_transD_size"<<'\t'<<"Has_5'_inverted_sequence?"<<'\t'<<"5'_inverted_seq_end"<<'\t'<<"5'_seq_start"<<endl;
    
    string sys_final_tsd_title = WD+output+"_TSD_reads.txt";
    char *syst_final_tsd_title = new char[sys_final_tsd_title.length()+1];
    strcpy(syst_final_tsd_title, sys_final_tsd_title.c_str());
    ofstream file31;
    file31.open(syst_final_tsd_title,ios::trunc);
    
    file31<<"cluster_id"<<'\t'<<"read_name.info"<<'\t'<<"5'_TSD"<<'\t'<<"3'_TSD"<<'\t'<<"Predicted_transD"<<'\t'<<"Unique_26mer_at_5'junction"<<'\t'<<"Whole_insertion_seq"<<endl;
    
    /*
    for(int i=0;i!=line_index;){
        file2>>chr;
        file2>>start;
        file2>>end;
        
        stringstream ss1, ss2;
        ss1 << start;
        string s_start =ss1.str();
        ss2 << end;
        string s_end =ss2.str();
        
        if(CHR=="ALL"){
            string sys_final="cat "+WD+chr+"_"+s_start+"_"+s_end+"/calls.txt >> "+sys_final_title;
            //cout<<sys_final<<endl;
            char *syst_final = new char[sys_final.length()+1];
            strcpy(syst_final, sys_final.c_str());
            system(syst_final);
            
            string sys_final_tsd="cat "+WD+chr+"_"+s_start+"_"+s_end+"/TSD_output.txt >> "+sys_final_tsd_title;
            //cout<<sys_final<<endl;
            char *syst_final_tsd = new char[sys_final_tsd.length()+1];
            strcpy(syst_final_tsd, sys_final_tsd.c_str());
            system(syst_final_tsd);
            ++i;
        }
        
        else if(CHR==chr&&flag_s==0&&flag_e==0){
            string sys_final="cat "+WD+chr+"_"+s_start+"_"+s_end+"/calls.txt >> "+sys_final_title;
            //cout<<sys_final<<endl;
            char *syst_final = new char[sys_final.length()+1];
            strcpy(syst_final, sys_final.c_str());
            system(syst_final);
            
            string sys_final_tsd="cat "+WD+chr+"_"+s_start+"_"+s_end+"/TSD_output.txt >> "+sys_final_tsd_title;
            //cout<<sys_final<<endl;
            char *syst_final_tsd = new char[sys_final_tsd.length()+1];
            strcpy(syst_final_tsd, sys_final_tsd.c_str());
            system(syst_final_tsd);
            ++i;
        }
        
        else if(CHR==chr&&flag_s==1&&flag_e==1){
            if(!((region_s>end)||(region_e<start))){
                string sys_final="cat "+WD+chr+"_"+s_start+"_"+s_end+"/calls.txt >> "+sys_final_title;
                //cout<<sys_final<<endl;
                char *syst_final = new char[sys_final.length()+1];
                strcpy(syst_final, sys_final.c_str());
                system(syst_final);
                
                string sys_final_tsd="cat "+WD+chr+"_"+s_start+"_"+s_end+"/TSD_output.txt >> "+sys_final_tsd_title;
                //cout<<sys_final<<endl;
                char *syst_final_tsd = new char[sys_final_tsd.length()+1];
                strcpy(syst_final_tsd, sys_final_tsd.c_str());
                system(syst_final_tsd);
                ++i;
            }
        }
    }*/
    
    for(const auto &region:regions){
        stringstream ss1, ss2;
        ss1 << region.start;
        string s_start =ss1.str();
        ss2 << region.end;
        string s_end =ss2.str();
        
        string sys_final="cat "+WD+region.chr+"_"+s_start+"_"+s_end+"/calls.txt >> "+sys_final_title;
        char *syst_final = new char[sys_final.length()+1];
        strcpy(syst_final, sys_final.c_str());
        system(syst_final);
        
        string sys_final_tsd="cat "+WD+region.chr+"_"+s_start+"_"+s_end+"/TSD_output.txt >> "+sys_final_tsd_title;
        char *syst_final_tsd = new char[sys_final_tsd.length()+1];
        strcpy(syst_final_tsd, sys_final_tsd.c_str());
        system(syst_final_tsd);
    }
    

    cout<<"Merging step completed."<<endl;

    string sys_final_vcf = WD+output+"_integrated.vcf";
    string command_line;
    for (int arg_i = 0; arg_i < argc; ++arg_i) {
        if (arg_i > 0) {
            command_line += " ";
        }
        command_line += argv[arg_i];
    }

    write_vcf_output(sys_final_title, sys_final_tsd_title, sys_final_vcf, T, ref_fa, command_line);

    //colapse for the redundant calls

    //kmer



    vector<string> intermediate_dirs;
    file2.close();
    file2.clear();
    file2.open(syst_region_index);

    if(file2.is_open()){
        for(int i=0;i!=line_index;){
            file2>>chr;
            file2>>start;
            file2>>end;

            int chr_index=0;
            if(CHR=="ALL"){
                chr_index=1;
            }
            else if(CHR==chr&&flag_s==0&&flag_e==0){
                chr_index=1;

            }
            else if(CHR==chr&&flag_s==1&&flag_e==1){
                if(!((region_s>end)||(region_e<start))){
                    chr_index=1;
                }
            }

            if(chr_index==1){
                ++i;
                stringstream ss1, ss2;
                ss1 << start;
                ss2 << end;
                string s_start = ss1.str();
                string s_end = ss2.str();
                string dir_path = WD+chr+"_"+s_start+"_"+s_end+"/";
                intermediate_dirs.push_back(dir_path);
            }
        }
    }

    file2.close();
    file2.clear();

    if(intermediate==0){
        cout<<"Removing intermediate subfolders. Use --intermediate 1 to retain them."<<endl;
        for(const auto& dir_path : intermediate_dirs){
            string remove_cmd = "rm -rf \"" + dir_path + "\"";
            char *syst_remove_cmd = new char[remove_cmd.length()+1];
            strcpy(syst_remove_cmd, remove_cmd.c_str());
            system(syst_remove_cmd);
            delete [] syst_remove_cmd;
        }
    } else {
        cout<<"Intermediate subfolders retained per --intermediate setting."<<endl;
    }

    cout<<"Final calls finished."<<endl;
    cout<<"Results are in "+WD+output<<endl;


    
}
