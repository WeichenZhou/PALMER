////copyright by ArthurZhou @ UMich&Fudan&HUST
#include "common.hpp"

int samtools(string working_dir, string input_bam, string chr, string start, string end, string fasta, string MAPQ){
    string output_path = working_dir + "region.sam";
    ofstream region_stream(output_path);
    if (!region_stream.is_open()) {
        cout << "CANNOT WRITE FILE, 'region.sam'" << endl;
        return 1;
    }

    vector<string> args = {"samtools", "view", "-q", MAPQ, "-F", "0x100", "-F", "0x200",
                            "-F", "0x400", "-T", fasta, input_bam, chr + ":" + start + "-" + end};

    bool success = stream_process_output(args, [&](const char *buffer, ssize_t count) {
        bool in_space = false;
        for (ssize_t i = 0; i < count; ++i) {
            char c = buffer[i];
            if (c == ' ') {
                if (!in_space) {
                    region_stream.put('_');
                    in_space = true;
                }
            } else {
                region_stream.put(c);
                in_space = false;
            }
        }
        return true;
    });

    region_stream.close();
    if (!success) {
        remove(output_path.c_str());
        return 1;
    }

    return 0;
}
