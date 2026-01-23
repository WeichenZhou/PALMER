////copyright by ArthurZhou @ UMich&Fudan&HUST
#include "common.hpp"
#include <htslib/kstring.h>
#include <htslib/sam.h>

int samtools(string working_dir, string input_bam, string chr, string start, string end, string fasta, string MAPQ){

    string region = chr + ":" + start + "-" + end;
    string output_path = working_dir + "region.sam";

    ofstream out_stream(output_path);
    if (!out_stream.is_open()) {
        return -1;
    }

    htsFile *in = sam_open(input_bam.c_str(), "r");
    if (!in) {
        return -1;
    }
    if (hts_set_fai_filename(in, fasta.c_str()) != 0) {
        sam_close(in);
        return -1;
    }

    bam_hdr_t *header = sam_hdr_read(in);
    if (!header) {
        sam_close(in);
        return -1;
    }

    hts_idx_t *idx = sam_index_load(in, input_bam.c_str());
    if (!idx) {
        bam_hdr_destroy(header);
        sam_close(in);
        return -1;
    }

    hts_itr_t *itr = sam_itr_querys(idx, header, region.c_str());
    if (!itr) {
        hts_idx_destroy(idx);
        bam_hdr_destroy(header);
        sam_close(in);
        return -1;
    }

    bam1_t *record = bam_init1();
    kstring_t formatted = {0, 0, nullptr};
    int mapq_threshold = atoi(MAPQ.c_str());

    while (sam_itr_next(in, itr, record) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) {
            continue;
        }

        if (record->core.qual < mapq_threshold) {
            continue;
        }

        if (sam_format1(header, record, &formatted) < 0) {
            break;
        }

        for (size_t i = 0; i < formatted.l; ++i) {
            if (formatted.s[i] == ' ') {
                formatted.s[i] = '_';
            }
        }

        out_stream.write(formatted.s, formatted.l);
        out_stream.put('\n');
    }

    free(formatted.s);
    bam_destroy1(record);
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(in);

    return out_stream ? 0 : -1;
}
