#include "samview.h"

Samview::Samview(/* args */)
{
    this->settings = {
        rghash : NULL,
        tvhash : NULL,
        min_mapQ : 0,
        flag_on : 0,
        flag_off : 0,
        flag_alloff : 0,
        min_qlen : 0,
        remove_B : 0,
        subsam_seed : 0,
        subsam_frac : -1.,
        library : NULL,
        bed : NULL,
        remove_aux_len : 0,
        remove_aux : NULL,
        multi_region : 0,
        tag : NULL
    };
    this->ga = {
        {unknown_category, unknown_format, 0, 0, no_compression, 0},
        {unknown_category, unknown_format, 0, 0, no_compression, 0},
        0,
        0,
        0};
}

Samview::~Samview()
{
    printf("this->headerChr.clear();\n");
    this->headerChr.clear();
    printf("this->headerLength.clear();\n");
    this->headerLength.clear();
    printf("this->regionLines.clear();\n");
    this->regionLines.clear();
    // close files, free and return
    printf("free(settings.library);\n");
    free(settings.library);
    if (settings.bed)
    {
        printf("bed_destroy(settings.bed);\n");
        bed_destroy(settings.bed);
    }
    if (settings.rghash)
    {
        printf("1 free((char *)kh_key(settings.rghash, k));\n");
        khint_t k;
        for (k = 0; k < kh_end(settings.rghash); ++k)
            if (kh_exist(settings.rghash, k))
                free((char *)kh_key(settings.rghash, k));
        printf("1 kh_destroy(rg, settings.rghash);\n");
        kh_destroy(rg, settings.rghash);
    }
    if (settings.tvhash)
    {
        printf("2 free((char *)kh_key(settings.tvhash, k));\n");
        khint_t k;
        for (k = 0; k < kh_end(settings.tvhash); ++k)
            if (kh_exist(settings.tvhash, k))
                free((char *)kh_key(settings.tvhash, k));
        printf("2 kh_destroy(tv, settings.tvhash);\n");
        kh_destroy(tv, settings.tvhash);
    }
    if (settings.remove_aux_len)
    {
        printf("free(settings.remove_aux);\n");
        free(settings.remove_aux);
    }
    if (settings.tag)
    {
        printf("free(settings.tag);\n");
        free(settings.tag);
    }
    printf("sam_global_args_free(&ga);\n");
    sam_global_args_free(&ga);
}

int Samview::check_sam_write1(const sam_hdr_t *h, const bam1_t *b)
{
    // NewCode
    kstring_t kstr = KS_INITIALIZE;
    printf("sam_format1(h, b, &kstr);\n");
    int ret = sam_format1(h, b, &kstr);
    SamLine sl = SamLine(kstr.s);
    regionLines.push_back(sl);
    // printf("sam_write:: %d\n", ret);

    if (kstr.s)
    {
        printf("free(kstr.s);\n");
        free(kstr.s);
    }

    return ret;
}

int Samview::process_aln_palmer(const sam_hdr_t *h, bam1_t *b, samview_settings_t *settings)
{
    if (settings->remove_B){
        printf("bam_remove_B(b);\n");
        bam_remove_B(b);
    }
    if (settings->min_qlen > 0)
    {
        int k, qlen = 0;
        uint32_t *cigar = bam_get_cigar(b);
        for (k = 0; k < b->core.n_cigar; ++k)
            if ((bam_cigar_type(bam_cigar_op(cigar[k])) & 1) || bam_cigar_op(cigar[k]) == BAM_CHARD_CLIP)
                qlen += bam_cigar_oplen(cigar[k]);
        if (qlen < settings->min_qlen)
            return 1;
    }
    if (b->core.qual < settings->min_mapQ || ((b->core.flag & settings->flag_on) != settings->flag_on) || (b->core.flag & settings->flag_off))
        return 1;
    if (settings->flag_alloff && ((b->core.flag & settings->flag_alloff) == settings->flag_alloff))
        return 1;
    if (!settings->multi_region && settings->bed && (b->core.tid < 0 || !bed_overlap(settings->bed, sam_hdr_tid2name(h, b->core.tid), b->core.pos, bam_endpos(b))))
        return 1;
    if (settings->subsam_frac > 0.)
    {
        uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(b)) ^ settings->subsam_seed);
        if ((double)(k & 0xffffff) / 0x1000000 >= settings->subsam_frac)
            return 1;
    }
    if (settings->rghash)
    {
        uint8_t *s = bam_aux_get(b, "RG");
        if (s)
        {
            khint_t k = kh_get(rg, settings->rghash, (char *)(s + 1));
            if (k == kh_end(settings->rghash))
                return 1;
        }
    }
    if (settings->tvhash && settings->tag)
    {
        uint8_t *s = bam_aux_get(b, settings->tag);
        if (s)
        {
            khint_t k = kh_get(tv, settings->tvhash, (char *)(s + 1));
            if (k == kh_end(settings->tvhash))
                return 1;
        }
        else
        {
            return 1;
        }
    }
    if (settings->library)
    {
        const char *p = bam_get_library((sam_hdr_t *)h, b);
        if (!p || strcmp(p, settings->library) != 0)
            return 1;
    }
    if (settings->remove_aux_len)
    {
        size_t i;
        for (i = 0; i < settings->remove_aux_len; ++i)
        {
            uint8_t *s = bam_aux_get(b, settings->remove_aux[i]);
            if (s)
            {
                bam_aux_del(b, s);
            }
        }
    }
    return 0;
}

int Samview::SamViewHeaderOnly(const char *inFileName)
{
    int ret = 0;

    sam_hdr_t *header = NULL;
    samFile *in = 0;
    char *fn_in = 0;

    fn_in = strdup(inFileName);

    if ((in = sam_open_format(fn_in, "r", &ga.in)) == 0)
    {
        print_error_errno("[SamView]", "failed to open \"%s\" for reading", fn_in);
        return 1;
    }
    if ((header = sam_hdr_read(in)) == 0)
    {
        fprintf(stderr, "[SamView] fail to read the header from \"%s\".\n", fn_in);
        return 1;
    }

    string str = string(header->text);
    size_t pos = 0;
    size_t pos_SQ = 0;
    size_t pos_lineEnd = 0;
    size_t i, j;

    this->headerChr.clear();
    this->headerLength.clear();
    while (true)
    {
        pos_SQ = str.find("@SQ", pos);
        if (pos_SQ == string::npos)
        {
            break;
        }

        pos_lineEnd = str.find("\n", pos_SQ);
        if (pos_lineEnd == string::npos)
        {
            pos_lineEnd = str.length();
        }
        pos = pos_lineEnd;

        string SQLine = str.substr(pos_SQ, pos_lineEnd - pos_SQ);
        i = SQLine.find("\tSN:");
        j = SQLine.find("\tLN:");
        this->headerChr.push_back(SQLine.substr(i + 4, j - (i + 4)));
        this->headerLength.push_back(std::stoi(SQLine.substr(j + 4, SQLine.length())));
    }

    if (in)
        check_sam_close("view", in, fn_in, "standard input", &ret);
    if (fn_in)
        free(fn_in);
    if (header)
        sam_hdr_destroy(header);

    return ret;
}

int Samview::SamViewCommand(int argc, char *argv[], const char *inFileName, int argMinMapQ, const char *argRegion)
// int Samview::SamViewCommand(const char *inFileName, int argMinMapQ, const char *argRegion)
{
    this->regionLines.clear();

    int ret = 0;
    int64_t count = 0;
    samFile *in = 0;
    sam_hdr_t *header = NULL;
    char *fn_in = 0;
    char *arg_list = NULL;

    /* fake argument parse */
    // -F
    settings.flag_off |= strtol("0x100", 0, 0);
    settings.flag_off |= strtol("0x200", 0, 0);
    settings.flag_off |= strtol("0x400", 0, 0);
    settings.flag_off |= strtol("0x800", 0, 0);
    // -q 10
    settings.min_mapQ = argMinMapQ;

    fn_in = strdup(inFileName);

    // open file handlers
    if ((in = sam_open_format(fn_in, "r", &ga.in)) == 0)
    {
        print_error_errno("view", "failed to open \"%s\" for reading", fn_in);
        ret = 1;
        goto view_end;
    }

    if ((header = sam_hdr_read(in)) == 0)
    {
        fprintf(stderr, "[samtools-view] fail to read the header from \"%s\".\n", fn_in);
        ret = 1;
        goto view_end;
    }
    if (settings.rghash)
    {
        sam_hdr_remove_lines(header, "RG", "ID", settings.rghash);
    }

    {
        if (!(arg_list = stringify_argv(argc + 1, argv - 1)))
        {
            print_error("view", "failed to create arg_list");
            ret = 1;
            goto view_end;
        }

        if (arg_list)
        {
            printf("############# %s \n", arg_list);
        }

        if (sam_hdr_add_pg(header, "samtools",
                           "VN", SAMTOOLS_VERSION,
                           arg_list ? "CL" : NULL,
                           arg_list ? arg_list : NULL,
                           NULL))
        {
            print_error("view", "failed to add PG line to the header");
            ret = 1;
            goto view_end;
        }
    }

    {
        // retrieve alignments in specified regions
        // int i;
        bam1_t *b;
        hts_idx_t *idx = NULL;
        {
            idx = sam_index_load(in, fn_in);
        }
        if (idx == 0)
        { // index is unavailable
            fprintf(stderr, "[samtools-view] random alignment retrieval only works for indexed BAM or CRAM files.\n");
            ret = 1;
            goto view_end;
        }
        b = bam_init1();
        {
            int result;
            hts_itr_t *iter = sam_itr_querys(idx, header, argRegion); // parse a region in the format like `chr2:100-200'
            if (iter == NULL)
            { // region invalid or reference name not found
                fprintf(stderr, "[samtools-view] region \"%s\" specifies an invalid region or unknown reference. Continue anyway.\n", argRegion);
            }
            // fetch alignments
            while ((result = sam_itr_next(in, iter, b)) >= 0)
            {
                if (!this->process_aln_palmer(header, b, &settings))
                {
                    {
                        if (this->check_sam_write1(header, b) < 0)
                            break;
                    }
                    count++;
                }
                else
                {
                }
            }
            hts_itr_destroy(iter);
            if (result < -1)
            {
                fprintf(stderr, "[samtools-view] retrieval of region \"%s\" failed due to truncated file or corrupt BAM index file\n", argRegion);
                ret = 1;
            }
        }
        bam_destroy1(b);
        hts_idx_destroy(idx); // destroy the BAM index
    }

view_end:
    // close files, free and return
    if (in)
        check_sam_close("view", in, fn_in, "standard input", &ret);
    if (fn_in)
        free(fn_in);

    if (header)
        sam_hdr_destroy(header);

    if (arg_list)
        free(arg_list);

    return ret;
}
