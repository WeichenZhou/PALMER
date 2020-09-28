
#ifndef PALMER_Extern_Sam_H
#define PALMER_Extern_Sam_H
extern "C"
{
#include "samtools.h"
#include "sam_opts.h"
#include "bedidx.h"
#include "version.h"

    KHASH_SET_INIT_STR(rg)
    KHASH_SET_INIT_STR(tv)

    typedef khash_t(rg) * rghash_t;
    typedef khash_t(tv) * tvhash_t;

    // This structure contains the settings for a samview run
    typedef struct samview_settings
    {
        rghash_t rghash;
        tvhash_t tvhash;
        int min_mapQ;
        int flag_on;
        int flag_off;
        int flag_alloff;
        int min_qlen;
        int remove_B;
        uint32_t subsam_seed;
        double subsam_frac;
        char *library;
        void *bed;
        size_t remove_aux_len;
        char **remove_aux;
        int multi_region;
        char *tag;
    } samview_settings_t;

    extern const char *bam_get_library(sam_hdr_t *header, const bam1_t *b);
    extern int bam_remove_B(bam1_t *b);
    extern char *samfaipath(const char *fn_ref);
}

#endif