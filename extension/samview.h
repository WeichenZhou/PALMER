#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/thread_pool.h"
#include "htslib/header.h"

#include "extension/extern-sam.h"
#include "extension/samline.h"

#ifndef PALMER_SamView_H
#define PALMER_SamView_H

class Samview
{
private:
    samview_settings_t settings;
    sam_global_args ga;

public:
    std::vector<SamLine> regionLines;
    std::vector<string> headerChr;
    std::vector<int> headerLength;

public:
    Samview(/* args */);
    ~Samview();

    int SamViewCommand(int argc, char *argv[], const char *inFileName, int argMinMapQ, const char *argRegion);
    // int SamViewCommand(const char *inFileName, int argMinMapQ, const char *argRegion);
    int SamViewHeaderOnly(const char *inFileName);

private:
    int check_sam_write1(const sam_hdr_t *h, const bam1_t *b);
    int process_aln_palmer(const sam_hdr_t *h, bam1_t *b, samview_settings_t *settings);
};

#endif