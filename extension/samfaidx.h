#ifndef PALMER_SamFaidx_H
#define PALMER_SamFaidx_H

#include <config.h>

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/hfile.h>
#include <htslib/kstring.h>

#include <samtools/samtools.h>

// #include "extension/extern-sam.h"

class SamFaidx
{
private:

public:

public:
    SamFaidx(/* args */);
    ~SamFaidx();

    int SamFaidxCommand(const char *faFileName, const char *regions, const char* output_file);
private:
};

#endif