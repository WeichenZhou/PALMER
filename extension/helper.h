#include <ctime>

#include "extension/extern-samview.h"

#ifndef PALMER_HELPER_H
#define PALMER_HELPER_H

int XPrintSettings(samview_settings_t *settings)
{
    char *prefix = "[XPrintSettings]";
    printf("%s settings->min_mapQ = %d\n", prefix, settings->min_mapQ);
    printf("%s settings->flag_on = %d\n", prefix, settings->flag_on);
    printf("%s settings->flag_off = %d\n", prefix, settings->flag_off);
    printf("%s settings->flag_alloff = %d\n", prefix, settings->flag_alloff);
    printf("%s settings->min_qlen = %d\n", prefix, settings->min_qlen);
    printf("%s settings->remove_B = %d\n", prefix, settings->remove_B);
    printf("%s settings->subsam_seed = %d\n", prefix, settings->subsam_seed);
    printf("%s settings->subsam_frac = %f\n", prefix, settings->subsam_frac);
    printf("%s settings->library = %s\n", prefix, settings->library);
    printf("%s settings->bed = %p\n", prefix, settings->bed);
    printf("%s settings->remove_aux_len = %ld\n", prefix, settings->remove_aux_len);
    printf("%s settings->multi_region = %d\n", prefix, settings->multi_region);
    printf("%s settings->tag = %s\n", prefix, settings->tag);

    return 0;
}

class Ticker
{
private:
    clock_t start;
    clock_t end;

public:
    Ticker(/* args */) { start = clock(); };
    ~Ticker(){};
    int Start()
    {
        start = clock();
        // printf("[Ticker] start = %ld s\n", start);
    };
    int End()
    {
        end = clock();
        // printf("[Ticker] end = %ld s\n", end);
        // printf("[Ticker] CLOCKS_PER_SEC = %ld s\n", CLOCKS_PER_SEC);
        double time_duration = end - start;
        printf("[Ticker] time_duration = %f s\n", time_duration / CLOCKS_PER_SEC);
    };
};

#endif