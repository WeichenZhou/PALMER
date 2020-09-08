////copyright by ArthurZhou @ UMich&FUDAN&HUST
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

#include "extension/samview.h"

#define ARG_SIZE 12

using namespace std;

int samtools(string working_dir, string input_bam, string chr, string start, string end, Samview *samview)
{
    int argc = ARG_SIZE;
    char *argv[ARG_SIZE] = {
        "samtools",
        "view",
        "-q",
        "10",
        "-F",
        "0x100",
        "-F",
        "0x200",
        "-F",
        "0x400",
    };
    string region = (chr + ":" + start + "-" + end);
    argv[ARG_SIZE - 2] = (char *)malloc(sizeof(input_bam));
    argv[ARG_SIZE - 1] = (char *)malloc(sizeof(region));
    strcpy(argv[ARG_SIZE - 2], input_bam.c_str());
    strcpy(argv[ARG_SIZE - 1], region.c_str());

    samview->SamViewCommand(argc - 1, argv + 1, input_bam.c_str(), 10, region.c_str());
    if (samview->regionLines.size() > 0)
    {
        cout << "samview.regionLines.size() = " << samview->regionLines.size() << endl;
    }

    return 0;
}
