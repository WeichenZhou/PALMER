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

    //cout<<"Samtools Step is now running."<<endl;
    //##########hard code warming########## -q -F ##########
    //std::ios::sync_with_stdio(false);
    //std::cin.tie(0);

    // string sys;
    // sys = "samtools view -q 10 -F 0x100 -F 0x200 -F 0x800 -F 0x400 " + input_bam + " " + chr + ":" + start + "-" + end + " |sed -e 's/[ ][ ]*/_/g'  > " + working_dir + "region.sam";

    // char *syst = new char[sys.length() + 1];
    // strcpy(syst, sys.c_str());
    // system(syst);

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

    free(argv[ARG_SIZE - 2]);
    free(argv[ARG_SIZE - 1]);

    //string sys_replace;
    //sys_replace="sed -e 's/[ ][ ]*/_/g' "+working_dir+"region.pre.sam"+" > "+working_dir+"region.sam";

    ///char *syst_replace = new char[sys_replace.length()+1];
    //strcpy(syst_replace, sys_replace.c_str());
    //system(syst_replace);

    return 0;
}
