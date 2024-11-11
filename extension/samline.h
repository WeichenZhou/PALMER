#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#ifndef PALMER_SamLine_H
#define PALMER_SamLine_H

using namespace std;

class SamLine
{
public:
    string QNAME; // SRR3658462.112531
    string LEN;   // length=23647
    string FLAG;  // 0
    string RNAME; // chr19
    string POS;   // 4977409
    string MAPQ;  // 60
    string CIGAR; // 8S10=1D5=1X18=1I7=1D33=...
    string RNEXT; // *
    string PNEXT; // 0
    string TLEN;  // 0
    string SEQ;   // ATGC
    string QUAL;  // !@#$%^&*()

    vector<string> UNKNOWNs;

    string QNAME_LEN; // SRR3658462.112531_length=23647
    int INT_FLAG;
    int INT_POS;
    int INT_MAPQ;
    int INT_PNEXT;
    int INT_TLEN;

public:
    SamLine(char *line)
    {
        istringstream issLine(line);
        issLine >> QNAME; //QNAME
        issLine >> LEN;   //LEN
        issLine >> FLAG;  //FLAG
        issLine >> RNAME; //RNAME
        issLine >> POS;   //POS
        issLine >> MAPQ;  //MAPQ
        issLine >> CIGAR; //CIGAR
        issLine >> RNEXT; //RNEXT
        issLine >> PNEXT; //PNEXT
        issLine >> TLEN;  //TLEN
        issLine >> SEQ;   //SEQ
        issLine >> QUAL;  //QUAL

        QNAME_LEN = QNAME + "_" + LEN;
        INT_FLAG = atoi(FLAG.c_str());
        INT_POS = atoi(POS.c_str());
        INT_MAPQ = atoi(MAPQ.c_str());
        INT_PNEXT = atoi(PNEXT.c_str());
        INT_TLEN = atoi(TLEN.c_str());

        while (!issLine.eof())
        {
            string tmp;
            issLine >> tmp;
            UNKNOWNs.push_back(tmp);
        }
        // printf("[SamLine] UNKNOWNs.size = %d\n", UNKNOWNs.size());
    };
    ~SamLine()
    {
        UNKNOWNs.clear();
    };
};

#endif