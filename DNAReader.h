#ifndef DNAREADER_H
#define DNAREADER_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;

class DNAReader {

  public:
    DNAReader();
    DNAReader(string fileName);

    int ComputeTotals();
    float ComputeMean();
    float ComputeVariance();
    float ComputeStddev();

    float ComputeRelA();
    float ComputeRelC();
    float ComputeRelT();
    float ComputeRelG();
    float ComputeRelAA();
    float ComputeRelAC();
    float ComputeRelAT();
    float ComputeRelAG();
    float ComputeRelCA();
    float ComputeRelCC();
    float ComputeRelCT();
    float ComputeRelCG();
    float ComputeRelTA();
    float ComputeRelTC();
    float ComputeRelTT();
    float ComputeRelTG();
    float ComputeRelGA();
    float ComputeRelGC();
    float ComputeRelGT();
    float ComputeRelGG();
    int ComputeRandomGaussian();
    string ProduceGaussianLine();

    void WriteAll();

  private:
    string ifilename;
    string ofilename;
    string nucleodump;

    int totalLines;
    int totalA;
    int totalC;
    int totalT;
    int totalG;
    int totalAA;
    int totalAC;
    int totalAT;
    int totalAG;
    int totalCA;
    int totalCC;
    int totalCT;
    int totalCG;
    int totalTA;
    int totalTC;
    int totalTT;
    int totalTG;
    int totalGA;
    int totalGC;
    int totalGT;
    int totalGG;

    float relA;
    float relC;
    float relT;
    float relG;
    float relAA;
    float relAC;
    float relAT;
    float relAG;
    float relCA;
    float relCC;
    float relCT;
    float relCG;
    float relTA;
    float relTC;
    float relTT;
    float relTG;
    float relGA;
    float relGC;
    float relGT;
    float relGG;

    int sum;
    int totalbigrams;
    float mean;
    float variance;
    float stddev;

};

#endif
