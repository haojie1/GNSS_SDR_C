#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <gsl/gsl_vector_int.h>

#define MAX_FILENAME_LENGTH 128

// Unit: m/s
#define SPEED_OF_LIGHT  (299792458)

// Unit: ms
#define START_OFF_START (68.802)

#define gpsPI   (3.1415926535898)
#define DATATYPE int8_t

struct ENUStructure {
    double E, N, U;
};

struct ECEFStructure {
    double X, Y, Z;
};

struct settings {
    int msToProcess;
    int numberOfChannel;
    int skipNumberOfByte;
    char fileName[MAX_FILENAME_LENGTH];
    double intermediatFreq;
    double samplingFreq;
    double codeFreqBasis;
    int codeLength; 
    gsl_vector_int *acqStatelliteList; 
    int acqSearchBand;
    double acqThreshold;
    double startOffset;

    // For DLL
    double dllDampingRatio;
    double dllNoiseBandWidth;
    double dllCorrelatorSpacing;

    // For PLL
    double pllDampingRatio;
    double pllNoiseBandWidth;

    // For Navigation Solution
    int navSolPeriod;
    int elevationMask;
    int useTropCorr;

    // For perior true position
    struct ENUStructure truePosition;
};



#endif
