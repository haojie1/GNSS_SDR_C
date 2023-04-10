#ifndef __TRACK_H__
#define __TRACK_H__

#include <stdint.h>
#include <math.h>
#include "config.h"
#include "acquisition.h"
#include "generateCACode.h"

struct trackResultStruct{
    //I path's output
    double * absoluteSamples;
    double * codeFreq;
    double * carrFreq;

    //Output From E P L Path
    double * I_E;
    double * I_P;
    double * I_L;
    double * Q_E;
    double * Q_P;
    double * Q_L;

    // Loop discriminators
    double * dllDiscr;
    double * dllDiscrFilt;
    double * pllDiscr;
    double * pllDiscrFilt;

    // flag determing whether this satellite is tracked
    uint8_t flag;
};

typedef struct trackResultStruct trackResult;

void trackResultInit(trackResult **trackResults, struct settings *receiverSetting) ;
void trackProcess(trackResult **trackResults, struct settings *receiverSetting, acquisitionResult* acqResults, FILE* fid) ;



#endif  //__TRACK_H__
