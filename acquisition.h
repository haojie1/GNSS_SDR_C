#ifndef __ACQUSITION_H__
#define __ACQUSITION_H__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_vector.h>
#include "tool.h"
#include "generateCACode.h"
#include "config.h"

// This structure used to store information about acquisition
struct acquisitionResultStruct {
    double freq;
    uint32_t codePhase;
    uint32_t flag;
};

typedef struct acquisitionResultStruct acquisitionResult;

acquisitionResult* acquisitionProcess(FILE* fid, const struct settings* receiverSetting) ;

#endif
