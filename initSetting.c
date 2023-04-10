#include <string.h>
#include <gsl/gsl_math.h>
#include "config.h"


void initSetting(struct settings *initSetting)
{
    //strcpy(initSetting->fileName, "../GNSS_signal_records/GPSdata-DiscreteComponents-fs38_192-if9_55.bin");     //The file containing signal
    strcpy(initSetting->fileName, "L1Ant1.dat");     //The file containing signal

    initSetting->msToProcess = 37000;                	//Total Time 37000ms
    initSetting->numberOfChannel = 8;               	//Total channel numbers
    initSetting->intermediatFreq = 1.42E6;          	//IF
    initSetting->samplingFreq    = 10E6;        	//Sampling Frequency
    initSetting->acqStatelliteList = gsl_vector_int_alloc(32);  //The list of satellite needed to be dealed with
    for (int i=0; i<32; i++) {
        gsl_vector_int_set(initSetting->acqStatelliteList, i, i+1);
    }
    initSetting->acqSearchBand = 14;                	//-7KHz - 7KHz
    initSetting->acqThreshold = 2;
    initSetting->codeFreqBasis = 1.023E6;
    initSetting->codeLength = 1023;
    initSetting->startOffset = 68.802;

    initSetting->dllCorrelatorSpacing = 0.5;        	//The space between early code and late code
    initSetting->dllDampingRatio = 0.707;
    initSetting->dllNoiseBandWidth = 2;

    initSetting->pllDampingRatio = 0.707;
    initSetting->pllNoiseBandWidth = 25;

    initSetting->elevationMask = 10;                	//Exclude the satellites at low elevation
    initSetting->navSolPeriod = 500;                	//Output the estimated position once 500ms
    initSetting->skipNumberOfByte = 0;

    initSetting->skipNumberOfByte = 0;

    initSetting->truePosition.E = GSL_NAN;
    initSetting->truePosition.N = GSL_NAN;
    initSetting->truePosition.U = GSL_NAN;          	//The true position of receiver

    initSetting->useTropCorr = 1;			//Use troposphere correction
}
