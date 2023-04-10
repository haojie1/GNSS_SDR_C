#include "track.h"
#include "initSetting.h"
#include "config.h"

// Init TrackResult Structure
void trackResultInit(trackResult **trackResults, struct settings *receiverSetting){
    uint8_t trackResultNumber;
    trackResultNumber = receiverSetting->acqStatelliteList->size;

    uint32_t processTimeMS = receiverSetting->msToProcess;
    *trackResults = (trackResult*)malloc(sizeof(trackResult)*trackResultNumber);

    for (int i=0; i<trackResultNumber; i++) {
        (*trackResults)[i].absoluteSamples = (double*)malloc(sizeof(double)*processTimeMS);
        (*trackResults)[i].carrFreq = (double*)malloc(sizeof(double)*processTimeMS);
        (*trackResults)[i].codeFreq = (double*)malloc(sizeof(double)*processTimeMS);

        (*trackResults)[i].dllDiscr = (double*)malloc(sizeof(double)*processTimeMS);
        (*trackResults)[i].dllDiscrFilt = (double*)malloc(sizeof(double)*processTimeMS);
        (*trackResults)[i].pllDiscr = (double*)malloc(sizeof(double)*processTimeMS);
        (*trackResults)[i].pllDiscrFilt = (double*)malloc(sizeof(double)*processTimeMS);

        (*trackResults)[i].I_E = (double*)malloc(sizeof(double)*processTimeMS);
        (*trackResults)[i].I_L = (double*)malloc(sizeof(double)*processTimeMS);
        (*trackResults)[i].I_P = (double*)malloc(sizeof(double)*processTimeMS);
        (*trackResults)[i].Q_E = (double*)malloc(sizeof(double)*processTimeMS);
        (*trackResults)[i].Q_L = (double*)malloc(sizeof(double)*processTimeMS);
        (*trackResults)[i].Q_P = (double*)malloc(sizeof(double)*processTimeMS);
        (*trackResults)[i].flag = 0;
    }
}


// Calculate the coefficient of 2nd system
void calcLoopCoef(double* tau1code, double* tau2code, double bandWidth, double dampingRatio, double loopGain) {
    double Wn = bandWidth*8*dampingRatio / (4*dampingRatio*dampingRatio+1);
    *tau1code = loopGain/(Wn*Wn);
    *tau2code = 2.0*dampingRatio/Wn;
}

extern uint8_t bytePerSample[3];

// Tracking Process
void trackProcess(trackResult **trackResults, struct settings *receiverSetting, acquisitionResult* acqResults, FILE* fid) {
    //how many codePeriods needed to be process
    uint32_t codePeriods = receiverSetting->msToProcess;    

    double samplingFreq = receiverSetting->samplingFreq;
    double earlyLateSpc = receiverSetting->dllCorrelatorSpacing;

    //Summation interval
    double PDIcode = 0.001;
    double PDIcarr = 0.001;

    //The coefficients for DLL and PLL
    double tau1code, tau2code;
    double tau1carr, tau2carr;
    calcLoopCoef(&tau1code, &tau2code, receiverSetting->dllNoiseBandWidth, receiverSetting->dllDampingRatio, 1);
    calcLoopCoef(&tau1carr, &tau2carr, receiverSetting->pllNoiseBandWidth, receiverSetting->pllDampingRatio, 0.25);

    //Allocate space for PRN code
    int8_t *caCodeStandard = (int8_t*)malloc(sizeof(int8_t)*receiverSetting->codeLength); 

    for (uint8_t i=0; i<receiverSetting->acqStatelliteList->size; i++) {

        if(acqResults[i].flag == 1) {	// This satellite has been acquised
            double codeFreq = receiverSetting->codeFreqBasis;
            double remCodePhase = 0;
            double carrFreqBasis = acqResults[i].freq;
            double carrFreq = carrFreqBasis;
            double remCarrPhase = 0;

            double oldCodeNco = 0;
            double oldCodeError = 0;
            double oldCarrNco = 0;
            double oldCarrError = 0;

            printf("Track %d\n", i+1);

            (*trackResults)[i].flag = 1;

            //Make the code phase alligned
            fseek(fid, acqResults[i].codePhase*sizeof(DATATYPE), SEEK_SET);

            for (uint32_t j=0; j<codePeriods; j++) {
                double codePhaseStep = codeFreq/samplingFreq;
                //Make sure that the last sample is the first point in the next bit
                uint32_t blksize = ceil((receiverSetting->codeLength-remCodePhase)/codePhaseStep);
                //Here need deal more carefully ***************************************************************************
                DATATYPE *rawSignal = (DATATYPE*)malloc(sizeof(DATATYPE)*blksize);
                //Read data from signal file
                size_t sampleRead = fread(rawSignal, sizeof(DATATYPE), blksize, fid);
                if(sampleRead != blksize) {
                    printf("Not able to read the specified number of samples for tracking, exiting!\n");
                    fclose(fid);
                    exit(-1);
                }

                //generate E P L code
                CACode(receiverSetting->acqStatelliteList->data[i], caCodeStandard);
                int8_t * earlyCode  = (int8_t*)malloc(sizeof(int8_t)*blksize);
                int8_t * promptCode = (int8_t*)malloc(sizeof(int8_t)*blksize);
                int8_t * lateCode   = (int8_t*)malloc(sizeof(int8_t)*blksize);

                for(uint32_t i=0; i<blksize; i++) {
                    int32_t index = codePhaseStep*i+remCodePhase;
                    // Prompt code
                    promptCode[i] = caCodeStandard[index%receiverSetting->codeLength];
                    // Early code
                    index = floor(codePhaseStep*i+remCodePhase-earlyLateSpc);
                    if (index<0) 
                        index = 1022;
                    earlyCode[i] = caCodeStandard[index%receiverSetting->codeLength];
                    // Late Code
                    index = floor(codePhaseStep*i+remCodePhase+earlyLateSpc);
                    lateCode[i] = caCodeStandard[index%receiverSetting->codeLength];
                }
                remCodePhase = codePhaseStep*(blksize)+remCodePhase - receiverSetting->codeLength;

                // Generate the carrier frequency to mix the signal to baseband
                double * carrCos, * carrSin;
                carrCos = (double*)malloc(sizeof(double)*blksize);
                carrSin = (double*)malloc(sizeof(double)*blksize);

                for(uint32_t i=0; i<blksize; i++) {
                    carrCos[i] = cos(2*M_PI*carrFreq/samplingFreq*i + remCarrPhase);
                    carrSin[i] = sin(2*M_PI*carrFreq/samplingFreq*i + remCarrPhase);
                }
                remCarrPhase = remainder(remCarrPhase+(blksize)*2*M_PI*carrFreq/samplingFreq, 2*M_PI);

                double* qBasebandSignal, *iBasebandSignal;
                qBasebandSignal = (double*)malloc(sizeof(double)*blksize);
                iBasebandSignal = (double*)malloc(sizeof(double)*blksize);

                // Not be necessary
                for (uint32_t i=0; i<blksize; i++) {
                    iBasebandSignal[i] = carrSin[i]*rawSignal[i];
                    qBasebandSignal[i] = carrCos[i]*rawSignal[i];
                }

                double I_E, Q_E, I_P, Q_P, I_L, Q_L;
                I_E = Q_E = I_P = Q_P = I_L = Q_L = 0;
                for (uint32_t i=0; i<blksize; i++) {
                    I_E = I_E + iBasebandSignal[i]*earlyCode[i];
                    Q_E = Q_E + qBasebandSignal[i]*earlyCode[i];

                    I_P = I_P + iBasebandSignal[i]*promptCode[i];
                    Q_P = Q_P + qBasebandSignal[i]*promptCode[i];

                    I_L = I_L + iBasebandSignal[i]*lateCode[i];
                    Q_L = Q_L + qBasebandSignal[i]*lateCode[i];
                }

                // Find PLL error and update carrier NCO
                double carrError = atan(Q_P/I_P) / (2*M_PI);

                double carrNco = oldCarrNco + (tau2carr/tau1carr) * (carrError-oldCarrError) + carrError *(PDIcarr/tau1carr);
                oldCarrNco = carrNco;
                oldCarrError = carrError;
                
                //upated carrier freq
                carrFreq = carrFreqBasis + carrNco;
                (*trackResults)[i].carrFreq[j] = carrFreq;


                //Find DLL error and update code NCO
                double codeError = (sqrt(I_E*I_E+Q_E*Q_E)-sqrt(I_L*I_L+Q_L*Q_L)) / (sqrt(I_E*I_E+Q_E*Q_E)+sqrt(I_L*I_L+Q_L*Q_L));
                double codeNco = oldCodeNco +(tau2code/tau1code)*(codeError-oldCodeError) + codeError *(PDIcode/tau1code);
                oldCodeNco = codeNco;
                oldCodeError = codeError;


                //updated code freq
                codeFreq = receiverSetting->codeFreqBasis - codeNco;
                (*trackResults)[i].codeFreq[j] = codeFreq;


                //store information
                (*trackResults)[i].absoluteSamples[j] = ftell(fid);
                (*trackResults)[i].dllDiscr[j] = codeError;
                (*trackResults)[i].dllDiscrFilt[j] = codeNco;
                (*trackResults)[i].pllDiscr[j] = carrError;
                (*trackResults)[i].pllDiscrFilt[j] = carrNco;
                (*trackResults)[i].I_E[j] = I_E;
                (*trackResults)[i].Q_E[j] = Q_E;
                (*trackResults)[i].I_P[j] = I_P;
                (*trackResults)[i].Q_P[j] = Q_P;
                (*trackResults)[i].I_L[j] = I_L;
                (*trackResults)[i].Q_L[j] = Q_L;


                // freee space
                free(rawSignal);
                free(earlyCode);
                free(promptCode);
                free(lateCode);
                free(carrCos);
                free(carrSin);
                free(iBasebandSignal);
                free(qBasebandSignal);
            }
        }
    }
    free(caCodeStandard);
}
