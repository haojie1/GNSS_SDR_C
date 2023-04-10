#include "acquisition.h"
#include "config.h"

// The Structure of acquisition result
struct acquisitionPoint {
    uint32_t freqBin;
    uint32_t codePhase;
    uint32_t value;
};

struct finedFreq {
    uint32_t findFreqBin;
    double value; 
};

// Acquisition Process
acquisitionResult* acquisitionProcess(FILE* fid, const struct settings* receiverSetting) {

    uint32_t sampleNumbersPerCode = receiverSetting->samplingFreq / 1E3;
    uint32_t sampleNumbersPerChip = round(receiverSetting->samplingFreq / 1E3 / receiverSetting->codeLength);
    // Allocate memory for every satellite acquisition
    acquisitionResult * acqResult = (acquisitionResult*)malloc(sizeof(acquisitionResult)*receiverSetting->acqStatelliteList->size);

    // Flag the acquisition: 0--negetive; 1--positive
    for (uint32_t i=0; i<receiverSetting->acqStatelliteList->size; i++) {
        acqResult[i].flag = 0;
    }

        // 1. Read 2ms Data From the signal File
        DATATYPE*v1 = (DATATYPE*)malloc(sizeof(DATATYPE)*sampleNumbersPerCode);
        DATATYPE*v2 = (DATATYPE*)malloc(sizeof(DATATYPE)*sampleNumbersPerCode);
        // 11ms data to estimate the fined frequency
        DATATYPE*fineFreqDat = (DATATYPE*)malloc(sizeof(DATATYPE)*sampleNumbersPerCode*11);

        fread(v1, sizeof(DATATYPE), sampleNumbersPerCode, fid);
        fread(v2, sizeof(DATATYPE), sampleNumbersPerCode, fid);
        fseek(fid, 0, SEEK_SET);
        fread(fineFreqDat, sizeof(DATATYPE), sampleNumbersPerCode*11, fid);

        // 2. Get prn code
        double* codeData = (double*)malloc(sizeof(double)*sampleNumbersPerCode*2);
        int8_t*  caCodeSampling = (int8_t*)malloc(sizeof(int8_t)*sampleNumbersPerCode);
        for (uint32_t i=0; i<(receiverSetting->acqStatelliteList)->size; i++) {
            caCodeAfterSampling(receiverSetting, sampleNumbersPerCode, caCodeSampling, receiverSetting->acqStatelliteList->data[i]);

            // pack int8_t to gsl_vector_complex
            for (uint32_t i=0; i<sampleNumbersPerCode; i++) {
		    REAL(codeData,i) = *(caCodeSampling+i);
		    IMAG(codeData,i) = 0;
            }

            // calculate the fft of C/A code
            gsl_fft_complex_wavetable* waveTable;
            gsl_fft_complex_workspace* workSpace;

            waveTable = gsl_fft_complex_wavetable_alloc(sampleNumbersPerCode);
            workSpace = gsl_fft_complex_workspace_alloc(sampleNumbersPerCode);

            gsl_fft_complex_forward(codeData, 1, sampleNumbersPerCode, waveTable, workSpace);

            // conjugat the code fft result
            for (uint32_t i=0; i<sampleNumbersPerCode; i++) {
                IMAG(codeData,i) = -IMAG(codeData,i);
            }

            // check every possible frequency bin, step is 0.5Khz
            double * sinData = (double*)malloc(sizeof(double)*sampleNumbersPerCode);
            double * cosData = (double*)malloc(sizeof(double)*sampleNumbersPerCode);
            double * i_data_base_1 = (double*)malloc(sizeof(double)*sampleNumbersPerCode);
            double * q_data_base_1 = (double*)malloc(sizeof(double)*sampleNumbersPerCode);
            double * complex_data_base_1 = (double*)malloc(sizeof(double)*2*sampleNumbersPerCode);
            double * i_data_base_2 = (double*)malloc(sizeof(double)*sampleNumbersPerCode);
            double * q_data_base_2 = (double*)malloc(sizeof(double)*sampleNumbersPerCode);
            double * complex_data_base_2 = (double*)malloc(sizeof(double)*2*sampleNumbersPerCode);
            double * multiplexed_fft_result_1 = (double*)malloc(sizeof(double)*2*sampleNumbersPerCode);
            double * multiplexed_fft_result_2 = (double*)malloc(sizeof(double)*2*sampleNumbersPerCode);

            double phasePerSamplePoint = 2*M_PI/receiverSetting->samplingFreq;
            double **twoDimResult = (double**)malloc(sizeof(double*)*receiverSetting->acqSearchBand*2);

            // store the acquisition result
            for (int i=0; i<receiverSetting->acqSearchBand*2; i++) {
                twoDimResult[i] = (double*)malloc(sizeof(double)*sampleNumbersPerCode);
            }

            for (int i=0; i<receiverSetting->acqSearchBand*2; i++) {
                // frequency for now
                double freqBin = receiverSetting->intermediatFreq - receiverSetting->acqSearchBand/2*1E3 + i*0.5*1E3;
                // local oscillator
                for (uint32_t i=0; i<sampleNumbersPerCode; i++) {
                    sinData[i] = sin(i*freqBin*phasePerSamplePoint);
                    cosData[i] = cos(i*freqBin*phasePerSamplePoint);
                }
                // downconversion to baseband
                for (uint32_t i=0; i<sampleNumbersPerCode; i++) {
                    i_data_base_1[i] = sinData[i] * v1[i];
                    q_data_base_1[i] = cosData[i] * v1[i];
                    REAL(complex_data_base_1, i) = i_data_base_1[i];
                    IMAG(complex_data_base_1, i) = q_data_base_1[i];

                    i_data_base_2[i] = sinData[i] * v2[i];
                    q_data_base_2[i] = cosData[i] * v2[i];
                    REAL(complex_data_base_2, i) = i_data_base_2[i];
                    IMAG(complex_data_base_2, i) = q_data_base_2[i];
                } 

                
                // Calculate the fft of the baseband signal
                gsl_fft_complex_forward(complex_data_base_1, 1, sampleNumbersPerCode, waveTable, workSpace);
                gsl_fft_complex_forward(complex_data_base_2, 1, sampleNumbersPerCode, waveTable, workSpace);

                // execute the multiplication of fft result
                for (uint32_t i=0; i<sampleNumbersPerCode; i++) {
                    REAL(multiplexed_fft_result_1, i) = REAL(codeData, i)*REAL(complex_data_base_1, i) - IMAG(codeData, i)*IMAG(complex_data_base_1, i);
                    IMAG(multiplexed_fft_result_1, i) = REAL(codeData, i)*IMAG(complex_data_base_1, i) + IMAG(codeData, i)*REAL(complex_data_base_1, i);
                    REAL(multiplexed_fft_result_2, i) = REAL(codeData, i)*REAL(complex_data_base_2, i) - IMAG(codeData, i)*IMAG(complex_data_base_2, i);
                    IMAG(multiplexed_fft_result_2, i) = REAL(codeData, i)*IMAG(complex_data_base_2, i) + IMAG(codeData, i)*REAL(complex_data_base_2, i);
                }


                // invert the fft result to get the corelation 
                gsl_fft_complex_backward(multiplexed_fft_result_1, 1, sampleNumbersPerCode, waveTable, workSpace);
                gsl_fft_complex_backward(multiplexed_fft_result_2, 1, sampleNumbersPerCode, waveTable, workSpace);

                // gsl's invert fft is portional to ture fft result
                for (uint32_t i=0; i<sampleNumbersPerCode; i++) {
                    REAL(multiplexed_fft_result_1, i) = REAL(multiplexed_fft_result_1, i)/sampleNumbersPerCode;
                    IMAG(multiplexed_fft_result_1, i) = IMAG(multiplexed_fft_result_1, i)/sampleNumbersPerCode;
                    REAL(multiplexed_fft_result_2, i) = REAL(multiplexed_fft_result_2, i)/sampleNumbersPerCode;
                    IMAG(multiplexed_fft_result_2, i) = IMAG(multiplexed_fft_result_2, i)/sampleNumbersPerCode;
                }

                // find the max power bwtween multiplexed_fft_reslt_1 and multiplexed_fft_result2
                double *abs_correlation_value_1, *abs_correlation_value_2;
                double abs_correlation_value_max_1, abs_correlation_value_max_2;
                abs_correlation_value_max_1 = abs_correlation_value_max_2 = 0;
                abs_correlation_value_1 = (double*)malloc(sizeof(double)*sampleNumbersPerCode);
                abs_correlation_value_2 = (double*)malloc(sizeof(double)*sampleNumbersPerCode);

                for (uint32_t i=0; i<sampleNumbersPerCode; i++) {
                    abs_correlation_value_1[i] = REAL(multiplexed_fft_result_1, i)*REAL(multiplexed_fft_result_1, i)+IMAG(multiplexed_fft_result_1, i)*IMAG(multiplexed_fft_result_1, i);
                    abs_correlation_value_max_1 = abs_correlation_value_max_1 >= abs_correlation_value_1[i] ? abs_correlation_value_max_1 : abs_correlation_value_1[i];
                    abs_correlation_value_2[i] = REAL(multiplexed_fft_result_2, i)*REAL(multiplexed_fft_result_2, i)+IMAG(multiplexed_fft_result_2, i)*IMAG(multiplexed_fft_result_2, i);
                    abs_correlation_value_max_2 = abs_correlation_value_max_2 >= abs_correlation_value_2[i] ? abs_correlation_value_max_2 : abs_correlation_value_2[i];
                }


                if (abs_correlation_value_max_1>abs_correlation_value_max_2) {
                    memcpy(twoDimResult[i], abs_correlation_value_1, sizeof(double)*sampleNumbersPerCode);
                } else {
                    memcpy(twoDimResult[i], abs_correlation_value_2, sizeof(double)*sampleNumbersPerCode);
                }
                free(abs_correlation_value_1);
                free(abs_correlation_value_2);

            }

            // Free space
            gsl_fft_complex_wavetable_free(waveTable);
            gsl_fft_complex_workspace_free(workSpace);

            // find the max point in the 2D plan
            double max = 0;
            struct acquisitionPoint acq;
            for (int i=0; i<receiverSetting->acqSearchBand*2; i++) {
                for (uint32_t j=0; j<sampleNumbersPerCode; j++) {
                    double temp = twoDimResult[i][j]; 
                    if (max<temp) {
                        max = temp;
                        acq.freqBin = i;
                        acq.codePhase = j;            
                        acq.value = max;
                    }
                }
            }

            // find the second max point in the 2D plan
            uint32_t freqBin;
            freqBin = acq.freqBin;

            double secondMax = 0;
            uint32_t startPoint, endPoint;
            startPoint = endPoint = 0;
            if (acq.codePhase <= sampleNumbersPerChip) {
                startPoint = acq.codePhase + sampleNumbersPerChip;
                endPoint = sampleNumbersPerCode - (sampleNumbersPerChip - acq.codePhase);
                for (uint32_t i=startPoint; i<endPoint; i++) {
                    if (secondMax<twoDimResult[freqBin][i]) {
                        secondMax = twoDimResult[freqBin][i];
                    }
                }
            } else if (acq.codePhase+sampleNumbersPerChip >= sampleNumbersPerCode) {
                startPoint = acq.codePhase+sampleNumbersPerChip - sampleNumbersPerCode;
                endPoint = acq.codePhase-sampleNumbersPerChip;
                for (uint32_t i=startPoint; i<endPoint; i++) {
                    if (secondMax<twoDimResult[freqBin][i]) {
                        secondMax = twoDimResult[freqBin][i];
                    }
                }
            } else {
                for (uint32_t i=0; i<acq.codePhase-sampleNumbersPerChip; i++) {
                    if (secondMax<twoDimResult[freqBin][i]) {
                        secondMax = twoDimResult[freqBin][i];
                    }
                }
                for (uint32_t i=acq.codePhase+sampleNumbersPerChip; i<sampleNumbersPerCode; i++) {
                    if (secondMax<twoDimResult[freqBin][i]) {
                        secondMax = twoDimResult[freqBin][i];
                    }
                }
            }

            if (max/secondMax >= receiverSetting->acqThreshold) {
                printf("satellite %d is Acqitioned\n", i);
                // fine the frequency
                double meanValue, sum;
                sum = meanValue = 0;
                for (uint32_t i=0; i<sampleNumbersPerCode*11; i++) {
                    sum += fineFreqDat[i];
                }
                meanValue = sum/sampleNumbersPerCode/11;

                uint32_t totalNumberForFFT = next2pow(sampleNumbersPerCode*10)*8;
                double * paddingDataFFT = (double*)malloc(sizeof(double)*totalNumberForFFT*2);
                double * absFFTValue = (double*)malloc(sizeof(double)*totalNumberForFFT);
                int8_t * caCodeSamplingFor10MS = (int8_t*)malloc(sizeof(int8_t)*sampleNumbersPerCode*10);
                // generate long CACode
                caCodeAfterSampling(receiverSetting, sampleNumbersPerCode*10, caCodeSamplingFor10MS, receiverSetting->acqStatelliteList->data[i]);

                memset(paddingDataFFT, 0, sizeof(double)*totalNumberForFFT*2);
                for(uint32_t i=0; i<sampleNumbersPerCode*10; i++) {
                    REAL(paddingDataFFT, i) = (fineFreqDat[i+acq.codePhase]-meanValue)*caCodeSamplingFor10MS[i];
                }


                //FILE *fid_temp = fopen("paddingDataFFT.dat", "w");
                //fwrite(paddingDataFFT, sizeof(double), totalNumberForFFT*2, fid_temp);
                //fclose(fid_temp);

                gsl_fft_complex_wavetable* fineWaveTable = gsl_fft_complex_wavetable_alloc(totalNumberForFFT);
                gsl_fft_complex_workspace* fineWorkSpace = gsl_fft_complex_workspace_alloc(totalNumberForFFT);

                gsl_fft_complex_forward(paddingDataFFT, 1, totalNumberForFFT, fineWaveTable, fineWorkSpace);

                for (uint32_t i=0; i<totalNumberForFFT; i++) {
                    absFFTValue[i] = sqrt(REAL(paddingDataFFT, i)*REAL(paddingDataFFT, i) + IMAG(paddingDataFFT, i)*IMAG(paddingDataFFT, i));
                }
                free(paddingDataFFT);
                gsl_fft_complex_wavetable_free(fineWaveTable);
                gsl_fft_complex_workspace_free(fineWorkSpace);

                uint32_t uniqueFreq = ceil((totalNumberForFFT+1.0)/2);

                struct finedFreq finedFreqStructure;
                finedFreqStructure.findFreqBin = 0;
                finedFreqStructure.value  = 0;
                for (uint32_t i=5; i<uniqueFreq-5; i++) {
                    if (finedFreqStructure.value<absFFTValue[i]) {
                        finedFreqStructure.value = absFFTValue[i];
                        finedFreqStructure.findFreqBin = i;
                    }
                }
                acqResult[i].flag = 1;
                acqResult[i].codePhase = acq.codePhase;
                acqResult[i].freq = receiverSetting->samplingFreq/totalNumberForFFT*(finedFreqStructure.findFreqBin+1);

                free(absFFTValue);
                free(caCodeSamplingFor10MS);
		for (int i=0; i<receiverSetting->acqSearchBand*2; i++) {
			free(twoDimResult[i]);
		}
                free(twoDimResult);
            } else {
                printf("satellite %d isnot Acqitioned\n", i);
            }
	    free(complex_data_base_1);
	    free(complex_data_base_2);
            free(i_data_base_1);
            free(q_data_base_1);
            free(i_data_base_2);
            free(q_data_base_2);
            free(sinData);
            free(cosData);
            free(multiplexed_fft_result_1);
            free(multiplexed_fft_result_2);

          }  
	    free(v1);
	    free(v2);
	    free(codeData);
	    free(caCodeSampling);
	    free(fineFreqDat);

    return acqResult;

}
