/**
 * @file main.c
 * @author haojie
 * @brief 
 *  Create A GPS L1 Receiver Base on GNSS_SDR
 * @version 0.1
 * @date 2023-03-05
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "tool.h"
#include "initSetting.h"
#include "generateCACode.h"
#include "acquisition.h"
#include "track.h"
#include "findPreamble.h"
#include "postProceeding.h"


#define     BUFSIZE     1024

struct settings receiverSetting;

acquisitionResult *acqResults;
trackResult *trackResults;

int8_t  buff[BUFSIZE];

int main(int argc, char *argv[]) 
{
    printf("Welcome to: GNSS_SDR_C\n\n");
    printf("Load Configuration...\n");

    // Init the receiver setting
    initSetting(&receiverSetting);

    // Read signal from file
    FILE *fid = fopen(receiverSetting.fileName, "r");

    // Acqusition 
    acqResults = acquisitionProcess(fid, &receiverSetting);

    printf("\n\nACQUISITION RESULT: \n");
    for (int i=0; i<32; i++) {
        if(acqResults[i].flag == 1) {
            printf("Satellite %d: %lf, %d\n", i+1, acqResults[i].freq, acqResults[i].codePhase+1);
        }
    }

    // Tracking
    time_t t;
    time(&t);
    char* timeString = ctime(&t);

    printf("tracking start: %s\n", timeString);

    trackResultInit(&trackResults, &receiverSetting);
    trackProcess(&trackResults, &receiverSetting, acqResults, fid);
    time(&t);
    timeString = ctime(&t);

    printf("tracking end: %s\n", timeString);

    // Stroe Tracking Results
    FILE *fid_IE_21 = fopen("IE.dat", "w"); 
    FILE *fid_QE_21 = fopen("QE.dat", "w"); 
    FILE *fid_IP_21 = fopen("IP.dat", "w"); 
    FILE *fid_QP_21 = fopen("QP.dat", "w"); 
    FILE *fid_IL_21 = fopen("IL.dat", "w"); 
    FILE *fid_QL_21 = fopen("QL.dat", "w"); 
    FILE *fid_absolutSamples = fopen("absoluteSamples.dat", "w");

    for (uint32_t i=0; i<receiverSetting.acqStatelliteList->size; i++) {
        if (acqResults[i].flag) {
            fwrite(trackResults[i].I_E, sizeof(double), receiverSetting.msToProcess, fid_IE_21);
            fwrite(trackResults[i].Q_E, sizeof(double), receiverSetting.msToProcess, fid_QE_21);
            fwrite(trackResults[i].I_P, sizeof(double), receiverSetting.msToProcess, fid_IP_21);
            fwrite(trackResults[i].Q_P, sizeof(double), receiverSetting.msToProcess, fid_QP_21);
            fwrite(trackResults[i].I_L, sizeof(double), receiverSetting.msToProcess, fid_IL_21);
            fwrite(trackResults[i].Q_L, sizeof(double), receiverSetting.msToProcess, fid_QL_21);
            fwrite(trackResults[i].absoluteSamples, sizeof(double), receiverSetting.msToProcess, fid_absolutSamples);
        }
    }

    fclose(fid_IE_21);
    fclose(fid_QE_21);
    fclose(fid_IP_21);
    fclose(fid_QP_21);
    fclose(fid_IL_21);
    fclose(fid_QL_21);
    fclose(fid_absolutSamples);

    fclose(fid);

    free(acqResults);
    free(trackResults);

    /*

    int32_t * pos = tool_find_preambles(trackResults, &receiverSetting);

    // the preamble position is index from 0, not 1
    for (uint32_t i = 0; i<receiverSetting.acqStatelliteList->size; i++) {
        printf("preamble pos: %d\n", pos[i]);
    }

    */
    
    // Get Receiver Position
   struct posLLH pos = getPos(receiverSetting); 

   printf("latitude: %f\n", pos.latitude);
   printf("longitude: %f\n", pos.longitude);
   printf("height: %f\n", pos.heigth);

    return 0;
}
