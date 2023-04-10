#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector_int.h>
#include "track.h"
#include "config.h"
#include "initSetting.h"
#include "findPreamble.h"
#include "ephemeris.h"
#include "satpos.h"
#include "common.h"
#include "cart2geo.h"

// Find the start position of frame and correct the invert problem
int8_t** assamble_1ms_to_bit(trackResult* trackResults, int32_t* pos, struct settings *receiverSetting) {
    int8_t** bits = (int8_t**)malloc(sizeof(int8_t*)*32);
    for (uint32_t i=0; i<receiverSetting->acqStatelliteList->size; i++) {
        bits[i] = (int8_t*)malloc(sizeof(int8_t)*(1500+1));
        memset(bits[i], -2, sizeof(int8_t)*(1500+1));
    }

    for (uint32_t i=0; i<receiverSetting->acqStatelliteList->size; i++) {
        if (trackResults[i].flag == 1) {
            // a complete frame and a last bit of the previous frame 
            int8_t counter = 0;
            int32_t j=0;
            double sum = 0;
            for (int32_t k=pos[i]-20; k<(pos[i]+(1500*20)); k++) {
                if(counter++ == 20) {
                    bits[i][j++] =  sum > 0 ? 1 : 0;
                    counter = 1;
                    sum = trackResults[i].I_P[k];
                } else {
                    sum += trackResults[i].I_P[k];
                }
            }

	   // Deal with the last bit
	   bits[i][j] = sum > 0 ? 1 : 0;
        }
    }
    return bits;
}

// use two uint8_t array to combind a unsigned integer
// most bit first
// array1[len-1] -> lowest bit
// array1[len-1] array1[len-2] ... array1[0] array2[len2-1] array2[len2-2] ... array2[0]
uint32_t intArray_to_bin(uint8_t *array1, uint8_t *array2, uint8_t len1, uint8_t len2) {
    uint32_t result = 0, temp = 0;
    if (array2 == NULL) {
        for (int8_t i=len1-1; i>-1; i--){
            temp = 1<<(len1-1-i);
            if (array1[i] == 1) {
                result += temp;
            }
        }
    } else {
        for (int8_t i=len1-1; i>-1; i--) {
            temp = 1<<(len1-1-i);
            if (array1[i] == 1) {
                result += temp;
            }
        }
        for (int8_t i=len2-1; i>-1; i--) {
            temp = 1<<(len2-1-i+len1);
            if (array2[i] == 1) {
                result += temp;
            }
        }
    }

    return result;
}

// use two uint8_t array to combind a signed integer
// most bit first
// array1[len-1] -> lowest bit
// array1[len-1] array1[len-2] ... array1[0] array2[len2-1] array2[len2-2] ... array2[0]
int32_t intArray_to_complement_bin(uint8_t *array1, uint8_t *array2, uint8_t len1, uint8_t len2) {
    int32_t result = 0, temp = 0;
    if (array2 == NULL) {
	// check if this is a negtive number
	if (array1[0] == 1) {	// negtive
        	for (uint8_t i=len1-1; i>0; i--){
        	    temp = 1<<(len1-1-i);
        	    if (array1[i] == 1) {
        	        result += temp;
        	    }
		}
		result = result - pow(2, len1-1);
	} else {
        	for (uint8_t i=len1-1; i>0; i--){
        	    temp = 1<<(len1-1-i);
        	    if (array1[i] == 1) {
        	        result += temp;
        	    }
		}
        }

    } else {
	if (array2[0] == 1) { //negtive
        	for (int8_t i=len1-1; i>-1; i--){
        	    temp = 1<<(len1-1-i);
        	    if (array1[i] == 1) {
        	        result += temp;
        	    }
        	}
        	for (uint8_t i=len2-1; i>0; i--) {
        	    temp = 1<<(len2-1-i+len1);
        	    if (array2[i] == 1) {
        	        result += temp;
        	    }
        	}

		result = result - pow(2, len1+len2-1);

	} else {
        	for (int8_t i=len1-1; i>-1; i--) {
        	    temp = 1<<(len1-1-i);
        	    if (array1[i] == 1) {
        	        result += temp;
        	    }
        	}
        	for (uint8_t i=len2-1; i>0; i--) {
        	    temp = 1<<(len2-1-i+len1);
        	    if (array2[i] == 1) {
        	        result += temp;
        	    }
        	}
	}
    }

    return result;
}

void checkPhase(uint8_t* word, int8_t d30Start) {
	if (d30Start == 1) {
		for (uint8_t i=0; i<24; i++) {
			word[i] = (word[i] == 0);
		}
	}
}

void get_ephemeris(struct ephemeris* ephemerisStructs, int8_t ** bits) {
    uint8_t subframe[300];
    uint8_t subframe_id;
    FILE* fid_eph = fopen("eph.dat", "w");

    for (uint32_t i=0; i<32; i++) {
        // this bits[i] is tracked
        if (bits[i][0] != -2) {
            for (uint8_t k=0; k<5; k++) {
		// throw the first bit
                memcpy(subframe, bits[i]+k*300+1, sizeof(int8_t)*300);
		int8_t d30Start = *(bits[i]);

		for (uint8_t i=0; i<10; i++) {
			checkPhase(subframe+i*30, d30Start);
			d30Start = subframe[30*(i+1)-1];
		}

                subframe_id = intArray_to_bin(subframe+49, NULL, 3,  0);

                switch (subframe_id)
                {
                case 1: 
			ephemerisStructs[i].weekNumber	= intArray_to_bin(subframe+60, NULL, 10, 0) + 1024;
			ephemerisStructs[i].accuracy	= intArray_to_bin(subframe+72, NULL ,4, 0);
			ephemerisStructs[i].health 	= intArray_to_bin(subframe+76, NULL,5, 0);
			ephemerisStructs[i].T_GD 		= intArray_to_complement_bin(subframe+196, NULL, 8, 0) * pow(2, -31);
			ephemerisStructs[i].IODC 		= intArray_to_bin(subframe+210, subframe+82, 8, 2);
			ephemerisStructs[i].t_oc 		= intArray_to_bin(subframe+218, NULL, 16, 0) * pow(2, 4);
			ephemerisStructs[i].a_f2 		= intArray_to_complement_bin(subframe+240, NULL, 8, 0) * pow(2, -55);
			ephemerisStructs[i].a_f1 		= intArray_to_complement_bin(subframe+248, NULL, 16, 0) * pow(2, -43);
			ephemerisStructs[i].a_f0 		= intArray_to_complement_bin(subframe+270, NULL, 22, 0) * pow(2, -31);
                    break;
                case 2:
		    	ephemerisStructs[i].IODE_sf2 	= intArray_to_bin(subframe+60, NULL, 8, 0);
			ephemerisStructs[i].C_rs 		= intArray_to_complement_bin(subframe+68, NULL, 16, 0) * pow(2, -5);
			ephemerisStructs[i].deltan 	= intArray_to_complement_bin(subframe+90, NULL, 16, 0) * pow(2, -43) * gpsPI;
			ephemerisStructs[i].M_0 		= intArray_to_complement_bin(subframe+120, subframe+106, 24, 8) * pow(2, -31) * gpsPI;
			ephemerisStructs[i].C_uc 		= intArray_to_complement_bin(subframe+150, NULL, 16, 0) * pow(2, -29);
			ephemerisStructs[i].e 		= intArray_to_complement_bin(subframe+180, subframe+166, 24, 8) * pow(2, -33);
			ephemerisStructs[i].C_us 		= intArray_to_complement_bin(subframe+210, NULL, 16, 0) * pow(2, -29);
			ephemerisStructs[i].sqrtA 	= intArray_to_bin(subframe+240, subframe+226, 24, 8) * pow(2, -19);
			ephemerisStructs[i].t_oe 		= intArray_to_bin(subframe+270, NULL, 16, 0) * pow(2, 4);
                    break;
                case 3:
		    	ephemerisStructs[i].C_ic 		= intArray_to_complement_bin(subframe+60, NULL, 16, 0) * pow(2, -29);
		    	ephemerisStructs[i].omega_0 		= intArray_to_complement_bin(subframe+90, subframe+76, 24, 8) * pow(2, -31) * gpsPI;
		    	ephemerisStructs[i].C_is 		= intArray_to_complement_bin(subframe+120, NULL, 16, 0) * pow(2, -29);
		    	ephemerisStructs[i].i_0 		= intArray_to_complement_bin(subframe+150, subframe+136, 24, 8) * pow(2, -31) * gpsPI;
		    	ephemerisStructs[i].C_rc 		= intArray_to_complement_bin(subframe+180, NULL, 16, 0) * pow(2, -5);
		    	ephemerisStructs[i].omega 		= intArray_to_complement_bin(subframe+210, subframe+196, 24, 8) * pow(2, -31) * gpsPI;
		    	ephemerisStructs[i].omegaDot 		= intArray_to_complement_bin(subframe+240, NULL, 24, 0) * pow(2, -43) * gpsPI;
		    	ephemerisStructs[i].IODE_sf3 		= intArray_to_bin(subframe+270, NULL, 8, 0);
		    	ephemerisStructs[i].iDot 		= intArray_to_complement_bin(subframe+278, NULL, 14, 0) * pow(2, -43) * gpsPI;
                    break;
                case 4:
                    break;
                case 5:
                    break;
                
                default:
                    break;
                }

		// The transmit time for the firset received subframe
		ephemerisStructs[i].TOW = intArray_to_bin(subframe+30, NULL, 17, 0) * 6 - 30;
            }

	    fprintf(fid_eph, "index: %d\n", i);
	    //subframe1
	    fprintf(fid_eph, "WeekNumber: 	%d\n", ephemerisStructs->weekNumber);
	    fprintf(fid_eph, "Accuracy: 	%d\n", ephemerisStructs->accuracy);
	    fprintf(fid_eph, "health: 		%d\n", ephemerisStructs->health);
	    fprintf(fid_eph, "T_GD: 		%e\n", ephemerisStructs->T_GD);
	    fprintf(fid_eph, "IODC: 		%d\n", ephemerisStructs->IODC);
	    fprintf(fid_eph, "t_oc: 		%d\n", ephemerisStructs->t_oc);
	    fprintf(fid_eph, "a_f2: 		%e\n", ephemerisStructs->a_f2);
	    fprintf(fid_eph, "a_f1: 		%e\n", ephemerisStructs->a_f1);
	    fprintf(fid_eph, "a_f0: 		%e\n", ephemerisStructs->a_f0);

	    //subframe2
	    fprintf(fid_eph, "IODE_sf2: 	%d\n", ephemerisStructs->IODE_sf2);
	    fprintf(fid_eph, "C_rs: 		%e\n", ephemerisStructs->C_rs);
	    fprintf(fid_eph, "deltan: 		%e\n", ephemerisStructs->deltan);
	    fprintf(fid_eph, "M_0: 		%e\n", ephemerisStructs->M_0);
	    fprintf(fid_eph, "C_uc: 		%e\n", ephemerisStructs->C_uc);
	    fprintf(fid_eph, "e: 		%e\n", ephemerisStructs->e);
	    fprintf(fid_eph, "C_us: 		%e\n", ephemerisStructs->C_us);
	    fprintf(fid_eph, "sqrtA: 		%e\n", ephemerisStructs->sqrtA);
	    fprintf(fid_eph, "t_oe: 		%d\n", ephemerisStructs->t_oe);

	    //subframe3
	    fprintf(fid_eph, "C_ic: 		%e\n", ephemerisStructs->C_ic);
	    fprintf(fid_eph, "omega_0: 		%e\n", ephemerisStructs->omega_0);
	    fprintf(fid_eph, "C_is: 		%e\n", ephemerisStructs->C_is);
	    fprintf(fid_eph, "i_0: 		%e\n", ephemerisStructs->i_0);
	    fprintf(fid_eph, "C_rc: 		%e\n", ephemerisStructs->C_rc);
	    fprintf(fid_eph, "omega: 		%e\n", ephemerisStructs->omega);
	    fprintf(fid_eph, "omegaDot: 	%e\n", ephemerisStructs->omegaDot);
	    fprintf(fid_eph, "IODE_sf3: 	%d\n", ephemerisStructs->IODE_sf3);
	    fprintf(fid_eph, "iDot: 		%e\n", ephemerisStructs->iDot);

	    fprintf(fid_eph, "\n\n");
        }
    }

    fclose(fid_eph);
}


// calculate the pseudorange 
// trackResults:	store the tracked result
// ms:			
// usdeSatellite:	the tracked satellite list
// receiverSetting:	the setting for receiver
void calculatePseudoranges(trackResult* trackResults, int32_t* ms, uint32_t* usedSatellite, struct settings* receiverSetting, double* pseudoranges, uint32_t usedSatelliteNumber) {
	double samplesPerCode = (receiverSetting->samplingFreq + 0.0) / ((receiverSetting->codeFreqBasis+0.0)/receiverSetting->codeLength);
	double* travleTime = (double*)malloc(sizeof(double)*usedSatelliteNumber);
	for (uint8_t i=0; i<usedSatelliteNumber; i++) {
		// calculate how many code periods have been passed
		travleTime[i] = (trackResults[usedSatellite[i]].absoluteSamples[ms[i]]+1.0) / samplesPerCode;
	}	
	double minimum = travleTime[0];
	for (uint8_t i=1; i<usedSatelliteNumber; i++) {
		if (minimum > travleTime[i]) 
			minimum = travleTime[i];
	}
	minimum = floor(minimum);

	for (uint8_t i=0; i<usedSatelliteNumber; i++) {
		pseudoranges[i] = (travleTime[i] - minimum + receiverSetting->startOffset) * (SPEED_OF_LIGHT/1000);
	}
	free(travleTime);
}

struct posLLH getPos(struct settings receiverSetting) 
{
    struct posLLH result = {0, 0, 0};
    trackResult *trackResults;
    trackResultInit(&trackResults, &receiverSetting);

    FILE *fid_IE = fopen("IE.dat", "r"); 
    FILE *fid_QE = fopen("QE.dat", "r"); 
    FILE *fid_IP = fopen("IP.dat", "r"); 
    FILE *fid_QP = fopen("QP.dat", "r"); 
    FILE *fid_IL = fopen("IL.dat", "r"); 
    FILE *fid_QL = fopen("QL.dat", "r"); 
    FILE *fid_absoluteSamples = fopen("absoluteSamples.dat", "r");

    uint32_t i;	//tracked satellite number
    for (i=0; i<receiverSetting.acqStatelliteList->size; i++) {
        int readNum = fread(trackResults[i].I_E, sizeof(double), receiverSetting.msToProcess, fid_IE);
        if (readNum != receiverSetting.msToProcess) {
            break;
        }
        readNum = fread(trackResults[i].Q_E, sizeof(double), receiverSetting.msToProcess, fid_QE);
        readNum = fread(trackResults[i].I_P, sizeof(double), receiverSetting.msToProcess, fid_IP);
        readNum = fread(trackResults[i].Q_P, sizeof(double), receiverSetting.msToProcess, fid_QP);
        readNum = fread(trackResults[i].I_L, sizeof(double), receiverSetting.msToProcess, fid_IL);
        readNum = fread(trackResults[i].Q_L, sizeof(double), receiverSetting.msToProcess, fid_QL);
        readNum = fread(trackResults[i].absoluteSamples, sizeof(double), receiverSetting.msToProcess, fid_absoluteSamples);
        trackResults[i].flag = 1;
    }
    uint32_t trackedSatelliteNumber = i;

    fclose(fid_IE);
    fclose(fid_QE);
    fclose(fid_IP);
    fclose(fid_QP);
    fclose(fid_IL);
    fclose(fid_QL);
    fclose(fid_absoluteSamples);


    int32_t *pos = tool_find_preambles(trackResults, &receiverSetting);

    // assamble 20 points to 1 bits

    int8_t** bits;

    // get actual bits of GPS signal
    bits = assamble_1ms_to_bit(trackResults, pos, &receiverSetting);

    // analy the ephemeris
    struct ephemeris ephemerisStructs[32];
    get_ephemeris(ephemerisStructs, bits);

    // every ms get one positions
    int maxSubframeStart = pos[0];
    for (uint8_t i=1; i<receiverSetting.acqStatelliteList->size; i++){
	if (maxSubframeStart < pos[i]) 
		maxSubframeStart = pos[i];
    }

    // total time for getting position
    uint32_t currMeasNr = floor(receiverSetting.msToProcess - (maxSubframeStart+1.0)) / receiverSetting.navSolPeriod;
    //uint32_t usedSatellite[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    uint32_t* usedSatellite = (uint32_t*)malloc(sizeof(uint32_t)*trackedSatelliteNumber);
    for (uint32_t i=0; i<trackedSatelliteNumber; i++) {
	    usedSatellite[i] = i;
    }

    double* pseudoranges = (double*)malloc(sizeof(double)*trackedSatelliteNumber);
    double satClkCorr[32];
    //struct ECEFStructure ECEFpos[8];
    struct ECEFStructure* ECEFpos = (struct ECEFStructure*)malloc(sizeof(struct ECEFStructure)*trackedSatelliteNumber);
    double transmitTime = ephemerisStructs[0].TOW;

    FILE* fid_pos = fopen("pos.dat", "w");

    for (uint32_t i=0; i<currMeasNr; i++) {
		calculatePseudoranges(trackResults, pos, usedSatellite, &receiverSetting, pseudoranges, trackedSatelliteNumber);	
	    	for (uint32_t i=0; i<trackedSatelliteNumber; i++) {
			pos[i] += receiverSetting.navSolPeriod;
		}
		satpos(transmitTime, usedSatellite, trackedSatelliteNumber, ephemerisStructs, &receiverSetting, ECEFpos, satClkCorr);


		// Correct the pseudorange
		for (uint32_t i=0; i<trackedSatelliteNumber; i++) {
			pseudoranges[i] += satClkCorr[i]*SPEED_OF_LIGHT;
		}

		// Cauculate the position of satellite
		gsl_vector* pseudorangesVec = gsl_vector_alloc(trackedSatelliteNumber);
		assign_value_to_vec(pseudorangesVec, pseudoranges);

		gsl_vector** pos_vec = (gsl_vector**)malloc(sizeof(gsl_vector*)*trackedSatelliteNumber);
		for (uint32_t i=0; i<trackedSatelliteNumber; i++) {
			pos_vec[i] = gsl_vector_alloc(3);
		}

		for (uint32_t i=0; i<trackedSatelliteNumber; i++) {
			pos_vec[i]->data[0] = ECEFpos[i].X;
			pos_vec[i]->data[1] = ECEFpos[i].Y;
			pos_vec[i]->data[2] = ECEFpos[i].Z;
		}

		gsl_vector* b = pseudorangesVec;
		gsl_matrix* A = gsl_matrix_alloc(trackedSatelliteNumber, 4);

		gsl_vector* pos_r = gsl_vector_alloc(4);
		initial_vector(pos_r, 0);

		gsl_vector* b0 = gsl_vector_alloc(trackedSatelliteNumber);
		gsl_vector* tmp = gsl_vector_alloc(trackedSatelliteNumber);
		initial_vector(tmp, 1);

		double* rho = (double*)malloc(sizeof(double)*trackedSatelliteNumber);

		assign_vector_to_matrix(A, tmp, 3, 0);

		
		for (int i=0; i<10; i++) {
			for (uint32_t i=0; i<trackedSatelliteNumber; i++) {
				rho[i] = norm2Diff(pos_vec[i], pos_r);
			}

			assign_value_to_vec(b0, rho);

			for (uint32_t i=0; i<3; i++) {
				for (uint32_t j=0; j<trackedSatelliteNumber; j++) {
					tmp->data[j] = (pos_r->data[i] - pos_vec[j]->data[i]) / rho[j];
				}
				assign_vector_to_matrix(A, tmp, i, 0);
			}

			gsl_vector* err = gsl_vector_alloc(trackedSatelliteNumber);
			gsl_vector_memcpy(err, b);

			gsl_vector* dtErr = gsl_vector_alloc(trackedSatelliteNumber);
			for (uint32_t i=0; i<trackedSatelliteNumber; i++) {
				dtErr->data[i] = pos_r->data[3];
			}
			gsl_vector_sub(err, b0);
			gsl_vector_sub(err, dtErr);

			gsl_matrix* p_inv_A = p_invert_a_matrix(A);

			gsl_vector* x_ref = multiplication_matrix_and_vec(p_inv_A, err);

			gsl_vector_add(pos_r, x_ref);

			gsl_vector_free(err);
			gsl_vector_free(dtErr);
			gsl_matrix_free(p_inv_A);
			gsl_vector_free(x_ref);
		}
			gsl_vector_free(b);
			gsl_vector_free(b0);
			gsl_vector_free(tmp);
			gsl_matrix_free(A);

		for (uint32_t i=0; i<trackedSatelliteNumber; i++) {
			gsl_vector_free(pos_vec[i]);
		}
			free(pos_vec);

		struct posLLH p = cart2geo(pos_r->data[0], pos_r->data[1], pos_r->data[2], 4);
		gsl_vector_free(pos_r);

		result.heigth += p.heigth;
		result.latitude += p.latitude;
		result.longitude += p.longitude;


		fprintf(fid_pos, "latitude : %f ", p.latitude);
		fprintf(fid_pos, "longitude: %f\n", p.longitude);
		fprintf(fid_pos, "heigtht  : %f\n", p.heigth);

		transmitTime += 1.0*receiverSetting.navSolPeriod/1E3;
    }

    free(pos);
    free(usedSatellite);
    free(pseudoranges);

    fclose(fid_pos);

    result.latitude /= currMeasNr;
    result.longitude /= currMeasNr;
    result.heigth /= currMeasNr;

    return result;
}

