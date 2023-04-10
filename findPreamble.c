#include "findPreamble.h"

void tool_calculate_correlation(int32_t* correlation_value, int32_t* bits, uint32_t len1, int32_t* preamble_bits_upsampling, uint32_t len2) {
   // this version doesn't use fft 
   // 1. check whether the size of the two sequences is same
   int32_t * padding_compared_sequence = (int32_t *)malloc(sizeof(int32_t)*len1);
   memset(padding_compared_sequence, 0, sizeof(int32_t)*len1);
   memcpy(padding_compared_sequence, preamble_bits_upsampling, sizeof(int32_t)*len2);

   //2. start cauculate the correlation value
   for (uint32_t i=0; i<(len1-len2); i++) {
    correlation_value[i] = 0;
    for (uint32_t j=0; j<len2; j++) {
        correlation_value[i] += bits[i+j]*padding_compared_sequence[j];
    }
   }

   for (uint32_t i=len1-len2, k=0; i<len1; i++, k++) {
    correlation_value[i] = 0;
    for (uint32_t j=0; j<len2-k; j++) {
        correlation_value[i] += bits[i+j]*padding_compared_sequence[j];
    } 
   }

   free(padding_compared_sequence);
}

int32_t* tool_find_preambles(trackResult* trackResults, struct settings* receiverSettings) {
   int32_t preamble_bits[8] = {1, -1, -1, -1, 1, -1, 1, 1}; 

   // here is fixed to GPS signal
   int32_t *preamble_bits_upsampling = (int32_t *)malloc(sizeof(int32_t)*8*20);
   int32_t *sub_frame_pos = (int32_t*)malloc(sizeof(int32_t)*32);

   // upsamping the preamble_bits to suit the coherent period
   for (uint32_t i=0; i<8*20; i++) {
    preamble_bits_upsampling[i] = preamble_bits[i/20];
   }

   int32_t* bits = (int32_t*)malloc(sizeof(int32_t)*receiverSettings->msToProcess);
   int32_t* correlation_value = (int32_t*)malloc(sizeof(int32_t)*receiverSettings->msToProcess);
   uint32_t *pos = (uint32_t*)malloc(sizeof(uint32_t)*100);

   // deal with every tracked satellite
   for (uint32_t i=0; i<receiverSettings->acqStatelliteList->size; i++) {
    int32_t sub_frame_pos_now = -1;

    if (trackResults[i].flag == 1) {
        // make the bits be -1 or 1
        for (int j=0; j<receiverSettings->msToProcess; j++) {
            if (trackResults[i].I_P[j] > 1) {
                bits[j] = 1;
            } else {
                bits[j] = -1;
            }
        }

        // calculate the correlation value
        tool_calculate_correlation(correlation_value, bits, receiverSettings->msToProcess, preamble_bits_upsampling, 20*8);

        // find the max value pos
        uint32_t potential_pos_num = 0;
        for (int i=0; i<receiverSettings->msToProcess; i++) {
            if (abs(correlation_value[i]) > 150) {
                pos[potential_pos_num++] = i;
            }
        }


        // find the true sub frame head pos

        for (uint32_t k=0; k<potential_pos_num; k++) {
            uint32_t num = 0;
            for (uint32_t j=k+1; j<potential_pos_num; j++) {
                if ((pos[j]-pos[k])%6000 == 0) {
                    num++;
                }
                if (num == 5) {
                    sub_frame_pos_now = pos[k];
                    goto next_loop;
                }
            }
        }
    }

next_loop:
    sub_frame_pos[i] = sub_frame_pos_now;
   }

    free(correlation_value);
    free(bits);
    free(pos);
    free(preamble_bits_upsampling);

   return sub_frame_pos;
}
