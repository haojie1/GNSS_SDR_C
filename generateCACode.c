#include "generateCACode.h"

const uint32_t g2s[32] = {
    5, 6, 7, 8, 17, 18, 139, 140,
    141, 251, 252, 254, 255, 256, 257, 258,
    469, 470, 471, 472, 473, 474, 509, 512,
    513, 514, 515, 516, 859, 860, 861, 862
};

int8_t g1[1023];
int8_t g2[1023];


int8_t reg[10] ;

void CACode(uint8_t prn, int8_t *caCode) {
    // Clear reg
    memset(reg, -1, sizeof(int8_t)*10);

    // -1 stands for 1
    reg[9] = -1;

    int8_t temp;
    for (int i=0; i<1023; i++) {
        g1[i] = reg[9];
        temp = reg[2]*reg[9]; 
        memmove(reg+1, reg, sizeof(int8_t)*9); 
        reg[0] = temp;
    }

//    for (int i=0; i<1023; i++) {
//        printf("%2d", g1[i]);
//    }

    // create G2 code
    memset(reg, -1, sizeof(int8_t)*10);
    for (int i=0; i<1023; i++) {
        g2[i] = reg[9];
        temp = reg[1]*reg[2]*reg[5]*reg[7]*reg[8]*reg[9]; 
        memmove(reg+1, reg, sizeof(int8_t)*9); 
        reg[0] = temp;
    }

    // shift G2 Code
    int8_t buffForG2[1023];
    memcpy(buffForG2, g2, sizeof(int8_t)*(1023-g2s[prn-1]));
    memmove(g2, g2+(1023-g2s[prn-1]), g2s[prn-1]);
    memcpy(g2+g2s[prn-1], buffForG2, sizeof(int8_t)*(1023-g2s[prn-1]));

    for (int i=0; i<1023; i++) {
        *(caCode+i) = -1*g1[i]*g2[i];
    }
}


void caCodeAfterSampling(const struct settings* receiverSetting, uint32_t number, int8_t *caCodeSamplings, uint8_t prn) {
    int8_t caCode[1023];
    uint32_t index;
    CACode(prn, caCode);

    //FILE* fid = fopen("index.dat", "w");
    for (uint32_t i=0; i<number; i++) {
        // Make the first sample is 1/fs not 0
        if (fmod((i+1.0)*receiverSetting->codeFreqBasis, receiverSetting->samplingFreq) == 0) {
            // if change the multilication and division may change the result
            index = (i+1.0)*receiverSetting->codeFreqBasis/receiverSetting->samplingFreq-1;
        } else {
            index = (i+1.0)*receiverSetting->codeFreqBasis/receiverSetting->samplingFreq;
        }
     //   fwrite(&index, sizeof(uint32_t), 1, fid);
	index = fmod(index, 1023);
        *(caCodeSamplings+i) = caCode[index]; 
    }
    //fclose(fid);

    // correct the last number
    //*(caCodeSamplings+number-1) = caCode[1022];
}
