#include "tool.h"

uint32_t next2pow(uint32_t data) {
    if ((data & (data-1)) == 0) {
        return data;
    } else {
        uint32_t top = 0x80000000;
        while ((top & data) == 0) {
            top >>= 1;
        }
        top <<= 1;

        return top;
    }
}

double* tool_padding_data_packed_double(double* data, uint32_t originalN, uint32_t N) {
    if (originalN > N) {
        printf("Error In tool_padding_data\n");
        exit(EXIT_FAILURE);
	return NULL;
    } else {
        double *newData = (double*)malloc(sizeof(double)*2*N);
        memcpy(newData, data, sizeof(double)*originalN*2);
        for (int i=originalN; i<N; i++) {
            REAL(newData, i) = 0;
            IMAG(newData, i) = 0;
        }
        return newData;
    }
}


