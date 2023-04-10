#ifndef __TOOL_H__
#define __TOOL_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "config.h"

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

uint32_t next2pow(uint32_t data) ;
double* tool_padding_data_packed_double(double* data, uint32_t originalN, uint32_t N) ;
uint32_t next2pow(uint32_t data) ;

#endif