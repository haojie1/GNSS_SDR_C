#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "initSetting.h"

void CACode(uint8_t prn, int8_t *caCode) ;
void caCodeAfterSampling(const struct settings* receiverSetting, uint32_t number, int8_t *caCodeSamplings, uint8_t prn) ;