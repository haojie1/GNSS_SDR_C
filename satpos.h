#ifndef __SATPOS_H__
#define __SATPOS_H__

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "config.h"
#include "ephemeris.h"

void satpos(double transmitTime, uint32_t* usedSatellite, uint8_t usedSatelliteNumber, struct ephemeris* ephemerisStructs, 
		struct settings* receiverSetting, struct ECEFStructure* ECEFpos, double* satClkCorr) ;

#endif
