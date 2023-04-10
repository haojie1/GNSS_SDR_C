#ifndef __CART2GEO_H__
#define __CART2GEO_H__

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

struct posLLH {
	double latitude, longitude, heigth;
};


struct posLLH cart2geo(double X, double Y, double Z, uint32_t i) ;

#endif
