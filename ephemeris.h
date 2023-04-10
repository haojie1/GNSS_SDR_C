#ifndef __EPHEMERIS_H__
#define __EPHEMERIS_H__

#include <stdint.h>

struct ephemeris {
    // subfram1 parameters
    int32_t weekNumber;
    int8_t accuracy;
    double T_GD;
    int16_t IODC;
    uint8_t health;
    uint32_t t_oc;
    double a_f2;
    double a_f1; 
    double a_f0;

    // subframe2 parameters
    uint8_t IODE_sf2;    
    double C_rs;
    double deltan;
    double M_0;
    double C_uc;
    double  e;
    double C_us;
    double sqrtA;
    uint32_t t_oe;

    // subframe3 parameters 
    double C_ic;
    double omega_0;
    double C_is;
    double i_0;
    double C_rc;
    double omega;
    double omegaDot;
    uint8_t IODE_sf3;
    double iDot;

    // time of the first subframe transmit
    int32_t TOW;
};



#endif 
