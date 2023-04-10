#include "satpos.h"
#define 	Omegae_dot 	(7.2921151467e-5)
#define 	GM		(3.986005e14)
#define 	F		(-4.442807633e-10)

double check_t(double t) {

	int32_t half_week = 302400;
	double corrTime = t;

	if (t > half_week) {
		corrTime = t - 2*half_week;
	} else if(t < -half_week) {
		corrTime = t + 2*half_week;
	}

	return corrTime;
}

void satpos(double transmitTime, uint32_t* usedSatellite, uint8_t usedSatelliteNumber, struct ephemeris* ephemerisStructs, 
		struct settings* receiverSetting, struct ECEFStructure* ECEFpos, double* satClkCorr) {
	for (uint8_t k=0; k<usedSatelliteNumber; k++) {
		uint8_t satelliteIndex = usedSatellite[k];
		int32_t dt = check_t(transmitTime - ephemerisStructs[satelliteIndex].t_oc);
		//calculate clock correction
		satClkCorr[satelliteIndex] = (ephemerisStructs[satelliteIndex].a_f2 * dt + ephemerisStructs[satelliteIndex].a_f1) * dt + 
				ephemerisStructs[satelliteIndex].a_f0 - ephemerisStructs[satelliteIndex].T_GD;
		double time = transmitTime - satClkCorr[satelliteIndex];


		//find satellite's position
		//Restore semi-major axis
		double a = ephemerisStructs[satelliteIndex].sqrtA * ephemerisStructs[satelliteIndex].sqrtA;
		//Time correction
		double tk  = check_t(time - ephemerisStructs[satelliteIndex].t_oe);
		//Initial mean motion
		double n0 = sqrt(GM/pow(a, 3));
		//Mean motion 
		double n = n0 + ephemerisStructs[satelliteIndex].deltan;
		//Mean anomaly
		double M = ephemerisStructs[satelliteIndex].M_0 + n * tk;
		//Reduce mean anomyly to between 0 and 360 deg
		M = fmod(M + 2*gpsPI, 2*gpsPI);
		//Initial guess of eccentric anomaly
		double E = M;

		//Iteratively compute eccentric anomaly
		double E_old, dE;
		for (uint8_t i=0; i<10; i++) {
			E_old = E;
			E = M + ephemerisStructs[satelliteIndex].e * sin(E);
			dE = fmod(E - E_old, 2*gpsPI);

			if (fabs(dE) < 1e-12){
				break;
			}
		}

		//Reduce eccentric anomaly to between 0 and 360 deg
		E = fmod(E+2*gpsPI, 2*gpsPI);

		//Compute relativistic correction term
		double dtr = F * ephemerisStructs[satelliteIndex].e * ephemerisStructs[satelliteIndex].sqrtA * sin(E);

		//Calculate the true anomaly
		double nu = atan2(sqrt(1-pow(ephemerisStructs[satelliteIndex].e, 2)) * sin(E), cos(E) - ephemerisStructs[satelliteIndex].e);

		//Compute angle phi
		double phi = nu + ephemerisStructs[satelliteIndex].omega;
		phi = fmod(phi, 2*gpsPI);

		//Correct argument of latitude
		double u = phi + 
				ephemerisStructs[satelliteIndex].C_uc * cos(2*phi) + 
				ephemerisStructs[satelliteIndex].C_us * sin(2*phi);
		//Correct radius
		double r = a * (1 - ephemerisStructs[satelliteIndex].e*cos(E)) + 
				ephemerisStructs[satelliteIndex].C_rc * cos(2*phi) + 
				ephemerisStructs[satelliteIndex].C_rs * sin(2*phi);
		//Correct inclination
		double i = ephemerisStructs[satelliteIndex].i_0 +  ephemerisStructs[satelliteIndex].iDot * tk + 
			  	ephemerisStructs[satelliteIndex].C_ic * cos(2*phi) + 
				ephemerisStructs[satelliteIndex].C_is * sin(2*phi);

		//Compute the angle between tha ascending node and the Greenwich meridian
		double Omega = ephemerisStructs[satelliteIndex].omega_0 + (ephemerisStructs[satelliteIndex].omegaDot- Omegae_dot) * tk - 
			Omegae_dot * ephemerisStructs[satelliteIndex].t_oe;

		//Reduce to between 0 and 360 deg
		Omega = fmod(Omega + 2*gpsPI, 2*gpsPI);

		//Compute satellite coordinates
		ECEFpos[satelliteIndex].X = cos(u) * r * cos(Omega) - sin(u) * r * cos(i) * sin(Omega);
		ECEFpos[satelliteIndex].Y = cos(u) * r * sin(Omega) + sin(u) * r * cos(i) * cos(Omega);
		ECEFpos[satelliteIndex].Z = sin(u) * r * sin(i);
		

		satClkCorr[satelliteIndex] = (ephemerisStructs[satelliteIndex].a_f2 * dt + ephemerisStructs[satelliteIndex].a_f1) * dt + 
						ephemerisStructs[satelliteIndex].a_f0 - 
						ephemerisStructs[satelliteIndex].T_GD + dtr;
		

	}

}
