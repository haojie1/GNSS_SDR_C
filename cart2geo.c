#include "cart2geo.h"

struct posLLH cart2geo(double X, double Y, double Z, uint32_t i) {
	struct posLLH pos;

	const double a[5] = {6378388, 6378160, 6378135,  6378137, 6378137};
	const double f[5] = {1.0/297, 1.0/298.247, 1.0/298.26, 1.0/298.257222101, 1.0/298.257223563};

	double lambda = atan2(Y, X);
	double ex2 = (2-f[i])*f[i] / (pow((1-f[i]), 2));
	double c = a[i]*sqrt(1+ex2);
	double phi = atan(Z/((sqrt(X*X+Y*Y)*(1-(2-f[i]))*f[i])));

	double h = 0.1, oldh = 0;

	while (fabs(h-oldh) > 1e-12) {
		oldh = h;
		double N = c/sqrt(1+ex2*cos(phi)*cos(phi));
		phi = atan(Z/((sqrt(X*X+Y*Y)*(1-(2-f[i])*f[i]*N/(N+h)))));
		h = sqrt(X*X+Y*Y)/cos(phi)-N;
	}

	phi = phi * 180 / M_PI;
	lambda = lambda * 180 / M_PI;

	pos.latitude = phi;
	pos.longitude = lambda;
	pos.heigth = h;

	return pos;
}
