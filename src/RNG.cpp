#include "../include/RNG.h"
#include "../include/Controller.h"

RNG::RNG()
{
	a = 16807;
	M = 2147483647;
	Q = 127773;
	R = 2836;
	Minv = (double)1. / M;
	seed = HoController::INPUT_SEED;
}

double RNG::rand_num()
{
	k = (seed / Q);
	seed = a*(seed % Q) - R*k;
	if (seed < 0) seed += M;
	return double(seed*Minv);
}