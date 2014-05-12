#include "../include/MMSFixedSource.h"


MMSFixedSource::MMSFixedSource()
{
	//Currently does nothing
}

MMSFixedSource::MMSFixedSource(double coeff, double x_coeff, double mu_coeff) :
_a(coeff),
_a_x(x_coeff),
_a_mu(mu_coeff)
{};

MMSFixedSource::~MMSFixedSource()
{
}

double MMSFixedSource::getValue(const std::vector<double> & coors) const
{
	return _a + coors[0] * _a_x + coors[1] * _a_mu;
}
