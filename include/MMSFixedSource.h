//This class could be improved by overloading the constructor to 
//take in a function pointer for more arbitrary evaluation
//of the source function.

#pragma once
#include "FixedSourceFunctor.h"

class MMSFixedSource : 	public FixedSourceFunctor
{
private:
	double _a;
	double _a_x;
	double _a_mu;

public:

	MMSFixedSource();
	MMSFixedSource(double scalar_coeff, double x_coeff, double mu_coeff);
	~MMSFixedSource();
	virtual double getValue(const std::vector<double> & coors) const; //evaluates function
};

