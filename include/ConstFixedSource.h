#pragma once
#include "FixedSourceFunctor.h"

class ConstFixedSource :
	public FixedSourceFunctor
{
protected:
	double _value;
public:
	ConstFixedSource(double value); //value must be in p/(sec-str-cm^3)
	~ConstFixedSource();
	virtual double getValue(const std::vector<double> & coordinates) const; //only returns value, doenst care about vector size
};

