#include "../include/ConstFixedSource.h"


ConstFixedSource::ConstFixedSource(double value) :
_value(value)
{}


ConstFixedSource::~ConstFixedSource()
{
	//nothing dynamic
}

double ConstFixedSource::getValue(const std::vector<double> & coordinates) const
{
	return _value;
}