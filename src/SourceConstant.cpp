#include "../include/SourceConstant.h"

SourceConstant::SourceConstant()
{
	//should never be called
}


SourceConstant::SourceConstant(double value) : Source()
{
	_value = value;
}

void SourceConstant::sampleSource()
{
	//For a constant source just need to sample which cell it is in
	int cell;
}