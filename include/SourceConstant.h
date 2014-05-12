#ifndef _SOURCECONSTANT_H
#define _SOURCECONSTANT_H

//THIS CLASS IS NOT ACTUALLY IMPLEMENTED YET

#include "Source.h"

class SourceConstant : public Source
{
private:

	double _value; //constant source value, particles per second, per unit volume
	double _total_source_strength; //constant source value, particles per second
	SourceConstant(); //default constructor, never use

public:

	SourceConstant(double value); //Constant source value constructor
	virtual void sampleSource();  //overwrite sample source function.  Sample is done to particle p

};


#endif