// This class is for a residual source using standard monte carlo sampling of which element
// you are born in, and then sampling from the residual PDF for location within the element.
// The sampling of which bin you are in is based on an alias sampler.
//
//  @ Project : Untitled
//  @ File Name : ResidualSource.h
//  @ Date : 2/8/2014
//  @ Author : 
//
//


#if !defined(_STANDARDRESIDUALSOURCE_H)
#define _STANDARDRESIDUALSOURCE_H

#include "ResidualSource.h"
#include "Source.h"
#include "AliasSampler.h"

class StandardResidualSource : public ResidualSource
{
protected:

	AliasSampler* _element_source;  //Sampler to determine which source you are in
	AliasSampler* _face_source;	    //Sampler to determine which face you are on

public:

	~StandardResidualSource();
	StandardResidualSource(Particle1D* particle, FixedSourceFunctor & q);
	virtual void sampleSourceParticle();

};







#endif // _STANDARDRESIDUALSOURCE_H
