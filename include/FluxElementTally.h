//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : FluxElementTally.h
//  @ Date : 1/31/2014
//  @ Author : 
//
//


#if !defined(_FLUXELEMENTTALLY_H)
#define _FLUXELEMENTTALLY_H

#include "ElementTally.h"

class FluxElementTally : public ElementTally
{
public:
	FluxElementTally();
	FluxElementTally(int n_angular_bins, int n_spatial_moments);

	virtual void incrementScore(double weight, double path_length, double normal_cosine,
		double volume, double normalized_position);
};

#endif  //_FLUXELEMENTTALLY_H
