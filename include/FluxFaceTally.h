//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : FluxFaceTally.h
//  @ Date : 1/31/2014
//  @ Author : 
//
//


#if !defined(_FLUXFACETALLY_H)
#define _FLUXFACETALLY_H

#include "FaceTally.h"

class FluxFaceTally : public FaceTally
{
public:
	FluxFaceTally();
	FluxFaceTally(int n_angular_bins);
	virtual void incrementScore(double weight, double normal_cosine, double surface_area);  //all information any face tally needs
};

#endif  //_FLUXFACETALLY_H
