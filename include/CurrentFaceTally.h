//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : CurrentFaceTally.h
//  @ Date : 1/31/2014
//  @ Author : 
//
//


#if !defined(_CURRENTFACETALLY_H)
#define _CURRENTFACETALLY_H

#include "FaceTally.h"

class CurrentFaceTally : public FaceTally
{
public:
	CurrentFaceTally();
	CurrentFaceTally(int n_angular_bins);
	virtual void incrementScore(double weight, double normal_cosine, double surface_area);

};

#endif  //_CURRENTFACETALLY_H
