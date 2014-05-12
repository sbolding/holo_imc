//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  High Order Solver
//  @ File Name : FaceTally.h
//  @ Date : 1/27/2014
//  @ Author : Simon R Bolding
//
// Currently it is hardcoded to have 2 angular bins and take the 1st and second
// spatial moments
//


#if !defined(_FACETALLY_H)
#define _FACETALLY_H

#include "Tally.h"
#include "GlobalConstants.h"

using HoConstants::COSINE_CUTOFF;
using HoConstants::COSINE_SUBSTITUTE_VALUE;

class FaceTally : public Tally
{
private:
	FaceTally(); //Never used default constructor

public:
	virtual void incrementScore(double weight, double normal_cosine, double surface_area) = 0; //Pure virtual, it must be implemented by current or flux tally
	FaceTally(int n_angular_bins);
};

#endif  //_FACETALLY_H
