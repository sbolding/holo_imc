//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  High Order Solver
//  @ File Name : ElementTally.h
//  @ Date : 1/27/2014
//  @ Author : Simon R Bolding
//
// For the element tallies the first spatial moment is taken in terms of a unitless
// basis function xi.  Thus the normalized_position is the distance from the left edge
// of a cell, i.e.: normalized_position = _position/_element_width;
//
// The path length and the volume need to be passed in units of cm.


#if !defined(_ELEMENTTALLY_H)
#define _ELEMENTTALLY_H

#include "Tally.h"

class ElementTally : public Tally
{
private:
	ElementTally(); //should never be called
	ElementTally(const ElementTally & element_tally); //should never be called

public:
	virtual void incrementScore(double weight, double path_length, double normal_cosine, 
		double volume, double normalized_position) = 0;	//Pure virtual, must be implemented as current or flux tally, normalized position is unitless
	ElementTally(int n_angular_bins, int n_spatial_moments);

};

#endif  //_ELEMENTTALLY_H
