//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : CurrentElementTally.cpp
//  @ Date : 1/31/2014
//  @ Author : 
//
//


#include <stdlib.h>
#include "../include/CurrentElementTally.h"
#include <iostream>

CurrentElementTally::CurrentElementTally(int n_angular_bins, int n_spatial_moments) : ElementTally(n_angular_bins, n_spatial_moments)
{}

void CurrentElementTally::incrementScore(double weight, double path_length,
	double normal_cosine, double volume, double normalized_position)
{
	double value;  //tally value
	int angular_bin;
	//Determine which bin you are in (currently it assumes only two)
	if (_bin_sums.size() == 2)
	{
		angular_bin = floor(normal_cosine) + 1; //either 0 or 1
	}
	else if(_bin_sums.size() == 1)
	{
		angular_bin = 0;
	}
	else
	{
		std::cerr << "Tallies only handle 1 and 2 angular directions currently \n";
		exit(1);
	}
	//loop over spatial moments, recursively, zeroth bin is 0th moment
	value = weight*normal_cosine*path_length / volume; //0th spatial moment
	if (path_length*normal_cosine < 0.)
	{
		std::cerr << "You have passed in a negative pathlength or cosine to the tally, this may not be correct for ECMC case, CurrentElementTally.cpp";
		system("pause");
		exit(1);
	}
	for (int i = 0; i < _bin_sums[0].size(); ++i)
	{
		_bin_sums[angular_bin][i] += value;
		_bin_sums_sq[angular_bin][i] += value*value;
	    value *= normalized_position; //get spatial moment for next iterat.
	}

}