//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : CurrentFaceTally.cpp
//  @ Date : 1/31/2014
//  @ Author : 
//
//


#include <stdlib.h>
#include "../include/CurrentFaceTally.h"

CurrentFaceTally::CurrentFaceTally(int n_angular_bins) : FaceTally(n_angular_bins)
{}

void CurrentFaceTally::incrementScore(double weight, double normal_cosine,
	double area)
{
	//Determine which bin you are in (currently it assumes only two)
	if (_bin_sums.size() == 2)
	{
		int angular_bin = floor(normal_cosine) + 1; //either 0 or 1
		double value = weight / area;
		_bin_sums[angular_bin][0] += value;
		_bin_sums_sq[angular_bin][0] += value*value;
	}
	else
	{
		std::cerr << "Tallies only handle two directions currently \n";
		exit(29);
	}
}