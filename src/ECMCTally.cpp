//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : ECMCTally.cpp
//  @ Date : 2/10/2014
//  @ Author : 
//
//


#include "../include/ECMCTally.h"


ECMCTally::ECMCTally() : _angular_moment(1, 1), _spatial_moments(1,2)
{}

ECMCTally::~ECMCTally()
{
	//do nothing
}

void ECMCTally::incrementScores(double weight, double path_length, double normalized_dir_cosine,
	double volume, double normalized_position)
{
	_angular_moment.incrementScore(weight, path_length, normalized_dir_cosine, volume, normalized_position);
	_spatial_moments.incrementScore(weight, path_length, normalized_dir_cosine, volume, normalized_position);
}

double ECMCTally::getAngularMoment(int n_histories)
{
	return _angular_moment.getScore(n_histories, 0,0); //there is only one angular bin
}

std::vector<double> ECMCTally::getSpatialMoments(int n_histories)
{
	return _spatial_moments.getScores(n_histories)[0]; //there is only one angular bin
}

void ECMCTally::reset()
{
	_angular_moment.reset();
	_spatial_moments.reset();
}