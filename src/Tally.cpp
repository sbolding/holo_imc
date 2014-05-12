//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  High Order Solver
//  @ File Name : Tally.cpp
//  @ Date : 1/27/2014
//  @ Author : Simon R Bolding
//



#include <stdlib.h>
#include "../include/Tally.h"

Tally::~Tally()
{
	//no dynamic members, no need to delete anything
}

double Tally::getScore(int n_histories, int angular_bin, int spatial_moment) const
{
	if (spatial_moment >= _bin_sums[0].size() || angular_bin >= _bin_sums.size())
	{
		std::cerr << "Trying to access tally member that is not available\n";
		exit(1);
	}

	return _bin_sums[angular_bin][spatial_moment] / (float)n_histories;
}

//For this function it is assumed that there is only the 0th spatial moment (for face tallies)
double Tally::getScore(int n_histories, int angular_bin) const
{
	if (angular_bin >= _bin_sums.size())
	{
		std::cerr << "Trying to access tally member that is not available\n";
		exit(1);
	}

	return _bin_sums[angular_bin][0] / (float)n_histories;
}

double Tally::getScoreAngularIntegrated(int n_histories, int spatial_moment) const
{
	if (spatial_moment >= _bin_sums[0].size())
	{
		std::cerr << "Trying to access tally member that is not available\n";
		exit(1);
	}
	double sum = 0.;
	for (int i = 0; i < _bin_sums.size(); ++i)
	{
		sum += _bin_sums[i][spatial_moment];
	}
		
	return sum / (float)n_histories;
}

double Tally::getScoreAngularIntegrated(int n_histories) const
{
	double sum = 0.;
	for (int i = 0; i < _bin_sums.size(); ++i)
	{
		sum += _bin_sums[i][0];
	}

	return sum / (float)n_histories;
}

std::vector<std::vector<double>> Tally::getScores(int n_histories) const
{
	std::vector<std::vector<double>> scores;
	std::vector<double> temp_scores;
	for (int i = 0; i < _bin_sums.size(); ++i) //loop over angular bins, from - to +
	{
		temp_scores.clear();
		for (int j = 0; j < _bin_sums[i].size(); ++j) //loop over spatial bins, from 0 to ...
		{
			temp_scores.push_back(getScore(n_histories,i,j));
		}
		scores.push_back(temp_scores);
	}

	return scores;
}

double Tally::getStdDev(int n_histories, int angular_bin, int spatial_moment) const
{
	if (spatial_moment >= _bin_sums[0].size() || angular_bin >= _bin_sums.size())
	{
		std::cerr << "Trying to access tally member that is not available\n";
		exit(1);
	}

	double variance = (_bin_sums_sq[angular_bin][spatial_moment] / (float)n_histories) -
		pow(getScore(n_histories, angular_bin, spatial_moment), 2);
	variance /= (n_histories - 1);
	
	return sqrt(variance);
}

double Tally::getStdDev(int n_histories, int angular_bin) const
{
	if (angular_bin >= _bin_sums.size())
	{
		std::cerr << "Trying to access tally member that is not available\n";
		exit(1);
	}

	return getStdDev(n_histories, angular_bin, 0);
}

std::vector<std::vector<double>> Tally::getStdDevs(int n_histories) const
{
	std::vector<std::vector<double>> std_devs;
	std::vector<double> temp_std_devs;
	for (int i = 0; i < _bin_sums.size(); ++i)
	{
		temp_std_devs.clear();
		for (int j = 0; j < _bin_sums[i].size(); ++j)
		{
			temp_std_devs.push_back(getStdDev(n_histories, i, j));
		}
		std_devs.push_back(temp_std_devs);
	}

	return std_devs;
}


Tally::Tally(int n_angle_bins, int n_spatial_moment_bins)
{
	//set all tally bins to zero
	_bin_sums.resize(n_angle_bins);
	_bin_sums_sq.resize(n_angle_bins);
	for (int i = 0; i < n_angle_bins; ++i)
	{
		_bin_sums[i].assign(n_spatial_moment_bins, 0.);
		_bin_sums_sq[i].assign(n_spatial_moment_bins, 0.);
	}
}

void Tally::reset()
{
	//set all tally bins to zero
	for (int i = 0; i < _bin_sums.size(); ++i)
	{
		std::fill(_bin_sums[i].begin(), _bin_sums[i].end(), 0.0);
		std::fill(_bin_sums_sq[i].begin(), _bin_sums_sq[i].end(), 0.0);
	}
}

Tally::Tally()
{
	//Default constructor
}

