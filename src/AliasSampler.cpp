#include <stdlib.h>
#include "../include/AliasSampler.h"
#include <vector>
#include <numeric>
#include <algorithm>
#include <limits>

AliasSampler::AliasSampler(std::vector<double> bin_probabilities, bool normalized)
{
	if (!normalized) //then normalize bin_probabilities
	{
		double sum = 0;
		sum = std::accumulate(bin_probabilities.begin(), bin_probabilities.end(), 0.);
		double inv_sum = 1. / sum; //invert sum for cheaper multiply
		for (int i = 0; i < bin_probabilities.size(); ++i)
		{
			if (bin_probabilities[i] < 0.0)
			{
				std::cerr << "You have passed in a negative probability to AliasSampler\n";
				system("pause");
				exit(1);
			}
			bin_probabilities[i] *= inv_sum;
		}
	} 

	//create the list of probabilities
	_n_bins = bin_probabilities.size();
	createAliasTable(bin_probabilities);
}

AliasSampler::AliasSampler()
{
	//default do nothing
}

void AliasSampler::createAliasTable(std::vector<double> & bin_probs)
{
	//local variables
	std::vector<unsigned int> lo_bins; //bins that have probability below avg
	std::vector<unsigned int> hi_bins; //bins that have probability above avg
	double p_avg = 1. / bin_probs.size();
	double inv_p_avg = (float)_n_bins; // = 1/p_avg = _n_bins

	_alias_data.resize(_n_bins); //initialize list

	//Check that the probability is not already uniform, if so do special case
	//This should be redundant, as the later algorithm should catch this case, but in cases of very high numbers of bins coudl be needed due to large roundoff
	bool all_bins_equal = true;
	double hundred_epsilon = 100. * std::numeric_limits<double>::epsilon();
	for (int i = 0; i < bin_probs.size(); ++i)
	{
		if (abs(p_avg - bin_probs[i]) > hundred_epsilon)
		{
			all_bins_equal = false;
			break;
		}
	}

	//Do case if all_bins_equal
	if (all_bins_equal)
	{
		for (int i = 0; i < bin_probs.size(); ++i)
		{
			_alias_data[i].non_alias_probability = 1.0001;
		}
		return;
	}


	//Sort the bins, this could in theory have trouble if one is approx the average, but shouldnt since
	//bins that are "at probability" go in to the higher bin
	for (size_t bin=0; bin < bin_probs.size(); bin++)
	{
		if (bin_probs[bin] < p_avg)
		{
			lo_bins.push_back(bin);
		}
		else
		{
			hi_bins.push_back(bin);
		}
	}

	//Start filling the table by looping over the small ones
	int n_small = lo_bins.size()-1; //these are just to save access time getting the size over and over
	unsigned int i_small, i_big; //index in the actual array of the current big and small bin
	int n_big = hi_bins.size()-1;

	while (n_small >= 0)
	{
		i_small = lo_bins[n_small]; 		//figure out what the current bins are
		i_big = hi_bins[n_big];

		_alias_data[i_small].alias_event = i_big;	//where the extra prob came from
		_alias_data[i_small].non_alias_probability = bin_probs[i_small] * inv_p_avg; // p_{i_small}/p_avg
		lo_bins.pop_back();	//get rid of last element since dont need it anymore
		n_small--;

		bin_probs[i_big] -= (p_avg - bin_probs[i_small]); //reduce probability of the big bin

		//Determine if the current hi_bin goes into the lo_bins
		if (bin_probs[i_big] < p_avg)
		{
			if (n_big < 1 ) //don't get rid of the last big bin, this is special case because of roundoff error
			{
				if (n_small < 0) //Still have more lo bins to initialize, again, because of roundoff
				{
					break;
				}
				else
				{
					continue;
				}
			}
			lo_bins.push_back(i_big); //move hi bin to lo bin
			hi_bins.pop_back();	//delete hi bin
			n_big--;
			n_small++;
		}
		//else, do nothing, index is correct
	}
	//Any ho_bins that are left need to be initialized
	for (int i = 0; i < (n_big+1); i++)
	{
		_alias_data[hi_bins[i]].non_alias_probability = 1.0001; //ensure alias bin is never accessed
		_alias_data[hi_bins[i]].alias_event = hi_bins[i];
	}

}

void HistogramUtilities::normalizeDiscreteDistribution(std::vector<double> & bin_probabilities)
{
	double sum = 0;
	sum = std::accumulate(bin_probabilities.begin(), bin_probabilities.end(), 0.);
	double inv_sum = 1. / sum; //invert sum for cheaper multiply
	for (int i = 0; i < bin_probabilities.size(); ++i)
	{
		bin_probabilities[i] *= inv_sum;
	}
}

void HistogramUtilities::scaleDiscreteDistribution(std::vector<double> & bin_probabilities, double scalar)
{
	for (int i = 0; i < bin_probabilities.size(); ++i)
	{
		bin_probabilities[i] *= scalar;
	}
}