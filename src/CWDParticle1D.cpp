#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include "../include/CWDParticle1D.h"


CWDParticle1D::CWDParticle1D(HoMesh* mesh, Source* src, RNG* rng, string method_str) : Particle1D(mesh, src, rng, method_str)
{
	_is_analog = false;
	_non_analog_mfp = HoConstants::MAX_STREAMING_MFP; 
	if (_method != HoMethods::HOLO_ECMC)
	{
		std::cerr << "You can not make a continuous weight deposition particle with standard monte carlo\n";
		exit(1);
	}
}


void CWDParticle1D::scoreElementTally(double path_start_mfp, double path_end_mfp)
{
	if (_is_analog)
	{
		//call the parent classes element tally score, no need to adjust weight
		Particle1D::scoreElementTally(path_start_mfp, path_end_mfp);
		return;
	}
	else
	{
		//Score Element tally, need to convert path_length and volume
		//to cm, rather than mfp
		double path_length_mfp = std::abs((path_start_mfp - path_end_mfp) / _mu);
		double path_length_cm = path_length_mfp*_mfp_tot;
		double normalized_position = (path_start_mfp + path_end_mfp) / _element_width_mfp; //location of the center of pathlength
		

		//calculate the new weight and the amount of weight deposited
		double weight_new = _weight*exp(-path_length_mfp);
		double weight_deposited = weight_new - _weight;

		//ScoreECMCTallies, using a normalized direction cosine to ensure positive tallies
		std::cout << "This does not work yet, in CWDParticle1D.cpp\n";
		exit(1);
		_current_element->incrementTallyScores(weight_deposited, _mfp_tot,
			_mu, normalized_position);

		//update the weight
		_weight = weight_new;
	}
}

void CWDParticle1D::runHistory()
{
	sampleSourceParticle(); //samples position, direction, and initalizes weight
	streamToNextEvent(_non_analog_mfp); //stream until leaked or destroyed
	if (_is_dead)
	{
		return;
	}
	else //switch to analog transport
	{
		_is_analog = true;
		double path_length_mfp = samplePathLengthMFP(); //sample a pathlength
		streamToNextEvent(path_length_mfp); //stream to leakage or absorption
	}
}
