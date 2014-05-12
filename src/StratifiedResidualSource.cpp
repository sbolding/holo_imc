#include "../include/StratifiedResidualSource.h"
#include "../include/Particle1D.h" //Because a friend class, MUST include this here and not in the header file
#include <cmath>
#include "../include/AliasSampler.h"

StratifiedResidualSource::StratifiedResidualSource(Particle1D* particle,
  FixedSourceFunctor & q, unsigned int & n_histories) : 
ResidualSource(particle, q)
{
	//determine how many particles to sample from each element
	unsigned int n_active_elements = _particle->_mesh->getNumActiveElements();
	_last_el_id = _particle->_mesh->getNumElems() - 1;
	while (n_histories % n_active_elements != 0)
	{
		n_histories++;
	}
	_n_samples_per_element = n_histories / (int)n_active_elements;

	if (_n_samples_per_element < 100)
	{
		std::cout << "DEBUG less than 100: " << _n_samples_per_element << std::endl;
		std::cout << "n_histories" << n_histories << " N_elements" << n_active_elements << std::endl;
		system("pause");
	}

	//initialize the element id, the next two lines will lead to starting on the 0th element the first time through
	_current_element_id = -1; //always start at the first element, if inactive it will be determined later
	_n_sampled_from_current_element = _n_samples_per_element;

	//normalize the face and element magnitudes to form a pdf of a particle being born in each element (either face or element)
	double norm_const = 1./(_face_src_total + _vol_src_total);
	HistogramUtilities::scaleDiscreteDistribution(_res_face_mags, norm_const);
	HistogramUtilities::scaleDiscreteDistribution(_res_element_mags, norm_const);

	//determine the stratified PDF's bin probability.  Here it is assumed uniform
	_stratified_probability = 1. / (double)n_active_elements;
}

StratifiedResidualSource::~StratifiedResidualSource()
{
	//no dynamic memory
}

void StratifiedResidualSource::sampleSourceParticle()
{
	//In the event of zero probability bins don't want to sample all of them
	if (_current_element_id > _last_el_id)
	{
		_particle->terminateHistory();
		return;
	}

	//Determine which element you are in
	if (_n_sampled_from_current_element == _n_samples_per_element)
	{
		//find new active element (has no children)
		_current_element_id++;
		while (true)
		{
			//If sampled all elements stop sampling, set particle to dead
			if (_current_element_id == _last_el_id + 1)
			{
				_particle->terminateHistory();
				return;
			}
			else if (!_particle->_mesh->getElement(_current_element_id)->hasChildren())
			{
				//if the residual for an element is zero, need to forego sampling
				if (_res_face_mags[_current_element_id] < 0.0001*GlobalConstants::RELATIVE_TOLERANCE) //source is normalized to one, so accurate 
				{
					if (_res_element_mags[_current_element_id] < 0.0001*GlobalConstants::RELATIVE_TOLERANCE)
					{
						_current_element_id++; 
						//DEBUG STUFF
						//std::cout << "Element has zero residual, ID =  " << _current_element_id - 1 << "FACE: " << _res_element_mags[_current_element_id-1]
							//<< " ELEMENT: " << _res_element_mags[_current_element_id-1] << std::endl;
						continue;					
					}
				}
				//found a working element
				_n_sampled_from_current_element = 0; //reset counter
				_curr_el_face_probability = _res_face_mags[_current_element_id];
				_curr_el_element_probability = _res_element_mags[_current_element_id];
				break;
			}
			else //parent element
			{
				_current_element_id++;
			}
		}
	}  //else, stay on the current element

	//set particle in cell
	_particle->_current_element_ID = _current_element_id;
	_particle->_current_element = _particle->_mesh->getElement(_current_element_id);
	_particle->updateElementProperties();

	//Sample if face or element source
	if (_rng->rand_num()*(_curr_el_face_probability + _curr_el_element_probability) < _curr_el_element_probability) //element source
	{
		//determine location and direction within element using rejection method
		sampleElementSource();
	}
	else //face source
	{
		//determine direction and put on correct face
		sampleFaceSource();
	}

	//adjust weight: w(x_i) = p(x_i)/p*(x_i)
	_particle->_weight *= (_curr_el_element_probability + _curr_el_face_probability) / _stratified_probability;

	//increment counter for the current element
	_n_sampled_from_current_element++;
}
