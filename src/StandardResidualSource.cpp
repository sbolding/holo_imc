#include "../include/StandardResidualSource.h"
#include "../include/Particle1D.h" //Because a friend class, MUST include this here and not in the header file
#include <cmath>

StandardResidualSource::StandardResidualSource(Particle1D* particle, FixedSourceFunctor & q) :
ResidualSource(particle, q)
{
	//Create sampler with alias sampling, let it normalize, delete unneccessary data
	_face_source = new AliasSampler(_res_face_mags, false);
	_element_source = new AliasSampler(_res_element_mags, false);

	//Get rid of the res_face and element magnitudes, they will never be needed again
	_res_face_mags.clear();
	_res_element_mags.clear();
	
	std::vector<double> temp_empty_vector1;
	std::vector<double> temp_empty_vector2;

	//make vectors swap memory with an empty vector that will free after this function call, avoids using c++11 function vector::shrink_to_fit
	_res_face_mags.swap(temp_empty_vector1);
	_res_element_mags.swap(temp_empty_vector2);
}

StandardResidualSource::~StandardResidualSource()
{
	delete _element_source;
	delete _face_source;
}

void StandardResidualSource::sampleSourceParticle()
{
	//Sample if face or element source
	if (_rng->rand_num() < (_vol_src_total / (_vol_src_total + _face_src_total))) //element source
	{
		_particle->_current_element_ID = _element_source->sampleBin(_rng->rand_num(), _rng->rand_num()); //sample bin location
		_particle->_current_element = _particle->_mesh->getElement(_particle->_current_element_ID); //update bin
		_particle->updateElementProperties();

		//determine location and direction within element using rejection method
		sampleElementSource();
	}
	else //face source
	{
		_particle->_current_element_ID = _face_source->sampleBin(_rng->rand_num(), _rng->rand_num()); //sample bin location
		_particle->_current_element = _particle->_mesh->getElement(_particle->_current_element_ID); //update bin
		_particle->updateElementProperties();

		//determine location and direction within element using rejection method
		sampleFaceSource();
	}
}
