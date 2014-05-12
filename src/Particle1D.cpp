//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  High Order Solver
//  @ File Name : Particle1D.cpp
//  @ Date : 1/27/2014
//  @ Author : Simon R Bolding
//
//

#include <stdlib.h>
#include "../include/Particle1D.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../include/Controller.h"
#include "../include/Source.h"
#include "../include/LinDiscSource.h"
#include "../include/ResidualSource.h"
#include "../include/StandardResidualSource.h"

Particle1D::Particle1D(HoMesh* mesh, Source* src, RNG* rng, string method_str) :
_source(src) //source may be passed a NULL, to be set later
{
	_rng = rng;
	_mesh = mesh;
	_position_mfp = -9999.; //initialize outside teh domain to check that source samplign works correctly
	_n_elements = _mesh->getNumElems();
	_current_element_ID = 0;	   //particle needs to be somewhere to initialize material properties
	_current_element = _mesh->getElement(_current_element_ID);
	_current_mat_ID = 0;
	updateElementProperties(); //intiialize to the material properties of the 0-th element, will likely change once sampling occurs, but must initialize
	_method = HoMethods::method_map.at(method_str);
	_is_dead = true;

	if (HoController::PARTICLE_BALANCE)
	{
		_n_abs = 0;
		_n_leak = 0;
		_n_scat = 0;
		_n_terminations=0;
	}
}

inline double Particle1D::samplePathLength()
{
	return -1.*log(_rng->rand_num())*_mfp_tot;
}

inline double Particle1D::samplePathLengthMFP()
{
	return -1.*log(_rng->rand_num());
}

void Particle1D::sampleCollision()
{
	if (_is_dead)
	{
		return;
	}
	else if ( (_method == HoMethods::HOLO_ECMC) || (_method == HoMethods::HOLO_STANDARD_MC) ) //then a pure absorber problem, end the history
	{
		terminateHistory();
		if (HoController::PARTICLE_BALANCE)
		{
			_n_abs++;
		}
	}
	else if (_method == HoMethods::STANDARD_MC) //usual MC sample which event, sample new direction if scattering
	{
		//determine if a scattering event
		if (_rng->rand_num() < _scat_ratio)
		{
			if (HoController::PARTICLE_BALANCE)
			{
				_n_scat++;
			}
			double mu_scat = sampleAngleIsotropic();
			//use angle addition to get the new scattered cosine
			_mu = _mu*mu_scat - sqrt(1. - mu_scat*mu_scat)*sqrt(1. - _mu*_mu);
		}
		else //Non scattering event
		{
			//TODO, if there were something besides absorption possible you would put it here
			terminateHistory();
			if (HoController::PARTICLE_BALANCE) //debug stuff
			{
				_n_abs++;
			}
		}
	}
	else	
	{
		std::cerr << "Input an incorrect mode of operation\n";
		exit(1);
	}

}

void Particle1D::resetParticleBalance()
{	//reset debug variables
	if (HoController::PARTICLE_BALANCE)
	{
		_n_abs = 0;
		_n_leak = 0;
		_n_scat = 0;
		_n_terminations = 0;
	}
}

void Particle1D::sampleSourceParticle()
{
	initializeHistory(); //Reset the basic parameters
	_source->sampleSourceParticle();
}

inline double Particle1D::sampleAngleIsotropic()
{
	return _rng->rand_num()*2.0 - 1.;
}

void Particle1D::leaveElement()
{
	//Contribute to tallies, determine which element you are entering
	scoreFaceTally();
	if (_current_element->getDownStreamElement() == NULL) //Leaked out of the problem
	{
		if (HoController::PARTICLE_BALANCE)
		{
			_n_leak++;
		}
		terminateHistory();
	}
	else //move to the new element and update properties
	{
		_current_element = _current_element->getDownStreamElement(); //move to the downstream element
		if (_current_element->hasChildren()) //new cell is more refined, need to find new cell
		{
			_current_element = _current_element->findChildEntered(_mu);
		}
		updateElementProperties();
		if (_mu > 0) //moving to the right
		{
			_position_mfp = 0.0;
		}
		else
		{
			_position_mfp = _element_width_mfp;
		}
	}
}

void Particle1D::runHistory()
{
	//start history
	sampleSourceParticle(); //samples position, direction, and initalizes weight
	if (_is_dead)
	{
		return; //have to check because of stratified sampling
	}

	double path_length_mfp;  

	while (true) //stream the particle until it is absorbed or leaks
	{
		path_length_mfp = samplePathLengthMFP(); //sample a pathlength
		streamToNextEvent(path_length_mfp); //stream to leakage or interaction
		sampleCollision(); //if particle is already dead, this will return so
		if (_is_dead)
		{
			return;
		}
	} //end of outer history while loop
}


//return random numbers for use by source or whoever needs one
inline double Particle1D::getRandNum()
{
	return _rng->rand_num();
}

//for when particle has entered a new cell need to update the properties to the current cell
void Particle1D::updateElementProperties()
{
	_element_width_mfp = _sigma_tot*_current_element->getElementDimensions()[0]; //always need to update the width in general
	_element_angular_width = _current_element->getElementDimensions().back(); //angular dimension is always last
	if (_spatial_element == _current_element->getSpatialElement())
	{
		return; //no need to update, materials etc. are still the same
	}
	else
	{
		//check to see if materials have changed
		_spatial_element = _current_element->getSpatialElement();
		if (_current_mat_ID == _spatial_element->getMaterialID())
		{
			return;
		}
		
		MaterialConstant mat;
		mat = _current_element->getSpatialElement()->getMaterial();
		_current_mat_ID = mat.getID();
		_sigma_tot = mat.getSigmaT();
		_mfp_tot = 1. / _sigma_tot;
		_sigma_scat = mat.getSigmaS();
		_sigma_abs = mat.getSigmaA();
		//Set the material data
		if (_method == HoMethods::HOLO_ECMC || _method == HoMethods::HOLO_STANDARD_MC) //pure absorber problem
		{
			_scat_ratio = 0.0;
		}
		else
		{
			_scat_ratio = _sigma_scat / _sigma_tot;
		}
	}

}

void Particle1D::scoreElementTally(double path_start_mfp, double path_end_mfp)
{
	//Score Element tally, need to convert path_length and volume
	//to cm, rather than mfp
	double path_length_cm = std::abs((path_start_mfp - path_end_mfp)*_mfp_tot/_mu);
	double normalized_position = 0.5*(path_start_mfp + path_end_mfp) / _element_width_mfp; //location of the center of pathlength

	//ScoreECMCTallies, using a normalized direction cosine to ensure positive tallies
	_current_element->incrementTallyScores(_weight, path_length_cm, 
		_mu, normalized_position);
}

inline void Particle1D::terminateHistory()
{
	//TODO may need to do other stuff here
	if (HoController::PARTICLE_BALANCE)
	{
		_n_terminations++;
	}
	_is_dead = true;
}

void Particle1D::printParticleBalance(int n_hist, bool reset_particle_balance)
{
	cout << "---------------------------------------------------\n"
		<< "               Particle balance\n"
		<< "---------------------------------------------------\n\n"
		<< "	 Number Created: " << n_hist << endl
		<< "	Number Absorbed: " << _n_abs << endl
		<< "      Number Leaked: " << _n_leak << endl
		<< "    Number Scatters: " << _n_scat << endl
		<< "  Number Terminated: " << _n_terminations << endl;	

	//optionally reset the particle balances, this is true by default
	if (reset_particle_balance)
	{
		resetParticleBalance();
	}
}

inline void Particle1D::initializeHistory()
{
	_weight = 1.0;
	_is_dead = false;
}

inline void Particle1D::scoreFaceTally()
{
	return;
	//Currently inactive, if you wanted to do standard MC you would need to add these back in and 
	//make sure you only score whe yo0u leave a spatial element most likely


	//The face tally is scored before you have left the cell, so everything is
	//based on the cell you are leaving, not the cell you are entering
	/*
	//determine face based on direction
	int face_id = 0; //Leaving to the left
	if (_mu >= 0.0) //Leaving to the right
	{
		face_id = 1;
	}

	//increment the correct tallies
	int face_index = _spatial_element->getID()+face_id;
	_current_face_tallies[face_index]->incrementScore(_weight, _mu, 1.0); //for one-d, per cm sq
	_flux_face_tallies[face_index]->incrementScore(_weight, _mu, 1.0);
	*/ 
}

void Particle1D::streamToNextEvent(double path_length_mfp)
{
	double displacement_mfp;
	double new_position_mfp;

	//determine horizontal displacement
	displacement_mfp = path_length_mfp*_mu;
	new_position_mfp = displacement_mfp + _position_mfp;

	//determine if the particle has left the current element, or not
	while ((new_position_mfp < 0.) || (new_position_mfp > _element_width_mfp)) //particle has left the cell
	{
		//tally variables, corresponding to path across current cell
		double path_start = _position_mfp;
		double path_end;

		//Determine the number of mean free paths remaining to stream after entering new cell
		if (_mu >= 0.0) //streaming to the right
		{
			displacement_mfp += (_position_mfp - _element_width_mfp);
			path_end = _element_width_mfp; //leaves to the right
		}
		else //streaming to the left
		{
			displacement_mfp += _position_mfp;
			path_end = 0.0; //leaves to the left
		}
		scoreElementTally(path_start, path_end);
		leaveElement(); //move to the next element

		if (_is_dead) //then particle has leaked, do not score anything else
		{
			return;
		}
		new_position_mfp = _position_mfp + displacement_mfp; //determine where the particle would be now
	}

	//Now position is within the current cell
	scoreElementTally(_position_mfp, new_position_mfp);
	_position_mfp = new_position_mfp;
}