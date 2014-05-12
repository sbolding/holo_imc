#include <stdlib.h>
#include "../include/Source.h"
#include "../include/LinDiscSource.h"
#include "../include/Particle1D.h"
#include "../include/Controller.h"

LinDiscSource::LinDiscSource(Particle1D* particle) : Source(particle)
{
	//need to initialize alias sampler
	//get the area of the source, and the total external source nodal values

	//local variables
	std::vector<ECMCElement1D*>* elements;
	elements = _particle->_mesh->getElements();
	std::vector<double> source_strength_per_cell;
	std::vector<ECMCElement1D*>::iterator it_el;          //element iterator
	it_el = elements->begin();			            //initialize iterator
	double ext_source_el;		             		//magnitude of external source of curren ECMC element, units of p/sec
	std::vector<double> total_src_nodal_values_el, ext_source_LD_values;
	bool only_one_spatial_element = true;

	Element* spatial_element; //spatial element corresponding to the current element
	Element* first_spat_element;
	first_spat_element = (*elements)[0]->getSpatialElement();

	for (; it_el != elements->end(); it_el++)
	{
		if ((*it_el)->hasChildren()) //This source is not implemented for cases with refinedmesh
		{
			std::cerr << "Standard LinDiscSource is not capable of operating on refined meshes, just use residual source\n";
			exit(1);
		}
		try
		{
			//Map external source (and scattering source if ECMC) onto each element
			spatial_element = (*it_el)->getSpatialElement();
			if (first_spat_element != spatial_element)
			{
				only_one_spatial_element = false;
			}
			mapExtSrcToElement(ext_source_LD_values, ext_source_el, spatial_element, *it_el); //determine external source LD values over each element (p/(s-str))
			convertMomentsToNodalValuesIsotropic(ext_source_LD_values, total_src_nodal_values_el); //determine integrated external source edge values (p/s)

			//Add element values to total array
			_total_src_nodal_values.push_back(total_src_nodal_values_el); //have units of particles per steradian, not exactly right, but since normalized doesnt matter
			source_strength_per_cell.push_back(ext_source_el);
			_vol_src_total += ext_source_el; //units of p / sec
		}
		catch (...)
		{
			std::cerr << "The HoSolver had trouble initialized because the Lo System was not properly initialized and not solved, in Particle1D" << std::endl;
			exit(1);
		}
	}

	//Create sampler with alias sampling, let it normalize, delete unneccessary data
	_alias_sampler = new AliasSampler(source_strength_per_cell, false);

	//Create boundary condition source info
	std::vector<DirichletBC1D*> bcs = _particle->_mesh->getDirichletBCs();
	std::vector<int> bc_elem_IDs = _particle->_mesh->findUpwindBoundaryCells();
	std::vector<double> bc_magnitudes; //magnitude of bc on each element
	ECMCElement1D* element;
	double incident_current, two_pi_incident_flux, direction;
	if (bcs.size() > 2)
	{
		std::cout << "Currently only implemented BCs for 1D problems, in LinDiscSource.cpp\n";
	}

	for (int i_bc=0; i_bc < bcs.size(); i_bc++)
	{
		incident_current =  bcs[i_bc]->getCurrent();	//p/sec entering the region
		_BC_src_total += incident_current;
		two_pi_incident_flux = 2.*incident_current; //p/sec-cm^2 This assumes incident flux is isotropic in halfspace

		//map flux to boundary elements
		for (int id=0; id < bc_elem_IDs.size(); id++)
		{
			int el_id = bc_elem_IDs[id]; //id of current boundary element
			element = _particle->_mesh->getElement(el_id);
			if (element->getSpatialElement() != bcs[i_bc]->getElement()) //on the wrong face
			{
				continue;
			}
			else if (only_one_spatial_element)
			{
				std::cout << "This function only works with more than one spatial element, in LinDiscSource.  Use Residual source instead\n";
				exit(1);
			}
			else
			{
				bc_magnitudes.push_back(incident_current* //fraction in this bin
					2.0*element->getAngularWidth()*abs(element->getAngularCoordinate())); //based on cosine law PDF
				_bc_element_map[bc_magnitudes.size() - 1] = el_id; //what element does the last bc value correspond to
			}
		}
	}
	_bc_element_sampler = new AliasSampler(bc_magnitudes);

	if (_vol_src_total + _BC_src_total < GlobalConstants::RELATIVE_TOLERANCE)
	{
		std::cerr << "Source is zero, need to specify boundary or volumetric source\n";
		system("pause");
		exit(1);
	}
}

LinDiscSource::~LinDiscSource()
{
	delete _alias_sampler;
}

void LinDiscSource::sampleSourceParticle()
{
	//Determine if it is volumetric source, or surface source (depending on the mode you are in, may sample scattering source as well)
	//Store the entire source (ext + scattering) into the other one and compute its area.  With the area you can easily determine if sample
	//is from isotropic source or if it is from boundary
	if (_rng->rand_num() < _vol_src_total / (_BC_src_total + _vol_src_total))
	{
		_particle->_current_element_ID = _alias_sampler->sampleBin(_rng->rand_num(), _rng->rand_num()); //sample bin location
		_particle->_current_element = _particle->_mesh->getElement(_particle->_current_element_ID); //update bin
		_particle->updateElementProperties(); 
		sampleLinDiscSource(_total_src_nodal_values[_particle->_current_element_ID]); //sample position in bin
		
		//determine direction within the element
		double mu_center = _particle->_current_element->getAngularCoordinate();
		double half_angular_width = 0.5*_particle->_current_element->getAngularWidth();

		sampleAngleIsotropic(mu_center - half_angular_width, mu_center + half_angular_width); //sample isotropically within the bin
		
	}
	else //Boundary Source
	{
		int alias_element = _bc_element_sampler->sampleBin(_rng->rand_num(), _rng->rand_num());
		_particle->_current_element_ID = _bc_element_map[alias_element];
		_particle->_current_element = _particle->_mesh->getElement(_particle->_current_element_ID); //update bin
		_particle->updateElementProperties();

		//sample direction
		double mu_center = _particle->_current_element->getAngularCoordinate();
		double half_angular_width = 0.5*_particle->_current_element->getAngularWidth();
		sampleAngleCosineLaw(mu_center - half_angular_width, mu_center + half_angular_width);

		//put particle on the correct face
		if (_particle->_mu < 0.0)
		{
			_particle->_position_mfp = _particle->_element_width_mfp;
		}
		else //moving to the right
		{
			_particle->_position_mfp = 0.0;
		}
	}

}

void LinDiscSource::sampleLinDiscSource(std::vector<double> nodal_values)
{
	//Sample the position based on the nodal values, should write a function to get the area of the source from the element somehow

	//If this routine is too slow, do a soft check to see if they are different first, then do the check below
	if (abs(nodal_values[0] - nodal_values[1]) / nodal_values[0] < 1.E-10) //then effectively a constant source, sampling is uniform across the cell
	{
		_particle->_position_mfp = _rng->rand_num()*_particle->_element_width_mfp;
	}
	else //need to sample from lin discontinuous source //THIS ROUTINE WORKS
	{
		double left_hat, right_hat; //Normalized nodal values, such that CDF is normalized
		left_hat = 2.0*nodal_values[0] / (nodal_values[1] + nodal_values[0]);
		right_hat = 2.0 - left_hat;
		//use direct inversion of CDF to sample position, based on quadratic formula
		_particle->_position_mfp = -left_hat + sqrt(left_hat*left_hat + 2 * _rng->rand_num()*(right_hat - left_hat));
		_particle->_position_mfp /= (right_hat - left_hat);
		_particle->_position_mfp *= _particle->_element_width_mfp; //convert to mfp
	}
}

double LinDiscSource::getTotalSourceStrength()
{
	return _BC_src_total + _vol_src_total;
}