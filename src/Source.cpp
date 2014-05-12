#include <stdlib.h>
#include "../include/Source.h"
#include "../include/Particle1D.h"
#include "../include/GlobalConstants.h"
#include <iostream>
#include "../include/FEMUtilities.h"

Source::Source()
{
	//do nothing;
}

Source::~Source()
{
	//No dynamic memory
}

Source::Source(Particle1D* particle) :
_particle(particle),
_ext_src_functor(NULL)
{
	_particle = particle;
	_rng = particle->_rng; 
	_vol_src_total = 0.0;
	_BC_src_total = 0.0;

	//make sure particle's pointer is to this source
	_particle->_source = this;
}

Source::Source(Particle1D* particle, FixedSourceFunctor* q) :
_particle(particle),
_ext_src_functor(q)
{
	_rng = particle->_rng;
	_vol_src_total = 0.0;
	_BC_src_total = 0.0;

	//make sure particle's pointer is to this source
	_particle->_source = this;
}

void Source::sampleAngleIsotropic()
{
	_particle->_mu = _rng->rand_num()*2.0 - 1.;
}

void Source::sampleAngleIsotropic(double min_cosine, double max_cosine, bool directional_bias)
{
	//This function works correctly for - or + numbers, as long as min_cosine < max_cosine, e.g., -1 < -0.5
	if (min_cosine > max_cosine)
	{
		std::cerr << "Passed a min cosine greater than max cosine to sampling routine" << std::endl;
		exit(1);
	}
	_particle->_mu = _rng->rand_num()*(max_cosine - min_cosine) + min_cosine;

	//Adjust weight to be unbiased if necessary
	if (!directional_bias)
	{
		_particle->_weight *= (max_cosine - min_cosine)*0.5;
	}
}

void Source::sampleAngleCosineLaw(int n, double min_cosine, double max_cosine)
{
	if (min_cosine > max_cosine || n < 0 || min_cosine*max_cosine < 0.0)
	{
		if (std::abs(min_cosine) > GlobalConstants::RELATIVE_TOLERANCE
			&& std::abs(max_cosine) > GlobalConstants::RELATIVE_TOLERANCE) //special case for roundoff error (0 can be negative)
		{
			std::cerr << "Passed a min cosine greater than max cosine to sampling routine, or cell across zero" << std::endl;
			exit(1);
		}
	}
	//sample from pdf f(mu) = norm_const*mu^n, mu in [min_cosine, max_cosine]
	double norm_const = (std::pow(max_cosine, n + 1) - std::pow(min_cosine, (n + 1)));
	if (min_cosine >= 0.0)
	{
		_particle->_mu = _rng->rand_num()*norm_const + pow(min_cosine, n + 1);
		_particle->_mu = pow(_particle->_mu, 1. / (n + 1.));
	}
	else //mu is negative and odd function so special case
	{
		_particle->_mu = _rng->rand_num()*abs(norm_const) + std::abs(pow(max_cosine, n + 1));
		_particle->_mu = -1.*pow(_particle->_mu, 1. / (n + 1.));
	}

}

void Source::sampleAngleCosineLaw(double min_cosine, double max_cosine)
{
	if (min_cosine > max_cosine || min_cosine*max_cosine < 0.0)
	{
		if (std::abs(min_cosine) > GlobalConstants::RELATIVE_TOLERANCE
			&& std::abs(max_cosine) > GlobalConstants::RELATIVE_TOLERANCE) //special case for roundoff error (0 can be negative)
		{
			std::cerr << "Passed a min cosine greater than max cosine to sampling routine, or cell across zero" << std::endl;
			exit(1);
		}
	}
	
	if (min_cosine < -1.*GlobalConstants::RELATIVE_TOLERANCE) //for negative case f(mu) = nrom_const*|mu|, mu in [min_cosine, max_cosine]
	{
		_particle->_mu = -1.*std::sqrt(_rng->rand_num()*(min_cosine*min_cosine-
			max_cosine*max_cosine) + max_cosine*max_cosine);
	}
	else //sample from pdf f(mu) = norm_const*mu, mu in [min_cosine, max_cosine]
	{
		_particle->_mu = std::sqrt(_rng->rand_num()*(max_cosine*max_cosine -
			min_cosine*min_cosine) + min_cosine*min_cosine);
	}
}

double Source::sampleLinDiscFunc1D(std::vector<double> nodal_values, double left_node_coor, double right_node_coor)
{
	std::cout << "This function has not been checked yet" << std::endl;
	system("pause");
	exit(1); 

	//Sample the position based on the nodal values, should write a function to get the area of the source from the element somehow
	double coordinate;
	double width = (right_node_coor - left_node_coor);
	if (width > 0.0)
	{
		std::cerr << "Passed in coordinates in reverse order to sampleLinDiscFunc, in source.cpp" << std::endl;
		exit(1);
	}

	//If this routine is too slow, do a soft check to see if they are different first, then do the check below
	if (abs(nodal_values[0] - nodal_values[1]) / nodal_values[0] < 1.E-10) //then effectively a constant source, sampling is uniform across the cell
	{
		coordinate = _rng->rand_num()* width + left_node_coor;
	}
	else //need to sample from lin discontinuous source //THIS ROUTINE WORKS
	{
		double left_hat, right_hat; //Normalized nodal values, such that CDF is normalized
		left_hat = 2.0*nodal_values[0] / (nodal_values[1] + nodal_values[0]);
		right_hat = 2.0 - left_hat;
		//use direct inversion of CDF to sample position, based on quadratic formula
		coordinate = -left_hat + sqrt(left_hat*left_hat + 2 * _rng->rand_num()*(right_hat - left_hat));
		coordinate /= (right_hat - left_hat);
		coordinate = coordinate*width + left_node_coor; //convert to mfp
	}
	return coordinate;
}

void Source::mapExtSrcToElement(std::vector<double> & ext_src_ld_dof, double & tot_src_strength,
	Element* spatial_element, ECMCElement1D* element)
{
	//Determine if the source is isotropic.
	bool is_isotropic = false;
	if (_ext_src_functor == NULL) is_isotropic = true;
	if (element->hasChildren())
	{
		std::cerr << "No abstract elements should reach this function, mapExtSrcToElement in Source.cpp\n";
		exit(1);
	}
	
	//initialize the LD moments of external source over the element to zero or scattering source (isotropic)
	ext_src_ld_dof.assign(3, 0.0);
	std::vector<double> spatial_x_coors;
	double x_left_el;
	double x_right_el;
	if (_particle->_method == HoMethods::HOLO_ECMC || _particle->_method == HoMethods::HOLO_STANDARD_MC) //append scattering source
	{
		double sigma_s_el = spatial_element->getMaterial().getSigmaS();
		std::vector<double> scat_src_nodal_values_spat = spatial_element->getScalarFluxNodalValues();//This should return 0 if LO system hasnt been solved yet
		
		//map the ext source strength nodal values on spatial element to nodal values on the current ECMCelement
		spatial_x_coors = spatial_element->getNodalCoordinates();
		x_left_el = element->getSpatialCoordinate() - 0.5*element->getSpatialWidth();
		x_right_el = x_left_el + element->getSpatialWidth();
		std::vector<double> scat_nodal_values_el(2);
		scat_nodal_values_el[0] = evalLinDiscFunc1D(scat_src_nodal_values_spat, spatial_x_coors, x_left_el);
		scat_nodal_values_el[1] = evalLinDiscFunc1D(scat_src_nodal_values_spat, spatial_x_coors, x_right_el);

		std::vector<double> scat_src_avg_slope; //convert edge values to average and slope
		FEMUtilities::convertEdgeValuesToAvgSlope1D(scat_nodal_values_el, scat_src_avg_slope); 
		
		//initialize values to isotropic external source moments (p/(sec str))
		for (int mom = 0; mom < scat_src_avg_slope.size(); mom++)
		{
			ext_src_ld_dof[mom] += scat_src_avg_slope[mom] * sigma_s_el/GlobalConstants::FOUR_PI; //phi*_sigma_s/4pi
		}
	}
	
	//handle isotropic case and functor case seperately
	if (is_isotropic) //use LD angular integrated values.  If eigenvalue problem, q_nodal_vals is actually \nu_sigma_f/k_eff*phi
	{
		std::vector<double> q_nodal_values_spat_el(spatial_element->getExtSourceNodalValues()); //initialize to ext source values, this has units of particles/sec-cm
		std::vector<double> q_moments_int; //integrated over angle

		//map the ext source strength nodal values on spatial element to nodal values on the current ECMCelement
		if (_particle->_method != HoMethods::HOLO_ECMC && _particle->_method != HoMethods::HOLO_STANDARD_MC) //save time
		{
			spatial_x_coors = spatial_element->getNodalCoordinates();
			x_left_el = element->getSpatialCoordinate() - 0.5*element->getSpatialWidth();
			x_right_el = x_left_el + element->getSpatialWidth();
		}
		std::vector<double> q_nodal_values_el(2);
		q_nodal_values_el[0] = evalLinDiscFunc1D(q_nodal_values_spat_el, spatial_x_coors, x_left_el);
		q_nodal_values_el[1] = evalLinDiscFunc1D(q_nodal_values_spat_el, spatial_x_coors, x_right_el);

		//convert ECMC element values
		FEMUtilities::convertEdgeValuesToAvgSlope1D(q_nodal_values_el, q_moments_int);

		//add the isotropic values to the totals
		for (int mom = 0; mom < q_moments_int.size(); mom++)
		{
			ext_src_ld_dof[mom] += q_moments_int[mom] / GlobalConstants::FOUR_PI; // (p/(sec-str))
		}
	}
	else //use q functor
	{
		std::vector<double> q_moments = _ext_src_functor->getHoMoments(
			element->getElementCoordinates(), element->getElementDimensions()); //moments in p/(sec-str)
		for (int mom = 0; mom < q_moments.size(); mom++)
		{
			ext_src_ld_dof[mom] += q_moments[mom];
		}
	}

	//compute total sorce strength \int_mu \int_x q(x,mu) dx dmu = h_x*h_mu*q_avg
	tot_src_strength = ext_src_ld_dof[0]*element->getSpatialWidth()*element->getAngularWidth(); //fraction in this angular element, units of p / (sec)
}

void Source::convertNodalValuesToMoments(std::vector<double> & nodal_values,
	std::vector<double> & ld_moments, bool nodal_values_isotropic)
{
	if (nodal_values.size() != 2)
	{
		std::cerr << "This method is only implemented for 1D, in Source::convertNodalValues\n";
		system("pause");
		exit(1);
	}
	//Set ld_moments to have one more DOF than nodal values
	ld_moments.clear();
	ld_moments.assign(nodal_values.size() + 1, 0.0);

	//compute average and spatial moment
	ld_moments[0] = 0.5*(nodal_values[0] + nodal_values[1]);
	ld_moments[1] = 0.5*(nodal_values[1] - nodal_values[0]);
	if (!nodal_values_isotropic)
	{
		std::cerr << "I have not implemented this method yet" << std::endl;
		system("pause");
		exit(1);
	}
	else //isotropic
	{
		//angular moment is zero because uniform value in angle, the nodal values have taken into account angular fraction already	
	}
}

void Source::convertMomentsToNodalValuesIsotropic(std::vector<double> & ld_moments, std::vector<double> & nodal_values)
{
	if (nodal_values.size() != 2)
	{
		std::cerr << "This method is only implemented for 1D, in Source::convertMomentsToNodalValues\n";
		exit(1);
	}
	//Set nodal_values to have one less DOF than moments
	nodal_values.assign(ld_moments.size() - 1, 0.0);

	//compute left and right edge values 
	nodal_values[0] =2.*(ld_moments[0] - ld_moments[1]);
	nodal_values[1] = 2.*(ld_moments[0] + ld_moments[1]);

	std::cerr << "I have not checked this yet. in source.cpp\n";
	exit(1);
}