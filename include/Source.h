#ifndef _SOURCE_H
#define _SOURCE_H

#include <vector>
#include "RNG.h"
#include "Element.h"
#include "Controller.h"
#include "GlobalConstants.h"
#include "ECMCElement1D.h"
#include "FixedSourceFunctor.h"

class Particle1D; //forward declaration

class Source
{
protected:

	Particle1D* _particle; //pointer to the particle class, source is a friend of particle class
	Source(); //Never use the default constructor
	RNG* _rng;	//pointer to particles random number generator
	double _vol_src_total; //The total volumetric source (p/sec)
	double _BC_src_total; //The total BC source (p/sec)
	FixedSourceFunctor* _ext_src_functor; //This is only used for computations in anisotropic sources

	//Angle sampling methods
	void sampleAngleIsotropic();
	void sampleAngleIsotropic(double min_cosine, double max_cosine, bool directional_bias=true);  //if false, will adjust weight to produce unbiased particle equivalent to source distribution uniform over sphere, i.e., 1/4PI 
	void sampleAngleCosineLaw(int order, double min_cosine, double max_cosine); //mu^(order) sampling law
	void sampleAngleCosineLaw(double min_cosine, double max_cosine); //mu^1 sampling law, more efficient implementation

	//Spatial sampling methods
	double sampleLinDiscFunc1D(std::vector<double>, double left_node_coor, double right_node_coor); // method samples a coordinate between left and right coors, from linear fn.

	//Useful tools for use with LD sources
	void mapExtSrcToElement(std::vector<double> & ext_src_ld_dof, double & tot_src_strength,
		Element* spatial_element, ECMCElement1D* element); //for computing the ext source ld moments over an element (p/str-cm^3);
	void convertNodalValuesToMoments(std::vector<double> & nodal_values, std::vector<double> 
		& ld_moments, bool nodal_values_isotropic = false); //convert LD nodal values to LD avg and moments, nodal_values has units of p/str-cm^3
	void convertMomentsToNodalValuesIsotropic(std::vector<double> & ld_moments, std::vector<double> & nodal_values); //converts moments for isotropic fluxes to nodal values
	void initializeSamplingSource(); //will create the total source vector, as well as initilize sampling routines, is a unique function because source can be reset between cycles

public:
	Source(Particle1D* particle);
	Source(Particle1D* particle, FixedSourceFunctor *q); //if desired can pass in a functor that is needed for mapping to elements
	virtual ~Source();

	virtual void sampleSourceParticle() = 0; //samples a source particle direction and location
	virtual double getTotalSourceStrength() = 0; //return the total source strenght, virtual because this will be uniqe for residual sources
	double getAreaLinDiscFunction(std::vector<double> nodal_values, double element_volume) //for computing total source magnitudes
	{
		return 0.5*(nodal_values[0] + nodal_values[1])*element_volume;
	};
	double evalLinDiscFunc2D(std::vector<double> dof, std::vector<double> dimensions, 
		std::vector<double>coors , double x, double mu)
	{
		return dof[0] + 2. / dimensions[0]*dof[1] * (x - coors[0])
			+ 2. / dimensions[1]*dof[2] * (mu - coors[1]);
	}
	double evalLinDiscFunc1D(std::vector<double> nodal_values, std::vector<double> x_nodes, double x)
	{
		double width = x_nodes[1] - x_nodes[0];
		if (width < 0.0)
		{
			std::cerr << "Error in evalLinDiscFunc1D, in Source.h\n";
			exit(1);
		}
		return nodal_values[0] * (x_nodes[1] - x) / width + nodal_values[1] * (x - x_nodes[0]) / width;
	}

};








#endif