//	File contains quadrature and other useful functions for computations on
//  2D LD elements.
//
//  @ File Name : FEMUtilies.h
//  @ Date : 3/18/2013
//  @ Author : Simon R Bolding
//
//

#ifndef _FEMUTILITIES_H
#define _FEMUTILITIES_H

#include <vector>
#include <iostream>

class GaussQuadrature
{
private:

	//Static variables, values set in constructor, currently only 2D
	std::vector<double> _weights; //weights corresponding to each gauss point
	std::vector<double> _gauss_points; //1D gauss points on reference element [-1,1]
	int _n_points;

public:

	GaussQuadrature(); //default 2 gauss points
	GaussQuadrature(int n_gauss_points); //currently only works for 2 and 3
	std::vector<double> getQuadraturePoints(double center, double width); //returns gauss points over the interval a to b
	std::vector<double> getQuadratureWeights(double center, double width); //note, for a higherDimension domain weights have to be multiplied together
	int getNumPoints();
};


namespace FEMUtilities
{
	inline double evalLinDiscFunc2D(std::vector<double> dof, std::vector<double> dimensions,
		std::vector<double>coors, double x, double mu) //evaluate a 2D LD function at x and mu
	{
		return dof[0] + 2. / dimensions[0] * dof[1] * (x - coors[0])
			+ 2. / dimensions[1] * dof[2] * (mu - coors[1]);
	} 
	inline double evalBilinearFunc2D(std::vector<double> dof, std::vector<double> dimensions,
		std::vector<double> coors, double x, double mu) //evaluate a 2D LD bilinear function at x and mu
	{
		return dof[0] + 2. / dimensions[0] * dof[1] * (x - coors[0])
			+ 2. / dimensions[1] * dof[2] * (mu - coors[1]) 
			+ 4 / (dimensions[0] * dimensions[1])*(x - coors[0])*(mu - coors[1]);
	}
	void convertMomentsToEdgeValues1D(const std::vector<double> & moment_dof, std::vector<double> & nodal_values); //convert moment dof to edge values
	void convertAvgSlopeToBasisMoments1D(const std::vector<double> & moment_dof, std::vector<double> & left_right_moments); //convert from average and slope moments to left and right basis function moments
	void convertEdgeValuesToAvgSlope1D(const std::vector<double> & nodal_values, std::vector<double> & moment_dof); //convert edge values to average slope moments
}








#endif // _FEMUTILITIES_H