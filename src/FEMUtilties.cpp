#include <stdlib.h>
#include "../include/FEMUtilities.h"

GaussQuadrature::GaussQuadrature():
_weights(), _gauss_points()
{
	_n_points = 2;
	_weights = { 1., 1. };
	_gauss_points = { -0.5773502691896257, 0.5773502691896257 };
}

GaussQuadrature::GaussQuadrature(int n_gauss_points):
_weights(), _gauss_points()
{
	if (n_gauss_points < 2 || n_gauss_points > 4)
	{
		std::cerr << "Only 2 and 3 point quadrature currently available\n";
		exit(1);
	}
	_n_points = n_gauss_points;
	switch (_n_points)
	{
	case(2):
		_weights = { 1., 1. };
		_gauss_points = { -0.5773502691896257, 0.5773502691896257 };
		break;
	case(3):
		_weights = { 0.5555555555555556, 0.8888888888888888, 0.5555555555555556 };
		_gauss_points = { -0.7745966692414834, 0.0L, 0.7745966692414834 };
		break;
	case(4):
		_weights = { 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538 };
		_gauss_points = { -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526 };
	}
}

int GaussQuadrature::getNumPoints()
{
	return _n_points;
}

std::vector<double> GaussQuadrature::getQuadraturePoints(double center, double width)
{
	std::vector<double> points(_n_points);
	if (width < 0.0)
	{
		std::cerr << "Passed in a negative width to quadrature, in FEMutitlites.cpp\n";
		exit(1);
	}
	double a = center - 0.5*width;
	double b = center + 0.5*width;
		
	for (int i = 0; i < _n_points; i++)
	{
		points[i] = a + 0.5*width + _gauss_points[i] * 0.5*width;
	}
	return points;
}

std::vector<double> GaussQuadrature::getQuadratureWeights(double center, double width)
{
	double a = center - 0.5*width;
	double b = center + 0.5*width;

	std::vector<double> wgts(_n_points);
	for (int i = 0; i < _n_points; ++i)
	{
		wgts[i] = _weights[i] * 0.5*width; //multiply by 0.5 because gauss quadrature is over intervale [-1,1]
	}
	return wgts;
}

void FEMUtilities::convertMomentsToEdgeValues1D(const std::vector<double> & moment_dof, std::vector<double> & nodal_values) //convert moment dof to edge values
{
	if (moment_dof.size() != 2)
	{
		std::cerr << "Passed in incorrect length of moment vector to FEMUtilties::convertMomentsToEdgeValues1D\n";
		exit(1);
	}
	nodal_values.resize(2);
	nodal_values[0] = moment_dof[0] - moment_dof[1];
	nodal_values[1] = moment_dof[0] + moment_dof[1];
}

void FEMUtilities::convertAvgSlopeToBasisMoments1D(const std::vector<double> & moment_dof, std::vector<double> & left_right_moments)
{
	if (moment_dof.size() != 2)
	{
		std::cerr << "Passed in incorrect length of moment vector to FEMUtilties::convertAvgSlopeToBasis\n";
		std::cerr << "Length = " << moment_dof.size() << std::endl;
		exit(1);
	}
	left_right_moments.resize(2);
	left_right_moments[0] = moment_dof[0] - moment_dof[1] / 3.; //left basis moment
	left_right_moments[1] = moment_dof[0] + moment_dof[1] / 3.; //right basis moment
}

void FEMUtilities::convertEdgeValuesToAvgSlope1D(const std::vector<double> & nodal_values, std::vector<double> & moment_dof)
{
	if (nodal_values.size() != 2)
	{
		std::cerr << "Passed in incorrect length of moment vector to FEMUtilties\n";
		std::cerr << "Length = " << moment_dof.size() << std::endl;
		exit(1);
	}
	moment_dof.resize(2);
	moment_dof[0] = 0.5*(nodal_values[0]+nodal_values[1]); //left basis moment
	moment_dof[1] = 0.5*(nodal_values[1] - nodal_values[0]); //right basis moment
}

