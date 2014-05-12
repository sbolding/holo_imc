#include <stdlib.h>
#include "../include/FixedSourceFunctor.h"
#include "../include/GlobalConstants.h"

FixedSourceFunctor::FixedSourceFunctor()
{
	//currently nothing to do
}

FixedSourceFunctor::~FixedSourceFunctor()
{
	//currently nothing to do
}

std::vector<double> FixedSourceFunctor::getHoMoments(const std::vector<double> & coors, const std::vector<double> & dimens) const
{
	//local variables
	GaussQuadrature quad(4);
	int n_qps = quad.getNumPoints();
	size_t dof_size = dimens.size() + 1;
	std::vector<double> x_pnts(n_qps), x_wgts(n_qps), 
		mu_pnts(n_qps), mu_wgts(n_qps);

	//Check length
	if ( (coors.size() < 1) || (coors.size() > 2) )
	{
		std::cerr << "Passed in too small of vector to FixedSourceFunctor::getHoMoments, also only implemented for 1D\n";
		exit(1);
	}

	//get quadrature information based on dimensions
	x_pnts = quad.getQuadraturePoints(coors[0], dimens[0]);
	x_wgts = quad.getQuadratureWeights(coors[0], dimens[0]);
	mu_pnts = quad.getQuadraturePoints(coors[1], dimens[1]);
	mu_wgts = quad.getQuadratureWeights(coors[1], dimens[1]);

	std::vector<double> ho_moments(dof_size,0.0); //output variable
	std::vector<double> qp_coors(coors.size()); //coordinates of the qp
	double Q_qp;
	for (int i_qp = 0; i_qp < n_qps; ++i_qp) //x_qps
	{
		qp_coors[0] = x_pnts[i_qp]; //update x coordinate 

		for (int j_qp = 0; j_qp < n_qps; ++j_qp) //mu_qps
		{
			qp_coors[1] = mu_pnts[j_qp]; //update mu coordinate

			//1/(h_x*h_mu)*q(x,mu)), where hx and hmu are for the half range element
			Q_qp = x_wgts[i_qp] * mu_wgts[j_qp] / (dimens[0] * dimens[1])* getValue(qp_coors); //compute contribution to cell average

			//compute the integrals for each basis function
			ho_moments[0] += Q_qp;
			ho_moments[1] += 6.0*Q_qp*(x_pnts[i_qp] - coors[0]) / dimens[0]; //these two lines could be replaced by a for loop if x,mu pnts in one vector
			ho_moments[2] += 6.0*Q_qp*(mu_pnts[j_qp] - coors[1]) / dimens[1]; 
		}
	}

	return ho_moments;
}

std::vector<double> FixedSourceFunctor::getLoNodalValues(const std::vector<double> & spatial_coors, const std::vector<double> & spatial_dimens,
	double angular_coordinate) const
{
	//Return Q_L and Q_R for 1D (edge values, not moments).  Easiest way to do it is by calling ho moments and adjusting
	//Determine if on positive or negative halfrange element
	double mu_center;
	if (angular_coordinate < 0.0)
	{
		mu_center = -0.5;
	}
	else
	{
		mu_center = 0.5;
	}

	std::vector<double> coors(spatial_coors.begin(), spatial_coors.end());
	coors.push_back(mu_center);  //Add Half range coordinate
	std::vector<double> dimens(spatial_dimens.begin(), spatial_dimens.end());
	dimens.push_back(1.0); //Add Half range angular width
	std::vector<double> ho_moments(getHoMoments(coors, dimens));

	//compute LO edge values.  The average in angle is equivalent to integral over 0,1 d\mu
	std::vector<double> lo_nodal_vals(2);
	lo_nodal_vals[0] = (ho_moments[0] - ho_moments[1]); //Q_L, 
	lo_nodal_vals[1] = (ho_moments[0] + ho_moments[1]); //Q_R
	return lo_nodal_vals;
}

