#include <stdlib.h>
#include "../include/ResidualSource.h"
#include "../include/Particle1D.h"
#include <cmath>
#include <algorithm>

ResidualSource::ResidualSource(Particle1D* particle, FixedSourceFunctor & q) : Source(particle, &q)
{
	//member variables initialize
	_face_src_total = 0.0;

	//local variables
	std::vector<ECMCElement1D*>* elements;
	std::vector<ECMCElement1D*>::iterator it_el;
	elements = _particle->_mesh->getElements();
	double res_element_mag_el;		         //mag of element residual
	double res_face_mag_el;					//magnitude of residaul on face

	//Local variables to prevent repeated accessing of particle class
	std::vector<double> res_LD_values_el; //the dof of the residual for a particular element
	std::vector<double> res_LD_values_face; //the dof of the residual for the delta source

	//for computing face residuals on boundary
	std::vector<int> boundary_cells = _particle->_mesh->findUpwindBoundaryCells();

	//initialize res vectors
	int n_elems = _particle->_mesh->getNumElems();
	_res_face_mags.resize(n_elems);
	_res_element_mags.resize(n_elems);
	_residual_element_LD_values.resize(n_elems);
	_residual_face_LD_values.resize(n_elems);
	std::vector<double> zeros(2, 0.0);

	for (it_el = elements->begin(); it_el != elements->end(); it_el++)
	{
		if ((*it_el)->hasChildren()) //children are in elements vector, so no need to explicitly add them here
		{
			//TODO, if replaced with a map, no need to add these values here, however,
			//currently zeros are stored for the abstract elements
			_residual_element_LD_values[(*it_el)->getID()] = zeros;
			_residual_face_LD_values[(*it_el)->getID()] = zeros;
			_res_face_mags[(*it_el)->getID()] = 0.0;
			_res_element_mags[(*it_el)->getID()] = 0.0;
			continue;
		}

		//reset magnitudes just in case
		res_face_mag_el = 0.0;
		res_element_mag_el = 0.0;

		//compute volume integrals, for cells that have no children
		computeElementResidual(*it_el, res_LD_values_el, res_element_mag_el);
		_residual_element_LD_values[(*it_el)->getID()] = res_LD_values_el;
		_res_element_mags[(*it_el)->getID()] = res_element_mag_el;

		if ((*it_el)->getDownStreamElement() != NULL) //in this case there is a down wind cell to add delta residual to
		{
			ECMCElement1D* ds_element = (*it_el)->getDownStreamElement();
			if (ds_element->hasChildren()) //ds element more refined
			{
				//need to compute a face residual term for both down stream children
				ECMCElement1D* child_minus = ds_element->findChildEntered(
					ds_element->getAngularCoordinate() - 0.25*ds_element->getAngularWidth()); //get bottom child
				computeFaceResidual(*it_el, child_minus, res_LD_values_face, res_face_mag_el, false);
				_residual_face_LD_values[child_minus->getID()] = res_LD_values_face;
				_res_face_mags[child_minus->getID()] = res_face_mag_el;

				ECMCElement1D* child_plus = ds_element->findChildEntered(
					ds_element->getAngularCoordinate() + 0.25*ds_element->getAngularWidth()); //get top child
				computeFaceResidual(*it_el, child_plus, res_LD_values_face, res_face_mag_el, false);
				_residual_face_LD_values[child_plus->getID()] = res_LD_values_face;
				_res_face_mags[child_plus->getID()] = res_face_mag_el;
			}
			else if (ds_element->getRefinementLevel() == (*it_el)->getRefinementLevel()) //same refinement
			{
				computeFaceResidual(*it_el, ds_element, res_LD_values_face, res_face_mag_el, false); //for non boundary cells
				_residual_face_LD_values[ds_element->getID()] = res_LD_values_face;
				_res_face_mags[ds_element->getID()] = res_face_mag_el;
			}
			else if (ds_element->getRefinementLevel() == (*it_el)->getRefinementLevel() - 1) //ds less refined
			{
				//store LD values for both downstream "ghost" cells.  The sampling routine later must handle the
				//case in which res_ld_values is larger than normal
				int ds_id = ds_element->getID();
				computeFaceResidual(*it_el, ds_element, res_LD_values_face, res_face_mag_el, false);
				if (_residual_face_LD_values[ds_id].size() == 0) //need to initialize
				{
					_res_face_mags[ds_id] = 0.0; //force initialization in case of dictionary
					_residual_face_LD_values[ds_id].resize(res_LD_values_face.size() * 2); //for bottom and top cells
				}
				_res_face_mags[ds_id] += res_face_mag_el; //construct total probability of being in ds cell

				//add the bottom first, then the top (should be an it_el that comes after the bottom one, but forced to add in correct order)
				int begin; //where to begin in vector, 0 for bottom cell, 0+dof.size() for top cell
				begin = ((*it_el)->getAngularCoordinate() < ds_element->getAngularCoordinate() ?
					0 : res_LD_values_face.size());
				for (int i_vec = 0; i_vec < res_LD_values_face.size(); i_vec++)
				{
					_residual_face_LD_values[ds_element->getID()][i_vec + begin] = res_LD_values_face[i_vec];
				}
			}
			else
			{
				std::cerr << "Non conforming mesh, in Residual Source constructor\n";
				exit(1);
			}
		}

		//add magnitudes to appropriate values,
		_vol_src_total += res_element_mag_el; //units of p / sec
		_face_src_total += res_face_mag_el; //units of p / sec
	}

	//compute the boundary condition LD representation over the face of each boundary element
	computeBCAngularFluxDof();

	int bc_element_ID;
	for (int i = 0; i < boundary_cells.size(); i++)
	{
		//add terms to the elements on the boundary
		res_face_mag_el = 0.0; //initialize
		bc_element_ID = boundary_cells[i];
		computeFaceResidual((*elements)[bc_element_ID],(*elements)[bc_element_ID], res_LD_values_face, res_face_mag_el, true); //in this case the d.s. element is the boundary element

		//add values to arrays
		_residual_face_LD_values[bc_element_ID] = res_LD_values_face;
		_res_face_mags[bc_element_ID] = res_face_mag_el;
		_face_src_total += res_face_mag_el; //units of p / sec
	}

	//No BC source, it is implicitly included in the face residual calculations
}

ResidualSource::~ResidualSource()
{
	//No dynamic memory
}

void ResidualSource::sampleFaceSource()
{
	ECMCElement1D* element = _particle->_current_element;
	//put particle on the upwind face
	if (element->getAngularCoordinate() > 0.0) //particle moving to the right
	{
		_particle->_position_mfp = 0.0; 
	}
	else //moving to the left, on right face
	{
		_particle->_position_mfp = _particle->_element_width_mfp;
	}

	std::vector<double> res_dof = _residual_face_LD_values[element->getID()]; //res_x is 0
	std::vector<double> coors = { element->getSpatialCoordinate(), element->getAngularCoordinate() };
	std::vector<double> dimensions = element->getElementDimensions();
	double & h_mu = dimensions[1]; //references to the dimensions vector
	double & dir_coor = coors[1];  //mu_m, direction coordinate

	//Handle special case where element is not as refined as its upstream elements
	if (res_dof.size() > element->getAngularFluxDOF().size()) //then ghost cells
	{
		//figure out which cell you are in
		h_mu *= 0.5;
		std::vector<double> res_dof_bot(res_dof.begin(), res_dof.begin() + 3);
		std::vector<double> res_dof_top(res_dof.begin() + 3, res_dof.end());
		double bot_coor = dir_coor - 0.5*h_mu;
		double top_coor = dir_coor + 0.5*h_mu;
		double mag_bottom = evalDeltaResidualIntegral(res_dof_bot, h_mu, bot_coor);
		double mag_top = evalDeltaResidualIntegral(res_dof_top, h_mu, top_coor);

		if (_rng->rand_num() > mag_bottom / (mag_top+mag_bottom)) //top cell
		{
			dir_coor = top_coor;
			res_dof = res_dof_top;
		}
		else //bottom cell
		{
			dir_coor = bot_coor;
			res_dof = res_dof_bot;
		}
		res_dof.shrink_to_fit();
	}
	
	//Sample angle using rejection method
	double & res_avg = res_dof[0];
	double & res_mu = res_dof[2];
	if (res_dof.size() != 3)
	{
		std::cerr << "Evaluating Lin Disc Function that is not 2D, in ResidualSource.cpp\n";
		system("pause");
		exit(1);
	}
	
	//Determine max value of residual pdf function "f"
	double mu_max = 0.5*dir_coor - 0.25*res_avg*h_mu / res_mu;
	double f_minus = std::abs((dir_coor - 0.5*h_mu)*(res_avg - res_mu)); //value on minus boundary
	double f_plus = std::abs((dir_coor + 0.5*h_mu)*(res_avg + res_mu)); // value on plus boundary
	double f_max = std::max(f_minus, f_plus);
	if (std::abs(mu_max - dir_coor) < 0.5*h_mu) //max of quadratic is within range
	{
		f_max = std::max(f_max, std::abs( mu_max*(evalLinDiscFunc2D(res_dof,dimensions,coors,0.0, mu_max))));
	}

	//Perform rejection sampling for direction
	double mu_new = 0;
	double f_new;
	while (true)
	{
		mu_new = dir_coor + h_mu*(_rng->rand_num() - 0.5); //pick new direction
		f_new = (mu_new*evalLinDiscFunc2D(res_dof, dimensions,coors,0.0,mu_new)); //evaluate residual at new mu TODO abs of mu is artifact, not sure why yet, jakes code has this to force to work, probably because deltas are defined going in the wrong direction
		if (_rng->rand_num()*f_max < std::abs(f_new)) //keep mu
		{
			_particle->_mu = mu_new;
			if (f_new < 0.0) //is residual negative here?
			{
				_particle->_weight *= -1.0;
			}
			break;
		}
	}
}

void ResidualSource::sampleElementSource()
{
	ECMCElement1D* element = _particle->_current_element;
	
	//local variables
	double h_mu = element->getAngularWidth();
	double h_x = element->getSpatialWidth();
	double dir_coor = element->getAngularCoordinate();  //mu_m, direction coordinate, cell center mu
	double pos_coor = element->getSpatialCoordinate();  //x_i, cell center
	std::vector<double> res_dof = _residual_element_LD_values[element->getID()]; //res_x is 0
	double & res_avg = res_dof[0];
	double & res_x = res_dof[1];
	double & res_mu = res_dof[2];

	//determine maximum of function over the element (one of the nodes since it is linear)
	double f_max = std::abs(res_avg) + std::abs(res_x) + std::abs(res_mu);

	//perform rejection sampling
	double f_new, mu_new, x_new; //sampled variables
	std::vector<double> coors = { element->getSpatialCoordinate(), element->getAngularCoordinate() };
	std::vector<double> dimensions = element->getElementDimensions();
	while (true)
	{
		x_new = pos_coor + h_x*(_rng->rand_num() - 0.5); 
		mu_new = dir_coor + h_mu*(_rng->rand_num() - 0.5); 
		f_new = evalLinDiscFunc2D(res_dof, dimensions, coors, x_new, mu_new);

		if ((_rng->rand_num()*f_max) < std::abs(f_new)) //keep sample
		{
			_particle->_mu = mu_new;
			_particle->_position_mfp = ((x_new - pos_coor) / h_x + 0.5)
				* _particle->_element_width_mfp; //normalized position in element
			
			if (f_new < 0.0) //residual was negative set weight
			{
				_particle->_weight *= -1.0;
			}
			break;
		}
	}
}

void ResidualSource::computeElementResidual(ECMCElement1D* element,
	std::vector<double> & res_LD_values_el, double & res_mag)
{
	//initialize variables
	res_mag = 0.0;
	res_LD_values_el.assign(3, 0.0);
	double h_mu = element->getAngularWidth();
	double dir_coor = element->getAngularCoordinate();
	double h_x = element->getSpatialWidth();
	std::vector<double> ang_flux_dof = element->getAngularFluxDOF();
	double & psi_avg = ang_flux_dof[0]; //references for better legibility
	double & psi_x = ang_flux_dof[1];
	double & psi_mu = ang_flux_dof[2];
	double & res_avg = res_LD_values_el[0]; //aliases for computing residual magnitude
	double & res_x = res_LD_values_el[1];
	double & res_mu = res_LD_values_el[2];

	//Local variables to prevent repeated accessing of particle class
	Element* spatial_element; //spatial element corresponding to the current element
	std::vector<double> ext_source_LD_values; //the external source
	double ext_source_mag;
	std::vector<double> total_src_nodal_values_el;   //external source nodal values

	//Map external source (and scattering source if ECMC) onto each element
	spatial_element = element->getSpatialElement();
	double sigma_tot = spatial_element->getMaterial().getSigmaT(); //total cross section
	mapExtSrcToElement(ext_source_LD_values, ext_source_mag, spatial_element, element); //determine external source nodal values over the element

	//compute residual LD values for the element
	res_LD_values_el[0] = ext_source_LD_values[0] - sigma_tot*psi_avg
		- dir_coor*2. / h_x*psi_x;  //average term
	res_LD_values_el[1] = ext_source_LD_values[1] - sigma_tot*psi_x; //x moment
	res_LD_values_el[2] = ext_source_LD_values[2] - sigma_tot*psi_mu
		- psi_x*h_mu / h_x;  //mu moment

	//compute four corner values to figure out which integral to use
	double r_right_plus = res_LD_values_el[0] + res_LD_values_el[1] + res_LD_values_el[2];
	double r_right_minus = res_LD_values_el[0] + res_LD_values_el[1] - res_LD_values_el[2];
	double r_left_plus = res_LD_values_el[0] - res_LD_values_el[1] + res_LD_values_el[2];
	double r_left_minus = res_LD_values_el[0] - res_LD_values_el[1] - res_LD_values_el[2];

	//this first block is for the special cases to eliminate divide by zero errors
	if (std::abs(res_x) <= std::abs(res_avg)*GlobalConstants::RELATIVE_TOLERANCE) //x_slope~0
	{
		if (std::abs(res_mu) <= std::abs(res_avg)*GlobalConstants::RELATIVE_TOLERANCE) //mu_slope~0 also
		{
			res_mag = std::abs(res_avg);
		}
		else //mu slope non-zero
		{	
			if (r_right_minus*r_right_plus > 0.0) //no sign change
			{
				res_mag = std::abs(res_avg);
			}
			else
			{
				res_mag = (res_mu*res_mu + res_avg*res_avg) / std::abs(2.*res_mu);
			}
		}
	}
	else if (std::abs(res_mu) <= std::abs(res_avg)*GlobalConstants::RELATIVE_TOLERANCE) //mu_slope~0
	{
		if (r_left_plus*r_right_plus > 0.0) //no sign change
		{
			res_mag = std::abs(res_avg);
		}
		else //sign change in x
		{
			res_mag = (res_x*res_x + res_avg*res_avg) / std::abs(2.*res_x);
		}
	}
	//---
	else if (r_left_plus*r_right_plus >= 0.0) //no sign change on top
	{
		if (r_left_minus*r_right_minus >= 0.0) //no change on bottom
		{
			if (r_left_plus*r_left_minus >= 0.0) //no change at all
			{
				res_mag = std::abs(res_LD_values_el[0]);
			}
			else if (r_left_plus*r_right_minus < 0.0) //change on left and right
			{
				res_mag = (res_avg*res_avg + res_x*res_x / 3. + res_mu*res_mu) / std::abs(2.*res_mu);
			}
			else
			{
				std::cerr << "You should not be here, error in integral logic, ResidualSource.cpp\n";
				system("pause");
				exit(1);
			}
		}
		else if (r_left_minus*r_right_minus < 0.0) //change on bottom (minus)
		{
			if (r_left_plus*r_right_minus < 0.0) // change on bottom  and right
			{
				res_mag = std::abs(res_avg + pow(r_right_minus, 3) / (12.*res_x*res_mu));
			}
			else if (r_right_plus*r_left_minus < 0.0) //change on bottom and left
			{
				res_mag = std::abs(res_avg + pow(r_left_minus, 3) / (12.*res_x*res_mu));
			}
			else
			{
				std::cerr << "You should not be here, error in integral logic, ResidualSource.cpp\n";
				system("pause");
				exit(1);
			}
		}
		else
		{
			std::cerr << "You should not be here, error in integral logic, ResidualSource.cpp\n";
			system("pause");
			exit(1);
		}
	}
	else if (r_left_plus*r_right_plus < 0.0) //change on top (plus)
	{
		if (r_left_minus*r_right_minus>=0.0) //no change on bottom
		{
			if (r_right_plus*r_right_minus <= 0.0) //change on top and right
			{
				res_mag = std::abs(res_avg + pow(r_right_plus, 3) / (12.*res_x*res_mu));
			}
			else if (r_left_plus*r_left_minus <= 0.0) //change on top and left
			{
				res_mag = std::abs(res_avg + pow(r_left_plus, 3) / (12.*res_x*res_mu));
			}
			else
			{
				std::cerr << "You should not be here, error in integral logic, ResidualSource.cpp\n";
				system("pause");
				exit(1);
			}
		}
		else // change on top and bottom
		{
			res_mag = (res_avg*res_avg + res_x*res_x + res_mu*res_mu / 3.) / (2 * std::abs(res_x));
		}
	}
	else
	{
		std::cerr << "You should not be here, error in integral logic, ResidualSource.cpp\n";
		system("pause");
		exit(1);
	}

	//multiply by area
	res_mag *= h_x*h_mu;
}

void ResidualSource::computeFaceResidual(ECMCElement1D* element, ECMCElement1D* ds_element, std::vector<double> & res_LD_values_face, double & res_mag,
	bool on_boundary)
{
	//This function projects to the finest cell (between element and ds_element) and does computation based on those dimensions

	if (ds_element == NULL) 
	{
		res_LD_values_face.clear();
		res_mag = 0.;
		return;
	}

	//initialize variables
	res_mag = 0.0;
	res_LD_values_face.assign(3, 0.0);
	double h_mu = element->getAngularWidth();
	double dir_coor = element->getAngularCoordinate(); //this may be changed in logic later on
	double & res_avg = res_LD_values_face[0];
	double & res_mu = res_LD_values_face[2]; 
	double mu_sgn = dir_coor / std::abs(dir_coor); //negative or positive direction

	std::vector<double> psi_down; //downstream values
	std::vector<double> psi_up; //upstream element values

	if (on_boundary) //for non isotropic boundary conditions this will not be true
	{
		psi_up = _bc_dof[_bc_element_to_dof_map[element->getID()]];
		psi_down = ds_element->getAngularFluxDOF();
	}
	else
	{	
		if (ds_element->getRefinementLevel() == element->getRefinementLevel()+1) 
		{
			//downstream element is more refined, map dof values to ds elem
			std::vector<double> psi_up_big = element->getAngularFluxDOF();
			dir_coor = ds_element->getAngularCoordinate();
			//need to change dimensions to be on the downstream element scale
			h_mu = ds_element->getAngularWidth();
			double top_or_bottom_sgn = (dir_coor > element->getAngularCoordinate() ? 1. : -1.);

			//map the upstream values to the ds scale
			psi_up.resize(psi_up_big.size());
			psi_up[1] = psi_up_big[1] * 0.5;
			psi_up[2] = psi_up_big[2] * 0.5;
			psi_up[0] = psi_up_big[0] + mu_sgn*psi_up[1] + top_or_bottom_sgn*psi_up[2]; //map to (top_or_bottom) (right: mu>0, left: mu<0);

			psi_down = ds_element->getAngularFluxDOF();//set the ds element
		}
		else if (ds_element->getRefinementLevel() == element->getRefinementLevel() - 1) 
		{
			//upstream element is more refined, map dof values to upstr elem
			std::vector<double> psi_down_big = ds_element->getAngularFluxDOF();
			double top_or_bottom_sgn = (dir_coor > ds_element->getAngularCoordinate() ? 1. : -1.);

			//map the downstream values to the upstream scale
			psi_down.resize(psi_down_big.size());
			psi_down[1] = psi_down_big[1] * 0.5;
			psi_down[2] = psi_down_big[2] * 0.5;
			psi_down[0] = psi_down_big[0] - mu_sgn*psi_down[1] + top_or_bottom_sgn*psi_down[2]; //map to (top_or_bottom) (left: mu>0, right: mu<0);

			psi_up = element->getAngularFluxDOF();
		}
		else if (ds_element->getRefinementLevel() == element->getRefinementLevel()) //same refinement
		{
			psi_up = element->getAngularFluxDOF(); //upstream element values	
			psi_down = ds_element->getAngularFluxDOF(); //downstream values
			if (ds_element != element->getDownStreamElement())
			{
				std::cerr << "The ds element passed in is not actually the down stream element, in ResidualSource::computeFaceResidual\n";
				exit(1);
			}
		}
		else
		{
			std::cerr << "Refinement is not regularized, one cell is more than 1 level refined than an adjacent one, in ResidualSource::computeFaceResidual\n";
			exit(1);
		}
	}

	//mu_sgn*psi[1] tells you if on left or right face
	//res_X is zero here, note that the residual here is -1.*dpsi/dx (NEGATIVE), and does not include mu product, that is added when the residual is sampeld
	res_avg = mu_sgn*((psi_up[0] + mu_sgn*psi_up[1]) - (psi_down[0] - mu_sgn*psi_down[1]));
	res_mu = mu_sgn*(psi_up[2] - psi_down[2]);

    //compute magnitude of integral
	res_mag = evalDeltaResidualIntegral(res_LD_values_face, h_mu, dir_coor);
}

double ResidualSource::evalDeltaResidualIntegral(std::vector<double> res_LD_face_values, double h_mu, double dir_coor)
{
	if (res_LD_face_values.size() != 3 || res_LD_face_values[1] > GlobalConstants::RELATIVE_TOLERANCE)
	{
		std::cerr << "Not using evalDeltaResidualIntegral properily, in ResidualSource.cpp\n";
		exit(1);
	}

	double & res_mu = res_LD_face_values[2];
	double & res_avg = res_LD_face_values[0];
	double res_mag;

	//compute magnitude of integral
	if (std::abs(res_mu / res_avg) < GlobalConstants::RELATIVE_TOLERANCE)
	{
		res_mag = h_mu*std::abs(res_avg*dir_coor);
	}
	if (std::abs(res_avg / res_mu) < 1) //sign change
	{
		double ratio = res_avg / res_mu; //see Jake Peterson A&M thesis for terms, computation has been checked
		res_mag = h_mu*std::abs(dir_coor*0.5*res_mu - h_mu*ratio*ratio*res_avg / 12. + dir_coor*0.5*ratio*ratio*res_mu
			+ 0.25*h_mu*res_avg);
	}
	else //no sign change
	{
		res_mag = h_mu*std::abs(dir_coor*res_avg + h_mu*res_mu / 6.);
	}
	return res_mag;
}

void ResidualSource::computeBCAngularFluxDof()
{
	//Calculate the magnitude of the boundary condition source over elements
	std::vector<DirichletBC1D*> bcs = _particle->_mesh->getDirichletBCs();
	std::vector<int> bc_elem_IDs = _particle->_mesh->findUpwindBoundaryCells();
	std::vector<double> bc_dof_el; //LD dof for each bc, for now isotropic so will only be the average term
	ECMCElement1D* element;
	double direction, bc_avg, bc_mu, psi_avg, psi_mu; //"bc" dof are for an equivalent half range sized element, "psi" dof are for the actual element of interest
	if (bcs.size() > 2)
	{
		std::cout << "Currently only implemented BCs for 1D problems, in ResidualSource.cpp\n";
	}

	for (int i_bc = 0; i_bc < bcs.size(); i_bc++)
	{
		bc_avg = bcs[i_bc]->getIncFluxAvg(); //average flux over half range
		bc_mu = bcs[i_bc]->getIncFluxAngMoment(); //psi_mu over half range
		
		//map flux to boundary elements
		for (int id = 0; id < bc_elem_IDs.size(); id++)
		{
			//get info about element
			int el_id = bc_elem_IDs[id]; //id of current boundary element
			element = _particle->_mesh->getElement(el_id);
			double h_mu = element->getAngularWidth();
			double mu_m = element->getAngularCoordinate(); //incident fluxes are defined over each appropriate half range - or + mu

			//determine the value of mu on the boundary
			double mu_bc = 0.5;
			if (mu_m < 0.0) mu_bc *= -1.;

			if (element->hasChildren())
			{
				std::cerr << "Error in computeBC in ResidualSource, only active elements should be passed from findUpwindBoundaryCells\n";
				exit(1);
			}
			if (element->getSpatialElement() != bcs[i_bc]->getElement()) //on the wrong face
			{
				continue;
			}
			else
			{
				//Need to compute the average and mu moment over the projected ECMC element
				bc_dof_el.assign(3, 0.0);
				bc_dof_el[0] = bc_avg + 2 * bc_mu*(mu_m - mu_bc);
				bc_dof_el[2] = bc_mu*h_mu; //simply scaled by h_mu/h_bc = h_mu/(1.0)
				_bc_dof.push_back(bc_dof_el);
				_bc_element_to_dof_map[el_id] = _bc_dof.size() - 1; //index in bc_dof the last element corresponds to
			}
		}
	}
}

double ResidualSource::getTotalSourceStrength()
{
	return _face_src_total + _vol_src_total;
}



