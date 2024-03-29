//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : holo
//  @ File Name : LoSolver1DEigen.cpp
//  @ Date : 5/4/2014
//  @ Author : srb
//
//


#include <stdlib.h>
#include "../include/LoSolver1DEigen.h"

double LoSolver1DEigen::updateFissionSource()
{
	using GlobalConstants::FOUR_PI;
	using EquationMaps1D::EQN_MAP;

	//Update the scalar flux before updating the fission source, this always must happen in general
	updateSystem();

	//set vector to fission source on right hand side of the equation (no 1/k factor)
	std::vector<Element*>::const_iterator it_el;  //element iterator
	std::vector<Element*>* elements = _mesh->getElements();				//initialize iterator
	it_el = elements->begin();
	int n_elem_nodes = (*it_el)->getNumNodes(); //Assuming one mesh type
	std::vector<int> eqn_nums; //What row in load vector each belongs to
	std::vector<double> temp_vec, phi; //LD scalar flux over element (left and right face values)
	double h; //width of element

	//For computing integral of fission source over domain (no 1/k factor)
	double src_integral = 0.0;

	for (; it_el != elements->end(); ++it_el)
	{
		//get equ nums and scalar flux LD values
		eqn_nums = (*it_el)->getEqnNumbers();
		(*it_el)->getScalarFluxHOClosure(phi, temp_vec);
		h = (*it_el)->getElementDimensions()[0];
		double nu_sig_f = (*it_el)->getMaterial().getNuSigmaF();

		//Compute left and right basis moments
		double fiss_source_left = h / FOUR_PI * nu_sig_f * (2. / 3.* phi[0] + 1. / 3.*phi[1]);
		double fiss_source_right = h / FOUR_PI * nu_sig_f * (1. / 3.* phi[0] + 2. / 3.*phi[1]);

		//compute addition to the integral over the fission source (integrated over space and angle)
		src_integral += h * nu_sig_f * 0.5*(phi[0] + phi[1]);

		//Add them to the appropriate equation numbers based on eqn_map ordering
		_system_vec->setCoeff(eqn_nums[EQN_MAP[0]], fiss_source_left); //<.>_L^+ mom eqn.
		_system_vec->setCoeff(eqn_nums[EQN_MAP[1]], fiss_source_right); //<.>_R^+ mom eqn., etc.
		_system_vec->setCoeff(eqn_nums[EQN_MAP[2]], fiss_source_left); // <.>_L^-
		_system_vec->setCoeff(eqn_nums[EQN_MAP[3]], fiss_source_right);// <.>_R^-
	}
	return src_integral;
}

LoSolver1DEigen::LoSolver1DEigen(Mesh* mesh, std::string solver_mode, double conv_tol) : LoSolver1D(mesh),
_eigen_method_str(solver_mode),
_convg_tol(conv_tol),
_k_eff(1.0), //guess fixed currently
_k_eff_previous(_k_eff),
_first_solve(true),
_fiss_src_integral(0.0),
_fiss_src_integral_prev(0.0),
_fiss_src(_system_vec)
{
	if (solver_mode.compare("arnoldi") == 0)
	{
		_eigen_method = LoMethods::ARNOLDI;
	}
	else if (solver_mode.compare("nka") == 0)
	{
		_eigen_method = LoMethods::NKA;
	}
	else //fixed point
	{
		_eigen_method = LoMethods::FIXED_POINT;
	}

	//Check for vaccuum boundary conditions
	std::vector<DirichletBC1D*> dirichlet_bcs = _mesh->getDirichletBCs();
	std::vector<DirichletBC1D *>::iterator it_bc; //Boundary condition iterator
	it_bc = dirichlet_bcs.begin();
	for (; it_bc != dirichlet_bcs.end(); it_bc++)
	{
		if (abs((*it_bc)->getIncFluxAvg()) != 0.0)
		{
			std::cerr << "Eigenvalue problems currently only work for vaccuum boundary conditions";
			exit(1);
		}
	}
}

LoSolver1DEigen::~LoSolver1DEigen()
{
	//all dynamic memory want to keep
}

void LoSolver1DEigen::fixedPointIteration(double tol, unsigned int n_max_iters)
{
	/* Standard fixed point iteration to 
	determine the eigenvalue of the system:
	L Phi = M/k Phi */

	//Check convergence based on fission source
	numVector* fiss_src_prev = new numVector(*_fiss_src); //fiss src vector has already been initialized
	(*_fiss_src) *= (1. / _k_eff_previous);

	//Perform power iterations until converged
	double error_k_eff;

	for (unsigned int iter = 0; iter < n_max_iters; iter++)
	{
		//Invert the system, i.e., phi^l+1 = L^-1*M/k phi
		solveLinearSystem();

		//Update the fission source and calculate integral over domain
		_fiss_src_integral = updateFissionSource();

		//Calculate a new eigenvalue
		_k_eff = _k_eff_previous* (_fiss_src_integral / _fiss_src_integral_prev);

		//Store the fission source for next iteration
		_fiss_src_integral_prev = _fiss_src_integral;
		(*fiss_src_prev) = (*_fiss_src);

		//divide fission source by k_eff for the next solve
		(*_system_vec) *= (1/_k_eff);
		//_system_vec->print(std::cout);

		//Compute _k_eff error and store for next iteration
		error_k_eff = std::fabs(_k_eff - _k_eff_previous) / (_k_eff);
		_k_eff_previous = _k_eff;

		if (LoController::WRITE_ITERATIONS)
		{
			std::cout << "Iteration " << iter+1 << ", new estimate for k_eff = " << _k_eff
				<< ", relative change of " << error_k_eff << std::endl;
		}

		if (error_k_eff < tol)
		{
			std::cout << "Fixed point iteration for k_eff converged on iteration: " << iter+1
				<< " to error " <<   error_k_eff << std::endl;
			break;
		}
	}
}

void LoSolver1DEigen::arnolidi(double tol, unsigned int n_restart_iter)
{
	std::cerr << "Arnoldi not implemented yet\n";
	exit(1);
}

void LoSolver1DEigen::nKA(double tol)
{
	std::cerr << "NKA not implemented yet\n";
	exit(1);
}

void LoSolver1DEigen::solveSystem()
{
	//Output
	std::cout << "Solving the low order system...\n";

	//Save a copy of the solution vec before assembleSystem() overwrites
	numVector* _sol_copy(NULL);
	if (!_first_solve)
	{
		_sol_copy = new numVector(*_sol_vec);
	}

	//Assemble the standard 1D system, load and sys vector will be zero, but initialized after this
	assembleSystem();

	//If first time solving, make an initial guess of phi = 1 
	if (_first_solve)
	{
		//Initialize scalar flux guess to 1.0 everywhere, and k_eff = k_eff_guess
		_sol_vec->assign(1./(GlobalConstants::FOUR_PI)); 
		_fiss_src_integral_prev = updateFissionSource();
		_first_solve = false;
	}
	else
	{
		//use old solution (stored in _sol_copy) from last time for guess everywhere
		(*_sol_vec) = (*_sol_copy);
		_fiss_src_integral_prev = updateFissionSource();
	}

	//Output to screen if desired
	if (LoController::WRITE_MATRIX)
	{
		std::cout << "\nSystem before solving:\n";
		printSystem(std::cout);
	}

	//Solve the system using desired solver type
	switch (_eigen_method)
	{
	case LoMethods::ARNOLDI:
		arnolidi(_convg_tol, 10);
		break;
	case LoMethods::NKA:
		nKA(_convg_tol);
		break;
	case LoMethods::FIXED_POINT:
		fixedPointIteration(_convg_tol, 1000);
		break;
	}

	//Update the system for use by ho solver (just in case)
	updateSystem();

	//Map the fission source to the external source nodal values for use by the
	//HoSolver
	_mesh->setExternalSource(_k_eff, GlobalMethods::EIGEN_VALUE); //this is an eigenvalue problem

	//Free up memory from matrix and load vector, but keep the solution vector
	deleteMatrixVector();
	delete _sol_copy; //delete copy
}


