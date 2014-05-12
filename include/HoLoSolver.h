// Created: May 2 2014
// Author: srb
// Project: holo
//
// *************************************************************
// This class is the main solver/interface with the user.  
// Most of the work is done in the high order solver.  This
// class sets most of the global constants throughout the 
// code by accessing static variables.  Input is read from the
// passed in xml file
// *************************************************************

#ifndef _HOLOSOLVER_H
#define _HOLOSOLVER_H

#include "Mesh.h"
#include "MaterialConstant.h"
#include <iostream>
#include "LoSolver1D.h"
#include "HoSolver.h"
#include <fstream>
#include "DataTransfer.h"
#include <cmath>
#include "ConstFixedSource.h"
#include "MMSFixedSource.h"
#include "numMatrixBanded.h"
#include "numVector.h" 
#include "LoSolver1DEigen.h"

class HoLoSolver
{
protected:

	/* ****** Problem Parameters ********** */
	//Main solvers and io
	LoSolver* _lo_solver;
	HoSolver* _ho_solver;
	ofstream _out_file;
	ifstream _input_file;
	MaterialConstant* _mat;
	Mesh* _mesh;
	FixedSourceFunctor* _ext_source; //not always needed

	//Data needed to reconstruct ho solver each time
	unsigned int _n_histories;
	int _n_batches; //maximum number of ECMC batches
	double _exp_convg_rate;
	string _solver_mode; //"standard-mc", "holo-ecmc", "holo-standard-mc"
	string _sampling_method; //"stratified", "standard"

	//Type of problem you are running (eigenvalue, fixed_source, etc.)
	GlobalMethods::ProblemType _prob_type;

	//Problem Data
	int _dimension;
	int _num_elems;
	int _n_ang_elements; //number angles in half ranges
	double _lo_tol;
	size_t _n_holo_solves;  //max number of outer sweeps
	/* ************************************ */

	//For constructor
	void parseInputFile();
public:
	HoLoSolver();
	HoLoSolver(ifstream & input_file);
	~HoLoSolver();
	void solveProblem(); //Run the HoLo solution
};

#endif 
