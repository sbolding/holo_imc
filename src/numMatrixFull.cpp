// numMatrix class definition

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include "../include/Define.h"
#include "../include/numMatrixFull.h"

using std::cout;
using std::cerr;
using std::endl;

// Constructors
numMatrixFull::numMatrixFull(int nr, int nc)
{
	using std::cout;
	_n_rows = nr;
	_n_cols = nc;

    //allocate pointers to rows
	_coeff.resize(_n_rows);
	
	// Allocate columns for each row
	for(int i=0; i<_n_rows; i++) {
		_coeff[i].resize(_n_cols);
	}
	zero();
}

//Destructor
numMatrixFull::~numMatrixFull()
{
}		

//Default constructor should never be called
numMatrixFull::numMatrixFull()
{
	std::cerr << "this should not be getting called, default constructor for numMatrixFull" << std::endl;
	exit(0);
}


// Functions
void numMatrixFull::setCoeff(int i, int j, double value)
{
	using std::cout;
	if((i<0) || (i>_n_rows)){
		cout << "Error in SetCoeff";
		exit(0);
	}
	if((j<0) || (j>_n_cols)){
		cout << "Error in SetCoeff";
		exit(0);
	}
	_coeff[i][j] = value;
}

void numMatrixFull::addCoeff(int i, int j, double value)
{
	using std::cout;
	if((i<0) || (i>_n_rows)){
		cout << "Error in AddCoeff";
		exit(0);
	}
	if((j<0) || (j>_n_cols)){
		cout << "Error in AddCoeff";
		exit(0);
	}
	_coeff[i][j] += value;
}

double numMatrixFull::getCoeff(int i, int j)
{
	using std::cout;
	if((i<0) || (i>_n_rows)){
		cout << "Error in GetCoeff";
		exit(0);
	}
	if((j<0) || (j>_n_cols)){
		cout << "Error in GetCoeff";
		exit(0);
	}
	return _coeff[i][j];
}

void numMatrixFull::mult(numMatrix *b, numMatrix *c)
{
	// Check for consistency.
	if(_n_cols != b->getNumRows()){
		cout << "Unable to multiply matrices in numMatrix.Mult";
		exit(0);
	}

	// Multiply the two matrices.
	for(int i=0; i<_n_rows; i++){
		for(int j=0; j<b->getNumCols(); j++){
			c->setCoeff(i, j, 0.);
			for(int k=0; k<_n_cols; k++){
				c->addCoeff(i, j, _coeff[i][k]*b->getCoeff(k,j));
			}
		}
	}
}

void numMatrixFull::mult(numVector *b, numVector *c)
{
	// Check for consistency.
	if(_n_cols != b->getNumRows()){
		cout << "Unable to multiply matrices in numMatrix.Mult";
		exit(0);
	}

	// Multiply the matrix times the vector.
	for(int i=0; i<_n_rows; i++){
		c->setCoeff(i, 0.);
		for(int k=0; k<_n_cols; k++){
			c->addCoeff(i, _coeff[i][k]*b->getCoeff(k));
		}
	}
}

void numMatrixFull::trans(numMatrix *a)
{
	// Check for consistency.
	if((_n_rows != a->getNumCols()) || 
	   (_n_cols != a->getNumRows())){
		cout << "Unable to transpose matrix in numMatrix.Trans";
		exit(0);
	}

	//transpose the matrix
	for(int i=0; i<a->getNumRows(); i++){
		for(int j=0; j<a->getNumCols(); j++){
			_coeff[j][i] = a->getCoeff(i, j);
		}
	}
}

void numMatrixFull::solve(numVector *b, numVector *x)
{
/****************************************************************************
**                                                                         **
**  This function solves a linear system of equations.  It uses a          **
**  uses a standard Gauss elimination scheme.  A is the coefficient        **
**  matrix, B is the right side vector, and NEQ is the number of           **
**  equations.  On return, data in A has been destroyed and the solution   **
**  is stored in the calling vector.                                       **
****************************************************************************/

	double pivot;
	double factor;
	int size, irow, ipv, icol; 
	size = _n_rows;
	for(irow=0; irow<size; irow++ ){
		for(ipv=0; ipv<size; ipv++){

			// Normalize the pivot row by dividing across by the pivot.
			pivot = _coeff[ipv][ipv];
			for(icol=0; icol<size; icol++){
                _coeff[ipv][icol] /= pivot;
//				a.SetCoeff(ipv, icol, a.GetCoeff(ipv, icol)/pivot);
			}
			b->setCoeff(ipv, b->getCoeff(ipv)/pivot);

            // Replace row(irow) with row(irow)-factor*row(pivot).  
			// This makes the column elements all zero except for the 
			// pivot row. 
			for(irow=0; irow<size; irow++){
				if(irow != ipv){
                    factor = _coeff[irow][ipv];
//					factor = a.GetCoeff(irow, ipv);
					for(icol=ipv; icol<size; icol++){
                        _coeff[irow][icol] -= factor*_coeff[ipv][icol];
//						a.AddCoeff(irow, icol, -factor*a.GetCoeff(ipv, icol));
					}
					b->addCoeff(irow, -factor*b->getCoeff(ipv));
				}
			}
		}
	}
	for(int i=0; i<size; i++){
		x->setCoeff(i, b->getCoeff(i));
	}
}
