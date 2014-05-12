// numMatrix class definition

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include "../include/Define.h"
#include "../include/numMatrixBanded.h"
#include <cmath>
#include <algorithm>



// Constructors
numMatrixBanded::numMatrixBanded(int nr, int band_width) :
_is_LU_decomped(false)
{
	using std::cout;
	_n_rows = nr;
	_n_cols = nr; //assuming square matrix
	_band_width = band_width;
	if (_band_width % 2 != 1)
	{
		std::cout << "Only works for odd banded matrices because of indexing\n";
		exit(1);
	}
	
	//create a sparse matrix
	_coeff.resize(_n_rows);
	
	// Allocate band for each row, need extra 50% of columns for pivoting.  
	// The pivot rows are just stored at the end during solving, everything
	// else is oblivious
	for(int i=0; i<_n_rows; i++) {
		_coeff[i].resize(3 * ((_band_width - 1) / 2) + 1);
		//_coeff[i].resize(_band_width);
	}
}

//Destructor
numMatrixBanded::~numMatrixBanded()
{
}		

//Default constructor should never be called
numMatrixBanded::numMatrixBanded()
{
	std::cerr << "this should not be getting called, default constructor for numMatrixBanded" << std::endl;
	exit(2222);
}


// Functions
void numMatrixBanded::setCoeff(int i, int j, double value)
{
	using std::cout;
	using std::cerr;
	using std::endl;

	if ((i<0) || (i>_n_rows)){
		cout << "Error in GetCoeff";
		exit(0);
	}
	//find col based on band size
	int col = getBandIdx(i, j);
	_coeff[i][col] = value;
}

void numMatrixBanded::addCoeff(int i, int j, double value)
{
	using std::cout;
	using std::cerr;
	using std::endl;

	if ((i<0) || (i>_n_rows)){
		cout << "Error in GetCoeff";
		exit(0);
	}
	j = getBandIdx(i, j);
	_coeff[i][j] += value;
}

double numMatrixBanded::getCoeff(int i, int j)
{
	using std::cout;
	if((i<0) || (i>_n_rows)){
		cout << "Error in GetCoeff";
		exit(0);
	}
	j = getBandIdx(i, j);
	return _coeff[i][j];
}

void numMatrixBanded::mult(numMatrix *b, numMatrix *c)
{
	using std::cout;
	using std::cerr;
	using std::endl;
	std::cout << "not implemented yet in numMatrixBanded.cpp\n";
	exit(1);
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

void numMatrixBanded::mult(numVector *b, numVector *c)
{
	using std::cout;
	using std::cerr;
	using std::endl;
	// Check for consistency.
	if(_n_cols != b->getNumRows()){
		cout << "Unable to multiply matrices in numMatrix.Mult";
		exit(0);
	}

	// Multiply the matrix times the vector.
	for(int i=0; i<_n_rows; i++){
		c->setCoeff(i, 0.);
		for(int k=0; k<_n_cols; k++){
			//check to make sure not a zero in full matrix
			if ( std::abs(k-i) <= (_band_width/2) )
			{
				c->addCoeff(i, _coeff[i][getBandIdx(i,k)] * b->getCoeff(k)); //multiplies by appropriate banded member
			}
			else
			{
				continue;
			}
		}
	}
}

void numMatrixBanded::trans(numMatrix *a)
{
	std::cout << "not implemented yet in numMatrixBanded.cpp\n";
	exit(1);
	// Check for consistency.
	if((_n_rows != a->getNumCols()) || 
	   (_n_cols != a->getNumRows())){
		std::cerr << "Unable to transpose matrix in numMatrix.Trans";
		exit(0);
	}

	//transpose the matrix
	for(int i=0; i<a->getNumRows(); i++){
		for(int j=0; j<a->getNumCols(); j++){
			_coeff[j][i] = a->getCoeff(i, j);
		}
	}
}

void numMatrixBanded::solve(numVector *b, numVector *x)
{

/****************************************************************************
**                                                                         **
**  Solves a banded system by an LU decomposition with partial pivoting
**  then inverting the system.  x is teh solution, b is the source vector.
**  although it is not const, b will not be changed here.
**
**  This will only do inversion first time, after that it will no longer invert
****************************************************************************/
	if (!_is_LU_decomped) //only invert if necessary
	{
		LUDecomposition();
		_is_LU_decomped = true; //save inverted system.  There are no risks with this
		//since inversion destroys original matrix anyways
	}
	invertLUDecomposition(*b, *x);
}

inline int numMatrixBanded::getBandIdx(int row, int col) const
{
	int idx;
	if (row >= 0)
	{
		idx = (col - row) + int(_band_width/2);
	}
	else
	{
		std::cerr << "Cannot use a negative in numMatrixBanded" << std::endl;
	}
	if ( (idx >= _band_width) || (idx < 0) )
	{
		std::cerr << "Accesssed member that is outside of the band" << std::endl;
		exit(1);
	}
	return idx;
}

inline int numMatrixBanded::getColIdx(int row, int band_idx) const
{
	int col;
	if (row >= 0)
	{
		col = band_idx + row - int(_band_width / 2);
	}
	else
	{
		std::cerr << "Cannot use a negative in numMatrixBanded" << std::endl;
		exit(1);
	}
	return col;	
}


void numMatrixBanded::LUDecomposition()
{
	//Column locations (uses integer math)
	int iBand = _band_width;
	int diag = (iBand - 1) / 2;
	int ColEnd = 3 * ((iBand - 1) / 2) + 1;

	//Create pivot vector. References rows.
	//If pivotVec[i] == -1, no pivoting is used.
	//Else, pivoting is used.
	int n = _coeff.size();
	_pivot_vec.resize(n, -1);

	//Clear Extra Matrix Elements
	for (int iRow = 0; iRow <n; iRow++) {
		for (int iCol = iBand; iCol <ColEnd; iCol++) {
			_coeff[iRow][iCol] = 0.0;
		}
	}

	/*! Loop for elimination */
	for (int iRow = 0; iRow < (n - 1); iRow++) {
		/*! Get number of column elements to be eliminated*/
		int iNumElem = std::min(diag, (n - 1) - iRow);

		/*! Use pivoting algorithm.
		*  Obtain the entry in column iRow which has highest absolute value.
		*  This algorithm is necessary because row operation involves divison by pivot.
		*  Computationally roundoff errors are reduced if the denominator is large
		*  Here obtain the index of maximum value for row interchange. */
		int iPivot = -1;
		double iMaxCol = fabs(_coeff[iRow][diag]);
		double iTempColVal;

		for (int iColOffset = 0; iColOffset < iNumElem; iColOffset++) {
			iTempColVal = fabs(_coeff[iRow + (iColOffset + 1)][diag - (iColOffset + 1)]);

			if (iTempColVal > iMaxCol) {
				iPivot = (iColOffset);
				iMaxCol = iTempColVal;
			}
		}

		/*! Put the interchange index in pivot row for updating b[matrix]
		*  in future */
		_pivot_vec[iRow] = iPivot;
		double dTemp1;

		if (iPivot != -1) {
			/*! Interchange the rows for row reduction */
			for (int iColOffset = diag; iColOffset < ColEnd; iColOffset++) {
				dTemp1 = _coeff[iRow][iColOffset];
				_coeff[iRow][iColOffset] = _coeff[iRow + (iPivot + 1)][iColOffset - (iPivot + 1)];
				_coeff[iRow + (iPivot + 1)][iColOffset - (iPivot + 1)] = dTemp1;
			}
		}

		/*! Calculate pivots for row reduction*/
		for (int iColOffset = 0; iColOffset < iNumElem; iColOffset++) {
			_coeff[iRow][diag - (iColOffset + 1)] =
				_coeff[iRow + (iColOffset + 1)][diag - (iColOffset + 1)] / _coeff[iRow][diag];
		}

		/*! Based on the pivots, calculate the new row elements*/
		for (int iRowOffset = 0; iRowOffset < iNumElem; iRowOffset++) {
			for (int iColOffset = (diag + 1); iColOffset < ColEnd; iColOffset++) {
				_coeff[iRow + (iRowOffset + 1)][iColOffset - (iRowOffset + 1)] =
					_coeff[iRow + (iRowOffset + 1)][iColOffset - (iRowOffset + 1)] -
					(_coeff[iRow][diag - (iRowOffset + 1)] * _coeff[iRow][iColOffset]);
			}
		}
	}

	/*! Decomposition of matrix is completed.*/
}

/**
Solves the (n-band) diagonal matrix using pivoting.
Based on Akansha's Math_Utils::solveLUDecomposition,
which was based on Dr. J Morel's Fortran code.
**/
void numMatrixBanded::invertLUDecomposition(const numVector &b, numVector &c)
{
	//copy elements of b to c
	c = b;

	//Column locations (uses integer math)
	int iBand = _band_width;
	int diag = (iBand - 1) / 2;
	int ColEnd = 3 * ((iBand - 1) / 2) + 1;

	int n = c.getNumRows();

	//Assumes LU decomposition of the matrix has been completed.

	/*! Forward Elimination */
	for (int iRow = 0; iRow < (n - 1); iRow++) {

		/*! Interchange row elements in right hand side constants vector*/
		int iRowOffset = (_pivot_vec[iRow]+1);

		if (iRowOffset != -1) {
			double dTemp = c[iRow];
			c[iRow] = c[iRow + iRowOffset];
			c[iRow + iRowOffset] = dTemp;
		}

		/*! Number of Elements in elimination step*/
		int iNumColElemElim = std::min(diag, ((n - 1) - iRow));

		/*! Perform Elimination*/
		for (int iColOffset = 0; iColOffset <iNumColElemElim; iColOffset++) {
			c[iRow + (iColOffset + 1)] = c[iRow + (iColOffset + 1)] -
				(_coeff[iRow][diag - (iColOffset + 1)] * c[iRow]);
		}


	}

	/*! Perform Back substitution */
	/*! Substitur Last element */
	c[n - 1] = c[n - 1] / _coeff[n - 1][diag];

	/*! Back substitute all the elements*/
	for (int iRow = (n - 2); iRow >= 0; iRow--) {
		/*! Calculate the number of elements in back elimination step*/
		int iNumElem = std::min(iBand - 1, ((n - 1) - iRow));
		/*! Perform back substitution*/
		double sum = 0.0;

		for (int iRowOffset = 0; iRowOffset < iNumElem; iRowOffset++) {
			sum += _coeff[iRow][diag + (iRowOffset + 1)] * c[iRow + (iRowOffset + 1)];
		}

		c[iRow] = (c[iRow] - sum) /
			_coeff[iRow][diag];
	}

}

void numMatrixBanded::print(std::ostream & out){
	using std::ostream;
	using std::endl;
	using std::cout;
	using std::ios;

	for (int i = 0; i<_n_rows; i++){
		out << "Row " << i << endl;
		for (int j = 0; j<_n_cols; j++){
			out.setf(ios::scientific);
			out.precision(4);
			out.width(13);
			if (std::abs(j - i) > (_band_width / 2))
			{
				out << 0.0;
			}
			else
			{
				out << _coeff[i][getBandIdx(i, j)];
			}
			//if ((j + 1) % 6 == 0){
			//	out << endl;
			
		}
		out << endl;
		out << endl;
	}
}
