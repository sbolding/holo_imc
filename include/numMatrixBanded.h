/* 
 numMatrixBanded

 This class currently just called Morel's sparse matrix class.
*/

#ifndef numMatrixBanded_h
#define numMatrixBanded_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "numMatrix.h"

class numMatrixBanded: public numMatrix{
  protected:
	int _band_width;
	bool _is_LU_decomped; //For solve() method, if it is already inverted, don't want to redo the decomposition
	int getBandIdx(int row, int col) const; //based on a row, column in full matrix, return index in the band for that row
	int getColIdx(int row, int band_idx) const; //based on a row and location in band, return column in the full matrix
    numMatrixBanded(const numMatrixBanded *matrix);

	//for solving the system
	void LUDecomposition(void); //decompose the matrix using LU decomposition
	void invertLUDecomposition(const numVector &b, numVector &c); //b is the right hand side, c is the solution of system Ac=b
	std::vector<int> _pivot_vec;

  public:
	// Constructors
	numMatrixBanded();
    numMatrixBanded(int nr, int band_width); //construct based on the number of rows and how wide of a band (number of non-zero entries across)
    //Destructor
    virtual ~numMatrixBanded();
    // Functions
	virtual void setCoeff(int i, int j, double value);
	virtual void addCoeff(int i, int j, double value);
	virtual double getCoeff(int i, int j);
	virtual void mult(numMatrix *b, numMatrix *c);
	virtual void mult(numVector *b, numVector *c);
    virtual void trans(numMatrix *a);
	virtual void solve(numVector *b, numVector *x);
	virtual void print(std::ostream & out);
};

#endif
