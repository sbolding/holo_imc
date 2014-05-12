/* 
 numMatrixFull

 This class stores all the coefficients of the matrix.
*/

#ifndef numMatrixFull_h
#define numMatrixFull_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "numMatrix.h"

class numMatrixFull: public numMatrix{
  protected:
    numMatrixFull(const numMatrixFull *matrix);
  public:
	// Constructors
	numMatrixFull();
    numMatrixFull(int nr, int nc);
    //Destructor
    virtual ~numMatrixFull();
    // Functions
	virtual void setCoeff(int i, int j, double value);
	virtual void addCoeff(int i, int j, double value);
	virtual double getCoeff(int i, int j);
	virtual void mult(numMatrix *b, numMatrix *c);
	virtual void mult(numVector *b, numVector *c);
    virtual void trans(numMatrix *a);
	virtual void solve(numVector *b, numVector *x);
};

#endif
