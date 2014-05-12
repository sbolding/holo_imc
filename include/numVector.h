/* 
 numVector

 This class implements a vector.
*/

#ifndef numVector_h
#define numVector_h

#include <stdio.h>
#include <stdlib.h>
#include <ostream>

class numVector{
  protected:
    double *_coeff;		// Coefficients of matrix
    int _n_rows;       		// Number of rows

  public:
	// Constructors
	numVector();
	numVector(int nr);
	//Shallow copy assignment (creates new memory, doesnt copy pointers)
	numVector(const numVector &vector);
	numVector& operator=(const numVector&);
	numVector& operator*(double scalar); //overload multiply operator
	numVector& operator*=(double scalar);

    //Destructor
    virtual ~numVector();
    // Functions
	void zero();
	void assign(double val); //set vector to constant value
	int getNumRows() const;
	void setCoeff(int i, double value);
	void addCoeff(int i, double value);
	double getCoeff(int i) const;
	double& operator[] (int i) const;
	double& operator[] (int i);
	void print(std::ostream &out) const;
};

#endif
