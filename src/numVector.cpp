// numVector class definition

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <ostream>
//#include "../include/Define.h"
#include "../include/numMatrix.h"
#include "../include/numVector.h"
using namespace std;

numVector::numVector()
{
	_n_rows = 0;
	_coeff = NULL;
}

numVector::numVector(int nr)
{
	_n_rows = nr;

	// Allocate pointer to rows
	_coeff = new double[_n_rows];
	if(_coeff == NULL) {
		cout << "allocation failure in numVector";
		exit(0);
	}
	zero();
}

//Destructor
numVector::~numVector()
{
	delete [] _coeff;
}		

// Functions
void numVector::zero()
{
	for(int i=0; i<_n_rows; i++){
		_coeff[i] = 0.;
	}
}

int numVector::getNumRows() const
{
	return _n_rows;
}

void numVector::setCoeff(int i, double value)
{
	if((i<0) || (i>_n_rows)){
		cout << "Error in SetCoeff";
		exit(0);
	}
	_coeff[i] = value;
}

void numVector::addCoeff(int i, double value)
{
	if((i<0) || (i>_n_rows)){
		cout << "Error in AddCoeff";
		exit(0);
	}
	_coeff[i] += value;
}

double numVector::getCoeff(int i) const
{
	if((i<0) || (i>_n_rows)){
		cout << "Error in GetCoeff";
		exit(0);
	}
	return _coeff[i];
}

void numVector::print(ostream &out) const
{
	for(int i=0; i<_n_rows; i++){
		out << "Row " << i << "  ";
		out.setf(ios::scientific);
		out.precision(3);
		out.width(12);
		out << _coeff[i] << endl;	
	}
    out << endl;
}

double& numVector::operator[](int i) const
{
	if ((i<0) || (i>_n_rows)){
		cout << "Error in [] vector";
		exit(0);
	}
	return _coeff[i];
}

double&  numVector::operator[](int i) 
{
	if ((i<0) || (i>_n_rows)){
		cout << "Error in [] vector";
		exit(0);
	}
	return _coeff[i];
}

void numVector::assign(double value)
{
	for (unsigned int i=0; i < _n_rows; ++i)
	{
		_coeff[i] = value;
	}
}

numVector::numVector(const numVector& rhs)
{
	_n_rows = rhs.getNumRows();

	// Allocate pointer to rows
	_coeff = new double[_n_rows];
	if (_coeff == NULL) {
		cout << "allocation failure in numVector";
		exit(0);
	}

	//copy members
	for (unsigned int i = 0; i < _n_rows; ++i)
	{
		_coeff[i] = rhs.getCoeff(i);
	}
}

numVector& numVector::operator=(const numVector& rhs)
{
	/*In the case of self assignment, this is necessary to ensure scenarios
	such as (*numVector) = (*numVector)*5. are done correctly. This is exception
	safe since _coeff initialization check done by hand*/
	if (this == &rhs) return *this;

	_n_rows = rhs.getNumRows();

	// Allocate pointer to rows
	delete [] _coeff;
	_coeff = new double[_n_rows];
	if (_coeff == NULL) {
		cout << "allocation failure in numVector";
		exit(0);
	}

	//copy members
	for (unsigned int i = 0; i < _n_rows; ++i)
	{
		_coeff[i] = rhs.getCoeff(i);
	}
	return *this;
}

numVector& numVector::operator*(double scalar)
{
	for (unsigned int i = 0; i < _n_rows; ++i)
	{
		_coeff[i] *= scalar;
	}
	return *this;
}

numVector& numVector::operator*=(double scalar)
{
	for (unsigned int i = 0; i < _n_rows; ++i)
	{
		_coeff[i] *= scalar;
	}
	return *this;
}
