// Dof class member functions
 
#include <stdlib.h>
#include "../include/Dof.h"
#include <stdlib.h>

using namespace std;

// Constructors
Dof::Dof(){
	_value = 0.;
	_eqn_number = 0;
	_active = true;
}

// Destructors
Dof::~Dof()
{
}

// Functions
void Dof::setNotActive()
{
	_active = false;
}

bool Dof::isActive() const
{
	return _active;
}

int Dof::setEqnNum(int currEqnNum)
{
	if(_active)
    {
		_eqn_number = currEqnNum;
		currEqnNum++;
	}
	else{
		std::cout << "Error, this case is not implemented yet" << std::endl;
		exit(0);
	}
	return currEqnNum;
}

int Dof::getEqnNum() const
{
	return _eqn_number;
}

void Dof::setValue(double value)
{
	_value = value;
}

double Dof::getValue() const
{
	return _value;
}

void Dof::print(ostream &out) const{
	out << " EqnNum = ";
	out.width(6);
	out << _eqn_number;
	if (_active)
	{
		out << " Active, ";
	}
	else 
	{
		out << " Not Active, ";
	}
	out << " Value = ";
	out.setf(ios::scientific);
	out.precision(3);
	out.width(11);
	out << _value;
}
