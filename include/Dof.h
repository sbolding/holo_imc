/* 
 Dof

 This class implements the degree of freedom.  It stores a value, an
 equation number (where it fits in the global matrix), and a flag whether
 it is active (to be solved).  The   assumption is that if the node is not active,
 then the value is an  essential boundary condition.  For the simple 1D case, all
 nodes should always be active because the DOFs are not on nodes, but BC are on
 nodes.  BC must be handled in a different class.

*/

#ifndef Dof_h
#define Dof_h

#include <iostream>

enum DOFType {T=0};

class Dof
{
  public:
	// Constructors
	Dof();
    //Destructor
    virtual ~Dof();
    // Functions
	void setNotActive();
	bool isActive() const;
	int setEqnNum(int currEqnNum);
	int getEqnNum() const;
	void setValue(double value);
	double getValue() const;
	void print(std::ostream &out) const;

  private:
	double _value;
	int _eqn_number;
	bool _active;

	//never used copy constructor
	Dof (Dof *dof);
};

#endif
