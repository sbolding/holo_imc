/*
	Class contains all of the information needed to handle a dirichlet BC.

	The boundary condition knows which element and node it belongs to.  The element
	has to figure out how to handle the boundary condition when it constructs the 
	local matrices.  This may be done seperately from the orginial matrix construction
	and just overwrite the values it has set.

*/

#ifndef _DIRICHLETBC1D_H
#define _DIRICHLETBC1D_H

#include "Element.h"
#include "Node.h"

class DirichletBC1D
{
private:

	//TODO you dont really need to know the node for 1D, but it will help in 2D
	Node *  _node;			//The node for which the BC is specified on
	Element * _element;		//The element ''						  "
	double _value_current;	//Value of boundary condition based on incident current
	double _inc_flux_avg, _inc_flux_mu; //angular flux avg and moment, constructors must compute these
	int _id;			    //BC id

	//Never used constructor and copiers
	DirichletBC1D();
	DirichletBC1D operator=(const DirichletBC1D &);
	DirichletBC1D(const DirichletBC1D &);

public:

	DirichletBC1D(int id, Element* element, Node* node, double incid_current); //Construction based on specified incoming current, assumes that incident FLUX is isotropic
	/*construct based on specified angular moments.  Useful for MMS solutions and anisotropic solutions.
	The angular flux dof are exact over their respective half ranges, but the current is defined to always
	be positive. Essentially, the BC angular flux DOF are treated as any other cell, and then current
	is defined as always positive*/
	DirichletBC1D(int id, Element* element, Node* node, std::vector<double> inc_flux_moments); 
	//Public access functions
	int getID() const;
	double getCurrent() const;
	std::vector<double> getAngularFluxDOF() const;
	double getIncFluxAvg() const;
	double getIncFluxAngMoment() const;
	int getElementID() const;
	Element* getElement() const;
	Node* getNode() const;
	int getNodeID() const;

};

#endif