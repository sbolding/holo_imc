#include "../include/DirichletBC1D.h"

//Constructor
DirichletBC1D::DirichletBC1D(int id, Element* element, Node* node, double val)
{
	_id = id;
	_element = element;
	_node = node;
	_value_current = val;
	_inc_flux_avg = 2.0*_value_current;
	_inc_flux_mu = 0.0L;
}

DirichletBC1D::DirichletBC1D(int id, Element* element, Node* node, std::vector<double> inc_flux_moments) :
_id(id),
_element(element),
_node(node)
{
	//set incident flux moments
	_inc_flux_avg = inc_flux_moments[0];
	_inc_flux_mu = inc_flux_moments[1];

	if (node->getID() == 0) //left boundary condition
	{
		//Calculate the half range current based on incident (\psi_inc = \psi_a + 2(mu+-0.5)\psi_mu)
		_value_current = 0.5*(_inc_flux_avg + _inc_flux_mu / 3.);
	}
	else //right boundary, inflow is opposite
	{
		_value_current = 0.5*(_inc_flux_avg - _inc_flux_mu / 3.);
	}
}

int DirichletBC1D::getElementID() const
{
	return _element->getID();
}

Element* DirichletBC1D::getElement() const
{
	return _element;
}

Node* DirichletBC1D::getNode() const
{
	return _node;
}

double DirichletBC1D::getCurrent() const 
{
	return _value_current;
}

int DirichletBC1D::getID() const
{
	return _id;
}

int DirichletBC1D::getNodeID() const
{
	return _node->getID();
}

std::vector<double> DirichletBC1D::getAngularFluxDOF() const
{
	std::vector<double> dof;
	dof.push_back(_inc_flux_avg);
	dof.push_back(_inc_flux_mu);
	return dof;
}

double DirichletBC1D::getIncFluxAvg() const
{
	return _inc_flux_avg;
}

double DirichletBC1D::getIncFluxAngMoment() const
{
	return _inc_flux_mu;
}

