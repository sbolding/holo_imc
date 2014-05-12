//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : Mesh.cpp
//  @ Date : 11/1/2013
//  @ Author : SRB
//
//


#include <stdlib.h>
#include "../include/Mesh.h"
#include "../include/Node1D.h"
#include "../include/Element1D.h"
#include "../include/DirichletBC1D.h"
#include <vector>
using std::cout;
using std::endl;

Mesh::Mesh()
{
	//Should never be called
	using namespace std;
	cerr << "The default constructor should never be called for Mesh object" << endl;
	exit(0);
}

Mesh::Mesh(int dim, int number_elements, double width, MaterialConstant* material, double* bc_values): //TODO currently constant material accross the problem
_ext_source_functor(NULL) //this will not always be used
{
	//set fstream to dummy file
	_dim = dim;
	Element1D* new_element;
	Node* new_node;
	_mat_list.push_back(material); //TODO this is hardcoded for a single material

	if (_dim != 1)
		exit(0); //not implemented beyond 1D yet

	if (dim == 1)
	{
		//Initially size variables
		_n_nodes = number_elements + 1;
		_n_edges = _n_nodes;
		_n_elems = number_elements;

		//Create Node objects from points
		double elem_width = width / _n_elems;
		double pnt;
		for (int i = 0; i < _n_nodes; ++i)
		{
			pnt = i*elem_width;
			new_node = new Node1D(i, pnt);
			_nodes.push_back(new_node);
		}

		//Create Connectivity Array, Straight forward for 1D (0: 0 1; 1: 1 2; 2: 2,3;... 
		std::vector<int> temp_vec;
		for (int i = 0; i < _n_elems; ++i)
		{	
			temp_vec.push_back(i);
			temp_vec.push_back(i + 1);			
			_connectivity_array.push_back(temp_vec);
			temp_vec.clear();
		}
				
		//Create elements, this code should work for higher dimensions as well			
		std::vector<Node* > element_nodes;
		for (int i = 0; i < _n_elems; ++i)
		{
			for (size_t j = 0; j < _connectivity_array[i].size(); ++j)
			{
				element_nodes.push_back(_nodes[_connectivity_array[i][j]]);
			}			
			new_element = new Element1D(i, material, element_nodes);
			_elements.push_back(new_element);
			element_nodes.clear(); //reset temp vector
		}

		//Set boundary conditions, currently only correct in 1D
		//*******************************************************************
		if (bc_values != NULL)
		{
			_dirichlet_bcs.push_back(new DirichletBC1D(0, _elements[0], _nodes[0], bc_values[0]));
			_dirichlet_bcs.push_back(new DirichletBC1D(1, _elements[_n_elems - 1], _nodes[_n_nodes - 1], bc_values[1]));
			_n_dirichlet_bc = _dirichlet_bcs.size();
		}
	}	
}

Mesh::Mesh(int dim, std::ifstream file)
{
	using namespace std;
	cerr << "Not implemented yet";
}

void Mesh::print(std::ostream & out) const
{
	using namespace std;
	//Print all nodes
	out << "----------------------------Mesh----------------------------" << endl;
	out << "============================Nodes===========================" << endl;
	for (unsigned int i = 0; i < _nodes.size(); i++) //Print out all nodes
		_nodes[i]->print(out); 	
	out << endl;
	out << "=====================Connectivity Array=====================" << endl;
	out << " (Elem. ID) : (Node1 ID) (Node2 ID) ..." <<  endl;
	for (unsigned int i = 0; i < _connectivity_array.size(); ++i)
	{
		out.width(8);
		out << i << "  :  ";
		out.width(4);
		for (unsigned int j = 0; j < _connectivity_array[i].size(); ++j)
		{
			out.width(4);
			out << _connectivity_array[i][j];
		}
		out << endl;
	}
	out << endl;
	out << "==========================Elements==========================" << endl;
	for (int i = 0; i < _n_elems; i++)
		_elements[i]->print(out);
	out << endl;
	out << "=========================Materials==========================" << endl;
	for (unsigned int i = 0; i < _mat_list.size(); ++i)
		_mat_list[i]->print(out);
	out << endl;
}


//-----------------Read and Set functions -----------------------------------------

int Mesh::getNumElems() const
{
	return _n_elems;
}

int Mesh::getNumNodes() const
{
	return _n_nodes;
}

int Mesh::getNumEdges() const
{
	return _n_edges;
}

std::vector<Element* >* Mesh::getElements(void)
{
	return & _elements; //pointer, be careful not to change these values
}

std::vector<DirichletBC1D* > Mesh::getDirichletBCs(void)
{
	return _dirichlet_bcs; //Copy, not a pointer
}

int Mesh::getNumDirichletBC() const
{
	return _n_dirichlet_bc;
}

void Mesh::setExternalSource(double ext_source_constant) 
{
	//Have the elements do the setting
	std::vector<Element*>::const_iterator it_el;  //element iterator
	it_el = _elements.begin();				//initialize iterator
	int n_elem_nodes = (*it_el)->getNumNodes(); //Assuming one mesh type
	std::vector<double> ext_source_nodal_values;
	ext_source_nodal_values.resize(n_elem_nodes);
	
	for (; it_el != _elements.end(); ++it_el)
	{
		for (int i_node = 0; i_node < n_elem_nodes; ++i_node)
		{
			ext_source_nodal_values[i_node] = ext_source_constant;
		}
		(*it_el)->setExtSourceNodalValues(ext_source_nodal_values);
	}
}

void Mesh::setExternalSource(double factor, GlobalMethods::ProblemType p)
{
	using namespace GlobalMethods;
	if (p == FIXED_SOURCE)
	{
		setExternalSource(factor); //set constant isotropic value over domain
	}
	else if (p == EIGEN_VALUE)
	{
		double k_eff = factor; //passed in factor is keff
		std::vector<double> phi_nodal_values;

		//Get the scalar flux nodal values
		getDiscScalarFluxVector(phi_nodal_values);

		//Loop over elements and set the external source values
		std::vector<Element*>::const_iterator it_el;  //element iterator
		it_el = _elements.begin();				//initialize iterator
		int n_elem_nodes = (*it_el)->getNumNodes(); //Assuming one mesh type
		std::vector<double> ext_source_nodal_values;
		ext_source_nodal_values.resize(n_elem_nodes);
		size_t node_idx=0; //which node you are on

		//Set the nodal values of the fission source to nu_sig_f/k phi
		for (; it_el != _elements.end(); ++it_el)
		{
			for (int i_node = 0; i_node < n_elem_nodes; ++i_node)
			{
				double nu_sig_f = (*it_el)->getMaterial().getNuSigmaF();
				ext_source_nodal_values[i_node] = (nu_sig_f/k_eff)*phi_nodal_values[node_idx];
				node_idx++;
			}
			(*it_el)->setExtSourceNodalValues(ext_source_nodal_values);
			std::fill(ext_source_nodal_values.begin(), ext_source_nodal_values.end(), 0.0); //reset
		}

	}
	else
	{
		std::cerr << "Have not implemented this problem type in Mesh.cpp\n";
		exit(1);
	}
}



void Mesh::setExternalSource(FixedSourceFunctor & q)
{
	//store the functor.  You could set the edge valeus if functor is isotropic,
	//but in general no need to do that because the functor case will be handled
	//correctly.
	_ext_source_functor = &q;
}

void Mesh::setBoundaryConditions(const std::vector<std::vector<double>> & ang_flux_moments)
{
	if (_dim != 1)
	{
		std::cerr << "Boundary conditions in Mesh.cpp not implemented yet for higher dimensions\n";
		exit(1);
	}
	_dirichlet_bcs.clear();   // need to clear incase already set
	_dirichlet_bcs.push_back(new DirichletBC1D(0, _elements[0], _nodes[0], ang_flux_moments[0]));
	_dirichlet_bcs.push_back(new DirichletBC1D(1, _elements[_n_elems - 1], _nodes[_n_nodes - 1], ang_flux_moments[1]));
	_n_dirichlet_bc = _dirichlet_bcs.size();
}

Element* Mesh::getElement(int element_id) const
{
	//Check to make sure array is ordered by id
	Element* element = _elements[element_id];
	if (element->getID() != element_id)
	{
		std::cerr << "The elements were not ordered properly" << std::endl;
	}
	return element;
}

int Mesh::getFaceIndex(int element_id, int face_id) const
{
	//use connectivity array to get the correct face index, used for tallying
	return _connectivity_array[element_id][face_id];
}

void Mesh::printLDScalarFluxValues(std::ostream &out) const
{
	//Have the elements do the printing
	using std::vector;
	std::vector<Element*>::const_iterator it_el;  //element iterator
	it_el = _elements.begin();				//initialize iterator

	out << "---------------------------------------------------\n";
	out << "LD scalar flux values: (coordinate), value" << endl;
	out << "---------------------------------------------------\n\n";

	for (; it_el != _elements.end(); ++it_el)
	{
		(*it_el)->printLDScalarFluxValues(out);
	}
}

int Mesh::getSpatialDimension() const
{
	return _dim;
}

void Mesh::getDiscScalarFluxVector(std::vector<double> & v) const
{
	//Have the elements do the getting
	std::vector<Element*>::const_iterator it_el;  //element iterator
	it_el = _elements.begin();				//initialize iterator
	std::vector<double> flux_edge_values, dummy;
	v.clear();

	for (; it_el != _elements.end(); ++it_el)
	{
		(*it_el)->getScalarFluxHOClosure(flux_edge_values, dummy);
		for (int i = 0; i < flux_edge_values.size(); ++i)
		{
			v.push_back(flux_edge_values[i]);
		}
	}
}

