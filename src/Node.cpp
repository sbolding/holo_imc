#include "../include/Node.h"

Node::Node(int id)
{
	_id = id;
}

double Node::getX() const
{
	return _x;
}

double Node::getY() const
{
	return _y;
}

int Node::getID() const
{
	return _id;
}