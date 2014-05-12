#ifndef _RNG_H
#define _RNG_H

#include<string>
#include<sstream>
#include<vector>
#include<iostream>
#include<iomanip>




using namespace std;

class RNG
{
private:
	//Constants
	int a, M, Q, R;
	double Minv; //Minv is inverse of M

	//Variables
	int seed, k;  //seed is set initially in initialize(), then overwritten with new integervalue

public:
	RNG(); //sets the values of the modulus, multiplier, and seed
	double rand_num(); //returns a random number between 0. and 1.
};

#endif
