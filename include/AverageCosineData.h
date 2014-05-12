//Contains structures for handling cosine data;  Create structure for mu's averaged on faces (at points in 1D) or volumes

#ifndef _AVERAGEDCOSINEDATA_
#define _AVERAGEDCOSINEDATA_

struct AveragedCosines
{
	double _mu_left_plus;   //Average value of mu on left face, in plus direction
	double _mu_right_plus;  //Average value of mu on right face, in minus direction
	double _mu_right_minus; //Averages in other direction
	double _mu_left_minus;
};

#endif