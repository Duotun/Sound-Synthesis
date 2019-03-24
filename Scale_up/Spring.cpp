#include "Spring.h"
#include <iostream>
#include <math.h>

spring::spring(double young, double area, pointmass *mass1, pointmass *mass2)
{
	restlen = (mass1->position.real() - mass2->position.real()).norm();
	k = young * area / restlen;
	m1 = mass1; m2 = mass2; dscale = 1; E = 0.0;
}

double spring::get_restlen()
{
	return restlen;
}
void spring::set_e_hessian(double E)
{
	this->E = E;
}
void spring::set_hooker(double young)
{
	k = young * restlen;    //1D Case
}
void spring::set_restlen(double restlength)
{
	restlen = restlength;
}

// unit measurement can be changed 
void spring::set_dscale(double scale)
{
	dscale = scale;
}
Eigen::Vector3cd spring::getforce(pointmass *pulse)   //need to be forced on the two ends
{
	if (pulse != m1 && pulse != m2)
	{
		return Eigen::Vector3cd(0, 0, 0);
	}

	double dis = (m1->position - m2->position).norm()*1e-2;
	double forcescale = dscale * k*(dis - restlen);

	Eigen::Vector3cd dir = (m1->position - m2->position) / dis;

	//don't calculate damping currently     d*velocity

	if (pulse == m1)
	{
		return forcescale * dir;
	}
	else
	{
		return -forcescale * dir;   //opposite toward two ends
	}

}
