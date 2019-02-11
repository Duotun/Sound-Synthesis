#pragma once
#ifndef  SPRING_H
#define SPRING_H
#include "Point.h"

class pointmass;

class spring
{

public:
	spring(double young, double restlength, pointmass *mass1, pointmass *mass2);
	void set_hooker(double young);
	void set_restlen(double restlength);
	void set_dscale(double scale);
	Eigen::Vector3cd getforce(pointmass *pulse);
	double k; //spring constant
	pointmass *m1;
	pointmass *m2;    //two end, make it public?
private:
	double restlen;   
	double dscale;  
	
};

#endif // ! SPRING_H

