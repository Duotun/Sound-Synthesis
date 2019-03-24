#pragma once
#ifndef POINTMASS_H
#define POINTMASS_H

#include<stdlib.h>
#include<iostream>
#include <vector>
#include<Eigen\Dense>
#include "Spring.h"
//frequency is going with point, 
//initial fs =-f_ext differs, currently fg is zero

class spring;
class pointmass {

public:
	pointmass(double mass, double x, double y, double z, int ind);

	void set_mass(double m);
	void set_position(double x, double y, double z);
	void set_fixedposition(bool fix);
	void set_distancescale(double scale);
	void add_spring(spring *s);

	Eigen::Vector3cd calculteforce(Eigen::Vector3cd initialforce);
	void update(double dt, Eigen::Vector3cd initialforce);

	double mass;
	double dscale;
	int index;
	Eigen::Vector3cd position;
	Eigen::Vector3cd velocity;
	Eigen::Vector3cd gravity;
	bool isfixedposition;   //mass fixed or not
	std::vector<spring*> springs; //all the springs 
	std::vector<double> frequencies; // the frequencies corresponding to the springs

};

#endif // !1
