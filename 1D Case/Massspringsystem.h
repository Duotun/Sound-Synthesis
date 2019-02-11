#pragma once
#ifndef  MASSSPRINGSYSTEM_H
#define  MASSSPRINGSYSTEM_H

#include<vector>
#include "Point.h"
#include "Spring.h"
#include <Eigen\Dense>
#include <iostream>

class massspringsystem
{
public:
	massspringsystem() {
		dscale = 1;   //meters, currently
	};
	~massspringsystem()
	{
		for (int i = 0; i < masses.size(); i++)
		{
			delete masses[i];
		}

		for (int i = 0; i < springs.size(); i++)
		{
			delete springs[i];
		}
	};
	void addmass(double mass, double x, double y, double z,int ind);
	pointmass * massget(int index);
	spring* addspring(double young, double restlength, pointmass *mass1, pointmass *mass2);

	void setdistancescale(double scale);
	void update(double dt, Eigen::Vector3cd initialforce);
	std::vector<spring*> getspring();
	std::vector<pointmass*> getmasses();
private:
	std::vector<pointmass*> masses;
	std::vector<spring*>springs;
	double dscale;    //corresponding scale

};


#endif //! MASSSPRINGSYSTEM_H
