#include "Massspringsystem.h"

void massspringsystem::addmass(double mass, double x, double y, double z, int ind)
{
	pointmass *m = new pointmass(mass, x, y, z,ind);
	m->set_distancescale(dscale);
	masses.push_back(m);
	//return m;
}

pointmass* massspringsystem::massget(int index)
{
	return masses[index];
}

spring* massspringsystem::addspring(double young, double restlength, pointmass *mass1, pointmass *mass2)
{
	spring* s = new spring(young, restlength, mass1, mass2);
	s->set_dscale(dscale);
	mass1->add_spring(s);
	mass2->add_spring(s);
	springs.push_back(s);
	return s;
}

void massspringsystem::setdistancescale(double scale)
{
	dscale = scale;

	for (int i = 0; i<masses.size(); i++) {
		masses[i]->set_distancescale(dscale);
	}

	for (int i = 0; i<springs.size(); i++) {
		springs[i]->set_dscale(dscale);
	}
}

void massspringsystem::update(double dt,Eigen::Vector3cd initialforce)
{
	for (int i = 0; i < masses.size(); i++)
	{
		masses[i]->update(dt,initialforce);
	}
}

std::vector<spring*> massspringsystem::getspring()
{
	return springs;
}
std::vector<pointmass*> massspringsystem::getmasses()
{
	return masses;
}