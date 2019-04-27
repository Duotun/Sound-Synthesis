#include "Massspringsystem.h"

void massspringsystem::addmass(double x, double y, double z, int ind)
{        // default 1.0f
	pointmass *m = new pointmass(0.0f, x, y, z,ind);
	m->set_distancescale(dscale);
	masses.push_back(m);
	//return m;
}

void massspringsystem::U_T_Compressed(Eigen::MatrixXd	&U_T, std::vector<double> &frequency)  // compress according to the frequencies
{

}
spring* massspringsystem::springget(int index)
{
	if(index>=0&&index<springs.size())
	return springs[index];
}

//vec
std::vector<double> massspringsystem::masses_3d()
{
	const std::vector<pointmass*> mass = this->getmasses();
	std::vector<double> dd;
	int cnt = 0;
	for (size_t i = 0; i < mass.size() * 3; i++)   //corresponding to the mesh size
	{
		if (i != 0 && i % 3 == 0) cnt++;
		dd.push_back(mass[cnt]->mass);    //originally built diagonal
	}
	return dd;

}
pointmass* massspringsystem::massget(int index)
{
	if (index >= 0 && index < masses.size())
	return masses[index];
}

spring* massspringsystem::addspring(double young, double area, pointmass *mass1, pointmass *mass2)
{
	spring* s = new spring(young, area, mass1, mass2);
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

void massspringsystem::reset_mass(double p, double area, pointmass *mass1, pointmass * mass2)
{
	double restlen = (mass1->position.real() - mass2->position.real()).norm();
	double mass = p * area*restlen/2;   //equally distribution
	mass1->set_mass(mass);    //add
	mass2->set_mass(mass);

}
std::vector<spring*> massspringsystem::getspring()
{
	return springs;
}
std::vector<pointmass*> massspringsystem::getmasses()
{
	return masses;
}