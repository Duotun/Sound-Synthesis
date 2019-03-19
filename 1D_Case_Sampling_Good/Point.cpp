#include "Point.h"

pointmass::pointmass(double m, double x, double y, double z, int ind)
{
	mass = m; dscale = 1;
	position = Eigen::Vector3cd(x, y, z);
	velocity = Eigen::Vector3cd(0.0, 0.0, 0.0);
	gravity = Eigen::Vector3cd(0.0, 0.0, 0.0);  //no gravity now
	index = ind;
	isfixedposition = false;
}
void pointmass::set_position(double x, double y, double z)
{
	position = Eigen::Vector3cd(x, y, z);
}

void pointmass::set_mass(double m)   //add mass
{
	mass += m;  //spring spreading
}

void pointmass::set_fixedposition(bool fix)
{
	isfixedposition = fix;
}
void pointmass::set_distancescale(double scale)
{
	dscale = scale;
}

//connect a spring to a point
void pointmass::add_spring(spring *s)
{
	springs.push_back(s);
	double frequency = sqrt(s->k / mass);
	frequencies.push_back(frequency);    //record frequency
}



//currently take no gravity into consideration
Eigen::Vector3cd pointmass::calculteforce(Eigen::Vector3cd initialforce)
{
	Eigen::Vector3cd fg = gravity * mass;  //currently zero
	Eigen::Vector3cd fs = -initialforce;

	for (int i = 0; i < springs.size(); i++)
	{
		fs += springs[i]->getforce(this);
	}

	Eigen::Vector3cd force = fg + fs;
	return force;
}

void pointmass::update(double dt, Eigen::Vector3cd initialforce)
{
	if (isfixedposition) { return; }
	
	Eigen::Vector3cd accelerate = calculteforce(initialforce) / mass;

	velocity = velocity + accelerate * dt;
	position = position + velocity * dscale*dt;
}
