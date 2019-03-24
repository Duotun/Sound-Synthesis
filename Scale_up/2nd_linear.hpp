#include <iostream>
#include <Eigen\Dense>
#include <stdlib.h>
#include "Massspringsystem.h"
#include <stdio.h>
#include <cmath>
// solve 1d equation
void solve_2nd_1d_no_damping(Eigen::MatrixXd &S, double &force, std::vector<double> &frequency, std::vector<double> & amp, massspringsystem &s, int index)
{
	double w_2;
	if (frequency[index] != 0.0)
		w_2 = abs(S(index, index));
	else
	{
		amp.push_back(0.0);   //if frequency is zero
		return;
	}
	double b = force;   //6 zero in 3d, 0 in 1d, u^t
	pointmass* m = s.massget(index);     // 3
										 // 1d position
										 //double position = m->position(0).real();
										 //std::cout << position << std::endl;
										 //std::cout << position << std::endl;
										 //double  amplitude = (position - b /w_2);
	double amplitude = b / m->mass / w_2;
	//std::cout << "amplitude: " <<amplitude << std::endl;
	amp.push_back(amplitude);  //one frequency is added
}

void solve_2nd_1d_damping(Eigen::MatrixXd &S, double &force, std::vector<double> &frequency, std::vector<double> & amp, massspringsystem &s, int index)
{
	double w_2;
	if (frequency[index] != 0.0)
		w_2 = abs(S(index, index));
	else
	{
		amp.push_back(0.0);   //if frequency is zero
		return;
	}
	//double epsilon = s.springget(index)->d / 2 / s.massget(index)->mass / w_2;
	//double expotential=exp(-epsilon*w_2);  //t£¿
	double b = force;
	pointmass* m = s.massget(index);
	// 1d position
	//double position = m->position(0).real();
	//std::cout << position << std::endl;
	//std::cout << position << std::endl;
	//double  amplitude = (position - b /w_2);
	double amplitude = b / m->mass / w_2;
	//std::cout << "amplitude: " <<amplitude << std::endl;
	amp.push_back(amplitude);  //one frequency is added
}

//amplitude_normalize
void amplitude_normalize(std::vector<double> &ampli, double maxamp)
{
	std::cout << "Amplitude of All Synthesized Sounds: " << std::endl;
	for (int i = 0; i < ampli.size(); i++)
	{
		double tmp = abs(ampli[i]) / maxamp;
		tmp = round(tmp * 1000) / 1000.0;    // 0.001 precision
		ampli[i] = tmp;
		if (tmp != 0.0)
			std::cout << tmp << std::endl;
	}
	std::cout << std::endl;
}
//y = kx+b
void linear_relations_for_mode_compression(double x1, double y1, double x2, double y2, double &k, double &b)  //0-2khz
{
	k = (y2 - y1) / (x2 - x1);
	b = y2 - k * x2;
}
