#include <iostream>
#include <Eigen\Dense>
#include <stdlib.h>
#include "Massspringsystem.h"
#include <stdio.h>

// solve 1d equation
void solve_2nd_1d(Eigen::MatrixXd &S,double &force, std::vector<double> &frequency, std::vector<double> & amp, massspringsystem &s, int index)
{
	double w_2;
	if (frequency[index] != 0.0)
		w_2 = abs(S(index, index));
	else
	{
		amp.push_back(0.0);   //if frequency is zero
		return;
	}
	double b = force;
	pointmass* m = s.massget(index);
    // 1d position
	//double position = m->position(0).real();
	//std::cout << position << std::endl;
	//std::cout << position << std::endl;
	//double  amplitude = (position - b /w_2);
	double amplitude = b / w_2;
	//std::cout << "amplitude: " <<amplitude << std::endl;
	amp.push_back(amplitude);  //one frequency is added
}