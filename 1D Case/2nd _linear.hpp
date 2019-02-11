#include <iostream>
#include <Eigen\Dense>
#include <stdlib.h>
#include "Massspringsystem.h"
#include <stdio.h>

// solve 1d equation
void solve_2nd_1d(Eigen::MatrixXd &S,Eigen::VectorXd &force, std::vector<double> & amp, massspringsystem &s, int index)
{
	double w_2 = S(index, index); double b = force(index);
	pointmass* m = s.massget(index);
    // 1d position
	double position = m->position(0).real();
	std::cout << position << std::endl;
	//std::cout << position << std::endl;
	//double  amplitude = (position - b /w_2);
	double amplitude = b / w_2;
	std::cout << "amplitude: " <<amplitude << std::endl;
	amp.push_back(amplitude);
}