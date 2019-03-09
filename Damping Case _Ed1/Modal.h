#pragma once
#ifndef MODAL_H
#define MODAL_H

#include<iostream>
#include<stdlib.h>
#include "Massspringsystem.h"

#endif // !MODAL_H
class modalmodel {
private:
	std::vector<double>eigenmodes;
	std::vector<double> omega_;
	std::vector<double> omega_d;   //damping frequency
	std::vector<double> c_;     //  rayleigh damping
	int nummodes;
	double beta, alpha;   // 1e-7 , alpha =0
public:
	modalmodel():beta(0.00000025),alpha(0){};   //currently don't consider outer input
	const std::vector<double> omega_all()const
	{
		return omega_;
	}
	double omega(int mid) const
	{
		assert(mid >= 0 && mid < nummodes);
		return omega_[mid];
	}
	int num_modes()const
	{
		return nummodes;
	}
	const std::vector<double> damped_omega()const
	{
		return omega_d;
	}

	const std::vector<double> eigenmodes_()const
	{
		return eigenmodes;
	}
	const std::vector<double> damping_co()const
	{
		return c_;
	}
	void load_eigenmodes(std::vector<double> &frequency, massspringsystem &s);
	void accm_modal_impulse();
};