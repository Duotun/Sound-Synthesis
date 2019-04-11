#pragma once
#ifndef MODALMODEL_H
#define MODALMODEL_H

#include <assert.h>
#include<stdlib.h>
#include<iostream>
#include<vector>
#include <fstream>
#endif
class modalmodel {
private:
	std::vector<double> eigenmodes;
	std::vector<double> eigenvec;
	std::vector<double> omega_;
	std::vector<double> omega_d;   //damping frequency
	std::vector<double> c_;     //  rayleigh damping
	std::vector<double> freq;
	int nummodes; int nrows; //length of each eig-vector
	double beta, alpha;   // 1e-7 , alpha =0
public:   //0.00000025
	int zero_num;
	modalmodel() :beta(0.00000025), alpha(0) {};   //currently don't consider outer input
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
	int nrows_get() const
	{
		return nrows;
	}
	const std::vector<double> damped_omega()const
	{
		return omega_d;
	}
    
	const std::vector<double> frequency()const
	{
		return freq;
	}

	const std::vector<double> eigenmodes_()const
	{
		return eigenmodes;
	}
	const std::vector<double> eigv()
	{
		return eigenvec;
	}
	const std::vector<double> damping_co()const
	{
		return c_;
	}
	int load_eigenmodes(const char* path);
	void accm_modal_impulse();
};
