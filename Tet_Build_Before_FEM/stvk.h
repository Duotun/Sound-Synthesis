#pragma once
#include "Tensor.hpp"
#include<iostream>
#include<cmath>

#ifndef  STVK_H
#define  STVK_H

class StVKMaterial
{
public:
	StVKMaterial() {};
	StVKMaterial(double youngs, double poisson) {
		lambda_ = 0.5*youngs / (1 + poisson);
		mu_ = youngs * poisson / (1 + poisson) / (1 - 2 * poisson);
	};
	StVKMaterial(const StVKMaterial &mat) :lambda_(mat.lambda_), mu_(mat.mu_) {}
	StVKMaterial(const StVKMaterial *mat) :lambda_(mat->lambda_), mu_(mat->mu_) {}

	void set_parameter(double l, double m)
	{
		lambda_ = l; mu_ = m;
	}
	void first_piola_kirchhoff(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const;
	void second_piola_kirchhoff(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const;
	void cauchy_stress_tensor(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const;

	//evaluate the stiffness matrix for a single tetrahedron
	void stiffness_matrix() const;

	//strain energy density in undeformed volume yes! to deformed
	double strain_energy_density(const Eigen::MatrixX3d &F) const;
	~StVKMaterial();

private:
	double lambda_;
	double mu_;
};

StVKMaterial::~StVKMaterial()  //currently non definition
{
}

#endif // ! STVK_H
