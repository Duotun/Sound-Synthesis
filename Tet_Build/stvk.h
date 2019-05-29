#pragma once
#include "Tensor.hpp"
#include<iostream>
#include<cmath>
#include <vector>
#include<algorithm>
#ifndef  STVK_H
#define  STVK_H


class StVKMaterial
{
public:
	StVKMaterial() {};
	~StVKMaterial() {};  //currently non definition
	StVKMaterial(double youngs, double poisson,double density) {
		mu_ = 0.5*youngs / (1 + poisson);
		lambda_ = youngs * poisson / (1 + poisson) / (1 - 2 * poisson);
		youngs_ = youngs;
		possion_ = possion_;
		density_ = density;
	};
	StVKMaterial(const StVKMaterial &mat) :lambda_(mat.lambda_), mu_(mat.mu_) {}
	StVKMaterial(const StVKMaterial *mat) :lambda_(mat->lambda_), mu_(mat->mu_) {}


	double get_mu() { return mu_; }
	void set_parameter(double l, double m)
	{
		lambda_ = l; mu_ = m;
	}
	void first_piola_kirchhoff(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const;
	void second_piola_kirchhoff(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const;
	void cauchy_stress_tensor(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const;

	//evaluate the stiffness matrix for a single tetrahedron
	void stiffness_matrix(const double density, const std::map<int, double> &mass, const std::vector<Eigen::Vector4i> &tetra, const std::vector<Eigen::Vector3d> &positions, Eigen::SparseMatrix<double, Eigen::RowMajor> &stiffmatrix, int pointsize) const;
	
	//strain energy density in undeformed volume yes! to deformed
	double strain_energy_density(const Eigen::MatrixX3d &F) const;   //check if how to decide? frobenius?
	double get_density() { return density_; }

private:
	double lambda_;
	double mu_;
	double youngs_;
	double possion_; 
	double density_;

};



#endif // ! STVK_H
