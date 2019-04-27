#include "stvk.h"
#include "Tensor.hpp"
void  StVKMaterial::first_piola_kirchhoff(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const
{
	second_piola_kirchhoff(F, out);
	out = F * out;    //P=F*S, S- Second piola_kirchhoff
}
//Saint Venant-Kirchhoff model for Hyperelastic Materials
void StVKMaterial::second_piola_kirchhoff(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const
{
	Eigen::Matrix3d C, E;
	right_cauchy_green_deformation_tensor(F, C);
	green_strain_tensor(C, E);

	out = Eigen::Matrix3d::Identity()*E.trace() + 2 * mu_*E;
	//out =\lambda*tr(E)*I+2*\MU*E

}
// cauchy=1/J *F*S*F^T, S- second piolar stress tensor, J is the determinant of F
void StVKMaterial::cauchy_stress_tensor(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const
{
	second_piola_kirchhoff(F, out);
	out = F * out*F.transpose();
	out /= F.determinant();

}
double StVKMaterial::strain_energy_density(const Eigen::MatrixX3d &F) const
{
	Eigen::Matrix3d C, E;
	right_cauchy_green_deformation_tensor(F, C);
	green_strain_tensor(C, E);

	//\psi = \mu *tr(E^T*E)+0.5*\lambda*tr(E)^2;
	Eigen::Matrix3d colon = E.transpose()*E;
	const double tr = E.trace()*E.trace();  //tr(E)^2
	return mu_ * colon.trace() + 0.5*lambda_*tr;
}