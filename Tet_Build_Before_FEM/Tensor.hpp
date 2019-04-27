#pragma once
#include<Eigen\Dense>
#include<Eigen\Core>

//double form
inline void right_cauchy_green_deformation_tensor(const Eigen::Matrix3d &F,Eigen::Matrix3d &C)
{
	C = F.transpose()*F;
}
//E=0.5*(F^T*F-I)
inline void green_strain_tensor(const Eigen::Matrix3d &C, Eigen::Matrix3d &E)
{
	E = (C - Eigen::Matrix3d::Identity())*0.5;
}