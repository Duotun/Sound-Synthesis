#pragma once
#include <iostream>
#include <Eigen\Dense>
#include "Massspringsystem.h"
#include <stdlib.h>
#include <math.h>



void calculate_gradient(massspringsystem &s, Eigen::Matrix3Xd &gradient)  //to identify the correctness of hessian
{
	std::vector<spring*> tmpspring = s.getspring();
	size_t numspring = tmpspring.size();
	//currently test first 3*3 matrix   
	gradient.resize(3, tmpspring.size());
	for (size_t i = 0; i < tmpspring.size(); i++) {
		pointmass *m1, *m2; m1 = tmpspring[i]->m1, m2 = tmpspring[i]->m2;
		Eigen::Vector3cd vary(1, 0, 0);
		Eigen::Vector3cd dis = (m1->position - m2->position + vary); double k = tmpspring[i]->k;
		double l0 = tmpspring[i]->get_restlen();
		double dist = dis.norm();  //distance ||x1-x2||
		Eigen::Vector3d tmpgradient = k * (dist - l0)*(dis / dist).real();
		//tmpgradient.addTo(gradient);
		gradient.col(i) = tmpgradient;
		//std::cout << gradient.col(i) << std::endl;
		//std::cout << "K is: " << k << std::endl;
		
	}
	std::cout << "Gradient is: " << std::endl << gradient << std::endl << std::endl;
}

void build_mass_matrix(massspringsystem &s, Eigen::MatrixXd &mass)   //diagonal,3d
{
	
	size_t numpoint = s.getmasses().size();
	mass.resize(numpoint*3, numpoint*3); mass.setZero();
	int cnt = 0;
	for (size_t i = 0; i < numpoint*3; i++)   //corresponding to the mesh size
	{
		if (i != 0 && i % 3 == 0) cnt++;
		mass(i,i) = s.massget(cnt)->mass;    //originally built diagonal
	    
	}
	//mass=mass.diagonal();
	std::cout << std::endl;
	std::cout << "Mass Matrix:" << std::endl;
	std::cout << mass<<std::endl;
}
void calculate_E(spring *sp, double dist, double a, double b, double c,Eigen::Matrix3d &matrix3)
{
	Eigen::Matrix3d tmpmatrix; tmpmatrix.setZero();
	double tmp = sp->get_restlen() / (pow(dist, 3)); 
	double d1, d2, d3;  //diagnal parameters
	d1 = 1.0 + a * a*tmp - sp->get_restlen() / dist;
	d2 = 1.0 + b * b*tmp - sp->get_restlen() / dist;
	d3 = 1.0 + c * c*tmp - sp->get_restlen() / dist;
	tmpmatrix(0, 0) = d1; tmpmatrix(1, 1) = d2; tmpmatrix(2, 2)= d3;
	tmpmatrix(0, 1) = tmpmatrix(1, 0) = a * b*tmp;
	tmpmatrix(0, 2) = tmpmatrix(2, 0) = a * c*tmp;
	tmpmatrix(1, 2) = tmpmatrix(2, 1) = a * c*tmp;
	matrix3 << tmpmatrix;
	//std::cout << "3*3 "<<matrix3 << std::endl<< std::endl;
}
// e -e e -e two points connected
void build_stiff_matrix(massspringsystem &s, Eigen::MatrixXd &stiff)
{
	// -k k+; /
	
	std::vector<spring*> tmpspring = s.getspring();
	size_t numpoint = s.getmasses().size();
	stiff.resize(numpoint, numpoint); stiff.setZero();
	for (size_t i = 0; i < tmpspring.size(); i++)
	{
		stiff(tmpspring[i]->m1->index, tmpspring[i]->m2->index) += -tmpspring[i]->k;
		stiff(tmpspring[i]->m2->index, tmpspring[i]->m1->index) += -tmpspring[i]->k;
		stiff(tmpspring[i]->m1->index, tmpspring[i]->m1->index) += tmpspring[i]->k;
		stiff(tmpspring[i]->m2->index, tmpspring[i]->m2->index) += tmpspring[i]->k;
	}
	stiff = -stiff;   //flip, 

	std::cout << std::endl;
	std::cout << "Stiff Matrix:" << std::endl;
	std::cout << stiff<<std::endl;
}

void stiffness_build_3d(massspringsystem &s,Eigen::MatrixXd &stiff)   //3d, E 2 differential
{
	stiff.resize(s.getmasses().size() * 3, s.getmasses().size() * 3); stiff.setZero();
	std::vector<spring*> tmpspring = s.getspring();
	size_t numpoint = s.getmasses().size(); int cnt = 0;
	for (size_t i = 0; i < tmpspring.size(); i++)   //*k
	{	
		pointmass *m1, *m2; double a, b, c; // a,b,c present three parts of dis vector
		m1 = tmpspring[i]->m1; m2 = tmpspring[i]->m2;  
		Eigen::Vector3cd dis = (m1->position - m2->position);
		//std::cout << dis << std::endl;
		double dist = dis.norm();
		a = dis(0).real(), b = dis(1).real(), c = dis(2).real();
		Eigen::Matrix3d tmpmatrix; tmpmatrix.setZero();calculate_E(tmpspring[i], dist, a, b, c,tmpmatrix);
		stiff.block<3, 3>(3*i, 3*i) += tmpspring[i]->k*tmpmatrix;
		stiff.block<3, 3>(3 * (i+1), 3 * (i+1)) += tmpspring[i]->k*tmpmatrix;
		stiff.block<3, 3>(3*i, 3*(i + 1)) += tmpspring[i]->k *-tmpmatrix;
		stiff.block<3, 3>(3*(i + 1), 3*i) += tmpspring[i]->k *-tmpmatrix;
		//std::cout << i << std::endl;		
	}
	std::cout << std::endl;
	std::cout << "Stiff Matrix:" << std::endl;
	std::cout << stiff << std::endl;
}

//Get D, G
void eigendecomposition(Eigen::MatrixXd &A, Eigen::MatrixXd &B, Eigen::MatrixXd & out, std::vector<double> &frequency, Eigen::MatrixXd &Ut)    //not general deomposition
{
	//real and part?
	//Eigen::EigenSolver<Eigen::MatrixXd>dec;
	Eigen::GeneralizedEigenSolver<Eigen::MatrixXd>dec;
	dec.compute(A,B);
	out.resize(A.rows(), A.cols()); 
	Eigen::VectorXd eigenvalueresults;
	eigenvalueresults.resize(dec.eigenvalues().real().size()); 
	//std::cout << dec.eigenvalues().real().size()<<std::endl <<std::endl;
	eigenvalueresults << dec.eigenvalues().real();    Ut.resize(dec.eigenvectors().rows(), dec.eigenvectors().cols()); Ut.setZero();
	//std::cout << "eigen_vectors" << dec.eigenvectors().cols()<<std::endl;
	size_t numpoint = eigenvalueresults.size();
	out.resize(numpoint, numpoint); out.setZero(); 
	//std::cout << numpoint << std::endl;
	//std::cout << "eigen_values: " << dec.eigenvalues().real() << std::endl << std::endl;
	for (size_t i = 0; i < numpoint; i++)   //corresponding to the mesh size
	{
		out(i, i) = eigenvalueresults(i);
		if(sqrt(abs(eigenvalueresults(i)))<=22000&&sqrt(abs(eigenvalueresults(i))>=2))
		       frequency.push_back(sqrt(abs(eigenvalueresults(i))));
		else   frequency.push_back(0.0);  //zero frequency - out of range
	}
	for (unsigned int i = 0; i < dec.eigenvectors().cols(); i++)
	{
		Ut.col(i) = dec.eigenvectors().real().col(i);
		//std::cout << dec.eigenvectors().real().col(i)<< std::endl<<std::endl;
	}
	//eigenvalueresults << dec.eigenvalues().real();
	//std::cout << out<<std::endl;
	//std::cout << frequency[0]<<std::endl;
	Ut = Ut.reverse();   //u^t
	//std::cout << "S "<<std::endl<<out << std::endl;
}


// get rid of frequencies that are out of the range: 2hz - 22000hz