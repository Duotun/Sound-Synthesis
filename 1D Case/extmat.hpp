#pragma once
#include <iostream>
#include <Eigen\Dense>
#include "Massspringsystem.h"
#include <stdlib.h>
#include <math.h>
void build_mass_matrix(massspringsystem &s, Eigen::MatrixXd &mass)   //diagonal
{

	int numpoint = s.getmasses().size();
	mass.resize(numpoint, numpoint); mass.setZero();
	for (int i = 0; i < numpoint; i++)   //corresponding to the mesh size
	{
		mass(i,i) = s.massget(i)->mass;    //originally built diagonal
	}
	//mass=mass.diagonal();
	//std::cout << mass;
}

void build_stiff_matrix(massspringsystem &s, Eigen::MatrixXd &stiff)
{
	// -k k+; /
	
	std::vector<spring*> tmpspring = s.getspring();
	int numpoint = s.getmasses().size();
	stiff.resize(numpoint, numpoint); stiff.setZero();
	for (int i = 0; i < tmpspring.size(); i++)
	{
		stiff(tmpspring[i]->m1->index, tmpspring[i]->m2->index) += -tmpspring[i]->k;
		stiff(tmpspring[i]->m2->index, tmpspring[i]->m1->index) += -tmpspring[i]->k;
		stiff(tmpspring[i]->m1->index, tmpspring[i]->m1->index) += tmpspring[i]->k;
		stiff(tmpspring[i]->m2->index, tmpspring[i]->m2->index) += tmpspring[i]->k;
	}
	stiff = -stiff;   //flip, 
	//std::cout << stiff<<std::endl;
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
	eigenvalueresults.resize(dec.eigenvalues().rows()); 
	eigenvalueresults << dec.eigenvalues().real();    Ut.resize(eigenvalueresults.size(), eigenvalueresults.size());
	int numpoint = eigenvalueresults.size();
	out.resize(numpoint, numpoint); out.setZero(); 
	for (int i = 0; i < numpoint; i++)   //corresponding to the mesh size
	{
		out(i, i) = eigenvalueresults(i);
		frequency.push_back(sqrt(abs(eigenvalueresults(i))));
	}
	//eigenvalueresults << dec.eigenvalues().real();
	//std::cout << out<<std::endl;
	//std::cout << frequency[0]<<std::endl;
	Ut<< dec.eigenvectors().real().inverse();
	//std::cout << Ut;
	//std::cout << out << std::endl;
}


// get rid of frequencies that are out of the range: 2hz - 22000hz