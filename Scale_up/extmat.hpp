#pragma once
#include <Eigen\Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include "Massspringsystem.h"
#include <stdlib.h>
#include <math.h>
#include <Eigen/SparseCore>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>

// change to sparse matrix edition
typedef Eigen::Triplet<double> T;
using namespace Spectra;
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

void build_mass_matrix_sparse(massspringsystem &s, Eigen::SparseMatrix<double> &m) //diagonal, 3d, sparse
{
	size_t numpoint = s.getmasses().size();
	m.resize(numpoint * 3, numpoint * 3); m.setZero();
	int cnt = 0;
	for (size_t i = 0; i < numpoint * 3; i++)   //corresponding to the mesh size
	{
		if (i != 0 && i % 3 == 0) cnt++;
		m.insert(i, i) = s.massget(cnt)->mass;    //originally built diagonal

	}
	m.makeCompressed();
	std::cout << "Mass Matrix: " << m << std::endl;
}

void build_mass_matrix(massspringsystem &s, Eigen::MatrixXd &mass)   //diagonal,3d
{

	size_t numpoint = s.getmasses().size();
	mass.resize(numpoint * 3, numpoint * 3); mass.setZero();
	int cnt = 0;
	for (size_t i = 0; i < numpoint * 3; i++)   //corresponding to the mesh size
	{
		if (i != 0 && i % 3 == 0) cnt++;
		mass(i, i) = s.massget(cnt)->mass;    //originally built diagonal

	}
	//mass=mass.diagonal();
	//std::cout << std::endl;
	//std::cout << "Mass Matrix:" << std::endl;
	//std::cout << mass<<std::endl;
}

void build_mass_matrix_1d(massspringsystem &s, Eigen::MatrixXd &mass)   //diagonal,3d
{

	size_t numpoint = s.getmasses().size();
	mass.resize(numpoint, numpoint); mass.setZero();
	int cnt = 0;
	for (size_t i = 0; i < numpoint; i++)   //corresponding to the mesh size
	{
		mass(i, i) = s.massget(i)->mass;    //originally built diagonal

	}
	//mass=mass.diagonal();
	std::cout << std::endl;
	std::cout << "Mass Matrix:" << std::endl;
	std::cout << mass << std::endl;
}
void calculate_E(spring *sp, double dist, double a, double b, double c, Eigen::Matrix3d &matrix3)
{
	Eigen::Matrix3d tmpmatrix; tmpmatrix.setZero();
	double tmp = sp->get_restlen() / (pow(dist, 3));
	double d1, d2, d3;  //diagnal parameters
	d1 = 1.0 + a * a*tmp - sp->get_restlen() / dist;
	d2 = 1.0 + b * b*tmp - sp->get_restlen() / dist;
	d3 = 1.0 + c * c*tmp - sp->get_restlen() / dist;
	tmpmatrix(0, 0) = d1; tmpmatrix(1, 1) = d2; tmpmatrix(2, 2) = d3;
	tmpmatrix(0, 1) = tmpmatrix(1, 0) = a * b*tmp;
	tmpmatrix(0, 2) = tmpmatrix(2, 0) = a * c*tmp;
	tmpmatrix(1, 2) = tmpmatrix(2, 1) = a * c*tmp;
	matrix3 << tmpmatrix;
	//std::cout << "3*3 "<<matrix3 << std::endl<< std::endl;
}

// first calculate the corresponding dense matrix and convert to triplet
void calculate_E_sparse(std::vector<T> &tripletlist, const Eigen::MatrixXd &condensematrix)
{
	for (int i = 0; i < condensematrix.rows(); i++)
	{
		for (int j = 0; j < condensematrix.cols(); j++)
		{
			if(condensematrix(i,j)!=0.)
			tripletlist.push_back(T(i, j, condensematrix(i, j)));   //build triplet
		}
	}

}
// e -e e -e two points connected
void build_stiff_matrix_1d(massspringsystem &s, Eigen::MatrixXd &stiff)
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
	//stiff = -stiff;   //flip, 

					  //std::cout << std::endl;
					  //std::cout << "Stiff Matrix:" << std::endl;
					  //std::cout << stiff<<std::endl;
}

void stiffness_build_3d(massspringsystem &s, Eigen::MatrixXd &stiff)   //3d, E 2 differential
{
	stiff.resize(s.getmasses().size() * 3, s.getmasses().size() * 3); stiff.setZero();
	std::vector<spring*> tmpspring = s.getspring();
	size_t numpoint = s.getmasses().size(); int cnt = 0;
	for (size_t i = 0; i < tmpspring.size(); i++)   //*k
	{
		Eigen::Triplet<double> s;
		pointmass *m1, *m2; double a, b, c; // a,b,c present three parts of dis vector
		m1 = tmpspring[i]->m1; m2 = tmpspring[i]->m2;
		Eigen::Vector3cd dis = (m1->position - m2->position);
		//std::cout << dis << std::endl;
		double dist = dis.norm();
		a = dis(0).real(), b = dis(1).real(), c = dis(2).real();
		Eigen::Matrix3d tmpmatrix; tmpmatrix.setZero(); calculate_E(tmpspring[i], dist, a, b, c, tmpmatrix);
		stiff.block<3, 3>(3 * i, 3 * i) += tmpspring[i]->k*tmpmatrix;
		stiff.block<3, 3>(3 * (i + 1), 3 * (i + 1)) += tmpspring[i]->k*tmpmatrix;
		stiff.block<3, 3>(3 * i, 3 * (i + 1)) += tmpspring[i]->k *-tmpmatrix;
		stiff.block<3, 3>(3 * (i + 1), 3 * i) += tmpspring[i]->k *-tmpmatrix;
		//std::cout << i << std::endl;		
	}
	//stiff = -stiff;
	//std::cout << std::endl;
	//std::cout << "Stiff Matrix:" << std::endl;
	//std::cout << stiff << std::endl;
}

void build_stiff_matrix_sparse_3d(massspringsystem &s, Eigen::SparseMatrix<double> &stiff)   //sparse
{
	Eigen::MatrixXd stifftmp; 
	stiff.resize(s.getmasses().size() * 3, s.getmasses().size() * 3); stiff.setZero();
	stifftmp.resize(s.getmasses().size() * 3, s.getmasses().size() * 3); stifftmp.setZero();
	std::vector<spring*> tmpspring = s.getspring();
	size_t numpoint = s.getmasses().size(); int cnt = 0;
	std::vector<T> tripletlist;
	for (size_t i = 0; i < tmpspring.size(); i++)   //*k
	{
		pointmass *m1, *m2; double a, b, c; // a,b,c present three parts of dis vector
		m1 = tmpspring[i]->m1; m2 = tmpspring[i]->m2;
		Eigen::Vector3cd dis = (m1->position - m2->position);
		//std::cout << dis << std::endl;
		double dist = dis.norm();
		a = dis(0).real(), b = dis(1).real(), c = dis(2).real();
		Eigen::Matrix3d tmpmatrix; tmpmatrix.setZero(); calculate_E(tmpspring[i], dist, a, b, c, tmpmatrix);
		stifftmp.block<3, 3>(3 * i, 3 * i) += tmpspring[i]->k*tmpmatrix;
		stifftmp.block<3, 3>(3 * (i + 1), 3 * (i + 1)) += tmpspring[i]->k*tmpmatrix;
		stifftmp.block<3, 3>(3 * i, 3 * (i + 1)) += tmpspring[i]->k *-tmpmatrix;
		stifftmp.block<3, 3>(3 * (i + 1), 3 * i) += tmpspring[i]->k *-tmpmatrix;
		//std::cout << i << std::endl;		
	}
	calculate_E_sparse(tripletlist, stifftmp);
	stiff.setFromTriplets(tripletlist.begin(), tripletlist.end());
	//stiff = -stiff;
	stiff.makeCompressed();
	std::cout << "Stiff Matrix: " << std::endl << stiff << std::endl;
}

//Get D, G
void eigendecomposition(Eigen::MatrixXd &A, Eigen::MatrixXd &B, Eigen::MatrixXd & out, std::vector<double> &frequency, Eigen::MatrixXd &Ut)    //not general deomposition
{
	//real and part?
	//Eigen::EigenSolver<Eigen::MatrixXd>dec;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> dec(A,B);
	//Eigen::GeneralizedEigenSolver<Eigen::MatrixXd>dec;
	//dec.compute(A, B);
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
		
		//out(i, i) = eigenvalueresults(i);
		if (sqrt(abs(eigenvalueresults(i))) <= 22000 * 2 * EIGEN_PI&&sqrt(abs(eigenvalueresults(i)) >= 2 * 2 * EIGEN_PI))
		{
			frequency.push_back(sqrt(abs(eigenvalueresults(i))));
			out(i, i) = eigenvalueresults(i);
		}
		else   frequency.push_back(0.0);  //zero frequency - out of range
	}
	//out.conservativeResize(frequency.size(), frequency.size());
	//Ut.resize(A.rows(), frequency.size());
	for (unsigned int i = 0; i <dec.eigenvectors().cols(); i++)
	{
	    if(frequency[i]!=0.)
		Ut.col(i) = dec.eigenvectors().real().col(i);
		//std::cout << dec.eigenvectors().real().col(i)<< std::endl<<std::endl;
	}
	Eigen::MatrixXd tmp = Ut.transpose();

	std::cout << "Test: " << std::endl << tmp*A*Ut << std::endl;
	std::cout << "S: " << std::endl << out << std::endl;
	Ut.transposeInPlace();   //u^t
	//Ut = Ut.reverse();
	std::cout << "U_T "<<std::endl<<Ut<< std::endl;
}

// use sparsematrix
void eigendecomposition_sparse(Eigen::SparseMatrix<double> &A, Eigen::SparseMatrix<double> &B, std::vector<double> &frequency, Eigen::MatrixXd &Ut)    //not general deomposition
{
	A.makeCompressed(); B.makeCompressed();
	Eigen::GeneralizedEigenSolver<Eigen::MatrixXd>dec;
	dec.compute(A, B);
	
	Eigen::VectorXd eigenvalueresults;
	eigenvalueresults.resize(dec.eigenvalues().real().size());
	eigenvalueresults << dec.eigenvalues().real();    Ut.resize(dec.eigenvectors().rows(), dec.eigenvectors().cols()); Ut.setZero();
	size_t numpoint = eigenvalueresults.size();
	
	Eigen::SparseMatrix<double> S(A.rows(), A.rows()); S.setZero();
	//std::cout << "eigen_values: " << dec.eigenvalues().real() << std::endl << std::endl;
	for (size_t i = 0; i < numpoint; i++)   //corresponding to the mesh size
	{
	    S.insert(i,i)= eigenvalueresults(i);
		if (sqrt(abs(eigenvalueresults(i))) <= 22000 * 2 * EIGEN_PI&&sqrt(abs(eigenvalueresults(i)) >= 2 * 2 * EIGEN_PI))
			frequency.push_back(sqrt(abs(eigenvalueresults(i))));
			//frequency.push_back(eigenvalueresults(i));
		else   frequency.push_back(0.0);  //zero frequency - out of range
	}

	for (unsigned int i = 0; i < dec.eigenvectors().cols(); i++)
	{
		Ut.col(i) = dec.eigenvectors().real().col(i);
		//std::cout << dec.eigenvectors().real().col(i)<< std::endl<<std::endl;
	}
	//std::cout << "Stiff Matrix (back): " << std::endl;
	//std::cout << A-B * Ut*S*Ut.transpose() << std::endl;

	Eigen::MatrixXd tmp = Ut.transpose();
    
	std::cout << "Test: " << std::endl << tmp*A*Ut << std::endl;
	std::cout << "S: " << std::endl << S << std::endl;
	Ut << Ut.transpose();  //u^t
	//Ut.reser
	//std::cout << "U_T "<<std::endl<<Ut << std::endl;

}

// K is symmetric
#pragma warning(disable : 4996)
void eigendecomposition_sparse_spectra(Eigen::SparseMatrix<double> &A, Eigen::SparseMatrix<double> &B, std::vector<double> &frequency, Eigen::MatrixXd &Ut)    //not general deomposition
{   
  Spectra::SparseSymMatProd<double> op(A); //stiff
  Spectra::SparseCholesky<double> Bop(B);   //mass
  Spectra::SymGEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEIGS_CHOLESKY> geigs(&op, &Bop, A.rows()-1,A.rows());
  
  Ut.resize(A.rows()-1, A.cols()); Ut.setZero();
 // Eigen::MatrixXd tmp(A.rows(), A.cols());
  geigs.init();
  int nconv = geigs.compute();
  Eigen::MatrixXd S(A.rows() - 1, A.rows() - 1); S.setZero();
  if (geigs.info() == SUCCESSFUL)
  {
	  //tmp = geigs.eigenvectors();
	  Ut = geigs.eigenvectors();
	  Ut.transposeInPlace();
	  //Ut = tmp.transpose();   //n-1 *n
	  //std::cout << Ut << std::endl;
	 
	  for (int i = 0; i < geigs.eigenvalues().size(); i++)
	  {
		  S(i, i) = geigs.eigenvalues()(i);
		 if(sqrt(abs(geigs.eigenvalues()(i))) <= 22000 * 2 * EIGEN_PI&&sqrt(abs(geigs.eigenvalues()(i)) >= 2 * 2 * EIGEN_PI))
			  frequency.push_back(sqrt(abs(geigs.eigenvalues()(i))));
		 else   frequency.push_back(0.0);  //zero frequency - out of range
		  //frequency.push_back(geigs.eigenvalues()(i));
		  //std::cout << frequency[i] << std::endl;
	  }
  }
  //std::cout << S << std::endl;   //能通过验证
  //std::cout <<"First column: "<<geigs.eigenvectors().row(0)<< std::endl;
}

#pragma warning(disable : 4996)
void eigendecomposition_dense_spectra(Eigen::MatrixXd &A, Eigen::MatrixXd &B, Eigen::MatrixXd & out, std::vector<double> &frequency, Eigen::MatrixXd &Ut)    //not general deomposition	
{
	//std::cout << "A:" << std::endl << A << std::endl<<std::endl;
	//std::cout << "B:" << std::endl << B << std::endl << std::endl;
	Spectra::DenseSymMatProd<double> op(A);
	Spectra::DenseCholesky<double> bop(B);
	int min_val = std::min(A.rows() - 2, 1);
	Spectra::SymGEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double>, DenseCholesky<double>, GEIGS_CHOLESKY>
		es(&op, &bop, min_val, A.rows());
	es.init();
	es.compute();
	if (es.info() == SUCCESSFUL)
	{
		//std::cout << "hehehe" << std::endl;
		size_t numpoint = es.eigenvalues().size();
		out.resize(numpoint, numpoint); out.setZero();
		Eigen::VectorXd tmp = es.eigenvalues();
		//std::cout << tmp << std::endl;
		
		
		for (size_t i = 0; i < numpoint; i++)   //corresponding to the mesh size
		{
			//out(i, i) = eigenvalueresults(i);
			if (sqrt(abs(es.eigenvalues()(i))) <= 22000 * 2 * EIGEN_PI&&sqrt(abs(es.eigenvalues()(i)) >= 2 * 2 * EIGEN_PI))
			{
				frequency.push_back(sqrt(abs(es.eigenvalues()(i))));
				out(i, i) = es.eigenvalues()(i);
			}
			else   frequency.push_back(0.0);  //zero frequency - out of range
		}
		//out.conservativeResize(frequency.size(), frequency.size());
		//Ut.resize(A.rows(), frequency.size());
		Ut = es.eigenvectors();
		Ut.transposeInPlace();
	}
	//std::cout << "U_T " << std::endl << Ut << std::endl;
}

