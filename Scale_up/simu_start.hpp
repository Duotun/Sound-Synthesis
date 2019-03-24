#pragma once
//#include <Eigen\Dense>
#include <iostream>
#include "Modal.h"
#include"Sound.h"
#include"extmat.hpp"
#include"export.hpp"
#include <string>
#include<stdio.h>
#include<vector>
//scale  - meters
using std::cout;

static double dt;   //wood density -400, pine wood young's modulus 9 gpa

//double area[4] = { 1e-4,1e-4,1e-4,1e-4 };   //1cm^2    ,areaûʲôӰ��
//double rho[4] = { 2.7e3,2.7e3,2.7e3,2.7e3 };  //5 points  gold and wood, right now wood  
//double youngs[4] = { 9e10,9e10,9e10,9e10 }; //4 springs
//Eigen::Vector3cd positions[5] = { Eigen::Vector3cd(1e-1,2e-1,3e-1),Eigen::Vector3cd(2e-1,3e-1,4e-1),Eigen::Vector3cd(3e-1,4e-1,5e-1) ,Eigen::Vector3cd(4e-1,5e-1,6e-1),Eigen::Vector3cd(5e-1,6e-1,7e-1) };  //m


double *area;
double *rho;
double *youngs;
Eigen::Vector3cd *positions;

void init(massspringsystem &s, int count)   //add point+spring+fix or not
{
	for (int i = 0; i < count; i++)
	{
		s.addmass(positions[i].x().real(), positions[i].y().real(), positions[i].z().real(), i);  //use position to get distance
	}
	for (int i = 0; i < count; i++)   //reset mass based on springs
	{
		if (i + 1 < count)   // increase with count
		{
			s.addspring(youngs[i], area[i], s.massget(i), s.massget(i + 1));
			s.reset_mass(rho[i], area[i], s.massget(i), s.massget(i + 1));    //3n*3n
		}
	}

}

//modeling
void init_parameters(int count)
{
	int i = 0;  //i for cycling
	area = new double[count];
	rho = new double[count];
	youngs = new double[count];
	positions = new Eigen::Vector3cd[count];
	for (int i = 0; i < count; i++)
	{
		area[i] = 1e-4;
		rho[i] = 2.7e3;
		youngs[i] = 9e10;
		positions[i] = Eigen::Vector3cd(i*1e-1, (i+1)*1e-1, (i+2)*1e-1);
	}
}

//currently not used
void init_3d(massspringsystem &s)   //add point+spring+fix or not   
{
	for (int i = 0; i < 5; i++)
	{
		s.addmass(positions[i].x().real(), positions[i].y().real(), positions[i].z().real(), i);  //use position to get distance
	}
	for (int i = 0; i < 5; i++)   //reset mass based on springs
	{
		if (i + 1 < 5)     //should be 3*n  masses
		{
			s.addspring(youngs[i], area[i], s.massget(i), s.massget(i + 1));
			s.reset_mass(rho[i], area[i], s.massget(i), s.massget(i + 1));    //3n*3n
		}
	}

}

void sound_1d_damping()  // difference lay in the stiffmatrix and mass matrix
{
	massspringsystem systems;   dt = 0.1f;
	init_parameters(5); //build_model
	init(systems,5);  //5 points
	cout << "Pre Computing......" << std::endl;
	Eigen::MatrixXd massmatrix;
	build_mass_matrix_1d(systems, massmatrix);
	Eigen::MatrixXd stiffmatrix;
	build_stiff_matrix_1d(systems, stiffmatrix);
	modalmodel mod;
	//Eigen::Matrix3Xd gradient;
	//calculate_gradient(systems, gradient);
	std::vector<double> frequencies;
	Eigen::MatrixXd S;  Eigen::MatrixXd U_T;
	std::vector<double> Force_Ud;
	eigendecomposition(stiffmatrix, massmatrix, S, frequencies, U_T);   // D = frequency^2, //general eigen decomposition
	Eigen::VectorXd force(systems.getmasses().size()); force.setZero();
	cout << "Successfully get all the matrixs!!" << std::endl << std::endl;

	cout << "All the Frequencies in the hearing range: " << std::endl;
	for (int i = 0; i < frequencies.size(); i++)
	{
		if (frequencies[i] != 0)
			cout << frequencies[i] << std::endl;
	}
	mod.load_eigenmodes(frequencies, systems);
	cout << std::endl;
	while (true)   //force input!
	{
		Force_Ud.resize(0);
		cout << "Force: (Assigned)" << std::endl;
		force.setZero();
		force(2) = 1000;   //first
		std::vector<double> amplitude;  double maxamp = 0;
		for (int i = 0; i < U_T.rows(); i++)   //  frequency ->amplitude, 15 dof
		{
			if (frequencies[i] != 0.) {   //cutting mode out of range
				double fi = U_T.row(i) * force;   //fi
												  //std::cout << "Size: " << Force_Ud.size() << std::endl;
				Force_Ud.push_back(fi);
			}

			//double tmp = amplitude[i];
			//if (abs(tmp) > maxamp)
			//maxamp = abs(tmp);
		}
		//cout << U_T << std::endl;
		//cout << Force_Ud.size() << std::endl << std::endl;
		cout << "Tap 'Enter' to trigger Sound or '0' to Exit" << std::endl << std::endl;
		std::string cc;
		std::getline(std::cin, cc);
		if (cc._Equal("0")) break;
		//cout << "max" << maxamp << std::endl;
		//amplitude_normalize(amplitude, maxamp);
		sound::soundplay_al_new(Force_Ud, mod, systems);
		//sound::soundplay_sfml(Force_Ud, mod, systems);
	}

}

void sound_1d_no_damping()
{
	massspringsystem systems;   dt = 0.1f;
	init(systems,5);
	cout << "Pre Computing......" << std::endl;
	Eigen::MatrixXd massmatrix;
	build_mass_matrix_1d(systems, massmatrix);
	Eigen::MatrixXd stiffmatrix;
	build_stiff_matrix_1d(systems, stiffmatrix);
	modalmodel mod;
	//Eigen::Matrix3Xd gradient;
	//calculate_gradient(systems, gradient);
	std::vector<double> frequencies;
	Eigen::MatrixXd S;  Eigen::MatrixXd U_T;
	std::vector<double> Force_Ud;
	eigendecomposition(stiffmatrix, massmatrix, S, frequencies, U_T);   // D = frequency^2, //general eigen decomposition
	Eigen::VectorXd force(systems.getmasses().size()); force.setZero();
	cout << "Successfully get all the matrixs!!" << std::endl << std::endl;

	cout << "All the Frequencies in the hearing range: " << std::endl;
	for (int i = 0; i < frequencies.size(); i++)
	{
		if (frequencies[i] != 0)
			cout << frequencies[i] << std::endl;
	}
	mod.load_eigenmodes(frequencies, systems);
	cout << std::endl;
	while (true)   //force input!
	{
		Force_Ud.resize(0);  //re-counting
		cout << "Force: (Assigned)" << std::endl;
		force.setZero();
		force(1) = 1000;   //first
		std::vector<double> amplitude;  double maxamp = 0;
		for (int i = 0; i < U_T.rows(); i++)   //  frequency ->amplitude, 15 dof
		{
			if (frequencies[i] != 0.) {  //cutting mode out of ranges
				double fi = U_T.row(i) * force;   //fi
				Force_Ud.push_back(fi);
			}

			//double tmp = amplitude[i];
			//if (abs(tmp) > maxamp)
			//maxamp = abs(tmp);
		}
		//cout << Force_Ud.size() << std::endl << std::endl;    //
		cout << "Tap 'Enter' to trigger Sound or '0' to Exit" << std::endl << std::endl;
		std::string cc;
		std::getline(std::cin, cc);
		if (cc._Equal("0")) break;
		//cout << "max" << maxamp << std::endl;
		//amplitude_normalize(amplitude, maxamp);
		sound::soundplay_al_new_nodamping(Force_Ud, mod, systems);
	}
}
// check if go in once
void damping_sound()   //need to be changed to only once
{
	    std::cout << "Input the size of Eige matrix: " << std::endl;
		int size_p;  std::cin >> size_p;
	    massspringsystem systems;   dt = 0.1f;
		init_parameters(size_p); //build_model		
		init(systems, size_p);
		cout << "Pre Computing......" << std::endl;
		//Eigen::SparseMatrix<double> massmatrix_sparse;
		//Eigen::SparseMatrix<double> stiffmatrix_sparse;
		//build_mass_matrix_sparse(systems, massmatrix_sparse);
		//build_stiff_matrix_sparse_3d(systems, stiffmatrix_sparse);
		Eigen::MatrixXd massmatrix;
		build_mass_matrix(systems, massmatrix);
		Eigen::MatrixXd stiffmatrix;
		stiffness_build_3d(systems, stiffmatrix);
		
		modalmodel mod;
		//Eigen::Matrix3Xd gradient;
		//calculate_gradient(systems, gradient);
		std::vector<double> frequencies;
		std::vector<double> eigenvalues;
		//std::vector<double> eigenvalues;
		//Eigen::SparseMatrix<double> ss;
		Eigen::MatrixXd S;  Eigen::MatrixXd U_T;
		std::vector<double> Force_Ud;
		output_tofile(stiffmatrix,massmatrix);

		std::cout << "Output to Matlab Successfully!" << std::endl;
		char ctmp2;   //wait for matlab to process
		std::cin >> ctmp2;
		import_eigenvalues(eigenvalues);   // record eigenvalues
		import_frequencies(frequencies);  //include zero
		import_U_matrix(U_T,size_p);
		//eigendecomposition_dense_spectra(stiffmatrix, massmatrix, S, frequencies, U_T);
		//eigendecomposition_sparse_spectra(stiffmatrix_sparse, massmatrix_sparse, frequencies, U_T);
		//eigendecomposition_sparse(stiffmatrix_sparse, massmatrix_sparse, frequencies, U_T);
		//eigendecomposition(stiffmatrix, massmatrix, S, frequencies, U_T);   // D = frequency^2, //general eigen decomposition
		Eigen::VectorXd force(systems.getmasses().size() * 3); force.setZero();
		cout << "Successfully get all the matrixs!!" << std::endl << std::endl;
		//std::cout << U_T.rows() << std::endl;
		
	cout << "All the Frequencies in the hearing range: " << std::endl;
	for (int i = 0; i < frequencies.size(); i++)
	{
		if (frequencies[i] != 0)
			cout << frequencies[i]<< std::endl;
	}

	//output to txt file eigenvalues
	output_tofile_non_zero(eigenvalues);
	//char ctmp;
	//std::cin >> ctmp;
	mod.load_eigenmodes(frequencies, systems);
	cout << std::endl;
	while (true)   //force input!
	{
		Force_Ud.resize(0);
		cout << "Force: (Assigned)" << std::endl;
		force.setZero();
		force(5) = 1000;   //first
		std::vector<double> amplitude;  double maxamp = 0;
		for (int i = 0; i < U_T.rows(); i++)   //  frequency ->amplitude, 15 dof
		{
			if (frequencies[i] != 0.) {
				double fi = U_T.row(i) * force;   //fi
				//std::cout << "F: " << fi << std::endl;
				Force_Ud.push_back(fi);
				//double tmp = amplitude[i];
				//if (abs(tmp) > maxamp)
				//maxamp = abs(tmp);
			}
		}
		//std::cout << U_T << std::endl;
		//std::cout << Force_Ud.size() << std::endl;
		cout << "Tap 'Enter' to trigger Sound or '0' to Exit" << std::endl << std::endl;
		std::string cc;
		std::getline(std::cin, cc);
		if (cc._Equal("0")) break;
		//cout << "max" << maxamp << std::endl;
		//amplitude_normalize(amplitude, maxamp);
		sound::soundplay_al_new(Force_Ud, mod, systems);
		//sound::soundplay_sfml(Force_Ud, mod, systems);
		//sound::soundplay_al_one_frequency(Force_Ud, mod, systems, 3);
		//sound::soundplay_al_one_frequency_nodamping(Force_Ud, mod, systems, 3);
		//sound::sound_test();
	}
	
}

void no_damping_sound()
{
	massspringsystem systems;   dt = 0.1f;
	init_parameters(5); //build_model	
	init(systems,5);
	
	cout << "Pre Computing......" << std::endl;
	Eigen::MatrixXd massmatrix;
	build_mass_matrix(systems, massmatrix);
	Eigen::MatrixXd stiffmatrix;
	stiffness_build_3d(systems, stiffmatrix);
	modalmodel mod;
	//Eigen::Matrix3Xd gradient;
	//calculate_gradient(systems, gradient);
	std::vector<double> frequencies;
	Eigen::MatrixXd S;  Eigen::MatrixXd U_T;
	std::vector<double> Force_Ud;
	
	eigendecomposition(stiffmatrix, massmatrix, S, frequencies, U_T);   // D = frequency^2, //general eigen decomposition
	Eigen::VectorXd force(systems.getmasses().size() * 3); force.setZero();
	cout << "Successfully get all the matrixs!!" << std::endl << std::endl;

	cout << "All the Frequencies in the hearing range: " << std::endl;
	for (int i = 0; i < frequencies.size(); i++)
	{
		if (frequencies[i] != 0)
			cout << frequencies[i] << std::endl;
	}
	mod.load_eigenmodes(frequencies, systems);
	cout << std::endl;
	while (true)   //force input!
	{
		Force_Ud.resize(0);
		cout << "Force: (Assigned)" << std::endl;
		force.setZero();
		force(0) = 1000;   //first
		std::vector<double> amplitude;  double maxamp = 0;
		for (int i = 0; i < U_T.rows(); i++)   //  frequency ->amplitude, 15 dof
		{
			if (frequencies[i] != 0.) {
				double fi = U_T.row(i) * force;   //fi
				Force_Ud.push_back(fi);
				//double tmp = amplitude[i];
				//if (abs(tmp) > maxamp)
				//maxamp = abs(tmp);
			}
		}
		//cout << Force_Ud[0] << std::endl << std::endl;
		cout << "Tap 'Enter' to trigger Sound or '0' to Exit" << std::endl << std::endl;
		std::string cc;
		std::getline(std::cin, cc);
		if (cc._Equal("0")) break;
		//cout << "max" << maxamp << std::endl;
		//amplitude_normalize(amplitude, maxamp);
		sound::soundplay_al_new_nodamping(Force_Ud, mod, systems);

	}
}