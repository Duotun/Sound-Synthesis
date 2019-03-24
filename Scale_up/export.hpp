#pragma once
#include <iostream>
#include <fstream>
#include <Eigen\Dense>
#include <stdlib.h>
#include <stdio.h>
#include <string>
void output_tofile(const Eigen::MatrixXd &K, const Eigen::MatrixXd &M)
{
	std::ofstream file_m;
	std::ofstream file_k;

	file_m.open("mass.txt", std::ios::out | std::ios::trunc);
	file_k.open("stiff.txt", std::ios::out | std::ios::trunc);
	if (file_m.is_open()&&file_k.is_open()) {
		//std::cout << "here" << std::endl;
		file_k << K;
		file_m << M;
	}
	file_m.close();
	file_k.close();


}

void import_eigenvalues(std::vector<double> &eigenvalues)
{
	std::ifstream myfile;
	myfile.open("eigenvalues.txt", std::ios::in);
	int tmpvalue;
	if (myfile.is_open())
	{
		double num = 0.0;
		std::string s;
		while (!std::getline(myfile,s).eof())
		{
			//if (num >= 5e-1)   //avoid too small
			//{
			num=std::stod(s);
			if(num>=1.)
			eigenvalues.push_back(num);
			//}
		}
	}
}

void import_frequencies(std::vector<double> &eigenvalues)
{
	std::ifstream myfile;
	myfile.open("eigenvalues.txt", std::ios::in);
	int tmpvalue;
	if (myfile.is_open())
	{
		double num = 0.0;
		std::string s;
		while (!std::getline(myfile, s).eof())
		{
			//if (num >= 5e-1)   //avoid too small
			//{
			num = std::stod(s); 
			if (num >= 1.)
			{
				num = sqrt(num);   //change into frequency
				eigenvalues.push_back(num);
			}
			else
			{
				eigenvalues.push_back(0.);
			}
			//}
		}
	}
}

void import_U_matrix(Eigen::MatrixXd &U, int size)
{
	std::ifstream myfile;
	myfile.open("U_matrix.txt", std::ios::in);
	int tmpvalue;
	U.resize(size*3, size*3 - 2); U.setZero();   //3d
	if (myfile.is_open())
	{
		double num = 0.0;
		std::string s;
        int cnt = 0; int col = 0;   //count for input
		while(!std::getline(myfile, s).eof())
		{
			num = std::stod(s);
			
			U(cnt++, col)=num;
			if (cnt == size*3)
			{
				cnt = 0;
				col++;   //next column
			}
		}
	}
	U.transposeInPlace();   //transpose itself
	//std::cout << "U_T: " << std::endl << U << std::endl;
}
void output_tofile_non_zero(const std::vector<double> &eigenvalues)  //need to be put non-zero vectors
{
	std::ofstream file_e;
	file_e.open("eigenvalues_nonzero.txt", std::ios::out | std::ios::trunc);
	Eigen::VectorXd tmp = Eigen::VectorXd::Map(eigenvalues.data(), eigenvalues.size());
	if (file_e.is_open())
	{
		file_e << tmp;
	}
	file_e.close();
}
