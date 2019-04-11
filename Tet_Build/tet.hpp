//for tet load and compute
#pragma once
#include<mkl.h>
#define EIGEN_USE_MKL_ALL
#include<Eigen\Core>
#include<Eigen\Dense>
#include<Eigen\Sparse>
#include<iostream>
#include <fstream>
std::string filepath=".\\Tetgen_exe\\";

const double density=2.7e3; 

double scale = 1.0e-1; //currently for 1e-1 meter

double tet_volume(const std::vector<Eigen::Vector3d> tet);  //declaration
double Calculation_V_L(const std::vector<Eigen::Vector3d> tet);

void tet_load_nodes(std::vector<Eigen::Vector3d> &positions)
{
	std::fstream file(filepath + "plate-surf.1.node", std::ios_base::in);
	std::string lines;

	// first 3 parameters
	int index; int index_tmp; //for every line
	double pos[3];
	
	if (file.is_open()) {   //opened check
		if (getline(file, lines))   //read index first
		{
			std::istringstream is(lines);
			is >> index; int tmp;
			while (is >> tmp);
		}
		//std::cout << "Index: " << index << std::endl;


		while (getline(file, lines))
		{
			if (lines[0] == '#') break;
			std::istringstream is(lines);
			is >> index_tmp;
			for (int i = 0; i < 3; i++)
			{
				is >> pos[i];
			}
			positions.push_back(Eigen::Vector3d(pos[0], pos[1], pos[2]));
		}
		//check logic
		if (positions.size() == index) std::cout << "Read Correctly!" << std::endl;
		file.close();
	}
	else
	{
		std::cout << "Error: " << std::endl;
	}
}

//map for volume
void tet_load_tetrahedron(std::vector<Eigen::Vector4i> &tetra,std::map<int,double> &mass, const std::vector<Eigen::Vector3d> &positions)
{
	
	std::fstream file(filepath + "plate-surf.1.ele", std::ios_base::in);
	std::string lines;

	// first 3 parameters
	int index; int index_tmp; //for every line
	int point_index[4];
	if (file.is_open()) {   //opened check
		if (getline(file, lines))   //read index first
		{
			std::istringstream is(lines);
			is >> index; int tmp;
			while (is >> tmp);
		}
		//std::cout << "Index: ele" << index << std::endl;
		int cnt=0; //for map
		while (getline(file, lines))
		{
			if (lines[0] == '#') break;
			std::vector<Eigen::Vector3d> forvolume;  //for volume computing
			std::istringstream is(lines);
			is >> index_tmp;  //line title 
			for (int i = 0; i < 4; i++)
			{
				is >> point_index[i];
				forvolume.push_back(positions[point_index[i]]);
				
			} //tetra input
		
			tetra.push_back(Eigen::Vector4i(point_index[0], point_index[1], point_index[2], point_index[3]));
			mass[cnt++] = tet_volume(forvolume)*density;  //get the mass of tetrahedron
			//std::cout << "Tetrahedron " << cnt << ": " << " V/L: " << Calculation_V_L(forvolume) << std::endl;
			//std::cout << "Mass: " << mass[cnt - 1] << std::endl;
		}
		if (tetra.size() == index) std::cout << "Read Correctly!" << std::endl;
		file.close();
	}
	else
	{
		std::cout << "Load Error!" << std::endl;
	}

}

// calculate v/l -area testing
double Calculation_V_L(const std::vector<Eigen::Vector3d> tet)
{
	double volume = tet_volume(tet);
	double length = 0; 
	for (int i = 0; i < tet.size(); i++)
	{
		for (int j = i + 1; j < tet.size(); j++)   //
		{
			length += (tet[j] - tet[i]).norm();   //add
			
		}
	}
	return volume / length;
}

//volume computation
double tet_volume(const std::vector<Eigen::Vector3d> tet)   //1/6 determinant of for vertices
{
	Eigen::Matrix4d determinant;   //4*4
	for (int i = 0; i < 4; i++) {
		Eigen::Vector4d tmp;                 //1*4
		tmp << tet[i], 1;
		determinant.col(i) << tmp;  
		
	}
	
	//std::cout <<std::endl<<"Volume Matrix:"<< std::endl<<determinant << std::endl;
	
	return abs(determinant.determinant());  //volume positive
}

void edge_check(Eigen::SparseMatrix<int,Eigen::RowMajor> &spmat_edge, int index_i, int index_j)
{
	//insert once
	spmat_edge.coeffRef(index_i, index_j) = 1;  

}

//get edges from tet, please avoid repeating i-j, j-i edge
void tet_load_edge(const std::vector<Eigen::Vector4i> &tet, std::vector<Eigen::Vector2i> &edge,int pointsize)
{
	Eigen::SparseMatrix<int,Eigen::RowMajor> spmat_edge;
	spmat_edge.resize(pointsize, pointsize);  //n*n matrix
	
	for (int i = 0; i < tet.size(); i++)   //matrix to store edge
	{
		edge_check(spmat_edge, tet[i](0), tet[i](1)); edge_check(spmat_edge, tet[i](1), tet[i](0));
		edge_check(spmat_edge, tet[i](0), tet[i](2)); edge_check(spmat_edge, tet[i](2), tet[i](0));
		edge_check(spmat_edge, tet[i](0), tet[i](3)); edge_check(spmat_edge, tet[i](3), tet[i](0));
		edge_check(spmat_edge, tet[i](1), tet[i](2)); edge_check(spmat_edge, tet[i](2), tet[i](1));
		edge_check(spmat_edge, tet[i](1), tet[i](3)); edge_check(spmat_edge, tet[i](3), tet[i](1));
		edge_check(spmat_edge, tet[i](2), tet[i](3)); edge_check(spmat_edge, tet[i](3), tet[i](2));
		//std::cout << "Here" << std::endl;
	}
	spmat_edge.makeCompressed();

	
	//check
	//Eigen::MatrixXi tmp = spmat_edge.toDense();
	//std::cout << "Matrix: " << std::endl << tmp;

	//only upper of triangle sparse matrix
	int cnt = 0;
	Eigen::SparseMatrix<int, 1> spupper = spmat_edge.triangularView<Eigen::Upper>();
	for (int k = 0; k< spupper.outerSize(); ++k)
		for (Eigen::SparseMatrix<int,1>::InnerIterator it(spupper, k); it; ++it)
		{
			edge.push_back(Eigen::Vector2i(it.row(),it.col()));
		}
	std::cout <<"Load Edges Successfully!" << std::endl;  //18 or 19
	
}
