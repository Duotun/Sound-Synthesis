#pragma once    //for constructing mass matrix and stiff matrix
#include<mkl.h>
#define EIGEN_USE_MKL_ALL
#include<Eigen\Core>
#include<Eigen\Dense>
#include <Eigen\Sparse>
#include <algorithm>

typedef Eigen::Triplet<double>T;
const double youngs = 9e10;
const double area = 1e-1;
void construct_mass_matrix(const std::vector<Eigen::Vector4i> &tetra,const std::map<int, double> &mass, Eigen::SparseMatrix<double, Eigen::RowMajor> &massmatrix, int pointsize)  //use sparse directly
{
	std::vector<double> sum_mass; sum_mass.resize(pointsize, 0);
	
	int index_1, index_2, index_3, index_4;  //tetra
	for (int i = 0; i < tetra.size(); i++)   //spreading on every point
	{
		index_1 = tetra[i](0); index_2 = tetra[i](1); 
		index_3 = tetra[i](2); index_4 = tetra[i](3);
		sum_mass[index_1] += mass.at(i)/4.0; sum_mass[index_2] += mass.at(i) / 4.0;
		sum_mass[index_3] += mass.at(i) / 4.0; sum_mass[index_4] += mass.at(i) / 4.0;
	}
	std::vector<T>coefficients;  //make diagonal sparse mass matrix
	for (int i = 0; i < pointsize; i++)
	{
		for(int j=0;j<3;j++)    //3n *3n
		coefficients.push_back(T(3*i+j, 3*i+j, sum_mass[i]));
	}
	massmatrix.resize(pointsize*3, pointsize*3);  //need to express size explicitly
	massmatrix.setFromTriplets(coefficients.begin(), coefficients.end());  //already compressed
	std::cout << "Construct Matrix Successfully!" <<std::endl;
}

//3*3
void Hessian_insert_two_point(std::vector<T> &triplelist,int index_i,int index_j,double d1,double d2, double d3, 
double a,double b, double c, double tmp,int tag,double k)  //k hook constant
{
	if (tag > 0) {  //positive
		triplelist.push_back(T(3 * index_i, 3 * index_j, d1*k));
		triplelist.push_back(T(3 * index_i + 1, 3 * index_j + 1, d2*k));
		triplelist.push_back(T(3 * index_i + 2, 3 * index_j + 2, d3*k));

		triplelist.push_back(T(3 * index_i, 3 * index_j + 1, a*b*tmp*k));
		triplelist.push_back(T(3 * index_j + 1, 3 * index_i, a*b*tmp*k));
		triplelist.push_back(T(3 * index_i, 3 * index_j + 2, a*c*tmp*k));
		triplelist.push_back(T(3 * index_j + 2, 3 * index_i, a*c*tmp*k));
		triplelist.push_back(T(3 * index_i + 1, 3 * index_j + 2, b*c*tmp*k));
		triplelist.push_back(T(3 * index_j + 2, 3 * index_i + 1, b*c*tmp*k));
	}
	else   // negative
	{
		triplelist.push_back(T(3 * index_i, 3 * index_j, -d1*k));
		triplelist.push_back(T(3 * index_i + 1, 3 * index_j + 1, -d2*k));
		triplelist.push_back(T(3 * index_i + 2, 3 * index_j + 2, -d3*k));

		triplelist.push_back(T(3 * index_i, 3 * index_j + 1, -a*b*tmp*k));
		triplelist.push_back(T(3 * index_j + 1, 3 * index_i, -a*b*tmp*k));
		triplelist.push_back(T(3 * index_i, 3 * index_j + 2, -a*c*tmp*k));
		triplelist.push_back(T(3 * index_j + 2, 3 * index_i, -a*c*tmp*k));
		triplelist.push_back(T(3 * index_i + 1, 3 * index_j + 2, -b*c*tmp*k));
		triplelist.push_back(T(3 * index_j + 2, 3 * index_i + 1, -b*c*tmp*k));
	}
}

void calculate_Hessian_Sparse(std::vector<T> &triplelist,int index_i, int index_j, double dist, double a, double b, double c,double k) 
{
	
	double tmp = 1 / (pow(dist, 2));
	double d1, d2, d3;   //get all parameters
	d1 = 1.0 + a * a*tmp - 1.0;
	d2 = 1.0 + b * b*tmp - 1.0;
	d3 = 1.0 + c * c*tmp - 1.0;
	
	//if (a < 0) std::cout << "A<0:" << std::endl; a->vector
	// currently manual add....
	//i i
	Hessian_insert_two_point(triplelist, index_i, index_i, d1, d2, d3, a, b, c, tmp,1,k);
	
	//j j
	Hessian_insert_two_point(triplelist, index_j, index_j, d1, d2, d3, a, b, c, tmp, 1,k);

	//i j
	Hessian_insert_two_point(triplelist, index_i, index_j, d1, d2, d3, a, b, c, tmp, -1,k);

	//j i
	Hessian_insert_two_point(triplelist, index_j, index_i, d1, d2, d3, a, b, c, tmp, -1,k);

	//Hessian.makeCompressed();
}

// positions, edges, areas, currently spring models
void construct_stiff_matrix_spring_model(const std::vector<Eigen::Vector2i> &edge,const std::vector<Eigen::Vector3d> &positions, Eigen::SparseMatrix<double, Eigen::RowMajor> &stiffmatrix,int pointsize)
{
	stiffmatrix.resize(positions.size() * 3, positions.size() * 3); stiffmatrix.setZero();
	std::cout << "Size: " << positions.size()*3 << std::endl;
	int edgenum = edge.size();
	std::vector<T> buildlist;  
	for (int i = 0; i < edgenum; i++)
	{
		Eigen::Vector3d p1, p2;  double a, b, c; //parameters for hessian
		p1 = positions[edge[i](0)];  p2 = positions[edge[i](1)];
		Eigen::Vector3d dis = p2 - p1; double dist = dis.norm();
		a = dis(0); b = dis(1); c = dis(2);
		Eigen::SparseMatrix<double, Eigen::RowMajor> tmpspmat; 
		//std::cout << "Edge: " << edge[i](0) << " " << edge[i](1) << std::endl;
		double k = youngs * area / dist;
		calculate_Hessian_Sparse(buildlist, edge[i](0), edge[i](1), dist, a, b, c,k);  //need to be changed
		
		
		//std::cout << "Value: " << stiffmatrix.coeffRef(0, 0) << std::endl;
	}
	stiffmatrix.setFromTriplets(buildlist.begin(), buildlist.end());
	//stiffmatrix.makeCompressed();  //compressed
	//stiffmatrix = -stiffmatrix;
	//Eigen::MatrixXd  t=stiffmatrix.toDense();
	std::cout << "Construct StiffMatrix Successfully!" <<std::endl; //<<t<<std::endl;
}

//utilize mkl libraries, nomalize vectors or not?
int eigen_decomposition_mkl_feast(std::vector<int> ia[2],
std::vector<int> ja[2], std::vector<double> data[2], int nrowk,int numeigv)
{ 
	
	int freqLow = 2; int freqHigh = 22000;
	int info = 3;   //check if eigendecomposition is successfull
	double epsout;   //error 
	int loop;
	const char UPLO = 'F'; //upper part or full part?
	MKL_INT feastparam[128];
	feastinit(feastparam);

	int M; //total numbers found in the interval
	
	const double emin = freqLow * 2.*EIGEN_PI*freqLow * 2.*EIGEN_PI;  //whether * density, currently I think is not
	const double emax = freqHigh * 2.*EIGEN_PI*freqHigh * 2.*EIGEN_PI;
	//std::cout << "Emax: " << emax << std::endl;
	int M0 = numeigv; 
	std::vector<double> eval(M0), evec(M0*nrowk), res(M0);
	// size?
	//std::cout << "Data: " << ia[0][24]<< std::endl;   
	//std::cout << "Data: " << ia[1].size() << std::endl;

	dfeast_scsrgv(&UPLO,&nrowk,data[1].data(),ia[1].data(),ja[1].data(),
		data[0].data(), ia[0].data(), ja[0].data(),feastparam,&epsout,&loop
	,&emin,&emax,&M0,eval.data(),evec.data(),&M,res.data(),&info);

	//dfeast_scsrev(&UPLO, &nrowk, data[1].data(), ia[1].data(), ja[1].data(),feastparam, &epsout, &loop
		//, &emin, &emax, &M0, eval.data(), evec.data(), &M, res.data(), &info);

	switch (info)
	{
	case 0: printf("Successfully run FEAST!\n"); for(int i=0;i<M;i++) std::cout << "Eval: " << eval[i] << std::endl;
		break;
	case 1:
		printf("Failed to run FEAST. [code = %d]\n", info);
		printf("no eigenvalue found in the search interval\n");
		return 1;
	case 3:
		printf("Failed to run FEAST. [code = %d]\n", info);
		printf("Size of M0 (=%d) is too small\n", numeigv);
		return 1;
	default:
		printf("Failed to run FEAST. [code = %d]\n", info);
		return 1;
	}
	

	//output eigenvalues and eigenvectors, rightnow directly to a matrix (Cube is small)
	std::cout <<std::endl<< "Number of Eigenvalues: " << M << std::endl;
	
	const char* path = "eig.csr";   //eig results;
	std::ofstream fout;
	fout.open(path, std::ios::binary);
    
	//M - founded eigenvalues, nrowk - size of matrix(eigvector)

	if (!fout.good())
	{
		std::cerr << "Fail to open file to write!" << std::endl;
		return 0;
	}

	fout.write((char*)&M, sizeof(int));    //first write size, here different from Dingzeyu, he first writes nrowk
	fout.write((char*)&nrowk, sizeof(int));  

	for (int vid = 0; vid < M; vid++)
	{
		fout.write((char*)&eval[vid], sizeof(double));
		//printf("ev#%3d:  %lf %lfHz\n", vid, eval[vid], sqrt(eval[vid])*0.5/EIGEN_PI);
	}

	for (int vid = 0; vid < M; vid++)  //need to be checked
	{
		fout.write((char *)&evec[vid*nrowk], sizeof(double)*nrowk);
	}

	fout.close();

}



