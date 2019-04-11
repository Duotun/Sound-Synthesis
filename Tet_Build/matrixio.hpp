#pragma once   //for outputing mass matrix and stiff matrix
#include<mkl.h>
#define EIGEN_USE_MKL_ALL
#include<Eigen\Core>
#include<Eigen\Dense>
#include<Eigen\Sparse>

// csr format  - data and idx:col are the same size, 

void write_csr_matrix(Eigen::SparseMatrix<double,1> &spmat,const char *file)  // file write in binary mode
{
	std::ofstream fout;
	fout.open(file, std::ios::binary);
	if (!fout.good())
	{
		std::cerr << "Fail to open file to write!" << std::endl;
		return;
	}
	//first write row, col, numzeros
	int rows = spmat.rows(), cols = spmat.cols(), n = spmat.nonZeros();
	fout.write((char*)&rows, sizeof(int));
	fout.write((char*)&cols, sizeof(int));
	fout.write((char*)&n, sizeof(int));

	//then csr format, size of ia is (row+1), the last one is size of row
	int* ia = spmat.outerIndexPtr();
	int* ja = spmat.innerIndexPtr();
	const double *data = spmat.valuePtr();
	
	// from 0-based to 1-based
	for (int i = 0; i < rows + 1; i++) ia[i] += 1;
	for (int i = 0; i < n; i++) ja[i] += 1;

	fout.write((const char*)ia, sizeof(int)*(rows+1));  
	fout.write((const char*)ja, sizeof(int)*n);
	fout.write((const char*)data, sizeof(double)*n);

	fout.close();

}
//read csr into three vectors
void read_csr_matrix(const char *file,std::vector<int> &ia, std::vector<int>&ja,std::vector<double> &data, int &rows, int &cols)
{
	std::ifstream fin;
	fin.open(file, std::ios::binary);
	if (fin.fail()) {
		std::cout << "read_csc_dmatrix:: Cannot open file [%s] to read\n";
		return;
	}
	int n; //for data
	fin.read((char *)&rows, sizeof(int));
	fin.read((char *)&cols, sizeof(int));
	fin.read((char *)&n, sizeof(int));

	ia.resize(rows + 1);
	ja.resize(n);
	data.resize(n);

	fin.read((char*)(ia.data()), sizeof(int)*(rows + 1));
	fin.read((char*)(ja.data()), sizeof(int)*(n));
	fin.read((char*)(data.data()), sizeof(double)*(n));

	fin.close();
	std::cout << "Read Matrix Successfully!" << std::endl;

}

void write_txt_matrix(Eigen::SparseMatrix<double, 1> &spmat, const char *file)
{
	std::ofstream fout;
	fout.open(file, std::ios_base::out);
	if (!fout.good())
	{
		std::cerr << "Fail to open file to write!" << std::endl;
		return;
	}
	//first write row, col, numzeros
	int rows = spmat.rows(), cols = spmat.cols(), n = spmat.nonZeros();
	fout <<rows<< " "; fout << cols << " "; fout << n << " "<<"\n";

	//Eigen::MatrixXd nma= spmat.toDense().triangularView<Eigen::Upper>();
	//Eigen::MatrixXd m= spmat.toDense().triangularView<Eigen::Lower>();
	//fout << m<<"\n"<<"\n";
	fout << spmat.toDense();
	fout.close();


}