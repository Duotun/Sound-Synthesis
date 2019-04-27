// MKL + Eigen are enough
#include "tetra.hpp"  
#include "extmat.hpp"
#include "matrixio.hpp"
#include "Modalmodel.h"
#include "sound.h"

#define min(x,y) (((x) < (y)) ? (x) : (y))   //min define, and can also define max

std::vector<Eigen::Vector3d> positions;
std::vector<int> ia[2];
std::vector<int> ja[2];
std::vector<double> data[2];
void simple_ui();
void construct_massmatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& massmatrix);
void construct_stiffmatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& stiffmatrix); 
void sound_synthesis();
void main()
{
	simple_ui();
}

void simple_ui()
{
	int cnt; 
	std::cout << "Simple UI - 0 to exit" << std::endl;
	const char* path_mass = "mass.csr"; const char * path_stiff = "stiff.csr";
	std::cin >> cnt;
	Eigen::SparseMatrix<double, Eigen::RowMajor> massmatrix;
	Eigen::SparseMatrix<double, Eigen::RowMajor> stiffmatrix;
	int numeigv;
	int rows_m = 0, cols_m = 0;   int rows_k = 0 , cols_k = 0;
	while (cnt)
	{
		switch(cnt)   //12 -io:mass , 34 io:stiff, 5-eigendecomposition, 6-sound synthesis
		{
		case 1: read_csr_matrix(path_mass, ia[0], ja[0], data[0], rows_m,cols_m); break;
		case 2: construct_massmatrix(massmatrix); break;
		case 3: read_csr_matrix(path_stiff, ia[1], ja[1], data[1], rows_k, cols_k); break;
		case 4: construct_stiffmatrix(stiffmatrix); break;
		case 5: numeigv = min(200, rows_k - 2); eigen_decomposition_mkl_feast(ia,ja,data,rows_k,numeigv); break;
		case 6: sound_synthesis(); break;
		}
		std::cin >> cnt;
	}
}
//lumped mass matrix
void construct_massmatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& massmatrix)
{
	if(positions.size()==0)
	tet_load_nodes(positions);
	std::vector<Eigen::Vector4i> tetra;  //index
	std::map<int, double> mass;   //corresponding mass
	tet_load_tetrahedron(tetra, mass, positions);   //get mass per tetra
	construct_mass_matrix(tetra, mass, massmatrix,positions.size()); //construct mass matrix

	const char* path = "mass.csr";
	write_csr_matrix(massmatrix, path);
	//write_txt_matrix(massmatrix, "mass.txt");
	std::cout << "Nonzeros of Mass Matrix: " << massmatrix.nonZeros() << std::endl;
}

void construct_stiffmatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& stiffmatrix)
{
	if (positions.size() == 0)
	tet_load_nodes(positions);
	std::vector<Eigen::Vector4i> tetra;  //index
	std::map<int, double> mass;   //corresponding mass
	tet_load_tetrahedron(tetra, mass, positions);   //here only need tetra
	std::vector<Eigen::Vector2i> edge;
	tet_load_edge(tetra, edge,positions.size());
	std::cout << "Size of Edge: " << edge.size() << std::endl;
	construct_stiff_matrix_spring_model(edge, positions, stiffmatrix, positions.size());
	
	const char* path = "stiff.csr";
	write_csr_matrix(stiffmatrix, path);
	//write_txt_matrix(stiffmatrix, "stiff.txt");
	std::cout << "Size of Stiff Matrix: " << stiffmatrix.rows() << std::endl;
}

void sound_synthesis()
{
	modalmodel m;  //read lambdas (eigenmodes)
	const char* path = "eig.csr";
	m.load_eigenmodes(path);
	//Eigen::VectorXd force(m.nrows_get()); force.setZero();
	std::vector<double> force; force.resize(m.nrows_get());
	double *out; force[3] = 1000;
	//force should build from impulse? maybe loop now
	//every 3 eigvec together 
	//cblas_dgemv(CblasRowMajor,CblasTrans,)

	int nrowk = m.nrows_get(); int numodes = m.num_modes();
	out = (double*)mkl_malloc(m.nrows_get() * sizeof(double), 64);
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, numodes, 1, nrowk, 1.0, m.eigv().data(), nrowk, force.data(), 1, 0.0, out, 1);
	//std::cout << "U_T*force: " << nrowk<<" "<<out[23] << std::endl;

	sound::soundplay_al(out, m, data[0],nrowk);

}