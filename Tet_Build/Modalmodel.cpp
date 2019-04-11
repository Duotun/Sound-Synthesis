#include<mkl.h>
#include "Modalmodel.h"
#include<Eigen\Core>
const double cuttingmode = 16000 * 2*EIGEN_PI;
int  modalmodel::load_eigenmodes(const char* path)
{
	std::ifstream fin(path, std::ios::binary);
	int M, nrowk;
	fin.read((char*)&M, sizeof(int));
	fin.read((char*)&nrowk, sizeof(int));
	nummodes = M; nrows = nrowk;
	if (fin.fail()) {
		printf("Cannot read file: %s\n", path);
		return 0;
	}
	eigenmodes.resize(M); //cutting try
	fin.read((char*)&eigenmodes[0], sizeof(double)*nummodes);
	int nmds = 0;
	for (; nmds < nummodes; ++nmds)
	{
		if (eigenmodes[nmds] > cuttingmode*cuttingmode) break;
	}
	nummodes = nmds;   //may reduce
	eigenvec.resize(nrowk*nummodes);   //current number of eigenvec
	fin.read((char *)&eigenvec[0], sizeof(double)*eigenvec.size());

	//compute omega right now
	omega_.resize(nummodes);
	omega_d.resize(nummodes);
	freq.resize(nummodes);
	c_.resize(nummodes);

	for (int i = 0; i < nummodes; i++)
	{
		omega_[i] = sqrt(eigenmodes[i]);
		freq[i] = omega_[i] * 0.5 / EIGEN_PI;
		c_[i] = beta * eigenmodes[i] + alpha;
		double xi = c_[i] / (2.*omega_[i]);

		if (xi >= 1.)
		{
			printf("c[%d] should always be in the range [0, 1]: %lf\n", i, xi);
			return 0;
		}
		omega_d[i] = omega_[i] * sqrt(1. - xi * xi);
	}
	printf("Load Eigenmodes Successfully!");
	return 1;
}

void modalmodel::accm_modal_impulse()   // compute force, // out += U'*imp/rho, rho?
{

}