#include "Modal.h"

void modalmodel::load_eigenmodes(std::vector<double> &frequency, massspringsystem &s)
{
	eigenmodes.resize(frequency.size());
	int nmd = 0; int zerocnt = 0;
	for (; nmd < frequency.size(); nmd++)
	{
		eigenmodes[nmd] = frequency[nmd];
		//std::cout <<eigenmodes[nmd] << std::endl;
		if (frequency[nmd] == 0.0)  zerocnt++;
	}
	nummodes = nmd; nmd = zerocnt;   //change nummodes later
	std::cout << "Non-Zerocnt: " << frequency.size()-zerocnt << std::endl;
	omega_.resize(nummodes - zerocnt);
	omega_d.resize(nummodes - zerocnt);
	c_.resize(nummodes - zerocnt);

	//avoid zero computing  
	zerocnt = 0;
	for (int i = 0; i < nummodes; i++)
	{

		if (eigenmodes[i] != 0.) {
			omega_[zerocnt] = eigenmodes[i];
			c_[zerocnt] = beta * eigenmodes[i] * eigenmodes[i] + alpha;
			double xi = c_[zerocnt] / (2 * omega_[zerocnt]);    //epsilon
																//std::cout << xi << std::endl;
			assert(xi <= 1);
			omega_d[zerocnt] = omega_[zerocnt] * sqrt(1.0 - xi * xi);
			zerocnt++;
			//std::cout << "Omega: " << omega_[zerocnt] << std::endl;
		}
	}
	nummodes -= nmd;
}