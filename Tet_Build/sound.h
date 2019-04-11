#pragma once
#ifndef  SOUND_H
#define SOUND_H

#include <Windows.h>
#include <vector>
#include <alc.h>
#include <al.h>
#include <math.h>
#include <mkl.h>
#define EIGEN_USE_MKL_ALL
#include <Eigen\Dense>
#include "Modalmodel.h"
namespace sound {
	//force after 
	void soundplay_al(double *force,modalmodel &modal, std::vector<double> &mass, int nrowk);

}
#endif