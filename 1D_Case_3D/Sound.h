#pragma once
#ifndef  SOUND_H
#define SOUND_H

#include <math.h>
#include <Eigen\Dense>
#include <iostream>
#include <Windows.h>
#include <vector>
#include <alc.h>
#include <al.h>

namespace sound {
#define two_PI 6.283185307

	//a single sine wave changes with time (one sound)
	short Sinewave(double time, double freq, double amplitude);

	void soundplay(Eigen::MatrixXd &S, std::vector<double>frequency, std::vector<double> amp, Eigen::VectorXd &force);
	void soundplay_al(Eigen::MatrixXd &S, std::vector<double>frequency, std::vector<double> amp);
	//void exemplt_frequency(std::vector<double>&frequency);
	// 2-22000hz
}
#endif //
