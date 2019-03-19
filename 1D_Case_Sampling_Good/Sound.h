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
#include "Modal.h"
namespace sound {
#define two_PI 6.283185307

	//a single sine wave changes with time (one sound)
	void soundplay_sfml(std::vector<double> const &force, modalmodel &modal, massspringsystem &sys);
	short Sinewave(double time, double freq, double amplitude);
	void soundplay_al_new_nodamping(std::vector<double> const &force, modalmodel &modal, massspringsystem &sys);
	//void soundplay(Eigen::MatrixXd &S, std::vector<double>frequency, std::vector<double> amp, Eigen::VectorXd &force);
	void soundplay_al(Eigen::MatrixXd &S, std::vector<double>frequency, std::vector<double> amp);
	void soundplay_al_new(std::vector<double> const &force, modalmodel &modal, massspringsystem &sys);
	void soundplay_al_one_frequency(std::vector<double> const &force, modalmodel &modal, massspringsystem &sys,int tag);
	void soundplay_al_one_frequency_nodamping(std::vector<double> const &force, modalmodel &modal, massspringsystem &sys, int tag);
	void sound_test();
	//void exemplt_frequency(std::vector<double>&frequency);
	// 2-22000hz
}
#endif //
