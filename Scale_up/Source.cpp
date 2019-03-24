#pragma once
#include<iostream>
#include<stdio.h>
#include <math.h>
#include "Massspringsystem.h"
#include<SFML/Graphics.hpp>
#include<SFML\Audio.hpp>
#include <alc.h>
#include <al.h>
#include <Windows.h>
#include <iomanip>
#include "simu_start.hpp"
#define maxforce 3000
//scale  - meters
// five Points, 4 Spring
// alloyed steel 2e-7 young's modulus
using std::cin; using std::cout;
//K= y/r_restlength   

//double masses[5] = { 400,400,400,400,400 };  //5 points  gold and wood  
//double youngs[4] = { 9e9,9e9,9e9,9e9}; //4 springs

//double restlength[4] = { 1,1,1,1};   //shold be very small


void pushforce(Eigen::Vector3cd force, massspringsystem &s);
//void getallfrequencies(std::vector<double>&frequency, massspringsystem &s);
void playsound(std::vector<double>frequency, std::vector<double> amp);
void simpleui();
void random_test();   //test eigen efficiency
void random_test_spectra();
//position <-> sampling  ??

int main()
{
	simpleui();
	//random_test();
	//random_test_spectra();
	return 0;
}

void random_test()
{
  Eigen::MatrixXd X = Eigen::MatrixXd::Random(500, 500);
  Eigen::MatrixXd A = X + X.transpose();
	//cout << "Here is a random symmetric matrix, A:" << std::endl << A << std::endl;
	X = Eigen::MatrixXd::Random(500, 500);
	Eigen::MatrixXd B = X * X.transpose();
	//cout << "and a random postive-definite matrix, B:" << std::endl << B << std::endl <<std::endl;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A, B);
	std::cout << "Completed" << std::endl << std::endl;
	//cout << "The eigenvalues of the pencil (A,B) are:" << std::endl << es.eigenvalues() << std::endl;
	//cout << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl << std::endl;

}

void random_test_spectra()
{
	Eigen::MatrixXd X = Eigen::MatrixXd::Random(3000, 3000);
	Eigen::MatrixXd A = X + X.transpose();
	X = Eigen::MatrixXd::Random(3000, 3000);
	Eigen::MatrixXd B = X * X.transpose();

	Spectra::DenseSymMatProd<double> op(A);
	Spectra::DenseCholesky<double> bop(B);
	Spectra::SymGEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double>, DenseCholesky<double>, GEIGS_CHOLESKY>
		es(&op, &bop, 498, 500);
	std::cout << "Completed" << std::endl << std::endl;
}
void simpleui()
{
	//int goincnt = 0;
	while (true)
	{
		std::cout << "Choose Modals:  " << std::endl << std::endl;
		std::cout << "0. Exit" << std::endl;
		std::cout << "1. 1D - No Damping" << std::endl;
		std::cout << "2. 1D - Damping" << std::endl;
		std::cout << "3. 3D - No Damping" << std::endl;
		std::cout << "4. 3D - Damping" << std::endl << std::endl;
		int cnt = 0;
		std::cin >> cnt;
		std::string cc;
		std::getline(std::cin, cc);   //delete enter
		switch (cnt)
		{
		case 1: sound_1d_no_damping(); break;
		case 2: sound_1d_damping(); break;
		case 3: no_damping_sound(); break;
		case 4: damping_sound(); break;
		default: break;

		}
		if (cnt <= 0) break;
		std::cout << std::endl << "Play the Sound Successfully!" << std::endl;
	}


}

void pushforce(Eigen::Vector3cd force, massspringsystem &s)
{
	for (double i = 0.0; i < 1.0; i += dt)
	{
		s.update(dt, force);
	}

}

/*void playsound(std::vector<double>frequency, std::vector<double> amp)
{
const unsigned SAMPLES = 44100;
const unsigned SAMPLE_RATE = 44100;
const unsigned AMPLITUDE = 33000;

sf::Int16 raw[SAMPLES];
const double TWO_PI = 6.28318;
double *increment = new double[frequency.size()];
for (int i = 0; i < frequency.size(); i++)
{
increment[i] = frequency[i] / 44100;
}
//const double increment2 = 350. / 44110;
double *x = new double[frequency.size()]{ 0 }; int cnt = 0;
while (cnt<frequency.size())
{
double tmp = 0.0;
for (unsigned i = 0; i < SAMPLES; i++)
{
tmp = cos(x[cnt] * TWO_PI);
x[cnt] += increment[cnt];

}
raw[cnt] += amp[cnt] * AMPLITUDE*tmp;
cnt++;
}

//cout << "psps" << std::endl;
sf::SoundBuffer Buffer;
if (!Buffer.loadFromSamples(raw, SAMPLES, 1, SAMPLE_RATE)) {
std::cerr << "Loading failed!" << std::endl;
return ;
}

sf::Sound Sound;
Sound.setBuffer(Buffer);
Sound.setLoop(false);
Sound.play();
sf::sleep(sf::milliseconds(150));
Sound.stop();
//while (1) {
//	sf::sleep(sf::milliseconds(100));
//}
//sound->play();
cout << "Play Sound Successfully!" << std::endl;
}
*/

void playsound(std::vector<double>frequency, std::vector<double> ampnow)    //combine sin waves not considering amplitudes
{
	const unsigned SAMPLES = 44100;
	const unsigned SAMPLE_RATE = 44100;
	const unsigned AMPLITUDE = 32767;

	sf::Int16 raw[SAMPLES];
	const double TWO_PI = 6.28318;
	double *increment = new double[frequency.size()];
	for (int i = 0; i < frequency.size(); i++)
	{
		increment[i] = frequency[i] / 44100;
		if (ampnow[i] > 32767)
			ampnow[i] = 32767;
		else if (ampnow[i] < -32767)
			ampnow[i] = 32767;
	}
	//const double increment2 = 350. / 44110;
	double *x = new double[frequency.size()]{ 0 };
	for (unsigned i = 0; i < SAMPLES; i++) {
		double tmp = 0.0;
		double amptmp = 0.0;
		for (int j = 0; j < frequency.size(); j++)
		{
			if ((tmp + sin(x[j] * TWO_PI)*ampnow[j]) <= 32767 && (tmp + sin(x[j] * TWO_PI)*ampnow[j]) >= -32767)
				tmp += sin(x[j] * TWO_PI)*abs(ampnow[j]);
			else if ((tmp + sin(x[j] * TWO_PI)*ampnow[j]) > 32767)
			{
				tmp = 32767;  cout << "ppp" << std::endl;
			}
			else if ((tmp + sin(x[j] * TWO_PI)*ampnow[j]) < -32767)
			{
				tmp = -32767;
			}
			x[j] += increment[j];
		}
		//tmp = sin(x[0] * TWO_PI);
		//x[0] += increment[0];
		//tmp /= frequency.size();
		raw[i] = tmp;

	}
	//cout << "psps" << std::endl;
	sf::SoundBuffer Buffer;
	if (!Buffer.loadFromSamples(raw, SAMPLES, 1, SAMPLE_RATE)) {
		std::cerr << "Loading failed!" << std::endl;
		return;
	}

	sf::Sound Sound;
	Sound.setBuffer(Buffer);
	Sound.setLoop(false);
	Sound.play();
	sf::sleep(sf::seconds(1));
	Sound.stop();
	cout << "Play Sound Successfully!" << std::endl;
}


/*
/*
ALCdevice *device;
ALCcontext *context;
ALshort data[BUFFER_LENGTH * 2];
ALuint buffer, source;
int i;

// Initialization
device = alcOpenDevice(NULL);
context = alcCreateContext(device, NULL);
alcMakeContextCurrent(context);
alGenBuffers(1, &buffer);

// Generate sine wave data
for (i = 0; i < BUFFER_LENGTH; ++i) {
data[i * 2] += sin(2 * EIGEN_PI * SOUND_HZ * i / BUFFER_LENGTH) * SHRT_MAX;
data[i * 2 + 1] += -1 * sin(2 * EIGEN_PI * SOUND_HZ * i / BUFFER_LENGTH) * SHRT_MAX; // antiphase
data[i * 2] += sin(2 * EIGEN_PI * SOUND_HZ *2* i / BUFFER_LENGTH) * SHRT_MAX;
data[i * 2 + 1] += -1 * sin(2 * EIGEN_PI *2* SOUND_HZ * i / BUFFER_LENGTH) * SHRT_MAX; // antiphase
}

// Output looping sine wave
alBufferData(buffer, AL_FORMAT_STEREO16, data, sizeof(data), BUFFER_LENGTH * 2);
alGenSources(1, &source);
alSourcei(source, AL_BUFFER, buffer);
alSourcei(source, AL_LOOPING, AL_TRUE);
alSourcePlay(source);

// Wait to exit
printf("Press any key to exit.");
getchar();

// Finalization
alSourceStop(source);
alDeleteSources(1, &source);
alDeleteBuffers(1, &buffer);
alcMakeContextCurrent(NULL);
alcDestroyContext(context);
alcCloseDevice(device);

return 0;
*/

