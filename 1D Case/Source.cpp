#include<iostream>
#include<stdio.h>
#include <math.h>
#include "Massspringsystem.h"
#include<SFML/Graphics.hpp>
#include<stdio.h>
#include<vector>
#include<SFML\Audio.hpp>
#include"Sound.h"
#include"extmat.hpp"
#include "2nd _linear.hpp"
#define maxforce 3000
//scale  - meters
// five Points, 4 Spring
// alloyed steel 2e-7 young's modulus
using std::cin; using std::cout;
//Wave is defined by the sqrt(k/m)
//K= y/r_restlength   

// better to use double
double dt;
double masses[5] = { 2.7e4,2.7e4,2.7e4,2.7e4,2.7e4 };  //5 points  gold   
double youngs[4] = { 7e10,7e10,7e10,7e10}; //4 springs
Eigen::Vector3cd positions[5] = {Eigen::Vector3cd(1,0,0),Eigen::Vector3cd(2,0,0),Eigen::Vector3cd(3,0,0) ,Eigen::Vector3cd(4,0,0),Eigen::Vector3cd(5,0,0)};
double restlength[4] = { 1,1,1,1};   //shold be very small


void init(massspringsystem &s);
void pushforce(Eigen::Vector3cd force, massspringsystem &s);
sf::SoundBuffer soundgeneration(double frequency, double amp);
//void getallfrequencies(std::vector<double>&frequency, massspringsystem &s);
void playsound(std::vector<double>frequency, std::vector<double> amp);

//position <-> sampling  ??
int main()
{

	massspringsystem systems;   dt = 0.1f;
	init(systems);
	cout << "Pre Computing......" << std::endl;
	Eigen::MatrixXd massmatrix;
	build_mass_matrix(systems, massmatrix);
	Eigen::MatrixXd stiffmatrix;
	build_stiff_matrix(systems, stiffmatrix);
	std::vector<double> frequencies;
	Eigen::MatrixXd D;  Eigen::MatrixXd U_T;
	eigendecomposition(stiffmatrix, massmatrix, D,frequencies,U_T);   // D = frequency^2, //general eigen decomposition
	Eigen::VectorXd force(systems.getmasses().size());
	cout << "Successfully get all the matrixs!!" << std::endl << std::endl;
	cout << frequencies[0] << std::endl; 
	cout << frequencies[1]<< std::endl; 
	cout << frequencies[2] << std::endl;
	cout << frequencies[3] << std::endl; 
	cout << frequencies[4] << std::endl;
	cout << "Input the force: " << std::endl;
	for (int i = 0; i < force.size(); i++)
	{
		
		cin >> force(i);
	}
	force = U_T * force;   //right side
	std::vector<double> amplitude;
	for (int i = 0; i < force.size(); i++)
	{
		solve_2nd_1d(D, force, amplitude, systems,i);
       
	}
	playsound(frequencies, amplitude);
	/*
	std::vector<double> frequencynow;
	getallfrequencies(frequencynow, systems);
	cout << std::endl;
	while (true) {
		cout << "Push Force: " << std::endl;
		double x, y, z;
		cin >> x; cin >> y; cin >> z;
		Eigen::Vector3cd force(x, y, z);

		pushforce(force, systems);
		cout << "Sound Generation: " << std::endl;
		sf::Sound sound[5];
		playsound(frequencynow);
		//std::vector<sf::SoundBuffer> buffer;
		
		//buffer.resize(frequencynow.size());  //size match with frequency
		//for (int i = 0; i<1; i++)
		//{
		//	buffer[i] = soundgeneration(frequencynow[i], 1.0);
		//}
		//

		
	}
	*/
	
	
	return 0;
}


void init(massspringsystem &s)   //add point+spring+fix or not
{
		for (int i = 0; i < 5; i++)
		{
			s.addmass(masses[i], positions[i].x().real(), positions[i].y().real(), positions[i].z().real(),i);
		}
		for (int i = 0; i < 5; i++)
		{
			if (i + 1 < 5)
			{
				s.addspring(youngs[i], restlength[i], s.massget(i), s.massget(i+1));
			}
		}
	
}

void pushforce(Eigen::Vector3cd force,massspringsystem &s)
{
	for (double i = 0.0; i < 1.0; i += dt)
	{
		s.update(dt, force);		
	}
	
}

sf::SoundBuffer soundgeneration(double frequency, double amp)  
{
	sf::SoundBuffer buffer;
	std::vector<sf::Int16> samples;

	for (int i = 0; i < 44100; i++)  //sampling 
	{
		samples.push_back(sound::Sinewave(i, 440, 0.9));
	}
	//one channel, sound per second
	buffer.loadFromSamples(&samples[0], samples.size(), 1, 44100);
	sf::Sound sound;
	sound.setBuffer(buffer);
	if (sound.getBuffer()!=NULL)
	{
		sound.play();
		cout << "sss" << std::endl;
	}
	return buffer;

	
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

void playsound(std::vector<double>frequency, std::vector<double> amp)    //combine sin waves not considering amplitudes
{
	const unsigned SAMPLES = 44100;
	const unsigned SAMPLE_RATE = 44100;
	const unsigned AMPLITUDE = 30000;

	sf::Int16 raw[SAMPLES];
	const double TWO_PI = 6.28318;
	double *increment = new double[frequency.size()];
	for (int i = 0; i < frequency.size(); i++)
	{
	increment[i] = frequency[i] / 44100;
	}
	//const double increment2 = 350. / 44110;
	double *x = new double[frequency.size()]{ 0 };
	for (unsigned i = 0; i < SAMPLES; i++) {
		double tmp = 0.0;
		for (int j = 0; j < frequency.size(); j++)
		{
			tmp += sin(x[j] * TWO_PI);
			x[j] += increment[j];
		}
		//tmp = sin(x[0] * TWO_PI);
		//x[0] += increment[0];
		tmp /= frequency.size();
		raw[i] = AMPLITUDE * tmp;

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
		//sf::sleep(sf::milliseconds(150));
		//Sound.stop();
		while (1) {
			sf::sleep(sf::milliseconds(100));
		}
		//sound->play();
		cout << "Play Sound Successfully!" << std::endl;
}


/*void getallfrequencies(std::vector<double>&frequency, massspringsystem &s)
{
	int num = 0;
	for (int i = 0; i < 5; i++)
	{
		pointmass *tmp = s.massget(i);
		
		for (int j = 0; j < tmp->frequencies.size(); j++)
		{
			double a = tmp->frequencies[j];
			frequency.resize(num + 1);
			frequency[num++] = a;
			cout << "Frequency: " << frequency[num - 1] << std::endl;
		}
	}
}
*/