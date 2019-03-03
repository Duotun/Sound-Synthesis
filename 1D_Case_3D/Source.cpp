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
#include <alc.h>
#include <al.h>
#include <Windows.h>
#include <iomanip>
#include <string>
#define maxforce 3000
//scale  - meters
// five Points, 4 Spring
// alloyed steel 2e-7 young's modulus
using std::cin; using std::cout;
//Wave is defined by the sqrt(k/m)
//K= y/r_restlength   

// better to use double
double dt;   //wood density -400, pine wood young's modulus 9 gpa
//double masses[5] = { 400,400,400,400,400 };  //5 points  gold and wood  
//double youngs[4] = { 9e9,9e9,9e9,9e9}; //4 springs
double area[4] = { 1e-4,1e-4,1e-4,1e-4 };   //1cm^2
double rho[4] = { 4e2,4e2,4e2,4e2 };  //5 points  gold and wood  
double youngs[4] = { 9e9,9e9,9e9,9e9 }; //4 springs
Eigen::Vector3cd positions[5] = {Eigen::Vector3cd(1,2,3),Eigen::Vector3cd(2,3,4),Eigen::Vector3cd(3,4,5) ,Eigen::Vector3cd(4,5,6),Eigen::Vector3cd(5,6,7)};  //m
//double restlength[4] = { 1,1,1,1};   //shold be very small


void init(massspringsystem &s);
void pushforce(Eigen::Vector3cd force, massspringsystem &s);
sf::SoundBuffer soundgeneration(double frequency, double amp);
//void getallfrequencies(std::vector<double>&frequency, massspringsystem &s);
void playsound(std::vector<double>frequency, std::vector<double> amp);
void amplitude_normalize(std::vector<double> &ampli, double maxamp);
//position <-> sampling  ??

#define SECOND 1
#define SAMPLING_HZ 44100
#define BUFFER_LENGTH (SECOND * SAMPLING_HZ)
#define SOUND_HZ 14706

int main()
{
	massspringsystem systems;   dt = 0.1f;
	init(systems);
	cout << "Pre Computing......" << std::endl;
	Eigen::MatrixXd massmatrix;
	build_mass_matrix(systems, massmatrix);
	Eigen::MatrixXd stiffmatrix;
	stiffness_build_3d(systems, stiffmatrix);
	//Eigen::Matrix3Xd gradient;
	//calculate_gradient(systems, gradient);
	std::vector<double> frequencies;
	Eigen::MatrixXd S;  Eigen::MatrixXd U_T;
	
	eigendecomposition(stiffmatrix, massmatrix, S,frequencies,U_T);   // D = frequency^2, //general eigen decomposition
	Eigen::VectorXd force(systems.getmasses().size()*3); force.setZero();
	cout << "Successfully get all the matrixs!!" << std::endl << std::endl;

	cout << "All the Frequencies in the hearing range: " << std::endl;
	for (int i = 0; i < frequencies.size(); i++)
	{
		if(frequencies[i]!=0)
		cout << frequencies[i] << std::endl;
	}
	cout << std::endl;
	while (true)   //force input!
	{
		cout << "Force: (Assigned)" << std::endl;
		force.setZero();
		force(1) = 1;   //first

		//for (int i = 0; i < force.cols(); i++)
		//{  
			//randomize right now
			//force<<1000,0,0,0,0
			//Eigen::VectorXd tmp; tmp.resize(systems.getmasses().size() * 3);
			//tmp.setOnes();  force.col(i) = tmp;
		//}
		//std::cout<<std::endl<< force << std::endl<<std::endl;
		   //right side
		std::vector<double> amplitude;  double maxamp = 0;
	
		for (int i = 0; i < U_T.rows(); i++)   //  frequency ->amplitude, 15 dof
		{
			double fi = U_T.row(i) * force;   //fi
			solve_2nd_1d(S, fi, frequencies,amplitude,systems,i);   // amplitude, frequency? dimension?
			double tmp = amplitude[i];
			if (abs(tmp) > maxamp)
				maxamp = abs(tmp);
		}

		cout << "Tap 'Enter' to trigger Sound " << std::endl << std::endl;
		std::string cc;
		std::getline(cin, cc);
		//cout << "max" << maxamp << std::endl;
		amplitude_normalize(amplitude, maxamp);
		sound::soundplay_al(S, frequencies, amplitude);
			
}
	
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
	memset(data, 0, sizeof(data));
	for (i = 0; i < BUFFER_LENGTH; ++i) {
		//data[i * 2] += sin(2 * EIGEN_PI * SOUND_HZ * i / BUFFER_LENGTH) * SHRT_MAX;
		//data[i * 2 + 1] += -1 * sin(2 * EIGEN_PI * SOUND_HZ * i / BUFFER_LENGTH) * SHRT_MAX; // antiphase
		data[i * 2] += sin(2 * EIGEN_PI * 300 * 2 * i / BUFFER_LENGTH) * SHRT_MAX;
		data[i * 2 + 1] += -1 * sin(2 * EIGEN_PI * 2 * 300 * i / BUFFER_LENGTH) * SHRT_MAX; // antiphase
	}

	// Output looping sine wave
	alBufferData(buffer, AL_FORMAT_STEREO16, data, sizeof(data), BUFFER_LENGTH * 2);
	alGenSources(1, &source);
	alSourcei(source, AL_BUFFER, buffer);
	alSourcei(source, AL_LOOPING, AL_FALSE);
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

	
	*/
	return 0;
}


void amplitude_normalize(std::vector<double> &ampli,double maxamp)
{
	cout << "Amplitude of All Synthesized Sounds: " << std::endl;
	for (int i = 0; i < ampli.size(); i++)
	{
		double tmp = abs(ampli[i]) / maxamp;
		tmp = round(tmp * 1000) / 1000.0;    // 0.001 precision
		ampli[i] = tmp;
		if(tmp!=0.0)
		cout << tmp << std::endl;		
	}
	cout << std::endl;
}

void init(massspringsystem &s)   //add point+spring+fix or not
{
		for (int i = 0; i < 5; i++)
		{
			s.addmass(positions[i].x().real(), positions[i].y().real(), positions[i].z().real(),i);  //use position to get distance
		}
		for (int i = 0; i < 5; i++)   //reset mass based on springs
		{
			if (i + 1 < 5)
			{
				s.addspring(youngs[i],area[i],s.massget(i), s.massget(i+1));
				s.reset_mass(rho[i], area[i], s.massget(i), s.massget(i + 1));    //3n*3n
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
			if((tmp+sin(x[j]*TWO_PI)*ampnow[j])<=32767&&(tmp + sin(x[j] * TWO_PI)*ampnow[j])>=-32767)
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
		return ;
		}

		sf::Sound Sound;
		Sound.setBuffer(Buffer);
		Sound.setLoop(false);
		Sound.play();
		sf::sleep(sf::seconds(1));
		Sound.stop();
		//while (1) {
		//	sf::sleep(sf::milliseconds(100));
		//}
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

