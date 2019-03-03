#include "Sound.h"
#include<SFML\Audio.hpp>
#include <Eigen\Dense>
#include <algorithm>


#define SECOND 1
#define SAMPLING_HZ 44100
#define BUFFER_LENGTH (SECOND * SAMPLING_HZ)
#define SOUND_HZ 14706

static const int SR = 44100;
static const int TS = 3; 
namespace sound {

	void soundplay(Eigen::MatrixXd &S, std::vector<double>frequency ,std::vector<double> amp, Eigen::VectorXd &force)
	{
		const unsigned SAMPLES = 44100;
		const unsigned SAMPLE_RATE = 44100;
		const unsigned AMPLITUDE = 32767;

		sf::Int16 raw[SR*TS] = { 0 };
		const int totticks = SR * TS;

		for (int i = 0; i < frequency.size(); i++)
		{
			const double SS = force[i] / S(i,i);
			for (int ti = 0; ti < totticks; ti++)
			{
				const double ts = static_cast<double>(ti) / static_cast<double>(SR);
				const double amp = 1.0;   //later consider damping
				raw[ti] += amp*SS*sin(frequency[i] * ts);   //3 s
			}
		}
		
	     double normalscale = 0.0;
		//normalize
		double AMP = 0;
		for (int ti = 0; ti < totticks; ++ti)
			AMP = std::fmax(AMP, abs(raw[ti]));
		normalscale = 0.6 / AMP;
		
		sf::SoundBuffer Buffer;
		for (int ti = 0; ti < totticks; ++ti)
			raw[ti] = raw[ti] * normalscale*32767;
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
		//while (1) {
		//	sf::sleep(sf::milliseconds(100));
		//}
	}

	short Sinewave(double time, double freq, double amplitude)
	{
		short result;
		double tpc = 44100 / freq;   //ticks per cycle for sampling  self-define sampling rate based on the frequency
		double cycles = time / tpc;
		double rad = two_PI * cycles;
		short tmpamplitude = 32767 * amplitude;   //amplitude can be 0-1

		result = tmpamplitude * sin(rad);

		return result;
	}
	// 
	void soundplay_al(Eigen::MatrixXd &S, std::vector<double>frequency, std::vector<double> ampli)
	{
		ALCdevice *device;
		ALCcontext *context;
		ALshort data[BUFFER_LENGTH * 2];
		ALuint buffer, source;

		//int i;


		// Initialization
		device = alcOpenDevice(NULL);
		if (device)
		{
			context = alcCreateContext(device, NULL);
			alcMakeContextCurrent(context);
		}
		alGetError();
		alGenBuffers(1, &buffer);

		if (alGetError() != AL_NO_ERROR)
		{
			std::cout << "error" << std::endl;
			return ;
		}
		memset(data, 0, sizeof(data));   //zero +
		// Generate sine wave data    //
		for (int j = 0; j < frequency.size(); j++) 
		{
			for (int i = 0; i < BUFFER_LENGTH; ++i) 
			{
				data[i * 2] += sin(2 * EIGEN_PI * frequency[j] * i / BUFFER_LENGTH) * ampli[j] * SHRT_MAX;   //normallized amplitude
				data[i * 2 + 1] += -1 * sin(2 * EIGEN_PI * frequency[j] * i / BUFFER_LENGTH) * ampli[j] * SHRT_MAX; // antiphase
			}
		}

		// Output looping sine wave
		//for (int i = 0; i < frequency.size(); i++)
	//	{
			alBufferData(buffer, AL_FORMAT_STEREO16, data, sizeof(data), BUFFER_LENGTH * 2);
			alGenSources(1, &source);
			alGetError();  //clear error code
			alSourcei(source, AL_BUFFER, buffer);
			alSourcei(source, AL_LOOPING, AL_FALSE);
			//alsourcef(source[i],AL_SEC_OFFSET)

		//}
		
		alSourcePlay(source);
		ALint s = AL_PLAYING;
		while (s == AL_PLAYING)
		{
			alGetSourcei(source, AL_SOURCE_STATE, &s);
			Sleep(50);   //50 ms
			alSourceStop(source);
		}
		//alSourcePlay(source2);
		// Wait to exit
		//printf("Press any key to exit.");
		//getchar();
		// Finalization

		alDeleteSources(1, &source);
		alDeleteBuffers(1, &buffer);
		alcMakeContextCurrent(NULL);
		alcDestroyContext(context);
		alcCloseDevice(device);

	}
	
}