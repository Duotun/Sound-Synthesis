#include "sound.h"
#include <algorithm>

#define SECOND 3
#define SAMPLING_HZ 44100
#define BUFFER_LENGTH (SECOND * SAMPLING_HZ)

static const int SR = 44100;
static const int TS = 3;

namespace sound
{
	//new method, currently old to produce sound first
	void soundplay_al(double *force, modalmodel &modal,std::vector<double> &mass,int nrowk) //csr (diag)
	{
		ALCdevice *device;
		ALCcontext *context;
		ALshort data[BUFFER_LENGTH];
		
		ALuint buffer, source;
		const std::vector<double> &omegaD = modal.damped_omega();
		const std::vector<double> &c = modal.damping_co();
		//const std::vector<double> mass = sys.masses_3d();     //masses

		const int totticks = SR * TS;   //sample the sound for 3 second 
        
		// Initialization
		device = alcOpenDevice(NULL);
		context = alcCreateContext(device, NULL);
		alcMakeContextCurrent(context);
		alGetError();
		alGenBuffers(1, &buffer);

		if (alGetError() != AL_NO_ERROR)
		{
			std::cout << "error" << std::endl;
			return;
		}
		//0 every time
		memset(data, 0, sizeof(data));   //zero +
										 // Generate sine wave data    //
		//std::cout << "Number of Modes: " << modal.num_modes()<<std::endl;
		//std::cout << modal.num_modes() << std::endl << std::endl;
		//std::cout << "force size: " << force.size() << std::endl;

		for (int j = 0; j < omegaD.size(); j++)  //all if force!= 0
		{
			
			if (round(force[j]) == 0) { continue; }
			//std::cout <<omegaD[j] << std::endl;
			//Currently not condiering distance between cam and object
			const double SS = force[j] / omegaD[j]; //need to be updated

																//ampvec.push_back(SS);  store amp

			//if (ampvec[j] <= 4.5e9) continue; //cut
			//zerocnt process
			//std::cout << "C: "<<c[j] << std::endl;
		    //std::cout << force[j]<< std::endl;
			std::cout << "Omega_D: " << omegaD[j] << " C: " << c[j] << std::endl;
			for (int i = 0; i < BUFFER_LENGTH; ++i)
			{
				const double ts = static_cast<double>(i) / static_cast<double>(SR);
				const double amp = exp(-c[j] * 0.5*ts);

				///check
				//if (ts > 2.)   //check too long sampling
				//{
				//	std::cout << ">2: " << "Omega_d: " << omegaD[j]<< " C: "<<c[j]<<std::endl; break;
				//}
				
				if (amp < 1e-5) {
					//std::cout << "Omega_d: " << omegaD[j] << "  T: " << ts <<" C: " << c[j] << std::endl;
					break;
				}  //cut amplitude
				   ///check end

				data[i] += sin(omegaD[j] * ts) * amp*SS*SHRT_MAX;   //normallized amplitude 
																	//data[2*i+1] += -sin(omegaD[j] * ts) * amp*SS*SHRT_MAX;   //normallized amplitude 
				if (data[i] > SHRT_MAX) std::cout << "Error Value" << std::endl;
			}
		}

		alBufferData(buffer, AL_FORMAT_MONO16, data, sizeof(data), SAMPLING_HZ);   //stereo shortens the sound (two channels)
		alGenSources(1, &source);
		alGetError();  //clear error code
		alSourcei(source, AL_BUFFER, buffer);
		alSourcei(source, AL_LOOPING, AL_FALSE);
		
		alSourcePlay(source);
		ALint s = AL_PLAYING;
		while (s == AL_PLAYING)
		{
			alGetSourcei(source, AL_SOURCE_STATE, &s);
			Sleep(10);   //50 ms
						 //alSourceStop(source);
		}
		alDeleteSources(1, &source);
		alDeleteBuffers(1, &buffer);
		alcMakeContextCurrent(NULL);
		alcDestroyContext(context);
		alcCloseDevice(device);
		
	}
}