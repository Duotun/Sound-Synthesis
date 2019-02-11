#pragma once
#ifndef  SOUND_H
#define SOUND_H

#include <math.h>

namespace sound {
#define two_PI 6.283185307

	//a single sine wave changes with time (one sound)
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


}
#endif //
