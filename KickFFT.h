/*
 FILENAME:	KickFFT.cpp
 AUTHOR:	Orlando S. Hoilett and Akio K. Fujita
 EMAIL:     orlandohoilett@gmail.com
 VERSION:	3.0.0
 
 
 AFFILIATIONS
 Linnes Lab, Weldon School of Biomedical Engineering,
 Purdue University, West Lafayette, IN 47907
 
 
 DESCRIPTION
 This is a templated, static class so funciton calls must be preceeded with
 "KickFFT<variable_type>::" where variable_type should be replaced with int16_t,
 int, float, etc.
 
 This class implements a Discrete Fourier Transform using look-up tables instead
 floating point math in order to optimize computation time. The library was
 developed to help refactor the KickDSP class.
 
 This class also computes the Power Spectral Density (PSD) function.
 
 This class also has a depedency on the TrigDef and KickMath classes.
 
 
 UPDATES
 Version 1.0.0
 2020/02/19:1200>
 			- Initial compilation.
 Version 1.0.1
 2020/03/23:1332>
			- Updated descriptions and comments.
			- Added future updates section.
 2020/06/28:1426>
 			- Updated comments a bit.
 Version 1.1.0
 2020/07/11:0655>
			- Added power spectral density (PSD) function.
			- Updated comments.
 Version 2.0.0
 2020/08/21:1555>
 			- Removed U2I function.
 			- Moved to templated class.
 Version 3.0.0
 2020/08/23:0406> (UTC-5)
 			- changed magnitude types to uint32_t to match isqrt function
 
 
 FUTURE UPDATES TO INCLUDE
 1. Refactor the code a bit to avoid redundancy
 2. Remove the need for multiple lookup tables of different sizes
 3. Expand to 1024 samples
 (CHECK) 4. Open-sourcing and moving to GitHub instead of the lab's private BitBucket.
 (CHECK) 5. Making this a templated class, meaning it will accept any data type for the
 data to be filtered.
 6. Implements inverse fourier transform
 (CHECK) 7. Consider making isqrt function to a uint32_t
 
 
 
 DISCLAIMER
 Linnes Lab code, firmware, and software is released under the
 MIT License (http://opensource.org/licenses/MIT).
 
 The MIT License (MIT)
 
 Copyright (c) 2020 Linnes Lab, Purdue University
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 of the Software, and to permit persons to whom the Software is furnished to do
 so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 
 */


#ifndef KickFFT_h
#define KickFFT_h


//Standard Arduino libraries
#include <Arduino.h>

//Custom Kick libraries
#include "KickMath.h"

//Custom Internal Libraries
#include "TrigDef.h"


template<typename Type>


class KickFFT
{
	
public:
	
	static void fft(uint16_t samples, const Type data[], uint32_t mag[]);
		
	static void fft(float fs, float f1, float f2, uint16_t samples, const Type data[], uint32_t mag[]);
	
	static void fft(float fs, float f1, float f2, uint16_t samples, const Type data[],
					uint32_t mag[], uint16_t &startIndex, uint16_t &endIndex);
	
	static void psd(float fs, float f1, float f2, uint16_t samples, const Type data[], uint32_t mag[]);
	
	static void psd(float fs, float f1, float f2, uint16_t samples, const Type data[],
					uint32_t mag[], uint16_t &startIndex, uint16_t &endIndex);

//	static void ifft(float fs, uint16_t samples, const Type mag[], Type output[]);

};



//void KickFFT<Type>::fft(uint16_t samples, const Type data[], int32_t mag[])
//samples	number of samples in input data array
//data		input data array
//mag		array to store calculated frequency magnitdues
//
//This method computes the Discrete Fourier Transform using the definition described
//here <https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definition>
//
//This method computes the entire DFT across 0 to Fs/2.
template<typename Type>
void KickFFT<Type>::fft(uint16_t samples, const Type data[], uint32_t mag[])
{
	uint16_t startIndex = 0;
	uint16_t endIndex = samples;
	
	for (uint16_t i = startIndex; i < endIndex; i++)
	{
		signed long int real = 0;
		signed long int imag = 0;
		
		
		//Euler's Identity
		for (uint16_t j = 0; j < samples; j++)
		{
			//uses lookup tables for trigonometric
			//functions to save compputing power
			switch(samples)
			{
				case 512:
					real += intcosine512[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine512[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 256:
					real += intcosine256[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine256[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 128:
					real += intcosine128[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine128[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 64:
					real += intcosine64[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine64[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 32:
					real += intcosine32[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine32[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				default:
					break;
			}
		}
		
		
		//dividing each number by 1000 to prevent each number from getting too large
		//Also adjusts for the fact that the trigonometric values were multiplied
		//by 1000 as well to make them integers instead of decimal values
		real = real/1000;
		imag = imag/1000;
		
		
		//calculating magnitude of the data by taking the square root of the
		//sum of the squares of the real and imaginary component of each signal
		mag[i] = KickMath<signed long int>::calcMagnitude(real, imag);
	}
}


//void KickFFT<Type>::fft(float fs, float f1, float f2, uint16_t samples, const Type data[], int32_t mag[])
//fs		sampling frequency for input data array Hertz (Hz)
//f1		low frequency bound to calculate frequency spectrum
//				- useful for limiting range and decreasing computation time
//				- set to 0 if no specific range is needed
//f2		high frequency bound to calculate frequency spectrum
//				- useful for limiting range and decreasing computation time
//				- set to fs if no specific range is needed
//samples	number of samples in input data array
//data		input data array
//mag		array to store calculated frequency magnitdues
//
//
//This method computes the Discrete Fourier Transform using the definition described
//here <https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definition>
//
//Computes the DFT across the range f1 to f2
template<typename Type>
void KickFFT<Type>::fft(float fs, float f1, float f2, uint16_t samples, const Type data[], uint32_t mag[])
{
	//changes f1 and f2 to indices
	//fs/samples gives the increments of frequency on the x-axis
	uint16_t startIndex = f1/(fs/samples);
	uint16_t endIndex = f2/(fs/samples);
	
	for (uint16_t i = startIndex; i < endIndex; i++)
	{
		signed long int real = 0;
		signed long int imag = 0;
		
		
		//Euler's Identity
		for (uint16_t j = 0; j < samples; j++)
		{
			//uses lookup tables for trigonometric
			//functions to save compputing power
			switch(samples)
			{
				case 512:
					real += intcosine512[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine512[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 256:
					real += intcosine256[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine256[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 128:
					real += intcosine128[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine128[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 64:
					real += intcosine64[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine64[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 32:
					real += intcosine32[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine32[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				default:
					break;
			}
		}
		
		
		//dividing each number by 1000 to prevent each number from getting too large
		//Also adjusts for the fact that the trigonometric values were multiplied
		//by 1000 as well to make them integers instead of decimal values
		real = real/1000;
		imag = imag/1000;
		
		
		//calculating magnitude of the data by taking the square root of the
		//sum of the squares of the real and imaginary component of each signal
		mag[i] = KickMath<signed long int>::calcMagnitude(real, imag);
	}
}


//void KickFFT<Type>::fft(float fs, float f1, float f2, uint16_t samples, const Type data[],
//						  int32_t mag[], uint16_t &startIndex, uint16_t &endIndex)
//
//fs			sampling frequency for input data array Hertz (Hz)
//f1			low frequency bound to calculate frequency spectrum
//					- useful for limiting range and decreasing computation time
//					- set to 0 if no specific range is needed
//f2			high frequency bound to calculate frequency spectrum
//					- useful for limiting range and decreasing computation time
//					- set to fs if no specific range is needed
//samples		number of samples in input data array
//data			input data array
//mag			array to store calculated frequency magnitdues
//startIndex	returns the index of the array that corresponds to f1
//endIndex		returns the index of the array that corresponds to f2
//
//
//This method computes the Fast Fourier Transform using the definition described
//here <https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definition>
template<typename Type>
void KickFFT<Type>::fft(float fs, float f1, float f2, uint16_t samples, const Type data[],
						uint32_t mag[], uint16_t &startIndex, uint16_t &endIndex)
{
	//changes f1 and f2 to indices
	//fs/samples gives the increments of frequency on the x-axis
	startIndex = f1/(fs/samples);
	endIndex = f2/(fs/samples);
	
	
	for (uint16_t i = startIndex; i < endIndex; i++)
	{
		signed long int real = 0;
		signed long int imag = 0;
		
		
		//Euler's Identity
		for (uint16_t j = 0; j < samples; j++)
		{
			//uses lookup tables for trigonometric
			//functions to save compputing power
			switch(samples)
			{
				case 512:
					real += intcosine512[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine512[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 256:
					real += intcosine256[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine256[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 128:
					real += intcosine128[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine128[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 64:
					real += intcosine64[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine64[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 32:
					real += intcosine32[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine32[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				default:
					break;
			}
		}
		
		
		//dividing each number by 1000 to prevent each number from getting too large
		//Also adjusts for the fact that the trigonometric values were multiplied
		//by 1000 as well to make them integers instead of decimal values
		real = real/1000;
		imag = imag/1000;
		
		
		//calculating magnitude of the data by taking the square root of the
		//sum of the squares of the real and imaginary component of each signal
		mag[i] = KickMath<signed long int>::calcMagnitude(real, imag);
	}
}


//void KickFFT<Type>::psd(float fs, float f1, float f2, uint16_t samples, const Type data[], int32_t mag[])
//fs		sampling frequency for input data array Hertz (Hz)
//f1		low frequency bound to calculate frequency spectrum
//				- useful for limiting range and decreasing computation time
//				- set to 0 if no specific range is needed
//f2		high frequency bound to calculate frequency spectrum
//				- useful for limiting range and decreasing computation time
//				- set to fs if no specific range is needed
//samples	number of samples in input data array
//data		input data array
//mag		array to store calculated frequency magnitdues
//
//
//This method computes the power spectral density of the signal. This is defined
//as the square of the absoluate values of the magnitude of the computed FFT.
//
//We intentionally do not call the fft function to save on processing power by
//avoiding having to scan through the arrays multiple times in the fft method
//and then again for the psd method.
template<typename Type>
void KickFFT<Type>::psd(float fs, float f1, float f2, uint16_t samples, const Type data[], uint32_t mag[])
{
	//changes f1 and f2 to indices
	//fs/samples gives the increments of frequency on the x-axis
	uint16_t startIndex = f1/(fs/samples);
	uint16_t endIndex = f2/(fs/samples);
	
	for (uint16_t i = startIndex; i < endIndex; i++)
	{
		signed long int real = 0;
		signed long int imag = 0;
		
		
		//Euler's Identity
		for (uint16_t j = 0; j < samples; j++)
		{
			//uses lookup tables for trigonometric
			//functions to save compputing power
			switch(samples)
			{
				case 512:
					real += intcosine512[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine512[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 256:
					real += intcosine256[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine256[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 128:
					real += intcosine128[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine128[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 64:
					real += intcosine64[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine64[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 32:
					real += intcosine32[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine32[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				default:
					break;
			}
		}
		
		
		//dividing each number by 1000 to prevent each number from getting too large
		//Also adjusts for the fact that the trigonometric values were multiplied
		//by 1000 as well to make them integers instead of decimal values
		real = real/1000;
		imag = imag/1000;
		
		
		//calculating magnitude of the data by taking the square root of the
		//sum of the squares of the real and imaginary component of each signal
		mag[i] = KickMath<signed long int>::calcMagnitude(real, imag);
		//calculates the PSD by squaring the magnitudes
		mag[i] = mag[i]*mag[i];
	}
}


//void KickFFT<Type>::psd(float fs, float f1, float f2, uint16_t samples, const Type data[],
//						  int32_t mag[], uint16_t &startIndex, uint16_t &endIndex)
//
//fs		sampling frequency for input data array Hertz (Hz)
//f1		low frequency bound to calculate frequency spectrum
//				- useful for limiting range and decreasing computation time
//				- set to 0 if no specific range is needed
//f2		high frequency bound to calculate frequency spectrum
//				- useful for limiting range and decreasing computation time
//				- set to fs if no specific range is needed
//startIndex	returns the index of the array that corresponds to f1
//endIndex		returns the index of the array that corresponds to f2
//
//
//This method computes the power spectral density of the signal. This is defined
//as the square of the absoluate values of the magnitude of the computed FFT.
//
//We intentionally do not call the fft function to save on processing power by
//avoiding having to scan through the arrays multiple times in the fft method
//and then again for the psd method.
template<typename Type>
void KickFFT<Type>::psd(float fs, float f1, float f2, uint16_t samples, const Type data[],
						uint32_t mag[], uint16_t &startIndex, uint16_t &endIndex)
{
	//changes f1 and f2 to indices
	//fs/samples gives the increments of frequency on the x-axis
	startIndex = f1/(fs/samples);
	endIndex = f2/(fs/samples);
	
	
	for (uint16_t i = startIndex; i < endIndex; i++)
	{
		signed long int real = 0;
		signed long int imag = 0;
		
		
		//Euler's Identity
		for (uint16_t j = 0; j < samples; j++)
		{
			//uses lookup tables for trigonometric
			//functions to save compputing power
			switch(samples)
			{
				case 512:
					real += intcosine512[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine512[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 256:
					real += intcosine256[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine256[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 128:
					real += intcosine128[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine128[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 64:
					real += intcosine64[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine64[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				case 32:
					real += intcosine32[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					imag += intsine32[ (i*j) - (samples*((i*j)/samples)) ] * data[j];
					break;
					
				default:
					break;
			}
		}
		
		
		//dividing each number by 1000 to prevent each number from getting too large
		//Also adjusts for the fact that the trigonometric values were multiplied
		//by 1000 as well to make them integers instead of decimal values
		real = real/1000;
		imag = imag/1000;
		
		
		//calculating magnitude of the data by taking the square root of the
		//sum of the squares of the real and imaginary component of each signal
		mag[i] = KickMath<signed long int>::calcMagnitude(real, imag);
		//calculates the PSD by squaring the magnitudes
		mag[i] = mag[i]*mag[i];
	}
}


//template<typename Type>
//void KickFFT<Type>::ifft(float fs, uint16_t samples, const Type mag[], Type output[])
//{
//	Type maxMagnitude = mag[0];
//	for(uint16_t i = 1; i < samples; i++)
//	{
//		if(mag[i] > maxMagnitude) maxMagnitude = mag[i];
//	}
//
//
//	for(uint16_t i = 0; i < samples; i++)
//	{
//		output[i] = 0;
//	}
//
//
//	for(uint16_t i = 0; i < samples; i++)
//	{
//		for(uint16_t j = 0; j < samples/2; j++)
//		{
//			output[i] += mag[j]*sin(2*PI*(j*(fs/samples))*(i*1/fs));
//
////			Serial.print(i);
////			Serial.print(",");
////			Serial.print(j);
////			Serial.println();
//		}
//	}
//}


#endif /* KickFFT_h */

