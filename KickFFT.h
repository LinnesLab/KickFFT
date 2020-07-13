/*
 FILENAME:	KickFFT.h
 AUTHOR:	Orlando S. Hoilett and Akio K. Fujita
 EMAIL:     orlandohoilett@gmail.com
 
 
 Please see .cpp file for extended descriptions, instructions, and version updates
 
 
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


class KickFFT
{
	
public:
	
	static void fft(uint16_t samples, const int16_t data[], uint16_t mag[]);
		
	static void fft(float fs, float f1, float f2, uint16_t samples, const int16_t data[], int32_t mag[]);
	
	static void fft(float fs, float f1, float f2, uint16_t samples, const int16_t data[],
					uint16_t mag[], uint16_t &startIndex, uint16_t &endIndex);
					
	static void U2I(uint16_t unsignedArray[], int16_t signedArray[], uint16_t samples);
	
	static void psd(float fs, float f1, float f2, uint16_t samples, const int16_t data[], int32_t mag[]);
	
	static void psd(float fs, float f1, float f2, uint16_t samples, const int16_t data[],
					int32_t mag[], uint16_t &startIndex, uint16_t &endIndex);

	//static void ifft(float fs, float f1, float f2, uint16_t samples, int16_t data[], const uint16_t mag[]);

};


#endif /* KickFFT_h */

