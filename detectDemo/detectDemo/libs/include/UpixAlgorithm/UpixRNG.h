#ifndef _UPIX_RNG_H_
#define _UPIX_RNG_H_

#include "UpixAlgorithmDllExport.h"
#include "UpixAlgorithmTypeDef.h"



class UPIX_ALGORITHM_CLASS UpixRNG
{
private:
	enum { UNIFORM=0, NORMAL=1 };

public:
	UpixRNG();
	UpixRNG(uint64 state);
	~UpixRNG();

public:
	// Get next number from random series
	unsigned Next();

	// Generate uniform random
	int Uniform(int a, int b);
	float Uniform(float a, float b);
	double Uniform(double a, double b);

	// Generate Gaussian random
	void Randn_0_1_32f( float *arr, int len, uint64 *state );
	double Gaussian(double sigma);

	operator uchar();
	operator ushort();
	operator short();
	operator unsigned();
	operator int();
	operator float();
	operator double();

	unsigned operator ()(unsigned N);
	unsigned operator ()();

private:
	void CreateGaussianTable();

public:
	uint64 state;

private:
	// For Gaussian
	static unsigned kn[128];
	static float wn[128];
	static float fn[128];
	static bool initialized;
};



#endif