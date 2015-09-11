/*
* Filename: UpixKMeans.h
* Author:   DuHaodong
* Date:     2013.10.11
*/
#ifndef _UPIX_KMEANS_H_
#define _UPIX_KMEANS_H_


#include "UpixMatrix.h"
#include "UpixArray.h"
#include "UpixPoint.h"
#include "UpixDebug.h"
#include "UpixRNG.h"



//kmeans parameters
class UPIX_ALGORITHM_CLASS UpixKmeansParams
{
public:
	enum { MAX_ITER = 1, EPS = 2 };

	enum { RANDOM_CENTERS=0, PP_CENTERS=1, USE_INITIAL_LABELS=2 };

	UpixKmeansParams() : attempts(3), flags(RANDOM_CENTERS), nType(0), nMaxIter(0), dEpsilon(0) {}

	UpixKmeansParams(int _attempts, int _flags, int _type, int _maxIter, double _epsilon)
		: attempts(_attempts), flags(_flags), nType(_type), nMaxIter(_maxIter), dEpsilon(_epsilon) {}

public:
	int attempts;		// Attempt times in K-Means
	int flags;			// Method of center initialization
	int nType;			// Iterator conditions
	int nMaxIter;		// Max iteration times
	double dEpsilon;	// Min iteration error
};



//kmeans
class UPIX_ALGORITHM_CLASS UpixKMeans
{
public:
	double KMeans(const UpixMatrix<float> &matData, int K, UpixArray<int> &arrayBestLabels, UpixKmeansParams &params, UpixMatrix<float> &_matCenters);

private:
	void UpixKMeansDistanceComput(double *pDistances, int *pLabels, const UpixMatrix<float> &matData, UpixMatrix<float> &matCenters);

	void UpixKMeansPPDistanceComput(float *tdist2, const float *pData, const float *dist, const int dims, const uint step, const uint stepci, int N);

	void UpixGenerateRandomCenter(const UpixMatrix<float> &matBox, float* pCenter, int Num, UpixRNG &rng);

	void UpixGenerateCentersPP(const UpixMatrix<float>& _matData, UpixMatrix<float>& _matOut_centers, UpixRNG &rng, int K, int trials);

};



#endif