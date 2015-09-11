#ifndef _PYR_YAPE_H_
#define _PYR_YAPE_H_


#include <stdlib.h>
#include "FernsDllExport.h"
#include "keypoint.h"


struct _IplImage;
typedef _IplImage IplImage;

class pyr_yape06;
class fine_gaussian_pyramid;


class FERNS_CLASS UpixPyrYape06
{
public:
	UpixPyrYape06(int nfeatures, int nlevels, int patchSize = 31);
	~UpixPyrYape06();

public:
	int Detect(IplImage *pGrayImage, IplImage *pMask = NULL);

	IplImage* DrawKeyPoints(IplImage *pGrayImage, keypoint *keypoints, int keyNum);

public:
	keypoint *keypoints;

private:
	pyr_yape06 *m_pDetector;
	fine_gaussian_pyramid *m_pPyramid;
	int m_patchSize;
	int m_nMaxNumOfKey;
};



class FERNS_CLASS UpixYape06
{
public:
	UpixYape06(int laplacian_threshold = 30, int min_eigenvalue_threshold = 25, int keyNum = 50000);
	~UpixYape06();

public:
	int Detect(IplImage *pGrayImage, IplImage *pMask = NULL);
	keypoint *GetKeyPoint();

//private:
	pyr_yape06 *m_pDetector;
};


#endif
