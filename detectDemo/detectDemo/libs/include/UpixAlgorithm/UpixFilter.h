/*
* Filename: UpixFilter.h
* Author:   Xiehong
* Date:     2013.9.5
*/

#ifndef _UPIX_FILTER_H_
#define _UPIX_FILTER_H_


#include "UpixVector.h"
#include "UpixArray.h"
#include "UpixMatrix.h"
#include "UpixPoint.h"



// Alpha-beta filter
class UPIX_ALGORITHM_TEMPLATE UpixAlphaBetaFilter
{
public:
	UpixAlphaBetaFilter();
	virtual ~UpixAlphaBetaFilter();

public:
	void Init(int nParamCount, float fDeltaTime, float alpha = 0.7f, float beta = 0.05f);

	void SetInitGuess(const UpixArrayf &vInitX, const UpixArrayf &vInitV = UpixArrayf());
	void SetInitGuess(float vInitX[], int nSize, float vInitV[] = NULL);

	void Predict(UpixArrayf &vPredicts);
	void Predict(float vPredicts[], int nSize);
	void Predict();

	void Correct(const UpixArrayf &vMeasures, UpixArrayf &vCorrects);
	void Correct(float vMeasures[], float vCorrects[], int nSize);

private:
	UpixVector<float> m_Xpre;
	UpixVector<float> m_Vpre;
	UpixVector<float> m_X;
	UpixVector<float> m_V;
	UpixVector<float> m_Xmeasure;
	UpixVector<float> m_Resid;
	float m_dT;
	float m_dTinv;
	float m_alpha;
	float m_beta;
	int m_nParamCount;
};




// Kalman filter
class UPIX_ALGORITHM_TEMPLATE UpixKalmanFilter
{
public:
	UpixKalmanFilter();
	virtual ~UpixKalmanFilter();

public:
	void Init(int dynamParams, int measureParams, int controlParams, 
			const UpixMatrix<float> &matTrans, const UpixMatrix<float> &matProcNoiseCov, 
			const UpixMatrix<float> &matControl = UpixMatrix<float>(), 
			const UpixMatrix<float> &matMeasure = UpixMatrix<float>(), 
			const UpixMatrix<float> &matMeasureNoiseCov = UpixMatrix<float>());
	void Init(int dynamParams, int measureParams, int controlParams, 
			float matTrans[], float matProcNoiseCov[], float matControl[] = NULL,
			float matMeasure[] = NULL, float matMeasureNoiseCov[] = NULL);

	void SetInitGuess(const UpixArrayf &vInitState);
	void SetInitGuess(float vInitState[], int nSize);

	// Return the predicted data
	void Predict(UpixArrayf &vPredicts, const UpixVector<float> &vControl = UpixVector<float>());
	void Predict(float vPredicts[], int nSize, const UpixVector<float> &vControl = UpixVector<float>());
	void Predict(const UpixVector<float> &vControl = UpixVector<float>());

	// The measurements will be modified to corrected data
	void Correct(const UpixArrayf &vMeasures, UpixArrayf &vCorrects);
	void Correct(float vMeasures[], float vCorrects[], int nSize);

private:
	void Allocate(int dynamParams, int measureParams, int controlParams);

private:
	int m_nStateNum;
	int m_nMeasureNum;
	int m_nControlNum;

	// Status
	UpixVector<float> m_vMeasurements;
	UpixVector<float> m_statePre;           //!< predicted state (x'(k)): x(k)=A*x(k-1)+B*u(k)
	UpixVector<float> m_statePost;          //!< corrected state (x(k)): x(k)=x'(k)+K(k)*(z(k)-H*x'(k))

	// Filter matrix
	UpixMatrix<float> m_transitionMatrix;   //!< state transition matrix (A)
	UpixMatrix<float> m_controlMatrix;      //!< control matrix (B) (not used if there is no control)
	UpixMatrix<float> m_measurementMatrix;  //!< measurement matrix (H)
	UpixMatrix<float> m_processNoiseCov;    //!< process noise covariance matrix (Q)
	UpixMatrix<float> m_measurementNoiseCov;//!< measurement noise covariance matrix (R)

	// error matrix
	UpixMatrix<float> m_errorCovPre;        //!< priori error estimate covariance matrix (P'(k)): P'(k)=A*P(k-1)*At + Q)*/
	UpixMatrix<float> m_gain;               //!< Kalman m_gain matrix (K(k)): K(k)=P'(k)*Ht*inv(H*P'(k)*Ht+R)
	UpixMatrix<float> m_errorCovPost;       //!< posteriori error estimate covariance matrix (P(k)): P(k)=(I-K(k)*H)*P'(k)

	// temporary matrices
	UpixMatrix<float> temp1;
	UpixMatrix<float> temp2;
	UpixMatrix<float> temp3;
	UpixMatrix<float> temp4;
	UpixVector<float> temp5;
};




// Kalman for 1-dimension movement
class UPIX_ALGORITHM_CLASS UpixKalmanFilter1D
{
public:
	UpixKalmanFilter1D();
	virtual ~UpixKalmanFilter1D();

public:
	void Init(float fDeltaTime, float sigmaProc = 1.0f, float sigmaMeasure = 1.0f);

	void SetInitGuess(float x);

	float Predict();

	float Correct(float x);

private:
	UpixKalmanFilter *m_pKalman;
};




// Kalman for 2-dimension movement
class UPIX_ALGORITHM_CLASS UpixKalmanFilter2D
{
public:
	UpixKalmanFilter2D();
	virtual ~UpixKalmanFilter2D();

public:
	void Init(float fDeltaTime, float sigmaProc = 1.0f, float sigmaMeasure = 1.0f);

	void SetInitGuess(const UpixPoint2f &pt);
	void SetInitGuess(float pt[], int nSize);

	void Predict(UpixPoint2f &pt);
	void Predict(float pt[], int nSize);
	void Predict();

	void Correct(UpixPoint2f &pt);
	void Correct(float pt[], int nSize);

private:
	UpixKalmanFilter *m_pKalman;
};




#endif