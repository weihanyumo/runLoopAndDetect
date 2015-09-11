/* class UpixLeastSquareNonlinear: base class for Levenberg-Marquardt optimization.

   This class is to solve least square problem using nonlinear method.
   Suppose {(x,y)} is a pair of measured sample, f = fun(x|p) is an estimation of y, with 
   parameter p, we have cost function such as:
	   e = ||f-y||^2 = ||fun(x|p)-y||^2
   Our purpose is to minimize the residual e of cost function, and find the right parameters 
   in mapping function fun() synchronously. Here we call measured y reference.

   In the abstract base class some basic options and methods were defined. But the mapping 
   function and Jacobian function and some other functions should be implemented in the 
   derived class.

   Class members:
    m_opts:      some options the Levenberg-Marquardt method may use.
    m_info:		 informations in computing and informations of result.
    m_param:	 parameters in mapping function.
    m_refer:	 reference data, just y which be measured
    m_szRefer:	 size of reference data.
    m_szParam:	 size of parameter.

   Override functions:
    constructor:       in constructor we should set measured data x, y and other constant 
					   data in use.
	destructor:		   release data of derive class.
	SetOptimizeParam:  parameterize, set structured parameters to m_param as an initial guess.
	GetOptimizeParam:  parse parameters, get structured parameters from m_param as a result.
	SetReferData:	   parameterize, set structured reference y to m_refer.
	fun:			   implement f = fun(x|p), f has the same size of y.
	Jacobian:		   implement jac = Jacobian(x|p), return Jacobian matrix.

   Example:
    class CLsqnl : public IUpixLeastSquareNonlinear; // Implemented
	class CParam : public IUpixOptimizeParams;		 // Implemented

	double a = a0;									 // Initial guess
	double b = b0;									 // Initial guess
	CParam p(a, b);									 // Construct parameters
	CLsqnl obj(double x[], double y[], void *pdata); // Construct LM optimizer
	obj.SetOptimizeParam(&p);
	obj.Optimize(100, true);						 // Optimize with Jacobian
	obj.GetOptimizeParam(&p);
	cout<<p.a<<" "<<p.b<<endl;						 // Print result

   Author: Xie Hong
   Data:   2010.9.25
*/

#ifndef _UPIX_LEAST_SQUARE_NONLINEAR_H_
#define _UPIX_LEAST_SQUARE_NONLINEAR_H_


#include "UpixAlgorithmTypeDef.h"




// Parameters to optimize. Each implement class that derived
// from UpixLeastSquareNonlinear should use an implement class 
// derived from this abstract interface.
class IUpixOptimizeParams
{
public:
	IUpixOptimizeParams() {}
	virtual ~IUpixOptimizeParams() {}
};



const int UPIX_LM_OPT_SIZE = 5;
const int UPIX_LM_INFO_SIZE = 10;


// The abstract interface for LM optimization
class UPIX_ALGORITHM_CLASS IUpixLeastSquareNonlinear
{
public:
	IUpixLeastSquareNonlinear(int param_num, int refer_num);
	virtual ~IUpixLeastSquareNonlinear();

public:
	// Parameterize to *param. Set initial guess to optimizer
	virtual void SetOptimizeParams(IUpixOptimizeParams *params) = 0;

	// Parse parameters from *param. Get the result of optimization.
	virtual void GetOptimizeParams(IUpixOptimizeParams *params) = 0;

	// Optimize the least square function
	double Optimize(int itmax = 100, bool bJacobian = false, double *work = NULL, double *covar = NULL);

	// Report result status
	void LogOptimizeStatus();

	// Get the size of optimize options
	int GetOptionSize();

	// Get the size of optimize information
	int GetInfoSize();

	// get optimize information
	double *GetInfo() { return m_info; }

	// get RMS error
	double GetResidual();

	// set optimize options
	void SetOptions(double optsValue, int optsIndex);

protected:
	// Parameterize to *refer
	virtual void SetReferData() = 0;

	// Mapping function: y = f(x|p). For measured y, e = ||f(x|p) - y||^2
	virtual void fun(double *param, double *f, int param_num, int refer_num, void *pdata) = 0;

	// Jacobian function
	virtual void Jacobian(double *param, double *jac, int param_num, int refer_num, void *pdata) {}

private:
	// Global cost function as the first parameter of levmar(). Make fun() could be passed into levmar().
	friend void gfun(double *param, double *f, int param_num, int refer_num, void *pdata);

	// Global Jacobian function
	friend void gJacfun(double *param, double *jac, int param_num, int refer_num, void *pdata);

	// Optimize by fun() that passed in, without Jacobian function
	int Optimize(void (*func)(double *param, double *f, int param_num, int refer_num, void *pdata), 
				int itmax = 100, void *pdata = NULL, double *work = NULL, double *covar = NULL);

	// Optimize by fun() that passed in, with Jacobian function
	int Optimize(void (*func)(double *param, double *f, int param_num, int refer_num, void *pdata), 
				void (*jacf)(double *param, double *jac, int param_num, int refer_num, void *pdata),
				int itmax = 100, void *pdata = NULL, double *work = NULL, double *covar = NULL);

	// allocate buffer
	void Allocate(int param_num, int refer_num);

	// release buffer
	void Release();

protected:
	double m_opts[UPIX_LM_OPT_SIZE];	// LM_OPTS_SZ = 5
	double m_info[UPIX_LM_INFO_SIZE];	// LM_INFO_SZ = 10
	double *m_param;
	double *m_refer; // just y which be measured
	int m_szRefer;
	int m_szParam;
};


#endif