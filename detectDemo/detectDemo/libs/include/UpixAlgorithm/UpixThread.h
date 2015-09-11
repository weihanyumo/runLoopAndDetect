/*
* Filename: UpixThread.h
* Author:   Xiehong
* Date:     2012.11.10
*/

#ifndef _UPIX_THREAD_H_
#define _UPIX_THREAD_H_


#include "UpixAlgorithmTypeDef.h"
#include "UpixMessageCenter.h"
#include "UpixException.h"
#include "UpixDebug.h"



#define UPIX_WAIT_INFINITE	-1

#define UPIX_SIGSPEND		 1
#define UPIX_SIGRESUME		 2

enum UpixThreadSatuts
{
	UPIX_THREAD_CANCEL = -1,
	UPIX_THREAD_FAILED,
	UPIX_THREAD_SUCCESS,
	UPIX_THREAD_USER = 1024
};


#ifdef _WIN32
typedef struct pthread_attr_t_ * UpixThreadAttr;
typedef struct pthread_mutex_t_ * UpixThreadMutex;
#else
typedef void * UpixThreadAttr;
typedef void * UpixThreadMutex;
#endif

typedef void * UpixHandleThread;



// Exceptions for thread
class UPIX_ALGORITHM_CLASS UpixThreadException : public UpixException
{
public:
	UpixThreadException(const char *strMsg) : UpixException(strMsg) {}
	UpixThreadException(int _code, const char *_err, const char *_func, const char *_file, int _line)
		: UpixException(_code, _err, _func, _file, _line)
	{

	}
};



// This class is thread safe, one object could start many threads.
class UPIX_ALGORITHM_CLASS IUpixThread : public IUpixMessageCenter
{
public:
	IUpixThread(IUpixMessageReceiver *pReceiver = NULL);
	virtual ~IUpixThread();

	// Start the thread
	bool Start(void *pParam);

	// Suspend the thread
	bool Suspend();

	// Resume the suspended thread
	bool Resume();

	// Stop the thread normally
	void Stop();

	//return the status of the thread
	bool IsStop();

	//wait for end
	UpixThreadSatuts WaitForEnd();

	//detach the thread 
	bool Detach();

	// Check thread if is active
	bool IsActive();

	// Get the thread id
	UpixHandleThread GetThreadID();

	// Sleep a time
	void Sleep(int nMilliSeconds);

	// Cancel a specified thread
	bool Cancel();

	//Kill a thread
	bool Kill();

	// The static thread function
	friend void* _ThreadProc(void *pData);

protected:
	// Thread callbacks
	// The really thread function, should be override
	virtual UpixThreadSatuts ThreadProc(void *pParam) = 0;

	// Set cancel type
	bool EnableCancel(bool Asychronous);
	bool DisableCancel();

	// Exit the thread with input data that returned by WaitForEnd
	void Exit(UpixThreadSatuts nRetVal);

	// Call back function for thread end
	virtual void OnEndThread();

	// Call back function for thread exception
	virtual void OnException(const UpixException *pExcp);

	// Call back function for cancel thread
	virtual void OnCancel();

	//Call back function for thread be killed
	virtual void OnKill();

	// Call back function for exit thread
	virtual void OnExit();

	//Call back function for Suspend
	virtual void OnSuspend();

	//Cal back function for Resume
	virtual void OnResume();

protected:
	void *m_pThreadParam;

private:
	bool m_bIsSuspend;
	bool m_bStop;
	bool m_bIsDetach;
	UpixHandleThread m_pThreadID;
	UpixThreadAttr m_pAttr;
};



// Class of lock
class UPIX_ALGORITHM_CLASS UpixLock
{ 
public: 
	UpixLock();
	virtual ~UpixLock();

	// Lock the resource. 
	void Lock();

	// Unlock the resource. 
	void UnLock();

protected: 
	UpixThreadMutex m_pMutex;
}; 


#endif