/*
* Filename: UpixTimer.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_TIME_COUNTER_H_
#define _UPIX_TIME_COUNTER_H_

#include "UpixAlgorithmTypeDef.h"
#include "UpixDebug.h"


//----------------------------------------------
// 高精度时钟
class UPIX_ALGORITHM_CLASS UpixTimer
{
public:
	// constructor & destructor
	UpixTimer();

	virtual ~UpixTimer() {}

public:
	void StartTimer();

	void ResetTimer();

	// Get elapsed time from the beginning (unit: ms)
	double GetElapsedTime();

	// Get elapsed time from last getting duration
	double GetElapsedTimeFromLast();

	void LogElapsedTime(const char *caption, const char *tag = NULL);

	void LogElapsedTimeFromLast(const char *caption, const char *tag = NULL);

protected:
	bool GetTime(double *pTime);

protected:
	double m_InvFrequency;
	double m_StartTime;
	double m_CurrentTime;
};

// 
class UPIX_ALGORITHM_CLASS UpixAutoTimer : public UpixTimer
{
public:
	UpixAutoTimer(const char *caption = "UpixAutoTimer: ", const char *tag = NULL);
	virtual ~UpixAutoTimer();

private:
	const char *m_caption;
	const char *m_tag;
};



//
class UPIX_ALGORITHM_CLASS UpixTimeCounter
{
public:
	UpixTimeCounter() : m_nStart(0), m_nCurrent(0) {}
	virtual ~UpixTimeCounter() {}

public:
	void StartTimer();
	void ResetTimer();
	double GetElapsedTime();
	double GetElapsedTimeFromLast();
	void LogElapsedTime(const char *caption, const char *tag = NULL);
	void LogElapsedTimeFromLast(const char *caption, const char *tag = NULL);

private:
	ulong GetTime();

protected:
	ulong m_nStart;
	ulong m_nCurrent;
};

// 
class UPIX_ALGORITHM_CLASS UpixAutoTimeCounter : public UpixTimeCounter
{
public:
	UpixAutoTimeCounter(const char *caption = "UpixAutoTimer: ", const char *tag = NULL);
	virtual ~UpixAutoTimeCounter();

private:
	const char *m_caption;
	const char *m_tag;
};


// Debug timer API
#define UpixDebugStartTimer(debug, t) \
if (debug) \
{ \
	t.StartTimer(); \
}

#define UpixDebugResetTimer(debug, t) \
if (debug) \
{ \
	t.ResetTimer(); \
}

#define UpixDebugLogElapsedTime(debug, t, caption, tag) \
if (debug) \
{ \
	t.LogElapsedTime(caption, tag); \
}

#define UpixDebugGetElapsedTime(debug, t, duration) \
if (debug) \
{ \
	duration = t.GetElapsedTime();  \
}



// Get current time
struct tm;

class UPIX_ALGORITHM_CLASS UpixLoacalTime
{
public:
	UpixLoacalTime();
	virtual ~UpixLoacalTime();

	void GetTime();
	int Year();
	int Month();
	int MonthDay();
	int WeekDay();
	int Hour();
	int Minute();
	int Second();

protected:
	tm *m_pTime;
};


#endif