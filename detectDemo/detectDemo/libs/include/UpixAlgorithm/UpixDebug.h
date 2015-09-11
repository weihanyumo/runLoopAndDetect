/*
* Filename: UpixDebug.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_DEBUG_H_
#define _UPIX_DEBUG_H_

#include "UpixAlgorithmTypeDef.h"


enum eLogPriority
{
	UPIX_LOG_DEBUG = 3,
	UPIX_LOG_INFO,
	UPIX_LOG_WARN,
	UPIX_LOG_ERROR,
};


class UPIX_ALGORITHM_CLASS UpixDebug
{
public:
	UpixDebug() {}
	~UpixDebug() {}

public:
	// 
	static void DebugLog(int prio, const char *tag, const char *format, ...);

	// 返回当前进程的内存使用量
	static ulong GetCurProcessMemoryUsage();
};



#define UpixLog(format, ...)		{ UpixDebug::DebugLog(UPIX_LOG_DEBUG, "UpixAlgorithm", format, ##__VA_ARGS__); }
#define UpixLogD(tag, format, ...)	{ UpixDebug::DebugLog(UPIX_LOG_DEBUG, tag, format, ##__VA_ARGS__); }
#define UpixLogI(tag, format, ...)	{ UpixDebug::DebugLog(UPIX_LOG_INFO, tag, format, ##__VA_ARGS__); }
#define UpixLogW(tag, format, ...)	{ UpixDebug::DebugLog(UPIX_LOG_WARN, tag, format, ##__VA_ARGS__); }
#define UpixLogE(tag, format, ...)	{ UpixDebug::DebugLog(UPIX_LOG_ERROR, tag, format, ##__VA_ARGS__); }


#define UpixDebugLog(debug, tag, format, ...) \
if (debug) \
{ \
	UpixLogD(tag, format, ##__VA_ARGS__); \
}


#define UpixAssert(condition) \
if(!(condition)) \
{ \
	UpixLog("UpixAssert: File: %s, Line: %s, Func: %s", __FILE__, __LINE__, __FUNCTION__); \
}


#endif
