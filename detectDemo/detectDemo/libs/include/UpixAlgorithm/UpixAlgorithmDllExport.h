/*
* Filename: UpixAlgorithmDllExport.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_ALGORITHM_DLL_EXPORT_H_
#define _UPIX_ALGORITHM_DLL_EXPORT_H_

#ifdef _WIN32
	#ifndef UPIX_ALGORITHM_EXPORTS
		#define UPIX_ALGORITHM_API		__declspec(dllimport)
		#define UPIX_ALGORITHM_CLASS	__declspec(dllimport)
		#define UPIX_ALGORITHM_TEMPLATE
	#else //EXPORT
		#define UPIX_ALGORITHM_API		__declspec(dllexport)
		#define UPIX_ALGORITHM_CLASS	__declspec(dllexport)
		#define UPIX_ALGORITHM_TEMPLATE __declspec(dllexport)
	#endif// UPIX_ALGORITHM_EXPORTS
#else//_WIN32
		#define UPIX_ALGORITHM_API
		#define UPIX_ALGORITHM_CLASS
		#define UPIX_ALGORITHM_TEMPLATE
#endif//_WIN32

#endif//_UPIX_ALGORITHM_DLL_EXPORT_H_
