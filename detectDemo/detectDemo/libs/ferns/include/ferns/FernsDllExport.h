/*
* Filename: FernsDllExport.h
* Author:   Xiehong
* Date:     2013.2.7
*/

#ifndef _FERNS_DLL_EXPORT_H_
#define _FERNS_DLL_EXPORT_H_


#ifdef _WIN32
#ifndef FERNS_EXPORTS
	#define FERNS_API		__declspec(dllimport)
	#define FERNS_CLASS		__declspec(dllimport)
	#define FERNS_TEMPLATE
#else //EXPORT
	#define FERNS_API		__declspec(dllexport)
	#define FERNS_CLASS		__declspec(dllexport)
	#define FERNS_TEMPLATE	__declspec(dllexport)
#endif
#else
	#define FERNS_API
	#define FERNS_CLASS
	#define FERNS_TEMPLATE
#endif


#endif