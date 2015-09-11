/*
* Filename: UpixAlgorithmTypeDef.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_ALGORITHM_TYPE_DEF_H_
#define _UPIX_ALGORITHM_TYPE_DEF_H_


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "UpixAlgorithmDllExport.h"


#if defined WIN32 || defined _WIN32
#	define UPIX_CDECL __cdecl
#	define UPIX_STDCALL __stdcall
#else
#	define UPIX_CDECL
#	define UPIX_STDCALL
#endif

#ifdef __cplusplus
#   define UPIX_EXTERN_C extern "C"
#	define UPIX_EXTERN_C_FUNCPTR(x) extern "C" { typedef x; }
#else
#	define UPIX_EXTERN_C
#   define UPIX_EXTERN_C_FUNCPTR(x) typedef x
#endif

#if defined __cplusplus
#   define UPIX_INLINE inline
#elif (defined WIN32 || defined _WIN32) && !defined __GNUC__
#   define UPIX_INLINE __inline
#else
#   define UPIX_INLINE static
#endif

#define UPIX_EXPORTS extern

#define UPIX_API(rettype) UPIX_EXTERN_C rettype UPIX_CDECL



// Constant in math
#define UPIX_INT_MAX		2147483647
#define UPIX_INT_MIN		(-2147483647 - 1)
#define UPIX_FLOAT_EPSILON	1.192092896e-07F
#define UPIX_DOUBLE_EPSILON 2.2204460492503131e-016
#define UPIX_REAL_MAX		1.7976931348623158e+308
#define UPIX_REAL_MIN		2.2250738585072014e-308
#define UPIX_FLOAT_MAX		3.402823466e+38F
#define UPIX_FLOAT_MIN		1.175494351e-38F
#define UPIX_PI				3.14159265358979323846
#define UPIX_PI_4			0.78539816339744831
#define UPIX_E				2.71828182845904523536
#define UPIX_SQRT2			1.41421356237309504880
#define UPIX_SQRT2_2		0.7071067811865475244
#define UPIX_INV_LOG10_2	3.32192809488736234787


// Useful functions
#define UpixMin(a, b)	((a) < (b) ? (a) : (b))
#define UpixMax(a, b)	((a) > (b) ? (a) : (b))
#define UpixAbs(x)		((x) < 0 ? -(x) : (x))
#define UpixRound(x)	( (int)((x) > 0 ? ((x) + 0.5) : ((x) - 0.5)) )
#define UpixCeil(x)		( UpixAbs((x) - (int)(x)) > UPIX_FLOAT_EPSILON ? ((int)((x) + 1)) : ((int)(x)) )
#define UpixRandom(n)	( rand()%(n) )	// generate a random number in the range of [0,n)
#define UpixClamp(x, low, high)	( UpixMax((low), UpixMin((x), (high))) )

// Linear interpolate, x0*(1 - p) + x1*p
#define UpixInterp(x0, x1, p)	( (x0) + (p)*((x1) - (x0)) )



#define UPIX_IN
#define UPIX_OUT



typedef enum _UpixBool
{
	UPIX_FALSE = 0,
	UPIX_TRUE = 1
} UpixBool;


// Short type define
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned long long ullong;
typedef long long llong;


#if defined WIN32 || defined _WIN32
typedef unsigned __int64 uint64;
#else
typedef unsigned long long uint64;
#endif


#endif