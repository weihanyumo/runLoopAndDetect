/*
* Filename: UpixMemory.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_MEMORY_H_
#define _UPIX_MEMORY_H_


#include "UpixAlgorithmTypeDef.h"
#include <memory.h>


template <class type, class ValType>
inline void UPIX_ALGORITHM_TEMPLATE Fill(type *p, int n, const ValType &val)
{
	for (int i=0; i<n; i++)
	{
		p[i] = val;
	}
}

template <class type, class ValType>
inline void UPIX_ALGORITHM_TEMPLATE Fill(type *p, int nBegin, int nEnd, const ValType &val)
{
	for (int i=nBegin; i<=nEnd; i++)
	{
		p[i] = val;
	}
}

template <class type>
inline void UPIX_ALGORITHM_TEMPLATE Fill(type *p, int n, type *pSrc)
{
	for (int i=0; i<n; i++)
	{
		p[i] = pSrc[i];
	}
}

template <class type>
inline void UPIX_ALGORITHM_TEMPLATE Fill(type *p, int n, int val)
{
	if (val == 0)
	{
		memset(p, val, n*sizeof(type));
	}
	else
	{
		type value = (type)val;
		for (int i=0; i<n; i++)
		{
			p[i] = value;
		}
	}
}

template <class type>
inline void UPIX_ALGORITHM_TEMPLATE Fill(type *p, int nBegin, int nEnd, int val)
{
	for (int i=nBegin; i<=nEnd; i++)
	{
		p[i] = val;
	}
}




typedef struct UPIX_ALGORITHM_CLASS _UpixMemory
{
	_UpixMemory() : pBuffer(NULL), nSize(0), nPosition(0) {}

	uchar *pBuffer;
	uint nSize;
	uint nPosition;
} UpixMemory;



uint UPIX_ALGORITHM_CLASS UpixMemoryRead(void *pDst, uint nElemSize, uint nCount, UpixMemory *pMem);




#endif