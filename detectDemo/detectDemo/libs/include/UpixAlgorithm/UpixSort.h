/*
* Filename: UpixSort.h
* Author:   Xiehong
* Date:     2012.11.20
*/

#ifndef _UPIX_SORT_H_
#define _UPIX_SORT_H_


#include "UpixAlgorithmTypeDef.h"
#include <math.h>


// Quick sort
template<class T>
void UPIX_ALGORITHM_TEMPLATE QuickSort(T *pData, int low, int high, bool bAscending = true)
{
	int i = low;
	int j = high;

	//int RandIndex = (low + high)/2;
	int RandIndex = rand()%(high - low +1) + low;
	T middle = pData[RandIndex];

	T temp;
	if (bAscending)	// Ascending
	{
		while(i <= j)
		{
			while((pData[i] < middle) && (i < high))
				i++;
			while((pData[j] > middle) && (j > low))
				j--;

			if (i <= j)
			{
				temp = pData[i];
				pData[i] = pData[j];
				pData[j] = temp;
				i++;
				j--;
			}
		}
	}
	else	// Descending
	{
		while(i <= j)
		{
			while((pData[i] > middle) && (i < high))
				i++;
			while((pData[j] < middle) && (j > low))
				j--;

			if (i <= j)
			{
				temp = pData[i];
				pData[i] = pData[j];
				pData[j] = temp;
				i++;
				j--;
			}
		}
	}

	if (j > low)
	{
		QuickSort(pData, low, j, bAscending);
	}
	if (i < high)
	{
		QuickSort(pData, i, high, bAscending);
	}
}

template<class T>
void UPIX_ALGORITHM_TEMPLATE QuickSortIndex(const T *pData, int low, int high, int *indices, bool bAscending = true)
{
	int i = low;
	int j = high;

	//int RandIndex = (low + high)/2;
	int RandIndex = rand()%(high - low +1) + low;
	T middle = pData[indices[RandIndex]];

	int temp = 0;
	if (bAscending)	// Ascending
	{
		while(i <= j)
		{
			while((pData[indices[i]] < middle) && (i < high))
				i++;
			while((pData[indices[j]] > middle) && (j > low))
				j--;

			if (i <= j)
			{
				temp = indices[i];
				indices[i] = indices[j];
				indices[j] = temp;
				i++;
				j--;
			}
		}
	}
	else	// Descending
	{
		while(i <= j)
		{
			while((pData[indices[i]] > middle) && (i < high))
				i++;
			while((pData[indices[j]] < middle) && (j > low))
				j--;

			if (i <= j)
			{
				temp = indices[i];
				indices[i] = indices[j];
				indices[j] = temp;
				i++;
				j--;
			}
		}
	}

	if (j > low)
	{
		QuickSortIndex(pData, low, j, indices, bAscending);
	}
	if (i < high)
	{
		QuickSortIndex(pData, i, high, indices, bAscending);
	}
}



// Heap sort
template <class T>
void HeapSink(T* array, int element, int max)
{
	while((element << 1) < max)
	{
		int j = (element << 1);
		if(j + 1 < max && array[j] < array[j + 1])
			j = j + 1;
		if(array[element] < array[j])
		{
			T t = array[j];
			array[j] = array[element];
			array[element] = t;
			element = j;
		}
		else
		{
			return;
		}
	}
}

template <class T>
void UPIX_ALGORITHM_TEMPLATE HeapSort(T* array, int size)
{
	T* virtualArray = array - 1;
	int virtualSize = size + 2;
	int i;
	for(i = ((size - 1) / 2); i >= 0; i--)
		HeapSink(virtualArray, i + 1, virtualSize - 1);

	for(i = size - 1; i >= 0; i--)
	{
		T t = array[0];
		array[0] = array[i];
		array[i] = t;
		HeapSink(virtualArray, 1, i + 1);
	}
}


#endif