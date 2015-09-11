/*
* Filename: UpixArray.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_ARRAY_H_
#define _UPIX_ARRAY_H_


#include "UpixAlgorithmTypeDef.h"
#include "UpixMemory.h"
#include "UpixSort.h"


template <class type>
class UPIX_ALGORITHM_TEMPLATE UpixArray
{
public:
	UpixArray() : m_Capacity(0), m_Size(0), m_pData(NULL) {}

	UpixArray(int nSize)
	{
		m_Capacity = nSize;
		m_Size = nSize;
		m_pData = new type[m_Capacity];
		type val = type();
		Fill(m_pData, m_Size, val);
	}

	UpixArray(int nSize, const type &val)
	{
		m_Capacity = nSize;
		m_Size = nSize;
		m_pData = new type[m_Capacity];
		Fill(m_pData, m_Size, val);
	}

	UpixArray(int nSize, type data[])
	{
		m_Capacity = nSize;
		m_Size = nSize;
		m_pData = new type[m_Capacity];
		Fill(m_pData, m_Size, data);
	}

	UpixArray(const UpixArray &other)
	{
		m_Capacity = other.m_Capacity;
		m_Size = other.m_Size;
		m_pData = new type[m_Capacity];
		Fill(m_pData, m_Size, other.m_pData);
	}

	~UpixArray() { Clear(); }

public:
	void Clear()
	{
		if (m_pData != NULL)
		{
			delete[] m_pData;
			m_pData = NULL;
		}
		m_Size = 0;
		m_Capacity = 0;
	}

	void Resize(int nSize)
	{
		Resize(nSize, type());
	}

	void Resize(int nSize, const type &val)
	{
		Clear();
		m_Capacity = nSize;
		m_Size = nSize;
		m_pData = new type[m_Capacity];
		Fill(m_pData, m_Size, val);
	}

	void Reserve(int nCapacity)
	{
		Clear();
		m_Capacity = nCapacity;
		m_pData = new type[m_Capacity];
	}

	void PushBack(const type &val)
	{
		if (m_Size + 1 > m_Capacity)	// If not enough
		{
			type *p = m_pData;

			m_Capacity = 2*m_Capacity + 1;
			m_pData = new type[m_Capacity];
			Fill(m_pData, m_Size, p);

			type e = type();
			Fill(m_pData, m_Size, m_Capacity - 1, e);

			delete[] p;
		}

		m_pData[m_Size] = val;
		m_Size++;
	}

	void Sort()
	{
		if (m_Size > 1)
		{
			HeapSort(m_pData, m_Size);
		}
	}

	// Must sort the array first
	int BinarySearch(const type &val) const
	{
		return BinarySearch(val, 0, m_Size - 1);
	}

	// Linear search for unsorted array
	int Search(const type &val)
	{
		for (int i=0; i<m_Size; i++)
		{
			if (val == m_pData[i])
			{
				return i;
			}
		}
		return -1;
	}

	inline int Capacity() const { return m_Capacity; }

	inline int Size() const { return m_Size; }

	void SetData(type data[])
	{
		Fill(m_pData, m_Size, data);
	}

	inline type *GetData() const { return m_pData; }

	inline bool IsEmpty() { return m_Size == 0; }

	UpixArray &operator=(const UpixArray &rhs)
	{
		if (this != &rhs)
		{
			if (m_Size != rhs.m_Size)
			{
				if (m_Capacity < rhs.m_Size)
				{
					Clear();
					m_Capacity = rhs.m_Capacity;
					m_pData = new type[m_Capacity];
				}
				m_Size = rhs.m_Size;
			}
			Fill(m_pData, m_Size, rhs.m_pData);
		}
		return *this;
	}

	inline type &operator[](int n) { return m_pData[n]; }
	inline const type &operator[](int n) const { return m_pData[n]; }

protected:
	int BinarySearch(const type &element, int left, int right) const
	{
		if(m_Size == 0)
		{
			return -1;
		}

		int m = 0;
		do 
		{
			m = (left + right) >> 1;
			if(element < m_pData[m])
				right = m - 1;
			else
				left = m + 1;

		} while ((element < m_pData[m] || m_pData[m] < element) && left <= right);

		if(!(element < m_pData[m]) && !(m_pData[m] < element))
			return m;
		return -1;
	}

protected:
	int m_Capacity;
	int m_Size;
	type *m_pData;
};


typedef UpixArray<uchar> UpixArrayuc;
typedef UpixArray<char> UpixArrayc;
typedef UpixArray<ushort> UpixArrayus;
typedef UpixArray<short> UpixArrays;
typedef UpixArray<uint> UpixArrayui;
typedef UpixArray<int> UpixArrayi;
typedef UpixArray<float> UpixArrayf;
typedef UpixArray<double> UpixArrayd;



#endif