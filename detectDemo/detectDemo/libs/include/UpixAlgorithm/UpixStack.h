/*
* Filename: UpixStack.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_STACK_H_
#define _UPIX_STACK_H_


#include "UpixAlgorithmTypeDef.h"



template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixStack
{
public:
	UpixStack()
	{
		m_pData = NULL;
		m_Size = 0;
		m_Capacity = 0;
	}

	UpixStack(int n)
	{
		m_pData = new T[n];
		m_Size = 0;
		m_Capacity = n;
	}

	~UpixStack()
	{
		Clear();
	}

public:
	void Resize(int n)
	{
		if (m_Capacity != n)
		{
			Clear();
			m_pData = new T[n];
			m_Capacity = n;
		}
		m_Size = 0;
	}

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

	bool Pop(T &val)
	{
		if (IsEmpty())
		{
			return false;
		}

		m_Size--;
		val = m_pData[m_Size];
		m_pData[m_Size] = T();
		return true;
	}

	bool Push(const T &val)
	{
		if (IsFull())
		{
			return false;
		}

		m_pData[m_Size] = val;
		m_Size++;
		return true;
	}

	int Size() { return m_Size; }

	int Capacity() { return m_Capacity; }

	bool IsEmpty() { return m_Size == 0; }

	bool IsFull() { return m_Size == m_Capacity; }

	T &operator[](int i)
	{
		assert(i>=0 && i<m_Size);
		return m_pData[nId];
	}

private:
	int m_Size;
	int m_Capacity;
	T *m_pData;
};


#endif