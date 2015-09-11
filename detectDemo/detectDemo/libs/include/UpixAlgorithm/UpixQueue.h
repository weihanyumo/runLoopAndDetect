/*
* Filename: UpixQueue.h
* Author:   Xiehong
* Date:     2013.7.6
*/
#ifndef _UPIX_QUEUE_H_
#define _UPIX_QUEUE_H_


#include "UpixAlgorithmTypeDef.h"
#include "UpixList.h"


// Linear queue
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixQueue
{
public:
	UpixQueue()
	{
		
	}

	~UpixQueue()
	{
		Clear();
	}

	//入队
	bool EnQueue(const T &val)
	{
		return m_list.AddTail(val);
	}
	//出队
	bool DeQueue(T &val)
	{
		UpixNode<T> *pNode = m_list.PopHead();
		if (pNode == NULL)
		{
			return false;
		}

		val = pNode->Data;
		delete pNode;
		pNode = NULL;
		return true;
	}

	void Clear()
	{
		m_list.Clear();
	}

	bool IsEmpty() { return m_list.IsEmpty(); }

	int Size()
	{
		return m_list.Length();
	}

private:
	UpixLinkedList<T> m_list;
};




// Circular queue
// The element in circular queue must overwrite its operator=()
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixCircularQueue
{
public:
	UpixCircularQueue(int nCapacity) : m_nCapacity(nCapacity)
	{
		m_pData = new T[m_nCapacity];
		m_nSize = 0;
		m_iHeadId = 0;
		m_iTailId = 0;
	}

	~UpixCircularQueue()
	{
		Clear();
	}

public:
	void Clear()
	{
		if (m_pData != NULL)
		{
			delete m_pData;//[]
			m_pData = NULL;
		}

		m_nCapacity = 0;
		m_nSize = 0;
		m_iHeadId = 0;
		m_iTailId = 0;
	}

	bool EnQueue(const T &val)
	{
		if (IsFull())
		{
			return false;
		}

		m_pData[m_iTailId++] = val;
		if (m_iTailId == m_nCapacity)
		{
			m_iTailId = 0;
		}

		m_nSize++;
		return true;
	}

	bool DeQueue(T &val)
	{
		if (IsEmpty())
		{
			return false;
		}

		val = m_pData[m_iHeadId];

		m_pData[m_iHeadId++] = T();
		if (m_iHeadId == m_nCapacity)
		{
			m_iHeadId = 0;
		}

		m_nSize--;
		return true;
	}

	int Size() { return m_nSize; }

	int Capacity() { return m_nCapacity; }

	bool IsEmpty() { return m_nSize == 0; }

	bool IsFull() { return m_nSize == m_nCapacity; }

	T &operator[](int i)
	{
		assert(i>=0 && i<m_nSize);
		int nId = m_iHeadId + i;
		if (nId >= m_nCapacity)
		{
			nId -= m_nCapacity;
		}
		return m_pData[nId];
	}

private:
	T *m_pData;
	int m_nCapacity;
	int m_nSize;
	int m_iHeadId;
	int m_iTailId;
};



#endif