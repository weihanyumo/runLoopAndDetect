/*
* Filename: UpixSortList.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_SORT_LIST_H_
#define _UPIX_SORT_LIST_H_


#include "UpixList.h"


template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixSortList : public UpixLinkedList<T>
{
public:
	UpixSortList() {}

	~UpixSortList() {}

	void InsertAscending(UpixNode<T> *pNode)
	{
		if (m_pHead == NULL || pNode->Data < m_pHead->Data)
		{
			Insert(NULL, pNode);
			return;
		}

		UpixNode<T> *pCurNode = m_pHead;
		while (pCurNode->pNext != NULL && pCurNode->pNext->Data < pNode->Data)
		{
			pCurNode = pCurNode->pNext;
		}
		Insert(pCurNode, pNode);
	}

	// Insert a value and return its pointer
	UpixNode<T> *InsertAscending(const T &val)
	{
		UpixNode<T> *pNode = new UpixNode<T>(val);
		InsertAscending(pNode);
		return pNode;
	}

	void InsertDescending(UpixNode<T> *pNode)
	{
		if (m_pHead == NULL || m_pHead->Data < pNode->Data)
		{
			Insert(NULL, pNode);
			return;
		}

		UpixNode<T> *pCurNode = m_pHead;
		while (pCurNode->pNext != NULL && pNode->Data < pCurNode->pNext->Data)
		{
			pCurNode = pCurNode->pNext;
		}
		Insert(pCurNode, pNode);
	}

	// Insert a value and return its pointer
	UpixNode<T> *InsertDescending(const T &val)
	{
		UpixNode<T> *pNode = new UpixNode<T>(val);
		InsertDescending(pNode);
		return pNode;
	}

	UpixNode<T> *SearchNotLess(const T &val)
	{
		UpixNode<T> *pNode = m_pHead;
		while (pNode != NULL && pNode->Data < val)
		{
			pNode = pNode->pNext;
		}
		return pNode;
	}

	UpixNode<T> *SearchNotGreater(const T &val)
	{
		UpixNode<T> *pNode = m_pHead;
		while (pNode != NULL && val < pNode->Data)
		{
			pNode = pNode->pNext;
		}
		return pNode;
	}
};



#endif