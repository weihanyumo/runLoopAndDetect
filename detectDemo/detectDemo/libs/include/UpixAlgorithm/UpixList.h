/*
* Filename: UpixList.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_LINKED_LIST_H_
#define _UPIX_LINKED_LIST_H_

#include "UpixAlgorithmTypeDef.h"
#include <iostream>
using namespace std;


// class of node
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixNode
{
public:
	// default constructor
	UpixNode(const T &val = T(), UpixNode *_pPrev = NULL, UpixNode *_pNext = NULL) 
	{
		Data = val;
		pPrev = _pPrev;
		pNext = _pNext;
	}

	// copy constructor
	UpixNode(const UpixNode &other)
	{
		Data = other.Data;
		pPrev = other.pPrev;
		pNext = other.pNext;
	}

	// operator =()
	UpixNode<T>& operator=(const UpixNode &rhs)
	{
		if (this != &rhs)
		{
			Data = rhs.Data;
			pPrev = rhs.pPrev;
			pNext = rhs.pNext;
		}
		return *this;
	}

	// destructor
	~UpixNode() {}

public:
	T Data;
	UpixNode *pPrev;
	UpixNode *pNext;
};


//---------------------------------------------------------------------
// class of linked list
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixLinkedList
{
public:
	UpixLinkedList(): m_pHead(NULL), m_pTail(NULL), m_nSize(0) {}

	UpixLinkedList(UpixNode<T> *pHead, UpixNode<T> *pTail, int nSize)
		: m_pHead(pHead), m_pTail(pTail), m_nSize(nSize) {}
	
	// destructor
	~UpixLinkedList() { Clear(); }

public:
	void Clear();

	void Drop() { m_pHead = NULL; m_pTail = NULL; m_nSize = 0; }

	// judge if is empty
	bool IsEmpty() { return m_pHead == NULL; }

	// length of linked list
	int Length() const { return m_nSize; }

	// get node at position
	UpixNode<T> *GetNode(int pos);

	// set data at position
	bool SetData(int pos, const T &val);

	// insert data after pNode
	bool Insert(UpixNode<T> *pNode, const T &val);

	// insert data after pos
	bool InsertAt(int pos, const T &val);

	// insert node after pNode
	bool Insert(UpixNode<T> *pNode, UpixNode<T> *pNewNode);

	// add data before the head
	bool AddHead(const T &val);

	bool AddHead(UpixNode<T> *pNode);

	// append data to the tail
	bool AddTail(const T &val);

	bool AddTail(UpixNode<T> *pNode);

	// add another list at tail
	bool AddList(UpixLinkedList<T> *pList);

	// pop the first node
	UpixNode<T> *PopHead();

	// pop the last node
	UpixNode<T> *PopTail();

	// remove a node and return the next
	UpixNode<T> *Remove(UpixNode<T> *pNode);

	// remove the node at position
	UpixNode<T> *Remove(int pos);

	// remove the first key word in list
	bool RemoveVal(const T &val);

	// flip the linked list
	void Flip();

	// Linear search the first key word in list
	UpixNode<T> *Search(const T &val);

	// get head of linked list
	UpixNode<T> *Head() const { return m_pHead; }

	UpixNode<T> *Tail() { return m_pTail; }

	// out put list from *pNode
	void Print(UpixNode<T> *pNode);

	// out put whole list
	void Print() { Print(m_pHead); }

protected:
	UpixNode<T> *m_pHead;
	UpixNode<T> *m_pTail;
	int m_nSize;
};


// implement of CLinkList

template<class T>
void UpixLinkedList<T>::Clear()
{
	if (m_pHead != NULL)
	{
		UpixNode<T> *pNode = m_pHead;
		UpixNode<T> *pTail = NULL;
		while (pNode != NULL)
		{
			pTail = pNode->pNext;
			delete pNode;
			pNode = pTail;
		}
		m_pHead = NULL;
		m_pTail = NULL;
		m_nSize = 0;
	}
}

template <class T>
UpixNode<T>* UpixLinkedList<T>::GetNode(int pos)
{
	if (pos >= m_nSize || pos < 0)
	{
		return NULL;
	}

	UpixNode<T> *pNode = NULL;
	if (pos <= m_nSize/2)
	{
		pNode = m_pHead;
		for (int i=0; i<pos; i++)
		{
			pNode = pNode->pNext;
		}
	}
	else
	{
		pNode = m_pTail;
		for (int i=m_nSize-1; i>pos; i--)
		{
			pNode = pNode->pPrev;
		}
	}

	return pNode;
}

template <class T>
bool UpixLinkedList<T>::SetData(int pos, const T &val)
{
	UpixNode<T> *pNode = GetNode(pos);
	if (pNode != NULL)
	{
		pNode->Data = val;
		return true;
	}
	return false;
}

template <class T>
bool UpixLinkedList<T>::Insert(UpixNode<T> *pNode, const T &val)
{
	UpixNode<T> *pNewNode = new UpixNode<T>(val);
	return Insert(pNode, pNewNode);
}

template <class T>
bool UpixLinkedList<T>::InsertAt(int pos, const T &val)
{
	if (pos == -1)
	{
		return Insert(NULL, val);
	}

	UpixNode<T> *pNode = GetNode(pos);
	if (pNode != NULL)
	{
		return Insert(pNode, val);;
	}
	return false;
}

template <class T>
bool UpixLinkedList<T>::Insert(UpixNode<T> *pNode, UpixNode<T> *pNewNode)
{
	if (pNode == NULL)
	{
		return AddHead(pNewNode);
	}

	if (pNode == m_pTail)
	{
		return AddTail(pNewNode);
	}

	pNewNode->pNext = pNode->pNext;
	pNewNode->pPrev = pNode;
	pNode->pNext->pPrev = pNewNode;
	pNode->pNext = pNewNode;
	m_nSize++;
	return true;
}

template <class T>
bool UpixLinkedList<T>::AddHead(UpixNode<T> *pNode)
{
	if (m_pHead == NULL)
	{
		m_pTail = pNode;
	}
	else
	{
		pNode->pNext = m_pHead;
		m_pHead->pPrev = pNode;
	}

	m_pHead = pNode;
	m_pHead->pPrev = NULL;
	m_nSize++;
	return true;
}

template <class T>
bool UpixLinkedList<T>::AddHead(const T &val)
{
	UpixNode<T> *pNode = new UpixNode<T>(val);
	return AddHead(pNode);
}

template <class T>
bool UpixLinkedList<T>::AddTail(UpixNode<T> *pNode)
{
	if (m_pTail == NULL)
	{
		m_pHead = pNode;
	}
	else
	{
		m_pTail->pNext = pNode;
		pNode->pPrev = m_pTail;
	}
	
	m_pTail = pNode;
	m_pTail->pNext = NULL;
	m_nSize++;
	return true;
}

template <class T>
bool UpixLinkedList<T>::AddTail(const T &val)
{
	UpixNode<T> *pNode = new UpixNode<T>(val);
	return AddTail(pNode);
}

template <class T>
bool UpixLinkedList<T>::AddList(UpixLinkedList<T> *pList)
{
	UpixNode<T> *pNode = pList->Head();
	if (m_pHead == NULL)
	{
		m_pHead = pNode;
	}
	else
	{
		m_pTail->pNext = pNode;
		pNode->pPrev = m_pTail;
	}
	m_pTail = pList->Tail();
	m_nSize += pList->Length();
	pList->Drop();

	return true;
}

template <class T>
UpixNode<T> *UpixLinkedList<T>::PopHead()
{
	UpixNode<T> *pNode = NULL;
	if (m_pHead != NULL)
	{
		pNode = m_pHead;
		m_pHead = m_pHead->pNext;
		pNode->pNext = NULL;
		if (m_pHead != NULL)
		{
			m_pHead->pPrev = NULL;
		}
		else
		{
			m_pTail = NULL;
		}
		m_nSize--;
	}
	return pNode;
}

template <class T>
UpixNode<T> *UpixLinkedList<T>::PopTail()
{
	UpixNode<T> *pNode = NULL;
	if (m_pTail != NULL)
	{
		pNode = m_pTail;
		m_pTail = m_pTail->pPrev;
		pNode->pPrev = NULL;
		if (m_pTail != NULL)
		{
			m_pTail->pNext = NULL;
		}
		else
		{
			m_pHead = NULL;
		}
		m_nSize--;
	}
	return pNode;
}

template <class T>
UpixNode<T> *UpixLinkedList<T>::Remove(UpixNode<T> *pNode)
{
	if (pNode == NULL || m_pHead == NULL)
	{
		return NULL;
	}

	UpixNode<T> *pNextNode = NULL;

	if (pNode == m_pHead)
	{
		PopHead();
		pNextNode = m_pHead;
	}
	else if (pNode == m_pTail)
	{
		PopTail();
	}
	else
	{
		pNextNode = pNode->pNext;
		pNode->pNext->pPrev = pNode->pPrev;
		pNode->pPrev->pNext = pNode->pNext;
		m_nSize--;
	}
	delete pNode;
	pNode = NULL;

	return pNextNode;
}

template <class T>
UpixNode<T> *UpixLinkedList<T>::Remove(int pos)
{
	UpixNode<T> *pNode = GetNode(pos);
	return Remove(pNode);
}

template <class T>
void UpixLinkedList<T>::Flip()
{
	UpixNode<T> *pNode = m_pHead;
	UpixNode<T> *pNextNode = NULL;
	UpixNode<T> *pTemp = NULL;
	while (pNode != NULL)
	{
		pNextNode = pNode->pNext;
		pTemp = pNode->pNext;
		pNode->pNext = pNode->pPrev;
		pNode->pPrev = pTemp;
		pNode = pNextNode;
	}
	pTemp = m_pHead;
	m_pHead = m_pTail;
	m_pTail = pTemp;
}

template <class T>
UpixNode<T>* UpixLinkedList<T>::Search(const T &val)
{
	UpixNode<T> *pNode = m_pHead;
	while (pNode != NULL)
	{
		if (pNode->Data == val)
		{
			break;
		}
		pNode = pNode->pNext;
	}
	return pNode;
}

template <class T>
bool UpixLinkedList<T>::RemoveVal(const T &val)
{
	bool bSearched = false;
	UpixNode<T> *pNode = Search(val);
	if (pNode != NULL)
	{
		Remove(pNode);
		bSearched = true;
	}
	return bSearched;
}

template <class T>
void UpixLinkedList<T>::Print(UpixNode<T> *pNode)
{
	while (pNode != NULL)
	{
		cout<<pNode->Data<<" ";
		pNode = pNode->pNext;
	}
	cout<<endl;
}


//----------------------------------------------------------------------
// class of circular linked list
template<class T>
class UPIX_ALGORITHM_TEMPLATE UpixCircularLinkedList
{
public:
	// default constructor
	UpixCircularLinkedList() { m_pHead = NULL; }

	// destructor
	~UpixCircularLinkedList() { Clear(); }

public:
	void Clear();

	// judge if empty
	bool IsEmpty() { return m_pHead == NULL; }

	// add data to tail
	bool Add(const T &val);

	// remove the node
	bool Remove(UpixNode<T> *pNode);

	// get head of the list
	UpixNode<T>* GetHead() { return m_pHead; }

	// out put list
	void Print();

private:
	UpixNode<T> *m_pHead;
};

template<class T>
void UpixCircularLinkedList<T>::Clear()
{
	if (m_pHead != NULL)
	{
		UpixNode<T> *pNode = m_pHead;
		UpixNode<T> *pNext = NULL;
		do
		{
			pNext = pNode->pNext;
			delete pNode;
			pNode = pNext;
		}while (pNode != m_pHead);
		m_pHead = NULL;
	}
}

template<class T>
bool UpixCircularLinkedList<T>::Add(const T &val)
{
	UpixNode<T> *pNode = new UpixNode<T>(val);

	if (m_pHead == NULL)
	{
		m_pHead = pNode;
		m_pHead->pNext = m_pHead;
		m_pHead->pPrev = m_pHead;
		return true;
	}

	UpixNode<T> *pTail = m_pHead->pPrev;
	pTail->pNext = pNode;
	pNode->pPrev = pTail;
	pNode->pNext = m_pHead;
	m_pHead->pPrev = pNode;
	return true;
}

template<class T>
bool UpixCircularLinkedList<T>::Remove(UpixNode<T> *pNode)
{
	if (pNode == NULL || m_pHead == NULL)
	{
		return false;
	}

	if (m_pHead == m_pHead->pNext && pNode == m_pHead)
	{
		m_pHead = NULL;
		delete pNode;
		return true;
	}

	if (pNode == m_pHead)
	{
		m_pHead = m_pHead->pNext;
	}
	
	pNode->pPrev->pNext = pNode->pNext;
	pNode->pNext->pPrev = pNode->pPrev;
	delete pNode;
	return true;
}

template<class T>
void UpixCircularLinkedList<T>::Print()
{
	if (m_pHead != NULL)
	{
		UpixNode<T> *pNode = m_pHead;
		do
		{
			cout<<pNode->Data<<" ";
			pNode = pNode->GetNext();
		}while (pNode != m_pHead);
		cout<<endl;
	}
}


// I/O stream
// for UpixNode
template <class T>
ostream& operator<<(ostream &os, UpixNode<T> &node)
{
	os << node.Data;
	return os;
}

// for CLinkList
template <class T>
ostream& operator<<(ostream &os, UpixLinkedList<T> &linklist)
{
	UpixNode<T> *pNode = linklist.GetHead();
	while (pNode != NULL)
	{
		os << pNode->Data<< " ";
		pNode = pNode->pNext;
	}
	return os;
}

// for CCircularList
template <class T>
ostream& operator<<(ostream &os, UpixCircularLinkedList<T> &linklist)
{
	UpixNode<T> *pHead = linklist.GetHead();
	if (pHead != NULL)
	{
		UpixNode<T> *pNode = pHead;
		do
		{
			os<<pNode->Data<<" ";
			pNode = pNode->pNext;
		}while (pNode != pHead);
	}
	return os;
}


#endif