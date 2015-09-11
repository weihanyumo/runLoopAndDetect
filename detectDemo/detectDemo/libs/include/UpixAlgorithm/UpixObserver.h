/*
* Filename: UpixObserver.h
* Author:   Xiehong
* Date:     2012.11.10
*/

#ifndef _UPIX_OBSERVER_H_
#define _UPIX_OBSERVER_H_


#include "UpixAlgorithmTypeDef.h"


template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixLinkedList;

class UPIX_ALGORITHM_CLASS IUpixSubject;



class  UPIX_ALGORITHM_CLASS IUpixObserver
{
public:
	IUpixObserver() {}
	virtual ~IUpixObserver() {}

	virtual int OnUpdate(IUpixSubject *pSubject, void *pData) = 0;
};



class UPIX_ALGORITHM_CLASS IUpixSubject
{
public:
	IUpixSubject();
	virtual ~IUpixSubject();

	// Register a observer
	void Attach(IUpixObserver *pObserver);

	// Remove a observer
	void Detach(IUpixObserver *pObserver);

	// Send a message. If observer is null, do broadcast
	int Notify(IUpixObserver *pObserver, void *pData);

protected:
	// Broadcast message
	void Notify(void *pData);

protected:
	UpixLinkedList<IUpixObserver *> *m_plistObserver;
};



#endif