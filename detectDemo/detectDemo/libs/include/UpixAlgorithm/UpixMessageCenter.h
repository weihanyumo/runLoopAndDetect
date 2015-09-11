/*
* Filename: UpixMessageCenter.h
* Author:   Xiehong
* Date:     2012.11.10
*/

#ifndef _UPIX_MESSAGE_CENTER_H_
#define _UPIX_MESSAGE_CENTER_H_


#include "UpixAlgorithmTypeDef.h"
#include "UpixObserver.h"



// Definition of messages
typedef enum _eMessageType
{
	// Messages for system reserved
	MSG_DEFAULT = -1,

	MSG_THREAD_END = 1024,
	MSG_THREAD_EXCEPTION,
	MSG_THREAD_CANCELED,
	MSG_THREAD_KILLED,
	MSG_THREAD_SUSPEND,
	MSG_THREAD_RESUME,
	MSG_THREAD_EXIT,

	// Messages for user
	MSG_USER_DEFAULT = 10000,
	
}eMessageType;



class UPIX_ALGORITHM_CLASS IUpixMessageReceiver : public IUpixObserver
{
public:
	IUpixMessageReceiver() {}
	virtual ~IUpixMessageReceiver() {}

	virtual int OnUpdate(IUpixSubject *pSubject, void *pData) { return -1; }

	virtual int OnMessage(int eMsg, long WParam, void *LParam)
	{
		return 0;
	}
};


class UPIX_ALGORITHM_CLASS IUpixMessageCenter : public IUpixSubject
{
public:
	IUpixMessageCenter() {}
	virtual ~IUpixMessageCenter() {}

	// Send a message. If receiver is null, do broadcast
	int SendMessages(IUpixMessageReceiver *pMsgReceiver,	// Receiver
					int eMsg,								// Message type
					long WParam = 0,						// Parameter 1
					void *LParam = NULL);					// Parameter 2

	// Post a message. If receiver is null, do broadcast
	int PostMessages(IUpixMessageReceiver *pMsgReceiver,	// Receiver
					int eMsg,								// Message type
					long WParam = 0,						// Parameter 1
					void *LParam = NULL);					// Parameter 2

protected:
	// Broadcast message
	void Notify(int eMsg, long WParam, void *LParam);
};



#endif
