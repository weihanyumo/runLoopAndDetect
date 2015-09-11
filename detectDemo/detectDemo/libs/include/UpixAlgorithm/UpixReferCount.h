/*
* Filename: UpixReferCount.h
* Author:   Xiehong
* Date:     2013.10.10
*/
#ifndef _UPIX_REFER_COUNT_H_
#define _UPIX_REFER_COUNT_H_


#include "UpixAlgorithmTypeDef.h"



// Abstract class of reference count object
class UPIX_ALGORITHM_CLASS IUpixReferCount
{
public:
	// Increase reference count
	void AddReference();

	// Decrease reference count. Delete itself when reference count reduce to 0.
	void ReduceReference();

protected:
	IUpixReferCount();
	IUpixReferCount(const IUpixReferCount &rhs);
	virtual ~IUpixReferCount() = 0;

protected:
	IUpixReferCount &operator=(const IUpixReferCount &rhs);
	
private:
	int m_nReferenceCount;
};



// Base class of smart pointer
template <class T>
class IUpixReferCountPtr
{
protected:
	IUpixReferCountPtr() : m_pHolder(NULL) {}

	IUpixReferCountPtr(const IUpixReferCountPtr &rhs)
		: m_pHolder(rhs.m_pHolder)
	{
		if (m_pHolder != NULL)
		{
			m_pHolder->AddReference();
		}
	}

public:
	virtual ~IUpixReferCountPtr()
	{
		if (m_pHolder != NULL)
		{
			m_pHolder->ReduceReference();
		}
	}

public:
	IUpixReferCountPtr &operator=(const IUpixReferCountPtr &rhs)
	{
		if (m_pHolder != rhs.m_pHolder)
		{
			if (m_pHolder != NULL)
			{
				m_pHolder->ReduceReference();
			}

			m_pHolder = rhs.m_pHolder;
			if (m_pHolder != NULL)
			{
				m_pHolder->AddReference();
			}
		}

		return *this;
	}

	bool operator==(const IUpixReferCountPtr &rhs) const
	{
		return m_pHolder == rhs.m_pHolder;
	}

	bool operator!=(const IUpixReferCountPtr &rhs) const
	{
		return m_pHolder != rhs.m_pHolder;
	}

	// Override operators for pointer
	inline T *operator->() const { return GetPtr(); }
	inline T &operator*() const { return GetObject(); }

	// Return raw pointer
	T *GetPtr() const
	{
		if (m_pHolder != NULL)
		{
			return m_pHolder->pRawPtr;
		}
		return NULL;
	}

	T &GetObject() const
	{
		assert(m_pHolder != NULL && m_pHolder->pRawPtr != NULL);
		return *(m_pHolder->pRawPtr);
	}

	void Clear()
	{
		if (m_pHolder != NULL)
		{
			m_pHolder->ReduceReference();
			m_pHolder = NULL;
		}
	}

protected:
	class IUpixReferCountHolder : public IUpixReferCount
	{
	public:
		IUpixReferCountHolder() : pRawPtr(NULL) {}
		virtual ~IUpixReferCountHolder() {}

	public:
		T *pRawPtr;
	};

protected:
	IUpixReferCountHolder *m_pHolder;
};




// Concrete class of smart object pointer
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixReferCountPtr : public IUpixReferCountPtr<T>
{
public:
	UpixReferCountPtr(T *ptr = NULL)
	{
		if (ptr != NULL)
		{
			this->m_pHolder = new UpixReferCountHolder;
			this->m_pHolder->pRawPtr = ptr;
			this->m_pHolder->AddReference();
		}
	}

	UpixReferCountPtr(const UpixReferCountPtr &rhs) : IUpixReferCountPtr<T>(rhs) {}

	virtual ~UpixReferCountPtr() {}

private:
	class UpixReferCountHolder : public IUpixReferCountPtr<T>::IUpixReferCountHolder
	{
	public:
		UpixReferCountHolder() {}

		virtual ~UpixReferCountHolder()
		{
			if (this->pRawPtr != NULL)
			{
				delete this->pRawPtr;
				this->pRawPtr = NULL;
			}
		}
	};
};




// Concrete class of smart array pointer
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixReferCountPtr_Array : public IUpixReferCountPtr<T>
{
public:
	UpixReferCountPtr_Array(T *ptr = NULL)
	{
		if (ptr != NULL)
		{
			this->m_pHolder = new UpixReferCountHolder;
			this->m_pHolder->pRawPtr = ptr;
			this->m_pHolder->AddReference();
		}
	}

	UpixReferCountPtr_Array(const UpixReferCountPtr_Array &rhs) : IUpixReferCountPtr<T>(rhs) {}

	virtual ~UpixReferCountPtr_Array() {}

private:
	class UpixReferCountHolder : public IUpixReferCountPtr<T>::IUpixReferCountHolder
	{
	public:
		UpixReferCountHolder() {}

		virtual ~UpixReferCountHolder()
		{
			if (this->pRawPtr != NULL)
			{
				delete[] this->pRawPtr;
				this->pRawPtr = NULL;
			}
		}
	};
};




#endif