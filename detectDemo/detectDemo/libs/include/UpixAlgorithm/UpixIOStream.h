/*
* Filename: UpixIOStream.h
* Author:   Xiehong
* Date:     2013.10.10
*/
#ifndef _UPIX_IO_STREAM_H_
#define _UPIX_IO_STREAM_H_


#include "UpixAlgorithmTypeDef.h"
#include "UpixException.h"



// Exceptions for I/O stream
class UPIX_ALGORITHM_CLASS UpixIOStreamException : public UpixException
{
public:
	UpixIOStreamException(const char *strMsg) : UpixException(strMsg) {}
	UpixIOStreamException(int _code, const char *_err, const char *_func, const char *_file, int _line)
		: UpixException(_code, _err, _func, _file, _line)
	{

	}
};



class UPIX_ALGORITHM_CLASS IUpixIOStream
{
public:
	IUpixIOStream() : m_bOpen(false) {}
	virtual ~IUpixIOStream() {}

public:
	bool IsOpen() { return m_bOpen; }
	virtual uint Read(void *pDst, uint nElemSize, uint nCount) = 0;
	virtual void Close() = 0;

protected:
	bool m_bOpen;
};




typedef struct _UpixMemory UpixMemory;

class UPIX_ALGORITHM_CLASS UpixMemoryStream : public IUpixIOStream
{
public:
	UpixMemoryStream(uchar *pBuffer, uint nSize);
	virtual ~UpixMemoryStream();

public:
	virtual uint Read(void *pDst, uint nElemSize, uint nCount);
	virtual void Close();

private:
	UpixMemory *m_pMem;
};




class UPIX_ALGORITHM_CLASS UpixFileStream : public IUpixIOStream
{
public:
	UpixFileStream(const char *strFile);
	virtual ~UpixFileStream();

public:
	virtual uint Read(void *pDst, uint nElemSize, uint nCount);
	virtual void Close();

private:
	FILE *m_pFile;
};




class UPIX_ALGORITHM_CLASS UpixFileStreamStream : public IUpixIOStream
{
public:
	UpixFileStreamStream(FILE *pFile);
	virtual ~UpixFileStreamStream();

public:
	virtual uint Read(void *pDst, uint nElemSize, uint nCount);
	virtual void Close();

private:
	FILE *m_pFile;
};




#endif