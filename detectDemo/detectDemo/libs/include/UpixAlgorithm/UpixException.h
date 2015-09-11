/*
* Filename: UpixException.h
* Author:   Xiehong
* Date:     2012.11.12
*/

#ifndef _UPIX_EXCEPTION_H_
#define _UPIX_EXCEPTION_H_


#include "UpixAlgorithmTypeDef.h"
#include <exception>
using namespace std;



class UPIX_ALGORITHM_CLASS UpixException : public exception
{
public:
	UpixException();
	UpixException(const char *strMsg);
	UpixException(int _code, const char *_err, const char *_func, const char *_file, int _line);
	virtual ~UpixException() throw() {}

public:
	virtual const char *What() const throw();
	void FormatMessage();

public:
	char msg[512];
	char err[256];
	char func[128];
	char file[128];
	int line;
	int code;
};



#endif