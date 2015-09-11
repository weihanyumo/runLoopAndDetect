/*
* Filename: UpixLine.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_LINE_H_
#define _UPIX_LINE_H_


#include "UpixAlgorithmTypeDef.h"
#include "UpixPoint.h"
#include <math.h>
#include <iostream>
using namespace std;



//-------------------------UpixLine-----------------------
template <class Type>
class UPIX_ALGORITHM_TEMPLATE UpixLine
{
public:
	Type a;
	Type b;
	Type c;

public: 
	// Not a line
	UpixLine() : a(0), b(0), c(0) {}

	// la*x + lb*y + lc = 0
	UpixLine(Type la, Type lb, Type lc): a(la), b(lb), c(lc) {}
	UpixLine(Type const vl[3]): a(vl[0]), b(vl[1]), c(vl[2]) {}

	// y = lk*x + lb
	UpixLine(Type lk, Type lb): a(lk), b(-1), c(lb) {}

	// (x - x0)(y1 - y0) = (y - y0)(x1 - x0)
	UpixLine(Type x0, Type y0, Type x1, Type y1): a(y0 - y1), b(x1 - x0), c(x0*y1 - x1*y0) {}
	UpixLine(const UpixPoint2D<Type> &p0, const UpixPoint2D<Type> &p1): a(p0.y - p1.y), b(p1.x - p0.x), c(p0.x*p1.y - p1.x*p0.y) {}

	template <class ValType>
	UpixLine(const UpixPoint2D<ValType> &p0, const UpixPoint2D<ValType> &p1): a(p0.y - p1.y), b(p1.x - p0.x), c(p0.x*p1.y - p1.x*p0.y) {}

	// y - y0 = lk*(x - x0)
	UpixLine(const UpixPoint2D<Type> &p0, Type lk): a(lk), b(-1), c(p0.y - lk*p0.x) {}

	// copy constructor
	UpixLine(const UpixLine<Type> &other): a(other.a), b(other.b), c(other.c) {}

	// destructor
	~UpixLine(){}

	// assignment
	UpixLine<Type>& operator=(const UpixLine<Type> &rhs)
	{
		if (&rhs != this)
		{
			a = rhs.a;
			b = rhs.b;
			c = rhs.c;
		}
		return *this;
	}

public:
	inline double Distance(Type x, Type y)
	{
		return UpixAbs(a*x + b*y + c)/sqrt((double)(a*a + b*b));
	}

	inline double Distance(const UpixPoint2D<Type> &p)
	{
		return Distance(p.x, p.y);
	}

	template <class ValType>
	inline ValType X(ValType y)
	{
		assert(a != 0);
		return -(b*y + c)/a;
	}

	template <class ValType>
	inline ValType Y(ValType x)
	{
		assert(b != 0);
		return -(a*x +c)/b;
	}

	inline Type Slope()
	{
		if (b == 0)
			return (Type)UPIX_REAL_MAX;
		return -a/b;
	}

	inline Type InterceptY()
	{
		if (b == 0)
			return (Type)UPIX_REAL_MAX;
		return -c/b;
	}

	inline Type InterceptX()
	{
		if (a == 0)
			return (Type)UPIX_REAL_MAX;
		return -c/a;
	}

	inline void Set(Type la, Type lb, Type lc) { a = la; b = lb; c = lc; }

	inline void Set(Type const vl[3]) { a = vl[0]; b = vl[1]; c = vl[2]; }

	inline bool IsIdeal(Type tol = (Type)0) const
	{
		return (UpixAbs(a)<=tol*UpixAbs(c)) && (UpixAbs(b)<=tol*UpixAbs(c) && c!=0);
	}
};



//  +-+-+ UpixLine simple I/O +-+-+

//: Write "<UpixLine a,b,c>" to stream
template <class Type>
ostream&  operator<<(ostream& s, UpixLine<Type> const& line)
{
	return s << line.a << " " << line.b << " " << line.c;
}

//: Read a b c from stream
template <class Type>
istream&  operator>>(istream& is,  UpixLine<Type>& line)
{
	is >> line.a >> line.b >> line.c;
	return is;
}

// Construct point from two lines.
template <class Type>
void UPIX_ALGORITHM_TEMPLATE InterceptPoint(UpixLine<Type> const &l0, UpixLine<Type> const &l1, UpixPoint2D<Type> &p)
{
	Type w = l0.a*l1.b - l0.b*l1.a;
	assert(w != 0);
	p.x = (l0.b*l1.c - l0.c*l1.b)/w;
	p.y = (l0.c*l1.a - l0.a*l1.c)/w;
}

template <class Type>
void UPIX_ALGORITHM_TEMPLATE InterceptPoint(UpixLine<Type> const &l0, UpixLine<Type> const &l1, UpixPointHomg2D<Type> &p)
{
	p.x = l0.b*l1.c - l0.c*l1.b;
	p.y = l0.c*l1.a - l0.a*l1.c;
	p.w = l0.a*l1.b - l0.b*l1.a;
}


#endif
