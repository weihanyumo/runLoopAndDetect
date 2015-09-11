/*
* Filename: UpixRect.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_RECT_H_
#define _UPIX_RECT_H_


#include "UpixAlgorithmTypeDef.h"
#include "UpixPoint.h"
#include <math.h>


template<class T>
class UPIX_ALGORITHM_TEMPLATE UpixRect
{
public:
	UpixRect(T _x = T(0), T _y = T(0), T _width = T(0), T _height = T(0))
		: x(_x), y(_y), width(_width), height(_height)
	{

	}

	UpixRect(const UpixPoint2D<T> &p1, const UpixPoint2D<T> &p2)
	{
		x = p1.x;
		y = p1.y;
		width = p2.x - p1.x;
		height = p2.y - p1.y;
	}

	UpixRect(const UpixRect &other)
	{
		x = other.x;
		y = other.y;
		width = other.width;
		height = other.height;
	}

	virtual ~UpixRect() {}

public:
	UpixRect &operator=(const UpixRect &rhs)
	{
		if (this != &rhs)
		{
			x = rhs.x;
			y = rhs.y;
			width = rhs.width;
			height = rhs.height;
		}
		return *this;
	}

	inline T Right() const
	{
		return x + width;
	}

	inline T Bottom() const
	{
		return y + height;
	}

	UpixPoint2D<T> Center() const
	{
		return UpixPoint2D<T>(x + width/2, y + height/2);
	}

	T Area()
	{
		return width*height;
	}

	template <class ValType>
	bool IsInRect(const UpixPoint2D<ValType> &pt)
	{
		return (x <= pt.x) && (pt.x <= x + width - 1) && 
				(y <= pt.y) && (pt.y <= y + height - 1);
	}

	bool IsOverlap(const UpixRect<T> &rc)
	{
		return (x <= rc.x + rc.width - 1) && (rc.x <= x + width - 1) && 
				(y <= rc.y + rc.height - 1) && (rc.y <= y + height - 1);
	}

public:
	T x;
	T y;
	T width;
	T height;
};


typedef UpixRect<int> UpixRecti;
typedef UpixRect<float> UpixRectf;
typedef UpixRect<double> UpixRectd;



template<class T>
class UPIX_ALGORITHM_TEMPLATE UpixRotRect
{
public:
	UpixRotRect(T _cx = T(0), T _cy = T(0), T _width = T(0), T _height = T(0), double _angle = T(0))
		: cx(_cx), cy(_cy), width(_width), height(_height), angle(_angle)
	{
		CalcRotate();
	}

	UpixRotRect(const UpixRotRect &other)
	{
		cx = other.cx;
		cy = other.cy;
		width = other.width;
		height = other.height;
		angle = other.angle;
		Fill(rot, 4, other.rot);
	}

	~UpixRotRect() {}

public:
	UpixRotRect &operator=(const UpixRotRect &rhs)
	{
		if (this != &rhs)
		{
			cx = rhs.cx;
			cy = rhs.cy;
			width = rhs.width;
			height = rhs.height;
			angle = rhs.angle;
			Fill(rot, 4, rhs.rot);
		}
		return *this;
	}

	T Area()
	{
		return width*height;
	}

	// Clock wise
	template <class ValType>
	void Vertex(UpixPoint2D<ValType> pts[4]) const
	{
		double halfWidth = width*0.5;
		double halfHeight = height*0.5;

		pts[0].x = (ValType)(cx - halfWidth*rot[0] - halfHeight*rot[1]);
		pts[0].y = (ValType)(cy - halfWidth*rot[2] - halfHeight*rot[3]);
		pts[1].x = (ValType)(cx + halfWidth*rot[0] - halfHeight*rot[1]);
		pts[1].y = (ValType)(cy + halfWidth*rot[2] - halfHeight*rot[3]);
		pts[2].x = (ValType)(cx + halfWidth*rot[0] + halfHeight*rot[1]);
		pts[2].y = (ValType)(cy + halfWidth*rot[2] + halfHeight*rot[3]);
		pts[3].x = (ValType)(cx - halfWidth*rot[0] + halfHeight*rot[1]);
		pts[3].y = (ValType)(cy - halfWidth*rot[2] + halfHeight*rot[3]);
	}

	template <class ValType>
	void GetRect(UpixRect<ValType> &rc) const
	{
		UpixPoint2D<ValType> pts[4];
		Vertex(pts);

		ValType x_max = pts[0].x;
		ValType y_max = pts[0].y;
		ValType x_min = pts[0].x;
		ValType y_min = pts[0].y;
		for (int i=1; i<4; i++)
		{
			UpixPoint2D<ValType> &p = pts[i];
			if (x_max < p.x)
			{
				x_max = p.x;
			}
			else if (x_min > p.x)
			{
				x_min = p.x;
			}

			if (y_max < p.y)
			{
				y_max = p.y;
			}
			else if (y_min > p.y)
			{
				y_min = p.y;
			}
		}

		rc.x = x_min;
		rc.y = y_min;
		rc.width = x_max - x_min + 1;
		rc.height = y_max - y_min + 1;
	}

	template <class ValType>
	bool IsInRect(const UpixPoint2D<ValType> &pt) const
	{
		UpixPoint2d p;
		p.x = (pt.x - cx)*rot[0] + (pt.y - cy)*rot[2];
		p.y = (pt.x - cx)*rot[1] + (pt.y - cy)*rot[3];

		double halfWidth = width*0.5;
		double halfHeight = height*0.5;
		return (-halfWidth <= p.x) && (p.x <= halfWidth) && 
				(-halfHeight <= p.y) && (p.y <= halfHeight);
	}

protected:
	void CalcRotate()
	{
		rot[0] = cos(angle);
		rot[1] = sin(angle);
		rot[2] = -rot[1];
		rot[3] = rot[0];
	}

public:
	T cx;
	T cy;
	T width;
	T height;
	double angle;
	double rot[4];
};


typedef UpixRotRect<int> UpixRotRecti;
typedef UpixRotRect<float> UpixRotRectf;
typedef UpixRotRect<double> UpixRotRectd;


#endif