/*
* Filename: UpixPoint.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_POINT_H_
#define _UPIX_POINT_H_


#include "UpixAlgorithmTypeDef.h"
#include "UpixArray.h"
#include <math.h>



#define UPIX_POINT_TOLERANCE	1e-6



//--------------------------------------UpixPoint2D----------------------------------------------
//: Represents a cartesian 2D point
template <class Type>
class UPIX_ALGORITHM_TEMPLATE UpixPoint2D
{
public:
	// the data associated with this point
	Type x;
	Type y;

public:
	// Constructors/Initializers/Destructor------------------------------------
	//: Construct from two Types.
	UpixPoint2D (Type px = Type(0), Type py = Type(0)) : x(px), y(py) {}

	//: Construct from 2-array.
	UpixPoint2D (const Type v[2]) : x(v[0]), y(v[1]) {}

	// Default copy constructor
	UpixPoint2D(const UpixPoint2D<Type> &p) : x(p.x), y(p.y) {}

	//: Destructor
	~UpixPoint2D () {}

	//: Assignment
	UpixPoint2D<Type> &operator=(const UpixPoint2D<Type> &p)
	{ 
		if (&p != this)
		{
			x = p.x; 
			y = p.y; 
		}
		return *this; 
	}

public:
	//: Test for equality
	bool operator==(const UpixPoint2D<Type> &p) const 
	{ 
		return (this == &p) || (UpixAbs(x - p.x) < UPIX_POINT_TOLERANCE && 
								UpixAbs(y - p.y) < UPIX_POINT_TOLERANCE); 
	}

	//: Set point
	inline void Set (Type px, Type py) { x = px; y = py; }

	inline void Set (const Type p[2]) { x = p[0]; y = p[1]; }

	//  +-+-+ point_2d arithmetic +-+-+

	// Dot product
	inline Type Dot(const UpixPoint2D<Type> &p)
	{
		return x*p.x + y*p.y;
	}

	//: The difference of two points
	UpixPoint2D<Type> operator-(const UpixPoint2D<Type> &p) const
	{ 
		return UpixPoint2D<Type>(x - p.x, y - p.y); 
	}

	//: Adding a vector to a point gives a new point at the end of that vector
	UpixPoint2D<Type> operator+(const UpixPoint2D<Type> &p) const
	{ 
		return UpixPoint2D<Type>(x + p.x, y + p.y); 
	}

	// Multiply a number
	UpixPoint2D<Type> operator*(const Type &val) const
	{
		return UpixPoint2D<Type>(x*val, y*val);
	}

	// Divide a number
	UpixPoint2D<Type> operator/(const Type &val) const
	{
		assert(val != 0);
		return UpixPoint2D<Type>(x/val, y/val);
	}

	//: Adding a vector to a point gives the point at the end of that vector
	UpixPoint2D<Type> operator+=(const UpixPoint2D<Type> &p)
	{
		x += p.x;
		y += p.y; 
		return *this; 
	}

	//: Subtracting a vector from a point is the same as adding the inverse vector
	UpixPoint2D<Type> operator-=(const UpixPoint2D<Type> &p)
	{ 
		x -= p.x;
		y -= p.y; 
		return *this; 
	}

	// Multiply a number from the point
	UpixPoint2D<Type> operator*=(const Type &val)
	{
		x *= val;
		y *= val;
		return *this;
	}

	// Divide a number from the point
	UpixPoint2D<Type> operator/=(const Type &val)
	{
		assert(val != 0);
		x /= val;
		y /= val;
		return *this;
	}

	UpixPoint2D<Type> operator-() const
	{
		return UpixPoint2D<Type>(-x, -y);
	}
};


typedef UpixPoint2D<int> UpixPoint2i;
typedef UpixPoint2D<float> UpixPoint2f;
typedef UpixPoint2D<double> UpixPoint2d;



//--------------------------------------UpixPoint3D------------------------------------------
template <class Type>
class UPIX_ALGORITHM_TEMPLATE UpixPoint3D
{
public:
	// the data associated with this point
	Type x;
	Type y;
	Type z;

public:
	// Constructors/Initializers/Destructor------------------------------------
	//: Construct from three Types.
	UpixPoint3D(Type px = Type(0), Type py = Type(0), Type pz = Type(0)) : x(px), y(py), z(pz) {}

	//: Construct from 3-array.
	UpixPoint3D(const Type v[3]) : x(v[0]), y(v[1]), z(v[2]) {}

	//: Copy constructor
	UpixPoint3D(const UpixPoint3D<Type> &p) : x(p.x), y(p.y), z(p.z) {}

	//: Destructor
	~UpixPoint3D () {}

	//: Assignment
	UpixPoint3D<Type> &operator=(const UpixPoint3D<Type> &p)
	{ 
		if (&p != this)
		{
			x = p.x;
			y = p.y;
			z = p.z;
		}
		return *this; 
	}

	//: Test for equality
	inline bool operator==(const UpixPoint3D<Type> &p) const
	{ 
		return (this == &p) || (UpixAbs(x - p.x) < UPIX_POINT_TOLERANCE &&
								UpixAbs(y - p.y) < UPIX_POINT_TOLERANCE &&
								UpixAbs(z - p.z) < UPIX_POINT_TOLERANCE); 
	}

	//: Set point
	inline void Set (Type px, Type py, Type pz) { x = px; y = py; z = pz; }

	inline void Set (Type const p[3]) { x = p[0]; y = p[1]; z = p[2]; }

	//  +-+-+ point_3d arithmetic +-+-+

	// Dot product
	inline Type Dot(const UpixPoint3D<Type> &p) const
	{
		return x*p.x + y*p.y + z*p.z;
	}

	// Cross product
	UpixPoint3D<Type> Cross(const UpixPoint3D<Type> &p) const
	{
		UpixPoint3D<Type> result;
		result.x = y*p.z - z*p.y;
		result.y = z*p.x - x*p.z;
		result.z = x*p.y - y*p.x;
		return result;
	}

	// Two normal
	Type TwoNorm() const
	{
		return sqrt(x*x + y*y + z*z);
	}

	// Normalize this vector
	void Normalize()
	{
		Type invd = (Type)(1.0/TwoNorm());
		x *= invd;
		y *= invd;
		z *= invd;
	}

	//: The difference of two points is the vector from second to first point
	inline UpixPoint3D<Type> operator-(const UpixPoint3D<Type> &p) const
	{ 
		return UpixPoint3D<Type>(x - p.x, y - p.y, z - p.z);
	}

	//: Adding a vector to a point gives a new point at the end of that vector
	inline UpixPoint3D<Type> operator+(const UpixPoint3D<Type> &p) const
	{ 
		return UpixPoint3D<Type>(x + p.x, y + p.y, z + p.z); 
	}

	// Multiply a number
	inline UpixPoint3D<Type> operator*(const Type &val) const
	{
		return UpixPoint3D<Type>(x*val, y*val, z*val);
	}

	// Divide a number
	inline UpixPoint3D<Type> operator/(const Type &val) const
	{
		assert(val != 0);
		Type s = 1/val;
		return UpixPoint3D<Type>(x*s, y*s, z*s);
	}

	//: Adding a vector to a point gives the point at the end of that vector
	inline UpixPoint3D<Type> operator+=(const UpixPoint3D<Type>& p)
	{ 
		x += p.x;
		y += p.y;
		z += p.z;
		return *this; 
	}

	//: Subtracting a vector from a point is the same as adding the inverse vector
	inline UpixPoint3D<Type> operator-=(const UpixPoint3D<Type> & p)
	{ 
		x -= p.x;
		y -= p.y;
		z -= p.z;
		return *this;
	}

	// Multiply a number from the point
	inline UpixPoint3D<Type> operator*=(const Type &val)
	{
		x *= val;
		y *= val;
		z *= val;
		return *this;
	}

	// Divide a number from the point
	inline UpixPoint3D<Type> operator/=(const Type &val)
	{
		assert(val != 0);
		Type s = 1/val;
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}

	UpixPoint3D<Type> operator-() const
	{
		return UpixPoint3D<Type>(-x, -y, -z);
	}
};


typedef UpixPoint3D<int> UpixPoint3i;
typedef UpixPoint3D<float> UpixPoint3f;
typedef UpixPoint3D<double> UpixPoint3d;



//-------------------------------------UpixPointHomg2D------------------------------------------
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixPointHomg2D
{
public:
	// the data associated with this point
	T x;
	T y;
	T w;

public:
	// Constructors/Initializers/Destructor------------------------------------
	//: Construct from two (nonhomogeneous) or three (homogeneous) Types.
	UpixPointHomg2D(T px = T(0), T py = T(0), T pw = T(1)) : x(px), y(py), w(pw) {}

	//: Construct from homogeneous 3-array.
	UpixPointHomg2D(const T v[3]) : x(v[0]), y(v[1]), w(v[2]) {}

	// Default copy constructor
	UpixPointHomg2D(const UpixPointHomg2D<T> &p) : x(p.x), y(p.y), w(p.w) {}

	// Destructor
	~UpixPointHomg2D() {}

	// Default assignment operator
	UpixPointHomg2D<T>& operator=(const UpixPointHomg2D<T> &p)
	{
		if (&p != this)
		{
			x = p.x;
			y = p.y;
			w = p.w;
		}
		return *this;
	}

	//: the comparison operator
	inline bool operator==(UpixPointHomg2D<T> const& p) const
	{
		return (this == &p) || (UpixAbs(x*p.w - w*p.x) < UPIX_POINT_TOLERANCE && 
								UpixAbs(y*p.w - w*p.y) < UPIX_POINT_TOLERANCE && 
								UpixAbs(y*p.x - x*p.y) < UPIX_POINT_TOLERANCE);
	}

	//: Set \a x,y,w
	inline void Set(T px, T py, T pw = T(1)) { x = px, y = py, w = pw; }

	inline void Set(T const p[3]) { x = p[0]; y = p[1]; w = p[2]; }

	bool Rescale(T new_w = T(1))
	{
		if (w == 0)
			return false;

		T s = new_w/w;
		x *= s;
		y *= s;
		w = new_w;

		return true;
	}

	//  +-+-+ homg_point_2d arithmetic +-+-+

	//: The difference of two points is the vector from second to first point
	// This function is only valid if the points are not at infinity.
	inline UpixPointHomg2D<T> operator-(UpixPointHomg2D<T> const& p)
	{
		assert(p.w && w);
		return UpixPointHomg2D<T>(x/w-p.x/p.w, y/w-p.y/p.w, 1);
	}

	inline UpixPointHomg2D<T> operator+(UpixPointHomg2D<T> const& p)
	{
		assert(p.w && w);
		return UpixPointHomg2D<T>(x/w+p.x/p.w, y/w+p.y/p.w, 1);
	}

	//: Adding a vector to a point gives the point at the end of that vector
	// If the point is at infinity, nothing happens.
	inline UpixPointHomg2D<T> operator+=(UpixPointHomg2D<T>& p)
	{
		assert(p.w && w);
		x = x/w + p.x/p.w;
		y = y/w + p.y/p.w;
		w = 1;
		return *this; 
	}

	//: Subtracting a vector from a point is the same as adding the inverse vector
	inline UpixPointHomg2D<T> operator-=(UpixPointHomg2D<T>& p)
	{ 
		assert(p.w && w);
		x = x/w - p.x/p.w;
		y = y/w - p.y/p.w;
		w = 1;
		return *this;  
	}

	//: Return true iff the point is at infinity (an ideal point).
	// The method checks whether |w| <= tol * max(|x|,|y|)
	inline bool IsIdeal(T tol = T(0)) const
	{
		return (UpixAbs(w) <= tol * UpixAbs(x)) || (UpixAbs(w) <= tol * UpixAbs(y));
	}
};


typedef UpixPointHomg2D<int> UpixPointHomg2i;
typedef UpixPointHomg2D<float> UpixPointHomg2f;
typedef UpixPointHomg2D<double> UpixPointHomg2d;




//----------------------------------------UpixPointHomg3D---------------------------------------
//: Represents a homogeneous 3D point
template <class Type>
class UPIX_ALGORITHM_TEMPLATE UpixPointHomg3D
{
public:
	// the data associated with this point
	Type x;
	Type y;
	Type z;
	Type w;

public:
	// Constructors/Initializers/Destructor------------------------------------
	//: Construct from three (nonhomogeneous) or four (homogeneous) Types.
	UpixPointHomg3D(Type px = Type(0), Type py = Type(0), Type pz = Type(0), Type pw = (Type)1) : x(px), y(py), z(pz), w(pw) {}

	//: Construct from homogeneous 4-array.
	UpixPointHomg3D(const Type v[4]) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) {}

	// Default copy constructor
	UpixPointHomg3D(const UpixPointHomg3D<Type>& that) : x(that.x), y(that.y), z(that.z), w(that.w) {}

	// Destructor
	~UpixPointHomg3D() {}

	// Default assignment operator
	UpixPointHomg3D<Type> &operator=(const UpixPointHomg3D<Type> &that)
	{
		if (&that != this)
		{
			x = that.x;
			y = that.y;
			z = that.z;
			w = that.w;
		}
		return *this;
	}

	//: the comparison operator
	bool operator==(UpixPointHomg3D<Type> const& p) const
	{
		return (this == &p) || (UpixAbs(x*p.y - y*p.x) < UPIX_POINT_TOLERANCE &&
								UpixAbs(x*p.z - z*p.x) < UPIX_POINT_TOLERANCE &&
								UpixAbs(x*p.w - w*p.x) < UPIX_POINT_TOLERANCE &&
								UpixAbs(y*p.z - z*p.y) < UPIX_POINT_TOLERANCE &&
								UpixAbs(y*p.w - w*p.y) < UPIX_POINT_TOLERANCE &&
								UpixAbs(z*p.w - w*p.z) < UPIX_POINT_TOLERANCE);
	}

	//: Set \a x,y,z,w
	inline void Set(Type px, Type py, Type pz, Type pw = (Type)1)
	{
		x = px; y = py; z = pz; w = pw; 
	}

	inline void Set(Type const p[4]) { x=p[0]; y=p[1]; z=p[2]; w=p[3]; }

	//: Return true iff the point is at infinity (an ideal point).
	// The method checks whether |w| <= tol * max(|x|,|y|,|z|)
	inline bool IsIdeal(Type tol = (Type)0) const
	{
		return UpixAbs(w) <= tol * UpixAbs(x) || UpixAbs(w) <= tol * UpixAbs(y) || UpixAbs(w) <= tol * UpixAbs(z);
	}

	inline bool Rescale(Type new_w = Type(1))
	{
		if (w == 0)
			return false;

		Type s = new_w / w;
		x *= s;
		y *= s;
		z *= s;
		w = new_w;

		return true;
	}

	void Normalize()
	{
		Type invd = (Type)(1/sqrt(x*x + y*y + z*z + w*w));
		x *= invd;
		y *= invd;
		z *= invd;
		w *= invd;
	}

	//  +-+-+ homg_point_3d arithmetic +-+-+

	//: The difference of two points is the vector from second to first point
	// This function is only valid if the points are not at infinity.
	inline UpixPointHomg3D<Type> operator-(UpixPointHomg3D<Type> const& p)
	{
		assert(p.w && w);
		return UpixPointHomg3D<Type>(x/w-p.x/p.w, y/w-p.y/p.w, z/w-p.z/p.w, 1);
	}

	inline UpixPointHomg3D<Type> operator+(UpixPointHomg3D<Type> const& p)
	{
		assert(p.w && w);
		return UpixPointHomg3D<Type>(x/w+p.x/p.w, y/w+p.y/p.w, z/w+p.z/p.w, 1);
	}

	//: Adding a vector to a point gives the point at the end of that vector
	// If the point is at infinity, nothing happens.
	inline UpixPointHomg3D<Type> operator+=(UpixPointHomg3D<Type>& p)
	{ 
		x = x/w + p.x/p.w;
		y = y/w + p.y/p.w;
		z = z/w + p.z/p.w;
		w = 1;
		return *this;
	}

	inline UpixPointHomg3D<Type> operator-=(UpixPointHomg3D<Type>& p)
	{ 
		x = x/w - p.x/p.w;
		y = y/w - p.y/p.w;
		z = z/w - p.z/p.w;
		w = 1;
		return *this;
	}

	inline UpixPointHomg3D<Type> operator*=(const Type &val)
	{
		x *= val;
		y *= val;
		z *= val;
		w *= val;
		return *this;
	}

	inline UpixPointHomg3D<Type> operator/=(const Type &val)
	{
		assert(val != 0);
		Type s = 1/val;
		x *= s;
		y *= s;
		z *= s;
		w *= s;
		return *this;
	}
};


typedef UpixPointHomg3D<int> UpixPointHomg3i;
typedef UpixPointHomg3D<float> UpixPointHomg3f;
typedef UpixPointHomg3D<double> UpixPointHomg3d;



//--------------------- type convert ------------------------
template<class T>
inline UpixPoint2D<T> UPIX_ALGORITHM_API Dehomogenized(const UpixPointHomg2D<T> &rhs)
{
	assert(rhs.w);
	return UpixPoint2D<T>(rhs.x/rhs.w, rhs.y/rhs.w);
}

template<class T>
inline UpixPoint3D<T> UPIX_ALGORITHM_API Dehomogenized(const UpixPointHomg3D<T> &rhs)
{
	assert(rhs.w);
	T s = 1/rhs.w;
	return UpixPoint3D<T>(rhs.x*s, rhs.y*s, rhs.z*s);
}

template<class T>
inline UpixPointHomg2D<T> UPIX_ALGORITHM_API Homogenized(const UpixPoint2D<T> &rhs)
{
	return UpixPointHomg2D<T>(rhs.x, rhs.y, T(1));
}

template<class T>
inline UpixPointHomg3D<T> UPIX_ALGORITHM_API Homogenized(const UpixPoint3D<T> &rhs)
{
	return UpixPointHomg3D<T>(rhs.x, rhs.y, rhs.z, T(1));
}




//------------Normalize Data-------------

// calculate the normalizing transformations for point sets:
// after the transformation the set will have the mass center at the coordinate origin
// and the RMS distance from the origin will be ~sqrt(2).
// T = [ sx    0   -sx*x0
//       0    sy   -sy*y0
//       0    0      1 ]
template <class T>
void UPIX_ALGORITHM_TEMPLATE GetNormalTransform(const UpixArray< UpixPoint2D<T> > &p, UpixPoint2D<T> &p0, UpixPoint2D<T> &scale)
{
	p0.x = 0;
	p0.y = 0;
	scale.x = 0;
	scale.y = 0;

	// compute centers and average distances for point sets
	int count  = p.Size();
	for( int i=0; i<count; i++ )
	{
		p0.x += p[i].x; 
		p0.y += p[i].y;
	}

	double t = 1./count;
	p0.x *= t; 
	p0.y *= t;

	for( int i=0; i<count; i++ )
	{
		scale.x += (p[i].x - p0.x)*(p[i].x - p0.x);
		scale.y += (p[i].y - p0.y)*(p[i].y - p0.y);
	}
	scale.x = sqrt(scale.x*t); // std(or RMS)
	scale.y = sqrt(scale.y*t);

	if (scale.x == 0)
	{
		scale.x = 1;
	}
	if (scale.y == 0)
	{
		scale.y = 1;
	}

	scale.x = UPIX_SQRT2/scale.x;
	scale.y = UPIX_SQRT2/scale.y;
}

template <class T>
void UPIX_ALGORITHM_TEMPLATE GetNormalTransform(const UpixArray< UpixPoint2D<T> > &p, const UpixArray<int> &indices,
												UpixPoint2D<T> &p0, UpixPoint2D<T> &scale)
{
	p0.x = 0;
	p0.y = 0;
	scale.x = 0;
	scale.y = 0;

	// compute centers and average distances for point sets
	int count = indices.Size();
	for( int i=0; i<count; i++ )
	{
		p0.x += p[indices[i]].x; 
		p0.y += p[indices[i]].y;
	}

	double t = 1./count;
	p0.x *= t; 
	p0.y *= t;

	for( int i=0; i<count; i++ )
	{
		const UpixPoint2D<T> &pt = p[indices[i]];
		scale.x += (pt.x - p0.x)*(pt.x - p0.x);
		scale.y += (pt.y - p0.y)*(pt.y - p0.y);
	}
	scale.x = sqrt(scale.x*t); // std(or RMS)
	scale.y = sqrt(scale.y*t);

	if (scale.x == 0)
	{
		scale.x = 1;
	}
	if (scale.y == 0)
	{
		scale.y = 1;
	}

	scale.x = UPIX_SQRT2/scale.x;
	scale.y = UPIX_SQRT2/scale.y;
}

// pn = T * p
// xn = sx(x - x0)
// yn = sy(y - y0)
template <class T>
void UPIX_ALGORITHM_TEMPLATE NormalizePointData(const UpixArray< UpixPoint2D<T> > &p, const UpixPoint2D<T> &p0, 
												const UpixPoint2D<T> &scale, UpixArray< UpixPoint2D<T> > &p_norm)
{
	int count = p.Size();
	if (p_norm.Size() != count)
	{
		p_norm.Resize(count);
	}

	for (int i=0; i<count; i++)
	{
		p_norm[i].x = scale.x*(p[i].x - p0.x);
		p_norm[i].y = scale.y*(p[i].y - p0.y);
	}
}

template <class T>
void UPIX_ALGORITHM_TEMPLATE NormalizePointData(const UpixArray< UpixPoint2D<T> > &p, const UpixArray<int> &indices, 
												const UpixPoint2D<T> &p0, const UpixPoint2D<T> &scale, 
												UpixArray< UpixPoint2D<T> > &p_norm)
{
	int count = indices.Size();
	if (p_norm.Size() != count)
	{
		p_norm.Resize(count);
	}

	for (int i=0; i<count; i++)
	{
		p_norm[i].x = scale.x*(p[indices[i]].x - p0.x);
		p_norm[i].y = scale.y*(p[indices[i]].y - p0.y);
	}
}


//------------- cross product for point 3d -------------
template <class T>
inline UpixPoint3D<T> UPIX_ALGORITHM_TEMPLATE Cross(const UpixPoint3D<T> &p1, const UpixPoint3D<T> &p2)
{
	return p1.Cross(p2);
}

// Dot product
template <class T>
inline T UPIX_ALGORITHM_TEMPLATE Dot(const UpixPoint2D<T> &p1, const UpixPoint2D<T> &p2)
{
	return p1.Dot(p2);
}

template <class T>
inline T UPIX_ALGORITHM_TEMPLATE Dot(const UpixPoint3D<T> &p1, const UpixPoint3D<T> &p2)
{
	return p1.Dot(p2);
}

// Distance between 2 points
template <class T>
double UPIX_ALGORITHM_TEMPLATE Distance(const UpixPoint2D<T> &p1, const UpixPoint2D<T> &p2)
{
	return sqrt((double)((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)));
}

template <class T>
double UPIX_ALGORITHM_TEMPLATE Distance(const UpixPoint3D<T> &p1, const UpixPoint3D<T> &p2)
{
	return sqrt((double)((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z)));
}

// Angle between 2 vector
template <class T>
double UPIX_ALGORITHM_TEMPLATE Angle(const UpixPoint3D<T> &p1, const UpixPoint3D<T> &p2)
{
	// cos(a) = (v1, v2)/(|v1||v2|)
	double cosA = p1.Dot(p2)/(p1.TwoNorm()*p2.TwoNorm());
	return acos(cosA)*180/UPIX_PI;
}


// ----------- Is point in triangle -----------
// If a point P is in triangle ABC, there must be equation as follow:
// P = aA + bB + cC, and to the non-negative a,b,c we have a + b + c = 1
// so P = (1 - b - c)A + bB + cC  =>  P = A +  b(B 每 A) + c(C - A) 
// let v0 = B 每 A, v1 = C 每 A, v2 = P 每 A
// we have v2 = bv0 + cv1, then:
// x2 = bx0 + cx1
// y2 = by0 + cy1
// from the above equations we can find:
// b = (x2*y1 - x1*y2) / (x0*y1 - x1*y0)
// c = (x0*y2 - x2*y0) / (x0*y1 - x1*y0)
template <class T>
bool IsPointInTriangle(const UpixPoint2D<T> &pA, const UpixPoint2D<T> &pB, const UpixPoint2D<T> &pC, const UpixPoint2D<T> &p)
{
	UpixPoint2D<T> v0 = pB - pA ;
	UpixPoint2D<T> v1 = pC - pA ;
	UpixPoint2D<T> v2 = p - pA ;
	double inverDeno = 1.0/(v0.x*v1.y - v1.x*v0.y);

	double b = (v2.x*v1.y - v1.x*v2.y)*inverDeno;
	if (b<0 || b>1)
	{
		return false;
	}

	double c = (v0.x*v2.y - v2.x*v0.y)*inverDeno;
	if (c<0 || c>1)
	{
		return false;
	}

	return b + c <= 1;
}

template <class T>
bool IsPointInConvex(const UpixArray< UpixPoint2D<T> > &p2d, const UpixPoint2D<T> &p)
{
	int nTriangleCount = p2d.Size() - 2;
	for (int i=0; i<nTriangleCount; i++)
	{
		if (IsPointInTriangle(p2d[0], p2d[i + 1], p2d[i + 2], p))
		{
			return true;
		}
	}
	return false;
}



// ---------- Interpolate ----------
template <class T>
UpixPoint2D<T> UPIX_ALGORITHM_TEMPLATE Interpolate(const UpixPoint2D<T> &p0, const UpixPoint2D<T> &p1, const T &percent)
{
	return UpixPoint2D<T>(UpixInterp(p0.x, p1.x, percent), UpixInterp(p0.y, p1.y, percent));
}

template <class T>
UpixPoint3D<T> UPIX_ALGORITHM_TEMPLATE Interpolate(const UpixPoint3D<T> &p0, const UpixPoint3D<T> &p1, const T &percent)
{
	return UpixPoint3D<T>(UpixInterp(p0.x, p1.x, percent), UpixInterp(p0.y, p1.y, percent), UpixInterp(p0.z, p1.z, percent));
}



#endif