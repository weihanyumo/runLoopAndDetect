/*
* Filename: UpixSize.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_SIZE_H_
#define _UPIX_SIZE_H_

#include "UpixAlgorithmTypeDef.h"


#define UPIX_SIZE_TOLERANCE	1e-6


// Specifies a 2 dimensional size.
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixSize
{
public:
	//! Default constructor for empty dimension
	UpixSize() : Width(0), Height(0) {}

	//! Constructor with width and height
	UpixSize(const T& width, const T& height) : Width(width), Height(height) {}

	//! Use this constructor only where you are sure that the conversion is valid.
	template <class U>
	explicit UpixSize(const UpixSize<U>& other) :
		Width((T)other.Width), Height((T)other.Height) {}

	template <class U>
	UpixSize<T>& operator=(const UpixSize<U>& other)
	{ 
		Width = (T) other.Width;
		Height = (T) other.Height;
		return *this;
	}

	//! Equality operator
	bool operator==(const UpixSize<T>& other) const
	{
		return UpixAbs(Width - other.Width)<UPIX_SIZE_TOLERANCE &&
				UpixAbs(Height - other.Height)<UPIX_SIZE_TOLERANCE;
	}

	//! Inequality operator
	bool operator!=(const UpixSize<T>& other) const
	{
		return ! (*this == other);
	}

	//! Set to new values
	UpixSize<T>& Set(const T& width, const T& height)
	{
		Width = width;
		Height = height;
		return *this;
	}

	//! Divide width and height by scalar
	UpixSize<T>& operator/=(const T& scale)
	{
		Width /= scale;
		Height /= scale;
		return *this;
	}

	//! Divide width and height by scalar
	UpixSize<T> operator/(const T& scale) const
	{
		return UpixSize<T>(Width/scale, Height/scale);
	}

	//! Multiply width and height by scalar
	UpixSize<T>& operator*=(const T& scale)
	{
		Width *= scale;
		Height *= scale;
		return *this;
	}

	//! Multiply width and height by scalar
	UpixSize<T> operator*(const T& scale) const
	{
		return UpixSize<T>(Width*scale, Height*scale);
	}

	//! Add another dimension to this one.
	UpixSize<T>& operator+=(const UpixSize<T>& other)
	{
		Width += other.Width;
		Height += other.Height;
		return *this;
	}

	//! Subtract a dimension from this one
	UpixSize<T>& operator-=(const UpixSize<T>& other)
	{
		Width -= other.Width;
		Height -= other.Height;
		return *this;
	}


	//! Add two dimensions
	UpixSize<T> operator+(const UpixSize<T>& other) const
	{
		return UpixSize<T>(Width+other.Width, Height+other.Height);
	}

	//! Get area
	T GetArea() const
	{
		return Width*Height;
	}

	//! Get the optimal size according to some properties
	/** This is a function often used for texture dimension
	calculations. The function returns the next larger or
	smaller dimension which is a power-of-two dimension
	(2^n,2^m) and/or square (Width=Height).
	\param requirePowerOfTwo Forces the result to use only
	powers of two as values.
	\param requireSquare Makes width==height in the result
	\param larger Choose whether the result is larger or
	smaller than the current dimension. If one dimension
	need not be changed it is kept with any value of larger.
	\param maxValue Maximum texturesize. if value > 0 size is
	clamped to maxValue
	\return The optimal dimension under the given
	constraints. */
	UpixSize<T> GetOptimalSize(
			bool requirePowerOfTwo = true,
			bool requireSquare = false,
			bool larger = true,
			uint maxValue = 0) const
	{
		uint i=1;
		uint j=1;
		if (requirePowerOfTwo)
		{
			while (i<(uint)Width)
				i<<=1;
			if (!larger && i!=1 && i!=(uint)Width)
				i>>=1;
			while (j<(uint)Height)
				j<<=1;
			if (!larger && j!=1 && j!=(uint)Height)
				j>>=1;
		}
		else
		{
			i=(uint)Width;
			j=(uint)Height;
		}

		if (requireSquare)
		{
			if ((larger && (i>j)) || (!larger && (i<j)))
				j=i;
			else
				i=j;
		}

		if ( maxValue > 0 && i > maxValue)
			i = maxValue;

		if ( maxValue > 0 && j > maxValue)
			j = maxValue;

		return UpixSize<T>((T)i, (T)j);
	}

	//! Get the interpolated dimension
	/** \param other Other dimension to interpolate with.
	\param d Value between 0.0f and 1.0f.
	\return Interpolated dimension. */
	UpixSize<T> GetInterpolated(const UpixSize<T>& other, float d) const
	{
		return UpixSize<T>(UpixInterp(other.Width, Width, d), UpixInterp(other.Height, Height, d));
	}

public:
	T Width;
	T Height;
};


typedef UpixSize<double> UpixSized;
typedef UpixSize<float> UpixSizef;
typedef UpixSize<uint> UpixSizeui;

//! Typedef for an integer dimension.
/** There are few cases where negative dimensions make sense. Please consider using
	AxSizeui instead. */
typedef UpixSize<int> UpixSizei;


#endif

