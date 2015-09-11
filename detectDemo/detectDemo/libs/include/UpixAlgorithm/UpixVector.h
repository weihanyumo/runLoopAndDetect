/*
* Filename: UpixVector.h
* Author:   Xiehong
* Date:     2012.9.29
*/
// This class is not a container, only base types are supported.

#ifndef _UPIX_VECTOR_H_
#define _UPIX_VECTOR_H_

#include "UpixAlgorithmTypeDef.h"
#include "UpixMemory.h"
#include <math.h>


template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixVector
{
public:
	//default constructor
	UpixVector() { m_pData = NULL; N = 0; }

	//construct by size
	UpixVector(uint n);

	//construct from a value
	UpixVector(uint n, const T &value);

	//construct from array
	UpixVector(uint n, T value[]);

	//construct from other
	UpixVector(const UpixVector &other);

	//operator =
	UpixVector& operator=(const UpixVector &rhs);

	// destructor
	~UpixVector() { Clear(); }

public:
	// get size
	uint Size() const { return N; }

	// change size
	void Resize(uint n);

	// clear buffer
	void Clear();

	// get data array
	T *GetData() const { return m_pData; }

	// set date by array
	void SetData(T *pData);

	// set all data to 0
	inline void SetZeros(uint nSize) { Resize(nSize); }

public:
	//the one normal
	T OneNorm();

	//the two normal
	T TwoNorm() const;

	//normalization(向量本身被修改)
	void Normalize();

	//normalization(结果给新向量)
	void Normalize(UpixVector &vec);

	//the sum
	T Sum();

	//mean value
	inline T Mean() { return Sum()/N; }

	//Standard deviation
	T StdDeviation();

	// Total deviation
	T Deviation();

	//minimum value
	T MinValue();

	//maximum value
	T MaxValue();

	//cross product
	UpixVector Cross(const UpixVector &other) const;

	//dot product
	T Dot(const UpixVector &other);

public:
	//operator []
	T& operator[](uint i) { return m_pData[i]; }

	const T& operator[](uint i) const { return m_pData[i]; }

	// operator ()
	T& operator()(uint i) { return m_pData[i]; }

	const T& operator()(uint i) const { return m_pData[i]; }

	//operator ==
	bool operator==(const UpixVector &other) const;

	//operator !=
	bool operator!=(const UpixVector &other) const { return !(*this == other); }

	//operator +
	UpixVector operator+(const UpixVector &other) const;

	// operator +
	UpixVector operator+(const T &value) const;

	//operator -
	UpixVector operator-(const UpixVector &other) const;

	// operator -
	UpixVector operator-(const T &value) const;

	//operator *
	UpixVector operator*(const T &value) const;

	//operator /
	UpixVector operator/(const T &value) const;

	//operator +=
	UpixVector &operator+=(const UpixVector &other);

	//operator -=
	UpixVector &operator-=(const UpixVector &other);

	// operator +=
	UpixVector<T> &operator+=(const T &value);

	// operator -=
	UpixVector<T> &operator-=(const T &value);

	//operator *=
	UpixVector &operator*=(const T &value);

	//operator /=
	UpixVector &operator/=(const T &value);

	//change sign
	UpixVector operator-();

protected:
	uint N;         //长度
	T* m_pData;     //数据缓冲区
};


//--------------------Implement of vector---------------------
template <class T>
UpixVector<T>::UpixVector(uint n)
{
	N = n;
	m_pData = new T[N];
	int val = 0;
	Fill(m_pData, N, val);
}

template <class T>
UpixVector<T>::UpixVector(uint n, const T &value)
{
	N = n;
	m_pData = new T[N];
	Fill(m_pData, N, value);
}

template <class T>
UpixVector<T>::UpixVector(uint n, T value[])
{
	N = n;
	m_pData = new T[N];
	Fill(m_pData, N, value);
}

template <class T>
UpixVector<T>::UpixVector(const UpixVector<T> &other)
{
	N = other.N;
	m_pData = new T[N];
	Fill(m_pData, N, other.m_pData);
}

template <class T>
UpixVector<T> &UpixVector<T>::operator=(const UpixVector<T> &rhs)
{
	if (&rhs != this)
	{
		if (N != rhs.N)
		{
			Clear();
			N = rhs.N;
			m_pData = new T[N];
		}
		Fill(m_pData, N, rhs.m_pData);
	}
	return *this;
}

template <class T>
void UpixVector<T>::Clear()
{
	if (m_pData != NULL)
	{
		delete [] m_pData;
		m_pData = NULL;
	}
	N = 0;
}

template <class T>
void UpixVector<T>::Resize(uint n)
{
	Clear();
	N = n;
	m_pData = new T[N];
	int val = 0;
	Fill(m_pData, N, val);
}

template <class T>
T UpixVector<T>::OneNorm()
{
	T tonenorm = 0;
	for (uint i=0; i<N; i++)
	{
		tonenorm += UpixAbs(m_pData[i]);
	}
	return tonenorm;
}

template <class T>
T UpixVector<T>::TwoNorm() const
{
	T ttwonorm = 0;
	for (uint i=0; i<N; i++)
	{
		ttwonorm += m_pData[i] * m_pData[i];
	}
	return sqrt(ttwonorm);
}

template <class T>
void UpixVector<T>::Normalize()
{
	T s = 1/TwoNorm();
	for (uint i=0; i<N; i++)
	{
		m_pData[i] *= s;
	}
}

template <class T>
inline void UpixVector<T>::Normalize(UpixVector<T> &vec)
{
	vec = *this;
	vec.Normalize();
}

template <class T>
T UpixVector<T>::Sum()
{
	T tsum = 0;
	for (uint i=0; i<N; i++)
	{
		tsum += m_pData[i];
	}
	return tsum;
}

template <class T>
T UpixVector<T>::StdDeviation()
{
	T aver = Mean();
	T deviation = 0;
	for (uint i=0; i<N; i++)
	{
		deviation += (m_pData[i] - aver)*(m_pData[i] - aver);
	}
	return sqrt(deviation/(N - 1));
}

template <class T>
T UpixVector<T>::Deviation()
{
	T aver = Mean();
	T deviation = 0;
	for (uint i=0; i<N; i++)
	{
		deviation += (m_pData[i] - aver)*(m_pData[i] - aver);
	}
	return sqrt(deviation/N);
}

template <class T>
void UpixVector<T>::SetData(T *pData)
{
	Fill(m_pData, N, pData);
}

template <class T>
T UpixVector<T>::MinValue()
{
	T min_data = m_pData[0];
	for (uint i=1; i<N; i++)
	{
		if (m_pData[i] < min_data)
		{
			min_data = m_pData[i];
		}
	}
	return min_data;
}

template <class T>
T UpixVector<T>::MaxValue()
{
	T max_data = m_pData[0];
	for (uint i=1; i<N; i++)
	{
		if (m_pData[i] > max_data)
		{
			max_data = m_pData[i];
		}
	}
	return max_data;
}

template <class T>
UpixVector<T> UpixVector<T>::Cross(const UpixVector<T> &other) const
{
	assert((N == 3) && (other.N == 3));
	UpixVector<T> result(3);
	T *pV = other.m_pData;
	result[0] = m_pData[1]*pV[2] - m_pData[2]*pV[1];
	result[1] = m_pData[2]*pV[0] - m_pData[0]*pV[2];
	result[2] = m_pData[0]*pV[1] - m_pData[1]*pV[0];
	return result;
}

template <class T>
T UpixVector<T>::Dot(const UpixVector<T> &other)
{
	assert(N == other.N);
	T result = 0;
	for (uint i=0; i<N; i++)
	{
		result += m_pData[i]*other.m_pData[i];
	}
	return result;
}

template <class T>
bool UpixVector<T>::operator==(const UpixVector<T> &other) const
{
	if (N != other.N) return false;
	for (uint i=0; i<N; i++)
		if (m_pData[i] != other.m_pData[i])
			return false;
	return true;
}

template <class T>
UpixVector<T> UpixVector<T>::operator+(const UpixVector<T> &other) const
{
	assert(N == other.N);
	UpixVector<T> result(N, m_pData);
	result += other;
	return result;
}

template <class T>
UpixVector<T> UpixVector<T>::operator+(const T &value) const
{
	UpixVector<T> result(N, m_pData);
	result += value;
	return result;
}

template <class T>
UpixVector<T> UpixVector<T>::operator-(const UpixVector<T> &other) const
{
	assert(N == other.N);
	UpixVector<T> result(N, m_pData);
	result -= other;
	return result;
}

template <class T>
UpixVector<T> UpixVector<T>::operator-(const T &value) const
{
	UpixVector<T> result(N, m_pData);
	result -= value;
	return result;
}

template <class T>
UpixVector<T> UpixVector<T>::operator*(const T &value) const
{
	UpixVector<T> result(N, m_pData);
	result *= value;
	return result;
}

template <class T>
UpixVector<T> UpixVector<T>::operator/(const T &value) const
{
	assert(value != 0);
	UpixVector<T> result(N, m_pData);
	result /= value;
	return result;
}

template <class T>
UpixVector<T> &UpixVector<T>::operator+=(const UpixVector<T> &other)
{
	assert(N == other.N);
	for (uint i=0; i<N; i++)
	{
		m_pData[i] += other.m_pData[i];
	}
	return *this;
}

template <class T>
UpixVector<T> &UpixVector<T>::operator+=(const T &value)
{
	for (uint i=0; i<N; i++)
	{
		m_pData[i] += value;
	}
	return *this;
}

template <class T>
UpixVector<T> &UpixVector<T>::operator-=(const UpixVector<T> &other)
{
	assert(N == other.N);
	for (uint i=0; i<N; i++)
	{
		m_pData[i] -= other.m_pData[i];
	}
	return *this;
}

template <class T>
UpixVector<T> &UpixVector<T>::operator-=(const T &value)
{
	for (uint i=0; i<N; i++)
	{
		m_pData[i] -= value;
	}
	return *this;
}

template <class T>
UpixVector<T> &UpixVector<T>::operator*=(const T &value)
{
	for (uint i=0; i<N; i++)
	{
		m_pData[i] *= value;
	}
	return *this;
}

template <class T>
UpixVector<T> &UpixVector<T>::operator/=(const T &value)
{
	assert(value != 0);
	T s = 1/value;
	for (uint i=0; i<N; i++)
	{
		m_pData[i] *= s;
	}
	return *this;
}

template <class T>
inline UpixVector<T> UpixVector<T>::operator-()
{
	T value = 0;
	return value - *this;
}


//-----------------------------Utilities--------------------------------
//向量左加数
template <class T>
inline UpixVector<T> operator+(const T &value, const UpixVector<T> &vec)
{
	return vec + value;
}

//向量左减数
template <class T>
UpixVector<T> operator-(const T &value, const UpixVector<T> &vec)
{
	uint n = vec.Size();
	UpixVector<T> result(n);
	for (uint i=0; i<n; i++)
	{
		result[i] = value - vec[i];
	}
	return result;
}

//向量左乘数
template <class T>
inline UpixVector<T> operator*(const T &value, const UpixVector<T> &vec)
{
	return vec*value;
}

//向量左除数
template <class T>
UpixVector<T> operator/(const T &value, const UpixVector<T> &vec)
{
	uint n = vec.Size();
	UpixVector<T> result(n);
	for (uint i=0; i<n; i++)
	{
		if (vec[i]==0)
		{
			assert(false);
			return result;
		}
		result[i] = value / vec[i];
	}
	return result;
}

template <class T>
inline UpixVector<T> Cross(const UpixVector<T> &vec1, const UpixVector<T> &vec2)
{
	assert((vec1.Size() == 3) && (vec2.Size() == 3));
	return vec1.Cross(vec2);
}

template <class T>
inline T Dot(const UpixVector<T> &vec1, const UpixVector<T> &vec2)
{
	assert(vec1.Size() == vec2.Size());
	return vec1.Dot(vec2);
}


#endif