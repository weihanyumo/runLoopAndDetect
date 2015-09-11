/*
* Filename: UpixMatrix4x4.h
* Author:   Xiehong
* Date:     2013.7.13
*/
#ifndef _UPIX_MATRIX_4X4_H_
#define _UPIX_MATRIX_4X4_H_


#include "UpixAlgorithmTypeDef.h"
#include "UpixPoint.h"


template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixMatrix4x4
{
public:
	UpixMatrix4x4()
	{
		SetValue(0);
	}

	UpixMatrix4x4(const T *pData)
	{
		SetData(pData);
	}

	UpixMatrix4x4(const T pData[4][4])
	{
		SetData(pData);
	}

	UpixMatrix4x4(const UpixMatrix4x4 &other)
	{
		SetData(other.m_pData);
	}

	~UpixMatrix4x4()
	{

	}

public:
	void SetValue(const T &val)
	{
		m_pData[0][0] = val; m_pData[0][1] = val; m_pData[0][2] = val; m_pData[0][3] = val;
		m_pData[1][0] = val; m_pData[1][1] = val; m_pData[1][2] = val; m_pData[1][3] = val;
		m_pData[2][0] = val; m_pData[2][1] = val; m_pData[2][2] = val; m_pData[2][3] = val;
		m_pData[3][0] = val; m_pData[3][1] = val; m_pData[3][2] = val; m_pData[3][3] = val;
	}

	void SetData(const T *pData)
	{
		m_pData[0][0] = pData[0]; m_pData[0][1] = pData[1]; m_pData[0][2] = pData[2]; m_pData[0][3] = pData[3];
		m_pData[1][0] = pData[4]; m_pData[1][1] = pData[5]; m_pData[1][2] = pData[6]; m_pData[1][3] = pData[7];
		m_pData[2][0] = pData[8]; m_pData[2][1] = pData[9]; m_pData[2][2] = pData[10]; m_pData[2][3] = pData[11];
		m_pData[3][0] = pData[12]; m_pData[3][1] = pData[13]; m_pData[3][2] = pData[14]; m_pData[3][3] = pData[15];
	}

	void SetData(const T pData[4][4])
	{
		T *p = (T*)pData;
		SetData(p);
	}

	void SetIdentity()
	{
		m_pData[0][0] = 1; m_pData[0][1] = 0; m_pData[0][2] = 0; m_pData[0][3] = 0;
		m_pData[1][0] = 0; m_pData[1][1] = 1; m_pData[1][2] = 0; m_pData[1][3] = 0;
		m_pData[2][0] = 0; m_pData[2][1] = 0; m_pData[2][2] = 1; m_pData[2][3] = 0;
		m_pData[3][0] = 0; m_pData[3][1] = 0; m_pData[3][2] = 0; m_pData[3][3] = 1;
	}

	inline void SetZeros()
	{
		SetValue(0);
	}

	void SetRow(int i, const T v[4])
	{
		assert(i>=0 && i<4);
		m_pData[i][0] = v[0];
		m_pData[i][1] = v[1];
		m_pData[i][2] = v[2];
		m_pData[i][3] = v[3];
	}

	void SetCol(int i, const T v[4])
	{
		assert(i>=0 && i<4);
		m_pData[0][i] = v[0];
		m_pData[1][i] = v[1];
		m_pData[2][i] = v[2];
		m_pData[3][i] = v[3];
	}

	void GetRow(int i, T v[4]) const
	{
		assert(i>=0 && i<4);
		v[0] = m_pData[i][0];
		v[1] = m_pData[i][1];
		v[2] = m_pData[i][2];
		v[3] = m_pData[i][3];
	}

	void GetCol(int i, T v[4]) const
	{
		assert(i>=0 && i<4);
		v[0] = m_pData[0][i];
		v[1] = m_pData[1][i];
		v[2] = m_pData[2][i];
		v[3] = m_pData[3][i];
	}

	T *GetData() const
	{
		return (T *)m_pData;
	}

	UpixMatrix4x4 &operator=(const UpixMatrix4x4 &rhs)
	{
		if (this != &rhs)
		{
			SetData(rhs.m_pData);
		}
		return *this;
	}

	T *operator[](int i)
	{
		return m_pData[i];
	}

	const T *operator[](int i) const
	{
		return m_pData[i];
	}

	T &operator()(int nRow, int nCol)
	{
		return m_pData[nRow][nCol];
	}

	const T &operator()(int nRow, int nCol) const
	{
		return m_pData[nRow][nCol];
	}

	UpixMatrix4x4<T> operator+(const UpixMatrix4x4<T> &other) const
	{
		UpixMatrix4x4<T> result(m_pData);
		result += other;
		return result;
	}

	UpixMatrix4x4<T> operator-(const UpixMatrix4x4<T> &other) const
	{
		UpixMatrix4x4<T> result(m_pData);
		result -= other;
		return result;
	}

	UpixMatrix4x4<T> operator*(const UpixMatrix4x4<T> &other) const
	{
		UpixMatrix4x4<T> result;
		T *pMr = (T *)result.m_pData;
		T *pMo = (T *)other.m_pData;

		pMr[0] = m_pData[0][0]*pMo[0] + m_pData[0][1]*pMo[4] + m_pData[0][2]*pMo[8] + m_pData[0][3]*pMo[12];
		pMr[1] = m_pData[0][0]*pMo[1] + m_pData[0][1]*pMo[5] + m_pData[0][2]*pMo[9] + m_pData[0][3]*pMo[13];
		pMr[2] = m_pData[0][0]*pMo[2] + m_pData[0][1]*pMo[6] + m_pData[0][2]*pMo[10] + m_pData[0][3]*pMo[14];
		pMr[3] = m_pData[0][0]*pMo[3] + m_pData[0][1]*pMo[7] + m_pData[0][2]*pMo[11] + m_pData[0][3]*pMo[15];

		pMr[4] = m_pData[1][0]*pMo[0] + m_pData[1][1]*pMo[4] + m_pData[1][2]*pMo[8] + m_pData[1][3]*pMo[12];
		pMr[5] = m_pData[1][0]*pMo[1] + m_pData[1][1]*pMo[5] + m_pData[1][2]*pMo[9] + m_pData[1][3]*pMo[13];
		pMr[6] = m_pData[1][0]*pMo[2] + m_pData[1][1]*pMo[6] + m_pData[1][2]*pMo[10] + m_pData[1][3]*pMo[14];
		pMr[7] = m_pData[1][0]*pMo[3] + m_pData[1][1]*pMo[7] + m_pData[1][2]*pMo[11] + m_pData[1][3]*pMo[15];

		pMr[8] = m_pData[2][0]*pMo[0] + m_pData[2][1]*pMo[4] + m_pData[2][2]*pMo[8] + m_pData[2][3]*pMo[12];
		pMr[9] = m_pData[2][0]*pMo[1] + m_pData[2][1]*pMo[5] + m_pData[2][2]*pMo[9] + m_pData[2][3]*pMo[13];
		pMr[10] = m_pData[2][0]*pMo[2] + m_pData[2][1]*pMo[6] + m_pData[2][2]*pMo[10] + m_pData[2][3]*pMo[14];
		pMr[11] = m_pData[2][0]*pMo[3] + m_pData[2][1]*pMo[7] + m_pData[2][2]*pMo[11] + m_pData[2][3]*pMo[15];

		pMr[12] = m_pData[3][0]*pMo[0] + m_pData[3][1]*pMo[4] + m_pData[3][2]*pMo[8] + m_pData[3][3]*pMo[12];
		pMr[13] = m_pData[3][0]*pMo[1] + m_pData[3][1]*pMo[5] + m_pData[3][2]*pMo[9] + m_pData[3][3]*pMo[13];
		pMr[14] = m_pData[3][0]*pMo[2] + m_pData[3][1]*pMo[6] + m_pData[3][2]*pMo[10] + m_pData[3][3]*pMo[14];
		pMr[15] = m_pData[3][0]*pMo[3] + m_pData[3][1]*pMo[7] + m_pData[3][2]*pMo[11] + m_pData[3][3]*pMo[15];

		return result;
	}

	UpixMatrix4x4<T> operator+(const T &value) const
	{
		UpixMatrix4x4<T> result(m_pData);
		result += value;
		return result;
	}

	UpixMatrix4x4<T> operator-(const T &value) const
	{
		UpixMatrix4x4<T> result(m_pData);
		result -= value;
		return result;
	}

	UpixMatrix4x4<T> operator*(const T &value) const
	{
		UpixMatrix4x4<T> result(m_pData);
		result *= value;
		return result;
	}

	UpixMatrix4x4<T> operator/(const T &value) const
	{
		UpixMatrix4x4<T> result(m_pData);
		result /= value;
		return result;
	}

	UpixPoint3D<T> operator*(const UpixPoint3D<T> &v) const
	{
		UpixPoint3D<T> result;
		T invw = 1/(m_pData[3][0]*v.x + m_pData[3][1]*v.y + m_pData[3][2]*v.z + m_pData[3][3]);
		result.x = (m_pData[0][0]*v.x + m_pData[0][1]*v.y + m_pData[0][2]*v.z + m_pData[0][3])*invw;
		result.y = (m_pData[1][0]*v.x + m_pData[1][1]*v.y + m_pData[1][2]*v.z + m_pData[1][3])*invw;
		result.z = (m_pData[2][0]*v.x + m_pData[2][1]*v.y + m_pData[2][2]*v.z + m_pData[2][3])*invw;
		return result;
	}

	UpixPointHomg3D<T> operator*(const UpixPointHomg3D<T> &v) const
	{
		UpixPointHomg3D<T> result;
		result.x = m_pData[0][0]*v.x + m_pData[0][1]*v.y + m_pData[0][2]*v.z + m_pData[0][3]*v.w;
		result.y = m_pData[1][0]*v.x + m_pData[1][1]*v.y + m_pData[1][2]*v.z + m_pData[1][3]*v.w;
		result.z = m_pData[2][0]*v.x + m_pData[2][1]*v.y + m_pData[2][2]*v.z + m_pData[2][3]*v.w;
		result.w = m_pData[3][0]*v.x + m_pData[3][1]*v.y + m_pData[3][2]*v.z + m_pData[3][3]*v.w;
		return result;
	}

	UpixMatrix4x4<T>& operator+=(const UpixMatrix4x4<T> &other)
	{
		T *pMo = (T *)other.m_pData;

		m_pData[0][0] += pMo[0]; m_pData[0][1] += pMo[1]; m_pData[0][2] += pMo[2]; m_pData[0][3] += pMo[3];
		m_pData[1][0] += pMo[4]; m_pData[1][1] += pMo[5]; m_pData[1][2] += pMo[6]; m_pData[1][3] += pMo[7];
		m_pData[2][0] += pMo[8]; m_pData[2][1] += pMo[9]; m_pData[2][2] += pMo[10]; m_pData[2][3] += pMo[11];
		m_pData[3][0] += pMo[12]; m_pData[3][1] += pMo[13]; m_pData[3][2] += pMo[14]; m_pData[3][3] += pMo[15];
		return *this;
	}

	UpixMatrix4x4<T>& operator-=(const UpixMatrix4x4<T> &other)
	{
		T *pMo = (T *)other.m_pData;

		m_pData[0][0] -= pMo[0]; m_pData[0][1] -= pMo[1]; m_pData[0][2] -= pMo[2]; m_pData[0][3] -= pMo[3];
		m_pData[1][0] -= pMo[4]; m_pData[1][1] -= pMo[5]; m_pData[1][2] -= pMo[6]; m_pData[1][3] -= pMo[7];
		m_pData[2][0] -= pMo[8]; m_pData[2][1] -= pMo[9]; m_pData[2][2] -= pMo[10]; m_pData[2][3] -= pMo[11];
		m_pData[3][0] -= pMo[12]; m_pData[3][1] -= pMo[13]; m_pData[3][2] -= pMo[14]; m_pData[3][3] -= pMo[15];
		return *this;
	}

	UpixMatrix4x4<T>& operator*=(const UpixMatrix4x4<T> &other)
	{
		UpixMatrix4x4<T> mat(m_pData);
		T *pMl = (T *)mat.m_pData;
		T *pMo = (T *)other.m_pData;
	
		m_pData[0][0] = pMl[0]*pMo[0] + pMl[1]*pMo[4] + pMl[2]*pMo[8] + pMl[3]*pMo[12];
		m_pData[0][1] = pMl[0]*pMo[1] + pMl[1]*pMo[5] + pMl[2]*pMo[9] + pMl[3]*pMo[13];
		m_pData[0][2] = pMl[0]*pMo[2] + pMl[1]*pMo[6] + pMl[2]*pMo[10] + pMl[3]*pMo[14];
		m_pData[0][3] = pMl[0]*pMo[3] + pMl[1]*pMo[7] + pMl[2]*pMo[11] + pMl[3]*pMo[15];

		m_pData[1][0] = pMl[4]*pMo[0] + pMl[5]*pMo[4] + pMl[6]*pMo[8] + pMl[7]*pMo[12];
		m_pData[1][1] = pMl[4]*pMo[1] + pMl[5]*pMo[5] + pMl[6]*pMo[9] + pMl[7]*pMo[13];
		m_pData[1][2] = pMl[4]*pMo[2] + pMl[5]*pMo[6] + pMl[6]*pMo[10] + pMl[7]*pMo[14];
		m_pData[1][3] = pMl[4]*pMo[3] + pMl[5]*pMo[7] + pMl[6]*pMo[11] + pMl[7]*pMo[15];

		m_pData[2][0] = pMl[8]*pMo[0] + pMl[9]*pMo[4] + pMl[10]*pMo[8] + pMl[11]*pMo[12];
		m_pData[2][1] = pMl[8]*pMo[1] + pMl[9]*pMo[5] + pMl[10]*pMo[9] + pMl[11]*pMo[13];
		m_pData[2][2] = pMl[8]*pMo[2] + pMl[9]*pMo[6] + pMl[10]*pMo[10] + pMl[11]*pMo[14];
		m_pData[2][3] = pMl[8]*pMo[3] + pMl[9]*pMo[7] + pMl[10]*pMo[11] + pMl[11]*pMo[15];

		m_pData[3][0] = pMl[12]*pMo[0] + pMl[13]*pMo[4] + pMl[14]*pMo[8] + pMl[15]*pMo[12];
		m_pData[3][1] = pMl[12]*pMo[1] + pMl[13]*pMo[5] + pMl[14]*pMo[9] + pMl[15]*pMo[13];
		m_pData[3][2] = pMl[12]*pMo[2] + pMl[13]*pMo[6] + pMl[14]*pMo[10] + pMl[15]*pMo[14];
		m_pData[3][3] = pMl[12]*pMo[3] + pMl[13]*pMo[7] + pMl[14]*pMo[11] + pMl[15]*pMo[15];

		return *this;
	}

	UpixMatrix4x4<T>& operator+=(const T &value)
	{
		m_pData[0][0] += value; m_pData[0][1] += value; m_pData[0][2] += value; m_pData[0][3] += value;
		m_pData[1][0] += value; m_pData[1][1] += value; m_pData[1][2] += value; m_pData[1][3] += value;
		m_pData[2][0] += value; m_pData[2][1] += value; m_pData[2][2] += value; m_pData[2][3] += value;
		m_pData[3][0] += value; m_pData[3][1] += value; m_pData[3][2] += value; m_pData[3][3] += value;

		return *this;
	}

	UpixMatrix4x4<T>& operator-=(const T &value)
	{
		m_pData[0][0] -= value; m_pData[0][1] -= value; m_pData[0][2] -= value; m_pData[0][3] -= value;
		m_pData[1][0] -= value; m_pData[1][1] -= value; m_pData[1][2] -= value; m_pData[1][3] -= value;
		m_pData[2][0] -= value; m_pData[2][1] -= value; m_pData[2][2] -= value; m_pData[2][3] -= value;
		m_pData[3][0] -= value; m_pData[3][1] -= value; m_pData[3][2] -= value; m_pData[3][3] -= value;

		return *this;
	}

	UpixMatrix4x4<T>& operator*=(const T &value)
	{
		m_pData[0][0] *= value; m_pData[0][1] *= value; m_pData[0][2] *= value; m_pData[0][3] *= value;
		m_pData[1][0] *= value; m_pData[1][1] *= value; m_pData[1][2] *= value; m_pData[1][3] *= value;
		m_pData[2][0] *= value; m_pData[2][1] *= value; m_pData[2][2] *= value; m_pData[2][3] *= value;
		m_pData[3][0] *= value; m_pData[3][1] *= value; m_pData[3][2] *= value; m_pData[3][3] *= value;

		return *this;
	}

	UpixMatrix4x4<T>& operator/=(const T &value)
	{
		T val = 1 / value;
		m_pData[0][0] *= val; m_pData[0][1] *= val; m_pData[0][2] *= val;
		m_pData[1][0] *= val; m_pData[1][1] *= val; m_pData[1][2] *= val;
		m_pData[2][0] *= val; m_pData[2][1] *= val; m_pData[2][2] *= val;

		return *this;
	}

	//change sign
	UpixMatrix4x4<T> operator-() const
	{
		UpixMatrix4x4<T> result;
		T *p = (T*)result.m_pData;

		p[0] = -m_pData[0][0]; p[1] = -m_pData[0][1]; p[2] = -m_pData[0][2]; p[3] = -m_pData[0][3];
		p[4] = -m_pData[1][0]; p[5] = -m_pData[1][1]; p[6] = -m_pData[1][2]; p[7] = -m_pData[1][3];
		p[8] = -m_pData[2][0]; p[9] = -m_pData[2][1]; p[10] = -m_pData[2][2]; p[11] = -m_pData[2][3];
		p[12] = -m_pData[3][0]; p[13] = -m_pData[3][1]; p[14] = -m_pData[3][2]; p[15] = -m_pData[3][3];

		return result;
	}

	void NormalizeRow()
	{
		T one(1);
		T invd = one/sqrt(m_pData[0][0]*m_pData[0][0] + m_pData[0][1]*m_pData[0][1] + m_pData[0][2]*m_pData[0][2] + m_pData[0][3]*m_pData[0][3]);
		m_pData[0][0] *= invd;
		m_pData[0][1] *= invd;
		m_pData[0][2] *= invd;
		m_pData[0][3] *= invd;

		invd = one/sqrt(m_pData[1][0]*m_pData[1][0] + m_pData[1][1]*m_pData[1][1] + m_pData[1][2]*m_pData[1][2] + m_pData[1][3]*m_pData[1][3]);
		m_pData[1][0] *= invd;
		m_pData[1][1] *= invd;
		m_pData[1][2] *= invd;
		m_pData[1][3] *= invd;

		invd = one/sqrt(m_pData[2][0]*m_pData[2][0] + m_pData[2][1]*m_pData[2][1] + m_pData[2][2]*m_pData[2][2] + m_pData[2][3]*m_pData[2][3]);
		m_pData[2][0] *= invd;
		m_pData[2][1] *= invd;
		m_pData[2][2] *= invd;
		m_pData[2][3] *= invd;

		invd = one/sqrt(m_pData[3][0]*m_pData[3][0] + m_pData[3][1]*m_pData[3][1] + m_pData[3][2]*m_pData[3][2] + m_pData[3][3]*m_pData[3][3]);
		m_pData[3][0] *= invd;
		m_pData[3][1] *= invd;
		m_pData[3][2] *= invd;
		m_pData[3][3] *= invd;
	}

	void NormalizeCol()
	{
		T one(1);
		T invd = one/sqrt(m_pData[0][0]*m_pData[0][0] + m_pData[1][0]*m_pData[1][0] + m_pData[2][0]*m_pData[2][0] + m_pData[3][0]*m_pData[3][0]);
		m_pData[0][0] *= invd;
		m_pData[1][0] *= invd;
		m_pData[2][0] *= invd;
		m_pData[3][0] *= invd;

		invd = one/sqrt(m_pData[0][1]*m_pData[0][1] + m_pData[1][1]*m_pData[1][1] + m_pData[2][1]*m_pData[2][1] + m_pData[3][1]*m_pData[3][1]);
		m_pData[0][1] *= invd;
		m_pData[1][1] *= invd;
		m_pData[2][1] *= invd;
		m_pData[3][1] *= invd;

		invd = one/sqrt(m_pData[0][2]*m_pData[0][2] + m_pData[1][2]*m_pData[1][2] + m_pData[2][2]*m_pData[2][2] + m_pData[3][2]*m_pData[3][2]);
		m_pData[0][2] *= invd;
		m_pData[1][2] *= invd;
		m_pData[2][2] *= invd;
		m_pData[3][2] *= invd;

		invd = one/sqrt(m_pData[0][3]*m_pData[0][3] + m_pData[1][3]*m_pData[1][3] + m_pData[2][3]*m_pData[2][3] + m_pData[3][3]*m_pData[3][3]);
		m_pData[0][3] *= invd;
		m_pData[1][3] *= invd;
		m_pData[2][3] *= invd;
		m_pData[3][3] *= invd;
	}

	T Trace() const
	{
		return (m_pData[0][0] + m_pData[1][1] + m_pData[2][2] + m_pData[3][3]);
	}

	// Determinant 
	T Det() const
	{
		return 
			(
			(m_pData[0][0] * m_pData[1][1] - m_pData[0][1] * m_pData[1][0]) * (m_pData[2][2] * m_pData[3][3] - m_pData[2][3] * m_pData[3][2]) +
			(m_pData[0][0] * m_pData[1][3] - m_pData[0][3] * m_pData[1][0]) * (m_pData[2][1] * m_pData[3][2] - m_pData[2][2] * m_pData[3][1]) + 
			(m_pData[0][1] * m_pData[1][2] - m_pData[0][2] * m_pData[1][1]) * (m_pData[2][0] * m_pData[3][3] - m_pData[2][3] * m_pData[3][0]) +
			(m_pData[0][2] * m_pData[1][0] - m_pData[0][0] * m_pData[1][2]) * (m_pData[2][1] * m_pData[3][3] - m_pData[2][3] * m_pData[3][1]) +
			(m_pData[0][2] * m_pData[1][3] - m_pData[0][3] * m_pData[1][2]) * (m_pData[2][0] * m_pData[3][1] - m_pData[2][1] * m_pData[3][0]) +
			(m_pData[0][3] * m_pData[1][1] - m_pData[0][1] * m_pData[1][3]) * (m_pData[2][0] * m_pData[3][2] - m_pData[2][2] * m_pData[3][0])
			);
	}

	UpixMatrix4x4<T> Transpose() const
	{
		UpixMatrix4x4<T> result;
		Transpose(result);
		return result;
	}

	inline void Transpose(UpixMatrix4x4<T> &ret) const
	{
		T *pMr = (T *)ret.m_pData;
		pMr[0] = m_pData[0][0]; pMr[1] = m_pData[1][0]; pMr[2] = m_pData[2][0]; pMr[3] = m_pData[3][0];
		pMr[4] = m_pData[0][1]; pMr[5] = m_pData[1][1]; pMr[6] = m_pData[2][1]; pMr[7] = m_pData[3][1];
		pMr[8] = m_pData[0][2]; pMr[9] = m_pData[1][2]; pMr[10] = m_pData[2][2]; pMr[11] = m_pData[3][2];
		pMr[12] = m_pData[0][3]; pMr[13] = m_pData[1][3]; pMr[14] = m_pData[2][3]; pMr[15] = m_pData[3][3];
	}

	UpixMatrix4x4<T> Inv() const
	{
		UpixMatrix4x4<T> ret;
		Inv(ret);
		return ret;
	}

	T Inv(UpixMatrix4x4<T> &m) const
	{
		T data1 = m_pData[2][2] * m_pData[3][3] - m_pData[2][3] * m_pData[3][2];
		T data2 = m_pData[2][1] * m_pData[3][3] - m_pData[2][3] * m_pData[3][1];
		T data3 = m_pData[2][1] * m_pData[3][2] - m_pData[2][2] * m_pData[3][1];
		T data4 = m_pData[2][3] * m_pData[3][0] - m_pData[2][0] * m_pData[3][3];
		T data5 = m_pData[2][0] * m_pData[3][2] - m_pData[2][2] * m_pData[3][0];
		T data6 = m_pData[2][0] * m_pData[3][1] - m_pData[2][1] * m_pData[3][0];

		T det = 
			(m_pData[0][0] * m_pData[1][1] - m_pData[0][1] * m_pData[1][0]) * data1 +
			(m_pData[0][0] * m_pData[1][3] - m_pData[0][3] * m_pData[1][0]) * data3 - 
			(m_pData[0][1] * m_pData[1][2] - m_pData[0][2] * m_pData[1][1]) * data4 +
			(m_pData[0][2] * m_pData[1][0] - m_pData[0][0] * m_pData[1][2]) * data2 +
			(m_pData[0][2] * m_pData[1][3] - m_pData[0][3] * m_pData[1][2]) * data6 +
			(m_pData[0][3] * m_pData[1][1] - m_pData[0][1] * m_pData[1][3]) * data5;

		if (UpixAbs(det) < UPIX_FLOAT_EPSILON)
		{
			return det;
		}

		T one(1);
		T dInvDet = one / det;

		T data7 = m_pData[1][2] * m_pData[3][3] - m_pData[1][3] * m_pData[3][2];
		T data8 = m_pData[1][3] * m_pData[3][1] - m_pData[1][1] * m_pData[3][3];
		T data9 = m_pData[1][1] * m_pData[3][2] - m_pData[1][2] * m_pData[3][1];
		T data10= m_pData[1][3] * m_pData[3][0] - m_pData[1][0] * m_pData[3][3];
		T data11= m_pData[1][0] * m_pData[3][2] - m_pData[1][2] * m_pData[3][0];
		T data12= m_pData[1][0] * m_pData[3][1] - m_pData[1][1] * m_pData[3][0];
		T data13= m_pData[1][2] * m_pData[2][3] - m_pData[1][3] * m_pData[2][2];
		T data14= m_pData[1][3] * m_pData[2][1] - m_pData[1][1] * m_pData[2][3];
		T data15= m_pData[1][1] * m_pData[2][2] - m_pData[1][2] * m_pData[2][1];
		T data16= m_pData[1][3] * m_pData[2][0] - m_pData[1][0] * m_pData[2][3];
		T data17= m_pData[1][0] * m_pData[2][2] - m_pData[1][2] * m_pData[2][0];
		T data18= m_pData[1][0] * m_pData[2][1] - m_pData[1][1] * m_pData[2][0];

		T * pm = m.GetData();

		pm[0] = (m_pData[1][1] * data1  - m_pData[1][2] * data2  + m_pData[1][3] * data3) * dInvDet;
		pm[1] =-(m_pData[0][1] * data1  - m_pData[0][2] * data2  + m_pData[0][3] * data3) * dInvDet;
		pm[2] = (m_pData[0][1] * data7  + m_pData[0][2] * data8  + m_pData[0][3] * data9) * dInvDet;
		pm[3] =-(m_pData[0][1] * data13 + m_pData[0][2] * data14 + m_pData[0][3] * data15) * dInvDet;
		pm[4] =-(m_pData[1][0] * data1  + m_pData[1][2] * data4  + m_pData[1][3] * data5) * dInvDet;
		pm[5] = (m_pData[0][0] * data1  + m_pData[0][2] * data4  + m_pData[0][3] * data5) * dInvDet;
		pm[6] =-(m_pData[0][0] * data7  + m_pData[0][2] * data10 + m_pData[0][3] * data11) * dInvDet;
		pm[7] = (m_pData[0][0] * data13 + m_pData[0][2] * data16 + m_pData[0][3] * data17) * dInvDet;
		pm[8] = (m_pData[1][0] * data2  + m_pData[1][1] * data4  + m_pData[1][3] * data6) * dInvDet;
		pm[9] =-(m_pData[0][0] * data2  + m_pData[0][1] * data4  + m_pData[0][3] * data6) * dInvDet;
		pm[10]= (m_pData[0][1] * data10 - m_pData[0][0] * data8  + m_pData[0][3] * data12) * dInvDet;
		pm[11]=-(m_pData[0][1] * data16 - m_pData[0][0] * data14 + m_pData[0][3] * data18) * dInvDet;
		pm[12]=-(m_pData[1][0] * data3  - m_pData[1][1] * data5  + m_pData[1][2] * data6) * dInvDet;
		pm[13]= (m_pData[0][0] * data3  - m_pData[0][1] * data5  + m_pData[0][2] * data6) * dInvDet;
		pm[14]=-(m_pData[0][0] * data9  - m_pData[0][1] * data11 + m_pData[0][2] * data12) * dInvDet;
		pm[15]= (m_pData[0][0] * data15 - m_pData[0][1] * data17 + m_pData[0][2] * data18) * dInvDet;

		return det;
	}

	// Get normal infinite, computes the maximum of the sums of absolute values over rows
	T InfiniteNorm() const
	{
		T sum1	= abs(m_pData[0][0]) + abs(m_pData[0][1]) + abs(m_pData[0][2]) + abs(m_pData[0][3]);
		T sum2	= abs(m_pData[1][0]) + abs(m_pData[1][1]) + abs(m_pData[1][2]) + abs(m_pData[1][3]);
		T sum3	= abs(m_pData[2][0]) + abs(m_pData[2][1]) + abs(m_pData[2][2]) + abs(m_pData[2][3]);
		T sum4	= abs(m_pData[3][0]) + abs(m_pData[3][1]) + abs(m_pData[3][2]) + abs(m_pData[3][3]);
		return UpixMax(UpixMax(UpixMax(sum1, sum2), sum3), sum4);
	}

	// Computes the matrix exponential of a matrix m by scaling m by 1/(powers of 2), 
	// using Taylor series and squaring again. The input matrix must be square
	UpixMatrix4x4<T> Exp() const
	{
		UpixMatrix4x4<T> result;
		Exp(result);
		return result;
	}

	void Exp(UpixMatrix4x4<T> &result) const
	{
		T level = log10(InfiniteNorm())*UPIX_INV_LOG10_2;	//const P l = log2(norm_inf(m));
		int scale = UpixMax(0, (int)ceil(level));

		exp_taylor(*this/(1<<scale), result);
		for(int i=0; i<scale; ++i)
		{
			result *= result;
		}
	}

private:
	// Exponentiate a matrix using a the Taylor series
	// This will not work if the norm of the matrix is too large.
	static void exp_taylor(const UpixMatrix4x4<T> &m, UpixMatrix4x4<T> &result)
	{
		result.SetValue(0);
		UpixMatrix4x4<T> f;
		T *pMf = f.GetData();
		pMf[0] = 1;
		pMf[5] = 1;
		pMf[10] = 1;
		pMf[15] = 1;

		T k = 1;
		while (f.InfiniteNorm() > 0)
		{
			result += f;
			f *= m;
			f /= k;
			k += 1;
		}
	}

private:
	T m_pData[4][4];
};



typedef UpixMatrix4x4<float> UpixMatrix4x4f;
typedef UpixMatrix4x4<double> UpixMatrix4x4d;




#endif
