/*
* Filename: UpixMatrix3x3.h
* Author:   Xiehong
* Date:     2013.7.7
*/

#ifndef _UPIX_MATRIX_3X3_H_
#define _UPIX_MATRIX_3X3_H_


#include "UpixAlgorithmTypeDef.h"
#include "UpixPoint.h"


template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixMatrix3x3
{
public:
	UpixMatrix3x3()
	{
		SetValue(0);
	}

	UpixMatrix3x3(const T *pData)
	{
		SetData(pData);
	}

	UpixMatrix3x3(const T pData[3][3])
	{
		SetData(pData);
	}

	UpixMatrix3x3(const UpixMatrix3x3 &other)
	{
		SetData(other.m_pData);
	}

	~UpixMatrix3x3()
	{

	}

public:
	void SetValue(const T &val)
	{
		m_pData[0][0] = val; m_pData[0][1] = val; m_pData[0][2] = val;
		m_pData[1][0] = val; m_pData[1][1] = val; m_pData[1][2] = val;
		m_pData[2][0] = val; m_pData[2][1] = val; m_pData[2][2] = val;
	}

	void SetData(const T *pData)
	{
		m_pData[0][0] = pData[0]; m_pData[0][1] = pData[1]; m_pData[0][2] = pData[2];
		m_pData[1][0] = pData[3]; m_pData[1][1] = pData[4]; m_pData[1][2] = pData[5];
		m_pData[2][0] = pData[6]; m_pData[2][1] = pData[7]; m_pData[2][2] = pData[8];
	}

	void SetData(const T pData[3][3])
	{
		T *p = (T*)pData;
		SetData(p);
	}

	void SetIdentity()
	{
		m_pData[0][0] = 1; m_pData[0][1] = 0; m_pData[0][2] = 0;
		m_pData[1][0] = 0; m_pData[1][1] = 1; m_pData[1][2] = 0;
		m_pData[2][0] = 0; m_pData[2][1] = 0; m_pData[2][2] = 1;
	}

	inline void SetZeros()
	{
		SetValue(0);
	}

	void SetRow(int i, const UpixPoint3D<T> &v)
	{
		assert(i>=0 && i<3);
		m_pData[i][0] = v.x;
		m_pData[i][1] = v.y;
		m_pData[i][2] = v.z;
	}

	void SetCol(int i, const UpixPoint3D<T> &v)
	{
		assert(i>=0 && i<3);
		m_pData[0][i] = v.x;
		m_pData[1][i] = v.y;
		m_pData[2][i] = v.z;
	}

	UpixPoint3D<T> GetRow(int i) const
	{
		UpixPoint3D<T> v;
		GetRow(i, v);
		return v;
	}

	void GetRow(int i, UpixPoint3D<T> &v) const
	{
		assert(i>=0 && i<3);
		v.x = m_pData[i][0];
		v.y = m_pData[i][1];
		v.z = m_pData[i][2];
	}

	UpixPoint3D<T> GetCol(int i) const
	{
		UpixPoint3D<T> v;
		GetCol(i, v);
		return v;
	}

	void GetCol(int i, UpixPoint3D<T> &v) const
	{
		assert(i>=0 && i<3);
		v.x = m_pData[0][i];
		v.y = m_pData[1][i];
		v.z = m_pData[2][i];
	}

	T *GetData() const
	{
		return (T *)m_pData;
	}

	UpixMatrix3x3 &operator=(const UpixMatrix3x3 &rhs)
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

	UpixMatrix3x3<T> operator+(const UpixMatrix3x3<T> &other) const
	{
		UpixMatrix3x3<T> result(m_pData);
		result += other;
		return result;
	}

	UpixMatrix3x3<T> operator-(const UpixMatrix3x3<T> &other) const
	{
		UpixMatrix3x3<T> result(m_pData);
		result -= other;
		return result;
	}

	UpixMatrix3x3<T> operator*(const UpixMatrix3x3<T> &other) const
	{
		UpixMatrix3x3<T> result;
		T *pMr = (T *)result.m_pData;
		T *pMo = (T *)other.m_pData;

		pMr[0] = m_pData[0][0] * pMo[0] + m_pData[0][1] * pMo[3] + m_pData[0][2] * pMo[6];
		pMr[1] = m_pData[0][0] * pMo[1] + m_pData[0][1] * pMo[4] + m_pData[0][2] * pMo[7];
		pMr[2] = m_pData[0][0] * pMo[2] + m_pData[0][1] * pMo[5] + m_pData[0][2] * pMo[8];

		pMr[3] = m_pData[1][0] * pMo[0] + m_pData[1][1] * pMo[3] + m_pData[1][2] * pMo[6];
		pMr[4] = m_pData[1][0] * pMo[1] + m_pData[1][1] * pMo[4] + m_pData[1][2] * pMo[7];
		pMr[5] = m_pData[1][0] * pMo[2] + m_pData[1][1] * pMo[5] + m_pData[1][2] * pMo[8];

		pMr[6] = m_pData[2][0] * pMo[0] + m_pData[2][1] * pMo[3] + m_pData[2][2] * pMo[6];
		pMr[7] = m_pData[2][0] * pMo[1] + m_pData[2][1] * pMo[4] + m_pData[2][2] * pMo[7];
		pMr[8] = m_pData[2][0] * pMo[2] + m_pData[2][1] * pMo[5] + m_pData[2][2] * pMo[8];

		return result;
	}

	UpixMatrix3x3<T> operator+(const T &value) const
	{
		UpixMatrix3x3<T> result(m_pData);
		result += value;
		return result;
	}

	UpixMatrix3x3<T> operator-(const T &value) const
	{
		UpixMatrix3x3<T> result(m_pData);
		result -= value;
		return result;
	}

	UpixMatrix3x3<T> operator*(const T &value) const
	{
		UpixMatrix3x3<T> result(m_pData);
		result *= value;
		return result;
	}

	UpixMatrix3x3<T> operator/(const T &value) const
	{
		UpixMatrix3x3<T> result(m_pData);
		result /= value;
		return result;
	}

	UpixPoint2D<T> operator*(const UpixPoint2D<T> &v) const
	{
		UpixPoint2D<T> result;
		T invw = 1/(m_pData[2][0] * v.x + m_pData[2][1] * v.y + m_pData[2][2]);
		result.x = (m_pData[0][0] * v.x + m_pData[0][1] * v.y + m_pData[0][2])*invw;
		result.y = (m_pData[1][0] * v.x + m_pData[1][1] * v.y + m_pData[1][2])*invw;
		return result;
	}

	UpixPoint3D<T> operator*(const UpixPointHomg2D<T> &v) const
	{
		UpixPointHomg2D<T> result;
		result.x = m_pData[0][0] * v.x + m_pData[0][1] * v.y + m_pData[0][2] * v.w;
		result.y = m_pData[1][0] * v.x + m_pData[1][1] * v.y + m_pData[1][2] * v.w;
		result.w = m_pData[2][0] * v.x + m_pData[2][1] * v.y + m_pData[2][2] * v.w;
		return result;
	}

	UpixPoint3D<T> operator*(const UpixPoint3D<T> &v) const
	{
		UpixPoint3D<T> result;
		result.x = m_pData[0][0] * v.x + m_pData[0][1] * v.y + m_pData[0][2] * v.z;
		result.y = m_pData[1][0] * v.x + m_pData[1][1] * v.y + m_pData[1][2] * v.z;
		result.z = m_pData[2][0] * v.x + m_pData[2][1] * v.y + m_pData[2][2] * v.z;
		return result;
	}

	UpixMatrix3x3<T>& operator+=(const UpixMatrix3x3<T> &other)
	{
		T *pMo = (T *)other.m_pData;
		
		m_pData[0][0] += pMo[0]; m_pData[0][1] += pMo[1]; m_pData[0][2] += pMo[2];
		m_pData[1][0] += pMo[3]; m_pData[1][1] += pMo[4]; m_pData[1][2] += pMo[5];
		m_pData[2][0] += pMo[6]; m_pData[2][1] += pMo[7]; m_pData[2][2] += pMo[8];

		return *this;
	}

	UpixMatrix3x3<T>& operator-=(const UpixMatrix3x3<T> &other)
	{
		T *pMo = (T *)other.m_pData;

		m_pData[0][0] -= pMo[0]; m_pData[0][1] -= pMo[1]; m_pData[0][2] -= pMo[2];
		m_pData[1][0] -= pMo[3]; m_pData[1][1] -= pMo[4]; m_pData[1][2] -= pMo[5];
		m_pData[2][0] -= pMo[6]; m_pData[2][1] -= pMo[7]; m_pData[2][2] -= pMo[8];

		return *this;
	}

	UpixMatrix3x3<T>& operator*=(const UpixMatrix3x3<T> &other)
	{
		UpixMatrix3x3<T> mat(m_pData);
		T *pMl = (T *)mat.m_pData;
		T *pMo = (T *)other.m_pData;

		m_pData[0][0] = pMl[0] * pMo[0] + pMl[1] * pMo[3] + pMl[2] * pMo[6];
		m_pData[0][1] = pMl[0] * pMo[1] + pMl[1] * pMo[4] + pMl[2] * pMo[7];
		m_pData[0][2] = pMl[0] * pMo[2] + pMl[1] * pMo[5] + pMl[2] * pMo[8];

		m_pData[1][0] = pMl[3] * pMo[0] + pMl[4] * pMo[3] + pMl[5] * pMo[6];
		m_pData[1][1] = pMl[3] * pMo[1] + pMl[4] * pMo[4] + pMl[5] * pMo[7];
		m_pData[1][2] = pMl[3] * pMo[2] + pMl[4] * pMo[5] + pMl[5] * pMo[8];

		m_pData[2][0] = pMl[6] * pMo[0] + pMl[7] * pMo[3] + pMl[8] * pMo[6];
		m_pData[2][1] = pMl[6] * pMo[1] + pMl[7] * pMo[4] + pMl[8] * pMo[7];
		m_pData[2][2] = pMl[6] * pMo[2] + pMl[7] * pMo[5] + pMl[8] * pMo[8];

		return *this;
	}

	UpixMatrix3x3<T>& operator+=(const T &value)
	{
		m_pData[0][0] += value; m_pData[0][1] += value; m_pData[0][2] += value;
		m_pData[1][0] += value; m_pData[1][1] += value; m_pData[1][2] += value;
		m_pData[2][0] += value; m_pData[2][1] += value; m_pData[2][2] += value;
		return *this;
	}

	UpixMatrix3x3<T>& operator-=(const T &value)
	{
		m_pData[0][0] -= value; m_pData[0][1] -= value; m_pData[0][2] -= value;
		m_pData[1][0] -= value; m_pData[1][1] -= value; m_pData[1][2] -= value;
		m_pData[2][0] -= value; m_pData[2][1] -= value; m_pData[2][2] -= value;
		return *this;
	}

	UpixMatrix3x3<T>& operator*=(const T &value)
	{
		m_pData[0][0] *= value; m_pData[0][1] *= value; m_pData[0][2] *= value;
		m_pData[1][0] *= value; m_pData[1][1] *= value; m_pData[1][2] *= value;
		m_pData[2][0] *= value; m_pData[2][1] *= value; m_pData[2][2] *= value;
		return *this;
	}

	UpixMatrix3x3<T>& operator/=(const T &value)
	{
		T val = 1 / value;
		m_pData[0][0] *= val; m_pData[0][1] *= val; m_pData[0][2] *= val;
		m_pData[1][0] *= val; m_pData[1][1] *= val; m_pData[1][2] *= val;
		m_pData[2][0] *= val; m_pData[2][1] *= val; m_pData[2][2] *= val;
		return *this;
	}

	//change sign
	UpixMatrix3x3<T> operator-() const
	{
		UpixMatrix3x3<T> result;
		T *p = (T*)result.m_pData;

		p[0] = -m_pData[0][0]; p[1] = -m_pData[0][1]; p[2] = -m_pData[0][2];
		p[3] = -m_pData[1][0]; p[4] = -m_pData[1][1]; p[5] = -m_pData[1][2];
		p[6] = -m_pData[2][0]; p[7] = -m_pData[2][1]; p[8] = -m_pData[2][2];

		return result;
	}

	void NormalizeRow()
	{
		T one(1);
		T invd = one/sqrt(m_pData[0][0]*m_pData[0][0] + m_pData[0][1]*m_pData[0][1] + m_pData[0][2]*m_pData[0][2]);
		m_pData[0][0] *= invd;
		m_pData[0][1] *= invd;
		m_pData[0][2] *= invd;

		invd = one/sqrt(m_pData[1][0]*m_pData[1][0] + m_pData[1][1]*m_pData[1][1] + m_pData[1][2]*m_pData[1][2]);
		m_pData[1][0] *= invd;
		m_pData[1][1] *= invd;
		m_pData[1][2] *= invd;

		invd = one/sqrt(m_pData[2][0]*m_pData[2][0] + m_pData[2][1]*m_pData[2][1] + m_pData[2][2]*m_pData[2][2]);
		m_pData[2][0] *= invd;
		m_pData[2][1] *= invd;
		m_pData[2][2] *= invd;
	}

	void NormalizeCol()
	{
		T one(1);
		T invd = one/sqrt(m_pData[0][0]*m_pData[0][0] + m_pData[1][0]*m_pData[1][0] + m_pData[2][0]*m_pData[2][0]);
		m_pData[0][0] *= invd;
		m_pData[1][0] *= invd;
		m_pData[2][0] *= invd;

		invd = one/sqrt(m_pData[0][1]*m_pData[0][1] + m_pData[1][1]*m_pData[1][1] + m_pData[2][1]*m_pData[2][1]);
		m_pData[0][1] *= invd;
		m_pData[1][1] *= invd;
		m_pData[2][1] *= invd;

		invd = one/sqrt(m_pData[0][2]*m_pData[0][2] + m_pData[1][2]*m_pData[1][2] + m_pData[2][2]*m_pData[2][2]);
		m_pData[0][2] *= invd;
		m_pData[1][2] *= invd;
		m_pData[2][2] *= invd;
	}

	T Trace() const
	{
		return (m_pData[0][0] + m_pData[1][1] + m_pData[2][2]);
	}

	// Determinant 
	T Det() const
	{
		return (m_pData[0][0]*(m_pData[1][1]*m_pData[2][2] - m_pData[1][2]*m_pData[2][1]) + 
				m_pData[0][1]*(m_pData[1][2]*m_pData[2][0] - m_pData[1][0]*m_pData[2][2]) + 
				m_pData[0][2]*(m_pData[1][0]*m_pData[2][1] - m_pData[1][1]*m_pData[2][0]));
	}

	int SVD(UpixMatrix3x3<T> &matU, UpixMatrix3x3<T> &matS, UpixMatrix3x3<T> &matVT, T eps) const
	{
		uint i;
		
		T u[9] = {0};
		T v[9] = {0};
		T a[9] = {0};
		T* pData = (T *)m_pData;
		for (i=0; i<9; i++)
		{
			a[i] = pData[i];
		}

		T s[4] = {0};
		T e[4] = {0};
		T w[4] = {0};

		uint it = 60; 
		const uint k = 2;
		const uint l= 1;
		uint ll = 2;

		T d = 0;
		uint kk, ix, iy, iz;
		for (kk=1; kk<=ll; kk++) // 1~2
		{ 
			d = 0.0f;
			for (i=kk; i<=3; i++) // (1~3)|(2~3)
			{
				ix = 3*i - 4 + kk; 
				d += a[ix]*a[ix];
			}
			s[kk-1] = sqrt(d);

			if (!IsEqualZero(s[kk-1]))
			{
				ix = 4*kk - 4;
				if (!IsEqualZero(a[ix]))
				{
					s[kk - 1] = fabs(s[kk - 1]);
					if (a[ix] < 0.0f)
					{
						s[kk - 1] = - s[kk - 1];
					}
				}
				for (i=kk; i<=3; i++)
				{
					iy = 3*i - 4 + kk;
					a[iy] = a[iy]/s[kk - 1];
				}
				a[ix] += 1.0f;
			}
			s[kk - 1] = -s[kk - 1];

			for (uint j=kk+1; j<=3; j++) // (2~3)~3
			{ 
				if ((kk<=k) && (!IsEqualZero(s[kk-1])))
				{ 
					d = 0.0f;
					for (i=kk; i<=3; i++)
					{
						ix = 3*i - 4 + kk;
						iy = 3*i - 4 + j;
						d += a[ix]*a[iy];
					}
					d /= -a[4*kk - 4];
					for (uint i=kk; i<=3; i++)
					{ 
						ix = 3*i - 4 + j;
						iy = 3*i - 4 + kk;
						a[ix] += d*a[iy];
					}
				}
				e[j - 1] = a[3*kk - 4 + j];
			}
		
			for (i=kk; i<=3; i++)
			{ 
				ix = 3*i - 4 + kk;
				u[ix] = a[ix];
			}

			if (kk <= l) // (1~2) <=1
			{
				d = 0.0f;
				for (i=kk+1; i<=3; i++)
				{
					d += e[i-1]*e[i-1];
				}
				e[kk-1] = sqrt(d);
				if (!IsEqualZero(e[kk-1]))
				{
					if (!IsEqualZero(e[kk]))
					{ 
						e[kk-1] = fabs(e[kk-1]);
						if (e[kk] < 0.0) 
						{
							e[kk-1] = -e[kk-1];
						}
					}
					for (i=kk+1; i<=3; i++)
					{
						e[i-1] /= e[kk-1];
					}
					e[kk] += 1.0f;
				}
				e[kk-1] = -e[kk-1];

				if ((kk+1<=3)&&(!IsEqualZero(e[kk-1])))
				{ 
					for (i=kk+1; i<=3; i++)
					{
						w[i-1] = 0.0f;
					}
					for (uint j=kk+1; j<=3; j++)
					{
						for (i=kk+1; i<=3; i++)
						{
							w[i-1] += e[j-1]*a[3*i - 4 + j];
						}
					}
					for (uint j=kk+1; j<=3; j++)
					{
						for (i=kk+1; i<=3; i++)
						{ 
							ix = 3*i - 4 + j;
							a[ix] -= w[i-1]*e[j-1]/e[kk];
						}
					}
				}
				for (i=kk+1; i<=3; i++)
				{
					v[3*i - 4 + kk] = e[i-1];
				}
			}
		}

		s[k] = a[k*4];
		e[l] = a[l*3 + 2];
		e[2] = 0.0f;
		
		for (i=1; i<=3; i++) // 1~3
		{
			u[3*i - 1] = 0.0f;
		}
		u[8] = 1.0f;
		
		for (ll=1; ll<=k; ll++) // 1~2
		{
			kk = k - ll + 1; // 1~2
			iz = 4*kk - 4;
			if (!IsEqualZero(s[kk-1]))
			{
				for (uint j=kk+1; j<=3; j++) // (2~3)~3
				{
					d = 0.0f;
					for (i=kk; i<=3; i++) // (1~2)~3
					{
						ix = 3*i - 4 + kk;
						iy = 3*i - 4 + j;
						d += u[ix]*u[iy]/u[iz];
					}
					d = -d;
					for (i=kk; i<=3; i++)
					{ 
						ix = 3*i - 4 + j;
						iy = 3*i - 4 + kk;
						u[ix] += d*u[iy];
					}
				}
	
				for (i=kk; i<=3; i++)
				{
					ix = 3*i - 4 + kk; 
					u[ix] = -u[ix];
				}
				u[iz] += 1.0f;
				if (kk-1 >= 1)
				{
					for (i=1; i<=kk-1; i++)
					{
						u[3*i - 4 + kk] = 0.0f;
					}
				}
			}
			else
			{ 
				for (i=1; i<=3; i++)
				{
					u[3*i - 4 + kk] = 0.0f;
				}
				u[4*kk - 4] = 1.0f;
			}
		}

		for (ll=1; ll<=3; ll++) // 1~3
		{ 
			kk = 4 - ll; // !~3
			iz = kk*4 - 1;
			if ((kk<=l) && (!IsEqualZero(e[kk-1])))
			{ 
				for (uint j=kk+1; j<=3; j++) // (2~4)~3
				{ 
					d = 0.0f;
					for (i=kk+1; i<=3; i++)
					{ 
						ix = 3*i - 4 + kk; 
						iy = 3*i - 4 + j;
						d += v[ix]*v[iy]/v[iz];
					}
					d = -d;
					for (i=kk+1; i<=3; i++)
					{
						ix = 3*i - 4 + j;
						iy = 3*i - 4 + kk;
						v[ix] += d*v[iy];
					}
				}
			}
			for (i=1; i<=3; i++)
			{
				v[3*i - 4 + kk] = 0.0f;
			}
			v[iz - 3] = 1.0f;
		}

		memset(a, 0, sizeof(T)*9);

		uint mm = 3;
		uint ks; 
		it = 60;
		T dd, t, sm, sm1, em1, sk, ek, b, c, shh;
		T fg[2], cs[2];
		while (1)
		{
			if (mm == 0)
			{
				ppp(a, e, s, v);
				matS.SetData(a);
				matU.SetData(u);
				matVT.SetData(v);
				return 1;
			}
			if (it == 0)
			{ 
				ppp(a, e, s, v);
				matS.SetData(a);
				matU.SetData(u);
				matVT.SetData(v);

				assert(false);//¹¤×÷Ê§°Ü
				return -1;
			}
			kk = mm - 1;
			while ((kk!=0) && (!IsEqualZero(e[kk-1])))
			{
				d = fabs(s[kk-1]) + fabs(s[kk]);
				dd = fabs(e[kk-1]);
				if (dd > eps*d)
				{
					kk--;
				}
				else
				{
					e[kk-1] = 0.0f;
				}
			}

			if (kk == mm - 1)
			{
				kk++;
				if (s[kk-1] < 0.0f)
				{
					s[kk-1] = -s[kk-1];
					for (i=1; i<=3; i++)
					{ 
						ix = 3*i - 4 + kk;
						v[ix] = -v[ix];
					}
				}
				while ((kk!=3) && (s[kk-1]<s[kk]))
				{
					d = s[kk-1]; s[kk-1] = s[kk]; s[kk] = d;
					if (kk < 3)
					{
						for (i=1; i<=3; i++)
						{ 
							ix = 3*i - 4 + kk; 
							iy = 3*i - 3 + kk;
							d = v[ix]; v[ix] = v[iy]; v[iy] = d;
						}
					}

					if (kk < 3)
					{
						for (i=1; i<=3; i++)
						{ 
							ix = 3*i - 4 + kk; 
							iy = 3*i - 3 + kk;
							d = u[ix]; u[ix] = u[iy]; u[iy] = d;
						}
					}
					kk++;
				}
				it = 60;
				mm--;
			}
			else
			{
				ks = mm;
				while ((ks>kk) && (!IsEqualZero(s[ks-1])))
				{
					d = 0.0f;
					if (ks != mm) 
					{
						d += fabs(e[ks-1]);
					}
					if (ks != kk+1)
					{
						d += fabs(e[ks-2]);
					}
					dd = fabs(s[ks-1]);
					if (dd > eps*d)
					{
						ks--;
					}
					else
					{
						s[ks-1] = 0.0f;
					}
				}
				if (ks == kk)
				{ 
					kk++;
					d = fabs(s[mm-1]);
					t = fabs(s[mm-2]);
					if (t > d)
					{
						d = t;
					}
					t = fabs(e[mm-2]);

					if (t > d)
					{
						d = t;
					}
					t = fabs(s[kk-1]);

					if (t > d)
					{
						d = t;
					}
					t = fabs(e[kk-1]);
					if (t > d) 
					{
						d = t;
					}
					sm = s[mm-1]/d; 
					sm1 = s[mm-2]/d;
					em1 = e[mm-2]/d;
					sk = s[kk-1]/d; 
					ek = e[kk-1]/d;
					b = ((sm1 + sm)*(sm1 - sm) + em1*em1)/2.0f;
					c = sm*em1; 
					c = c*c;
					shh = 0.0f;
					if ( ( !IsEqualZero(b) ) || ( !IsEqualZero(c) ) )
					{
						shh = sqrt(b*b + c);
						if (b < 0.0f)
						{
							shh = -shh;
						}
						shh = c/(b + shh);
					}
					fg[0] = (sk + sm)*(sk - sm) - shh;
					fg[1] = sk*ek;
					for (i=kk; i<=mm-1; i++)
					{ 
						sss(fg, cs);
						if (i != kk) 
						{
							e[i-2] = fg[0];
						}
						fg[0] = cs[0]*s[i-1] + cs[1]*e[i-1];
						e[i-1] = cs[0]*e[i-1] - cs[1]*s[i-1];
						fg[1] = cs[1]*s[i];
						s[i] = cs[0]*s[i];

						if ((!IsEqualZero(cs[0]-1.0f)) || (!IsEqualZero(cs[1])))
						{
							for (uint j=1; j<=3; j++)
							{
								ix = 3*j - 4 + i;
								iy = 3*j - 3 + i;
								d = cs[0]*v[ix] + cs[1]*v[iy];
								v[iy] = -cs[1]*v[ix] + cs[0]*v[iy];
								v[ix] = d;
							}
						}
						sss(fg, cs);
						s[i-1] = fg[0];
						fg[0] = cs[0]*e[i-1] + cs[1]*s[i];
						s[i] = -cs[1]*e[i-1] + cs[0]*s[i];
						fg[1] = cs[1]*e[i];
						e[i] = cs[0]*e[i];

						if (i < 3)
						{
							if ((!IsEqualZero(cs[0]-1.0f)) || (!IsEqualZero(cs[1])))
							{
								for (uint j=1; j<=3; j++)
								{
									ix = 3*j - 4 + i;
									iy = 3*j - 3 + i;
									d = cs[0]*u[ix] + cs[1]*u[iy];
									u[iy] = -cs[1]*u[ix] + cs[0]*u[iy];
									u[ix] = d;
								}
							}
						}
					}
					e[mm-2] = fg[0];
					it--;
				}
				else
				{ 
					if (ks == mm)
					{ 
						kk++;
						fg[1] = e[mm-2]; 
						e[mm-2] = 0.0f;
						for (ll=kk; ll<=mm-1; ll++)
						{ 
							i = mm + kk - ll - 1;
							fg[0] = s[i-1];
							sss(fg,cs);
							s[i-1] = fg[0];

							if (i != kk)
							{ 
								fg[1] = -cs[1]*e[i-2];
								e[i-2] = cs[0]*e[i-2];
							}
							if ((!IsEqualZero(cs[0]-1.0f)) || (!IsEqualZero(cs[1])))
							{
								for (uint j=1; j<=3; j++)
								{
									ix = 3*j - 4 + i;
									iy = 3*j - 4 + mm;
									d = cs[0]*v[ix] + cs[1]*v[iy];
									v[iy] = -cs[1]*v[ix] + cs[0]*v[iy];
									v[ix] = d;
								}
							}
						}
					}
					else
					{ 
						kk = ks + 1;
						fg[1] = e[kk-2];
						e[kk-2] = 0.0f;

						for (i=kk; i<=mm; i++)
						{
							fg[0] = s[i-1];
							sss(fg, cs);
							s[i-1] = fg[0];
							fg[1] = -cs[1]*e[i-1];
							e[i-1] = cs[0]*e[i-1];

							if ((!IsEqualZero(cs[0]-1.0f)) || (!IsEqualZero(cs[1])))
							{
								for (uint j=1; j<=3; j++)
								{ 
									ix = 3*j - 4 + i;
									iy = 3*j - 5 + kk;
									d = cs[0]*u[ix] + cs[1]*u[iy];
									u[iy] = -cs[1]*u[ix] + cs[0]*u[iy];
									u[ix] = d;
								}
							}
						}
					}
				}
			}
		}
	}

	UpixMatrix3x3<T> Transpose() const
	{
		UpixMatrix3x3<T> result;
		Transpose(result);
		return result;
	}

	inline void Transpose(UpixMatrix3x3<T> &ret) const
	{
		T *pMr = (T *)ret.m_pData;
		pMr[0] = m_pData[0][0]; pMr[1] = m_pData[1][0]; pMr[2] = m_pData[2][0];
		pMr[3] = m_pData[0][1]; pMr[4] = m_pData[1][1]; pMr[5] = m_pData[2][1];
		pMr[6] = m_pData[0][2]; pMr[7] = m_pData[1][2]; pMr[8] = m_pData[2][2];
	}

	UpixMatrix3x3<T> Inv() const
	{
		UpixMatrix3x3<T> ret;
		Inv(ret);
		return ret;
	}

	T Inv(UpixMatrix3x3<T> &m) const
	{
		T *pm = (T *)m.m_pData;
		
		pm[0] = m_pData[1][1] * m_pData[2][2] - m_pData[1][2] * m_pData[2][1];
		pm[1] = m_pData[0][2] * m_pData[2][1] - m_pData[0][1] * m_pData[2][2];
		pm[2] = m_pData[0][1] * m_pData[1][2] - m_pData[0][2] * m_pData[1][1];

		T det = m_pData[0][0] * pm[0] + m_pData[1][0] * pm[1] + m_pData[2][0] * pm[2];
		if (UpixAbs(det) < UPIX_FLOAT_EPSILON)
		{
			return det;
		}
		T one(1);
		T dInvDet = one / det;

		pm[0] *= dInvDet; 
		pm[1] *= dInvDet; 
		pm[2] *= dInvDet;
		pm[3] = (m_pData[1][2] * m_pData[2][0] - m_pData[1][0] * m_pData[2][2]) * dInvDet;
		pm[4] = (m_pData[0][0] * m_pData[2][2] - m_pData[0][2] * m_pData[2][0]) * dInvDet;
		pm[5] = (m_pData[0][2] * m_pData[1][0] - m_pData[0][0] * m_pData[1][2]) * dInvDet;
		pm[6] = (m_pData[1][0] * m_pData[2][1] - m_pData[1][1] * m_pData[2][0]) * dInvDet;
		pm[7] = (m_pData[0][1] * m_pData[2][0] - m_pData[0][0] * m_pData[2][1]) * dInvDet;
		pm[8] = (m_pData[0][0] * m_pData[1][1] - m_pData[0][1] * m_pData[1][0]) * dInvDet;

		return det;
	}

	// Get normal infinite, computes the maximum of the sums of absolute values over rows
	T InfiniteNorm() const
	{
		T sum1	= abs(m_pData[0][0]) + abs(m_pData[0][1]) + abs(m_pData[0][2]);
		T sum2	= abs(m_pData[1][0]) + abs(m_pData[1][1]) + abs(m_pData[1][2]);
		T sum3	= abs(m_pData[2][0]) + abs(m_pData[2][1]) + abs(m_pData[2][2]);
		return UpixMax(UpixMax(sum1, sum2), sum3);
	}

	// Computes the matrix exponential of a matrix m by scaling m by 1/(powers of 2), 
	// using Taylor series and squaring again. The input matrix must be square
	UpixMatrix3x3<T> Exp() const
	{
		UpixMatrix3x3<T> result;
		Exp(result);
		return result;
	}

	void Exp(UpixMatrix3x3<T> &result) const
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
	static void exp_taylor(const UpixMatrix3x3<T> &m, UpixMatrix3x3<T> &result)
	{
		result.SetValue(0);
		UpixMatrix3x3<T> f;
		T *pMf = f.GetData();
		pMf[0] = 1;
		pMf[4] = 1;
		pMf[8] = 1;

		T k = 1;
		while (f.InfiniteNorm() > 0)
		{
			result += f;
			f *= m;
			f /= k;
			k += 1;
		}
	}

	inline bool IsEqualZero(T _num) const
	{
		double eps = 1e-16f;
		return ( (_num > - eps) && (_num < eps) );
	}

	void ppp(T a[9], T e[4], T s[4], T v[9]) const
	{
		uint i = 3, j, p, q;
		T d;

		for (j=1; j<=i-1; j++)
		{ 
			a[4*j - 4] = s[j - 1];
			a[4*j - 3] = e[j-1];
		}
		a[4*i - 4] = s[i - 1];

		for (i=1; i<=2; i++)
		{
			for (j=i+1; j<=3; j++)
			{
				p = 3*i - 4 + j; 
				q = 3*j - 4 + i;
				d = v[p]; 
				v[p] = v[q];
				v[q] = d;
			}
		}
	}

	void sss(T fg[2], T cs[2]) const
	{
		T r,d;
		T temp1, temp2;
		temp1 = fabs(fg[0]);
		temp2 = fabs(fg[1]);
		if (IsEqualZero(temp1 + temp2))
		{
			cs[0] = 1.0f;
			cs[1] = 0.0f;
			d = 0.0f;
		}
		else 
		{ 
			d = sqrt(fg[0]*fg[0] + fg[1]*fg[1]);

			if (temp1 > temp2)
			{ 
				d = fabs(d);

				if (fg[0] < 0.0f)
				{
					d = -d;
				}
			}
			else
			{
				d = fabs(d);
				if (fg[1] < 0.0f)
				{
					d = -d;
				}
			}

			cs[0] = fg[0]/d; cs[1] = fg[1]/d;
		}

		r = 1.0f;
		if (temp1 > temp2) 
		{
			r = cs[1];
		}
		else
		{
			if (!IsEqualZero(cs[0]))
			{
				r = 1.0f/cs[0];
			}
		}
		fg[0] = d; fg[1] = r;
	}

private:
	T m_pData[3][3];
};



typedef UpixMatrix3x3<float> UpixMatrix3x3f;
typedef UpixMatrix3x3<double> UpixMatrix3x3d;



template <class T>
void UPIX_ALGORITHM_TEMPLATE MultVectorAsMatrix(const UpixPoint3D<T> &v1, const UpixPoint3D<T> &v2, UpixMatrix3x3<T> &ret)
{
	ret(0, 0) = v1.x*v2.x; ret(0, 1) = v1.x*v2.y; ret(0, 2) = v1.x*v2.z;
	ret(1, 0) = v1.y*v2.x; ret(1, 1) = v1.y*v2.y; ret(1, 2) = v1.y*v2.z;
	ret(2, 0) = v1.z*v2.x; ret(2, 1) = v1.z*v2.y; ret(2, 2) = v1.z*v2.z;
}




#endif