/*
* Filename: UpixMatrix.h
* Author:   Xiehong
* Date:     2012.9.29
*/
// This class is not a container, only base types are supported.

#ifndef _UPIX_MATRIX_H_
#define _UPIX_MATRIX_H_

#include "UpixVector.h"
#include <math.h>
#include "UpixMemory.h"


template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixMatrix
{
public:
	//default constructor
	UpixMatrix() { m_pData = NULL; m_nRow = 0; m_nCol = 0; }

	//construct by size
	UpixMatrix(uint nRow, uint nCol)
	{
		m_nRow = nRow;
		m_nCol = nCol;
		uint n = m_nRow*m_nCol;
		m_pData = new T[n];
		int val = 0;
		Fill(m_pData, n, val);
	}

	UpixMatrix(uint nRow, uint nCol, const T &val)
	{
		m_nRow = nRow;
		m_nCol = nCol;
		uint n = m_nRow*m_nCol;
		m_pData = new T[n];
		Fill(m_pData, n, val);
	}

	//construct from array
	UpixMatrix(uint nRow, uint nCol, T value[]); 

	//construct from vector
	UpixMatrix(uint nRow, uint nCol, const UpixVector<T> &vec);

	//construct by n*n size
	UpixMatrix(uint nSize)
	{
		m_nRow = nSize;
		m_nCol = nSize;
		uint n = m_nRow*m_nCol;
		m_pData = new T[n];
		int val = 0;
		Fill(m_pData, n, val);
	}

	//construct from array, n*n size
	UpixMatrix(uint nSize, T value[]);

	// Construct from vector in the Lie algebra
	UpixMatrix(uint nSize, const UpixVector<T> &v)
	{
		uint dim = nSize*nSize - 1;
		assert(v.Size() == dim);

		UpixMatrix<T> t(nSize, nSize, 0);
		UpixMatrix<T> base(nSize);
		for(int i=0; i<dim; ++i)
		{
			Generator(i, base);
			t += base * v[i];
		}
		*this = Exp(t);
	}

	//construct from other
	UpixMatrix(const UpixMatrix<T> &other);

	//operator = 
	UpixMatrix<T>& operator=(const UpixMatrix<T> &rhs);

	//default destructor
	~UpixMatrix() { Clear(); }

public:
	//transpose
	UpixMatrix<T> Transpose() const
	{
		UpixMatrix<T> result(m_nCol, m_nRow);
		Transpose(result);
		return result;
	}

	void Transpose(UpixMatrix<T> &m) const
	{
		assert(m.m_nRow == m_nCol && m.m_nCol == m_nRow);

		T *pDst = m.m_pData;
		T *pSrc = NULL;
		for(uint i=0; i<m_nCol; i++)
		{
			pSrc = m_pData + i;
			for(uint j=0; j<m_nRow; j++)
			{
				*pDst++ = *pSrc;
				pSrc += m_nCol;
			}
		}
	}

	//change size
	void Resize(uint nRow, uint nCol);

	//set size to n*n identity matrix
	void SetIdentity(uint nSize);

	// set zeros, n*n
	void SetZeros(uint nSize) { return SetZeros(nSize, nSize); }

	// set zeros, row*col
	void SetZeros(uint nRow, uint nCol);

	// set diag from vector
	void SetDiag(const UpixVector<T> &vec)
	{
		uint n = vec.Size();
		SetZeros(n);
		T *p = vec.GetData();
		for (uint i=0; i<n; i++)
		{
			m_pData[i*m_nCol+i] = p[i];
		}
	}

	//set data from array
	void SetData(uint nRow, uint nCol, T value[]);

	//set matrix from other
	void SetMatrix(const UpixMatrix<T> &M);

	//set element at (r,c)
	void SetElement(uint nRow, uint nCol, const T &value);

	// set column data
	void SetColumn(uint nCol, const UpixVector<T> &vec);

	//set row data
	void SetRow(uint nRow, const UpixVector<T> &vec)
	{
		SetRow(nRow, vec.GetData(), vec.Size());
	}

	void SetRow(uint nRow, T *pData, uint len)
	{
		assert((nRow >= 0) && (nRow < m_nRow) && m_nCol == len);
		Fill(m_pData + nRow*m_nCol, len, pData);
	}

	//get element at (r,c)
	T GetElement(uint nRow, uint nCol) const;

	//get column size
	inline uint GetColNum() const { return m_nCol; }

	//get row size
	inline uint GetRowNum() const { return m_nRow; }

	//get row data
	UpixVector<T> GetRowVector(uint nRow) const
	{
		assert(nRow >= 0 && nRow < m_nRow);
		return UpixVector<T>(m_nCol, m_pData + nRow*m_nCol);
	}

	void GetRowVector(uint nRow, UpixVector<T> &v) const
	{
		assert(nRow >= 0 && nRow < m_nRow);
		v.SetData(m_nCol, m_pData + nRow*m_nCol);
	}

	//get column data
	UpixVector<T> GetColVector(uint nCol) const
	{
		UpixVector<T> vec(m_nRow);
		GetColVector(nCol, vec);
		return vec;
	}

	void GetColVector(uint nCol, UpixVector<T> &v) const
	{
		assert(nCol >= 0 && nCol < m_nCol);
		T *p = m_pData + nCol;
		for(uint j = 0; j< m_nRow; ++j)
		{
			v[j] = *p;
			p += m_nCol;
		}
	}

	// get diag element
	UpixVector<T> GetDiag() const
	{
		UpixVector<T> vec;
		GetDiag(vec);
		return vec;
	}

	void GetDiag(UpixVector<T> &vec) const
	{
		uint n = UpixMin(m_nRow, m_nCol);
		if (vec.Size() != n)
		{
			vec.Resize(n);
		}

		T *p = vec.GetData();
		for (uint i=0; i<n; i++)
		{
			p[i] = m_pData[i*m_nCol+i];
		}
	}

	// get data buffer
	inline T* GetData() const { return m_pData; }

	//clear buffer
	void Clear();

	//normalization by row
	void NormalizeRow();

	//normalization by column
	void NormalizeCol();

	// Trace
	T Trace() const
	{
		uint n = UpixMin(m_nRow, m_nCol);

		T tr = 0;
		for (uint i=0; i<n; i++)
		{
			tr += m_pData[i*m_nCol+i];
		}

		return tr;
	}

	// Determinant 
	double Det() const;

	// Inverse
	UpixMatrix<T> Inv() const
	{
		UpixMatrix<T> ret(m_nRow);
		Inv(ret);
		return ret;	
	}

	bool Inv(UpixMatrix<T> &ret) const
	{
		assert(m_nRow == m_nCol);

		UpixMatrix<T> U, D, VT;
		if (SVD(U, D, VT, (T)1e-6) == -1)
		{
			return false;
		}

		UpixVector<T> r;
		D.GetDiag(r);
		T *p = r.GetData();
		uint n = r.Size();
		uint iCount = 0;
		while (iCount < n && *p != 0)
		{
			*p = 1/(*p++);
			iCount++;
		}
		D.SetDiag(r);

		// ret = VT.Transpose()*D*U.Transpose()
		MultMatrix(VT.Transpose()*D, U.Transpose(), ret);
		return true;
	}

	/*
	全选主元高斯-约旦法。
	高斯-约旦法（全选主元）求逆的步骤如下：
	首先，对于 k 从 0 到 n - 1 作如下几步：

	1. 从第 k 行、第 k 列开始的右下角子阵中选取绝对值最大的元素，并记住次元素所在的行号和列号，在通过行交换和列交换将它交换到主元素位置上。这一步称为全选主元。
	2. m(k, k) = 1 / m(k, k)
	3. m(k, j) = m(k, j) * m(k, k)，j = 0, 1, ..., n-1；j != k
	4. m(i, j) = m(i, j) - m(i, k) * m(k, j)，i, j = 0, 1, ..., n-1；i, j != k
	5. m(i, k) = -m(i, k) * m(k, k)，i = 0, 1, ..., n-1；i != k
	6. 最后，根据在全选主元过程中所记录的行、列交换的信息进行恢复，恢复的原则如下：在全选主元过程中，先交换的行（列）后进行恢复；原来的行（列）交换用列（行）交换来恢复。
	*/
	UpixMatrix<T> InvGaussian() const
	{
		UpixMatrix<T> ret(m_nRow);
		InvGaussian(ret);
		return ret;
	}

	T InvGaussian(UpixMatrix<T> &m) const
	{
		assert(m_nRow == m_nCol);

		m = *this;

		int *js = new int[m_nRow*2];
		int *is = js + m_nRow;
		// 	int *js = new int[m_nRow];
		// 	int *is = new int[m_nRow];

		double fDet = 1.0f;
		int f = 1;

		for(uint k=0; k<m_nRow; k++)
		{
			// 1. 全选主元
			double fMax = 0.0f;
			for (uint i=k; i<m_nRow; i++)
			{
				for (uint j=k; j<m_nRow; j++)
				{
					const double f = abs(m(i, j));
					if (f > fMax)
					{
						fMax	= f;
						is[k]	= i;
						js[k]	= j;
					}
				}
			}
#ifdef __APPLE__
			if (fabs(fMax) < UPIX_DOUBLE_EPSILON)
#else
            if (abs(fMax) < UPIX_DOUBLE_EPSILON)
#endif
			{
				delete[] js;
				return 0;
			}

			if (is[k] != k)
			{
				f = -f;
				for (uint i=0; i<m_nRow; i++)
				{
					Swap(m(k, i), m(is[k], i));
				}
			}

			if (js[k] != k)
			{
				f = -f;
				for (uint i=0; i<m_nRow; i++)
				{
					Swap(m(i, k), m(i, js[k]));
				}
			}

			// 计算行列值
			fDet *= m(k, k);

			// 计算逆矩阵

			// 2. 
			m(k, k) = 1.0f / m(k, k);

			// 3.
			for(uint j=0; j<m_nRow; j++)
			{
				if (j != k)
				{
					m(k, j) *= m(k, k);
				}
			}

			// 4. 
			for (uint i=0; i<m_nRow; i++)
			{
				if (i != k)
				{
					for (uint j=0; j<m_nRow; j++)
					{
						if (j != k)
						{
							m(i, j) = m(i, j) - m(i, k) * m(k, j);
						}
					}
				}
			}

			// 5.
			for (uint i=0; i<m_nRow; i++)
			{
				if (i != k)
				{
					m(i, k) *= -m(k, k);
				}
			}
		}

		// 6.
		for (int k=((int)m_nRow-1); k>=0; k--)
		{
			if (js[k] != k)
			{
				for (uint j=0; j<m_nRow; j++)
				{
					Swap(m(k, j), m(js[k], j));
				}
			}

			if (is[k] != k)
			{
				for (uint i=0; i<m_nRow; i++)
				{
					Swap(m(i, k), m(i, is[k]));
				}
			}
		}

		/*	delete[] is;*/
		delete[] js;
		return fDet * f;
	}

	// get rank
	int Rank();

	// SVD decompose                                                                   
	// u:   存放左奇异向量U        (大小m*m)                                                          
	// s:   对角线存放奇异值       (大小m*n)
	// vt:  存放右奇异向量V的转置V'(大小n*n)                                                        
	// eps: 给定的精度要求         (一般0.000001)
	int SVD(UpixMatrix<T>& matU, UpixMatrix<T>& matS, UpixMatrix<T>& matVT, T eps = 1e-6) const;

	//swap mat
	bool Swap(UpixMatrix<T> &M);

	// compute gradient
	void Gradient(UpixMatrix<double> &Fx, UpixMatrix<double> &Fy);

	// Get normal infinite, computes the maximum of the sums of absolute values over rows
	T InfiniteNorm() const
	{
		T norm = 0;
		T sum = 0;
		T *p = m_pData;
		for(uint i=0; i<m_nRow; ++i)
		{
			sum = 0;
			for(uint j=0; j<m_nCol; ++j)
			{
				sum += abs(*p);
				p++;
			}
			norm = UpixMax(norm, sum);
		}
		return norm;
	}

	// Frobenius (root of sum of squares) norm of matrix
	T FrobeniusNorm() const
	{
		T norm = 0;
		T *p = m_pData;
		for(uint i=0; i<m_nRow; ++i)
		{
			for(uint j=0; j<m_nCol; ++j)
			{
				norm += p[0]*p[0];
				p++;
			}
		}
		return sqrt(norm);
	}

	// Computes the matrix exponential of a matrix m by scaling m by 1/(powers of 2), 
	// using Taylor series and squaring again. The input matrix must be square
	UpixMatrix<T> Exp() const
	{
		UpixMatrix<T> result;
		Exp(result);
		return result;
	}

	void Exp(UpixMatrix<T> &result) const
	{
		assert(m_nRow == m_nCol);

		T level = log10(InfiniteNorm())*UPIX_INV_LOG10_2;	//const P l = log2(norm_inf(m));
		int scale = UpixMax(0, (int)ceil(level));

		exp_taylor((*this)/(1<<scale), result);
		for(int i=0; i<scale; ++i)
		{
			result *= result;
		}
	}

protected:
	bool IsEqualZero(T _num) const;
	void sss(T fg[2], T cs[2]) const;
	void ppp(T a[], T e[], T s[], T v[], uint m, uint n) const;

	// Exponentiate a matrix using a the Taylor series
	// This will not work if the norm of the matrix is too large.
	static void exp_taylor(const UpixMatrix<T> &m, UpixMatrix<T> &result)
	{
		int nRows = m.m_nRow;
		int nCols = m.m_nCol;
		assert(nRows == nCols);

		result.SetZeros(nRows);

		UpixMatrix<T> f;
		f.SetIdentity(nRows);

		T k = 1;
		while(f.InfiniteNorm() > 0)
		{
			result += f;
			f *= m;
			f /= k;
			k += 1;
		}
	}

	// returns one generator of the SL group.
	void Generator(int i, UpixMatrix<T> &matBase)
	{
		int dim = m_nRow*m_nCol - 1;
		assert(m_nRow == m_nCol && i > -1 && i < dim);
		matBase.SetZeros(m_nRow);

		int DIAG_LIMIT = m_nRow - 1;
		int SYMM_LIMIT = (dim - DIAG_LIMIT)/2 + DIAG_LIMIT;

		if(i < DIAG_LIMIT) // first ones are the diagonal ones
		{		
			matBase(i, i) = 1;
			matBase(i+1, i+1) = -1;
		}
		else if(i < SYMM_LIMIT)	// then the symmetric ones
		{
			int row = 0, col = i - DIAG_LIMIT + 1;
			while(col > (m_nRow - row - 1))
			{
				col -= (m_nRow - row - 1); 
				++row;
			}
			col += row;
			matBase(row, col) = matBase(col, row) = 1;
		} 
		else	// finally the antisymmetric ones
		{
			int row = 0, col = i - SYMM_LIMIT + 1;
			while(col > m_nRow - row - 1){
				col -= m_nRow - row - 1; 
				++row;
			}
			col += row;
			matBase(row, col) = -1;
			matBase(col, row) = 1;
		}
	}

	inline void Swap(T &a, T &b)
	{
		T c;
		c = a;
		a = b;
		b = c;
	}


public:
	//operator [] const
	T* operator[](uint i) const { return m_pData + i*m_nCol; }

	//operator ()
	T& operator()(uint row, uint col)
	{
		assert(row>=0 && row<m_nRow && col>=0 && col<m_nCol);
		return m_pData[row*m_nCol + col];
	}

	const T& operator()(uint row, uint col) const
	{
		assert(row>=0 && row<m_nRow && col>=0 && col<m_nCol);
		return m_pData[row*m_nCol + col];
	}

	//operator ==
	bool operator==(const UpixMatrix<T>& other) const;

	//operator +
	UpixMatrix<T> operator+(const UpixMatrix<T>& other) const;

	//operator -
	UpixMatrix<T> operator-(const UpixMatrix<T>& other) const;

	//operator +
	UpixMatrix<T> operator+(const T &value) const;

	//operator -
	UpixMatrix<T> operator-(const T &value) const;

	//operator *
	UpixMatrix<T> operator*(const T &value) const;

	//operator *
	UpixMatrix<T> operator*(const UpixMatrix<T>& other) const;

	//operator *
	UpixVector<T> operator*(const UpixVector<T> &other) const; 

	//operator / 
	UpixMatrix<T> operator/(const T &value) const; 

	//operator ~, transpose
	UpixMatrix<T> operator~(void) const { return Transpose(); }

	//operator +=
	UpixMatrix<T> &operator+=(const UpixMatrix<T> &other);

	//operator -=
	UpixMatrix<T> &operator-=(const UpixMatrix<T> &other);

	// operator +=
	UpixMatrix<T> &operator+=(const T &value);

	// operator -=
	UpixMatrix<T> &operator-=(const T &value);

	//operator *= 
	UpixMatrix<T> &operator*=(const T &value);

	//operator *=
	UpixMatrix<T> &operator*=(const UpixMatrix<T> &other)
	{
		SetMatrix((*this) * other);
		return *this;
	}

	//operator /= 
	UpixMatrix<T> &operator/=(const T &value);

	//change sign
	UpixMatrix<T> operator-() const
	{
		return T(0) - *this;
	}

protected:
	uint m_nRow;    //行数
	uint m_nCol;    //列数
	T *m_pData;
};


//----------------------implement-----------------------
template <class T>
UpixMatrix<T>::UpixMatrix(uint nRow, uint nCol, T value[])
{
	m_nRow = nRow;
	m_nCol = nCol;
	uint n = m_nRow*m_nCol;
	m_pData = new T[n];
	Fill(m_pData, n, value);
}

template <class T>
UpixMatrix<T>::UpixMatrix(uint nRow, uint nCol, const UpixVector<T> &vec)
{
	uint n = nRow*nCol;
	assert(n == vec.Size());

	m_nRow = nRow;
	m_nCol = nCol;
	m_pData = new T[n];
	T *pV = vec.GetData();
	Fill(m_pData, n, pV);
}

template <class T>
UpixMatrix<T>::UpixMatrix(uint nSize, T value[])
{
	m_nRow = nSize;
	m_nCol = nSize;
	uint n = m_nRow*m_nCol;
	m_pData = new T[n];
	Fill(m_pData, n, value);
}

template <class T>
UpixMatrix<T>::UpixMatrix(const UpixMatrix<T> &other)
{
	m_nRow = other.m_nRow;
	m_nCol = other.m_nCol;
	uint n = m_nRow*m_nCol;
	m_pData = new T[n];
	Fill(m_pData, n, other.m_pData);
}

template <class T>
UpixMatrix<T>& UpixMatrix<T>::operator=(const UpixMatrix<T> &rhs)
{
	if(&rhs != this)
	{
		uint n = rhs.m_nRow*rhs.m_nCol;
		if (m_nRow != rhs.m_nRow || m_nCol != rhs.m_nCol)
		{
			Clear();
			m_nRow = rhs.m_nRow;
			m_nCol = rhs.m_nCol;
			m_pData = new T[n];
		}
		Fill(m_pData, n, rhs.m_pData);
	}
	return *this;
}

template <class T>
void UpixMatrix<T>::SetMatrix(const UpixMatrix<T> &M)
{
	uint n = M.m_nRow*M.m_nCol;
	if (m_nRow != M.m_nRow || m_nCol != M.m_nCol)
	{
		Clear();
		m_nRow = M.m_nRow;
		m_nCol = M.m_nCol;
		m_pData = new T[n];
	}
	Fill(m_pData, n, M.m_pData);
}

template <class T>
void UpixMatrix<T>::SetColumn(uint nCol, const UpixVector<T> &vec)
{
	assert((nCol >= 0) && (nCol < m_nCol) );
	uint len = vec.Size();
	assert(m_nRow == len);
	T *pData = m_pData + nCol;
	T *pV = vec.GetData();
	for (uint j=0; j<len; j++)
	{
		*pData = pV[j];
		pData += m_nCol;
	}
}

template <class T>
void UpixMatrix<T>::Resize(uint nRow, uint nCol)
{
	Clear();
	m_nRow = nRow;
	m_nCol = nCol;
	uint n = m_nRow*m_nCol;
	m_pData = new T[n];
	int val = 0;
	Fill(m_pData, n, val);
}

template <class T>
void UpixMatrix<T>::Clear()
{
	if (m_pData != NULL)
	{
		delete [] m_pData;
		m_pData = NULL;
	}
	m_nRow = 0;
	m_nCol = 0;
}

template <class T>
void UpixMatrix<T>::SetIdentity(uint nSize)
{
	SetZeros(nSize, nSize);
	for (uint i=0; i<nSize; i++)
	{
		m_pData[i*m_nCol + i] = 1;
	}
}

template <class T>
void UpixMatrix<T>::SetZeros(uint nRow, uint nCol)
{
	uint n = nRow*nCol;
	if (m_nRow != nRow || m_nCol != nCol)
	{
		Clear();
		m_nRow = nRow;
		m_nCol = nCol;
		m_pData = new T[n];
	}
	int val = 0;
	Fill(m_pData, n, val);
}

template <class T>
void UpixMatrix<T>::SetData(uint nRow, uint nCol, T value[])
{
	uint n = nRow*nCol;
	if (m_nRow!=nRow || m_nCol!=nCol)
	{
		Clear();
		m_nRow = nRow;
		m_nCol = nCol;
		m_pData = new T[n];
	}
	Fill(m_pData, n, value);
}

template <class T>
void UpixMatrix<T>::SetElement(uint nRow, uint nCol, const T &value)
{
	assert((nCol >= 0 && nCol < m_nCol && nRow >= 0 && nRow < m_nRow));
	m_pData[nRow*m_nCol + nCol] = value;
}

template <class T>
T UpixMatrix<T>::GetElement(uint nRow, uint nCol) const
{
	assert(nCol>=0 && nCol<m_nCol && nRow>=0 && nRow<m_nRow);
	return m_pData[nRow*m_nCol + nCol];
}

template <class T>
bool UpixMatrix<T>::operator==(const UpixMatrix<T>& other) const
{
	if(m_nRow != other.m_nRow || m_nCol != other.m_nCol)
	{
		return false;
	}
	T *pMo = other.m_pData;
	uint n = m_nCol*m_nRow;
	for (uint i=0; i<n; i++)
	{
		if (m_pData[i] != pMo[i])
		{
			return false;
		}
	}
	return true;
}

template <class T>
UpixMatrix<T> UpixMatrix<T>::operator+(const UpixMatrix<T>& other) const
{
	assert(m_nRow == other.m_nRow && m_nCol == other.m_nCol);
	UpixMatrix<T> result(m_nRow, m_nCol, m_pData);
	result += other;
	return result;
}

template <class T>
UpixMatrix<T> UpixMatrix<T>::operator-(const UpixMatrix<T>& other) const
{
	assert(m_nRow == other.m_nRow && m_nCol == other.m_nCol);
	UpixMatrix<T> result(m_nRow, m_nCol, m_pData);
	result -= other;
	return result;
}

template <class T>
UpixMatrix<T> UpixMatrix<T>::operator *(const UpixMatrix<T> &other) const
{
	assert(m_nCol == other.m_nRow);
	uint otherColNum = other.m_nCol;
	UpixMatrix<T> result(m_nRow, otherColNum);

	T *pRes = result.m_pData;
	T *pRow = m_pData;
	T *pData = NULL;
	T *pCol = NULL;
	for(uint i=0; i<m_nRow; ++i)
	{
		for(uint j=0; j<otherColNum; ++j)
		{
			*pRes = 0;
			pData = pRow;
			pCol = other.m_pData + j;
			for(uint k=0; k<m_nCol; ++k)
			{
				*pRes += (*pData++) * (*pCol);
				pCol += otherColNum;
			}
			pRes++;
		}
		pRow += m_nCol;
	}
	return result;
}

template <class T>
UpixMatrix<T> UpixMatrix<T>::operator+(const T &value) const
{
	UpixMatrix<T> result(m_nRow, m_nCol, m_pData);
	result += value;
	return result;
}

template <class T>
UpixMatrix<T> UpixMatrix<T>::operator-(const T &value) const
{
	UpixMatrix<T> result(m_nRow, m_nCol, m_pData);
	result -= value;
	return result;
}

template <class T>
UpixMatrix<T> UpixMatrix<T>::operator *(const T &value) const
{
	UpixMatrix<T> result(m_nRow, m_nCol, m_pData);
	result *= value;
	return result;
}

template <class T>
UpixMatrix<T> UpixMatrix<T>::operator /(const T &value) const
{
	UpixMatrix<T> result(m_nRow, m_nCol, m_pData);
	result /= value;
	return result;
}

template <class T>
UpixVector<T> UpixMatrix<T>::operator*(const UpixVector<T> &other) const
{
	assert(m_nCol == other.Size());
	UpixVector<T> result(m_nRow);

	T *pData = m_pData;
	T *pV = other.GetData();
	T *pRes = result.GetData();

	for (uint j=0; j<m_nRow; j++)
	{
		*pRes = 0;
		for (uint i=0; i<m_nCol; i++)
		{
			*pRes += (*pData) * pV[i];
			pData++;
		}
		pRes++;
	}
	return result;
}

template <class T>
UpixMatrix<T> &UpixMatrix<T>::operator+=(const UpixMatrix<T> &other)
{
	assert(m_nRow == other.m_nRow && m_nCol == other.m_nCol);
	T *pMo = other.m_pData;
	uint n = m_nCol*m_nRow;
	for (uint i=0; i<n; i++)
	{
		m_pData[i] += pMo[i];
	}
	return *this;
}

template <class T>
UpixMatrix<T> &UpixMatrix<T>::operator-=(const UpixMatrix<T> &other)
{
	assert(m_nRow == other.m_nRow && m_nCol == other.m_nCol);
	T *pMo = other.m_pData;
	uint n = m_nCol*m_nRow;
	for (uint i=0; i<n; i++)
	{
		m_pData[i] -= pMo[i];
	}
	return *this;
}

template <class T>
UpixMatrix<T> &UpixMatrix<T>::operator+=(const T &value)
{
	uint n = m_nRow*m_nCol;
	for(uint j=0; j<n; j++)
	{
		m_pData[j] += value;
	}
	return *this;
}

template <class T>
UpixMatrix<T> &UpixMatrix<T>::operator-=(const T &value)
{
	uint n = m_nRow*m_nCol;
	for(uint j=0; j<n; j++)
	{
		m_pData[j] -= value;
	}
	return *this;
}

template <class T>
UpixMatrix<T> &UpixMatrix<T>::operator*=(const T &value)
{
	uint n = m_nRow*m_nCol;
	for(uint j=0; j<n; j++)
	{
		m_pData[j] *= value;
	}
	return *this;
}

template <class T>
UpixMatrix<T> &UpixMatrix<T>::operator/=(const T &value)
{
	assert(value != 0);
	T invVal = 1/value;
	int n = m_nRow*m_nCol;
	for(int j=0; j<n; j++)
	{
		m_pData[j] *= invVal;
	}
	return *this;
}


//矩阵左乘向量
template <class T>
UpixVector<T> operator*(const UpixVector<T>& vec, const UpixMatrix<T>& mat)
{
	uint nCol = mat.GetColNum();
	uint nRow = mat.GetRowNum();
	assert(vec.Size() == nRow);
	UpixVector<T> result(nCol);
	for (uint j=0; j<nCol; j++)
	{
		for (uint i=0; i<nRow; i++)
		{
			result[j] += vec[i]*mat(i,j);
		}
	}
	return result;
}

//矩阵左乘数
template <class T>
UpixMatrix<T> operator*(const T &value, const UpixMatrix<T>& mat)
{
	return mat * value;
}

//矩阵左加数
template <class T>
UpixMatrix<T> operator+(const T &value, const UpixMatrix<T>& mat)
{
	return mat + value;
}

//矩阵左减数
template <class T>
UpixMatrix<T> operator-(const T &value, const UpixMatrix<T>& mat)
{
	uint nCol = mat.GetColNum();
	uint nRow = mat.GetRowNum();
	UpixMatrix<T> result(nRow, nCol);
	T *pMr = result.GetData();
	T *pMm = mat.GetData();
	uint n = nCol*nRow;
	for (uint i=0; i<n; i++)
	{
		pMr[i] = value - pMm[i];
	}
	return result;
}

//矩阵左除数
template <class T>
UpixMatrix<T> operator/(const T &value, const UpixMatrix<T> &mat)
{
	uint nCol = mat.GetColNum();
	uint nRow = mat.GetRowNum();
	UpixMatrix<T> result(nRow, nCol);
	T *pMr = result.GetData();
	T *pMm = mat.GetData();
	uint n = nCol*nRow;
	for (uint i=0; i<n; i++)
	{
		if (pMm == 0)
		{
			assert(false);
			return result;
		}
		pMr[i] = value / pMm[i];
	}
	return result;
}

template <class T>
void UpixMatrix<T>::NormalizeRow()
{
	UpixVector<T> vec(m_nCol);
	for (uint j=0; j<m_nRow; j++)
	{
		vec = GetRowVector(j);
		vec.Normalize();
		SetRow(j, vec);
	}
}

template <class T>
void UpixMatrix<T>::NormalizeCol()
{
	UpixVector<T> vec(m_nRow);
	for (uint j=0; j<m_nCol; j++)
	{
		vec = GetColVector(j);
		vec.Normalize();
		SetColumn(j, vec);
	}
}

template <class T>
double UpixMatrix<T>::Det() const
{
	uint i,j,k,is,js,l,u,v;
	double f,det,q,d;
	f = 1.0; 
	det = 1.0;

	if (m_nRow != m_nCol)
	{
		assert(false);
		return det;
	}

	uint n = m_nCol;

	uint buff_size = m_nCol*m_nRow;
	T* a = new T[buff_size];
	Fill(a, buff_size, m_pData);

	for (k=0; k<=n-2; k++)
	{ 
		q=0.0;
		for (i=k; i<=n-1; i++)
			for (j=k; j<=n-1; j++)
			{ 
				l=i*n+j;
				d=fabs(a[l]);
				if (d>q)
				{ 
					q=d;
					is=i;
					js=j;
				}
			}

			if (q+1.0==1.0)
			{ 
				det=0.0;
				delete [] a;
				return(det);
			}

			if (is!=k)
			{ 
				f=-f;
				for (j=k; j<=n-1; j++)
				{ 
					u=k*n+j;
					v=is*n+j;
					d=a[u];
					a[u]=a[v];
					a[v]=d;
				}
			}

			if (js!=k)
			{ 
				f=-f;
				for (i=k; i<=n-1; i++)
				{ 
					u=i*n+js;
					v=i*n+k;
					d=a[u];
					a[u]=a[v];
					a[v]=d;
				}
			}

			l=k*n+k;

			det=det*a[l];

			for (i=k+1; i<=n-1; i++)
			{ 
				d=a[i*n+k]/a[l];
				for (j=k+1; j<=n-1; j++)
				{ 
					u=i*n+j;
					a[u]=a[u]-d*a[k*n+j];
				}
			}
	}

	det=f*det*a[n*n-1];

	delete [] a;

	return(det);
}

template <class T>
int UpixMatrix<T>::Rank()
{
	uint i,j,k,nn,is,js,l,ll,u,v;
	double q,d;

	uint m = m_nRow;
	uint n = m_nCol;

	uint buff_size = m*n;
	T* a = new T[buff_size];
	Fill(a, buff_size, m_pData);

	nn = m;

	if (m >= n)
	{
		nn = n;
	}

	k=0;
	for (l=0; l<=nn-1; l++)
	{ 
		q=0.0;
		for (i=l; i<=m-1; i++)
		{
			for (j=l; j<=n-1; j++)
			{ 
				ll=i*n+j;
				d=fabs(a[ll]);
				if (d>q)
				{ 
					q=d; 
					is=i;
					js=j;
				}
			}
		}

		if (q+1.0==1.0)
		{
			return k;
		}

		k=k+1;
		if (is!=l)
		{ 
			for (j=l; j<=n-1; j++)
			{ 
				u=l*n+j;
				v=is*n+j;
				d=a[u];
				a[u]=a[v];
				a[v]=d;
			}
		}

		if (js!=l)
		{ 
			for (i=l; i<=m-1; i++)
			{ 
				u=i*n+js;
				v=i*n+l;
				d=a[u];
				a[u]=a[v];
				a[v]=d;
			}
		}

		ll=l*n+l;

		for (i=l+1; i<=n-1; i++)
		{ 
			d=a[i*n+l]/a[ll];
			for (j=l+1; j<=n-1; j++)
			{ 
				u=i*n+j;
				a[u]=a[u]-d*a[l*n+j];
			}
		}
	}

	delete [] a;

	return k;
}

template <class T>
int UpixMatrix<T>::SVD(UpixMatrix<T>& matU, UpixMatrix<T>& matS, UpixMatrix<T>& matVT, T eps) const
{
	uint i;
	uint m = m_nRow;
	uint n = m_nCol;

	if (matU.m_nRow != m || matU.m_nCol != m)
	{
		matU.Resize(m, m);
	}

	if (matS.m_nRow != m || matS.m_nCol != n)
	{
		matS.Resize(m, n);
	}

	if (matVT.m_nRow != n || matVT.m_nCol != n)
	{
		matVT.Resize(n, n);
	}

	uint n_size = m*m;
	T* u = new T[n_size];
	memset(u, 0, sizeof(T)*n_size);

	n_size = n*n;
	T* v = new T[n_size];
	memset(v, 0, sizeof(T)*n_size);

	n_size = m*n;
	T* a = new T[n_size];
	Fill(a, n_size, m_pData);

	uint ka = UpixMax(m, n) + 1;
	T *s = new T[ka];
	T *e = new T[ka];
	T *w = new T[ka];

	uint it = 60; 
	uint k = n;
	if (m-1 < n) k = m - 1;
	uint l = m;
	if (n-2 < m) l = n - 2;
	if (l < 0) l = 0;
	uint ll = k;
	if (l > k) ll = l;

	T d;
	uint kk, ix, iy, iz;
	if (ll >= 1)
	{
		for (kk=1; kk<=ll; kk++)
		{ 
			if (kk <= k)
			{
				d = 0.0f;
				for (i=kk; i<=m; i++)
				{
					ix = (i - 1)*n + kk - 1; 
					d = d + a[ix]*a[ix];
				}
				s[kk-1] = sqrt(d);

				if (!IsEqualZero(s[kk-1]))
				{
					ix = (kk - 1)*n + kk - 1;
					if (!IsEqualZero(a[ix]))
					{
						s[kk-1] = fabs(s[kk-1]);
						if (a[ix] < 0.0f)
						{
							s[kk-1] = - s[kk-1];
						}
					}
					for (i=kk; i<=m; i++)
					{
						iy = (i - 1)*n + kk - 1;
						a[iy] = a[iy]/s[kk-1];
					}
					a[ix] = 1.0f + a[ix];
				}
				s[kk-1] = -s[kk-1];
			} 
			if (n >= kk+1)
			{
				for (uint j=kk+1; j<=n; j++)
				{ 
					if ((kk<=k) && (!IsEqualZero(s[kk-1])))
					{ 
						d = 0.0f;
						for (i=kk; i<=m; i++)
						{
							ix = (i - 1)*n + kk - 1;
							iy = (i - 1)*n + j - 1;
							d = d + a[ix]*a[iy];
						}
						d = -d/a[(kk - 1)*n + kk - 1];
						for (uint i=kk; i<=m; i++)
						{ 
							ix = (i - 1)*n + j - 1;
							iy = (i - 1)*n + kk - 1;
							a[ix] += d*a[iy];
						}
					}
					e[j-1] = a[(kk - 1)*n + j - 1];
				}
			}
			if (kk <= k)
			{
				for (i=kk; i<=m; i++)
				{ 
					ix = (i - 1)*m + kk - 1;
					iy = (i - 1)*n + kk - 1;
					u[ix] = a[iy];
				}
			}
			if (kk <= l)
			{
				d = 0.0f;
				for (i=kk+1; i<=n; i++)
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
					for (i=kk+1; i<=n; i++)
					{
						e[i-1] = e[i-1]/e[kk-1];
					}
					e[kk] += 1.0f;
				}
				e[kk-1] = -e[kk-1];
				if ((kk+1<=m)&&(!IsEqualZero(e[kk-1])))
				{ 
					for (i=kk+1; i<=m; i++)
					{
						w[i-1] = 0.0f;
					}
					for (uint j=kk+1; j<=n; j++)
					{
						for (i=kk+1; i<=m; i++)
						{
							w[i-1] += e[j-1]*a[(i - 1)*n + j - 1];
						}
					}
					for (uint j=kk+1; j<=n; j++)
					{
						for (i=kk+1; i<=m; i++)
						{ 
							ix = (i - 1)*n + j - 1;
							a[ix] = a[ix] - w[i-1]*e[j-1]/e[kk];
						}
					}
				}
				for (i=kk+1; i<=n; i++)
				{
					v[(i - 1)*n + kk - 1] = e[i-1];
				}
			}
		}
	}

	uint mm = n;
	if (m+1 < n) 
	{
		mm = m + 1;
	}
	if (k < n)
	{
		s[k] = a[k*n + k];
	}
	if (m < mm)
	{
		s[mm-1] = 0.0f;
	}
	if (l+1 < mm)
	{
		e[l] = a[l*n + mm - 1];
	}
	e[mm-1] = 0.0f;
	uint nn = m;
	if (m > n) 
	{
		nn = n;
	}
	if (nn >= k+1)
	{ 
		for (uint j=k+1; j<=nn; j++)
		{ 
			for (i=1; i<=m; i++)
			{
				u[(i - 1)*m + j - 1] = 0.0f;
			}
			u[(j - 1)*m + j - 1] = 1.0f;
		}
	}
	if (k >= 1)
	{ 
		for (ll=1; ll<=k; ll++)
		{
			kk = k - ll + 1; 
			iz = (kk - 1)*m + kk - 1;
			if (!IsEqualZero(s[kk-1]))
			{
				if (nn >= kk+1)
				{
					for (uint j=kk+1; j<=nn; j++)
					{
						d = 0.0f;
						for (i=kk; i<=m; i++)
						{
							ix = (i - 1)*m + kk - 1;
							iy = (i - 1)*m + j - 1;
							d = d + u[ix]*u[iy]/u[iz];
						}
						d = -d;
						for (i=kk; i<=m; i++)
						{ 
							ix = (i - 1)*m + j - 1;
							iy = (i - 1)*m + kk - 1;
							u[ix] = u[ix] + d*u[iy];
						}
					}
				}
				for (i=kk; i<=m; i++)
				{
					ix = (i - 1)*m + kk - 1; 
					u[ix] = -u[ix];
				}
				u[iz] += 1.0f;
				if (kk-1 >= 1)
				{
					for (i=1; i<=kk-1; i++)
					{
						u[(i - 1)*m + kk - 1] = 0.0f;
					}
				}
			}
			else
			{ 
				for (i=1; i<=m; i++)
				{
					u[(i - 1)*m + kk - 1] = 0.0f;
				}
				u[(kk - 1)*m + kk - 1] = 1.0f;
			}
		}
	}
	for (ll=1; ll<=n; ll++)
	{ 
		kk = n - ll + 1; 
		iz = kk*n + kk - 1;
		if ((kk<=l) && (!IsEqualZero(e[kk-1])))
		{ 
			for (uint j=kk+1; j<=n; j++)
			{ 
				d = 0.0f;
				for (i=kk+1; i<=n; i++)
				{ 
					ix = (i - 1)*n + kk - 1; 
					iy = (i - 1)*n + j - 1;
					d += v[ix]*v[iy]/v[iz];
				}
				d = -d;
				for (i=kk+1; i<=n; i++)
				{
					ix = (i - 1)*n + j - 1;
					iy = (i - 1)*n + kk - 1;
					v[ix] += d*v[iy];
				}
			}
		}
		for (i=1; i<=n; i++)
		{
			v[(i - 1)*n + kk - 1] = 0.0f;
		}
		v[iz - n] = 1.0f;
	}
	
	memset(a, 0, sizeof(T)*m*n);

	uint m1 = mm, ks; 
	it = 60;
	T dd, t, sm, sm1, em1, sk, ek, b, c, shh;
	T fg[2], cs[2];
	while (1)
	{
		if (mm == 0)
		{
			ppp(a, e, s, v, m, n);
			matS.SetData(m, n, a);
			matU.SetData(m, m, u);
			matVT.SetData(n, n, v);

			delete [] s;
			delete [] e;
			delete [] w;
			delete [] a;
			delete [] u;
			delete [] v;
			return 1;
		}
		if (it == 0)
		{ 
			ppp(a, e, s, v, m, n);
			matS.SetData(m, n, a);
			matU.SetData(m, m, u);
			matVT.SetData(n, n, v);

			delete [] s;
			delete [] e;
			delete [] w;
			delete [] a;
			delete [] u;
			delete [] v;

			assert(false);//工作失败

			return -1;
		}
		kk = mm - 1;
		while ((kk!=0) && (!IsEqualZero(e[kk-1])))
		{
			d = fabs(s[kk-1]) + fabs(s[kk]);
			dd = fabs(e[kk-1]);
			if (dd > eps*d)
			{
				kk = kk - 1;
			}
			else
			{
				e[kk-1] = 0.0f;
			}
		}
		if (kk == mm-1)
		{
			kk = kk + 1;
			if (s[kk-1] < 0.0f)
			{
				s[kk-1] = -s[kk-1];
				for (i=1; i<=n; i++)
				{ 
					ix = (i - 1)*n + kk - 1;
					v[ix] = -v[ix];
				}
			}
			while ((kk!=m1) && (s[kk-1]<s[kk]))
			{
				d = s[kk-1]; s[kk-1] = s[kk]; s[kk] = d;
				if (kk < n)
				{
					for (i=1; i<=n; i++)
					{ 
						ix = (i - 1)*n + kk - 1; 
						iy = (i - 1)*n + kk;
						d = v[ix]; v[ix] = v[iy]; v[iy] = d;
					}
				}

				if (kk < m)
				{
					for (i=1; i<=m; i++)
					{ 
						ix = (i - 1)*m + kk - 1; 
						iy = (i - 1)*m + kk;
						d = u[ix]; u[ix] = u[iy]; u[iy] = d;
					}
				}
				kk = kk + 1;
			}
			it = 60;
			mm = mm - 1;
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
					ks = ks-1;
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
						for (uint j=1; j<=n; j++)
						{
							ix = (j - 1)*n + i - 1;
							iy = (j - 1)*n + i;
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

					if (i < m)
					{
						if ((!IsEqualZero(cs[0]-1.0f)) || (!IsEqualZero(cs[1])))
						{
							for (uint j=1; j<=m; j++)
							{
								ix = (j - 1)*m + i - 1;
								iy = (j - 1)*m + i;
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
							for (uint j=1; j<=n; j++)
							{
								ix = (j - 1)*n + i - 1;
								iy = (j - 1)*n + mm - 1;
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
							for (uint j=1; j<=m; j++)
							{ 
								ix = (j - 1)*m + i - 1;
								iy = (j - 1)*m + kk - 2;
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

template <class T>
void UpixMatrix<T>::ppp(T a[], T e[], T s[], T v[], uint m, uint n) const
{
	uint i, j, p, q;
	T d;
	if (m >= n) i = n;
	else i = m;
	for (j=1; j<=i-1; j++)
	{ 
		a[(j-1)*n+j-1] = s[j-1];
		a[(j-1)*n+j]   = e[j-1];
	}
	a[(i-1)*n+i-1]=s[i-1];

	if (m<n) a[(i-1)*n+i]=e[i-1];

	for (i=1; i<=n-1; i++)
		for (j=i+1; j<=n; j++)
		{
			p=(i-1)*n+j-1; 
			q=(j-1)*n+i-1;
			d=v[p]; 
			v[p]=v[q];
			v[q]=d;
		}
}

template <class T>
void UpixMatrix<T>::sss(T fg[2], T cs[2]) const
{
	T r,d;
	T temp1, temp2;
	temp1 = fabs(fg[0]);
	temp2 = fabs(fg[1]);
	if (IsEqualZero(temp1 + temp2))
	{
		cs[0]=1.0f;
		cs[1]=0.0f;
		d=0.0f;
	}
	else 
	{ 
		d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);

		if (temp1 > temp2)
		{ 
			d=fabs(d);

			if (fg[0]<0.0f)
			{
				d=-d;
			}
		}
		else
		{
			d=fabs(d);
			if (fg[1]<0.0f)
			{
				d=-d;
			}
		}

		cs[0]=fg[0]/d; cs[1]=fg[1]/d;
	}
	r=1.0f;
	if (temp1 > temp2) 
	{
		r=cs[1];
	}
	else
	{
		if (!IsEqualZero(cs[0]))
		{
			r=1.0f/cs[0];
		}
	}
	fg[0]=d; fg[1]=r;
}

template <class T>
inline bool UpixMatrix<T>::IsEqualZero(T _num) const
{
	double eps = 1e-16f;
	return ( (_num > - eps) && (_num < eps) );
}

template <class T>
void UpixMatrix<T>::Gradient(UpixMatrix<double> &Fx, UpixMatrix<double> &Fy)
{
	if (Fx.GetColNum() != m_nCol || Fx.GetRowNum() != m_nRow)
	{
		Fx.Resize(m_nRow, m_nCol);
	}
	if (Fy.GetColNum() != m_nCol || Fy.GetRowNum() != m_nRow)
	{
		Fy.Resize(m_nRow, m_nCol);
	}

	for (uint j=0; j<m_nRow; j++)
	{
		Fx(j,0) = (double)(m_pData[j*m_nCol+1] - m_pData[j*m_nCol]);
		for (uint i=1; i<m_nCol-1; i++)
		{
			Fx(j,i) = (m_pData[j*m_nCol+(i+1)] - m_pData[j*m_nCol+(i-1)])*0.5;
		}
		Fx(j,m_nCol-1) = (double)(m_pData[j*m_nCol+(m_nCol-1)] - m_pData[j*m_nCol+(m_nCol-2)]);
	}

	for (uint i=0; i<m_nCol; i++)
	{
		Fy(0,i) = (double)(m_pData[m_nCol+i] - m_pData[i]);
		for (uint j=1; j<m_nRow-1; j++)
		{
			Fy(j,i) = (m_pData[(j+1)*m_nCol+i] - m_pData[(j-1)*m_nCol+i])*0.5;
		}
		Fy(m_nRow-1,i) = (double)(m_pData[(m_nRow-1)*m_nCol+i] - m_pData[(m_nRow-2)*m_nCol+i]);
	}
}

template <class T>
bool UpixMatrix<T>::Swap(UpixMatrix<T> &M)
{
	uint tmpCol = this->m_nCol;
	uint tmpRow = this->m_nRow;
	T*tmpPdata = this->m_pData;

	this->m_nCol = M.m_nCol;
	this->m_nRow = M.m_nRow;
	this->m_pData = M.m_pData;

	M.m_nCol = tmpCol;
	M.m_nRow = tmpRow;
	M.m_pData = tmpPdata;

	return true;
}

template <class T>
void UPIX_ALGORITHM_TEMPLATE MultVectorAsMatrix(const UpixVector<T> &v1, const UpixVector<T> &v2, UpixMatrix<T> &ret)
{
	uint nRows = v1.Size();
	uint nCols = v2.Size();
	assert(ret.GetRowNum() == nRows && ret.GetColNum() == nCols);

	T *pV1 = v1.GetData();
	T *pV2 = v2.GetData();
	T *pM = ret.GetData();

	for (uint i=0; i<nRows; i++)
	{
		for (uint j=0; j<nCols; j++)
		{
			*pM = (*pV1)*pV2[j];
			pM++;
		}
		pV1++;
	}
}

template <class T>
void UPIX_ALGORITHM_TEMPLATE MultVector(const UpixMatrix<T> &mat, const UpixVector<T> &vec, UpixVector<T> &ret)
{
	uint nRows = mat.GetRowNum();
	uint nCols = mat.GetColNum();
	assert(vec.Size() == nCols);

	T *pM = mat.GetData();
	T *pV = vec.GetData();
	T *p = ret.GetData();

	for (uint i=0; i<nRows; i++)
	{
		*p = 0;
		for (uint j=0; j<nCols; j++)
		{
			*p += (*pM) * pV[j];
			pM++;
		}
		p++;
	}
}

template <class T>
void UPIX_ALGORITHM_TEMPLATE MultMatrix(const UpixMatrix<T> &m1, const UpixMatrix<T> &m2, UpixMatrix<T> &ret)
{
	uint nRow = m1.GetRowNum();
	uint nCol = m2.GetColNum();
	uint nSharedRow = m1.GetColNum();

	assert(nSharedRow==m2.GetRowNum() && ret.GetRowNum()==nRow && ret.GetColNum()==nCol);

	T *pRes = ret.GetData();
	T *pRow = m1.GetData();
	T *pM2 = m2.GetData();
	T *pData = NULL;
	T *pCol = NULL;
	for(uint i=0; i<nRow; ++i)
	{
		for(uint j=0; j<nCol; ++j)
		{
			*pRes = 0;
			pData = pRow;
			pCol = pM2 + j;
			for(uint k=0; k<nSharedRow; ++k)
			{
				*pRes += (*pData++) * (*pCol);
				pCol += nCol;
			}
			pRes++;
		}
		pRow += nSharedRow;
	}
}



#endif