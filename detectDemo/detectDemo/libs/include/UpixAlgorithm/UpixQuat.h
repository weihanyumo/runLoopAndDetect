/*
* Filename: UpixQuat.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_QUAT_H_
#define _UPIX_QUAT_H_

#include <math.h>
#include "UpixMatrix.h"
#include "UpixMatrix3x3.h"


#define UPIX_QUAT_TOLERANCE	1e-6



template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixQuat
{
public:
	//默认构造，元素置(0,0,0,1)
	UpixQuat() : x(0), y(0), z(0), w(1) {}

	//元素赋值构造
	UpixQuat(const T &_x, const T &_y, const T &_z, const T &_w) : x(_x), y(_y), z(_z), w(_w) {}

	//转轴和转角构造
	UpixQuat(const UpixVector<T> &vec, const T &rad)
	{
		assert(vec.Size() == 3);
		MakeRot(UpixPoint3D<T>(vec[0], vec[1], vec[2]), rad);
	}

	UpixQuat(const UpixPoint3D<T> &vec, const T &rad)
	{
		MakeRot(vec, rad);
	}

	//旋转矩阵构造
	UpixQuat(const UpixMatrix<T> &rot)
	{
		assert(rot.GetRowNum()==3 && rot.GetColNum()==3);

		UpixMatrix3x3<T> mat(rot.GetData());
		MakeRot(mat);
	}

	UpixQuat(const UpixMatrix3x3<T> &rot)
	{
		MakeRot(rot);
	}

	//三方向旋转轴和转角构造
	UpixQuat(const UpixVector<T> &vec_1, const T &rad_1, 
			const UpixVector<T> &vec_2, const T &rad_2, 
			const UpixVector<T> &vec_3, const T &rad_3)
	{
		UpixQuat<T> q1(vec_1, rad_1);
		UpixQuat<T> q2(vec_2, rad_2);
		UpixQuat<T> q3(vec_3, rad_3);
		*this = q1*q2*q3;
	}

	//欧拉角构造
	UpixQuat(const T &theta_X, const T &theta_Y, const T &theta_Z)
	{
		UpixQuat<T> Rx(sin(theta_X/2), 0, 0, cos(theta_X/2));
		UpixQuat<T> Ry(0, sin(theta_Y/2), 0, cos(theta_Y/2));
		UpixQuat<T> Rz(0, 0, sin(theta_Z/2), cos(theta_Z/2));
		*this = Rx * Ry * Rz;
	}

	~UpixQuat() {}

	//Copy constructor 四元数构造
	UpixQuat(const UpixQuat<T> &q) : x(q.x), y(q.y), z(q.z), w(q.w)  {}

	//重载赋值 = 
	UpixQuat<T> &operator=(const UpixQuat<T> &q)
	{
		if (this != &q)
		{
			x = q.x;
			y = q.y;
			z = q.z;
			w = q.w;
		}
		return *this;
	}

public:
	//按指定轴和转角生成四元数
	void MakeRot(const UpixPoint3D<T> &vec, const T &rad)
	{
		if (rad == 0)
		{
			x = 0;
			y = 0;
			z = 0;
			w = 1;
		}
		else
		{
			T half(0.5);
			T sinHalfAngleDivideNorm = sin(half*rad)/vec.TwoNorm();
			x = vec.x * sinHalfAngleDivideNorm;
			y = vec.y * sinHalfAngleDivideNorm;
			z = vec.z * sinHalfAngleDivideNorm;
			w = cos(half*rad);
		}
	}

	void MakeRot(const UpixMatrix3x3<T> &rot)
	{
		UpixMatrix3x3<T> mat(rot);
		mat.NormalizeCol();	//set rot matrix normalize

		T s;
		T tr = mat.Trace();
		T one(1);
		T half(0.5);

		// check the diagonal
		if (tr > 0.0)
		{
			s = sqrt(tr + one);
			w = s * half;
			s = half / s;
			x = (mat[2][1] - mat[1][2]) * s;
			y = (mat[0][2] - mat[2][0]) * s;
			z = (mat[1][0] - mat[0][1]) * s;
		}
		else	// diagonal is negative
		{
			T tq[3];
			static const int nxt[3] = {1, 2, 0};

			int i = 0;
			if (mat[1][1] > mat[0][0])
			{
				i = 1;
			}

			if (mat[2][2] > mat[i][i])
			{
				i = 2;
			}

			int j = nxt[i];
			int k = nxt[j];
			
			s = sqrt((mat[i][i] - (mat[j][j] + mat[k][k])) + one);

			tq[i] = s * half;

			if (s != 0.0)
			{
				s = half / s;
			}

			tq[j] = (mat[j][i] + mat[i][j]) * s;
			tq[k] = (mat[k][i] + mat[i][k]) * s;

			w = (mat[k][j] - mat[j][k]) * s;
			x = tq[0];
			y = tq[1];
			z = tq[2];
		}
	}

	//赋值
	void Set(const T &_x, const T &_y, const T &_z, const T &_w)
	{
		x = _x;
		y = _y;
		z = _z;
		w = _w;
	}

	//重载恒等运算符 ==，若两四元数所有元素值相同返回true
	bool operator==(const UpixQuat<T> &other) const
	{
		return ((UpixAbs(x - other.x) < UPIX_QUAT_TOLERANCE) &&
			(UpixAbs(y - other.y) < UPIX_QUAT_TOLERANCE) &&
			(UpixAbs(z - other.z) < UPIX_QUAT_TOLERANCE) &&
			(UpixAbs(w - other.w) < UPIX_QUAT_TOLERANCE));
	}

	//重载运算符 !=
	bool operator!=(const UpixQuat<T>& other) const
	{
		return !(*this == other);
	}

	//重载运算符 + 右加四元数
	UpixQuat<T> operator+(const UpixQuat<T>& other) const
	{
		return UpixQuat<T>(x + other.x, y + other.y, z + other.z, w + other.w);
	}

	//重载运算符 - 右减四元数
	UpixQuat<T> operator-(const UpixQuat<T>& other) const
	{
		return UpixQuat<T>(x - other.x, y - other.y, z - other.z, w - other.w);
	}

	//重载运算符 * 右乘四元数
	UpixQuat<T> operator*(const UpixQuat<T>& other) const
	{
		UpixQuat<T> result(y*other.z - z*other.y + w*other.x + x*other.w,
						z*other.x - x*other.z + w*other.y + y*other.w,
						x*other.y - y*other.x + w*other.z + z*other.w,
						w*other.w - x*other.x - y*other.y - z*other.z);

		result.Normalize();
		return result;
	}

	//重载运算符 / 右除四元数
	UpixQuat<T> operator/(const UpixQuat<T>& other) const
	{
		return (*this) * (other.Inverse());
	}

	//重载运算符 += 加等四元数
	UpixQuat<T>& operator+=(const UpixQuat<T>& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		w += other.w;
		return *this;
	}
 
	//重载运算符 -= 减等四元数
	UpixQuat<T>& operator-=(const UpixQuat<T>& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		w -= other.w;
		return *this;
	}

	//重载运算符 *= 乘等四元数
	UpixQuat<T>& operator*=(const UpixQuat<T>& other)
	{
		UpixQuat<T> result(y*other.z - z*other.y + w*other.x + x*other.w,
						z*other.x - x*other.z + w*other.y + y*other.w,
						x*other.y - y*other.x + w*other.z + z*other.w,
						w*other.w - x*other.x - y*other.y - z*other.z);
		result.Normalize();
		x = result.x;
		y = result.y;
		z = result.z;
		w = result.w;
		return *this;
	}

	//重载运算符 /= 除等四元数
	UpixQuat<T>& operator/=(const UpixQuat<T>& other)
	{
		*this *= other.Inverse();
		return *this;
	}

	//重载运算符 + 右加数
	UpixQuat<T> operator+(const T &value) const
	{
		return UpixQuat<T>(x + value, y + value, z + value, w + value);
	}

	//重载运算符 - 右减数
	UpixQuat<T> operator-(const T &value) const
	{
		return UpixQuat<T>(x - value, y - value, z - value, w - value);
	}

	//重载运算符 * 右乘数
	UpixQuat<T> operator*(const T &value) const
	{
		return UpixQuat<T>(x*value, y*value, z*value, w*value);
	}

	//重载运算符 / 右除数
	UpixQuat<T> operator/(const T &value) const
	{
		assert(value != 0);
		T one(1);
		T invd = one/value;
		return UpixQuat<T>(x*invd, y*invd, z*invd, w*invd);
	}

	//重载运算符 += 加等数
	UpixQuat<T>& operator+=(const T &value)
	{
		x += value;
		y += value;
		z += value;
		w += value;
		return *this;
	}

	//重载运算符 -= 减等数
	UpixQuat<T>& operator-=(const T &value)
	{
		x -= value;
		y -= value;
		z -= value;
		w -= value;
		return *this;
	}

	//重载运算符 *= 乘等数
	UpixQuat<T>& operator*=(const T &value)
	{
		x *= value;
		y *= value;
		z *= value;
		w *= value;
		return *this;
	}

	//重载运算符 /= 除等数
	UpixQuat<T>& operator/=(const T &value)
	{
		assert(value != 0);
		T one(1);
		T invd = one/value;
		x *= invd;
		y *= invd;
		z *= invd;
		w *= invd;
		return *this; 
	}

	//所有元素取相反数
	UpixQuat<T> &operator-()
	{
		x = -x;
		y = -y;
		z = -z;
		w = -w;
		return *this;
	}

	UpixQuat<T> operator-() const
	{
		return UpixQuat<T>(-x, -y, -z, -w);
	}

	// Rotate a vector
	UpixPoint3D<T> operator*(const UpixPoint3D<T> &v) const
	{
		T Two(2);

		UpixPoint3D<T> uv, uuv;
		UpixPoint3D<T> QuatVector(x, y, z);
		uv = Cross(QuatVector, v);
		uuv = Cross(QuatVector, uv);
		uv *= (Two * w); 
		uuv *= Two; 

		return v + uv + uuv;
	}

	//正则化
	void Normalize() 
	{
		T one(1);
		T invd = one/TwoNrom();
		x *= invd;
		y *= invd;
		z *= invd;
		w *= invd;
	}

	//一阶模
	T OneNorm() const
	{ 
		return UpixAbs(x) + UpixAbs(y) + UpixAbs(z) + UpixAbs(w); 
	}

	//二阶模
	T TwoNrom() const
	{
		return sqrt(x*x + y*y + z*z + w*w);
	}

	//二阶模的平方
	T TwoNrom2() const
	{
		return x*x + y*y + z*z + w*w;
	}

	//共轭
	UpixQuat<T> Conj() const
	{
		return UpixQuat<T>(-x, -y, -z, w);
	}

	//乘法逆
	UpixQuat<T> Inverse() const
	{
		return Conj()/TwoNrom2();
	}

	//获取欧拉角
	void GetEulerAngles(T &theta_X, T &theta_Y, T &theta_Z)
	{
		Normalize();
		T two(2);

		T x2 = x*x;
		T y2 = y*y;
		T z2 = z*z;
		T w2 = w*w;

		T sin_theta = 2 * (w*y - x*z);
		if (UpixAbs(sin_theta) > 0.99999)
		{
			theta_X = 0;
			theta_Y = asin(sin_theta);
			theta_Z = atan2(two*(w*z - x*y), w2 - x2 + y2 - z2);
		}
		else
		{
			theta_X = atan2(two*(y*z + w*x), w2 - x2 - y2 + z2);
			theta_Y = asin(sin_theta);
			theta_Z = atan2(two*(x*y + w*z), w2 + x2 - y2 - z2);
		}
	}

	void GetEulerAngles(T &theta_X, T &theta_Y, T &theta_Z) const
	{
		UpixQuat<T> nq(x, y, z, w);
		nq.GetEulerAngles(theta_X, theta_Y, theta_Z);
	}

	//获取转轴和转角
	bool GetRotate(UpixVector<T> &vec, T &angle) const
	{
		if (vec.Size()!=3)
		{
			vec.Resize(3);
		}

		UpixPoint3D<T> pt;
		GetRotate(pt, angle);
		vec[0] = pt.x;
		vec[1] = pt.y;
		vec[2] = pt.z;
	}

	bool GetRotate(UpixPoint3D<T> &vec, T &angle) const
	{
		T sinhalfangle = sqrt(x*x + y*y + z*z);
		angle = 2 * atan2(sinhalfangle, w);

		if(UpixAbs(sinhalfangle) > UPIX_QUAT_TOLERANCE)
		{
			T one(1);
			T invd = one/sinhalfangle;
			vec.x = x * invd;
			vec.y = y * invd;
			vec.z = z * invd;
		} 
		else 
		{
			vec.x = 0;
			vec.y = 0;
			vec.z = 0;
		}
		return true;
	}

	//由旋转矩阵生成四元数
	void SetMatrix(const UpixMatrix<T> &m)
	{
		UpixMatrix3x3<T> mat(m.GetData());
		MakeRot(mat);
	}

	void SetMatrix(const UpixMatrix3x3<T> &m)
	{
		MakeRot(m);
	}

	//获取旋转矩阵R(3*3)
	void GetMatrix(UpixMatrix<T> &m) const
	{
		UpixMatrix3x3<T> rot;
		GetMatrix(rot);
		m.SetData(3, 3, rot.GetData());
	}

	void GetMatrix(UpixMatrix3x3<T> &m)
	{
		/* unit quaternion to rotation matrix
		| 1-2*y*y-2*z*z    2*x*y-2*w*z    2*x*z+2*w*y   |
		| 2*x*y+2*w*z      1-2*x*x-2*z*z  2*y*z-2*w*x   |
		| 2*x*z-2*w*y      2*y*z+2*w*x    1-2*x*x-2*y*y |
		*/

		Normalize();

		T wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;

		x2 = x + x;
		y2 = y + y;
		z2 = z + z;

		xx = x * x2;
		xy = x * y2;
		xz = x * z2;

		yy = y * y2;
		yz = y * z2;
		zz = z * z2;

		wx = w * x2;
		wy = w * y2;
		wz = w * z2;

		m[0][0] = 1 - (yy + zz);
		m[0][1] = xy - wz;
		m[0][2] = xz + wy;

		m[1][0] = xy + wz;
		m[1][1] = 1 - (xx + zz);
		m[1][2] = yz - wx;

		m[2][0] = xz - wy;
		m[2][1] = yz + wx;
		m[2][2] = 1 - (xx + yy);
	}

	void GetMatrix(UpixMatrix3x3<T> &m) const
	{
		UpixQuat<T> nq(x, y, z, w);
		nq.GetMatrix(m);
	}

public:
	T x;    //数据缓冲区
	T y;
	T z;
	T w;
};


typedef UpixQuat<float> UpixQuatf;
typedef UpixQuat<double> UpixQuatd;



// Left multiply a value
template <class T>
UpixQuat<T> operator*(double value, const UpixQuat<T> &q)
{
	return q*value;
}


// Dot multiply
template <class T> 
T Dot(const UpixQuat<T> &q1, const UpixQuat<T> &q2)
{
	return q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;
}



// Interpolate the quaternion
template <class T>
UpixQuat<T> Interpolate(const UpixQuat<T> &x, const UpixQuat<T> &y, const T &percent)
{
	T one(1);
	T zero(0);

	if(percent <= zero) return x;
	if(percent >= one) return y;

	float fCos = Dot(x, y);
	UpixQuat<T> y2(y);
	if(fCos < zero)
	{
		y2 = -y;
		fCos = -fCos;
	}

	//if(fCos > 1.0f) // problem
	float k0, k1;
	if(fCos > T(0.9999))
	{
		k0 = one - percent;
		k1 = zero + percent; //BUG!!! 1.0f + a;
	}
	else
	{
		T fSin = sqrt(one - fCos * fCos);
		T fAngle = atan2(fSin, fCos);
		T fOneOverSin = one / fSin;
		k0 = sin((one - percent) * fAngle) * fOneOverSin;
		k1 = sin((zero + percent) * fAngle) * fOneOverSin;
	}

	return UpixQuat<T>(
					k0 * x.x + k1 * y2.x,
					k0 * x.y + k1 * y2.y,
					k0 * x.z + k1 * y2.z,
					k0 * x.w + k1 * y2.w);
}



#endif