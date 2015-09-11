/*
* Filename: UpixRotation.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_ROTATION_H_
#define _UPIX_ROTATION_H_

#include "UpixQuat.h"



template <class T>
class  UPIX_ALGORITHM_TEMPLATE UpixRotation
{
public:
	// Constructors:-------------------------------------

	//: Construct the identity rotation
	UpixRotation() : q_(0,0,0,1) {}

	//: Construct from a quaternion.
	UpixRotation( const UpixQuat<T> &q ) : q_(q) { q_.Normalize(); }

	//: Construct from Euler angles.
	UpixRotation( const T &rx, const T &ry, const T &rz ) : q_(rx, ry, rz) {}

	//: Construct from a Rodrigues vector.
	UpixRotation( const UpixVector<T> &rvector )
	{
		T vec_norm = rvector.TwoNorm();
		if (vec_norm == 0)
		{
			q_ = UpixQuat<T>(rvector, 0);
		}
		else
		{
			q_ = UpixQuat<T>(rvector/vec_norm, vec_norm);
		}
	}

	UpixRotation( const UpixPoint3D<T> &rvector )
	{
		T vec_norm = rvector.TwoNorm();
		if (vec_norm == 0)
		{
			q_.MakeRot(rvector, 0);
		}
		else
		{
			q_.MakeRot(rvector/vec_norm, vec_norm);
		}
	}

	//: Construct from a 3x3 rotation matrix
	UpixRotation( const UpixMatrix<T> &matrix ) : q_(matrix) {}

	UpixRotation(const UpixMatrix3x3<T> &matrix) : q_(matrix) {}


	// Conversions:--------------------------------------

	//: Output unit quaternion.
	UpixQuat<T> AsQuaternion() const
	{
		return q_;
	}

	//: Output Euler angles.
	//  The first element is the rotation about the x-axis
	//  The second element is the rotation about the y-axis
	//  The third element is the rotation about the z-axis
	//  The total rotation is a composition of these rotations in this order
	UpixPoint3D<T> AsEulerAngles() const
	{
		UpixPoint3D<T> vAngles;
		q_.GetEulerAngles(vAngles.x, vAngles.y, vAngles.z);
		return vAngles;
	}

	void AsEulerAngles(UpixVector<T> &vAngles) const
	{
		T rx, ry, rz;
		q_.GetEulerAngles(rx, ry, rz);
		vAngles[0] = rx;
		vAngles[1] = ry;
		vAngles[2] = rz;
	}

	void AsEulerAngles(UpixPoint3D<T> &vAngles) const
	{
		q_.GetEulerAngles(vAngles.x, vAngles.y, vAngles.z);
	}

	//: Output Rodrigues vector.
	//  The direction of this vector is the axis of rotation
	//  The length of this vector is the angle of rotation in radians
	UpixPoint3D<T> AsRodrigues() const
	{
		UpixPoint3D<T> vec;
		T angle;
		q_.GetRotate(vec, angle);
		return vec*angle;
	}

	void AsRodrigues(UpixVector<T> &vec) const
	{
		T angle;
		q_.GetRotate(vec, angle);
		vec *= angle;
	}

	void AsRodrigues(UpixPoint3D<T> &vec) const
	{
		T angle;
		q_.GetRotate(vec, angle);
		vec *= angle;
	}

	//: Output the matrix representation of this rotation in 3x3 form.
	UpixMatrix3x3<T> AsMatrix() const
	{
		UpixMatrix3x3<T> m;
		q_.GetMatrix(m);
		return m;
	}

	void AsMatrix(UpixMatrix<T> &m) const
	{
		q_.GetMatrix(m);
	}

	void AsMatrix(UpixMatrix3x3<T> &m) const
	{
		q_.GetMatrix(m);
	}

	//: Returns the axis of rotation (unit vector) an the magnitude of the angle of rotation
	void AsAxisAngle(UpixVector<T> &vec, T &angle) const
	{
		q_.GetRotate(vec, angle);
	}

	void AsAxisAngle(UpixPoint3D<T> &vec, T &angle) const
	{
		q_.GetRotate(vec, angle);
	}

	// Operations:----------------------------------------

	//: Make the rotation the identity (i.e. no rotation)
	void SetIdentity() { q_.Set(0, 0, 0, 1); }

	//: The inverse rotation
	UpixRotation<T> Inverse() const { return UpixRotation<T>(q_.Conj()); }

	//: Composition of two rotations.
	UpixRotation<T> operator*( const UpixRotation<T>& first_rotation ) const
	{
		return UpixRotation<T>( q_*first_rotation.q_ );
	}

	//: Rotate a homogeneous point UpixVector(3).
	UpixVector<T> operator*( const UpixVector<T>& v ) const
	{
		assert(v.Size()==3);
		UpixVector<T> res(3);
		Rotate(v, res);
		return res;
	}

	UpixPoint3D<T> operator*( const UpixPoint3D<T>& v ) const
	{
		UpixPoint3D<T> res;
		Rotate(v, res);
		return res;
	}

	void Rotate( const UpixVector<T> &v, UpixVector<T> &res )
	{
		assert(v.Size() == 3);
		UpixMatrix3x3<T> R;
		AsMatrix(R);
		if (res.Size() != 3)
		{
			res.Resize(3);
		}

		res[0] = R[0][0]*v[0] + R[0][1]*v[1] + R[0][2]*v[2];
		res[1] = R[1][0]*v[0] + R[1][1]*v[1] + R[1][2]*v[2];
		res[2] = R[2][0]*v[0] + R[2][1]*v[1] + R[2][2]*v[2];
	}

	void Rotate( const UpixPoint3D<T> &v, UpixPoint3D<T> &res )
	{
		UpixMatrix3x3<T> R;
		AsMatrix(R);
		
		res.x = R[0][0]*v.x + R[0][1]*v.y + R[0][2]*v.z;
		res.y = R[1][0]*v.x + R[1][1]*v.y + R[1][2]*v.z;
		res.z = R[2][0]*v.x + R[2][1]*v.y + R[2][2]*v.z;
	}

protected:
	//: The internal representation of the rotation is a quaternion.
	UpixQuat<T> q_;
};



// External methods for constructing rotation matrix
// ----------------------------------------------------------------
template <class T>
void RotateX(UpixMatrix<T> &Rx, const T &fPitch)
{
	T cosPitch = cos(fPitch);
	T sinPitch = sin(fPitch);

	Rx.SetIdentity(3);
	Rx[1][1] = cosPitch;
	Rx[1][2] = - sinPitch;
	Rx[2][1] = sinPitch;
	Rx[2][2] = cosPitch;
}

template <class T>
void RotateX(UpixMatrix3x3<T> &Rx, const T &fPitch)
{
	T cosPitch = cos(fPitch);
	T sinPitch = sin(fPitch);

	Rx.SetIdentity(3);
	Rx[1][1] = cosPitch;
	Rx[1][2] = - sinPitch;
	Rx[2][1] = sinPitch;
	Rx[2][2] = cosPitch;
}

template <class T>
void RotateY(UpixMatrix<T> &Ry, const T &fPan)
{
	T cosPan = cos(fPan);
	T sinPan = sin(fPan);

	Ry.SetIdentity(3);
	Ry[0][0]  = cosPan;
	Ry[0][2]  = sinPan;
	Ry[2][0]  = -sinPan;
	Ry[2][2]  = cosPan;
}

template <class T>
void RotateY(UpixMatrix3x3<T> &Ry, const T &fPan)
{
	T cosPan = cos(fPan);
	T sinPan = sin(fPan);

	Ry.SetIdentity(3);
	Ry[0][0]  = cosPan;
	Ry[0][2]  = sinPan;
	Ry[2][0]  = -sinPan;
	Ry[2][2]  = cosPan;
}

template <class T>
void RotateZ(UpixMatrix<T> &Rz, const T &fRoll)
{
	T cosRoll = cos(fRoll);
	T sinRoll = sin(fRoll);

	Rz.SetIdentity(3);
	Rz[0][0] = cosRoll;
	Rz[0][1] = -sinRoll;
	Rz[1][0] = sinRoll;
	Rz[1][1] = cosRoll;
}

template <class T>
void RotateZ(UpixMatrix3x3<T> &Rz, const T &fRoll)
{
	T cosRoll = cos(fRoll);
	T sinRoll = sin(fRoll);

	Rz.SetIdentity(3);
	Rz[0][0] = cosRoll;
	Rz[0][1] = -sinRoll;
	Rz[1][0] = sinRoll;
	Rz[1][1] = cosRoll;
}

template <class T>
void MatrixRotateX4(UpixMatrix<T> &Rx, const T &fPitch)
{
	T cosPitch = cos(fPitch);
	T sinPitch = sin(fPitch);

	Rx.SetIdentity(4);
	Rx[1][1] = cosPitch;
	Rx[1][2] = - sinPitch;
	Rx[2][1] = sinPitch;
	Rx[2][2] = cosPitch;
}

template <class T>
void MatrixRotateY4(UpixMatrix<T> &Ry, const T &fPan)
{
	T cosPan = cos(fPan);
	T sinPan = sin(fPan);

	Ry.SetIdentity(4);
	Ry[0][0]  = cosPan;
	Ry[0][2]  = sinPan;
	Ry[2][0]  = -sinPan;
	Ry[2][2]  = cosPan;
}

template <class T>
void MatrixRotateZ4(UpixMatrix<T> &Rz, const T &fRoll)
{
	T cosRoll = cos(fRoll);
	T sinRoll = sin(fRoll);

	Rz.SetIdentity(4);
	Rz[0][0] = cosRoll;
	Rz[0][1] = -sinRoll;
	Rz[1][0] = sinRoll;
	Rz[1][1] = cosRoll;
}

#endif // _rotation_h_
