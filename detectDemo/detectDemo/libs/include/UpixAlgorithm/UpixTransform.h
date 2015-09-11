/*
* Filename: UpixTransform.h
* Author:   Xiehong
* Date:     2012.12.28
*/

#ifndef _UPIX_TRANSFORM_H_
#define _UPIX_TRANSFORM_H_


#include "UpixMatrix.h"
#include "UpixMatrix3x3.h"
#include "UpixPoint.h"



template <class Type>
void UPIX_ALGORITHM_TEMPLATE Homography(const UpixMatrix<Type> &H, const UpixPoint2D<Type> &p, UpixPoint2D<Type> &xp)
{
	Type *pH = H.GetData();
	Type invw = 1/(pH[6]*p.x + pH[7]*p.y + pH[8]);
	xp.x = (pH[0]*p.x + pH[1]*p.y + pH[2]) * invw;
	xp.y = (pH[3]*p.x + pH[4]*p.y + pH[5]) * invw;
}

template <class Type>
void UPIX_ALGORITHM_TEMPLATE Homography(const UpixMatrix3x3<Type> &H, const UpixPoint2D<Type> &p, UpixPoint2D<Type> &xp)
{
	Type *pH = H.GetData();
	Type invw = 1/(pH[6]*p.x + pH[7]*p.y + pH[8]);
	xp.x = (pH[0]*p.x + pH[1]*p.y + pH[2]) * invw;
	xp.y = (pH[3]*p.x + pH[4]*p.y + pH[5]) * invw;
}

template <class Type>
void UPIX_ALGORITHM_TEMPLATE Homography(const UpixMatrix3x3<Type> &H, const UpixPointHomg2D<Type> &p, UpixPointHomg2D<Type> &xp)
{
	Type *pH = H.GetData();
	xp.x = pH[0]*p.x + pH[1]*p.y + pH[2]*p.w;
	xp.y = pH[3]*p.x + pH[4]*p.y + pH[5]*p.w;
	xp.w = pH[6]*p.x + pH[7]*p.y + pH[8]*p.w;
}

template <class Type>
void UPIX_ALGORITHM_TEMPLATE ProjectPoint3D(const UpixMatrix<Type> &P, const UpixPoint3D<Type> &point_3d, UpixPoint2D<Type> &point_2d)
{
	Type *p = P.GetData();
	Type invw = 1/(p[8]*point_3d.x + p[9]*point_3d.y + p[10]*point_3d.z + p[11]);
	point_2d.x = (p[0]*point_3d.x + p[1]*point_3d.y + p[2]*point_3d.z + p[3]) * invw;
	point_2d.y = (p[4]*point_3d.x + p[5]*point_3d.y + p[6]*point_3d.z + p[7]) * invw;
}

template <class Type>
void UPIX_ALGORITHM_TEMPLATE ProjectPoint3D(const UpixMatrix<Type> &P, const UpixPointHomg3D<Type> &point_3d, UpixPointHomg2D<Type> &point_2d)
{
	Type *p = P.GetData();
	point_2d.x = p[0]*point_3d.x + p[1]*point_3d.y + p[2]*point_3d.z + p[3]*point_3d.w;
	point_2d.y = p[4]*point_3d.x + p[5]*point_3d.y + p[6]*point_3d.z + p[7]*point_3d.w;
	point_2d.w = p[8]*point_3d.x + p[9]*point_3d.y + p[10]*point_3d.z + p[11]*point_3d.w;
}



#endif