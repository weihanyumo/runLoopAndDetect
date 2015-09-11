/*
* Filename: UpixDBSCAN.h
* Author:   Xiehong
* Date:     2013.8.19
*/
#ifndef _UPIX_DBSCAN_H_
#define _UPIX_DBSCAN_H_


#include "UpixMatrix.h"
#include "UpixArray.h"
#include "UpixPoint.h"



class UPIX_ALGORITHM_TEMPLATE UpixDBSCANCluster 
{
public:
	UpixDBSCANCluster() : nClusterId(-1) {}
	~UpixDBSCANCluster() {}

public:
	bool operator<(const UpixDBSCANCluster &rhs) const
	{
		return vIndex.Size() < rhs.vIndex.Size();
	}

	bool operator>(const UpixDBSCANCluster &rhs) const
	{
		return vIndex.Size() > rhs.vIndex.Size();
	}

public:
	UpixArray<uint> vIndex; // 类簇中的样本点
	int nClusterId;
};




// Base class of cluster data set
class UPIX_ALGORITHM_TEMPLATE IUpixDBSCAN
{
public:
	IUpixDBSCAN();
	virtual ~IUpixDBSCAN();

public:
	bool Init(const UpixMatrix<float> &matData, float radius, int minPts, bool bNormalize = true);
	virtual int ClusterAnalysis(UpixArray<UpixDBSCANCluster> &vClusters) = 0;

	static IUpixDBSCAN *Create(bool bRecursive);
	static void Release(IUpixDBSCAN *&pClusterer);

protected:
	void NormalizeData();
	void CreateSqrDistanceTable(UpixMatrix<float> &matDistanceTable);
	float GetSqrDistance(uint nId1, uint nId2);

protected:
	UpixMatrix<float> m_matDataSet;
	uint m_nDimension;						//维度
	uint m_nDataCount;						//数据数量
	uint m_nClusterCount;					// The number of cluster count
	uint m_nMinPointCount;					//邻域最小数据个数
	float m_fRadius;						//半径
};




#endif