/*
* Filename: UpixGetData.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_GETDATA_H_
#define _UPIX_GETDATA_H_

#include "UpixAlgorithmTypeDef.h"

#include <fstream>
#include <string>
#include <sstream> 
#include <vector>
using namespace std;


// Get a line data from file
template<class T>
int UPIX_ALGORITHM_TEMPLATE _getData(ifstream &ifs, vector<T> &DataArray)
{
	DataArray.clear();
	int DataNum = 0;
	T tmp;

	string line;
	stringstream line_stream;

	if (ifs.eof())
	{
		return 0;
	}

	getline(ifs, line);
	line_stream<<line;   

	while( line_stream>>tmp )
	{
		DataArray.push_back(tmp);
		DataNum++;
	}

	line_stream.str() = "";
	line_stream.clear();

	return DataNum;
}

UPIX_INLINE int UPIX_CDECL getData(ifstream &ifs, vector<double> &DataArray)
		{ return _getData(ifs, DataArray); }

UPIX_INLINE int UPIX_CDECL getData(ifstream &ifs, vector<int> &DataArray)
		{ return _getData(ifs, DataArray); }

UPIX_INLINE int UPIX_CDECL getData(ifstream &ifs, vector<float> &DataArray)
		{ return _getData(ifs, DataArray); }


#endif