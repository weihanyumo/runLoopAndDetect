/*
* Filename: UpixFile.h
* Author:   Xiehong
* Date:     2012.10.10
*/

#ifndef _UPIX_FILE_H_
#define _UPIX_FILE_H_


#include "UpixAlgorithmTypeDef.h"


#ifdef _WIN32
#include <io.h>
#else
struct _finddata64i32_t {
	unsigned    attrib;
	long long  time_create;    /* -1 for FAT file systems */
	long long  time_access;    /* -1 for FAT file systems */
	long long  time_write;
	unsigned long    size;
	char        name[260];
};

typedef _finddata64i32_t _finddata_t;
#endif//_WIN32


#ifndef _WIN32
#define _MAX_PATH   260 /* max. length of full pathname */
#define _MAX_DRIVE  3   /* max. length of drive component */
#define _MAX_DIR    256 /* max. length of path component */
#define _MAX_FNAME  256 /* max. length of file name component */
#define _MAX_EXT    256 /* max. length of extension component */
#endif



// Return only extension
void UPIX_ALGORITHM_API UpixGetFileExtension(const char *strPathName, char *strExt);

// Return only file title
void UPIX_ALGORITHM_API UpixGetFileTitle(const char *strPathName, char *strFileTitle);

// Return file title with extension
void UPIX_ALGORITHM_API UpixGetFileName(const char *strPathName, char *strFileName);

// Return full path
void UPIX_ALGORITHM_API UpixGetFolderPath(const char *strPathName, char *strFolderPath);

// Check directory exist
bool UPIX_ALGORITHM_API UpixIsFolderExist(const char *strDir);

// Check file exist
bool UPIX_ALGORITHM_API UpixIsFileExist(const char *strPathName);

// Rename file
bool UPIX_ALGORITHM_API UpixRenameFile(const char *strOldName, const char *strNewName);

// Move file
bool UPIX_ALGORITHM_API UpixMoveFile(const char *strPathName, const char *strNewPath);

// Copy file
bool UPIX_ALGORITHM_API UpixCopyFile(const char *strPathName, const char *strNewPath);


// Find file in dir
#define UPIX_FIND_FILE_FAILED -1

typedef _finddata_t UpixFindData;

long UPIX_ALGORITHM_API UpixFindFirstFile(const char *strFileName, UpixFindData *pFindData);
bool UPIX_ALGORITHM_API UpixFindNextFile(long hFind, UpixFindData *pFindData);
bool UPIX_ALGORITHM_API UpixFindClose(long hFind);



#endif