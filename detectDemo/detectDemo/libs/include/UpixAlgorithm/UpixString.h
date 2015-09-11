/*
* Filename: UpixString.h
* Author:   Xiehong
* Date:     2012.11.20
*/

#ifndef _UPIX_STRING_H_
#define _UPIX_STRING_H_


#include "UpixAlgorithmTypeDef.h"


// Convert character to lower/upper
UPIX_INLINE char UPIX_ALGORITHM_API UpixCharLower(char x)
{
	return (x >= 'A' && x <= 'Z') ? (x + 0x20) : x;
}

UPIX_INLINE char UPIX_ALGORITHM_API UpixCharUpper(char x)
{
	return (x >= 'a' && x <= 'z') ? (x - 0x20) : x;
}


// Convert string to lower
void UPIX_ALGORITHM_API UpixStringLower(const char *src, char *dst);

// Convert string to upper
void UPIX_ALGORITHM_API UpixStringUpper(const char *src, char *dst);

// Compare strings
bool UPIX_ALGORITHM_API UpixIsStringEqual(const char *str1, const char *str2);



// Class for string
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixString
{
public:
	UpixString() : Array(NULL), Allocated(1), Used(1)
	{
		Array = new T[1];
		Array[0] = 0x0;
	}

	UpixString(const UpixString<T> &other) : Array(0), Allocated(0), Used(0)
	{
		*this = other;
	}

	template <class B>
	UpixString(const UpixString<B> &other) : Array(0), Allocated(0), Used(0)
	{
		*this = other;
	}

	explicit UpixString(const double number, const int decimalLen = 2) : Array(0), Allocated(0), Used(0)
	{
		char tmpbuf[256];
		char format[10];
		sprintf(format, "%%0.%df", decimalLen);
		sprintf(tmpbuf, format, number);
		*this = tmpbuf;
	}

	explicit UpixString(int number) : Array(0), Allocated(0), Used(0)
	{
		bool negative = false;
		if(number < 0)
		{
			number *= -1;
			negative = true;
		}

		char tmpbuf[16];
		tmpbuf[15] = 0;
		uint idx = 15;
		if(number == 0)
		{
			tmpbuf[14] = '0';
			*this = &tmpbuf[14];
			return;
		}

		while(number && idx)
		{
			idx--;
			tmpbuf[idx] = (char)('0' + (number % 10));
			number /= 10;
		}

		if(negative)
		{
			idx--;
			tmpbuf[idx] = '-';
		}
		*this = &tmpbuf[idx];
	}

	explicit UpixString(unsigned int number) : Array(0), Allocated(0), Used(0)
	{
		char tmpbuf[16];
		tmpbuf[15] = 0;
		uint idx = 15;
		if(number == 0)
		{
			tmpbuf[14] = '0';
			*this = &tmpbuf[14];
			return;
		}

		while(number && idx)
		{
			idx--;
			tmpbuf[idx] = (char)('0' + (number % 10));
			number /= 10;
		}
		*this = &tmpbuf[idx];
	}

	template <class B>
	UpixString(const B* const c, uint length) : Array(0), Allocated(0), Used(0)
	{
		if(c == NULL)
		{
			*this = "";
			return;
		}

		Allocated = Used = length + 1;
		Array = new T[Used];
		for(uint l=0; l<length; l++)
		{
			Array[l] = (T)c[l];
		}
		Array[length] = 0;
	}

	template <class B>
	UpixString(const B* const c) : Array(0), Allocated(0), Used(0)
	{
		*this = c;
	}

	~UpixString()
	{
		if (Array != NULL)
		{
			delete[] Array;
			Array = NULL;
		}
	}

	uint Size() const
	{
		return Used - 1;
	}

	const T* C_str() const
	{
		return Array;
	}

	UpixString<T>& operator = (const UpixString<T> &other)
	{
		if(this != &other)
		{
			if (Array != NULL)
			{
				delete[] Array;
				Array = NULL;
			}
		
			Allocated = Used = other.Used;
			Array = new T[Used];
			const T *p = other.Array;
			for(uint i=0; i<Used; i++, p++)
			{
				Array[i] = *p;
			}
		}
		return *this;
	}

	template <class B>
	UpixString<T>& operator = (const UpixString<B> &other)
	{
		*this = other.Array;
		return *this;
	}

	template <class B>
	UpixString<T>& operator = (const B* const c)
	{
		if(c == NULL)
		{
			if(Array == NULL)
			{
				Array = new T[1];
				Allocated = 1;
			}
			Used = 1;
			Array[0] = 0x0;
			return *this;
		}

		if((void*)c == (void*)Array)
			return *this;

		uint len = 0;
		const B *p = c;
		while(*p)
		{
			len++;
			p++;
		}
		len++;

		T *oldArray = Array;
		Allocated = Used = len;
		Array = new T[Used];
		for(uint l=0; l<len; l++)
			Array[l] = (T)c[l];

		delete[] oldArray;
		return *this;
	}

	void Append(T character)
	{
		if(Used + 1 > Allocated)
			ReAllocate(Used + 1);
		Used++;
		Array[Used - 2] = character;
		Array[Used - 1] = 0;
	}

	void Append(const T* const other)
	{
		if(other == NULL)
			return;

		uint len = 0;
		const T *p = other;
		while(*p)
		{
			len++;
			p++;
		}

		if(Used + len > Allocated)
			ReAllocate(Used + len);
		Used--;
		len++;

		for(uint l=0; l<len; l++)
		{
			Array[l + Used] = other[l];
		}
		Used += len;
	}

	void Append(const UpixString<T> &other)
	{
		Used--;
		uint len = other.Used;
		if(Used + len > Allocated)
			ReAllocate(Used + len);

		for(uint l=0; l<len; l++)
			Array[Used + l] = other[l];

		Used += len;
	}

	void Append(const UpixString<T> &other, uint length)
	{
		if(other.Size() < length)
		{
			Append(other);
			return;
		}
		if(Used + length > Allocated)
			ReAllocate(Used + length);
		Used--;

		for(uint l=0; l<length; l++)
			Array[l + Used] = other[l];

		Used += length;
		Array[Used] = 0;
		Used++;
	}

	void Reserve(uint count)
	{
		if(count < Allocated)
			return;
		ReAllocate(count);
	}

	UpixString<T> operator + (const UpixString<T> &other) const
	{
		UpixString<T> str(*this);
		str.Append(other);
		return str;
	}

	template <class B>
	UpixString<T> operator + (const B* const c) const
	{
		UpixString<T> str(*this);
		str.Append(c);
		return str;
	}

	T& operator [] (const uint index)
	{
		assert(index < Used);
		return Array[index];
	}

	const T& operator [] (const uint index) const
	{
		assert(index < Used);
		return Array[index];
	}

	bool operator == (const T* const str) const
	{
		if(str == NULL)
			return false;

		uint i;
		for(i=0; Array[i] && str[i]; i++)
			if(Array[i] != str[i])
				return false;

		return !Array[i] && !str[i];
	}

	bool operator == (const UpixString<T> &other) const
	{
		for(uint i=0; Array[i] && other.Array[i]; i++)
			if(Array[i] != other.Array[i])
				return false;

		return Used == other.Used;
	}

	bool operator < (const UpixString<T> &other) const
	{
		for(uint i=0; Array[i] && other.Array[i]; i++)
		{
			int diff = Array[i] - other.Array[i];
			if(diff)
				return diff < 0;
		}
		return Used < other.Used;
	}

	bool operator != (const T* const str) const
	{
		return !(*this == str);
	}

	bool operator != (const UpixString<T>& other) const
	{
		return !(*this == other);
	}

	void MakeLower()
	{
		for(uint i=0; i<Used; i++)
			Array[i] = AnsiLower(Array[i]);
	}

	void MakeUpper()
	{
		const T a = (T)'a';
		const T z = (T)'z';
		const T diff = (T)'A' - a;
		for(uint i=0; i<Used; i++)
		{
			if(Array[i] >= a && Array[i] <= z)
				Array[i] += diff;
		}
	}

	bool EqualsIgnoreCase(const UpixString<T> &other) const
	{
		for(uint i=0; Array[i] && other[i]; i++)
		{
			if(AnsiLower(Array[i]) != AnsiLower(other[i]))
				return false;
		}
		return Used == other.Used;
	}

	bool LowerIgnoreCase(const UpixString<T> &other) const
	{
		for(uint i=0; Array[i] && other.Array[i]; i++)
		{
			int diff = (int)AnsiLower(Array[i]) - (int)AnsiLower(other.Array[i]);
			if(diff < 0)
				return true;
		}
		return Used < other.Used;
	}

	bool EqualsN(const UpixString<T> &other, uint n) const
	{
		uint i;
		for(i=0; Array[i] && other[i] && i < n; i++)
			if(Array[i] != other[i])
				return false;
		return (i == n) || (Used == other.Used);
	}

	bool EqualsN(const T* const str, uint n) const
	{
		if(str == NULL)
			return false;

		uint i;
		for(i=0; Array[i] && str[i] && i < n; i++)
			if(Array[i] != str[i])
				return false;
		return (i == n) || (Array[i] == 0 && str[i] == 0);
	}

	int FindFirst(T c) const
	{
		for(uint i=0; i<Used; i++)
			if(Array[i] == c)
				return i;
		return -1;
	}

	int FindFirstChar(const T* const c, uint count) const
	{
		if(c == NULL)
			return -1;

		for(uint i=0; i<Used; i++)
		{
			for(uint j=0; j<count; j++)
				if(Array[i] == c[j])
					return i;
		}
		return -1;
	}

	template <class B>
	int FindFirstCharNotInList(const B* const c, uint count) const
	{
		for(uint i=0; i<Used - 1; i++)
		{
			uint j;
			for(j=0; j<count; j++)
				if(Array[i] == c[j])
					break;
			if(j == count)
				return i;
		}
		return -1;
	}

	template <class B>
	int FindLastCharNotInList(const B* const c, uint count) const
	{
		for(int i = (int)(Used - 2); i >= 0; i--)
		{
			uint j;
			for(j=0; j<count; j++)
				if(Array[i] == c[j])
					break;
			if(j == count)
				return i;
		}
		return -1;
	}

	int FindNext(T c, uint startPos) const
	{
		for(uint i=startPos; i<Used; i++)
			if(Array[i] == c)
				return i;
		return -1;
	}

	int FindLast(T c, int start = -1) const
	{
		start = UpixClamp(start < 0 ? (int)(Used) - 1 : start, 0, (int)(Used) - 1);
		for(int i=start; i>=0; i--)
			if(Array[i] == c)
				return i;
		return -1;
	}

	int FindLastChar(const T* const c, uint count) const
	{
		if(c == NULL)
			return -1;

		for(int i=Used - 1; i>=0; i--)
			for(uint j=0; j<count; j++)
				if(Array[i] == c[j])
					return i;
		return -1;
	}

	template <class B>
	int Find(const B* const str) const
	{
		if(str && *str)
		{
			uint len = 0;
			while(str[len])
				len++;
			if(len > Used - 1)
				return -1;
			for(uint i=0; i<Used - len; i++)
			{
				uint j = 0;
				while(str[j] && Array[i + j] == str[j])
					j++;
				if(!str[j])
					return i;
			}
		}
		return -1;
	}

	UpixString<T> SubString(uint begin, int length) const
	{
		if(length <= 0 || begin >= Size())
			return UpixString<T>("");
		if((length + begin) > Size())
			length = Size() - begin;
		UpixString<T> o;
		o.Reserve(length + 1);
		for(int i = 0; i < length; i++)
			o.Array[i] = Array[i + begin];
		o.Array[length] = 0;
		o.Used = o.Allocated;
		return o;
	}

	UpixString<T>& operator += (T c)
	{
		Append(c);
		return *this;
	}

	UpixString<T>& operator += (const T* const c)
	{
		Append(c);
		return *this;
	}

	UpixString<T>& operator += (const UpixString<T> &other)
	{
		Append(other);
		return *this;
	}

	UpixString<T>& operator += (const int i)
	{
		Append(UpixString<T>(i));
		return *this;
	}

	UpixString<T>& operator += (const unsigned int i)
	{
		Append(UpixString<T>(i));
		return *this;
	}

	UpixString<T>& operator += (const long i)
	{
		Append(UpixString<T>(i));
		return *this;
	}

	UpixString<T>& operator += (const unsigned long& i)
	{
		Append(UpixString<T>(i));
		return *this;
	}

	UpixString<T>& operator += (const double i)
	{
		Append(UpixString<T>(i));
		return *this;
	}

	UpixString<T>& operator += (const float i)
	{
		Append(UpixString<T>(i));
		return *this;
	}

	void Replace(T toReplace, T replaceWith)
	{
		for(uint i=0; i<Used; i++)
			if(Array[i] == toReplace)
				Array[i] = replaceWith;
	}

	UpixString<T>& Trim(const UpixString<T> &whitespace = " \t\n\r")
	{
		const int begin = FindFirstCharNotInList(whitespace.C_str(), whitespace.Used);
		if(begin == -1)
			return (*this = "");
		const int end = FindLastCharNotInList(whitespace.C_str(), whitespace.Used);
		return (*this = SubString(begin, (end + 1) - begin));
	}

	void Erase(uint index)
	{
		assert(index < Used);
		for(uint i=index + 1; i<Used; i++)
			Array[i - 1] = Array[i];
		Used--;
	}

private:
	inline T AnsiLower(uint x) const
	{
		return x >= 'A' && x <= 'Z' ? (T) x + 0x20 : (T) x;
	}

	void ReAllocate(uint new_size)
	{
		T* old_array = Array;
		Array = new T[new_size];
		Allocated = new_size;
		uint amount = Used < new_size ? Used : new_size;

		for(uint i=0; i<amount; i++)
			Array[i] = old_array[i];

		if(Allocated < Used)
			Used = Allocated;
		delete[] old_array;
	}

	T* Array;
	uint Allocated;
	uint Used;
};


typedef UpixString<char> UpixStringc;
typedef UpixString<wchar_t> UpixStringw;



#endif