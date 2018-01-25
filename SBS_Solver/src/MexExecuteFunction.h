// Create by Amir Michaelis.
// Modify $Id: MexExecuteFunction.h 335 2017-03-06 17:29:04Z amirm $

#pragma once

// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the MEXEXECUTEFUNCTION_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// MEXEXECUTEFUNCTION_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.

#ifdef MEX_EXECUTE_FUNCTION_EXPORTS
#define MEX_EXECUTE_FUNCTION_API __declspec(dllexport)
#else
#define MEX_EXECUTE_FUNCTION_API __declspec(dllimport)
#endif

#define DLL_EXPORT_SYM MEX_EXECUTE_FUNCTION_API

/* ****************** Mex Mediator and Converters ****************** */

#include <cstdint>
#include <type_traits>
#include <typeinfo>
#include <string>
#include <map>
#include <vector>
#include <list>

#include "mex.h"
#include "matrix.h"

/* ****************** Mex Atom Class ****************** */

class MxStruct {
	mxArray * ptrMxStruct;
	unsigned long structSize[2];

protected:
	const unsigned long sub2ind(const unsigned long row = 0, const unsigned long column = 0) const;

public:
	MxStruct();
	MxStruct(const std::vector<std::string> & fieldNames, const unsigned long row, const unsigned long column = 1);
	MxStruct(const unsigned long row, const unsigned long column = 1);
	MxStruct(mxArray * initPtr);
	MxStruct(const MxStruct & src);
	MxStruct & operator=(const MxStruct & src);

	mxArray * GetMxArray(void);
	const mxArray * GetMxArray(void) const;

	mxArray * operator() (const char * fieldName, const unsigned long row = 0, const unsigned long column = 0);
	mxArray * operator() (const std::string & fieldName, const unsigned long row = 0, const unsigned long column = 0);
	const mxArray * operator() (const char * fieldName, const unsigned long row = 0, const unsigned long column = 0) const;
	const mxArray * operator() (const std::string & fieldName, const unsigned long row = 0, const unsigned long column = 0) const;
	const std::vector<std::string> GetFieldsName(void) const;
	const unsigned long * GetSize(void) const;
	const unsigned long GetLength(void) const;
	const unsigned long GetNumElements(void) const;
	MxStruct & SetSize(const unsigned long row, const unsigned long column = 1);
	MxStruct & Set(const unsigned long row = 0, const unsigned long column = 1);
	MxStruct & Set(mxArray * value, const std::string & fieldName, const unsigned long row = 0, const unsigned long column = 0);
	MxStruct & Set(mxArray * value, const char * fieldName, const unsigned long row = 0, const unsigned long column = 0);
};

class MxCell {
	mxArray * ptrMxCell;
	unsigned long cellSize[2];

protected:
	const unsigned long sub2ind(const unsigned long row = 0, const unsigned long column = 0) const;

public:
	MxCell();
	MxCell(const unsigned long row, const unsigned long column = 1);
	MxCell(mxArray * initPtr);
	MxCell(const MxCell & src);
	MxCell & operator= (const MxCell & src);

	mxArray * GetMxArray(void);
	const mxArray * GetMxArray(void) const;

	mxArray * operator() (const unsigned long row, const unsigned long column = 0);
	const mxArray * operator() (const unsigned long row, const unsigned long column = 0) const;
	const unsigned long * GetSize(void) const;
	const unsigned long GetLength(void) const;
	const unsigned long GetNumElements(void) const;
	MxCell & SetSize(const unsigned long row, const unsigned long column = 1);
	MxCell & Set(const unsigned long row = 0, const unsigned long column = 1);
	MxCell & Set(mxArray * value, const unsigned long row, const unsigned long column = 0);
};

class MxScalar {
	mxArray * ptrMxScalar;
public:
	MxScalar(double iniValue);
	MxScalar(mxArray * iniPtr = nullptr);

	mxArray * GetMxArray(void);
	const mxArray * GetMxArray(void) const;

	double & operator() (void);
	const double & operator() (void) const;
	MxScalar & Set(const double value = 0);
};

class MxString {
	mxArray * ptrMxString;
	std::string stringFieldData;
public:
	MxString(const char *);
	MxString(const std::string & iniStr);
	MxString(mxArray * iniPtr = nullptr);
	~MxString();

	mxArray * GetMxArray(void);
	const mxArray * GetMxArray(void) const;

	std::string & operator() (void);
	const char * operator() (void) const;
	MxString & Set(const char * value = nullptr);
	MxString & Set(const std::string & value);
};

class MxMatrix {
	mxArray * ptrMxMatrix;
	unsigned long matrixSize[2];

protected:
	const unsigned long sub2ind(const unsigned long row = 0, const unsigned long column = 0) const;

public:
	MxMatrix();
	MxMatrix(const unsigned long row, const unsigned long column = 1);
	MxMatrix(mxArray * initPtr);
	MxMatrix(const MxMatrix & src);
	MxMatrix & operator= (const MxMatrix & src);

	mxArray * GetMxArray(void);
	const mxArray * GetMxArray(void) const;

	double * operator() (void);
	const double * operator() (void) const;
	double & operator() (const unsigned long row, const unsigned long column = 0);
	const double & operator() (const unsigned long row, const unsigned long column = 0) const;
	const unsigned long * GetSize(void) const;
	const unsigned long GetLength(void) const;
	const unsigned long GetNumElements(void) const;
	MxMatrix & SetSize(const unsigned long row, const unsigned long column = 1);
	MxMatrix & Set(const unsigned long row = 0, const unsigned long column = 1);
	MxMatrix & Set(const double value, const unsigned long row, const unsigned long column = 0);
};

/* ****************** Code Segments ****************** */
/*
MxAtomTypes::MxAtomTypes()
{

}

MxAtomTypes::MxAtomTypes(MxAtomTypes & init)
{

}

MxAtomTypes & MxAtomTypes::operator= (const MxAtomTypes &)
{
return (*this)
}

void MxAtomTypes::SetMxArray(mxArray * iniMxArray)
{
ptrMxArray = iniMxArray;
}

mxArray * & MxAtomTypes::SetMxArray(void)
{
return (ptrMxArray);
}

const mxArray * MxAtomTypes::GetMxArray(void) const
{
return(ptrMxArray);
}

struct MxAtomTypes {

mxArray * ptrMxArray;
public:
MxAtomTypes();
MxAtomTypes(MxAtomTypes &);
MxAtomTypes & operator= (const MxAtomTypes &);
~MxAtomTypes();

void SetMxArray (mxArray *);
mxArray * & SetMxArray(void);
const mxArray * GetMxArray (void) const;

MxStruct Struct;
MxCell Cell;
MxScalar Scalar;
MxString String;
MxMatrix Vector;
MxMatrix Matrix;
};

class MxVector {
mxArray * ptrMxVector;
public:
MxVector(mxArray * initPtr);

mxArray * GetMxArray(void);
const mxArray * GetMxArray(void) const;

double * operator() (void);
const double * operator() (void) const;
double & operator() (const unsigned long index);
const double & operator() (const unsigned long index) const;
const unsigned long GetSize(void) const;
MxVector & SetSize(unsigned long length);
MxVector & Set(const unsigned long index, const double value);
};

template<class T>
class MexScalar {
	mxArray * mxArrayPtr;
public:
	T Value(void) const;
	T & operator() (void);
	MexScalar & operator =(const T &);
};

class MexStruct {

public:

};

class MexCell {

public:

};

class MexVector {

public:

};

class MexMatrix {


};


template<class T>
T MexScalar<T>::Value(void) const
{

}


template<class T>
class MxMediator {
const mxArray * ptrMxArray;
T * ptrT;
unsigned long nElementsT;
public:
MxMediator();
MxMediator(const MxMediator &);
MxMediator & operator = (const MxMediator &);

MxMediator(const mxArray *);
~MxMediator();
T * GetArray(void);
unsigned long GetArraySize(void);
T * GetStructFiled(const std::string &);
T * GetStructFiled(const char *);
T * operator()(const std::string &);
T * operator()(const char *);
T * operator()(void);
};

// MexExecuteFunction.cpp.
//

template<class T>
MxMediator<T>::MxMediator()
{
ptrMxArray = nullptr;
ptrT = nullptr;
nElementsT = 0;
}

template<class T>
MxMediator<T>::MxMediator(const MxMediator<T> & initMxMediator)
{
ptrMxArray = initMxMediator.ptrMxArray;
ptrT = initMxMediator.ptrT;
nElementsT = initMxMediator.nElementsT;
}

template<class T>
MxMediator<T> & MxMediator<T>::operator =(const MxMediator<T> & src)
{
if (&src == this)
return (*this);

ptrMxArray = src.ptrMxArray;
ptrT = src.ptrT;
nElementsT = src.nElementsT;
return (*this);
}

template<class T>
MxMediator<T>::MxMediator(const mxArray * initMxArray)
{
ptrMxArray = initMxArray;
ptrT = nullptr;
nElementsT = 0;
}

template<class T>
MxMediator<T>::~MxMediator()
{
delete ptrT;
nElementsT = 0;
}

template<class T>
T * MxMediator<T>::GetArray(void)
{
if (!mxIsNumeric(ptrMxArray))
return nullptr;

// std::string mxTypeName( mxGetClassName(ptrMxArray) );
mxClassID category = mxGetClassID(ptrMxArray);

void * mxData = mxGetData(ptrMxArray);
unsigned long rows = mxGetM(ptrMxArray);
unsigned long columns = mxGetN(ptrMxArray);
unsigned long nElements = rows*columns;
ptrT = new T[nElements];
nElementsT = nElements;
for (unsigned int iData = 0; iData < nElements; iData++)
switch (category) {
case mxLOGICAL_CLASS:
*(ptrT + iData) = T(*(static_cast<bool *>(mxData) + iData));
break;
case mxCHAR_CLASS:
*(ptrT + iData) = T(*(static_cast<char *>(mxData) + iData));
break;
case mxDOUBLE_CLASS:
*(ptrT + iData) = T(*(static_cast<double *>(mxData) + iData));
break;
case mxSINGLE_CLASS:
*(ptrT + iData) = T(*(static_cast<float *>(mxData) + iData));
break;
case mxINT8_CLASS:
*(ptrT + iData) = T(*(static_cast<std::int8_t *>(mxData) + iData));
break;
case mxUINT8_CLASS:
*(ptrT + iData) = T(*(static_cast<std::uint8_t *>(mxData) + iData));
break;
case mxINT16_CLASS:
*(ptrT + iData) = T(*(static_cast<std::int16_t *>(mxData) + iData));
case mxUINT16_CLASS:
*(ptrT + iData) = T(*(static_cast<std::uint16_t *>(mxData) + iData));
break;
case mxINT32_CLASS:
*(ptrT + iData) = T(*(static_cast<std::int32_t *>(mxData) + iData));
break;
case mxUINT32_CLASS:
*(ptrT + iData) = T(*(static_cast<std::uint32_t *>(mxData) + iData));
break;
case mxINT64_CLASS:
*(ptrT + iData) = T(*(static_cast<std::int64_t *>(mxData) + iData));
break;
case mxUINT64_CLASS:
*(ptrT + iData) = T(*(static_cast<std::uint64_t *>(mxData) + iData));
break;
default:
mexErrMsgIdAndTxt("MATLAB:MexExecuteFunction:inputClassMismatch",
"Input class type currently not supported.");
break;
};

return(ptrT);
}

template<class T>
unsigned long MxMediator<T>::GetArraySize(void)
{
return(nElementsT);
}

template<class T>
T * MxMediator<T>::GetStructFiled(const std::string & fieldName)
{
if (!mxIsStruct(ptrMxArray))
return nullptr;

mwSize nStruct = mxGetNumberOfElements(ptrMxArray);
mwSize iStruct = 0;
if (nStruct != 1)
mexErrMsgIdAndTxt("MATLAB:MexExecuteFunction:MxMediator:sizeMismatch",
"Field numel must be one otherwise class currently not supported.");

mxArray * ptrFieldName = mxGetField(ptrMxArray, iStruct, fieldName.c_str());
ptrT = new T(ptrFieldName);
nElementsT = 1;

return(ptrT);
}

template<class T>
T * MxMediator<T>::GetStructFiled(const char * fieldName)
{
return (GetStructFiled(std::string(fieldName)));
}

template<class T>
T * MxMediator<T>::operator()(const std::string & fieldName)
{
return (GetStructFiled(fieldName));
}

template<class T>
T * MxMediator<T>::operator()(const char * fieldName)
{
return (GetStructFiled(fieldName));
}

template<class T>
T * MxMediator<T>::operator()(void)
{
return (GetArray());
}



*/

/* ****************** MxMediator ****************** */
/*
template<class T>
class MxMediator {
	const mxArray * ptrMxArray;
	T * ptrT;
	unsigned long nElementsT;
public:
	MxMediator();
	MxMediator(const MxMediator &);
	MxMediator & operator = (const MxMediator &);

	MxMediator(const mxArray *);
	~MxMediator();
	T * GetArray(void);
	unsigned long GetArraySize(void);
	T * operator()(void);
};

template<class T>
MxMediator<T>::MxMediator()
{
	ptrMxArray = nullptr;
	ptrT = nullptr;
	nElementsT = 0;
}

template<class T>
MxMediator<T>::MxMediator(const MxMediator<T> & initMxMediator)
{
	ptrMxArray = initMxMediator.ptrMxArray;
	ptrT = initMxMediator.ptrT;
	nElementsT = initMxMediator.nElementsT;
}

template<class T>
MxMediator<T> & MxMediator<T>::operator =(const MxMediator<T> & src)
{
	if (&src == this)
		return (*this);

	delete ptrT;

	ptrMxArray = src.ptrMxArray;
	ptrT = src.ptrT;
	nElementsT = src.nElementsT;
	return (*this);
}

template<class T>
MxMediator<T>::MxMediator(const mxArray * initMxArray)
{
	ptrMxArray = initMxArray;
	ptrT = nullptr;
	nElementsT = 0;
}

template<class T>
MxMediator<T>::~MxMediator()
{
	delete ptrT;
	nElementsT = 0;
}

template<class T>
T * MxMediator<T>::GetArray(void)
{
	if (!mxIsNumeric(ptrMxArray))
		return nullptr;

	// std::string mxTypeName( mxGetClassName(ptrMxArray) );
	mxClassID category = mxGetClassID(ptrMxArray);

	void * mxData = mxGetData(ptrMxArray);
	unsigned long rows = mxGetM(ptrMxArray);
	unsigned long columns = mxGetN(ptrMxArray);
	unsigned long nElements = rows*columns;
	ptrT = new T[nElements];
	nElementsT = nElements;
	for (unsigned int iData = 0; iData < nElements; iData++)
		switch (category) {
		case mxLOGICAL_CLASS:
			*(ptrT + iData) = T(*(static_cast<bool *>(mxData) + iData));
			break;
		case mxCHAR_CLASS:
			*(ptrT + iData) = T(*(static_cast<char *>(mxData) + iData));
			break;
		case mxDOUBLE_CLASS:
			*(ptrT + iData) = T(*(static_cast<double *>(mxData) + iData));
			break;
		case mxSINGLE_CLASS:
			*(ptrT + iData) = T(*(static_cast<float *>(mxData) + iData));
			break;
		case mxINT8_CLASS:
			*(ptrT + iData) = T(*(static_cast<std::int8_t *>(mxData) + iData));
			break;
		case mxUINT8_CLASS:
			*(ptrT + iData) = T(*(static_cast<std::uint8_t *>(mxData) + iData));
			break;
		case mxINT16_CLASS:
			*(ptrT + iData) = T(*(static_cast<std::int16_t *>(mxData) + iData));
		case mxUINT16_CLASS:
			*(ptrT + iData) = T(*(static_cast<std::uint16_t *>(mxData) + iData));
			break;
		case mxINT32_CLASS:
			*(ptrT + iData) = T(*(static_cast<std::int32_t *>(mxData) + iData));
			break;
		case mxUINT32_CLASS:
			*(ptrT + iData) = T(*(static_cast<std::uint32_t *>(mxData) + iData));
			break;
		case mxINT64_CLASS:
			*(ptrT + iData) = T(*(static_cast<std::int64_t *>(mxData) + iData));
			break;
		case mxUINT64_CLASS:
			*(ptrT + iData) = T(*(static_cast<std::uint64_t *>(mxData) + iData));
			break;
		default:
			mexErrMsgIdAndTxt("MATLAB:MexExecuteFunction:inputClassMismatch",
				"Input class type currently not supported.");
			break;
		};

	return(ptrT);
}

template<class T>
unsigned long MxMediator<T>::GetArraySize(void)
{
	return(nElementsT);
}

template<class T>
T * MxMediator<T>::operator()(void)
{
	return (GetArray());
}
*/

/* ****************** DLL Example ****************** */
/*
// This class is exported from the MexExecuteFunction.dll
class MEX_EXECUTE_FUNCTION_API CMexExecuteFunction {
public:
	CMexExecuteFunction(void);
	// TODO: add your methods here.
};

extern MEX_EXECUTE_FUNCTION_API int nMexExecuteFunction;

MEX_EXECUTE_FUNCTION_API int fnMexExecuteFunction(void);

// This is an example of an exported variable
MEX_EXECUTE_FUNCTION_API int nMexExecuteFunction = 0;

// This is an example of an exported function.
MEX_EXECUTE_FUNCTION_API int fnMexExecuteFunction(void)
{
	return 42;
}

// This is the constructor of a class that has been exported.
// see MexExecuteFunction.h for the class definition
CMexExecuteFunction::CMexExecuteFunction()
{
	return;
}

*/

/* ****************** CPP Example ****************** */
/*
// MexExecuteFunction.cpp : Defines the exported functions for the DLL application.
//

#include <string>
#include <map>
#include <vector>

#include "MexExecuteFunction.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
// varargout = cell(nargout, 1);
std::string mexAction;
int used_nlhs = 0;

// Examine input (right-hand-side) arguments.
if (nrhs <= 0 || !mxIsChar(prhs[0]))
mexErrMsgIdAndTxt("MATLAB:MexExecuteFunction:inputOutputMismatch",
	"Usage: MexExecuteFunction('Action', varargin).\n");
char * mexActionStr = mxArrayToString(prhs[0]);
mexAction = mexActionStr;
mxFree(mexActionStr);

if (mexAction == "Init") {
	mexPrintf("Do Init.\n");

	if (nrhs != 1)
		mexErrMsgIdAndTxt("MATLAB:MexExecuteFunction:Init:inputOutputMismatch",
			"Usage: MexExecuteFunction('Init').\n");
	// Create Class And Init.

	if (nlhs > used_nlhs)
		plhs[used_nlhs++] = mxCreateLogicalScalar(true);
};

if (mexAction == "Calc") {
	mexPrintf("Do Calc.\n");

	if (nrhs != 4)
		mexErrMsgIdAndTxt("MATLAB:MexExecuteFunction:inputOutputMismatch",
			"Usage: message = MexExecuteFunction('Calc', SignalVec, TimeVec, messageCols).\n");

	// Run the class finc.

	if (nlhs > used_nlhs)
		plhs[used_nlhs++] = mxCreateString("return message board");
};

if (mexAction == "Get Internal Signals") {
	mexPrintf("Get Internal Signals.\n");
	if (nrhs != 1)
		mexErrMsgIdAndTxt("MATLAB:MexExecuteFunction:inputOutputMismatch",
			"Usage: outputSignalStruct = MexExecuteFunction('Get Internal Signals').\n");
	// Return internal signals.

	if (nlhs > used_nlhs)
		plhs[used_nlhs++] = mxCreateString("return output Signal Struct");
};

for (mwSize ilhs = used_nlhs; ilhs < nlhs; ilhs++)
	plhs[ilhs] = mxCreateCellMatrix(1, 1);

}
*/
