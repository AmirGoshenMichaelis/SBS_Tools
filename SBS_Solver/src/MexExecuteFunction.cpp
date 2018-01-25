// Create by Amir Michaelis.
// Modify $Id: MexExecuteFunction.cpp 335 2017-03-06 17:29:04Z amirm $

#include "MexExecuteFunction.h"

// TODO
// Make this as MexApi.dll;

/* ****************** Mex Struct Class ****************** */

MxStruct::MxStruct() : ptrMxStruct(nullptr)
{
	structSize[0] = 0;
	structSize[1] = 0;
}

MxStruct::MxStruct(const std::vector<std::string> & fieldNames, const unsigned long row, const unsigned long column)
{
	structSize[0] = row;
	structSize[1] = column;
	unsigned int nFieldNames = static_cast<unsigned int>(fieldNames.size());

	ptrMxStruct = mxCreateStructMatrix(row, column, 0, nullptr);
	for (unsigned long iField = 0; iField < nFieldNames; ++iField)
		if (mxAddField(ptrMxStruct, fieldNames[iField].c_str()) == -1)
			throw "Can't create field name";
}

MxStruct::MxStruct(const unsigned long row, const unsigned long column) : ptrMxStruct(nullptr)
{
	structSize[0] = row;
	structSize[1] = column;
	ptrMxStruct = mxCreateStructMatrix(row, column, 0, nullptr);
}

MxStruct::MxStruct(mxArray * initPtr) : ptrMxStruct(nullptr), structSize{ 0, 0 }
{
	if (initPtr!=nullptr && mxIsStruct(initPtr)) {
		ptrMxStruct = initPtr;
		structSize[0] = static_cast<unsigned long>(mxGetM(ptrMxStruct));
		structSize[1] = static_cast<unsigned long>(mxGetN(ptrMxStruct));
	};
}

MxStruct::MxStruct(const MxStruct & src) : ptrMxStruct(nullptr), structSize{ 0, 0 }
{
	(*this) = src;
}

MxStruct & MxStruct::operator=(const MxStruct & src)
{
	if (this == &src)
		return (*this);

	structSize[0] = src.structSize[0];
	structSize[1] = src.structSize[1];

	mxDestroyArray(ptrMxStruct);
	if (src.ptrMxStruct == nullptr)
		ptrMxStruct = nullptr;
	else
		ptrMxStruct = mxDuplicateArray(src.ptrMxStruct);

	return (*this);
}

mxArray * MxStruct::GetMxArray(void)
{
	return(ptrMxStruct);
}

const mxArray * MxStruct::GetMxArray(void) const
{
	return(ptrMxStruct);
}

const unsigned long MxStruct::sub2ind(const unsigned long row, const unsigned long column) const
{
	if (ptrMxStruct == nullptr)
		return(0);

	mwSize subs[2];
	mwSize nsubs = 2;

	subs[0] = row;
	subs[1] = column;

	unsigned long index = mxCalcSingleSubscript(ptrMxStruct, nsubs, subs);
	return (index);
}

mxArray * MxStruct::operator() (const char * fieldName, const unsigned long row, const unsigned long column)
{
	if (ptrMxStruct == nullptr)
		return (nullptr);

	bool isExistFieldName = mxGetFieldNumber(ptrMxStruct, fieldName) != -1;

	if (!isExistFieldName)
		if (mxAddField(ptrMxStruct, fieldName) == -1)
			throw "Can't create field name";

	mxArray * fieldMxArray = mxGetField(ptrMxStruct, sub2ind(row, column), fieldName);

	return (fieldMxArray);
}

mxArray * MxStruct::operator() (const std::string & fieldName, const unsigned long row, const unsigned long column)
{
	return((*this)(fieldName.c_str(), row, column));
}

const mxArray * MxStruct::operator() (const char * fieldName, const unsigned long row, const unsigned long column) const
{
	if (ptrMxStruct == nullptr)
		return (nullptr);

	mxArray * fieldMxArray = mxGetField(ptrMxStruct, sub2ind(row, column), fieldName);

	return (fieldMxArray);
}

const mxArray * MxStruct::operator() (const std::string & fieldName, const unsigned long row, const unsigned long column) const
{
	return((*this)(fieldName.c_str(), row, column));
}

const std::vector<std::string> MxStruct::GetFieldsName(void) const
{
	if (ptrMxStruct == nullptr)
		return ( std::vector<std::string>() );

	unsigned long nFields = mxGetNumberOfFields(ptrMxStruct);
	std::vector<std::string> fieldsName(nFields);

	for (unsigned long iField = 0; iField < nFields; ++iField)
		fieldsName[iField] = mxGetFieldNameByNumber(ptrMxStruct, iField);

	return (fieldsName);
}

const unsigned long * MxStruct::GetSize(void) const
{
	return (structSize);
}

const unsigned long MxStruct::GetLength(void) const
{
	return (GetNumElements());
}

const unsigned long MxStruct::GetNumElements(void) const
{
	return (structSize[0] * structSize[1]);
}

MxStruct & MxStruct::SetSize(const unsigned long row, const unsigned long column)
{
	mxDestroyArray(ptrMxStruct);
	structSize[0] = row;
	structSize[1] = column;
	ptrMxStruct = mxCreateStructMatrix(row, column, 0, nullptr);

	return (*this);
}

MxStruct & MxStruct::Set(const unsigned long row, const unsigned long column)
{
	return(SetSize(row, column));
}

MxStruct & MxStruct::Set(mxArray * value, const std::string & fieldName, const unsigned long row, const unsigned long column)
{
	return(Set(value, fieldName.c_str(), row, column));
}

MxStruct & MxStruct::Set(mxArray * value, const char * fieldName, const unsigned long row, const unsigned long column)
{
	if (ptrMxStruct == nullptr)
		return(*this);

	unsigned long index = sub2ind(row, column);
	mxArray * currentFieldData = mxGetField(ptrMxStruct, index, fieldName);

	if (currentFieldData == value)
		return(*this);

	mxDestroyArray(currentFieldData);

	bool isExistFieldName = mxGetFieldNumber(ptrMxStruct, fieldName) != -1;
	if (!isExistFieldName)
		if (mxAddField(ptrMxStruct, fieldName) == -1)
			throw "Can't create field name";

	mxSetField(ptrMxStruct, index, fieldName, value);
	return(*this);
}

/* ****************** Mex Cell Class ****************** */

MxCell::MxCell() : ptrMxCell(nullptr)
{
	cellSize[0] = 0;
	cellSize[1] = 0;
}

MxCell::MxCell(const unsigned long row, const unsigned long column) : ptrMxCell(nullptr)
{
	SetSize(row, column);
}

MxCell::MxCell(mxArray * initPtr) : ptrMxCell(nullptr), cellSize{ 0, 0 }
{
	if (initPtr!=nullptr && mxIsCell(initPtr)) {
		ptrMxCell = initPtr;
		cellSize[0] = static_cast<unsigned long>(mxGetM(ptrMxCell));
		cellSize[1] = static_cast<unsigned long>(mxGetN(ptrMxCell));
	};
}

MxCell::MxCell(const MxCell & src) : ptrMxCell(nullptr), cellSize{ 0, 0 }
{
	*(this) = src;
}

MxCell & MxCell::operator= (const MxCell & src)
{
	if (this == &src)
		return (*this);

	cellSize[0] = src.cellSize[0];
	cellSize[1] = src.cellSize[1];

	mxDestroyArray(ptrMxCell);
	if (src.ptrMxCell == nullptr)
		ptrMxCell = nullptr;
	else
		ptrMxCell = mxDuplicateArray(src.ptrMxCell);

	return (*this);
}

mxArray * MxCell::GetMxArray(void)
{
	return(ptrMxCell);
}

const mxArray * MxCell::GetMxArray(void) const
{
	return(ptrMxCell);
}

const unsigned long MxCell::sub2ind(const unsigned long row, const unsigned long column) const
{
	if (ptrMxCell == nullptr)
		return(0);

	mwSize subs[2];
	mwSize nsubs = 2;

	subs[0] = row;
	subs[1] = column;

	unsigned long index = mxCalcSingleSubscript(ptrMxCell, nsubs, subs);
	return (index);
}

mxArray * MxCell::operator() (const unsigned long row, const unsigned long column)
{
	if (ptrMxCell == nullptr)
		return(nullptr);
	return (mxGetCell(ptrMxCell, sub2ind(row, column)));
}

const mxArray * MxCell::operator() (const unsigned long row, const unsigned long column) const
{
	if (ptrMxCell == nullptr)
		return(nullptr);
	return (mxGetCell(ptrMxCell, sub2ind(row, column)));
}

const unsigned long * MxCell::GetSize(void) const
{
	return (cellSize);
}

const unsigned long MxCell::GetLength(void) const
{
	return (GetNumElements());
}

const unsigned long MxCell::GetNumElements(void) const
{
	return (cellSize[0] * cellSize[1]);
}

MxCell & MxCell::SetSize(const unsigned long row, const unsigned long column)
{
	mxDestroyArray(ptrMxCell);
	cellSize[0] = row;
	cellSize[1] = column;
	ptrMxCell = mxCreateCellMatrix(row, column);

	return (*this);
}

MxCell & MxCell::Set(const unsigned long row, const unsigned long column)
{
	return(SetSize(row, column));
}

MxCell & MxCell::Set(mxArray * value, const unsigned long row, const unsigned long column)
{
	if (ptrMxCell == nullptr)
		return(*this);

	unsigned long index = sub2ind(row, column);
	mxArray * currentCell = mxGetCell(ptrMxCell, index);

	if (currentCell == value)
		return (*this);

	mxDestroyArray(currentCell);
	mxSetCell(ptrMxCell, index, value);
	return (*this);
}

/* ****************** Mex Scalar Class ****************** */

MxScalar::MxScalar(double iniValue) : ptrMxScalar(nullptr)
{
	ptrMxScalar = mxCreateDoubleScalar(iniValue);
}

MxScalar::MxScalar(mxArray * iniPtr) : ptrMxScalar(nullptr)
{

	if (iniPtr!=nullptr && mxIsDouble(iniPtr) && mxGetM(iniPtr) == 1 && mxGetN(iniPtr) == 1)
		ptrMxScalar = iniPtr;
}

mxArray * MxScalar::GetMxArray(void)
{
	return (ptrMxScalar);
}

const mxArray * MxScalar::GetMxArray(void) const
{
	return (ptrMxScalar);
}

double & MxScalar::operator() (void)
{
	if (ptrMxScalar == nullptr)
		ptrMxScalar = mxCreateDoubleScalar(0);

	return(*mxGetPr(ptrMxScalar));
}

const double & MxScalar::operator() (void) const
{
	if (ptrMxScalar == nullptr)
		return(0);
	return(*mxGetPr(ptrMxScalar));
}

MxScalar & MxScalar::Set(const double value)
{
	if (ptrMxScalar == nullptr)
		ptrMxScalar = mxCreateDoubleScalar(0);
	*mxGetPr(ptrMxScalar) = value;
	return(*this);
}

/* ****************** Mex String Class ****************** */

MxString::MxString(const char * iniStr) : ptrMxString(nullptr)
{
	Set(iniStr);
}

MxString::MxString(const std::string & iniStr) : ptrMxString(nullptr)
{
	Set(iniStr);
}

MxString::MxString(mxArray * iniPtr) : ptrMxString(nullptr)
{
	if (iniPtr != nullptr && mxIsChar(iniPtr)) {
		ptrMxString = iniPtr;
		char * mxString = mxArrayToString(ptrMxString);
		stringFieldData = mxString;
		mxFree(mxString);
	};
}

mxArray * MxString::GetMxArray(void)
{
	return(ptrMxString);
}

const mxArray * MxString::GetMxArray(void) const
{
	return(ptrMxString);
}

std::string & MxString::operator() (void)
{
	if (ptrMxString == nullptr)
		ptrMxString = mxCreateString(stringFieldData.c_str());

	return (stringFieldData);
}

const char * MxString::operator() (void) const
{
	return (stringFieldData.c_str());
}

MxString & MxString::Set(const char * value)
{

	if (value != nullptr)
		stringFieldData = value;

	if (ptrMxString != nullptr) {
		char * mxString = mxArrayToString(ptrMxString);
		if (stringFieldData != mxString) {
			mxDestroyArray(ptrMxString);
			ptrMxString = nullptr;
		};
		mxFree(mxString);
	};

	if (ptrMxString == nullptr)
		ptrMxString = mxCreateString(value);

	return (*this);
}

MxString & MxString::Set(const std::string & value)
{
	return (Set(value.c_str()));
}

MxString::~MxString()
{
	if(ptrMxString==nullptr)
		return;
	char * mxString = mxArrayToString(ptrMxString);
	if (stringFieldData != mxString) {

		mxArray * mxSrc = mxCreateString(stringFieldData.c_str());
		mxChar * mxCharSrc = mxGetChars(mxSrc);

		mxChar * mxCharDst = mxGetChars(ptrMxString);
		mxFree(mxCharDst);

		mxCharDst = (mxChar *)mxCalloc(mxGetNumberOfElements(mxSrc), mxGetElementSize(mxSrc));

		for (mwSize i = 0; i < mxGetNumberOfElements(mxSrc); ++i)
			mxCharDst[i] = mxCharSrc[i];

		mxSetData(ptrMxString, mxCharDst);
		mxSetDimensions(ptrMxString, mxGetDimensions(mxSrc), mxGetNumberOfDimensions(mxSrc));

		mxDestroyArray(mxSrc);
	};

	mxFree(mxString);
}

/* ****************** Mex Matrix Class ****************** */

MxMatrix::MxMatrix() : ptrMxMatrix(nullptr)
{
	matrixSize[0] = 0;
	matrixSize[1] = 0;
}

MxMatrix::MxMatrix(const unsigned long row, const unsigned long column)
{
	matrixSize[0] = row;
	matrixSize[1] = column;
	ptrMxMatrix = mxCreateDoubleMatrix(row, column, mxREAL);
}

MxMatrix::MxMatrix(mxArray * initPtr) : ptrMxMatrix(nullptr), matrixSize{ 0, 0 }
{
	if (initPtr!=nullptr && mxIsDouble(initPtr)) {
		ptrMxMatrix = initPtr;
		matrixSize[0] = static_cast<unsigned long>(mxGetM(ptrMxMatrix));
		matrixSize[1] = static_cast<unsigned long>(mxGetN(ptrMxMatrix));
	};
}

MxMatrix::MxMatrix(const MxMatrix & src) : ptrMxMatrix(nullptr), matrixSize{ 0, 0 }
{
	(*this) = src;
}

MxMatrix & MxMatrix::operator=(const MxMatrix & src)
{
	if (this == &src)
		return (*this);

	SetSize(src.matrixSize[0], src.matrixSize[1]);
	unsigned long nElements = GetNumElements();
	double * dst = (*this)();
	const double * srcData = src();
	for (unsigned long ind = 0; ind < nElements; ++ind)
		dst[ind] = srcData[ind];
	return (*this);
}

mxArray * MxMatrix::GetMxArray(void)
{
	return (ptrMxMatrix);
}

const mxArray * MxMatrix::GetMxArray(void) const
{
	return (ptrMxMatrix);
}

const unsigned long MxMatrix::sub2ind(const unsigned long row, const unsigned long column) const
{
	if (ptrMxMatrix == nullptr)
		return(0);

	mwSize subs[2];
	mwSize nsubs = 2;

	subs[0] = row;
	subs[1] = column;

	unsigned long index = mxCalcSingleSubscript(ptrMxMatrix, nsubs, subs);
	return (index);
}

double * MxMatrix::operator() (void)
{
	if (ptrMxMatrix == nullptr)
		return(nullptr);
	return (mxGetPr(ptrMxMatrix));
}

const double * MxMatrix::operator() (void) const
{
	if (ptrMxMatrix == nullptr)
		return(nullptr);
	return (mxGetPr(ptrMxMatrix));
}

double & MxMatrix::operator() (const unsigned long row, const unsigned long column)
{
	if (ptrMxMatrix == nullptr)
		throw "Memory Access Violation. Trying to access uninitialized matrix.";

	double * ptrValue = mxGetPr(ptrMxMatrix) + sub2ind(row, column);
	return (*ptrValue);
}

const double & MxMatrix::operator() (const unsigned long row, const unsigned long column) const
{
	if (ptrMxMatrix == nullptr)
		throw "Memory Access Violation. Trying to access uninitialized matrix.";

	double * ptrValue = mxGetPr(ptrMxMatrix) + sub2ind(row, column);
	return (*ptrValue);
}

const unsigned long * MxMatrix::GetSize(void) const
{
	return (matrixSize);
}

const unsigned long MxMatrix::GetLength(void) const
{
	return(GetNumElements());
}

const unsigned long MxMatrix::GetNumElements(void) const
{
	return (matrixSize[0] * matrixSize[1]);
}

MxMatrix & MxMatrix::SetSize(const unsigned long row, const unsigned long column)
{
	matrixSize[0] = row;
	matrixSize[1] = column;

	if (ptrMxMatrix == nullptr)
		ptrMxMatrix = mxCreateDoubleMatrix(row, column, mxREAL);
	else {
		mxFree(mxGetPr(ptrMxMatrix));

		double * mxDoubleMatrix = (double *)mxCalloc(row*column, sizeof(double));
		mxSetPr(ptrMxMatrix, mxDoubleMatrix);
		mxSetM(ptrMxMatrix, matrixSize[0]);
		mxSetN(ptrMxMatrix, matrixSize[1]);
	};

	return (*this);
}

MxMatrix & MxMatrix::Set(const unsigned long row, const unsigned long column)
{
	return(SetSize(row, column));
}

MxMatrix & MxMatrix::Set(const double value, const unsigned long row, const unsigned long column)
{
	if (ptrMxMatrix == nullptr)
		return(*this);

	double * ptrValue = mxGetPr(ptrMxMatrix) + sub2ind(row, column);
	*ptrValue = value;

	return (*this);
}
