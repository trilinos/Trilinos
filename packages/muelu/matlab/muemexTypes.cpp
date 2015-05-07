#include "muemexTypes.h"

//Debug this file
//#define VERBOSE_OUTPUT

using namespace std;
using namespace Teuchos;

/* Stuff for MATLAB R2006b vs. previous versions */
#if(defined(MX_API_VER) && MX_API_VER >= 0x07030000)
#else
typedef int mwIndex;
#endif

//Flag set to true if MATLAB's CSC matrix index type is not int (usually false)
bool rewrap_ints = sizeof(int) != sizeof(mwIndex);

int* mwIndex_to_int(int N, mwIndex* mwi_array)
{
	int* rv;
	rv = new int[N];
	for(int i = 0; i < N; i++)
		rv[i] = (int) mwi_array[i];
	return rv;
}

/* ******************************* */
/* Begin MuemexData implementation */
/* ******************************* */

MuemexArg::MuemexArg(MUEMEX_TYPE dataType)
{
	type = dataType;
}

//Fully generic methods
template<typename T>
MuemexData<T>::MuemexData(T& dataToCopy, MUEMEX_TYPE dataType) : MuemexArg(dataType)
{
	data = dataToCopy;
}

template<typename T>
MuemexData<T>::~MuemexData() {}

template<typename T>
T& MuemexData<T>::getData()
{
	return data;
}

template<typename T>
void MuemexData<T>::setData(T& newData)
{
	this->data = data;
}

//string specializations
template<>
MuemexData<string>::MuemexData(const mxArray* mxa) : MuemexArg(STRING)
{
	data = "";
	if(!mxGetClassID(mxa) != mxCHAR_CLASS)
	{
		throw runtime_error("Can't construct string from anything but a char array.");
	}
	data = string(mxArrayToString(mxa));
}

template<>
mxArray* MuemexData<string>::convertToMatlab()
{
	return mxCreateString(data.c_str());
}

//int
template<>
MuemexData<int>::MuemexData(const mxArray* mxa) : MuemexArg(INT)
{
	data = parseInt(mxa);
}

template<>
mxArray* MuemexData<int>::convertToMatlab()
{
	mxArray* output = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
	int* ptr = (int*) mxGetData(output);
	*ptr = data;
	return output;
}

//double
template<>
MuemexData<double>::MuemexData(const mxArray* mxa) : MuemexArg(DOUBLE)
{
	data = *((double*) mxGetPr(mxa));
}

template<>
mxArray* MuemexData<double>::convertToMatlab()
{
	return mxCreateDoubleScalar(data);
}

//complex scalar
template<>
MuemexData<complex_t>::MuemexData(const mxArray* mxa) : MuemexArg(COMPLEX)
{
	double* realPart = mxGetPr(mxa);
	double* imagPart = mxGetPi(mxa);
	data = complex_t(*realPart, *imagPart);
}

template<>
mxArray* MuemexData<complex_t>::convertToMatlab()
{
	mxArray* output = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
	double* realPart = mxGetPr(output);
	double* imagPart = mxGetPi(output);
	*realPart = std::real<double>(data);
	*imagPart = std::imag<double>(data);
	return output;
}

//Epetra_Crs
template<>
MuemexData<RCP<Epetra_CrsMatrix>>::MuemexData(const mxArray* mxa) : MuemexArg(EPETRA_CRSMATRIX)
{
	data = epetraLoadMatrix(mxa);
}

template<>
mxArray* MuemexData<RCP<Epetra_CrsMatrix>>::convertToMatlab()
{
	return saveEpetraMatrix(data);
}

//Tpetra_Crs double
template<>
MuemexData<RCP<Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(TPETRA_MATRIX_DOUBLE)
{
	data = tpetraLoadMatrix<double>(mxa);
}

template<>
mxArray* MuemexData<RCP<Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
	return saveMatrixToMatlab<double>(MueLu::TpetraCrs_To_XpetraMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(data));
}

//Tpetra_Crs complex
template<>
MuemexData<RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(TPETRA_MATRIX_COMPLEX)
{
	data = tpetraLoadMatrix<complex_t>(mxa);
}

template<>
mxArray* MuemexData<RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
	return saveMatrixToMatlab<complex_t>(MueLu::TpetraCrs_To_XpetraMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(data));
}

//Xpetra matrix double
template<>
MuemexData<RCP<Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(XPETRA_MATRIX_DOUBLE)
{
	data = xpetraLoadMatrix<double>(mxa);
}

template<>
mxArray* MuemexData<RCP<Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
	return saveMatrixToMatlab<double>(data);
}

//Xpetra matrix complex
template<>
MuemexData<RCP<Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(XPETRA_MATRIX_COMPLEX)
{
	data = xpetraLoadMatrix<complex_t>(mxa);
}

template<>
mxArray* MuemexData<RCP<Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
	return saveMatrixToMatlab<complex_t>(data);
}

//Epetra MV
template<>
MuemexData<RCP<Epetra_MultiVector>>::MuemexData(const mxArray* mxa) : MuemexArg(EPETRA_MULTIVECTOR)
{
	data = loadEpetraMV(mxa);
}

template<>
mxArray* MuemexData<RCP<Epetra_MultiVector>>::convertToMatlab()
{
	return saveEpetraMV(data);
}

//Tpetra MV double
template<>
MuemexData<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(TPETRA_MULTIVECTOR_DOUBLE)
{
	data = loadTpetraMV<double>(mxa);
}

template<>
mxArray* MuemexData<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
	return saveTpetraMV<double>(data);
}

//Tpetra MV complex
template<>
MuemexData<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(TPETRA_MULTIVECTOR_COMPLEX)
{
	data = loadTpetraMV<complex_t>(mxa);
}

template<>
mxArray* MuemexData<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
	return saveTpetraMV<complex_t>(data);
}

//Xpetra ordinal vector
template<>
MuemexData<RCP<Xpetra_ordinal_vector>>::MuemexData(const mxArray* mxa) : MuemexArg(XPETRA_ORDINAL_VECTOR)
{
	data = loadLOVector(mxa);
}

template<>
mxArray* MuemexData<RCP<Xpetra_ordinal_vector>>::convertToMatlab()
{
	return createMatlabLOVector(data);
}

//Xpetra multivector double
template<>
MuemexData<RCP<Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(XPETRA_MULTIVECTOR_DOUBLE)
{
	data = loadXpetraMV<double>(mxa);
}

template<>
mxArray* MuemexData<RCP<Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
	return saveMultiVectorToMatlab<double>(data);
}

//Xpetra multivector complex
template<>
MuemexData<RCP<Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(XPETRA_MULTIVECTOR_COMPLEX)
{
	data = loadXpetraMV<complex_t>(mxa);
}

template<>
mxArray* MuemexData<RCP<Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
	return saveMultiVectorToMatlab<complex_t>(data);
}

/* ***************************** */
/* End MuemexData implementation */
/* ***************************** */

template<typename Scalar = double>
mxArray* saveMatrixToMatlab(RCP<Xpetra::Matrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mat)
{
	int nr = mat->getGlobalNumRows();
	int nc = mat->getGlobalNumCols();
	int nnz = mat->getGlobalNumEntries();
	#ifdef VERBOSE_OUTPUT
	RCP<Teuchos::FancyOStream> fancyStream = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
	mat->describe(*fancyStream, Teuchos::VERB_EXTREME);
	#endif
	mxArray* mxa = createMatlabSparse<Scalar>(nr, nc, nnz);
	mwIndex* ir = mxGetIr(mxa);
	mwIndex* jc = mxGetJc(mxa);
	for(int i = 0; i < nc + 1; i++)
	{
		jc[i] = 0;
	}
	size_t maxEntriesPerRow = mat->getGlobalMaxNumRowEntries();
	int* rowProgress = new int[nc];
	//The array that will be copied to Pr and (if complex) Pi later
	Scalar* sparseVals = new Scalar[nnz];
	size_t numEntries;
	if(mat->isLocallyIndexed())
	{
		Scalar* rowValArray = new Scalar[maxEntriesPerRow];
		ArrayView<Scalar> rowVals(rowValArray, maxEntriesPerRow);
		mm_LocalOrd* rowIndicesArray = new mm_LocalOrd[maxEntriesPerRow];
		ArrayView<mm_LocalOrd> rowIndices(rowIndicesArray, maxEntriesPerRow);
		for(mm_LocalOrd m = 0; m < nr; m++)	//All rows in the Xpetra matrix
		{
			mat->getLocalRowCopy(m, rowIndices, rowVals, numEntries);	//Get the row
			for(mm_LocalOrd entry = 0; entry < int(numEntries); entry++)	//All entries in row
			{
				jc[rowIndices[entry] + 1]++; //for each entry, increase jc for the entry's column
			}
		}
		//now jc holds the number of elements in each column, but needs cumulative sum over all previous columns also
		int entriesAccum = 0;
		for(int n = 0; n <= nc; n++)
		{
			int temp = entriesAccum;
			entriesAccum += jc[n];
			jc[n] += temp;
		}
		//Jc now populated with colptrs
		for(int i = 0; i < nc; i++)
		{
			rowProgress[i] = 0;
		}
		//Row progress values like jc but keep track as the MATLAB matrix is being filled in
		for(mm_LocalOrd m = 0; m < nr; m++)	//rows
		{
			mat->getLocalRowCopy(m, rowIndices, rowVals, numEntries);
			for(mm_LocalOrd i = 0; i < int(numEntries); i++)	//entries in row m (NOT columns)
			{
				//row is m, col is rowIndices[i], val is rowVals[i]
				mm_LocalOrd col = rowIndices[i];
				sparseVals[jc[col] + rowProgress[col]] = rowVals[i];	//Set value
				ir[jc[col] + rowProgress[col]] = m;						//Set row at which value occurs
				rowProgress[col]++;
			}
		}
		delete[] rowIndicesArray;
	}
	else
	{
		ArrayView<const mm_GlobalOrd> rowIndices;
		ArrayView<const Scalar> rowVals;
		for(mm_GlobalOrd m = 0; m < nr; m++)
		{
			mat->getGlobalRowView(m, rowIndices, rowVals);
			for(mm_GlobalOrd n = 0; n < rowIndices.size(); n++)
			{
				jc[rowIndices[n] + 1]++;
			}
		}
		//Last element of jc is just nnz
		jc[nc] = nnz;
		//Jc now populated with colptrs
		for(int i = 0; i < nc; i++)
		{
			rowProgress[i] = 0;
		}
		int entriesAccum = 0;
		for(int n = 0; n <= nc; n++)
		{
			int temp = entriesAccum;
			entriesAccum += jc[n];
			jc[n] += temp;
		}
		//Row progress values like jc but keep track as the MATLAB matrix is being filled in
		for(mm_GlobalOrd m = 0; m < nr; m++)	//rows
		{
			mat->getGlobalRowView(m, rowIndices, rowVals);
			for(mm_LocalOrd i = 0; i < rowIndices.size(); i++)	//entries in row m (NOT == columns)
			{
				//row is m, col is rowIndices[i], val is rowVals[i]
				mm_GlobalOrd col = rowIndices[i];
				sparseVals[jc[col] + rowProgress[col]] = rowVals[i];	//Set value
				ir[jc[col] + rowProgress[col]] = m;						//Set row at which value occurs
				rowProgress[col]++;
			}
		}
	}
	//finally, copy sparseVals into pr (and pi, if complex)
	fillMatlabArray<Scalar>(sparseVals, mxa, nnz);
	delete[] sparseVals;
	delete[] rowProgress;
	return mxa;
}

template<> mxArray* createMatlabSparse<double>(int numRows, int numCols, int nnz)
{
	return mxCreateSparse(numRows, numCols, nnz, mxREAL);
}

template<> mxArray* createMatlabSparse<complex_t>(int numRows, int numCols, int nnz)
{
	return mxCreateSparse(numRows, numCols, nnz, mxCOMPLEX);
}

template<> void fillMatlabArray<double>(double* array, const mxArray* mxa, int n)
{
	memcpy(mxGetPr(mxa), array, n * sizeof(double));
}

template<> void fillMatlabArray<complex_t>(complex_t* array, const mxArray* mxa, int n)
{
	double* pr = mxGetPr(mxa);
	double* pi = mxGetPi(mxa);
	for(int i = 0; i < n; i++)
	{
		pr[i] = real<double>(array[i]);
		pi[i] = imag<double>(array[i]);
	}
}

RCP<Epetra_CrsMatrix> epetraLoadMatrix(const mxArray* mxa)
{
	RCP<Epetra_CrsMatrix> matrix;
	try
	{
		int* colptr;
		int* rowind;
		double* vals = mxGetPr(mxa);
		int nr = mxGetM(mxa);
		int nc = mxGetN(mxa);
		if(rewrap_ints)
		{
			colptr = mwIndex_to_int(nc + 1, mxGetJc(mxa));
			rowind = mwIndex_to_int(colptr[nc], mxGetIr(mxa));
		}
		else
		{
			rowind = (int*) mxGetIr(mxa);
			colptr = (int*) mxGetJc(mxa);
		}
		Epetra_SerialComm Comm;
		Epetra_Map RangeMap(nr, 0, Comm);
		Epetra_Map DomainMap(nc, 0, Comm);
		matrix = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, RangeMap, DomainMap, 0));
		/* Do the matrix assembly */
		for(int i = 0; i < nc; i++)
		{
			for(int j = colptr[i]; j < colptr[i + 1]; j++)
			{
				//		global row, # of entries, value array, column indices array
				matrix->InsertGlobalValues(rowind[j], 1, &vals[j], &i);
			}
		}
		matrix->FillComplete(DomainMap, RangeMap);
		if(rewrap_ints)
		{
			delete [] rowind;
			delete [] colptr;
		}
	}
	catch(exception& e)
	{
		mexPrintf("An error occurred while setting up an Epetra matrix:\n");
		cout << e.what() << endl;
	}
	return matrix;
}

template<>
RCP<Tpetra_CrsMatrix_double> tpetraLoadMatrix<double>(const mxArray* mxa)
{
	bool success = false;
	RCP<Tpetra_CrsMatrix_double> A;
	try
	{
		RCP<const Teuchos::Comm<int>> comm = rcp(new Teuchos::SerialComm<int>());
		//numGlobalIndices is just the number of rows in the matrix	
		const Tpetra::global_size_t numGlobalIndices = mxGetM(mxa);
		const mm_GlobalOrd indexBase = 0;
		RCP<const muemex_map_type> rowMap = rcp(new muemex_map_type(numGlobalIndices, indexBase, comm));
		RCP<const muemex_map_type> domainMap = rcp(new muemex_map_type(mxGetN(mxa), indexBase, comm));
		A = Tpetra::createCrsMatrix<double, mm_GlobalOrd, mm_LocalOrd, mm_node_t>(rowMap);
		double* valueArray = mxGetPr(mxa);
		int* colptr;
		int* rowind;
		//int nr = mxGetM(mxa);
		int nc = mxGetN(mxa);
		if(rewrap_ints)
		{
			//mwIndex_to_int allocates memory so must delete[] later
			colptr = mwIndex_to_int(nc + 1, mxGetJc(mxa));
			rowind = mwIndex_to_int(colptr[nc], mxGetIr(mxa));
		}
		else
		{
			rowind = (int*) mxGetIr(mxa);
			colptr = (int*) mxGetJc(mxa);
		}
		for(int i = 0; i < nc; i++)
		{
			for(int j = colptr[i]; j < colptr[i + 1]; j++)
			{
				//'array' of 1 element, containing column (in global matrix).
				ArrayView<mm_GlobalOrd> cols = ArrayView<mm_GlobalOrd>(&i, 1);
				//'array' of 1 element, containing value
				ArrayView<double> vals = ArrayView<double>(&valueArray[j], 1);
				A->insertGlobalValues(rowind[j], cols, vals);
			}
		}
		A->fillComplete(domainMap, rowMap);
		if(rewrap_ints)
		{
			delete[] rowind;
			delete[] colptr;
		}
		success = true;
	}
	catch(exception& e)
	{
		mexPrintf("Error while constructing Tpetra matrix:\n");
		cout << e.what() << endl;
	}
	if(!success)
		mexErrMsgTxt("An error occurred while setting up a Tpetra matrix.\n");
	return A;
}

template<>
RCP<Tpetra_CrsMatrix_complex> tpetraLoadMatrix<complex_t>(const mxArray* mxa)
{
	RCP<Tpetra_CrsMatrix_complex> A;
	//Create a map in order to create the matrix (taken from muelu basic example - complex)
	try
	{
		RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
		const Tpetra::global_size_t numGlobalIndices = mxGetM(mxa);
		const mm_GlobalOrd indexBase = 0;
		RCP<const muemex_map_type> map = rcp(new muemex_map_type(numGlobalIndices, indexBase, comm));
		A = rcp(new Tpetra_CrsMatrix_complex(map, 0));
		double* realArray = mxGetPr(mxa);
		double* imagArray = mxGetPi(mxa);
		int* colptr;
		int* rowind;
		int nc = mxGetN(mxa);
		if(rewrap_ints)
		{
			//mwIndex_to_int allocates memory so must delete[] later
			colptr = mwIndex_to_int(nc + 1, mxGetJc(mxa));
			rowind = mwIndex_to_int(colptr[nc], mxGetIr(mxa));
		}
		else
		{
			rowind = (int*) mxGetIr(mxa);
			colptr = (int*) mxGetJc(mxa);
		}
		for(int i = 0; i < nc; i++)
		{
			for(int j = colptr[i]; j < colptr[i + 1]; j++)
			{
				//here assuming that complex_t will always be defined as std::complex<double>
				//use 'value' over and over again with ArrayViews to insert into matrix
				complex_t value = std::complex<double>(realArray[j], imagArray[j]);
				ArrayView<mm_GlobalOrd> cols = ArrayView<mm_GlobalOrd>(&i, 1);
				ArrayView<complex_t> vals = ArrayView<complex_t>(&value, 1);
				A->insertGlobalValues(rowind[j], cols, vals);
			}
		}
		A->fillComplete();
		if(rewrap_ints)
		{
			delete[] rowind;
			delete[] colptr;
		}
	}
	catch(exception& e)
	{
		mexPrintf("Error while constructing tpetra matrix:\n");
		cout << e.what() << endl;
	}
	return A;
}

mxArray* createMatlabLOVector(RCP<Xpetra_ordinal_vector> vec)
{
	//this value might be a 64 bit int but it should never overflow a 32
	mwSize len = vec->getGlobalLength();
	//create a single column vector
	mwSize dimensions[] = {len, 1};
	return mxCreateNumericArray(2, dimensions, mxINT32_CLASS, mxREAL);
}

template<>
RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loadTpetraMV<double>(const mxArray* mxa)
{
	RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv;
	try
	{
		int nr = mxGetM(mxa);
		int nc = mxGetN(mxa);
		double* pr = mxGetPr(mxa);
		RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
		//numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
		RCP<const muemex_map_type> map = rcp(new muemex_map_type(nr, (mm_GlobalOrd) 0, comm));
		//Allocate a new array of complex values to use with the multivector
		ArrayView<const double> arrView(pr, nr * nc);
		mv = rcp(new Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(map, arrView, size_t(nr), size_t(nc)));
	}
	catch(exception& e)
	{
		mexPrintf("Error constructing Tpetra MultiVector.\n");
		cout << e.what() << endl;
	}
	return mv;
}

template<>
RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loadTpetraMV<complex_t>(const mxArray* mxa)
{
	RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv;
	try
	{
		int nr = mxGetM(mxa);
		int nc = mxGetN(mxa);
		double* pr = mxGetPr(mxa);
		double* pi = mxGetPi(mxa);
		RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
		//numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
		RCP<const muemex_map_type> map = rcp(new muemex_map_type(nr, (mm_GlobalOrd) 0, comm));
		//Allocate a new array of complex values to use with the multivector
		complex_t* myArr = new complex_t[nr * nc];
		for(int n = 0; n < nc; n++)
		{
			for(int m = 0; m < nr; m++)
			{
				myArr[n * nr + m] = complex_t(pr[n * nr + m], pi[n * nr + m]);
			}
		}
		ArrayView<complex_t> arrView(myArr, nr * nc);
		mv = rcp(new Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(map, arrView, nr, nc));
	}
	catch(exception& e)
	{
		mexPrintf("Error constructing Tpetra MultiVector.\n");
		cout << e.what() << endl;
	}
	return mv;
}

RCP<Epetra_MultiVector> loadEpetraMV(const mxArray* mxa)
{
	int nr = mxGetM(mxa);
	int nc = mxGetN(mxa);
	Epetra_SerialComm Comm;
	Epetra_BlockMap map(nr * nc, 1, 0, Comm);
	return rcp(new Epetra_MultiVector(Epetra_DataAccess::Copy, map, mxGetPr(mxa), nr, nc));
}

template<typename Scalar>
RCP<Xpetra::Matrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> xpetraLoadMatrix(const mxArray* mxa)
{
	RCP<Tpetra::CrsMatrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> tpetraMat = tpetraLoadMatrix<Scalar>(mxa);
	return MueLu::TpetraCrs_To_XpetraMatrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(tpetraMat);
}

template<> mxArray* createMatlabMultiVector<double>(int numRows, int numCols)
{
	return mxCreateDoubleMatrix(numRows, numCols, mxREAL);
}

template<> mxArray* createMatlabMultiVector<complex_t>(int numRows, int numCols)
{
	return mxCreateDoubleMatrix(numRows, numCols, mxCOMPLEX);
}

template<typename Scalar>
mxArray* saveTpetraMV(RCP<Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv)
{
	//Precondition: Memory has already been allocated by MATLAB for the array.
	int nr = mv->getGlobalLength();
	int nc = mv->getNumVectors();
	mxArray* output = createMatlabMultiVector<Scalar>(nr, nc);
	Scalar* data = new Scalar[nr * nc];
	for(int col = 0; col < nc; col++)
	{
		ArrayRCP<const Scalar> colData = mv->getData(col);
		for(int row = 0; row < nr; row++)
		{
			data[col * nr + row] = colData[row];
		}
	}
	fillMatlabArray<Scalar>(data, output, nc * nr);
	return output;
}

mxArray* saveEpetraMV(RCP<Epetra_MultiVector> mv)
{
	mxArray* output = mxCreateDoubleMatrix(mv->GlobalLength(), mv->NumVectors(), mxREAL);
	double* dataPtr = mxGetPr(output);
	mv->ExtractCopy(dataPtr, mv->GlobalLength());
	return output;
}

template<typename Scalar>
mxArray* saveMultiVectorToMatlab(RCP<Xpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv)
{
	//Precondition: Memory has already been allocated by MATLAB for the array.
	int nr = mv->getGlobalLength();
	int nc = mv->getNumVectors();
	mxArray* output = createMatlabMultiVector<Scalar>(nr, nc);
	Scalar* data = new Scalar[nr * nc];
	for(int col = 0; col < nc; col++)
	{
		ArrayRCP<const Scalar> colData = mv->getData(col);
		for(int row = 0; row < nr; row++)
		{
			data[col * nr + row] = colData[row];
		}
	}
	fillMatlabArray<Scalar>(data, output, nc * nr);
	return output;
}

int parseInt(const mxArray* mxa)
{
	mxClassID probIDtype = mxGetClassID(mxa);
	int rv;
	if(probIDtype == mxINT32_CLASS)
	{
		rv = *((int*) mxGetData(mxa));
	}
	else if(probIDtype == mxDOUBLE_CLASS)
	{
		rv = (int) *((double*) mxGetData(mxa));
	}
	else if(probIDtype == mxUINT32_CLASS)
	{
		rv = (int) *((unsigned int*) mxGetData(mxa));
	}
	else
	{
		rv = -1;
		throw runtime_error("Error: Unrecognized numerical type.");
	}
	return rv;
}

mxArray* saveEpetraMatrix(RCP<Epetra_CrsMatrix> mat)
{
	return saveMatrixToMatlab<double>(MueLu::EpetraCrs_To_XpetraMatrix(mat));
}

RCP<Xpetra_ordinal_vector> loadLOVector(const mxArray* mxa)
{
	RCP<const Teuchos::Comm<int>> comm = rcp(new Teuchos::SerialComm<int>());
	const Tpetra::global_size_t numGlobalIndices = mxGetM(mxa);
	RCP<const muemex_map_type> rowMap = rcp(new muemex_map_type(numGlobalIndices, (mm_GlobalOrd) 0, comm));
	if(mxGetClassID(mxa) != mxINT32_CLASS || mxGetN(mxa) != 1)
		throw runtime_error("Can only construct LOVector with int32 single vector.");
	int* array = (int*) mxGetData(mxa);
	ArrayView<int> dataView(array, mxGetM(mxa));
	RCP<Tpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loVec = rcp(new Tpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(rowMap, dataView));
	return Xpetra::toXpetra<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(loVec);
}

template<typename Scalar>
RCP<Xpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loadXpetraMV(const mxArray* mxa);
{
	RCP<Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> tmv = loadTpetraMV<Scalar>(mxa);
	return MueLu::TpetraMultiVector_To_XpetraMultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(tmv);
}
