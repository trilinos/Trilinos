// @HEADER
//
// ***********************************************************************
//
//		  MueLu: A package for multigrid based preconditioning
//					Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//					  Jonathan Hu		(jhu@sandia.gov)
//					  Andrey Prokopenko (aprokop@sandia.gov)
//					  Ray Tuminaro		(rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include "muemex.h"
#define IS_FALSE 0
#define IS_TRUE 1
#define MUEMEX_ERROR -1

//Do not compile MueMex if any of these aren't available
#if !defined HAVE_MUELU_EPETRA || !defined HAVE_MUELU_TPETRA || !defined HAVE_MUELU_MATLAB
#error "MueMex requires Epetra, Tpetra and MATLAB."
#endif

using namespace std;
using namespace Teuchos;

extern void _main();

/* MUEMEX Teuchos Parameters*/
#define MUEMEX_INTERFACE "Linear Algebra"

/* Default values */
#define MUEMEX_DEFAULT_LEVELS 10
#define MUEMEX_DEFAULT_NUMPDES 1
#define MUEMEX_DEFAULT_ADAPTIVEVECS 0
#define MUEMEX_DEFAULT_USEDEFAULTNS true
#define MMABS(x)   ((x)>0?(x):(-(x)))
#define MMISINT(x) ((x)==0?(((x-(int)(x))<1e-15)?true:false):(((x-(int)(x))<1e-15*MMABS(x))?true:false))

/* Debugging */
//#define VERBOSE_OUTPUT

/* Stuff for MATLAB R2006b vs. previous versions */
#if(defined(MX_API_VER) && MX_API_VER >= 0x07030000)
#else
typedef int mwIndex;
#endif

//Declare and call default constructor for data_pack_list vector (starts empty)
vector<RCP<MuemexSystem>> MuemexSystemList::list;
int MuemexSystemList::nextID = 0;

//Need a global flag to keep track of Epetra vs. Tpetra for constructing multivectors for param lists
bool useEpetra = false;
//Flag set to true if MATLAB's CSC matrix index type is not int (usually false)
bool rewrap_ints = false;

/* Epetra utility functions */

/* mwIndex_to_int - does a data copy and wraps mwIndex's to ints, in the case
   where they're not the same size.	 This routine allocates memory
   WARNING: This does not address overflow.
   Parameters:
   N		 - Number of unknowns in array [I]
   mwi_array - Array of mwIndex objects [I]
   Return value: mwIndex objects cast down to ints
*/

int* mwIndex_to_int(int N, mwIndex* mwi_array)
{
	int* rv;
	rv = new int[N];
	for(int i = 0; i < N; i++)
		rv[i] = (int) mwi_array[i];
	return rv;
}

//Parse a string to get each bit of Belos verbosity
int strToMsgType(const char* str)
{
	if(str == NULL)
		return Belos::Errors;
	else if(strstr(str, "Warnings") != NULL)
		return Belos::Warnings;
	else if(strstr(str, "IterationDetails") != NULL)
		return Belos::IterationDetails;
	else if(strstr(str, "OrthoDetails") != NULL)
		return Belos::OrthoDetails;
	else if(strstr(str, "FinalSummary") != NULL)
		return Belos::FinalSummary;
	else if(strstr(str, "TimingDetails") != NULL)
		return Belos::TimingDetails;
	else if(strstr(str, "StatusTestDetails") != NULL)
		return Belos::StatusTestDetails;
	else if(strstr(str, "Debug") != NULL)
		return Belos::Debug;
	//This has no effect when added/OR'd with flags
	return Belos::Errors;
}

HierAttribType strToHierAttribType(const char* str)
{
	if(strcmp(str, "A") == 0)
		return MATRIX;
	if(strcmp(str, "P") == 0)
		return MATRIX;
	if(strcmp(str, "R") == 0)
		return MATRIX;
	if(strcmp(str, "Nullspace") == 0)
		return MULTIVECTOR;
	if(strcmp(str, "Coordinates") == 0)
		return MULTIVECTOR;
	//TODO: Add eigenvalue scalar type here - what's it called?
	if(strcmp(str, "Aggregates") == 0)
		return LOVECTOR;
	return UNKNOWN;
}

//Parse a string to get Belos output style (Brief is default)
int strToOutputStyle(const char* str)
{
	if(strstr("General", str) != NULL)
		return Belos::General;
	else
		return Belos::Brief;
}

//Get Belos "Verbosity" setting for its ParameterList
int getBelosVerbosity(const char* input)
{
	int result = 0;
	char* str = (char*) input;
	char* pch;
	pch = strtok(str, " +,");
	if(pch == NULL)
		return result;
	result |= strToMsgType(pch);
	while(pch != NULL)
	{
		pch = strtok(NULL, " +,");
		if(pch == NULL)
			return result;
		result |= strToMsgType(pch);
	}
	return result;
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
#ifdef HAVE_COMPLEX_SCALARS
template<> mxArray* createMatlabSparse<complex_t>(int numRows, int numCols, int nnz)
{
	return mxCreateSparse(numRows, numCols, nnz, mxCOMPLEX);
}
#endif
template<> void fillMatlabArray<double>(double* array, const mxArray* mxa, int n)
{
	memcpy(mxGetPr(mxa), array, n * sizeof(double));
}
#ifdef HAVE_COMPLEX_SCALARS
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
#endif
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
#ifdef HAVE_COMPLEX_SCALARS
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
#endif
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
#ifdef HAVE_COMPLEX_SCALARS
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
#endif
RCP<Epetra_MultiVector> loadEpetraMV(const mxArray* mxa)
{
	int nr = mxGetM(mxa);
	int nc = mxGetN(mxa);
	Epetra_SerialComm Comm;
	Epetra_BlockMap map(nr * nc, 1, 0, Comm);
	return rcp(new Epetra_MultiVector(Epetra_DataAccess::Copy, map, mxGetPr(mxa), nr, nc));
}

template<typename Scalar>
mxArray* TpetraSystem<Scalar>::solve(RCP<ParameterList> params, RCP<Tpetra::CrsMatrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> matrix, const mxArray* b, int& iters)
{
	mxArray* output;
	try
	{
		int matSize = A->getGlobalNumRows();
		//Define Tpetra vector/multivector types for convenience
		typedef Tpetra::Vector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Vector;
		typedef Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_MultiVector;
		typedef Tpetra::Operator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Operator;
		RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
		//numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
		RCP<const muemex_map_type> map = rcp(new muemex_map_type(matSize, (mm_GlobalOrd) 0, comm));
		RCP<Tpetra_MultiVector> rhs = loadTpetraMV<Scalar>(b);
		RCP<Tpetra_MultiVector> lhs = rcp(new Tpetra_MultiVector(map, rhs->getNumVectors()));
		//rhs is initialized, lhs is not
		iters = 0;
		#ifdef VERBOSE_OUTPUT
		params->get("Verbosity", Belos::Errors | Belos::Warnings | Belos::Debug | Belos::FinalSummary | Belos::IterationDetails | Belos::OrthoDetails | Belos::TimingDetails | Belos::StatusTestDetails);
		params->get("Output Frequency", 1);
		params->get("Output Style", Belos::Brief);
		#else
		params->get("Verbosity", Belos::Errors + Belos::Warnings);
		#endif
		RCP<Belos::LinearProblem<Scalar, Tpetra_MultiVector, Tpetra_Operator>> problem = rcp(new Belos::LinearProblem<Scalar, Tpetra_MultiVector, Tpetra_Operator>(matrix, lhs, rhs));
		problem->setRightPrec(prec);
		bool set = problem->setProblem();
		TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error, "Linear Problem failed to set up correctly!");
		Belos::SolverFactory<Scalar, Tpetra_MultiVector, Tpetra_Operator> factory;
		string solverName = params->get("solver", "GMRES");
		RCP<Belos::SolverManager<Scalar, Tpetra_MultiVector, Tpetra_Operator>> solver = factory.create(solverName, params);
		solver->setProblem(problem);
		Belos::ReturnType ret = solver->solve();
		if(ret == Belos::Converged)
		{
			mexPrintf("Success, Belos converged!\n");
			iters = solver->getNumIters();
			output = saveTpetraMV<Scalar>(lhs);
		}
		else
		{
			mexPrintf("Belos failed to converge.\n");
			iters = 0;
			output = mxCreateDoubleScalar(0);
		}
	}
	catch(exception& e)
	{
		mexPrintf("Error occurred while running Belos solver:\n");
		cout << e.what() << endl;
		output = mxCreateDoubleScalar(0);
	}
	return output;
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

template<> RCP<Hierarchy_double> getDatapackHierarchy<double>(MuemexSystem* dp)
{
	RCP<MueLu::Hierarchy<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> hier;
	switch(dp->type)
	{
		case EPETRA:
		{
			EpetraSystem* pack = (EpetraSystem*) dp;
			hier = pack->getHierarchy();
			break;
		}
		case TPETRA:
		{
			TpetraSystem<double>* pack = (TpetraSystem<double>*) dp;
			hier = pack->getHierarchy();
			break;
		}
		default:
		{
			throw runtime_error("Got unexpected linear system type for real-valued functions.");
		}
	}
	return hier;
}

template<> RCP<Hierarchy_complex> getDatapackHierarchy<complex_t>(MuemexSystem* dp)
{
	return ((TpetraSystem<complex_t>*) dp)->getHierarchy();
}

//data pack base class implementation

MuemexSystem::MuemexSystem(DataPackType probType) : id(MUEMEX_ERROR), type(probType) {}
MuemexSystem::~MuemexSystem() {}

mxArray* MuemexSystem::getHierarchyData(string dataName, HierAttribType dataType, int levelID)
{
	mxArray* output = NULL;
	try
	{
		switch(dataType)
		{
			case MATRIX:
			{
				switch(this->type)	//datapack type (EPETRA TPETRA or TPETRA_COMPLEX)
				{
					//get real matrix, put into output
					case EPETRA:
					case TPETRA:
					{
						RCP<Hierarchy_double> hier = getDatapackHierarchy<double>(this);
						RCP<MueLu::Level> level = hier->GetLevel(levelID);
						RCP<Xpetra_Matrix_double> mat;
						level->Get(dataName, mat);
						output = saveMatrixToMatlab<double>(mat);
						break;
					}
					case TPETRA_COMPLEX:
					{
						RCP<Hierarchy_complex> hier = getDatapackHierarchy<complex_t>(this);
						RCP<MueLu::Level> level = hier->GetLevel(levelID);
						RCP<Xpetra_Matrix_complex> mat;
						level->Get(dataName, mat);
						output = saveMatrixToMatlab<complex_t>(mat);
						break;
					}
				}
				break;
			}
			case MULTIVECTOR:
			{
				switch(this->type)
				{
					case EPETRA:
					case TPETRA:
					{
						RCP<Hierarchy_double> hier = getDatapackHierarchy<double>(this);
						RCP<MueLu::Level> level = hier->GetLevel(levelID);
						RCP<Xpetra_MultiVector_double> mv;
						level->Get(dataName, mv);
						output = saveMultiVectorToMatlab<double>(mv);
						break;
					}
					case TPETRA_COMPLEX:
					{
						RCP<Hierarchy_complex> hier = getDatapackHierarchy<complex_t>(this);
						RCP<MueLu::Level> level = hier->GetLevel(levelID);
						RCP<Xpetra_MultiVector_complex> mv;
						level->Get(dataName, mv);
						output = saveMultiVectorToMatlab<complex_t>(mv);
						break;
					}
				}
				break;
			}
			case LOVECTOR:
			{
				RCP<Xpetra_ordinal_vector> loVec;
				switch(this->type)
				{
					case EPETRA:
					case TPETRA:
					{
						RCP<Hierarchy_double> hier = getDatapackHierarchy<double>(this);
						RCP<MueLu::Level> level = hier->GetLevel(levelID);
						#ifdef VERBOSE_OUTPUT
						level->print(cout, MueLu::MsgType::Extreme);
						#endif
						level->Get(dataName, loVec);
						break;
					}
					case TPETRA_COMPLEX:
					{
						RCP<Hierarchy_complex> hier = getDatapackHierarchy<complex_t>(this);
						RCP<MueLu::Level> level = hier->GetLevel(levelID);
						level->Get(dataName, loVec);
						break;
					}
				}
				output = createMatlabLOVector(loVec);
				break;
			}
			case SCALAR:
			{
				switch(this->type)
				{
					case EPETRA:
					case TPETRA:
					{
						double value;
						RCP<Hierarchy_double> hier = getDatapackHierarchy<double>(this);
						RCP<MueLu::Level> level = hier->GetLevel(levelID);
						level->Get(dataName, value);
						output = mxCreateDoubleScalar(value);
						break;
					}
					case TPETRA_COMPLEX:
					{
						complex_t value;
						RCP<Hierarchy_complex> hier = getDatapackHierarchy<complex_t>(this);
						RCP<MueLu::Level> level = hier->GetLevel(levelID);
						level->Get(dataName, value);
						output = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
						double* realPart = mxGetPr(output);
						*realPart = real<double>(value);
						double* imagPart = mxGetPi(output);
						*imagPart = imag<double>(value);
						break;
					}
				}
				break;
			}
			default:
			{
				throw runtime_error("getHierarchyData can't determine the type of the data requested.");
			}
		}
		if(output == NULL)
		{
			throw runtime_error("mxArray pointer was never initialized. Check data type and name.");
		}
	}
	catch(exception& e)
	{
		mexPrintf("Error occurred while getting hierarchy data.\n");
		cout << e.what() << endl;
	}
	return output;
}

//EpetraSystem impl

EpetraSystem::EpetraSystem() : MuemexSystem(EPETRA) {}
EpetraSystem::~EpetraSystem() {}

int EpetraSystem::status()
{
	mexPrintf("**** Problem ID %d [MueLu_Epetra] ****\n", id);
	if(!A.is_null())
		mexPrintf("Matrix: %dx%d w/ %d nnz\n", A->NumGlobalRows(), A->NumGlobalCols(), A->NumMyNonzeros());
	mexPrintf("Operator Complexity: %f\n", operatorComplexity);
	if(!List.is_null())
	{
		mexPrintf("Parameter List:\n");
		List->print();
	}
	mexPrintf("\n");
	return IS_TRUE;
}/*end status*/

int EpetraSystem::setup(const mxArray* mxa)
{
	bool success = false;
	try
	{
		/* Matrix Fill */
		A = epetraLoadMatrix(mxa);
		prec = MueLu::CreateEpetraPreconditioner(A, *List);
		//underlying the Epetra_Operator prec is a MueLu::EpetraOperator
		RCP<MueLu::EpetraOperator> meo = rcp_static_cast<MueLu::EpetraOperator, Epetra_Operator>(prec);
		operatorComplexity = meo->GetHierarchy()->GetOperatorComplexity();
		success = true;
	}
	catch(exception& e)
	{
		mexPrintf("Error occurred while setting up epetra problem:\n");
		cout << e.what() << endl;
	}
	return success ? IS_TRUE : IS_FALSE;
}/*end setup*/

/* EpetraSystem::solve - Given two Teuchos lists, one in the EpetraSystem, and one of
   solve-time options, this routine calls the relevant solver and returns the solution.
   TPL	   - Teuchos list of solve-time options [I]
   A	   - The matrix to solve with (may not be the one the preconditioned was used for)
   b	   - RHS vector [I]
   x	   - solution vector [O]
   iters   - number of iterations taken [O]
   Returns: IS_TRUE if solve was succesful, IS_FALSE otherwise
*/
mxArray* EpetraSystem::solve(RCP<ParameterList> TPL, RCP<Epetra_CrsMatrix> matrix, const mxArray* b, int& iters)
{
	mxArray* output;
	try
	{
		//Set up X and B
		Epetra_Map map = matrix->DomainMap();
		RCP<Epetra_MultiVector> rhs = loadEpetraMV(b);
		RCP<Epetra_MultiVector> lhs = rcp(new Epetra_MultiVector(map, rhs->NumVectors(), true));
		#ifdef VERBOSE_OUTPUT
		TPL->get("Verbosity", Belos::Errors | Belos::Warnings | Belos::Debug | Belos::FinalSummary | Belos::IterationDetails | Belos::OrthoDetails | Belos::TimingDetails | Belos::StatusTestDetails);
		TPL->get("Output Frequency", 1);
		TPL->get("Output Style", Belos::Brief);
		#else
		TPL->get("Verbosity", Belos::Errors | Belos::Warnings);
		#endif
		RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>> problem = rcp(new    Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>(matrix, lhs, rhs));
		RCP<Belos::EpetraPrecOp> epo = rcp(new Belos::EpetraPrecOp(prec));
		problem->setRightPrec(epo);
		bool set = problem->setProblem();
		TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error, "Linear Problem failed to set up correctly!");
		Belos::SolverFactory<double, Epetra_MultiVector, Epetra_Operator> factory;
		//Get the solver name from the parameter list, default to PseudoBlockGmres if none specified by user
		string solverName = TPL->get("solver", "GMRES");
		RCP<Belos::SolverManager<double, Epetra_MultiVector, Epetra_Operator>> solver = factory.create(solverName, TPL);
		solver->setProblem(problem);
		Belos::ReturnType ret = solver->solve();
		if(ret == Belos::Converged)
		{
			mexPrintf("Success, Belos converged!\n");
			iters = solver->getNumIters();
			output = saveEpetraMV(lhs);
		}
		else
		{
			mexPrintf("Belos failed to converge.\n");
			iters = 0;
			output = mxCreateDoubleScalar(0);
		}
		output = saveEpetraMV(lhs);
	}
	catch(exception& e)
	{
		mexPrintf("Error occurred during Belos solve:\n");
		cout << e.what() << endl;
		output = mxCreateDoubleScalar(0);
	}
	return output;
}

RCP<Hierarchy_double> EpetraSystem::getHierarchy()
{
	RCP<MueLu::EpetraOperator> meo = rcp_static_cast<MueLu::EpetraOperator, Epetra_Operator>(prec);
	return meo->GetHierarchy();
}

//tpetra_double_data_pack implementation

template<> TpetraSystem<double>::TpetraSystem() : MuemexSystem(TPETRA) {}
template<> TpetraSystem<double>::~TpetraSystem() {}

template<typename Scalar>
int TpetraSystem<Scalar>::setup(const mxArray* mxa)
{
	bool success = false;
	try
	{
		A = tpetraLoadMatrix<Scalar>(mxa);
		RCP<MueLu::TpetraOperator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mop = MueLu::CreateTpetraPreconditioner<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(A, *List);
		prec = rcp_implicit_cast<Tpetra::Operator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(mop);
		operatorComplexity = mop->GetHierarchy()->GetOperatorComplexity();
		success = true;
	}
	catch(exception& e)
	{
		mexPrintf("An error occurred while setting up tpetra problem:\n");
		cout << e.what() << endl;
	}
	return success ? IS_TRUE : IS_FALSE;
}

template<>
int TpetraSystem<double>::status()
{
	mexPrintf("**** Problem ID %d [MueLu_Tpetra] ****\n", id);
	if(!A.is_null())
		mexPrintf("Matrix: %dx%d w/ %d nnz\n", A->getGlobalNumRows(), A->getGlobalNumCols(), A->getGlobalNumEntries());
	mexPrintf("Operator Complexity: %f\n", operatorComplexity);
	if(!List.is_null())
	{
		mexPrintf("Parameter List:\n");
		List->print();
	}
	mexPrintf("\n");
	return IS_TRUE;
}

template<typename Scalar>
RCP<MueLu::Hierarchy<Scalar, mm_GlobalOrd, mm_LocalOrd, mm_node_t>> TpetraSystem<Scalar>::getHierarchy()
{
	RCP<MueLu::TpetraOperator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mueluOp = rcp_static_cast<MueLu::TpetraOperator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>, Tpetra::Operator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(prec);
	return mueluOp->GetHierarchy();
}

//tpetra_complex_data_pack implementation

#ifdef HAVE_COMPLEX_SCALARS
template<> TpetraSystem<complex_t>::TpetraSystem() : MuemexSystem(TPETRA_COMPLEX) {}
template<> TpetraSystem<complex_t>::~TpetraSystem() {}

template<>
int TpetraSystem<complex_t>::setup(const mxArray* mxa)
{
	bool success = false;
	try
	{
		A = tpetraLoadMatrix<complex_t>(mxa);
		RCP<MueLu::TpetraOperator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mop = MueLu::CreateTpetraPreconditioner<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(A, *List);
		prec = rcp_implicit_cast<Tpetra::Operator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(mop);
		operatorComplexity = mop->GetHierarchy()->GetOperatorComplexity();
		success = true;
	}
	catch(exception& e)
	{
		mexPrintf("An error occurred while setting up a Tpetra problem:\n");
		cout << e.what() << endl;
	}
	return success ? IS_TRUE : IS_FALSE;
}

template<>
int TpetraSystem<complex_t>::status()
{
	mexPrintf("**** Problem ID %d [MueLu_Tpetra (Complex Scalars)] ****\n", id);
	if(!A.is_null())
		mexPrintf("Matrix: %dx%d w/ %d nnz\n", A->getGlobalNumRows(), A->getGlobalNumCols(), A->getGlobalNumEntries());
	mexPrintf("Operator Complexity: %f\n", operatorComplexity);
	if(!List.is_null())
	{
		mexPrintf("Parameter List:\n");
		List->print();
	}
	mexPrintf("\n");
	return IS_TRUE;
}
#endif	//HAVE_COMPLEX_SCALARS

//MuemexSystemList namespace implementation

void MuemexSystemList::clearAll()
{
	//When items are cleared, RCPs will auto-delete the datapacks
	list.clear();
}

/* add - Adds an MuemexSystem to the list.
   Parameters:
   D	   - The MuemexSystem. [I]
   Returns: problem id number of D
*/
int MuemexSystemList::add(RCP<MuemexSystem> D)
{
	TEUCHOS_ASSERT(!D.is_null());
	D->id = nextID;
	nextID++;
	list.push_back(D);
	return D->id;
}	/*end add*/

/* find - Finds problem by id
   Parameters:
   id	   - ID number [I]
   Returns: pointer to MuemexSystem matching 'id', if found, NULL if not
   found.
*/
RCP<MuemexSystem> MuemexSystemList::find(int id)
{
	if(isInList(id))
	{
		for(auto problem : list)
		{
			if(problem->id == id)
				return problem;
		}
	}
	//auto-inits to NULL
	RCP<MuemexSystem> notFound;
	return notFound;
}/*end find*/

/* remove - Removes problem by id
   Parameters:
   id	   - ID number [I]
   Returns: IS_TRUE if remove was succesful, IS_FALSE otherwise
*/
int MuemexSystemList::remove(int id)
{
	int index = -1;
	for(int i = 0; i < int(list.size()); i++)
	{
		if(list[i]->id == id)
		{
			index = i;
			break;
		}
	}
	if(index == -1)
	{
		mexErrMsgTxt("Error: Tried to clean up a problem that doesn't exist.");
		return IS_FALSE;
	}
	list.erase(list.begin() + index);
	return IS_TRUE;
}/*end remove*/

/* size - Number of stored problems */
int MuemexSystemList::size()
{
	return list.size();
}

/* Returns the status of all members of the list
   Returns IS_TRUE
*/
int MuemexSystemList::status_all()
{
	//This prints all the existing problems in ascending order by ID
	for(int i = 0; i < nextID; i++)
	{
		for(auto problem : list)
		{
			if(problem->id == i)
			{
				problem->status();
				break;
			}
		}
	}
	return IS_TRUE;
}/*end status_all */

bool MuemexSystemList::isInList(int id)
{
	bool rv = false;
	for(auto problem : list)
	{
		if(problem->id == id)
		{
			rv = true;
			break;
		}
	}
	return rv;
}

/* sanity_check - sanity checks the first couple of arguements and returns the
   program mode.
   Parameters:
   nrhs	   - Number of program inputs [I]
   prhs	   - The problem inputs [I]
   Return value: Which mode to run the program in.
*/

MODE_TYPE sanity_check(int nrhs, const mxArray *prhs[])
{
	MODE_TYPE rv = MODE_ERROR;
	double *modes;
	/* Check for mode */
	if(nrhs == 0)
		mexErrMsgTxt("Error: Invalid Inputs\n");
	/* Pull mode data from 1st Input */
	modes = mxGetPr(prhs[0]);
	switch ((MODE_TYPE) modes[0])
	{
		case MODE_SETUP:
			if(nrhs > 1 && mxIsSparse(prhs[1]))
			{
				if(nrhs > 3 && mxIsSparse(prhs[2]) && mxIsSparse(prhs[3]))
					rv = MODE_ERROR;
				else
					rv = MODE_SETUP;
			}
			else
			{
				mexErrMsgTxt("Error: Invalid input for setup\n");
			}
			break;
		case MODE_SOLVE:
			//problem ID and matrix or rhs must be numeric
			if(nrhs >= 2 && mxIsNumeric(prhs[1]) && mxIsNumeric(prhs[2]))
				rv = MODE_SOLVE;
			else
				mexErrMsgTxt("Error: Invalid input for solve\n");
			break;
		case MODE_CLEANUP:
			if(nrhs == 1 || nrhs == 2)
				rv = MODE_CLEANUP;
			else
				mexErrMsgTxt("Error: Extraneous args for cleanup\n");
			break;
		case MODE_STATUS:
			if(nrhs == 1 || nrhs == 2)
				rv = MODE_STATUS;
			else
				mexErrMsgTxt("Error: Extraneous args for status\n");
			break;
		case MODE_AGGREGATE:
			if(nrhs > 1 && mxIsSparse(prhs[1]))
				//Uncomment the next line and remove one after when implementing aggregate mode
				//rv = MODE_AGGREGATE;
				rv = MODE_ERROR;
			else
				mexErrMsgTxt("Error: Invalid input for aggregate\n");
			break;
		case MODE_GET:
			if(nrhs < 4 || nrhs > 5)
				mexErrMsgTxt("Error: Wrong number of args for get\n");
			else
				rv = MODE_GET;
			break;
		default:
			  printf("Mode number = %d\n", (int) modes[0]);
			  mexErrMsgTxt("Error: Invalid input mode\n");
	};
	return rv;
}	
/*end sanity_check*/

void csc_print(int n, int* rowind, int* colptr, double* vals)
{
	int i, j;
	for(i = 0; i < n; i++)
	{
		for(j = colptr[i]; j < colptr[i + 1]; j++)
		{
			mexPrintf("%d %d %20.16e\n", rowind[j], i, vals[j]);
		}
	}
}

void parse_list_item(RCP<ParameterList> List, char *option_name, const mxArray *prhs)
{
	//List shouldn't be NULL but if it is, initialize here
	if(List.is_null())
	{
		List = rcp(new ParameterList);
	}
	mxClassID cid;
	int i, M, N, *opt_int;
	char *opt_char;
	double *opt_float;
	string opt_str;
	RCP<ParameterList> sublist = rcp(new ParameterList);
	mxArray *cell1, *cell2;
	/* Pull relevant info the the option value */
	cid = mxGetClassID(prhs);
	M = mxGetM(prhs);
	N = mxGetN(prhs);
	/* Add to the Teuchos list */
	switch(cid)
	{
		case mxCHAR_CLASS:
			// String
			opt_char = mxArrayToString(prhs);
			opt_str = opt_char;
			List->set(option_name, opt_str);
			if(strcmp(option_name, MUEMEX_INTERFACE) == 0)
			{
				if(strcmp(opt_str.c_str(), "epetra") == 0)
					useEpetra = true;
				else if(strcmp(opt_str.c_str(), "tpetra") == 0)
					useEpetra = false;
			}
			mxFree(opt_char);
			break;
		case mxDOUBLE_CLASS:
		case mxSINGLE_CLASS:
			// Single or double, real or complex
			if(mxIsComplex(prhs))
			{
				opt_float = mxGetPr(prhs);
				double* opt_float_imag = mxGetPi(prhs);
				//assuming user wants std::complex<double> here...
				if(M == 1 && N == 1)
				{
					List->set(option_name, complex_t(*opt_float, *opt_float_imag));
				}
				else if(M == 0 || N == 0)
				{
					List->set(option_name, (complex_t*) NULL);
				}
				else
				{
					List->set(option_name, xpetraLoadMatrix<complex_t>(prhs));
				}
			}
			else
			{
				opt_float = mxGetPr(prhs);
				if(M == 1 && N == 1 && MMISINT(opt_float[0]))
				{
					List->set(option_name, (int) opt_float[0]);
				}
				else if(M == 1 && N == 1)
				{
					List->set(option_name, opt_float[0]);
				}
				else if(M == 0 || N == 0)
				{
					List->set(option_name, (double*) NULL);
				}
				else
				{
					List->set(option_name, xpetraLoadMatrix<double>(prhs));
				}
			}
			break;
		case mxLOGICAL_CLASS:
		// Bool
		if(M == 1 && N == 1)
			List->set(option_name, mxIsLogicalScalarTrue(prhs));
		else
			List->set(option_name, mxGetLogicals(prhs));
			//NTS: The else probably doesn't work.
		break;
		case mxINT8_CLASS:
		case mxUINT8_CLASS:
		case mxINT16_CLASS:
		case mxUINT16_CLASS:
		case mxINT32_CLASS:
		case mxUINT32_CLASS:
			// Integer
			opt_int = (int*) mxGetData(prhs);
			if(M == 1 && N == 1)
				List->set(option_name, opt_int[0]);
			else
				List->set(option_name, opt_int);
			break;
			// NTS: 64-bit ints will break on a 32-bit machine.	 We
			// should probably detect machine type, or somthing, but that would
			// involve a non-trivial quantity of autoconf kung fu.
		case mxCELL_CLASS:
			// Interpret a cell list as a nested teuchos list.
			// NTS: Assuming that it's a 1D row ordered array
			for(i = 0; i < N; i += 2)
			{
				cell1 = mxGetCell(prhs, i);
				cell2 = mxGetCell(prhs, i + 1);
				if(!mxIsChar(cell1))
					mexErrMsgTxt("Error: Input options are not in ['parameter',value] format!\n");
				opt_char = mxArrayToString(cell1);
				parse_list_item(sublist, opt_char, cell2);
				List->set(option_name, *sublist);
				mxFree(opt_char);
			}
			break;
		case mxINT64_CLASS:
		case mxUINT64_CLASS:
		case mxFUNCTION_CLASS:
		case mxUNKNOWN_CLASS:
		case mxSTRUCT_CLASS:
		default:
			mexPrintf("Error parsing input option: %s [type=%d]\n", option_name, cid);
			mexErrMsgTxt("Error: An input option is invalid!\n");
	};
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
/* build_teuchos_list - takes the inputs (barring the solver mode and
  matrix/rhs) and turns them into a Teuchos list.
   Parameters:
   nrhs	   - Number of program inputs [I]
   prhs	   - The problem inputs [I]
   Return value: Teuchos list containing all parameters passed in by the user.
*/
RCP<ParameterList> build_teuchos_list(int nrhs, const mxArray *prhs[])
{
	printf("Building teuchos list with %.1f parameters.", 1.0 * nrhs / 2);
	RCP<ParameterList> TPL = rcp(new ParameterList);
	char* option_name;
	for(int i = 0; i < nrhs; i += 2)
	{
		if(i == nrhs - 1 || !mxIsChar(prhs[i]))
			mexErrMsgTxt("Error: Input options are not in ['parameter',value] format!\n");
		/* What option are we setting? */
		option_name = mxArrayToString(prhs[i]);
		cout << "Parsing parameter with name " << string(option_name) << endl;
		/* Parse */
		parse_list_item(TPL, option_name, prhs[i + 1]);
		/* Free memory */
		mxFree(option_name);
	}
	TPL->print(std::cout);
	return TPL;
}
/*end build_teuchos_list*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(sizeof(int) != sizeof(mwIndex))
		rewrap_ints = true;
	double* id;
	int rv;
	//Arrays representing vectors
	string intf;
	RCP<ParameterList> List;
	MODE_TYPE mode;
	RCP<MuemexSystem> D;
	/* Sanity Check Input */
	mode = sanity_check(nrhs, prhs);
	switch(mode)
	{
		case MODE_SETUP:
		{
			try
			{
				double oc = 0;
				if(nrhs > 2)
					List = build_teuchos_list(nrhs - 2, &(prhs[2]));
				else
					List = rcp(new ParameterList);
				if(mxIsComplex(prhs[1]))
				{
					//Abort if input is complex but complex isn't supported
					#ifndef HAVE_COMPLEX_SCALARS
					mexPrintf("Error: Complex scalars unsupported by this build of Trilinos.\n");
					throw runtime_error("Complex scalars not supported.");
					#endif
				}
				intf = List->get(MUEMEX_INTERFACE, "tpetra");
				List->remove(MUEMEX_INTERFACE);		//no longer need this parameter
				if(intf == "epetra")
				{
					if(mxIsComplex(prhs[1]))
					{
						mexPrintf("Error: Attempting to use complex-valued matrix with Epetra, which is unsupported.\n");
						mexPrintf("Use Tpetra with complex matrices instead.\n");
						throw runtime_error("Tried to use complex matrix with Epetra");
					}
					RCP<EpetraSystem> dp = rcp(new EpetraSystem());
					dp->List = List;
					dp->setup(prhs[1]);
					oc = dp->operatorComplexity;
					D = rcp_implicit_cast<MuemexSystem>(dp);
				}
				else if(intf == "tpetra")
				{
					//infer scalar type from prhs (can be double or std::complex<double>)
					if(mxIsComplex(prhs[1]))
					{
						#ifdef HAVE_COMPLEX_SCALARS
						RCP<TpetraSystem<complex_t>> dp = rcp(new TpetraSystem<complex_t>());
						dp->List = List;
						dp->setup(prhs[1]);
						oc = dp->operatorComplexity;
						D = rcp_implicit_cast<MuemexSystem>(dp);
						#else
						throw runtime_error("Complex scalars not supported.");
						#endif
					}
					else
					{
						RCP<TpetraSystem<double>> dp = rcp(new TpetraSystem<double>());
						dp->List = List;
						dp->setup(prhs[1]);
						oc = dp->operatorComplexity;
						D = rcp_implicit_cast<MuemexSystem>(dp);
					}
				}
				rv = MuemexSystemList::add(D);
				mexPrintf("Set up problem #%d\n", rv);
				if(nlhs > 0)
				{
					plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
					*((int*) mxGetData(plhs[0])) = rv;
					//output OC also
					if(nlhs > 1)
					{
						plhs[1] = mxCreateDoubleScalar(oc);
					}
				}
				mexLock();
			}
			catch(exception& e)
			{
				mexPrintf("An error occurred during setup routine:\n");
				cout << e.what() << endl;
				if(nlhs > 0)
				{
					plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
					*((int*) mxGetData(plhs[0])) = -1;	//output something that is obviously an error value
				}
				if(nlhs > 1)
				{
					plhs[1] = mxCreateDoubleScalar(0);
				}
			}
			break;
		}
		case MODE_SOLVE:
		{
			int iters;
			try
			{
				bool reuse;
				//MODE_SOLVE, probID, matrix, bVec, params OR
				//MODE_SOLVE, probID, bVec, params
				//prhs[0] holds the MODE_SOLVE enum value
				//if reusing, either nrhs == 3 or prhs[3] is not a string
				//(because bVec can't be a string)
				if(MuemexSystemList::size() == 0)
					throw runtime_error("No linear systems are set up.");
				if(nrhs == 3 || (nrhs > 3 && mxGetClassID(prhs[3]) == mxCHAR_CLASS))
				{
					//No matrix supplied as argument, use one from setup
					if(nrhs > 3)
						List = build_teuchos_list(nrhs - 3, &prhs[3]);
					else
						List = rcp(new ParameterList);
					reuse = true;
				}
				else
				{
					if(nrhs > 4)
						List = build_teuchos_list(nrhs - 4, &prhs[4]);
					else
						List = rcp(new ParameterList);
					reuse = false;
				}
				if(List->isType<string>("Output Style"))
				{
					int type = strToOutputStyle(List->get("Output Style", "Belos::Brief").c_str());
					List->remove("Output Style");
					//Reset the ParameterList entry to be of type int instead of string
					List->set("Output Style", type);
				}
				//Convert Belos msg type string to int in List
				//Note: if the parameter value is already an int, don't touch it.
				if(List->isType<string>("Verbosity"))
				{
					//Errors + Warnings is already the default Belos verbosity setting
					int verb = getBelosVerbosity(List->get("Verbosity", "Belos::Errors + Belos::Warnings").c_str());
					List->remove("Verbosity");
					List->set("Verbosity", verb);
				}
				int probID = parseInt(prhs[1]);
				RCP<MuemexSystem> dp = MuemexSystemList::find(probID);
				if(dp.is_null())
					throw runtime_error("Problem handle not allocated.");
				//get pointer to MATLAB array that will be "B" or "rhs" multivector
				const mxArray* rhs = reuse ? prhs[2] : prhs[3];
				switch(dp->type)
				{
					case EPETRA:
					{
						RCP<EpetraSystem> esys = rcp_static_cast<EpetraSystem, MuemexSystem>(dp);
						RCP<Epetra_CrsMatrix> matrix;
						if(reuse)
							matrix = esys->GetMatrix();
						else
							matrix = epetraLoadMatrix(prhs[2]);
						plhs[0] = esys->solve(List, matrix, rhs, iters);
						break;
					}
					case TPETRA:
					{
						RCP<TpetraSystem<double>> tsys = rcp_static_cast<TpetraSystem<double>, MuemexSystem>(dp);
						RCP<Tpetra_CrsMatrix_double> matrix;
						if(reuse)
							matrix = tsys->GetMatrix();
						else
							matrix = tpetraLoadMatrix<double>(prhs[2]);
						plhs[0] = tsys->solve(List, matrix, rhs, iters);
						break;
					}
					case TPETRA_COMPLEX:
					{
						RCP<TpetraSystem<complex_t>> tsys = rcp_static_cast<TpetraSystem<complex_t>, MuemexSystem>(dp);
						RCP<Tpetra_CrsMatrix_complex> matrix;
						if(reuse)
							matrix = tsys->GetMatrix();
						else
							matrix = tpetraLoadMatrix<complex_t>(prhs[2]);
						plhs[0] = tsys->solve(List, matrix, rhs, iters);
						break;
					}
				}
				if(nlhs > 1)
				{
					plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
					int* itersOutput = (int*) mxGetData(plhs[1]);
					*itersOutput = iters;
				}
			}
			catch(exception& e)
			{
				mexPrintf("An error occurred during the solve routine:\n");
				cout << e.what() << endl;
			}
			break;
		}
		case MODE_CLEANUP:
		{
			try
			{
				mexPrintf("MueMex in cleanup mode.\n");
				if(MuemexSystemList::size() > 0 && nrhs == 1)
				{
					/* Cleanup all problems */
					for(int i = 0; i < MuemexSystemList::size(); i++)
						mexUnlock();
					MuemexSystemList::clearAll();
					rv = 1;
				}
				else if(MuemexSystemList::size() > 0 && nrhs == 2)
				{
					/* Cleanup one problem */
					int probID = (int) *((double*) mxGetData(prhs[1]));
					mexPrintf("Cleaning up problem #%d\n", probID);
					rv = MuemexSystemList::remove(probID);
					if(rv)
						mexUnlock();
				}	/*end elseif*/
				else
				{
					rv = 0;
				}
				/* Set return value */
				plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
				id = (double*) mxGetData(plhs[0]);
				*id = double(rv);
			}
			catch(exception& e)
			{
				mexPrintf("An error occurred during the cleanup routine:\n");
				cout << e.what() << endl;
				if(nlhs > 0)
				{
					plhs[0] = mxCreateNumericMatrix(1,1, mxINT32_CLASS, mxREAL);
					*(int*) mxGetData(plhs[0]) = 0;
				}
			}
			break;
		}
		case MODE_STATUS:
		{
			try
			{
				//mexPrintf("MueMex in status checking mode.\n");
				if(MuemexSystemList::size() > 0 && nrhs == 1)
				{
				  /* Status check on all problems */
					rv = MuemexSystemList::status_all();
				}/*end if*/
				else if(MuemexSystemList::size() > 0 && nrhs == 2)
				{
					/* Status check one problem */
					int probID = parseInt(prhs[1]);
					D = MuemexSystemList::find(probID);
					if(D.is_null())
						throw runtime_error("Error: Problem handle not allocated.\n");
					rv = D->status();
				}/*end elseif*/
				else
					mexPrintf("No problems set up.\n");
					rv = 0;
				/* Set return value */
				plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
				id = (double*) mxGetData(plhs[0]);
				*id = double(rv);
			}
			catch(exception& e)
			{
				mexPrintf("An error occurred during the status routine:\n");
				cout << e.what() << endl;
				if(nlhs > 0)
					plhs[0] = mxCreateDoubleScalar(0);
			}
			break;
		}
		case MODE_GET:
		{
			try
			{
				int probID = parseInt(prhs[1]);
				int levelID = parseInt(prhs[2]);
				char* dataName = mxArrayToString(prhs[3]);
				HierAttribType outputType;
				RCP<MuemexSystem> dp = MuemexSystemList::find(probID);
				if(dp.is_null())
				{
					throw runtime_error("Problem handle not allocated.");
				}
				//See if typeHint was given
				//Note: nrhs is 4 or 5; already been sanity checked
				if(nrhs == 4)
				{
					outputType = strToHierAttribType(dataName);
					if(outputType == UNKNOWN)
						throw runtime_error("Unknown data type for hierarchy attribute. \
											Try passing type name manually to muemex.");
				}
				else if(nrhs == 5)
				{
					//Case insensitive compare with type names
					char* typeName = mxArrayToString(prhs[4]);
					char* iter = typeName;
					while(*iter != '\0')
					{
						*iter = (char) tolower((int) *iter);
						iter++;
					}
					if(strcmp(typeName, "matrix") == 0)
						outputType = MATRIX;
					else if(strcmp(typeName, "multivector") == 0)
						outputType = MULTIVECTOR;
					else if(strcmp(typeName, "lovector") == 0)
						outputType = LOVECTOR;
					else if(strcmp(typeName, "scalar") == 0)
						outputType = SCALAR;
					else
						throw runtime_error("Unknown data type for hierarchy attribute. \
											Must be one of 'matrix', 'multivector', 'lovector' or 'scalar'.");
				}
				plhs[0] = dp->getHierarchyData(string(dataName), outputType, levelID);
			}
			catch(exception& e)
			{
				mexPrintf("An error occurred during the get routine:\n");
				cout << e.what() << endl;
				plhs[0] = mxCreateDoubleScalar(0);
			}
			break;
		}
		case MODE_ERROR:
			mexPrintf("MueMex error.");
			break;
		//TODO: Will implement these modes later.
		case MODE_AGGREGATE:
		default:
			mexPrintf("Mode not supported yet.");
	}
}
