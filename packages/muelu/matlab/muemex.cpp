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
#ifdef HAVE_MUELU_MATLAB
#define IS_FALSE 0
#define IS_TRUE 1
#define MUEMEX_ERROR -1

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
#define VERBOSE_OUTPUT

/* Stuff for MATLAB R2006b vs. previous versions */
#if(defined(MX_API_VER) && MX_API_VER >= 0x07030000)
#else
typedef int mwIndex;
#endif

//Declare and call default constructor for data_pack_list vector (starts empty)
vector<RCP<muelu_data_pack>> muelu_data_pack_list::list;
int muelu_data_pack_list::nextID = 0;

//Need a global flag to keep track of Epetra vs. Tpetra for constructing multivectors for param lists
bool rewrap_ints = false;

/**************************************************************/
/**************************************************************/
/**************************************************************/
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

RCP<Epetra_CrsMatrix> epetra_setup(int Nrows, int Ncols, int* rowind, int* colptr, double* vals)
{
	RCP<Epetra_CrsMatrix> A;
	try
	{
		Epetra_SerialComm Comm;
		Epetra_Map RangeMap(Nrows, 0, Comm);
		Epetra_Map DomainMap(Ncols, 0, Comm);
		A = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, RangeMap, DomainMap, 0));
		/* Do the matrix assembly */
		for(int i = 0; i < Ncols; i++)
		{
			for(int j = colptr[i]; j < colptr[i + 1]; j++)
			{
				//		global row, # of entries, value array, column indices array
				A->InsertGlobalValues(rowind[j], 1, &vals[j], &i);
			}
		}
		A->FillComplete(DomainMap, RangeMap);
	}
	catch(exception& e)
	{
		mexPrintf("An error occurred while constructing Epetra matrix:\n");
		cout << e.what() << endl;
	}
	return A;
}

RCP<Epetra_CrsMatrix> epetra_setup_from_prhs(const mxArray* mxa)
{
	RCP<Epetra_CrsMatrix> A;
	try
	{
		int* colptr, *rowind;
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
		A = epetra_setup(nr, nc, rowind, colptr, vals);
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
	return A;
}

//Construct an RCP<Epetra_MultiVector> from a MATLAB array
RCP<Epetra_MultiVector> epetra_setup_multivector(const mxArray* mxa)
{
	RCP<Epetra_MultiVector> mv;
	try
	{
		int nr = mxGetM(mxa);
		int nc = mxGetN(mxa);
		Epetra_SerialComm comm;
		Epetra_Map map(nr, 0, comm);
		mv = rcp(new Epetra_MultiVector(Epetra_DataAccess::Copy, map, mxGetPr(mxa), nr, nc));
	}
	catch(exception& e)
	{
		mexPrintf("Error while constructing Epetra_MultiVector.\n");
		cout << e.what() << endl;
	}
	return mv;
}

RCP<Tpetra_CrsMatrix_double> tpetra_setup_real_prhs(const mxArray* mxa)
{
	bool success = false;
	RCP<Tpetra_CrsMatrix_double> A;
	try
	{
		RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
		//numGlobalIndices is just the number of rows in the matrix	
		const Tpetra::global_size_t numGlobalIndices = mxGetM(mxa);
		const mm_GlobalOrd indexBase = 0;
		RCP<const muemex_map_type> map = rcp(new muemex_map_type(numGlobalIndices, indexBase, comm));
		A = Tpetra::createCrsMatrix<double, mm_GlobalOrd, mm_LocalOrd, mm_node_t>(map);
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
		A->fillComplete();
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

RCP<Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> xpetra_setup_real(const mxArray* mxa)
{
	int nr = mxGetM(mxa);
	int nc = mxGetN(mxa);
	RCP<const Teuchos::Comm<int>> comm;
	typedef Xpetra::Map<mm_LocalOrd, mm_GlobalOrd, mm_node_t> XMap;
	typedef Xpetra::CrsMatrixWrap<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_CrsMatrixWrap;
	RCP<XMap> map = rcp_implicit_cast<XMap>(rcp(new Xpetra::TpetraMap<mm_LocalOrd, mm_GlobalOrd, mm_node_t>(nr, (mm_GlobalOrd) 0, comm)));
	RCP<Xpetra_CrsMatrixWrap> wrapMat = rcp(new Xpetra_CrsMatrixWrap(map,
}

RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> tpetra_setup_real_multivector(const mxArray* mxa)
{
	RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv;
	try
	{
		int nr = mxGetM(mxa);
		int nc = mxGetN(mxa);
		RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
		//numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
		RCP<const muemex_map_type> map = rcp(new muemex_map_type(nr, (mm_GlobalOrd) 0, comm));
		ArrayView<double> arrView(mxGetPr(mxa), nr * nc);
		mv = rcp(new Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(map, arrView, nr, nc));
	}
	catch(exception& e)
	{
		mexPrintf("Error constructing Tpetra MultiVector.\n");
		cout << e.what() << endl;
	}
	return mv;
}

#ifdef HAVE_COMPLEX_SCALARS
RCP<Tpetra_CrsMatrix_complex> tpetra_setup_complex_prhs(const mxArray* mxa)
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

RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> tpetra_setup_complex_multivector(const mxArray* mxa)
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
#endif	//HAVE_COMPLEX_SCALARS

/*******************************/
//Use Belos (unpreconditioned) to solve matrix
int epetra_unprec_solve(RCP<ParameterList> SetupList, RCP<ParameterList> TPL, RCP<Epetra_CrsMatrix> A, double* b, double* x, int numVecs, int &iters)
{
	int rv = IS_FALSE;
	try
	{
		int matSize = A->NumGlobalRows();
		Epetra_Map map = A->DomainMap();
		RCP<Epetra_MultiVector> lhs = rcp(new Epetra_MultiVector(map, numVecs, true));
		//Matlab 2-d array is column-by-column, which is perfect for pulling out multiple vectors into MultiVector
		RCP<Epetra_MultiVector> rhs = rcp(new Epetra_MultiVector(Epetra_DataAccess::Copy, map, b, matSize, numVecs));
		RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>> problem = rcp(new Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>(A, lhs, rhs));
		bool set = problem->setProblem();
		TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error, "Linear Problem failed to set up correctly!");
		#ifdef VERBOSE_OUTPUT
		TPL->get("Verbosity", Belos::Errors | Belos::Warnings | Belos::Debug | Belos::FinalSummary | Belos::IterationDetails | Belos::OrthoDetails | Belos::TimingDetails | Belos::StatusTestDetails);
		TPL->get("Output Frequency", 1);
		TPL->get("Output Style", Belos::Brief);
		#else
		TPL->get("Verbosity", Belos::Errors | Belos::Warnings);
		#endif
		string solverName = TPL->get("solver", "GMRES");
		Belos::SolverFactory<double, Epetra_MultiVector, Epetra_Operator> factory;
		RCP<Belos::SolverManager<double, Epetra_MultiVector, Epetra_Operator>> solver = factory.create(solverName, TPL);
		solver->setProblem(problem);
		Belos::ReturnType ret = solver->solve();
		if(ret == Belos::Converged)
		{
			mexPrintf("Success, Belos converged!\n");
			iters = solver->getNumIters();
			rv = IS_TRUE;
		}
		else
		{
			mexPrintf("Belos failed to converge.\n");
			iters = 0;
		}
		lhs->ExtractCopy(x, matSize);
	}
	catch(exception& e)
	{
		mexPrintf("An error occurred while setting up/running Belos solver:\n");
		cout << e.what() << endl;
	}
	return rv;
}

//same as above, but with MueLu-generated preconditioner used on right
int epetra_solve(RCP<ParameterList> SetupList, RCP<ParameterList> TPL, RCP<Epetra_CrsMatrix> A, RCP<Epetra_Operator> prec, double* b, double* x, int numVecs, int &iters)
{
	int rv = IS_FALSE;
	try
	{
		int matSize = A->NumGlobalRows();
		//Set up X and B
		Epetra_Map map = A->DomainMap();
		RCP<Epetra_MultiVector> lhs = rcp(new Epetra_MultiVector(map, numVecs, true));
		RCP<Epetra_MultiVector> rhs = rcp(new Epetra_MultiVector(Epetra_DataAccess::Copy, map, b, matSize, numVecs));
		#ifdef VERBOSE_OUTPUT
		TPL->get("Verbosity", Belos::Errors | Belos::Warnings | Belos::Debug | Belos::FinalSummary | Belos::IterationDetails | Belos::OrthoDetails | Belos::TimingDetails | Belos::StatusTestDetails);
		TPL->get("Output Frequency", 1);
		TPL->get("Output Style", Belos::Brief);
		#else
		TPL->get("Verbosity", Belos::Errors | Belos::Warnings);
		#endif
		RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>> problem = rcp(new    Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>(A, lhs, rhs));
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
			rv = IS_TRUE;
		}
		else
		{
			mexPrintf("Belos failed to converge.\n");
			iters = 0;
		}
		lhs->ExtractCopy(x, matSize);
	}
	catch(exception& e)
	{
		mexPrintf("Error occurred during Belos solve:\n");
		cout << e.what() << endl;
	}
	return rv;
}

int tpetra_double_solve(RCP<ParameterList> SetupList, RCP<ParameterList> TPL, RCP<Tpetra_CrsMatrix_double> A, RCP<Tpetra::Operator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> prec, double* b, double* x, int numVecs, int& iters)
{
	int rv = IS_FALSE;
	try
	{
		int matSize = A->getGlobalNumRows();
		//Define Tpetra vector/multivector types for convenience
		typedef Tpetra::Vector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Vector;
		typedef Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_MultiVector;
		typedef Tpetra::Operator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Operator;
		RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
		//numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
		RCP<const muemex_map_type> map = rcp(new muemex_map_type(matSize, (mm_GlobalOrd) 0, comm));
		ArrayView<double> xArrView(x, matSize * numVecs);
		ArrayView<double> bArrView(b, matSize * numVecs);
		//initially zero-out lhs
		RCP<Tpetra_MultiVector> lhs = rcp(new Tpetra_MultiVector(map, numVecs, true));
		//...and copy *b 2d array into rhs
		RCP<Tpetra_MultiVector> rhs = rcp(new Tpetra_MultiVector(map, bArrView, matSize, numVecs));
		//Set iters to 0 in case an error prevents it from being set later
		iters = 0;
		#ifdef VERBOSE_OUTPUT
		TPL->get("Verbosity", Belos::Errors | Belos::Warnings | Belos::Debug | Belos::FinalSummary | Belos::IterationDetails | Belos::OrthoDetails | Belos::TimingDetails | Belos::StatusTestDetails);
		TPL->get("Output Frequency", 1);
		TPL->get("Output Style", Belos::Brief);
		#else
		TPL->get("Verbosity", Belos::Errors + Belos::Warnings);
		#endif
		RCP<Belos::LinearProblem<double, Tpetra_MultiVector, Tpetra_Operator>> problem = rcp(new Belos::LinearProblem<double, Tpetra_MultiVector, Tpetra_Operator>(A, lhs, rhs));
		problem->setRightPrec(prec);
		bool set = problem->setProblem();
		TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error, "Linear Problem failed to set up correctly!");
		Belos::SolverFactory<double, Tpetra_MultiVector, Tpetra_Operator> factory;
		string solverName = TPL->get("solver", "GMRES");
		RCP<Belos::SolverManager<double, Tpetra_MultiVector, Tpetra_Operator>> solver = factory.create(solverName, TPL);
		solver->setProblem(problem);
		Belos::ReturnType ret = solver->solve();
		if(ret == Belos::Converged)
		{
			mexPrintf("Success, Belos converged!\n");
			iters = solver->getNumIters();
			rv = IS_TRUE;
		}
		else
		{
			mexPrintf("Belos failed to converge.\n");
			iters = 0;
		}
		//Copy lhs to *x
		lhs->get1dCopy(xArrView, matSize);
	}
	catch(exception& e)
	{
		mexPrintf("Error occurred while running Belos solver:\n");
		cout << e.what() << endl;
	}
	return rv;
}

#ifdef HAVE_COMPLEX_SCALARS
int tpetra_complex_solve(RCP<ParameterList> SetupList, RCP<ParameterList> TPL,
RCP<Tpetra_CrsMatrix_complex> A, RCP<Tpetra::Operator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> prec, complex_t* b, complex_t* x, int numVecs, int& iters)
{
	int rv = IS_FALSE;
	try
	{
		int matSize = A->getGlobalNumRows();
		//Define Tpetra vector/multivector types for convenience
		typedef Tpetra::Vector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Vector;
		typedef Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_MultiVector;
		typedef Tpetra::Operator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Operator;
		RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
		//numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
		RCP<const muemex_map_type> map = rcp(new muemex_map_type(matSize, (mm_GlobalOrd) 0, comm));
		ArrayView<complex_t> xArrView(x, matSize * numVecs);
		ArrayView<complex_t> bArrView(b, matSize * numVecs);
		RCP<Tpetra_MultiVector> lhs = rcp(new Tpetra_MultiVector(map, numVecs, true));
		RCP<Tpetra_MultiVector> rhs = rcp(new Tpetra_MultiVector(map, bArrView, matSize, numVecs));
		//Set iters to 0 in case an error prevents it from being set later
		iters = 0;
		#ifdef VERBOSE_OUTPUT
		TPL->get("Verbosity", Belos::Errors | Belos::Warnings | Belos::Debug | Belos::FinalSummary | Belos::IterationDetails | Belos::OrthoDetails | Belos::TimingDetails | Belos::StatusTestDetails);
		TPL->get("Output Frequency", 1);
		TPL->get("Output Style", Belos::General);
		#else
		TPL->get("Verbosity", Belos::Errors | Belos::Warnings);
		#endif
		RCP<Belos::LinearProblem<complex_t, Tpetra_MultiVector, Tpetra_Operator>> problem = rcp(new Belos::LinearProblem<complex_t, Tpetra_MultiVector, Tpetra_Operator>(A, lhs, rhs));
		problem->setRightPrec(prec);
		bool set = problem->setProblem();
		TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error, "Linear Problem failed to set up correctly!");
		string solverName = TPL->get("solver", "GMRES");
		Belos::SolverFactory<complex_t, Tpetra_MultiVector, Tpetra_Operator> factory;
		RCP<Belos::SolverManager<complex_t, Tpetra_MultiVector, Tpetra_Operator>> solver = factory.create(solverName, TPL);
		solver->setProblem(problem);
		Belos::ReturnType ret = solver->solve();
		if(ret == Belos::Converged)
		{
			mexPrintf("Success, Belos converged!\n");
			iters = solver->getNumIters();
			rv = IS_TRUE;
		}
		else
		{
			mexPrintf("Belos failed to converge.\n");
			iters = 0;
		}
		lhs->get1dCopy(xArrView, matSize);
	}
	catch(exception& e)
	{
		mexPrintf("An error occurred while running Belos solver:\n");
		cout << e.what() << endl;
	}
	return rv;
}
#endif //HAVE_COMPLEX_SCALARS

//data pack base class implementation

muelu_data_pack::muelu_data_pack(DataPackType probType) : id(MUEMEX_ERROR), type(probType) {}
muelu_data_pack::~muelu_data_pack() {}

//muelu_epetra_unprec implementation

muelu_epetra_unprec_data_pack::muelu_epetra_unprec_data_pack() : muelu_data_pack(EPETRA_UNPREC) {}
muelu_epetra_unprec_data_pack::~muelu_epetra_unprec_data_pack() {}

int muelu_epetra_unprec_data_pack::status()
{
	mexPrintf("**** Problem ID %d [Epetra (Unpreconditioned)] ****\n", id);
	if(!A.is_null())
		mexPrintf("Matrix: %dx%d w/ %d nnz\n", A->NumGlobalRows(), A->NumGlobalCols(), A->NumMyNonzeros());
	if(!List.is_null())
	{
		mexPrintf("Parameter List:\n");
		List->print();
	}
	mexPrintf("\n");
	return IS_TRUE;
}

int muelu_epetra_unprec_data_pack::setup(const mxArray* mxa)
{
	//Just set up matrix, no prec
	A = epetra_setup_from_prhs(mxa);
	return IS_TRUE;
}

int muelu_epetra_unprec_data_pack::solve(RCP<ParameterList> TPL, RCP<Epetra_CrsMatrix> Amat, double* b, double* x, int numVecs, int &iters)
{
	return epetra_unprec_solve(List, TPL, Amat, b, x, numVecs, iters);
}

//muelu_epetra_prec implementation

muelu_epetra_data_pack::muelu_epetra_data_pack() : muelu_data_pack(EPETRA) {}
muelu_epetra_data_pack::~muelu_epetra_data_pack() {}

int muelu_epetra_data_pack::status()
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

int muelu_epetra_data_pack::setup(const mxArray* mxa)
{
	bool success = false;
	try
	{
		/* Matrix Fill */
		A = epetra_setup_from_prhs(mxa);
		prec = MueLu::CreateEpetraPreconditioner(A, *List);
		//underlying the Epetra_Opreator prec is a MueLu::EpetraOperator
		RCP<MueLu::EpetraOperator> meo = rcp_static_cast<MueLu::EpetraOperator>(prec);
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

/* muelu_epetra_data_pack::solve - Given two Teuchos lists, one in the muelu_epetra_data_pack, and one of
   solve-time options, this routine calls the relevant solver and returns the solution.
   TPL	   - Teuchos list of solve-time options [I]
   A	   - The matrix to solve with (may not be the one the preconditioned was used for)
   b	   - RHS vector [I]
   x	   - solution vector [O]
   iters   - number of iterations taken [O]
   Returns: IS_TRUE if solve was succesful, IS_FALSE otherwise
*/
int muelu_epetra_data_pack::solve(RCP<ParameterList> TPL, RCP<Epetra_CrsMatrix> Amat, double* b, double* x, int numVecs, int& iters)
{
	return epetra_solve(List, TPL, Amat, prec, b, x, numVecs, iters);
}

//tpetra_double_data_pack implementation

muelu_tpetra_double_data_pack::muelu_tpetra_double_data_pack() : muelu_data_pack(TPETRA) {}
muelu_tpetra_double_data_pack::~muelu_tpetra_double_data_pack() {}

int muelu_tpetra_double_data_pack::setup(const mxArray* mxa)
{
	bool success = false;
	try
	{
		A = tpetra_setup_real_prhs(mxa);
		RCP<MueLu::TpetraOperator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mop = MueLu::CreateTpetraPreconditioner<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(A, *List);
		prec = rcp_implicit_cast<Tpetra::Operator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(mop);
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

int muelu_tpetra_double_data_pack::status()
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

int muelu_tpetra_double_data_pack::solve(Teuchos::RCP<Teuchos::ParameterList> TPL, Teuchos::RCP<Tpetra_CrsMatrix_double> Amat, double* b, double* x, int numVecs, int &iters)
{
	return tpetra_double_solve(List, TPL, Amat, prec, b, x, numVecs, iters);
}

//tpetra_complex_data_pack implementation

#ifdef HAVE_COMPLEX_SCALARS
muelu_tpetra_complex_data_pack::muelu_tpetra_complex_data_pack() : muelu_data_pack(TPETRA_COMPLEX) {}
muelu_tpetra_complex_data_pack::~muelu_tpetra_complex_data_pack() {}

int muelu_tpetra_complex_data_pack::setup(const mxArray* mxa)
{
	bool success = false;
	try
	{
		A = tpetra_setup_complex_prhs(mxa);
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

int muelu_tpetra_complex_data_pack::status()
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

int muelu_tpetra_complex_data_pack::solve(Teuchos::RCP<Teuchos::ParameterList> TPL, Teuchos::RCP<Tpetra_CrsMatrix_complex> Amat, complex_t* b, complex_t* x, int numVecs, int &iters)
{
	return tpetra_complex_solve(List, TPL, Amat, prec, b, x, numVecs, iters);
}
#endif	//HAVE_COMPLEX_SCALARS

//muelu_data_pack_list namespace implementation

void muelu_data_pack_list::clearAll()
{
	//When items are cleared, RCPs will auto-delete the datapacks
	list.clear();
}

/* add - Adds an muelu_data_pack to the list.
   Parameters:
   D	   - The muelu_data_pack. [I]
   Returns: problem id number of D
*/
int muelu_data_pack_list::add(RCP<muelu_data_pack> D)
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
   Returns: pointer to muelu_data_pack matching 'id', if found, NULL if not
   found.
*/
RCP<muelu_data_pack> muelu_data_pack_list::find(int id)
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
	RCP<muelu_data_pack> notFound;
	return notFound;
}/*end find*/

/* remove - Removes problem by id
   Parameters:
   id	   - ID number [I]
   Returns: IS_TRUE if remove was succesful, IS_FALSE otherwise
*/
int muelu_data_pack_list::remove(int id)
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
int muelu_data_pack_list::size()
{
	return list.size();
}

/* Returns the status of all members of the list
   Returns IS_TRUE
*/
int muelu_data_pack_list::status_all()
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

bool muelu_data_pack_list::isInList(int id)
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
			if(nrhs > 2 && mxIsNumeric(prhs[1]) && !mxIsSparse(prhs[2]) && mxIsNumeric(prhs[2]))
				rv = MODE_SOLVE;
			else
				mexErrMsgTxt("Error: Invalid input for solve\n");
			break;
		case MODE_SOLVE_NEWMATRIX:
			if(nrhs > 3 && mxIsNumeric(prhs[1]) && mxIsSparse(prhs[2]) && mxIsNumeric(prhs[3]))
				rv = MODE_SOLVE_NEWMATRIX;
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
		for(j = colptr[i]; j < colptr[i + 1]; j++)
			mexPrintf("%d %d %20.16e\n", rowind[j], i, vals[j]);
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
					RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> tmat = tpetra_setup_complex_prhs(prhs);
					List->set(option_name, tmat);
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
					#ifdef VERBOSE_OUTPUT
					mexPrintf("Creating a Epetra_CrsMatrix in TPL from MATLAB parameters.\n");
					#endif
					RCP<Tpetra_CrsMatrix_double> tcmat = tpetra_setup_real_prhs(prhs);
					typedef Xpetra::TpetraCrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_Tpetra_CrsMatrix;
					typedef Xpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_CrsMatrix;
					typedef Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_Matrix;
					RCP<Xpetra_CrsMatrix> xcmat = rcp_implicit_cast<Xpetra_CrsMatrix>(rcp(new Xpetra_Tpetra_CrsMatrix(tcmat)));
					List->set(option_name, Xpetra::CrsMatrixWrap<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(xcmat));
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
	RCP<ParameterList> TPL = rcp(new ParameterList);
	char* option_name;
	for(int i = 0; i < nrhs; i += 2)
	{
		if(i == nrhs - 1 || !mxIsChar(prhs[i]))
			mexErrMsgTxt("Error: Input options are not in ['parameter',value] format!\n");
		/* What option are we setting? */
		option_name = mxArrayToString(prhs[i]);
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
	int nr;
	int iters;
	string intf;
	RCP<ParameterList> List;
	MODE_TYPE mode;
	RCP<muelu_data_pack> D;
	/* Sanity Check Input */
	mode = sanity_check(nrhs, prhs);
	int res;
	switch(mode)
	{
		case MODE_SETUP:
		{
			try
			{
				bool hasOC = false;
				double oc = 0;
				nr = mxGetM(prhs[1]);
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
					intf = List->get(MUEMEX_INTERFACE, "tpetra");
				}
				else
				{
					intf = List->get(MUEMEX_INTERFACE, "epetra");
				}
				List->remove(MUEMEX_INTERFACE);
				if(intf == "epetra unprec")
				{
					if(mxIsComplex(prhs[1]))
					{
						throw runtime_error("Tried to use complex matrix with Epetra");
					}
					RCP<muelu_epetra_unprec_data_pack> dp = rcp(new muelu_epetra_unprec_data_pack());
					dp->List = List;
					dp->setup(prhs[1]);
					D = rcp_implicit_cast<muelu_data_pack>(dp);
				}
				else if(intf == "epetra")
				{
					if(mxIsComplex(prhs[1]))
					{
						mexPrintf("Error: Attempting to use complex-valued matrix with Epetra, which is unsupported.\n");
						mexPrintf("Use Tpetra with complex matrices instead.\n");
						throw runtime_error("Tried to use complex matrix with Epetra");
					}
					RCP<muelu_epetra_data_pack> dp = rcp(new muelu_epetra_data_pack());
					dp->List = List;
					dp->setup(prhs[1]);
					hasOC = true;
					oc = dp->operatorComplexity;
					D = rcp_implicit_cast<muelu_data_pack>(dp);
				}
				else if(intf == "tpetra")
				{
					//infer scalar type from prhs (can be double or std::complex<double>)
					if(mxIsComplex(prhs[1]))
					{
						#ifdef HAVE_COMPLEX_SCALARS
						RCP<muelu_tpetra_complex_data_pack> dp = rcp(new muelu_tpetra_complex_data_pack());
						dp->List = List;
						dp->setup(prhs[1]);
						hasOC = true;
						oc = dp->operatorComplexity;
						D = rcp_implicit_cast<muelu_data_pack>(dp);
						#else
						throw runtime_error("Complex scalars not supported.");
						#endif
					}
					else
					{
						RCP<muelu_tpetra_double_data_pack> dp = rcp(new muelu_tpetra_double_data_pack());
						dp->List = List;
						dp->setup(prhs[1]);
						hasOC = true;
						oc = dp->operatorComplexity;
						D = rcp_implicit_cast<muelu_data_pack>(dp);
					}
				}
				rv = muelu_data_pack_list::add(D);
				mexPrintf("Set up problem #%d\n", rv);
				if(nlhs > 0)
				{
					plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
					*((int*) mxGetData(plhs[0])) = rv;
					//output OC as well, if applicable
					if(nlhs > 1)
					{
						plhs[1] = mxCreateDoubleScalar(oc);
						if(!hasOC)
						{
							mexPrintf("Warning: Operator complexity not applicable for unpreconditioned problem.\n");
						}
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
					*((int*) mxGetData(plhs[0])) = -1;
				}
				if(nlhs > 1)
				{
					plhs[0] = mxCreateDoubleScalar(0);
				}
			}
			break;
		}
		case MODE_SOLVE_NEWMATRIX:
		{
			try
			{
				//Left hand is (x[, iters])
				//Right hand is (id, A, b[, params])
				/* Are there problems set up? */
				if(muelu_data_pack_list::size() == 0)
				{
					throw runtime_error("No problems set up, cannot perform a solve.");
				}
				/* Get the Problem Handle */
				int probID = parseInt(prhs[1]);
				D = muelu_data_pack_list::find(probID);
				if(D.is_null())
				{
					throw runtime_error("Problem handle not allocated.");
				}
				/* Pull Problem Size (from b vector) */
				nr = mxGetM(prhs[3]);
				//Sanity check: make sure A is square and b is the right size
				if(nr != int(mxGetM(prhs[2])) || mxGetM(prhs[2]) != mxGetN(prhs[2]))
				{
					mexPrintf("Error: Size Mismatch in Input\n");
					plhs[0] = mxCreateDoubleScalar(0);
					if(nlhs == 2)
						plhs[1] = mxCreateDoubleScalar(0);
					throw runtime_error("Size mismatch in input.");
				}
				/* Teuchos List*/
				if(nrhs > 5)
					List = build_teuchos_list(nrhs - 4, &(prhs[4]));
				else
					List = rcp(new ParameterList);
				if(List->isType<string>("Output Style"))
				{
					int type = strToOutputStyle(List->get("Output Style", "Belos::Brief").c_str());
					List->remove("Output Style");
					//Reset the ParameterList entry to be of type int instead of string
					List->set("Output Style", type);
				}
				//Convert Belos msg type string to int in List
				if(List->isType<string>("Verbosity"))
				{
					//Errors + Warnings is already the default Belos verbosity setting
					int verb = getBelosVerbosity(List->get("Verbosity", "Belos::Errors + Belos::Warnings").c_str());
					List->remove("Verbosity");
					List->set("Verbosity", verb);
				}
				switch(D->type)
				{
					//Different matrix/vector types, must static_cast datapack
					//to access type-specific functionality
					case EPETRA_UNPREC:
					{
						RCP<muelu_epetra_unprec_data_pack> dp = rcp_static_cast<muelu_epetra_unprec_data_pack>(D);
						RCP<Epetra_CrsMatrix> A = epetra_setup_from_prhs(prhs[2]);
						int numVecs = mxGetN(prhs[2]);
						double* b = mxGetPr(prhs[3]);
						//output solution vector (just single, contiguous array for now)
						plhs[0] = mxCreateDoubleMatrix(nr, numVecs, mxREAL);
						double* x = mxGetPr(plhs[0]);
						//Should copy sol'n multivec directly to plhs[0]
						res = dp->solve(List, A, b, x, numVecs, iters);
						break;
					}
					case EPETRA:
					{
						RCP<muelu_epetra_data_pack> dp = rcp_static_cast<muelu_epetra_data_pack>(D);
						RCP<Epetra_CrsMatrix> A = epetra_setup_from_prhs(prhs[2]);
						int numVecs = mxGetN(prhs[2]);
						double* b = mxGetPr(prhs[3]);
						plhs[0] = mxCreateDoubleMatrix(nr, numVecs, mxREAL);
						double* x = mxGetPr(plhs[0]);
						res = dp->solve(List, A, b, x, numVecs, iters);
						break;
					}
					case TPETRA:
					{
						RCP<muelu_tpetra_double_data_pack> dp = rcp_static_cast<muelu_tpetra_double_data_pack>(D);
						RCP<Tpetra_CrsMatrix_double> A = tpetra_setup_real_prhs(prhs[2]);
						int numVecs = mxGetN(prhs[2]);					
						double* b = mxGetPr(prhs[3]);
						plhs[0] = mxCreateDoubleMatrix(nr, numVecs, mxREAL);
						double* x = mxGetPr(plhs[0]);
						res = dp->solve(List, A, b, x, numVecs, iters);
						break;
					}
					case TPETRA_COMPLEX:
					{
						#ifdef HAVE_COMPLEX_SCALARS
						RCP<muelu_tpetra_complex_data_pack> dp = rcp_static_cast<muelu_tpetra_complex_data_pack>(D);
						RCP<Tpetra_CrsMatrix_complex> A = tpetra_setup_complex_prhs(prhs[2]);
						int numVecs = mxGetN(prhs[3]);
						bool complexB = true;
						if(!mxIsComplex(prhs[3]))
						{
							mexPrintf("Warning: Trying to use real vector with complex matrix.\n");
							mexPrintf("Will treat imaginary part of vector as 0.\n");
							complexB = false;
						}
						//Allocate contiguous arrays of complex to pass to tpetra_solve()
						//Will delete within this function
						complex_t* bArr = new complex_t[nr * numVecs];
						complex_t* xArr = new complex_t[nr * numVecs];
						//Allocate solution space in Matlab
						plhs[0] = mxCreateDoubleMatrix(nr, numVecs, mxCOMPLEX);
						double* br = mxGetPr(prhs[3]);
						if(complexB)
						{
							double* bi = mxGetPi(prhs[3]);
							for(int n = 0; n < numVecs; n++)
							{
								for(int m = 0; m < nr; m++)
								{
									bArr[n * nr + m] = complex_t(br[n * nr + m], bi[n * nr + m]);
								}
							}
						}
						else
						{
							for(int n = 0; n < numVecs; n++)
							{
								for(int m = 0; m < nr; m++)
								{
									bArr[n * nr + m] = complex_t(br[n * nr + m], 0);
								}
							}
						}
						res = dp->solve(List, A, bArr, xArr, numVecs, iters);
						//Copy the solution into plhs[0]
						double* solR = mxGetPr(plhs[0]);
						double* solI = mxGetPi(plhs[0]);
						for(int n = 0; n < numVecs; n++)
						{
							for(int m = 0; m < nr; m++)
							{
								solR[n * nr + m] = real<double>(xArr[n * nr + m]);
								solI[n * nr + m] = imag<double>(xArr[n * nr + m]);
							}
						}
						delete[] bArr;
						delete[] xArr;
						break;
						#else
						throw runtime_error("Complex scalars not supported.");
						#endif
					}
				}
				if(nlhs == 2)
					plhs[1] = mxCreateDoubleScalar(iters);
				mexPrintf("Belos solver returned %d\n", res);
			}
			catch(exception& e)
			{
				mexPrintf("An error occurred during the solve routine (newmatrix):\n");
				cout << e.what() << endl;
				if(nlhs > 0)
					plhs[0] = mxCreateDoubleScalar(0);
				if(nlhs == 2)
					plhs[1] = mxCreateDoubleScalar(0);
			}
			break;
		}
		case MODE_SOLVE:
		{
			try
			{
				//Left hand is (x[, iters])
				//Right hand is (id, b[, params])
				/* Are there problems set up? */
				if(muelu_data_pack_list::size() == 0)
				{
					throw runtime_error("No problems set up, can't solve.");
				}
				/* Get the Problem Handle */
				int probID = parseInt(prhs[1]);
				mexPrintf("Solving problem #%d", probID);
				D = muelu_data_pack_list::find(probID);
				if(D.is_null())
				{
					throw runtime_error("Problem handle not allocated.");
				}
				/* Pull Problem Size (number of rows in rhs multivec) */
				nr = mxGetM(prhs[2]);
				/* Teuchos List*/
				if(nrhs > 4)
					List = build_teuchos_list(nrhs - 3, &(prhs[3]));
				else
					List = rcp(new ParameterList);
				if(List->isType<string>("Output Style"))
				{
					int type = strToOutputStyle(List->get("Output Style", "Belos::Brief").c_str());
					List->remove("Output Style");
					//Reset the ParameterList entry to be of type int instead of string
					List->set("Output Style", type);
				}
				//Convert Belos msg type string to int in List
				if(List->isType<string>("Verbosity"))
				{
					int verb = getBelosVerbosity(List->get("Verbosity", "Belos::Errors + Belos::Warnings").c_str());
					List->remove("Verbosity");
					List->set("Verbosity", verb);
				}
				switch(D->type)
				{
					//Different matrix/vector types, must static_cast datapack
					//to access type-specific functionality
					case EPETRA_UNPREC:
					{
						RCP<muelu_epetra_unprec_data_pack> dp = rcp_static_cast<muelu_epetra_unprec_data_pack>(D);
						RCP<Epetra_CrsMatrix> A = dp->GetMatrix();
						//grab # of rhs (b) vectors to construct multivector
						int numVecs = mxGetN(prhs[2]);					
						//Sanity check: make sure A is square and b is the right size
						if(nr != A->NumMyRows() || A->NumMyRows() != A->NumMyCols())
						{
							plhs[0] = mxCreateDoubleScalar(0);
							if(nlhs == 2)
								plhs[1] = mxCreateDoubleScalar(0);
							throw runtime_error("Size mismatch in input to solve routine.");
						}
						double* b = mxGetPr(prhs[2]);
						plhs[0] = mxCreateDoubleMatrix(nr, numVecs, mxREAL);
						double* x = mxGetPr(plhs[0]);
						res = dp->solve(List, A, b, x, numVecs, iters);
						break;
					}
					case EPETRA:
					{
						RCP<muelu_epetra_data_pack> dp = rcp_static_cast<muelu_epetra_data_pack>(D);
						RCP<Epetra_CrsMatrix> A = dp->GetMatrix();
						int numVecs = mxGetN(prhs[2]);
						//Sanity check: make sure A is square and b is the right size
						if(nr != A->NumMyRows() || A->NumMyRows() != A->NumMyCols())
						{
							plhs[0] = mxCreateDoubleScalar(0);
							if(nlhs == 2)
								plhs[1] = mxCreateDoubleScalar(0);
							throw runtime_error("Size mismatch in input to solve routine.");
						}
						double* b = mxGetPr(prhs[2]);
						plhs[0] = mxCreateDoubleMatrix(nr, numVecs, mxREAL);
						double* x = mxGetPr(plhs[0]);
						res = dp->solve(List, A, b, x, numVecs, iters);
						break;
					}
					case TPETRA:
					{
						RCP<muelu_tpetra_double_data_pack> dp = rcp_static_cast<muelu_tpetra_double_data_pack>(D);
						RCP<Tpetra_CrsMatrix_double> A = dp->GetMatrix();
						int numVecs = mxGetN(prhs[2]);
						//Sanity check: make sure A is square and b i T & 	get (const std::string &name, T def_value)s the right size
						if(nr != int(A->getGlobalNumRows()) || A->getGlobalNumRows() != A->getGlobalNumCols())
						{
							mexPrintf("Error: Size Mismatch in Input\n");
							plhs[0] = mxCreateDoubleScalar(0);
							if(nlhs == 2)
								plhs[1] = mxCreateDoubleScalar(0);
							throw runtime_error("Size mismatch in input to solve routine.");
						}
						double* b = mxGetPr(prhs[2]);
						plhs[0] = mxCreateDoubleMatrix(nr, numVecs, mxREAL);
						double* x = mxGetPr(plhs[0]);
						res = dp->solve(List, A, b, x, numVecs, iters);
						break;
					}
					case TPETRA_COMPLEX:
					{
						#ifdef HAVE_COMPLEX_SCALARS
						RCP<muelu_tpetra_complex_data_pack> dp = rcp_static_cast<muelu_tpetra_complex_data_pack>(D);
						RCP<Tpetra_CrsMatrix_complex> A = dp->GetMatrix();
						int numVecs = mxGetN(prhs[2]);
						if(nr != int(A->getGlobalNumRows()) || A->getGlobalNumRows() != A->getGlobalNumCols())
						{
							mexPrintf("Error: Size Mismatch in Input\n");
							plhs[0] = mxCreateDoubleScalar(0);
							if(nlhs == 2)
								plhs[1] = mxCreateDoubleScalar(0);
							throw runtime_error("Size mismatch in input to solve routine.");
						}
						bool complexB = true;
						if(!mxIsComplex(prhs[2]))
						{
							mexPrintf("Warning: Trying to use real vector with complex matrix.\n");
							mexPrintf("Will treat imaginary part of vector as 0.\n");
							complexB = false;
						}
						//Allocate contiguous arrays of complex to pass to tpetra_solve()
						//Will delete within this function
						complex_t* bArr = new complex_t[nr * numVecs];
						complex_t* xArr = new complex_t[nr * numVecs];
						//Allocate solution space in Matlab
						plhs[0] = mxCreateDoubleMatrix(nr, numVecs, mxCOMPLEX);
						double* br = mxGetPr(prhs[2]);
						if(complexB)
						{
							double* bi = mxGetPi(prhs[2]);
							for(int n = 0; n < numVecs; n++)
							{
								for(int m = 0; m < nr; m++)
								{
									bArr[n * nr + m] = complex_t(br[n * nr + m], bi[n * nr + m]);
								}
							}
						}
						else
						{
							for(int n = 0; n < numVecs; n++)
							{
								for(int m = 0; m < nr; m++)
								{
									bArr[n * nr + m] = complex_t(br[n * nr + m], 0);
								}
							}
						}
						res = dp->solve(List, A, bArr, xArr, numVecs, iters);
						//Copy the solution into plhs[0]
						double* solR = mxGetPr(plhs[0]);
						double* solI = mxGetPi(plhs[0]);
						for(int n = 0; n < numVecs; n++)
						{
							for(int m = 0; m < nr; m++)
							{
								solR[n * nr + m] = real<double>(xArr[n * nr + m]);
								solI[n * nr + m] = imag<double>(xArr[n * nr + m]);
							}
						}
						delete[] bArr;
						delete[] xArr;
						#else
						throw runtime_error("Complex scalars are not supported.");
						#endif
						break;
					}
				}
				if(nlhs == 2)
					plhs[1] = mxCreateDoubleScalar(iters);
				mexPrintf("Belos solver returned %d\n", res);
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
				if(muelu_data_pack_list::size() > 0 && nrhs == 1)
				{
					/* Cleanup all problems */
					for(int i = 0; i < muelu_data_pack_list::size(); i++)
						mexUnlock();
					muelu_data_pack_list::clearAll();
					rv = 1;
				}
				else if(muelu_data_pack_list::size() > 0 && nrhs == 2)
				{
					/* Cleanup one problem */
					int probID = (int) *((double*) mxGetData(prhs[1]));
					mexPrintf("Cleaning up problem #%d\n", probID);
					rv = muelu_data_pack_list::remove(probID);
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
				if(muelu_data_pack_list::size() > 0 && nrhs == 1)
				{
				  /* Status check on all problems */
					rv = muelu_data_pack_list::status_all();
				}/*end if*/
				else if(muelu_data_pack_list::size() > 0 && nrhs == 2)
				{
					/* Status check one problem */
					int probID = parseInt(prhs[1]);
					D = muelu_data_pack_list::find(probID);
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
		case MODE_ERROR:
			mexPrintf("MueMex error.");
			break;
		//TODO: Will implement these modes later.
		case MODE_AGGREGATE:
		case MODE_SETUP_MAXWELL:
		default:
			mexPrintf("Mode not supported yet.");
	}
}

#else
#error "Do not have MATLAB"
#endif
/*

1-D problem (laplacianfun([x])   x = 3000
Matrix P which is 3 1s in the first column, then shifted down 3 othre cols

muelu(A,'level1',{'P',myP})
Make sure setup doesn't build P (we're providing this)

*/
