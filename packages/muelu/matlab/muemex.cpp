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
#define MUEMEX_INTERFACE "Problem Type"

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
vector<RCP<muelu_data_pack>> muelu_data_pack_list::list;
int muelu_data_pack_list::nextID = 0;

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
	int i, *rv = new int[N];
	for(i = 0; i < N; i++)
		rv[i] = (int) mwi_array[i];
	return rv;
}

RCP<Epetra_CrsMatrix> epetra_setup(int Nrows, int Ncols, int* rowind, int* colptr, double* vals)
{
	Epetra_SerialComm Comm;
	Epetra_Map RangeMap(Nrows, 0, Comm);
	Epetra_Map DomainMap(Ncols, 0, Comm);
	RCP<Epetra_CrsMatrix> A = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, RangeMap, DomainMap, 0));
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
	return A;
}

RCP<Epetra_CrsMatrix> epetra_setup_from_prhs(const mxArray* mxa, bool rewrap_ints)
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
	RCP<Epetra_CrsMatrix> A = epetra_setup(nr, nc, rowind, colptr, vals);
	if(rewrap_ints)
	{
		delete [] rowind;
		delete [] colptr;
	}
	return A;
}

RCP<Tpetra_CrsMatrix_double> tpetra_setup_real_prhs(const mxArray* mxa, bool rewrap_ints)
{
	//Create a map in order to create the matrix (taken from muelu basic example - complex)
	RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
	//numGlobalIndices is just the number of rows in the matrix	
	const Tpetra::global_size_t numGlobalIndices = mxGetM(mxa);
	const mm_GlobalOrd indexBase = 0;
	RCP<const muemex_map_type> map = rcp(new muemex_map_type(numGlobalIndices, indexBase, comm));
	RCP<Tpetra_CrsMatrix_double> A = Tpetra::createCrsMatrix<double, mm_GlobalOrd, mm_LocalOrd, mm_node_t>(map);
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
	return A;
}

RCP<Tpetra_CrsMatrix_complex> tpetra_setup_complex_prhs(const mxArray* mxa, bool rewrap_ints)
{
	//Create a map in order to create the matrix (taken from muelu basic example - complex)
	RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
	const Tpetra::global_size_t numGlobalIndices = mxGetM(mxa);
	const mm_GlobalOrd indexBase = 0;
	RCP<const muemex_map_type> map = rcp(new muemex_map_type(numGlobalIndices, indexBase, comm));
	RCP<Tpetra_CrsMatrix_complex> A = rcp(new Tpetra_CrsMatrix_complex(map, 0));
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
	return A;
}

/*******************************/
//Use Belos (unpreconditioned) to solve matrix
int epetra_unprec_solve(RCP<ParameterList> SetupList, RCP<ParameterList> TPL, RCP<Epetra_CrsMatrix> A, double* b, double* x, int &iters)
{
	//Set up X and B
	Epetra_Map map = A->DomainMap();
	RCP<Epetra_Vector> xVec = rcp(new Epetra_Vector(map));
	RCP<Epetra_Vector> bVec = rcp(new Epetra_Vector(Epetra_DataAccess::Copy, map, b));
	RCP<Epetra_MultiVector> lhs = rcp_implicit_cast<Epetra_MultiVector>(xVec);
	RCP<Epetra_MultiVector> rhs = rcp_implicit_cast<Epetra_MultiVector>(bVec);
	#ifdef VERBOSE_OUTPUT
	int matSize = A->NumGlobalRows();
	mexPrintf("lhs vec:\n");
	for(int i = 0; i < matSize; i++)
	{
		if(i % 10 == 0)
			mexPrintf("\n");
		mexPrintf("%f ", (*lhs)[0][i]);
	}
	mexPrintf("\n\nrhs vec:\n");
	for(int i = 0; i < matSize; i++)
	{
		if(i % 10 == 0)
			mexPrintf("\n");
		mexPrintf("%f ", (*rhs)[0][i]);
	}
	mexPrintf("\n\n");
	#endif
	RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>> problem = rcp(new Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>(A, lhs, rhs));
	bool set = problem->setProblem();
	TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error, "Linear Problem failed to set up correctly!");
	#ifdef VERBOSE_OUTPUT
	TPL->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::Debug + Belos::FinalSummary + Belos::IterationDetails + Belos::OrthoDetails + Belos::TimingDetails + Belos::StatusTestDetails);
	TPL->set("Output Frequency", 1);
	TPL->set("Output Style", Belos::Brief);
	#else
	TPL->set("Verbosity", Belos::Errors + Belos::Warnings);
	#endif
	string solverName = TPL->get("solver", "GMRES");
	Belos::SolverFactory<double, Epetra_MultiVector, Epetra_Operator> factory;
	RCP<Belos::SolverManager<double, Epetra_MultiVector, Epetra_Operator>> solver = factory.create(solverName, TPL);
	solver->setProblem(problem);
	Belos::ReturnType ret = solver->solve();
	int rv;
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
		rv = IS_FALSE;	  
	}
	xVec->ExtractCopy(x);
	return rv;
}	/*end solve*/

//same as above, but with MueLu-generated preconditioner used on right
int epetra_solve(RCP<ParameterList> SetupList, RCP<ParameterList> TPL, RCP<Epetra_CrsMatrix> A, RCP<Epetra_Operator> prec, double* b, double* x, int &iters)
{
	//Set up X and B
	Epetra_Map map = A->DomainMap();
	RCP<Epetra_Vector> xVec = rcp(new Epetra_Vector(map));
	RCP<Epetra_Vector> bVec = rcp(new Epetra_Vector(Epetra_DataAccess::Copy, map, b));
	RCP<Epetra_MultiVector> lhs = rcp_implicit_cast<Epetra_MultiVector>(xVec);
	RCP<Epetra_MultiVector> rhs = rcp_implicit_cast<Epetra_MultiVector>(bVec);
	#ifdef VERBOSE_OUTPUT
	TPL->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::Debug + Belos::FinalSummary + Belos::IterationDetails + Belos::OrthoDetails + Belos::TimingDetails + Belos::StatusTestDetails);
	TPL->set("Output Frequency", 1);
	TPL->set("Output Style", Belos::Brief);
	#else
	TPL->set("Verbosity", Belos::Errors + Belos::Warnings);
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
	int rv;
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
		rv = IS_FALSE;
	}
	xVec->ExtractCopy(x);
	return rv;
}

int tpetra_double_solve(RCP<ParameterList> SetupList, RCP<ParameterList> TPL, RCP<Tpetra_CrsMatrix_double> A, RCP<Tpetra::Operator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> prec, double* b, double* x, int& iters)
{
	int matSize = A->getGlobalNumRows();
	//Define Tpetra vector/multivector types for convenience
	typedef Tpetra::Vector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Vector;
	typedef Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_MultiVector;
	typedef Tpetra::Operator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Operator;
	RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
	//numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?	
	RCP<const muemex_map_type> map = rcp(new muemex_map_type(matSize, (mm_GlobalOrd) 0, comm));
	//Populate x with all 0z initially
	for(int i = 0; i < matSize; i++)
	{
		x[i] = 0;
	}
	ArrayView<double> xArrView(x, matSize);
	RCP<Tpetra_Vector> xVec = rcp(new Tpetra_Vector(map));
	xVec->putScalar(0);
	ArrayView<double> bArrView(b, matSize);
	RCP<Tpetra_Vector> bVec = rcp(new Tpetra_Vector(map, bArrView));
	//cast up to MV for use with Belos
	RCP<Tpetra_MultiVector> lhs = rcp_implicit_cast<Tpetra_MultiVector>(xVec);
	RCP<Tpetra_MultiVector> rhs = rcp_implicit_cast<Tpetra_MultiVector>(bVec);
	//Set iters to 0 in case an error prevents it from being set later
	iters = 0;
	#ifdef VERBOSE_OUTPUT
	TPL->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::Debug + Belos::FinalSummary + Belos::IterationDetails + Belos::OrthoDetails + Belos::TimingDetails + Belos::StatusTestDetails);
	TPL->set("Output Frequency", 1);
	TPL->set("Output Style", Belos::Brief);
	#else
	TPL->set("Verbosity", Belos::Errors + Belos::Warnings);
	#endif
	RCP<Belos::LinearProblem<double, Tpetra_MultiVector, Tpetra_Operator>> problem = rcp(new Belos::LinearProblem<double, Tpetra_MultiVector, Tpetra_Operator>(A, lhs, rhs));
//	RCP<MueLu::TpetraOperator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> prec = MueLu::CreateTpetraPreconditioner<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(A, *SetupList);
	problem->setRightPrec(prec);
	bool set = problem->setProblem();
	TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error, "Linear Problem failed to set up correctly!");
	Belos::SolverFactory<double, Tpetra_MultiVector, Tpetra_Operator> factory;
	string solverName = TPL->get("solver", "GMRES");
	RCP<Belos::SolverManager<double, Tpetra_MultiVector, Tpetra_Operator>> solver = factory.create(solverName, TPL);
	solver->setProblem(problem);
	Belos::ReturnType ret = solver->solve();
	int rv;
	if(ret == Belos::Converged)
	{
		mexPrintf("Success, Belos converged!\n");
		iters = solver->getNumIters();		 
		//Access a raw pointer to the array of doubles in the lhs multivector
		ArrayRCP<const double> solnView = lhs->getData(0);
		for(int i = 0; i < matSize; i++)
		{
			x[i] = solnView[i];
		}
		rv = IS_TRUE;
	}
	else
	{
		mexPrintf("Belos failed to converge.\n");
		iters = 0;
		//x array still has all 0s
		rv = IS_FALSE;
	}
	return rv;
}

//Note: b and x are contiguous arrays of std::complex<double>,
//so they must be allocated & populated when using data from Matlab
int tpetra_complex_solve(RCP<ParameterList> SetupList, RCP<ParameterList> TPL,
RCP<Tpetra_CrsMatrix_complex> A, RCP<Tpetra::Operator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> prec, complex_t* b, complex_t* x, int& iters)
{
	int matSize = A->getGlobalNumRows();
	//Define Tpetra vector/multivector types for convenience
	typedef Tpetra::Vector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Vector;
	typedef Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_MultiVector;
	typedef Tpetra::Operator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Operator;
	RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
	//numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
	RCP<const muemex_map_type> map = rcp(new muemex_map_type(matSize, (mm_GlobalOrd) 0, comm));
	//Fill x with (0 + 0i) initially
	for(int i = 0; i < matSize; i++)
	{
		x[i] = complex_t(0, 0);
	}
	ArrayView<complex_t> xArrView(x, matSize);
	RCP<Tpetra_Vector> xVec = rcp(new Tpetra_Vector(map));
	xVec->putScalar(0);
	ArrayView<complex_t> bArrView(b, matSize);
	RCP<Tpetra_Vector> bVec = rcp(new Tpetra_Vector(map, bArrView));
	//cast up to MV for use with Belos
	RCP<Tpetra_MultiVector> lhs = rcp_implicit_cast<Tpetra_MultiVector>(xVec);
	RCP<Tpetra_MultiVector> rhs = rcp_implicit_cast<Tpetra_MultiVector>(bVec);
	//Set iters to 0 in case an error prevents it from being set later
	iters = 0;
	#ifdef VERBOSE_OUTPUT
	TPL->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::Debug + Belos::FinalSummary + Belos::IterationDetails + Belos::OrthoDetails + Belos::TimingDetails + Belos::StatusTestDetails);
	TPL->set("Output Frequency", 1);
	TPL->set("Output Style", Belos::Brief);
	#else
	TPL->set("Verbosity", Belos::Errors + Belos::Warnings);
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
	int rv;
	if(ret == Belos::Converged)
	{
		mexPrintf("Success, Belos converged!\n");
		iters = solver->getNumIters();		 
		//Access a raw pointer to the array of doubles in the lhs multivector
		ArrayRCP<const complex_t> solnView = lhs->getData(0);
		for(int i = 0; i < matSize; i++)
		{
			x[i] = solnView[i];
		}
		rv = IS_TRUE;
	}
	else
	{
		mexPrintf("Belos failed to converge.\n");
		iters = 0;
		//x array still has all 0s
		rv = IS_FALSE;
	}
	return rv;
}

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

int muelu_epetra_unprec_data_pack::setup(const mxArray* mxa, bool rewrap_ints)
{
	//Just set up matrix, no prec
	A = epetra_setup_from_prhs(mxa, rewrap_ints);
	return IS_TRUE;
}

int muelu_epetra_unprec_data_pack::solve(RCP<ParameterList> TPL, RCP<Epetra_CrsMatrix> Amat, double* b, double* x, int &iters)
{
	return epetra_unprec_solve(List, TPL, Amat, b, x, iters);
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

int muelu_epetra_data_pack::setup(const mxArray* mxa, bool rewrap_ints)
{
	/* Matrix Fill */
	A = epetra_setup_from_prhs(mxa, rewrap_ints);
	prec = MueLu::CreateEpetraPreconditioner(A, *List);
	//underlying the Epetra_Opreator prec is a MueLu::EpetraOperator
	RCP<MueLu::EpetraOperator> meo = rcp_static_cast<MueLu::EpetraOperator>(prec);
	operatorComplexity = meo->GetHierarchy()->GetOperatorComplexity();
	return IS_TRUE;
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
int muelu_epetra_data_pack::solve(RCP<ParameterList> TPL, RCP<Epetra_CrsMatrix> Amat, double* b, double* x, int &iters)
{
	return epetra_solve(List, TPL, Amat, prec, b, x, iters);
}

//tpetra_double_data_pack implementation

muelu_tpetra_double_data_pack::muelu_tpetra_double_data_pack() : muelu_data_pack(TPETRA) {}
muelu_tpetra_double_data_pack::~muelu_tpetra_double_data_pack() {}

int muelu_tpetra_double_data_pack::setup(const mxArray* mxa, bool rewrap_ints)
{
	A = tpetra_setup_real_prhs(mxa, rewrap_ints);
	RCP<MueLu::TpetraOperator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mop = MueLu::CreateTpetraPreconditioner<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(A, *List);
	prec = rcp_implicit_cast<Tpetra::Operator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(mop);
	operatorComplexity = mop->GetHierarchy()->GetOperatorComplexity();
	return IS_TRUE;
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

int muelu_tpetra_double_data_pack::solve(Teuchos::RCP<Teuchos::ParameterList> TPL, Teuchos::RCP<Tpetra_CrsMatrix_double> Amat, double* b, double* x, int &iters)
{
	return tpetra_double_solve(List, TPL, Amat, prec, b, x, iters);
}

//tpetra_complex_data_pack implementation

muelu_tpetra_complex_data_pack::muelu_tpetra_complex_data_pack() : muelu_data_pack(TPETRA_COMPLEX) {}
muelu_tpetra_complex_data_pack::~muelu_tpetra_complex_data_pack() {}

int muelu_tpetra_complex_data_pack::setup(const mxArray* mxa, bool rewrap_ints)
{
	A = tpetra_setup_complex_prhs(mxa, rewrap_ints);
	RCP<MueLu::TpetraOperator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mop = MueLu::CreateTpetraPreconditioner<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(A, *List);
	prec = rcp_implicit_cast<Tpetra::Operator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(mop);
	operatorComplexity = mop->GetHierarchy()->GetOperatorComplexity();
	return IS_TRUE;
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

int muelu_tpetra_complex_data_pack::solve(Teuchos::RCP<Teuchos::ParameterList> TPL, Teuchos::RCP<Tpetra_CrsMatrix_complex> Amat, complex_t* b, complex_t* x, int &iters)
{
	return tpetra_complex_solve(List, TPL, Amat, prec, b, x, iters);
}

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
			// Single or double
			//NTS: Does not deal with complex args
			opt_float = mxGetPr(prhs);
			if(M == 1 && N == 1 && MMISINT(opt_float[0]))
			{
				List->set(option_name, (int) opt_float[0]);
			} /*end if*/
			else if(M == 1 && N == 1)
			{
				List->set(option_name, opt_float[0]);
			}/*end if*/
			else if(M == 0 || N == 0)
			{
				List->set(option_name, (double*) NULL);
			}
			else
			{
				List->set(option_name, opt_float);
			}/*end else*/
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
  matrix/rhs) and turns them into a Teuchos list for use by MLAPI.
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
	return TPL;
}
/*end build_teuchos_list*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double* id;
	int rv;
	//Arrays representing vectors
	int nr;
	int iters;
	string intf;
	RCP<ParameterList> List;
	MODE_TYPE mode;
	bool rewrap_ints = false;
	RCP<muelu_data_pack> D;
	/* Sanity Check Input */
	mode = sanity_check(nrhs, prhs);
	/* Set flag if mwIndex and int are not the same size */
	/* NTS: This can be an issue on 64 bit architectures */
	if(sizeof(int) != sizeof(mwIndex))
		rewrap_ints = true;
	int res;
	switch(mode)
	{
		case MODE_SETUP:
		{
			bool hasOC = false;
			double oc = 0;
			nr = mxGetM(prhs[1]);
			if(nrhs > 2)
				List = build_teuchos_list(nrhs - 2, &(prhs[2]));
			else
				List = rcp(new ParameterList);
			//have epetra (prec) be the default mode for now			
			intf = List->get(MUEMEX_INTERFACE, "epetra");
			//intf holds the problem type, now remove from TPL
			List->remove(MUEMEX_INTERFACE);
			if(intf == "epetra unprec")
			{
				RCP<muelu_epetra_unprec_data_pack> dp = rcp(new muelu_epetra_unprec_data_pack());
				dp->List = List;
				dp->setup(prhs[1], rewrap_ints);
				D = rcp_implicit_cast<muelu_data_pack>(dp);
			}
			else if(intf == "epetra")
			{
				RCP<muelu_epetra_data_pack> dp = rcp(new muelu_epetra_data_pack());
				dp->List = List;
				dp->setup(prhs[1], rewrap_ints);
				hasOC = true;
				oc = dp->operatorComplexity;
				D = rcp_implicit_cast<muelu_data_pack>(dp);
			}
			else if(intf == "tpetra")
			{
				//infer scalar type from prhs (can be double or std::complex<double>)
				if(mxIsComplex(prhs[1]))
				{
					RCP<muelu_tpetra_complex_data_pack> dp = rcp(new muelu_tpetra_complex_data_pack());
					dp->List = List;
					dp->setup(prhs[1], rewrap_ints);
					hasOC = true;
					oc = dp->operatorComplexity;
					D = rcp_implicit_cast<muelu_data_pack>(dp);
				}
				else
				{
					RCP<muelu_tpetra_double_data_pack> dp = rcp(new muelu_tpetra_double_data_pack());
					dp->List = List;
					dp->setup(prhs[1], rewrap_ints);
					hasOC = true;
					oc = dp->operatorComplexity;
					D = rcp_implicit_cast<muelu_data_pack>(dp);
				}
			}
			rv = muelu_data_pack_list::add(D);
			mexPrintf("Set up problem #%d\n", rv);
			plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
			*((int*) mxGetData(plhs[0])) = rv;
			//output OC as well, if applicable
			if(nlhs > 1)
			{
				plhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
				*((double*) mxGetData(plhs[1])) = oc;
				if(!hasOC)
				{
					mexPrintf("Operator complexity not applicable for unpreconditioned problem.\n");
				}
			}
			mexLock();
			break;
		}
		case MODE_SOLVE_NEWMATRIX:
		{
			//Left hand is (x[, iters])
			//Right hand is (id, A, b[, params])
			/* Are there problems set up? */
			if(muelu_data_pack_list::size() == 0)
			{
				mexErrMsgTxt("Error: No problems set up, cannot perform a solve.\n");
				break;
			}
			/* Get the Problem Handle */
			id = (double*) mxGetData(prhs[1]);
			D = muelu_data_pack_list::find(int(*id));
			if(D.is_null())
			{
				mexErrMsgTxt("Error: Problem handle not allocated.\n");
				break;
			}
			/* Pull Problem Size (from b vector) */
			nr = mxGetM(prhs[3]);
			//Sanity check: make sure A is square and b is the right size
			if(nr != int(mxGetM(prhs[2])) || mxGetM(prhs[2]) != mxGetN(prhs[2]))
			{
				mexErrMsgTxt("Error: Size Mismatch in Input\n");
				break;
			}
			/* Teuchos List*/
			if(nrhs > 5)
				List = build_teuchos_list(nrhs - 4, &(prhs[4]));
			else
				List = rcp(new ParameterList);
			switch(D->type)
			{
				//Different matrix/vector types, must static_cast datapack
				//to access type-specific functionality
				case EPETRA_UNPREC:
				{
					RCP<muelu_epetra_unprec_data_pack> dp = rcp_static_cast<muelu_epetra_unprec_data_pack>(D);
					RCP<Epetra_CrsMatrix> A = epetra_setup_from_prhs(prhs[2], rewrap_ints);
					double* b = mxGetPr(prhs[3]);
					//output solution vector (just single, contiguous array for now)
					plhs[0] = mxCreateDoubleMatrix(nr, 1, mxREAL);
					double* x = mxGetPr(plhs[0]);
					res = dp->solve(List, A, b, x, iters);
					break;
				}
				case EPETRA:
				{
					RCP<muelu_epetra_data_pack> dp = rcp_static_cast<muelu_epetra_data_pack>(D);
					RCP<Epetra_CrsMatrix> A = epetra_setup_from_prhs(prhs[2], rewrap_ints);
					double* b = mxGetPr(prhs[3]);
					plhs[0] = mxCreateDoubleMatrix(nr, 1, mxREAL);
					double* x = mxGetPr(plhs[0]);
					res = dp->solve(List, A, b, x, iters);
					break;
				}
				case TPETRA:
				{
					RCP<muelu_tpetra_double_data_pack> dp = rcp_static_cast<muelu_tpetra_double_data_pack>(D);
					RCP<Tpetra_CrsMatrix_double> A = tpetra_setup_real_prhs(prhs[2], rewrap_ints);
					double* b = mxGetPr(prhs[3]);
					plhs[0] = mxCreateDoubleMatrix(nr, 1, mxREAL);
					double* x = mxGetPr(plhs[0]);
					res = dp->solve(List, A, b, x, iters);
					break;
				}
				case TPETRA_COMPLEX:
				{
					RCP<muelu_tpetra_complex_data_pack> dp = rcp_static_cast<muelu_tpetra_complex_data_pack>(D);
					RCP<Tpetra_CrsMatrix_complex> A = tpetra_setup_complex_prhs(prhs[2], rewrap_ints);
					bool complexB = true;
					if(!mxIsComplex(prhs[3]))
					{
						mexPrintf("Warning: Trying to use real vector with complex matrix.\n");
						mexPrintf("Will treat imaginary part of vector as 0.\n");
						complexB = false;
					}
					//Allocate contiguous arrays of complex to pass to tpetra_solve()
					//Will delete within this function
					complex_t* bArr = new complex_t[nr];
					complex_t* xArr = new complex_t[nr];
					//Allocate solution space in Matlab
					plhs[0] = mxCreateDoubleMatrix(nr, 1, mxCOMPLEX);
					double* br = mxGetPr(prhs[3]);
					if(complexB)
					{
						double* bi = mxGetPi(prhs[3]);
						for(int i = 0; i < nr; i++)
						{
							bArr[i] = complex_t(br[i], bi[i]);
						}
					}
					else
					{
						for(int i = 0; i < nr; i++)
						{
							bArr[i] = complex_t(br[i], 0);
						}
					}
					res = dp->solve(List, A, bArr, xArr, iters);
					//Copy the solution into plhs[0]
					double* solR = mxGetPr(plhs[0]);
					double* solI = mxGetPi(plhs[0]);
					for(int i = 0; i < nr; i++)
					{
						//std::complex methods for getting real, imag parts
						solR[i] = real<double>(xArr[i]);
						solI[i] = imag<double>(xArr[i]);
					}
					delete[] bArr;
					delete[] xArr;
					break;
				}
			}
			if(nlhs == 2)
				plhs[1] = mxCreateDoubleScalar(iters);
			mexPrintf("Belos solver returned %d\n", res);
			break;
		}
		case MODE_SOLVE:
		{
			//Left hand is (x[, iters])
			//Right hand is (id, b[, params])
			/* Are there problems set up? */
			if(muelu_data_pack_list::size() == 0)
			{
				mexErrMsgTxt("Error: No problems set up, cannot perform a solve.\n");
				break;
			}
			/* Get the Problem Handle */
			id = (double*) mxGetData(prhs[1]);
			D = muelu_data_pack_list::find(int(*id));
			if(D.is_null())
			{
				mexErrMsgTxt("Error: Problem handle not allocated.\n");
				break;
			}
			/* Pull Problem Size */
			nr = mxGetM(prhs[2]);
			/* Teuchos List*/
			if(nrhs > 4)
				List = build_teuchos_list(nrhs - 3, &(prhs[3]));
			else
				List = rcp(new ParameterList);
			switch(D->type)
			{
				//Different matrix/vector types, must static_cast datapack
				//to access type-specific functionality
				case EPETRA_UNPREC:
				{
					RCP<muelu_epetra_unprec_data_pack> dp = rcp_static_cast<muelu_epetra_unprec_data_pack>(D);
					RCP<Epetra_CrsMatrix> A = dp->GetMatrix();
					//Sanity check: make sure A is square and b is the right size
					if(nr != A->NumMyRows() || A->NumMyRows() != A->NumMyCols())
					{
						mexErrMsgTxt("Error: Size Mismatch in Input\n");
						break;
					}
					double* b = mxGetPr(prhs[2]);
					plhs[0] = mxCreateDoubleMatrix(nr, 1, mxREAL);
					double* x = mxGetPr(plhs[0]);
					res = dp->solve(List, A, b, x, iters);
					break;
				}
				case EPETRA:
				{
					RCP<muelu_epetra_data_pack> dp = rcp_static_cast<muelu_epetra_data_pack>(D);
					RCP<Epetra_CrsMatrix> A = dp->GetMatrix();
					//Sanity check: make sure A is square and b is the right size
					if(nr != A->NumMyRows() || A->NumMyRows() != A->NumMyCols())
					{
						mexErrMsgTxt("Error: Size Mismatch in Input\n");
						break;
					}
					double* b = mxGetPr(prhs[2]);
					plhs[0] = mxCreateDoubleMatrix(nr, 1, mxREAL);
					double* x = mxGetPr(plhs[0]);
					res = dp->solve(List, A, b, x, iters);
					break;
				}
				case TPETRA:
				{
					RCP<muelu_tpetra_double_data_pack> dp = rcp_static_cast<muelu_tpetra_double_data_pack>(D);
					RCP<Tpetra_CrsMatrix_double> A = dp->GetMatrix();
					//Sanity check: make sure A is square and b is the right size
					if(nr != int(A->getGlobalNumRows()) || A->getGlobalNumRows() != A->getGlobalNumCols())
					{
						mexErrMsgTxt("Error: Size Mismatch in Input\n");
						break;
					}
					double* b = mxGetPr(prhs[2]);
					plhs[0] = mxCreateDoubleMatrix(nr, 1, mxREAL);
					double* x = mxGetPr(plhs[0]);
					res = dp->solve(List, A, b, x, iters);
					break;
				}
				case TPETRA_COMPLEX:
				{
					RCP<muelu_tpetra_complex_data_pack> dp = rcp_static_cast<muelu_tpetra_complex_data_pack>(D);
					RCP<Tpetra_CrsMatrix_complex> A = dp->GetMatrix();
					//Sanity check: make sure A is square and b is the right size
					if(nr != int(A->getGlobalNumRows()) || A->getGlobalNumRows() != A->getGlobalNumCols())
					{
						mexErrMsgTxt("Error: Size Mismatch in Input\n");
						break;
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
					complex_t* bArr = new complex_t[nr];
					complex_t* xArr = new complex_t[nr];
					//Allocate solution space in Matlab
					plhs[0] = mxCreateDoubleMatrix(nr, 1, mxCOMPLEX);
					double* br = mxGetPr(prhs[2]);
					if(complexB)
					{
						double* bi = mxGetPi(prhs[2]);
						for(int i = 0; i < nr; i++)
						{
							bArr[i] = complex_t(br[i], bi[i]);
						}
					}
					else
					{
						for(int i = 0; i < nr; i++)
						{
							bArr[i] = complex_t(br[i], 0);
						}
					}
					res = dp->solve(List, A, bArr, xArr, iters);
					//Copy the solution into plhs[0]
					double* solR = mxGetPr(plhs[0]);
					double* solI = mxGetPi(plhs[0]);
					for(int i = 0; i < nr; i++)
					{
						//std::complex methods for getting real, imag parts
						solR[i] = real<double>(xArr[i]);
						solI[i] = imag<double>(xArr[i]);
					}
					delete[] bArr;
					delete[] xArr;
					break;
				}
			}
			if(nlhs == 2)
				plhs[1] = mxCreateDoubleScalar(iters);
			mexPrintf("Belos solver returned %d\n", res);
			break;
		}
		case MODE_CLEANUP:
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
			break;
		}
		case MODE_STATUS:
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
				int probID = *((double*) mxGetData(prhs[1]));
				D = muelu_data_pack_list::find(probID);
				if(D.is_null())
					mexErrMsgTxt("Error: Problem handle not allocated.\n");
				rv = D->status();
			}/*end elseif*/
			else
				mexPrintf("No problems set up.\n");
				rv = 0;
			/* Set return value */
			plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
			id = (double*) mxGetData(plhs[0]);
			*id = double(rv);
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
