/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* This file provides a simple MEX interface to use ML in MATLAB.  This code is
   was originally based off of AdaptiveSA.cpp by Marzio Sala (although the two
   are radically different at this point), and is intended to be used with its
   MATLAB front-end ml.m.

   By: Chris Siefert <csiefer@sandia.gov>
   Version History
   10/28/2013 - Bug fixes for Maxwell2; adding some sanity checks.
   09/27/2013 - Bug fixes for Maxwell2.
   05/23/2013 - Adding thresholding control to Maxwell2
   05/01/2013 - Adding reuse support.
   04/26/2013 - Adding Maxwell2 support
   10/04/2012 - Some minor code cleanup and the memory bugfix.
   04/26/2011 - Bug fixed for error tolerance for ISINT checks. Removing gratuitous prints.
   07/31/2010 - Code cleanup, adding ability to get out AztecOO iteration counts.
   07/22/2010 - Adding ability to handle nested parameter lists via cell arrays
                (e.g. for sending in an ifpack list).
   05/22/2007 - Added fix for when mwIndex and int are not the same size data
                structure (this occurs on 64 bit architectures).
   10/12/2006 - Bug fixed for error tolerance for ISINT checks.
   10/09/2006 - Bug fixed where specifying coordinates for aggregation would not
                work in 2D.
   10/04/2006 - Added aggregate mode, patches for Matlab R2006b + mwArray.
   09/08/2006 - Memory leak fixed.
   08/30/2006 - Added ML_epetra interface functionality.
   08/15/2006 - Added operator complexity handling.
   08/08/2006 - Moved a few variable declarations to allow sun's less forgiving
                compiler to compile the code.
   08/07/2006 - Added status checking functionality.
   08/03/2006 - Added functionality to allow ML to have multiple systems ready
                to be solved at once.  This entailed adding the system handle stuff.
   05/22/2006 - Teuchos-friendly version
   05/16/2006 - Initial Version
*/

#include <stdio.h>
#include <string.h>

#include "ml_config.h"
#include "ml_common.h"
#include "mlmex.h"


#ifdef HAVE_ML_MATLAB
/* Needed for MEX */
using namespace Teuchos;
using namespace ML_Epetra;
using namespace MLAPI;
using namespace std;

extern void _main();

/* Macros */
#define ABS(x)   ((x)>0?(x):(-(x)))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define FLOOR(x) ((int)(x))
#define ISINT(x) ((x)==0?(((x-(int)(x))<1e-15)?true:false):(((x-(int)(x))<1e-15*ABS(x))?true:false))
#define IS_FALSE 0
#define IS_TRUE  1
#define MLMEX_ERROR -1

/* Mode Info */
typedef enum {MODE_SETUP=0, MODE_SOLVE, MODE_CLEANUP, MODE_STATUS, MODE_AGGREGATE, MODE_SETUP_MAXWELL, MODE_SOLVE_NEWMATRIX, MODE_ERROR} MODE_TYPE;

/* MLMEX Teuchos Parameters*/
#define MLMEX_INTERFACE "mlmex: interface"

/* Default values */
#define MLMEX_DEFAULT_LEVELS 10
#define MLMEX_DEFAULT_NUMPDES 1
#define MLMEX_DEFAULT_ADAPTIVEVECS 0
#define MLMEX_DEFAULT_USEDEFAULTNS true

/* Debugging */
//#define VERBOSE_OUTPUT

/* Compatibiliy with older version of MATLAB */
#ifndef MWSIZE_MAX
typedef int mwIndex;
#endif


/**************************************************************/
/**************************************************************/
/**************************************************************/

#ifdef TESTING_FUNCTIONS_ONLY
/* Printout function for testing */
void csc_print(int n,int *rowind,int* colptr, double* vals){
  int i,j;
  for(i=0;i<n;i++)
    for(j=colptr[i];j<colptr[i+1];j++)
      mexPrintf("%d %d %20.16e\n",rowind[j],i,vals[j]);
}/*end if*/
#endif


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* Epetra utility functions */

/* mwIndex_to_int - does a data copy and wraps mwIndex's to ints, in the case
   where they're not the same size.  This routine allocates memory
   WARNING: This does not address overflow.
   Parameters:
   N         - Number of unknowns in array [I]
   mwi_array - Array of mwIndex objects [I]
   Return value: mwIndex objects cast down to ints
*/
int* mwIndex_to_int(int N, mwIndex* mwi_array){
  int i,*rv = new int[N];
  for(i=0;i<N;i++)
    rv[i] = (int)mwi_array[i];
  return rv;
}

/*******************************/
Epetra_CrsMatrix* epetra_setup(int Nrows, int Ncols, int* rowind,int* colptr, double* vals) {
  Epetra_SerialComm Comm;
  Epetra_Map RangeMap(Nrows,0,Comm);
  Epetra_Map DomainMap(Ncols,0,Comm);

  Epetra_CrsMatrix *A=new Epetra_CrsMatrix(Copy,RangeMap,DomainMap,0);

  /* Do the matrix assembly */
  for(int i=0;i<Ncols;i++)
    for(int j=colptr[i];j<colptr[i+1];j++){
      A->InsertGlobalValues(rowind[j],1,&vals[j],&i);
    }
  A->FillComplete(DomainMap,RangeMap);

  return A;
}

/*******************************/
Epetra_CrsMatrix * epetra_setup_from_prhs(const mxArray * mxa, bool rewrap_ints){
  int *colptr, *rowind;
  double * vals=mxGetPr(mxa);

  int nr=mxGetM(mxa);
  int nc=mxGetN(mxa);

  if(rewrap_ints){
    colptr=mwIndex_to_int(nc+1,mxGetJc(mxa));
    rowind=mwIndex_to_int(colptr[nc],mxGetIr(mxa));
  }
  else{
    rowind=(int*)mxGetIr(mxa);
    colptr=(int*)mxGetJc(mxa);
  }

  Epetra_CrsMatrix * A=epetra_setup(nr,nc,rowind,colptr,vals);

  if(rewrap_ints){
    delete [] rowind;
    delete [] colptr;
  }
  return A;
}

/*******************************/
int epetra_solve(Teuchos::ParameterList *SetupList, Teuchos::ParameterList *TPL, Epetra_CrsMatrix *A, Epetra_Operator * Prec,
			       double*b, double*x,int &iters){

  int i, N=A->NumMyRows();
  Epetra_Vector LHS(A->RowMap());
  Epetra_Vector RHS(A->RowMap());

  /* Fill LHS/RHS */
  LHS.PutScalar(0);
  for(i=0;i<N;i++) RHS[i]=b[i];

  /* Create the new Teuchos List */
  Teuchos::ParameterList tmp=*SetupList;
  tmp.setParameters(*TPL);

  /* Define Problem / Solver */
  Epetra_LinearProblem Problem(A, &LHS, &RHS);
  AztecOO solver(Problem);
  solver.SetPrecOperator(Prec);

  /* Get solver options from Teuchos list */
  int    NumIters = tmp.get("krylov: max iterations", 1550);
  double Tol      = tmp.get("krylov: tolerance", 1e-9);
  string type     = tmp.get("krylov: type", "gmres");
  int    output   = tmp.get("krylov: output level",10);
  int    kspace   = tmp.get("krylov: kspace",30);
  string conv     = tmp.get("krylov: conv", "r0");

  /* Set solver options - Solver type*/
  if (type == "cg") solver.SetAztecOption(AZ_solver, AZ_cg);
  else if (type == "cg_condnum") solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  else if (type == "gmres") solver.SetAztecOption(AZ_solver, AZ_gmres);
  else if (type == "gmres_condnum") solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
  else if (type == "fixed point") solver.SetAztecOption(AZ_solver, AZ_fixed_pt);
  //  else ML_THROW("krylov: type has incorrect value (" + type + ")", -1);

  /* Set solver options - Convergence Criterion*/
  if(conv == "r0") solver.SetAztecOption(AZ_conv,AZ_r0);
  else if(conv == "rhs") solver.SetAztecOption(AZ_conv,AZ_rhs);
  else if(conv == "Anorm") solver.SetAztecOption(AZ_conv,AZ_Anorm);
  else if(conv == "noscaled") solver.SetAztecOption(AZ_conv,AZ_noscaled);
  else if(conv == "sol") solver.SetAztecOption(AZ_conv,AZ_sol);
  //  else ML_THROW("krylov: conv has incorrect value (" + conv + ")",-1);

  /* Set solver options - other */
  solver.SetAztecOption(AZ_output, output);
  solver.SetAztecOption(AZ_kspace, kspace);

  /* Do the solve */
  solver.Iterate(NumIters, Tol);

  /* Fill Solution */
  for(i=0;i<N;i++) x[i]=LHS[i];

  /* Solver Output */
  iters=solver.NumIters();
  return IS_TRUE;
}/*end solve*/





/**************************************************************/
/**************************************************************/
/**************************************************************/
/* ML data pack list */
class ml_data_pack_list;
class ml_data_pack_list{
public:
  ml_data_pack_list();
  ~ml_data_pack_list();

  /* add - Adds an ML_DATA_PACK to the list.
     Parameters:
     D       - The ML_DATA_PACK. [I]
     Returns: problem id number of D
  */
  int add(ml_data_pack *D);

  /* find - Finds problem by id
     Parameters:
     id      - ID number [I]
     Returns: pointer to ML_DATA_PACK matching 'id', if found, NULL if not
     found.
  */
  ml_data_pack* find(int id);

  /* remove - Removes problem by id
     Parameters:
     id      - ID number [I]
     Returns: IS_TRUE if remove was succesful, IS_FALSE otherwise
  */
  int remove(int id);

  /* size - Number of stored problems
     Returns: num_probs
  */
  int size();

  /* Returns the status of all members of the list
     Returns IS_TRUE
  */
  int status_all();

protected:
  int num_probs;
  /* Note: This list is sorted */
  ml_data_pack *L;
};


/**************************************************************/
/**************************************************************/
/**************** ml_data_pack class functions ****************/
/**************************************************************/
/**************************************************************/
ml_data_pack::ml_data_pack():id(MLMEX_ERROR),List(NULL),operator_complexity(MLMEX_ERROR),next(NULL){}
ml_data_pack::~ml_data_pack(){if(List) delete List;}


/**************************************************************/
/**************************************************************/
/*************** mlapi_data_pack class functions **************/
/**************************************************************/
/**************************************************************/
mlapi_data_pack::mlapi_data_pack():ml_data_pack(),FineSpace(NULL),A(NULL),Prec(NULL){}

mlapi_data_pack::~mlapi_data_pack(){
   if(FineSpace) delete FineSpace;
   if(A)         delete A;
   if(Prec)      delete Prec;
}/*end destructor*/


/* mlapi_data_pack::status - This function does a status query on the
   MLAPI_DATA_PACK passed in.
   Returns: IS_TRUE
*/
int mlapi_data_pack::status(){
  mexPrintf("**** Problem ID %d [MLAPI] ****\n",id);
  if(A) mexPrintf("Matrix: %dx%d w/ %d nnz\n",A->NumGlobalRows(),A->NumGlobalCols(),A->NumGlobalNonzeros());
  mexPrintf(" Operator complexity = %e\n",operator_complexity);
  if(List){mexPrintf("Parameter List:\n");List->print(cout,1);}
  mexPrintf("\n");
  return IS_TRUE;
}/*end status*/

/**************************************************************/
/**************************************************************/
/**************************************************************/
/* mlapi_data_pack::setup - This function does the setup phase for MLAPI, pulling
   key parameters from the Teuchos list, and calling the aggregation routines
   Parameters:
   N       - Number of unknowns [I]
   rowind  - Row indices of matrix (CSC format) [I]
   colptr  - Column indices of matrix (CSC format) [I]
   vals    - Nonzero values of matrix (CSC format) [I]
   Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
*/
int mlapi_data_pack::setup(int N,int* rowind,int* colptr, double* vals){

  int i,j;
  /* Initialize the workspace and set the output level */
  Init();

  /* Define the space for fine level vectors and operators */
  FineSpace=new Space(N);

  /* Do the matrix assembly --- Someone should really come up with a better way
     of doing this.  Argh! */
  A=new DistributedMatrix(*FineSpace,*FineSpace);
  for(i=0;i<N;i++)
    for(j=colptr[i];j<colptr[i+1];j++)
      (*A)(rowind[j],i)=vals[j];
  A->FillComplete();

  /* Pull key options from Teuchos */
  int numpdes=List->get("PDE equations",MLMEX_DEFAULT_NUMPDES);
  int cutoff=List->get("max levels",MLMEX_DEFAULT_LEVELS);
  int adaptivevecs=List->get("additional candidates",MLMEX_DEFAULT_ADAPTIVEVECS);
  bool UseDefaultNullSpace=List->get("use default null space", MLMEX_DEFAULT_USEDEFAULTNS);

  /* Allocate Memory */
  Prec=new MultiLevelAdaptiveSA(*A,*List,numpdes,cutoff);

  /* Build the Heirarchy */
  if(adaptivevecs>0) Prec->AdaptCompute(UseDefaultNullSpace,adaptivevecs);
  else Prec->Compute();

  operator_complexity=Prec->GetComplexity();
  printf("Smoothed Aggregation: operator complexity = %e\n",operator_complexity);

  Finalize();
  return IS_TRUE;
}/*end setup*/


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* mlapi_data_pack::solve - Given two Teuchos lists, one in the MLAPI_DATA_PACK, and one of
   solve-time options, this routine calls the relevant solver and returns the solution.
   Parameters:
   TPL     - Teuchos list of solve-time options [I]
   N       - Number of unknowns [I]
   b       - RHS vector [I]
   x       - solution vector [O]
   iters   - number of iterations taken [O] (NOT IMPLEMENTED)
   Returns: IS_TRUE if solve was succesful, IS_FALSE otherwise
*/
int mlapi_data_pack::solve(Teuchos::ParameterList *TPL, int N, double*b, double*x,int &iters){
  int i;
  Init();
  MultiVector LHS(A->GetDomainSpace());
  MultiVector RHS(A->GetRangeSpace());

  /* Create the new Teuchos List */
  Teuchos::ParameterList tmp=*List;
  tmp.setParameters(*TPL);

  /* Fill LHS/RHS */
  LHS=0;
  for(i=0;i<N;i++) RHS(i)=b[i];

#ifdef VERBOSE_OUTPUT
  tmp.print(cout,0,true);
#endif

  /* Do the solve */
  Krylov(*A,LHS,RHS,*Prec,tmp);

  /* Fill Solution */
  for(i=0;i<N;i++) x[i]=LHS(i);

  Finalize();

  iters=-1;

  return IS_TRUE;
}/*end solve*/


/**************************************************************/
/**************************************************************/
/************* ml_epetra_data_pack class functions ************/
/**************************************************************/
/**************************************************************/
ml_epetra_data_pack::ml_epetra_data_pack():ml_data_pack(),A(NULL),Prec(NULL){}
ml_epetra_data_pack::~ml_epetra_data_pack(){
  if(A)    delete A;
  if(Prec) delete Prec;
}/*end destructor*/

/* ml_epetra_data_pack_status - This function does a status query on the
   ML_EPETRA_DATA_PACK passed in.
   Returns: IS_TRUE
*/
int ml_epetra_data_pack::status(){
  mexPrintf("**** Problem ID %d [ML_Epetra] ****\n",id);
  if(A) mexPrintf("Matrix: %dx%d w/ %d nnz\n",A->NumGlobalRows(),A->NumGlobalCols(),A->NumMyNonzeros());
  mexPrintf(" Operator complexity = %e\n",operator_complexity);
  if(List){mexPrintf("Parameter List:\n");List->print(cout,1);}
  mexPrintf("\n");
  return IS_TRUE;
}/*end status*/


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* ml_epetra_data_pack::setup - This function does the setup phase for ML_Epetra, pulling
   key parameters from the Teuchos list, and calling the aggregation routines
   Parameters:
   N       - Number of unknowns [I]
   rowind  - Row indices of matrix (CSC format) [I]
   colptr  - Column indices of matrix (CSC format) [I]
   vals    - Nonzero values of matrix (CSC format) [I]
   Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
*/
int ml_epetra_data_pack::setup(int N,int* rowind,int* colptr, double* vals){
  /* Matrix Fill */
  A=epetra_setup(N,N,rowind,colptr,vals);

  /* Allocate Memory */
  Prec=new MultiLevelPreconditioner(*A, *List,false);

  /* Build the Heirarchy */
  Prec->ComputePreconditioner();

  /* Pull Operator Complexity */
  operator_complexity = Prec->GetML_Aggregate()->operator_complexity / Prec->GetML_Aggregate()->fine_complexity;

  return IS_TRUE;
}/*end setup*/


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* ml_epetra_data_pack::solve - Given two Teuchos lists, one in the ml_epetra_data_pack, and one of
   solve-time options, this routine calls the relevant solver and returns the solution.
   Parameters:
   TPL     - Teuchos list of solve-time options [I]
   A       - The matrix to solve with (may not be the one the preconditioned was used for)
   b       - RHS vector [I]
   x       - solution vector [O]
   iters   - number of iterations taken [O]
   Returns: IS_TRUE if solve was succesful, IS_FALSE otherwise
*/
int ml_epetra_data_pack::solve(Teuchos::ParameterList *TPL, Epetra_CrsMatrix *Amat, double*b, double*x,int &iters){
  return epetra_solve(List,TPL,Amat,Prec,b,x,iters);
}/*end solve*/


/**************************************************************/
/**************************************************************/
/************* ml_maxwell_data_pack class functions ************/
/**************************************************************/
/**************************************************************/
ml_maxwell_data_pack::ml_maxwell_data_pack():ml_data_pack(),EdgeMatrix(NULL),GradMatrix(NULL),NodeMatrix(NULL),DummyMatrix(NULL),DummyMatrix2(NULL),Prec(NULL){}

ml_maxwell_data_pack::~ml_maxwell_data_pack(){
  delete Prec;
  delete EdgeMatrix;
  delete NodeMatrix;
  delete GradMatrix;
  delete DummyMatrix;
  delete DummyMatrix2;
}/*end destructor*/

/* ml_epetra_data_pack_status - This function does a status query on the
   ML_EPETRA_DATA_PACK passed in.
   Returns: IS_TRUE
*/
int ml_maxwell_data_pack::status(){
  mexPrintf("**** Problem ID %d [RefMaxwell] ****\n",id);
  if(EdgeMatrix) mexPrintf("Matrix: %dx%d w/ %d nnz\n",EdgeMatrix->NumGlobalRows(),EdgeMatrix->NumGlobalCols(),EdgeMatrix->NumMyNonzeros());
  mexPrintf(" Operator complexity = %e\n",operator_complexity);
  if(List){mexPrintf("Parameter List:\n");List->print(cout,1);}
  mexPrintf("\n");
  return IS_TRUE;
}/*end status*/



/**************************************************************/
/**************************************************************/
/**************************************************************/
int ml_maxwell_data_pack::setup_matrix(const char * name, const mxArray * mxa, bool rewrap_ints){
  if(!strcmp("EdgeMatrix",name)) EdgeMatrix=epetra_setup_from_prhs(mxa,rewrap_ints);
  else if(!strcmp("GradMatrix",name)) GradMatrix=epetra_setup_from_prhs(mxa,rewrap_ints);
  else if(!strcmp("NodeMatrix",name)) NodeMatrix=epetra_setup_from_prhs(mxa,rewrap_ints);
  else return IS_FALSE;
  return IS_TRUE;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
int ml_maxwell_data_pack::setup_preconditioner(){
  // Don't overwrite what the user sent in.
  SetDefaultsRefMaxwell(*List,false);

  // Copy over the overall ML output, since cell arrays don't handle mixed lists well.
  Teuchos::ParameterList & List11_=List->sublist("refmaxwell: 11list");
  Teuchos::ParameterList & List22_=List->sublist("refmaxwell: 22list");
  Teuchos::ParameterList & List11c_=List11_.sublist("edge matrix free: coarse");
  List11_.set("ML output",List->get("ML output",0));
  List22_.set("ML output",List->get("ML output",0));
  List11c_.set("ML output",List->get("ML output",0));

  // Ditto for thresholding
  List11_.set("aggregation: threshold",List->get("aggregation: threshold",0.0));
  List22_.set("aggregation: threshold",List->get("aggregation: threshold",0.0));
  List11c_.set("aggregation: threshold",List->get("aggregation: threshold",0.0));

  // To avoid having to grab other inputs...
  List->set("refmaxwell: disable addon",true);

  if(EdgeMatrix->NumGlobalRows() != GradMatrix->NumGlobalRows() ||
     NodeMatrix->NumGlobalRows() != GradMatrix->NumGlobalCols()){
    printf("Global: EdgeMatrix = %dx%d[%d] NodeMatrix=%dx%d[%d] GradMatrix=%dx%d[%d]\n",
	   EdgeMatrix->NumGlobalRows(),EdgeMatrix->NumGlobalCols(),EdgeMatrix->NumGlobalNonzeros(),
	   NodeMatrix->NumGlobalRows(),NodeMatrix->NumGlobalCols(),NodeMatrix->NumGlobalNonzeros(),
	   GradMatrix->NumGlobalRows(),GradMatrix->NumGlobalCols(),GradMatrix->NumGlobalNonzeros());
    mexErrMsgTxt("Error: Matrix mismatch\n");
    operator_complexity=-1;
    return IS_FALSE;
  }
  /* Avoid aliasing */
  DummyMatrix = new Epetra_CrsMatrix(Copy,EdgeMatrix->RowMap(),0);
  for(int i=0; i<DummyMatrix->NumMyRows(); i++) {
    int GID=DummyMatrix->GRID(i);
    double one=1.0;
    DummyMatrix->InsertGlobalValues(GID,1,&one,&GID);
  }
  DummyMatrix->FillComplete();


  /* Build Hierarchy */
#ifdef ENABLE_MS_MATRIX
  DummyMatrix2=new Epetra_CrsMatrix(*EdgeMatrix);
  Prec=new RefMaxwellPreconditioner(*EdgeMatrix, *GradMatrix,*DummyMatrix2,*NodeMatrix,*DummyMatrix,*List);
#else
  Prec=new RefMaxwellPreconditioner(*EdgeMatrix, *GradMatrix,*NodeMatrix,*DummyMatrix,*List);
#endif

  /* Pull Operator Complexity */
  double nnz;
  Prec->Complexities(operator_complexity,nnz);

  return IS_TRUE;
}


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* ml_epetra_data_pack::solve - Given two Teuchos lists, one in the ml_epetra_data_pack, and one of
   solve-time options, this routine calls the relevant solver and returns the solution.
   Parameters:
   TPL     - Teuchos list of solve-time options [I]

   A       - The matrix to solve with (may not be the one the preconditioned was used for)
   b       - RHS vector [I]
   x       - solution vector [O]
   iters   - number of iterations taken [O]
   Returns: IS_TRUE if solve was succesful, IS_FALSE otherwise
*/
int ml_maxwell_data_pack::solve(Teuchos::ParameterList *TPL, Epetra_CrsMatrix * Amat, double*b, double*x,int &iters){
  return epetra_solve(List,TPL,Amat,Prec,b,x,iters);
}/*end solve*/




/**************************************************************/
/**************************************************************/
/***********  ml_data_pack list class functions ***************/
/**************************************************************/
/**************************************************************/

ml_data_pack_list::ml_data_pack_list():num_probs(0),L(NULL){}
ml_data_pack_list::~ml_data_pack_list(){
  ml_data_pack *old;
  while(L!=NULL){
    old=L; L=L->next;
    delete old;
  }/*end while*/
}/*end destructor*/


/* add - Adds an ml_data_pack to the list.
   Parameters:
   D       - The ml_data_pack. [I]
   Returns: problem id number of D
*/
int ml_data_pack_list::add(ml_data_pack *D){
  int idx_prev=0;
  ml_data_pack *iprev, *icurr;

  if(!L){L=D;D->next=NULL;D->id=1;}
  else {
    /* Find the first numbering gap + add in D appropriately */
    for(iprev=NULL,icurr=L; icurr && icurr->id-idx_prev==1; iprev=icurr,idx_prev=iprev->id,icurr=icurr->next);
    if(!iprev) L=D;
    else iprev->next=D;
    D->id=idx_prev+1;
    D->next=icurr;
    num_probs++;
  }/*end if*/
  return D->id;
}/*end add*/


/* find - Finds problem by id
   Parameters:
   id      - ID number [I]
   Returns: pointer to ml_data_pack matching 'id', if found, NULL if not
   found.
*/
ml_data_pack* ml_data_pack_list::find(int id){
  ml_data_pack *rv;
  for(rv=L;rv && rv->id<id;rv=rv->next);
  if(rv && rv->id==id) return rv;
  else return NULL;
}/*end find*/


/* remove - Removes problem by id
   Parameters:
   id      - ID number [I]
   Returns: IS_TRUE if remove was succesful, IS_FALSE otherwise
*/
int ml_data_pack_list::remove(int id){
  ml_data_pack *iprev, *icurr;
  if(!L) return IS_FALSE;
  for(iprev=NULL,icurr=L; icurr && icurr->id<id; iprev=icurr,icurr=icurr->next);

  if(!icurr || icurr->id!=id) return IS_FALSE;
  else{
    if(!iprev) L=icurr->next;
    else iprev->next=icurr->next;
    delete icurr;
    return IS_TRUE;
  }/*end else*/
}/*end remove*/


/* size - Number of stored problems
   Returns: num_probs
*/
int ml_data_pack_list::size(){return num_probs;}


/* Returns the status of all members of the list
   Returns IS_TRUE
*/
int ml_data_pack_list::status_all(){
  ml_data_pack *curr;
  for(curr=L;curr;curr=curr->next)
    curr->status();
  return IS_TRUE;
}/*end status_all */


/**************************************************************/
/**************************************************************/
/******************* Aggregation Functions ********************/
/**************************************************************/
/**************************************************************/
/* mlmex_aggregate -interface to ML's aggregation routines.
   Parameters:
   N       - Number of unknowns [I]
   rowind  - Row indices of matrix (CSC format) [I]
   colptr  - Column indices of matrix (CSC format) [I]
   vals    - Nonzero values of matrix (CSC format) [I]
   List    - Teuchos parameter list [I]
   agg     - allocated vector which will hold aggregates on return [O]
*/
void mlmex_aggregate(int N,int *colptr, int* rowind, double*vals, Teuchos::ParameterList* List, int* agg){
  ml_epetra_data_pack DPK;
  MultiLevelPreconditioner *Prec;

  /* Minimize work*/
  List->set("max levels",2);

  /* Setup a datapack */
  DPK.List=List;
  DPK.setup(N,rowind,colptr,vals);

  /* Pull the aggregates (by reference) */
  Prec=DPK.GetPreconditioner();

  /* Copy aggregate info over */
  memcpy(agg,Prec->GetML_Aggregate()->aggr_info[0],N*sizeof(int));

  /* Prevent the destructor from eating List */
  DPK.List=NULL;
}/*end mlmex_aggregate*/


/**************************************************************/
/**************************************************************/
/********************* Utility Functions **********************/
/**************************************************************/
/**************************************************************/
/* sanity_check - sanity checks the first couple of arguements and returns the
   program mode.
   Parameters:
   nrhs    - Number of program inputs [I]
   prhs    - The problem inputs [I]
   Return value: Which mode to run the program in.
*/

MODE_TYPE sanity_check(int nrhs, const mxArray *prhs[]){
  MODE_TYPE rv=MODE_ERROR;
  double *modes;
  /* Check for mode */
  if(nrhs==0)
    mexErrMsgTxt("Error: Invalid Inputs\n");

  /* Pull mode data from 1st Input */
  modes=mxGetPr(prhs[0]);

  switch ((MODE_TYPE)modes[0]){
  case MODE_SETUP:
    if(nrhs>1&&mxIsSparse(prhs[1])) {
      if(nrhs>3 && mxIsSparse(prhs[2]) && mxIsSparse(prhs[3]))
	rv=MODE_SETUP_MAXWELL;
      else
	rv=MODE_SETUP;
    }
    else mexErrMsgTxt("Error: Invalid input for setup\n");
    break;
  case MODE_SOLVE:
    if(nrhs>2&&mxIsNumeric(prhs[1])&&!mxIsSparse(prhs[2])&&mxIsNumeric(prhs[2])) rv=MODE_SOLVE;
    else mexErrMsgTxt("Error: Invalid input for solve\n");
    break;
  case MODE_SOLVE_NEWMATRIX:
    if(nrhs>3&&mxIsNumeric(prhs[1])&&mxIsSparse(prhs[2])&&mxIsNumeric(prhs[3])) rv=MODE_SOLVE_NEWMATRIX;
    else mexErrMsgTxt("Error: Invalid input for solve\n");
    break;
  case MODE_CLEANUP:
    if(nrhs==1 || nrhs==2) rv=MODE_CLEANUP;
    else mexErrMsgTxt("Error: Extraneous args for cleanup\n");
    break;
  case MODE_STATUS:
    if(nrhs==1 || nrhs==2) rv=MODE_STATUS;
    else mexErrMsgTxt("Error: Extraneous args for status\n");
    break;
  case MODE_AGGREGATE:
    if(nrhs>1&&mxIsSparse(prhs[1])) rv=MODE_AGGREGATE;
    else mexErrMsgTxt("Error: Invalid input for aggregate\n");
    break;
  default:
    printf("Mode number = %d\n",(int)modes[0]);
    mexErrMsgTxt("Error: Invalid input mode\n");
  };
  return rv;
}/*end sanity_check*/


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* build_teuchos_list - takes the inputs (barring the solver mode and
  matrix/rhs) and turns them into a Teuchos list for use by MLAPI.
   Parameters:
   nrhs    - Number of program inputs [I]
   prhs    - The problem inputs [I]
   Return value: Teuchos list containing all parameters passed in by the user.
*/
void parse_list_item(Teuchos::ParameterList & List,char *option_name,const mxArray *prhs){
  mxClassID cid;
  int i,M,N, *opt_int;
  char *opt_char;
  double *opt_float;
  string opt_str;
  Teuchos::ParameterList sublist;
  mxArray *cell1,*cell2;

  /* Pull relevant info the the option value */
  cid=mxGetClassID(prhs);
  M=mxGetM(prhs);
  N=mxGetN(prhs);

  /* Add to the Teuchos list */
  switch(cid){
  case mxCHAR_CLASS:
    // String
    opt_char=mxArrayToString(prhs);
    opt_str=opt_char;
    List.set(option_name, opt_str);
    mxFree(opt_char);
    break;
  case mxDOUBLE_CLASS:
  case mxSINGLE_CLASS:
    // Single or double
    //NTS: Does not deal with complex args
    opt_float=mxGetPr(prhs);
    if(M==1 && N==1 && ISINT(opt_float[0])) {
      List.set(option_name, (int)opt_float[0]);
    }/*end if*/
    else if(M==1 && N==1){
      List.set(option_name, opt_float[0]);
    }/*end if*/
    else if(M==0 || N==0){
      List.set(option_name,(double*)NULL);
    }
    else{
      List.set(option_name, opt_float);
    }/*end else*/
    break;
  case mxLOGICAL_CLASS:
    // Bool
    if(M==1 && N==1) List.set(option_name, mxIsLogicalScalarTrue(prhs));
    else List.set(option_name,mxGetLogicals(prhs));
    //NTS: The else probably doesn't work.
    break;
  case mxINT8_CLASS:
  case mxUINT8_CLASS:
  case mxINT16_CLASS:
  case mxUINT16_CLASS:
  case mxINT32_CLASS:
  case mxUINT32_CLASS:
    // Integer
    opt_int=(int*)mxGetData(prhs);
    if(M==1 && N==1) List.set(option_name, opt_int[0]);
    else List.set(option_name, opt_int);
    break;
    // NTS: 64-bit ints will break on a 32-bit machine.  We
    // should probably detect machine type, or somthing, but that would
    // involve a non-trivial quantity of autoconf kung fu.
  case mxCELL_CLASS:
    // Interpret a cell list as a nested teuchos list.
    // NTS: Assuming that it's a 1D row ordered array
    for(i=0;i<N;i+=2){
      cell1=mxGetCell(prhs,i);
      cell2=mxGetCell(prhs,i+1);
      if(!mxIsChar(cell1))
	mexErrMsgTxt("Error: Input options are not in ['parameter',value] format!\n");
      opt_char=mxArrayToString(cell1);
      parse_list_item(sublist,opt_char,cell2);
      List.set(option_name,sublist);
      mxFree(opt_char);
    }
    break;
  case mxINT64_CLASS:
  case mxUINT64_CLASS:
  case mxFUNCTION_CLASS:
  case mxUNKNOWN_CLASS:
  case mxSTRUCT_CLASS:
  default:
    mexPrintf("Error parsing input option: %s [type=%d]\n",option_name,cid);
    mexErrMsgTxt("Error: An input option is invalid!\n");
  };
}


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* build_teuchos_list - takes the inputs (barring the solver mode and
  matrix/rhs) and turns them into a Teuchos list for use by MLAPI.
   Parameters:
   nrhs    - Number of program inputs [I]
   prhs    - The problem inputs [I]
   Return value: Teuchos list containing all parameters passed in by the user.
*/
Teuchos::ParameterList* build_teuchos_list(int nrhs,const mxArray *prhs[]){
  Teuchos::ParameterList* TPL=new Teuchos::ParameterList;
  char * option_name;

  for(int i=0;i<nrhs;i+=2){
    if(i==nrhs-1 || !mxIsChar(prhs[i]))
      mexErrMsgTxt("Error: Input options are not in ['parameter',value] format!\n");

    /* What option are we setting? */
    option_name=mxArrayToString(prhs[i]);

    /* Parse */
    parse_list_item(*TPL,option_name,prhs[i+1]);

    /* Free memory */
    mxFree(option_name);
  }/*end for*/

  return TPL;
}/*end build_teuchos_list*/




/**************************************************************/
/**************************************************************/
/**************************************************************/
/* Calling syntax is done correctly by ml.m, please use that.
   The first rhs argument is always the program mode.  The other args are the
   ml.m parameter list in order (barring the mode arg).
*/

/* mexFunction is the gateway routine for the MEX-file. */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
  int i,*id,sz,rv,*rowind,*colptr, *agg;
  double *b, *x, *vals;
  int nr, nc,iters;
  string intf;
  Teuchos::ParameterList* List;
  MODE_TYPE mode;
  static ml_data_pack_list* PROBS=0;
  ml_data_pack *D=0;
  bool rewrap_ints=false;
  ml_maxwell_data_pack * Dhat=0;
  Epetra_CrsMatrix * A;

  /* Sanity Check Input */
  mode=sanity_check(nrhs,prhs);

  /* Set flag if mwIndex and int are not the same size */
  /* NTS: This can be an issue on 64 bit architectures */
  if(sizeof(int)!=sizeof(mwIndex)) rewrap_ints=true;


  switch(mode){
  case MODE_SETUP:
    if(!PROBS) PROBS=new ml_data_pack_list;

    /* Pull Problem Size */
    nr=mxGetM(prhs[1]);
    nc=mxGetN(prhs[1]);

    /* Teuchos List*/
    if(nrhs>2) List=build_teuchos_list(nrhs-2,&(prhs[2]));
    else List=new Teuchos::ParameterList;

    /* Pick the Interface Type */
    intf=List->get(MLMEX_INTERFACE,"epetra");
    if(intf=="mlapi") D=new mlapi_data_pack();
    else D=new ml_epetra_data_pack();
    D->List=List;

    /* Pull Matrix - CSC Format */
    vals=mxGetPr(prhs[1]);
    if(rewrap_ints){
      colptr=mwIndex_to_int(nc+1,mxGetJc(prhs[1]));
      rowind=mwIndex_to_int(colptr[nc],mxGetIr(prhs[1]));
    }
    else{
      rowind=(int*)mxGetIr(prhs[1]);
      colptr=(int*)mxGetJc(prhs[1]);
    }

    /* Construct the Heirarchy */
    D->setup(nr,rowind,colptr,vals);

    /* Add this problem to the list */
    rv=PROBS->add(D);

    /* Set return value(s) */
    plhs[0]=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    id=(int*)mxGetData(plhs[0]);id[0]=rv;
    if(nlhs>1) plhs[1]=mxCreateDoubleScalar(D->operator_complexity);

    /* Cleanup */
    if(rewrap_ints){
      delete [] rowind;
      delete [] colptr;
    }

    /* Lock so we can keep the memory for the heirarchy */
    mexLock();
    break;
  case MODE_SETUP_MAXWELL:
    if(!PROBS) PROBS=new ml_data_pack_list;

    /* Teuchos List*/
    if(nrhs>4) List=build_teuchos_list(nrhs-4,&(prhs[4]));
    else List=new Teuchos::ParameterList;

    /* Setup matrices*/
    Dhat=new ml_maxwell_data_pack();
    Dhat->List=List;

    Dhat->setup_matrix("EdgeMatrix",prhs[1],rewrap_ints);
    Dhat->setup_matrix("NodeMatrix",prhs[2],rewrap_ints);
    Dhat->setup_matrix("GradMatrix",prhs[3],rewrap_ints);

    Dhat->setup_preconditioner();

    /* Add this problem to the list */
    rv=PROBS->add(Dhat);

    /* Set return value(s) */
    plhs[0]=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    id=(int*)mxGetData(plhs[0]);id[0]=rv;
    if(nlhs>1) plhs[1]=mxCreateDoubleScalar(Dhat->operator_complexity);


    /* Lock so we can keep the memory for the heirarchy */
    mexLock();
    break;
  case MODE_SOLVE:
    /* Are there problems set up? */
    if(!PROBS) mexErrMsgTxt("Error: No problems set up, cannot solve.\n");

    /* Get the Problem Handle */
    id=(int*)mxGetData(prhs[1]);
    D=PROBS->find(id[0]);
    if(!D) mexErrMsgTxt("Error: Problem handle not allocated.\n");

    /* Pull Problem Size */
    nr=mxGetM(prhs[2]);
    A = D->GetMatrix();

    /* Pull RHS */
    b=mxGetPr(prhs[2]);

    /* Teuchos List*/
    if(nrhs>4) List=build_teuchos_list(nrhs-3,&(prhs[3]));
    else List=new Teuchos::ParameterList;

    /* Allocate Solution Space */
    plhs[0]=mxCreateDoubleMatrix(nr,1,mxREAL);
    x=mxGetPr(plhs[0]);

    /* Sanity Check Matrix / RHS */
    if(nr!=A->NumMyRows() || A->NumMyRows() != A->NumMyCols())
      mexErrMsgTxt("Error: Size Mismatch in Input\n");

    /* Run Solver */
    D->solve(List,A,b,x,iters);

    /* Output Iteration Count */
    if(nlhs>1){
      plhs[1]=mxCreateDoubleScalar((double)iters);
    }

    /* Cleanup */
    A=0; // we don't own this
    delete List;

    break;
  case MODE_SOLVE_NEWMATRIX:
    /* Are there problems set up? */
    if(!PROBS) mexErrMsgTxt("Error: No problems set up, cannot solve.\n");

    /* Get the Problem Handle */
    id=(int*)mxGetData(prhs[1]);
    D=PROBS->find(id[0]);
    if(!D) mexErrMsgTxt("Error: Problem handle not allocated.\n");

    /* Pull Problem Size */
    nr=mxGetM(prhs[2]);
    nc=mxGetN(prhs[2]);
    if(nr!=D->NumMyRows() && nc!=D->NumMyCols()) mexErrMsgTxt("Error: Problem size mismatch.\n");

    /* Pull RHS */
    b=mxGetPr(prhs[3]);

    /* Teuchos List*/
    if(nrhs>4) List=build_teuchos_list(nrhs-4,&(prhs[4]));
    else List=new Teuchos::ParameterList;

    /* Allocate Solution Space */
    plhs[0]=mxCreateDoubleMatrix(nr,1,mxREAL);
    x=mxGetPr(plhs[0]);

    /* Sanity Check Matrix / RHS */
    if(nr != nc || nr != (int)mxGetM(prhs[2]))
      mexErrMsgTxt("Error: Size Mismatch in Input\n");

    // Grab the input matrix
    A=epetra_setup_from_prhs(prhs[2],rewrap_ints);

    /* Run Solver */
    D->solve(List,A,b,x,iters);

    /* Output Iteration Count */
    if(nlhs>1){
      plhs[1]=mxCreateDoubleScalar((double)iters);
    }

    /* Cleanup */
    delete List;
    delete A;
    break;
  case MODE_CLEANUP:
    if(PROBS && nrhs==1){
      /* Cleanup all problems */
      sz=PROBS->size();
      for(i=0;i<sz;i++) mexUnlock();
      delete PROBS;PROBS=NULL;
      rv=1;
    }/*end if*/
    else if(PROBS && nrhs==2){
      /* Cleanup one problem */
      id=(int*)mxGetData(prhs[1]);
      rv=PROBS->remove(id[0]);
      if(rv) mexUnlock();
    }/*end elseif*/
    else rv=0;

    /* Set return value */
    plhs[0]=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    id=(int*)mxGetData(plhs[0]);id[0]=rv;
    break;

  case MODE_STATUS:
    if(PROBS && nrhs==1){
      /* Status check on all problems */
      rv=PROBS->status_all();
    }/*end if*/
    else if(PROBS && nrhs==2){
      /* Status check one problem */
      id=(int*)mxGetData(prhs[1]);
      D=PROBS->find(id[0]);
      if(!D) mexErrMsgTxt("Error: Problem handle not allocated.\n");
      rv=D->status();
    }/*end elseif*/
    else rv=0;

    /* Set return value */
    plhs[0]=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    id=(int*)mxGetData(plhs[0]);id[0]=rv;
    break;

  case MODE_AGGREGATE:
    /* Pull Problem Size */
    nr=mxGetM(prhs[1]);
    nc=mxGetN(prhs[1]);

    /* Teuchos List */
    if(nrhs>2) List=build_teuchos_list(nrhs-2,&(prhs[2]));
    else List=new Teuchos::ParameterList;

    /* Pull Matrix - CSC Format */
    vals=mxGetPr(prhs[1]);
    if(rewrap_ints){
      colptr=mwIndex_to_int(nc+1,mxGetJc(prhs[1]));
      rowind=mwIndex_to_int(colptr[nc],mxGetIr(prhs[1]));
    }
    else{
      rowind=(int*)mxGetIr(prhs[1]);
      colptr=(int*)mxGetJc(prhs[1]);
    }

    /* Allocate space for aggregate / return value */
    plhs[0]=mxCreateNumericMatrix(nr,1,mxINT32_CLASS,mxREAL);
    agg=(int*)mxGetData(plhs[0]);

    /* Do aggregation only */
    mlmex_aggregate(nr,colptr,rowind,vals,List,agg);

    /* Cleanup */
    delete List;
    if(rewrap_ints){
      delete [] rowind;
      delete [] colptr;
    }

    break;
  default:
    mexErrMsgTxt("Error: Generic\n");
  };
}/*end mexFunction*/

#else
#error "Do not have MATLAB"
#endif

