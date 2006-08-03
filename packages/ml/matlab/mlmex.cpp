//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/* This file provides a simple MEX interface to use ML in MATLAB.  This code is
   based off of AdaptiveSA.cpp by Marzio Sala, and is intended to be used with
   it's MATLAB front-end ml.m.

   By: Chris Siefert <csiefer@sandia.gov>
   Version History
   05/16/2006 - Initial Version 
   05/22/2006 - Teuchos-friendly version 
   08/03/2006 - Added functionality to allow ML to have multiple systems ready
                to be solved at once.  This entailed adding the system handle stuff.
*/

#include <stdio.h>
#include <string.h>

/* Needed for ML - Horked from the examples/AdaptiveSA.cpp.  Thanks, Marzio!*/
#include "ml_config.h"
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_MultiLevelAdaptiveSA.h"
#include "MLAPI_DistributedMatrix.h"
#include "MLAPI_Krylov.h"

#ifdef HAVE_ML_MATLAB
/* Needed for MEX */
#include "mex.h"

using namespace Teuchos;
using namespace MLAPI;
using namespace std;

extern void _main();

/* Macros */
#define MAX(x,y) ((x)>(y)?(x):(y))
#define FLOOR(x) ((int)(x))
#define ISINT(x) (((x-(int)(x))<1e-15*(x))?true:false)

#define IS_FALSE 0
#define IS_TRUE  1

/* Mode Info */
typedef enum {MODE_SETUP=0, MODE_SOLVE, MODE_CLEANUP} MODE_TYPE;

/* Default values */
#define MLMEX_DEFAULT_LEVELS 10
#define MLMEX_DEFAULT_NUMPDES 1
#define MLMEX_DEFAULT_ADAPTIVEVECS 0
#define MLMEX_DEFAULT_USEDEFAULTNS true

/* Debugging */
//#define VERBOSE_OUTPUT


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* The big struct of data */
struct MLAPI_DATA_PACK;
struct MLAPI_DATA_PACK{
  int id;
  Space * FineSpace;
  DistributedMatrix * A;
  Teuchos::ParameterList *List;
  MultiLevelAdaptiveSA *Prec;  
  struct MLAPI_DATA_PACK *next;
};


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
/* mlapi_data_pack_setup - This function does the setup phase for MLAPI, pulling
   key parameters from the Teuchos list, and calling the aggregation routines
   Parameters:
   D       - The MLAPI_DATA_PACK, with Teuchos list set. [I/O]
   N       - Number of unknowns [I]
   rowind  - Row indices of matrix (CSC format) [I]
   colptr  - Column indices of matrix (CSC format) [I]
   vals    - Nonzero values of matrix (CSC format) [I]
*/
void mlapi_data_pack_setup(MLAPI_DATA_PACK *D,int N,int* rowind,int* colptr, double* vals){

  int i,j;
  /* Initialize the workspace and set the output level */
  Init();  

  /* Define the space for fine level vectors and operators */
  D->FineSpace=new Space(N);
  
  /* Do the matrix assembly --- Someone should really come up with a better way
     of doing this.  Argh! */
  D->A=new DistributedMatrix(*(D->FineSpace),*(D->FineSpace)); 
  for(i=0;i<N;i++)
    for(j=colptr[i];j<colptr[i+1];j++)
      (*(D->A))(rowind[j],i)=vals[j];
  D->A->FillComplete();

  /* Pull key options from Teuchos */
  int numpdes=D->List->get("PDE equations",MLMEX_DEFAULT_NUMPDES);
  int cutoff=D->List->get("max levels",MLMEX_DEFAULT_LEVELS);
  int adaptivevecs=D->List->get("additional candidates",MLMEX_DEFAULT_ADAPTIVEVECS);
  bool UseDefaultNullSpace=D->List->get("use default null space", MLMEX_DEFAULT_USEDEFAULTNS);

  /* Allocate Memory */
  D->Prec=new MultiLevelAdaptiveSA(*(D->A), *(D->List),numpdes,cutoff);
    
  /* Build the Heirarchy */
  if(adaptivevecs>0) D->Prec->AdaptCompute(UseDefaultNullSpace,adaptivevecs);    
  else D->Prec->Compute();
  Finalize();
}/*end mlapi_data_pack_setup*/


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* mlapi_data_pack_cleanup - This function cleans up all the memory used during
   the setup phase.
*/
void mlapi_data_pack_cleanup(MLAPI_DATA_PACK *D){
  delete D->FineSpace;
  delete D->A;
  delete D->List;
  delete D->Prec;
}/*end mlapi_data_pack_cleanup*/


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* MLAPI data pack list */
class mlapi_data_pack_list;
class mlapi_data_pack_list{
public:
  mlapi_data_pack_list();
  ~mlapi_data_pack_list();

  /* add - Adds an MLAPI_DATA_PACK to the list.
     Parameters:
     D       - The MLAPI_DATA_PACK. [I]
     Returns: problem id number of D
  */
  int add(struct MLAPI_DATA_PACK *D);

  /* find - Finds problem by id
     Parameters:
     id      - ID number [I]
     Returns: pointer to MLAPI_DATA_PACK matching 'id', if found, NULL if not
     found.
  */
  struct MLAPI_DATA_PACK* find(int id);

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
  
protected:
  int num_probs;
  /* Note: This list is sorted */
  struct MLAPI_DATA_PACK *L;
};


mlapi_data_pack_list::mlapi_data_pack_list():num_probs(0),L(NULL){}
mlapi_data_pack_list::~mlapi_data_pack_list(){
  struct MLAPI_DATA_PACK *old;
  while(L!=NULL){
    old=L;
    mlapi_data_pack_cleanup(L);L=L->next;
    delete old;
  }/*end while*/
}/*end destructor*/

/* add - Adds an MLAPI_DATA_PACK to the list.
   Parameters:
   D       - The MLAPI_DATA_PACK. [I]
   Returns: problem id number of D
*/
int mlapi_data_pack_list::add(struct MLAPI_DATA_PACK *D){
  int idx_prev=0;
  struct MLAPI_DATA_PACK *iprev, *icurr;
  
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
   Returns: pointer to MLAPI_DATA_PACK matching 'id', if found, NULL if not
   found.
*/
struct MLAPI_DATA_PACK*  mlapi_data_pack_list::find(int id){
  struct MLAPI_DATA_PACK *rv;
  for(rv=L;rv && rv->id<id;rv=rv->next);
  if(rv && rv->id==id) return rv;
  else return NULL;  
}/*end find*/


/* remove - Removes problem by id
   Parameters:
   id      - ID number [I]
   Returns: IS_TRUE if remove was succesful, IS_FALSE otherwise
*/
int mlapi_data_pack_list::remove(int id){
  struct MLAPI_DATA_PACK *iprev, *icurr;
  if(!L) return IS_FALSE;
  for(iprev=NULL,icurr=L; icurr && icurr->id<id; iprev=icurr,icurr=icurr->next);

  if(!icurr || icurr->id!=id) return IS_FALSE;
  else{
    if(!iprev) L=icurr->next;
    else iprev->next=icurr->next;
    mlapi_data_pack_cleanup(icurr);
    free(icurr);
    return IS_TRUE;
  }/*end else*/    
}/*end remove*/


/* size - Number of stored problems
   Returns: num_probs
*/
int mlapi_data_pack_list::size(){return num_probs;}
/**************************************************************/
/**************************************************************/
/**************************************************************/
/* mlapi_solve - Given two Teuchos lists, one in the MLAPI_DATA_PACK, and one of
   solve-time options, this routine calls the relevant solver and returns the solution.
   Parameters:
   D       - The MLAPI_DATA_PACK, with Teuchos list set. [I]
   TPL     - Teuchos list of solve-time options [I]
   N       - Number of unknowns [I]
   b       - RHS vector [I]
   x       - solution vector [O]
*/
void mlapi_solve(MLAPI_DATA_PACK *D,Teuchos::ParameterList *TPL, int N, double*b, double*x){
  int i;
  Init();   
  MultiVector LHS(D->A->GetDomainSpace());
  MultiVector RHS(D->A->GetRangeSpace());

  /* Create the new Teuchos List */
  Teuchos::ParameterList tmp=*(D->List);
  tmp.setParameters(*TPL);
  
  /* Fill LHS/RHS */
  LHS=0;
  for(i=0;i<N;i++) RHS(i)=b[i];

#ifdef VERBOSE_OUTPUT
  tmp.print(cout,0,true);
#endif
  
  /* Do the solve */
  Krylov(*(D->A), LHS, RHS, *(D->Prec), tmp);

  /* Fill Solution */
  for(i=0;i<N;i++) x[i]=LHS(i);
  
  Finalize();
} /*end mlapi_solve*/




/**************************************************************/
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
  MODE_TYPE rv;
  double *modes;
  /* Check for mode */
  if(nrhs==0)
    mexErrMsgTxt("Error: Invalid Inputs\n");
  
  /* Pull mode data from 1st Input */
  modes=mxGetPr(prhs[0]);

  switch ((MODE_TYPE)modes[0]){
  case MODE_SETUP:
    if(nrhs>1&&mxIsSparse(prhs[1])) rv=MODE_SETUP;
    else mexErrMsgTxt("Error: Invalid input for setup\n");    
    break;
  case MODE_SOLVE:
    if(nrhs>3&&mxIsNumeric(prhs[1])&&mxIsSparse(prhs[2])&&mxIsNumeric(prhs[3])) rv=MODE_SOLVE;
    else mexErrMsgTxt("Error: Invalid input for solve\n");
    break;
  case MODE_CLEANUP:
    if(nrhs==1 || nrhs==2) rv=MODE_CLEANUP;
    else mexErrMsgTxt("Error: Extraneous args for cleanup\n");
    break;
  default:
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
Teuchos::ParameterList* build_teuchos_list(int nrhs,const mxArray *prhs[]){
  Teuchos::ParameterList* TPL=new Teuchos::ParameterList;
  mxClassID cid;
  int i,M,N;
  char * option_name;
  string opt_str;  
  for(i=0;i<nrhs;i+=2){
    if(i==nrhs-1 || !mxIsChar(prhs[i]))
      mexErrMsgTxt("Error: Input options are not in ['parameter',value] format!\n");

    /* What option are we setting? */
    option_name=mxArrayToString(prhs[i]);

    /* Pull relevant info the the option value */
    cid=mxGetClassID(prhs[i+1]);
    M=mxGetM(prhs[i+1]);
    N=mxGetN(prhs[i+1]);
#ifdef VERBOSE_OUTPUT
    mexPrintf("[%d] M=%d N=%d\n",i,M,N);
#endif

    /* Add to the Teuchos list */
    switch(cid){
    case mxCHAR_CLASS:
      char *opt_char=mxArrayToString(prhs[i+1]);
      opt_str=opt_char;
#ifdef VERBOSE_OUTPUT
      mexPrintf("[%s] String Found: %s\n",option_name,opt_char);
#endif
      TPL->set(option_name, opt_str);
      mxFree(opt_char);
      break;
    case mxDOUBLE_CLASS:
    case mxSINGLE_CLASS:
      //NTS: Does not deal with complex args
      double *opt_float=mxGetPr(prhs[i+1]);
      if(M==1 && N==1 && ISINT(opt_float[0])) {     
#ifdef VERBOSE_OUTPUT
        mexPrintf("[%s] Float(Int) Found!\n",option_name);
#endif
        TPL->set(option_name, (int)opt_float[0]);
      }/*end if*/
      else if(M==1 && N==1){
#ifdef VERBOSE_OUTPUT
        mexPrintf("[%s] Float Found!\n",option_name);
#endif
        TPL->set(option_name, opt_float[0]);
      }/*end if*/
      else{
#ifdef VERBOSE_OUTPUT
        mexPrintf("[%s] Float Found!\n",option_name);
#endif
        TPL->set(option_name, opt_float);
      }/*end else*/
      break;
    case mxLOGICAL_CLASS:
#ifdef VERBOSE_OUTPUT      
      mexPrintf("[%s] Logical Found!\n",option_name);
#endif
      if(M==1 && N==1) TPL->set(option_name, mxIsLogicalScalarTrue(prhs[i+1]));
      else TPL->set(option_name,mxGetLogicals(prhs[i+1]));
      //NTS: The else probably doesn't work.
      break;
    case mxINT8_CLASS:
    case mxUINT8_CLASS:
    case mxINT16_CLASS:
    case mxUINT16_CLASS:
    case mxINT32_CLASS:
    case mxUINT32_CLASS:

#ifdef VERBOSE_OUTPUT
      mexPrintf("[%s] Int Found!\n",option_name);
#endif
      int *opt_int=(int*)mxGetData(prhs[i+1]);
      if(M==1 && N==1) TPL->set(option_name, opt_int[0]);      
      else TPL->set(option_name, opt_int);      
      break;
      // NTS: 64-bit ints will break on a 32-bit machine.  We
      // should probably detect machine type, or somthing, but that would
      // involve a non-trivial quantity of autoconf kung fu.
    case mxINT64_CLASS:
    case mxUINT64_CLASS:      
    case mxFUNCTION_CLASS:
    case mxUNKNOWN_CLASS:
    case mxCELL_CLASS:
    case mxSTRUCT_CLASS:      
    default:
      mexPrintf("Error parsing input option #%d: %s [type=%d]\n",i,option_name,cid);
      mexErrMsgTxt("Error: An input option is invalid!\n");      
    };        

    /* Free memory */
    mxFree(option_name);   
  }/*end for*/  

  /* Return */
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
  int i,*id,sz,rv,*rowind,*colptr;
  double *b, *x, *vals,*opts_array;
  int nr, nc,no;
  bool UseDefaultNullSpace;
  double tol,agg_thresh;
  Teuchos::ParameterList* SolveOpts;
  mxClassID outclass;
  MODE_TYPE mode;
  static mlapi_data_pack_list* PROBS=NULL;
  MLAPI_DATA_PACK *D=NULL;
  
  /* Sanity Check Input */
  mode=sanity_check(nrhs,prhs);
  
  switch(mode){
  case MODE_SETUP:
    if(!PROBS) PROBS=new mlapi_data_pack_list;
    D=new MLAPI_DATA_PACK;

    /* Pull Problem Size */
    nr=mxGetM(prhs[1]);
    nc=mxGetN(prhs[1]);

    /* Teuchos List*/
    if(nrhs>2) D->List=build_teuchos_list(nrhs-2,&(prhs[2]));
    else D->List=new Teuchos::ParameterList;
    
    /* Pull Matrix - CSC Format */
    vals=mxGetPr(prhs[1]);
    rowind=mxGetIr(prhs[1]);
    colptr=mxGetJc(prhs[1]);

    /* Construct the Heirarchy */
    mlapi_data_pack_setup(D,nr,rowind,colptr,vals);

    /* Add this problem to the list */
    rv=PROBS->add(D);

    /* Set return value */
    plhs[0]=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    id=(int*)mxGetData(plhs[0]);id[0]=rv;
    
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
    nc=mxGetN(prhs[2]);
    
    /* Pull RHS */
    b=mxGetPr(prhs[3]);

    /* Teuchos List*/
    if(nrhs>4) SolveOpts=build_teuchos_list(nrhs-4,&(prhs[4]));
    else SolveOpts=new Teuchos::ParameterList;
    
    /* Allocate Solution Space */
    plhs[0]=mxCreateDoubleMatrix(nr,1,mxREAL);
    x=mxGetPr(plhs[0]);
    
    /* Sanity Check Matrix / RHS */
    if(nr != nc || nr != mxGetM(prhs[2]))
      mexErrMsgTxt("Error: Size Mismatch in Input\n");

    /* Run Solver */  
    mlapi_solve(D,SolveOpts,nr,b,x);

    /* Cleanup */
    delete SolveOpts;
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
  default:
    mexErrMsgTxt("Error: Generic\n");
  };
  
}/*end mexFunction*/

#else
#error "Do not have MATLAB"
#endif
#else
#error "Do not have MLAPI"
#endif
