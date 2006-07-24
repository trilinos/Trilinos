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
   Initial Version 05/16/2006
   Teuchos-friendly version 05/22/2006
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
struct MLAPI_DATA_PACK{
  Space * FineSpace;
  DistributedMatrix * A;
  Teuchos::ParameterList *List;
  MultiLevelAdaptiveSA *Prec;  
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
  int numpdes=D->List->get("mlmex: numpdes",MLMEX_DEFAULT_NUMPDES);
  int cutoff=D->List->get("mlmex: levels",MLMEX_DEFAULT_LEVELS);
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
    if(nrhs>2&&mxIsSparse(prhs[1])&&mxIsNumeric(prhs[2])) rv=MODE_SOLVE;
    else mexErrMsgTxt("Error: Invalid input for solve\n");
    break;
  case MODE_CLEANUP:
    if(nrhs==1) rv=MODE_CLEANUP;
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
    case mxINT64_CLASS:
    case mxUINT64_CLASS:
#ifdef VERBOSE_OUTPUT
      mexPrintf("[%s] Int Found!\n",option_name);
#endif
      int *opt_int=(int*)mxGetData(prhs[i+1]);
      if(M==1 && N==1) TPL->set(option_name, opt_int[0]);      
      else TPL->set(option_name, opt_int);      
      break;      
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
   The first rhs arguement is always the program mode.  The other args are the
   ml.m parameter list in order (barring the mode arg).
*/

/* mexFunction is the gateway routine for the MEX-file. */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
  int *rowind,*colptr;
  double *b, *x, *vals,*opts_array;
  int nr, nc,no,iters,numpdes,cutoff,solver,adaptivevecs,ssweeps;
  bool UseDefaultNullSpace;
  double tol,agg_thresh;
  Teuchos::ParameterList* SolveOpts;
  mxClassID outclass;
  MODE_TYPE mode;
  static MLAPI_DATA_PACK *D;
  
  /* Sanity Check Input */
  mode=sanity_check(nrhs,prhs);
  
  switch(mode){
  case MODE_SETUP:
    D=new MLAPI_DATA_PACK;

    /* Pull Problem Size */
    nr=mxGetM(prhs[1]);
    nc=mxGetN(prhs[1]);

    if(nrhs>2 && mxIsChar(prhs[2]))
      /* Teuchos List*/
      D->List=build_teuchos_list(nrhs-2,&(prhs[2]));
    else{
      numpdes=MLMEX_DEFAULT_NUMPDES;
      cutoff=MLMEX_DEFAULT_LEVELS;
      adaptivevecs=MLMEX_DEFAULT_ADAPTIVEVECS;
      UseDefaultNullSpace=MLMEX_DEFAULT_USEDEFAULTNS;
      agg_thresh=0;
      ssweeps=3;
      /* Old-School Options List */
      if(nrhs==3 && mxIsNumeric(prhs[2])) {
        opts_array=mxGetPr(prhs[2]);
        no=MAX(mxGetM(prhs[2]),mxGetN(prhs[2]));
        numpdes=(int)opts_array[0];
        if(no>1) cutoff=(int)opts_array[1];
        if(no>2) adaptivevecs=(int)opts_array[2];
        if(no>3) agg_thresh=opts_array[3];
        if(no>4) ssweeps=(int)opts_array[4];
      }
      /* Build the teuchos list from old school / default */
      D->List=new Teuchos::ParameterList;
      D->List->set("additional candidates", adaptivevecs);
      D->List->set("smoother: sweeps",ssweeps);
      D->List->set("use default null space", UseDefaultNullSpace);
      D->List->set("aggregation: threshold", agg_thresh);
    }/*end else*/
    
    /* Pull Matrix - CSC Format */
    vals=mxGetPr(prhs[1]);
    rowind=mxGetIr(prhs[1]);
    colptr=mxGetJc(prhs[1]);

    /* Construct the Heirarchy */
    mlapi_data_pack_setup(D,nr,rowind,colptr,vals);
    
    /* Lock so we can keep the memory for the heirarchy */
    mexLock();
    break;

  case MODE_SOLVE:
    /* Pull Problem Size */
    nr=mxGetM(prhs[1]);
    nc=mxGetN(prhs[1]);
    
    /* Pull RHS */
    b=mxGetPr(prhs[2]);
    
    if(nrhs>3 && mxIsChar(prhs[3]))
      /* Teuchos List*/
      SolveOpts=build_teuchos_list(nrhs-3,&(prhs[3]));
    else{
      iters=5;tol=1e-13;solver=1;
      if(nrhs==4 && mxIsNumeric(prhs[3])) {
        opts_array=mxGetPr(prhs[3]);
        no=MAX(mxGetM(prhs[3]),mxGetN(prhs[3]));
        iters=(int)opts_array[0];
        if(no>1) tol=opts_array[1];
        if(no>2) solver=(int)opts_array[2];
      }/*end if*/
      
      /* Update the Teuchos List */
      SolveOpts=new Teuchos::ParameterList;
      SolveOpts->set("krylov: tolerance", tol);
      SolveOpts->set("krylov: max iterations", iters);
      if(solver==1) SolveOpts->set("krylov: type", "fixed point");
      else if(solver==2) SolveOpts->set("krylov: type", "cg");
      else if(solver==3) SolveOpts->set("krylov: type", "gmres");
    }/*end else*/

       
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
    /* Cleanup the memory mess */
    mlapi_data_pack_cleanup(D);

    /* Unlock and let stuff go away */
    mexUnlock();
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
