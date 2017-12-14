#include <mex.h>
#include <matrix.h>
#include <stdlib.h>
#include <string.h>

#include "basker_decl.hpp"
#include "basker_def.hpp"

/* TODO: Need a better way to do this with matlab.*/
#define Int mwIndex
//#define mwIndex long
/*
 * 
 * Finds the L and U such that A=L*U.
 * 
 */

void mexFunction
(
    int nlhs,
    mxArray *plhs [],
    int nrhs,
    const mxArray *prhs []
)
{
    /* Compressed column form of L, x and b */
    mwIndex *Lp, *Ap, *Ai, *Up, *pp, *pi ;
    Int *Li, *Ui;
    mwIndex *Lp1, *Li1, *Up1, *Ui1 ;
    double *Lx, *Ax, *Ux, *px ;
    double *Lx1, *Ux1 ;
    double *t1, *t2 ;

    mwIndex anrow ;
    mwIndex ancol ;
    mwIndex lnnz ;                       
    mwIndex unnz ;                       
    mwIndex i ;
    mwIndex j ;
    mwIndex app_xnnz ;
    mwIndex memsize ;
    mwIndex result ;

    mwIndex *ws;  /* workspace */
    double *X ;     /* space to scatter x */
    mwIndex *pinv ; /* column permutation inverse */

    if ( nlhs != 3 || nrhs < 3 )
    {
      //mexErrMsgTxt (" Incorrect number of arguments to sproductmex \n") ;
    }

    Ap = mxGetJc(prhs[0]) ;
    Ai = mxGetIr(prhs[0]) ;
    Ax = mxGetPr(prhs[0]) ;

    anrow = mxGetM (prhs[0]) ;
    ancol = mxGetN (prhs[0]) ;

    t1 = mxGetPr (prhs[1]) ;
    t2 = mxGetPr (prhs[2]) ;

    lnnz = (mwIndex)*t1 ;
    unnz = (mwIndex)*t2 ;
    /*printf("lnnz=%d unnz=%d\n", lnnz, unnz);*/

   
    BaskerClassicNS::BaskerClassic <mwIndex, double> mybasker;
    mybasker.factor(anrow, ancol, lnnz, Ap, Ai, Ax);
    mybasker.returnL(&anrow, &lnnz, &Lp, &Li, &Lx);
    mybasker.returnU(&anrow, &unnz, &Up, &Ui, &Ux);
    mybasker.returnP(&pp);
    

    plhs[0] = mxCreateSparse (anrow, ancol, lnnz+1, mxREAL) ;
    Lp1 = mxGetJc (plhs[0]) ;
    Li1 = mxGetIr (plhs[0]) ;
    Lx1 = mxGetPr (plhs[0]) ;

    plhs[1] = mxCreateSparse (anrow, ancol, unnz, mxREAL) ;
    Up1 = mxGetJc (plhs[1]) ;
    Ui1 = mxGetIr (plhs[1]) ;
    Ux1 = mxGetPr (plhs[1]) ;

    
    mwIndex *pp1, *pp2;
    double *ppx;
    plhs[2] = mxCreateSparse (ancol, ancol, ancol, mxREAL);
    pp1 = mxGetJc (plhs[2]);
    pp2 = mxGetIr (plhs[2]);
    ppx = mxGetPr (plhs[2]);
    
    
    Lp1[0] = Lp[0];
    for ( i = 0 ; i < ancol ; i++)
    {
        Lp1[i+1] = Lp[i+1];
        for ( j = Lp[i] ; j < Lp[i+1] ; j++ )
        {
            Li1[j] = Li[j];
            Lx1[j] = Lx[j];
        }
    } 

    Up1[0] = Up[0];
    for ( i = 0 ; i < ancol ; i++)
    {
        Up1[i+1] = Up[i+1];
        for ( j = Up[i] ; j < Up[i+1] ; j++ )
        {
            Ui1[j] = Ui[j];
            Ux1[j] = Ux[j];
        }
    } 
  
    
    //mexPrintf("Perm \n");
    for ( i = 0; i < ancol; i++)
      {

	//mexPrintf("%d ", pp[i]);
	pp1[i] = i;
	//pp2[i] = i;
	pp2[pp[i]] = i ;
	ppx[i] = 1;
      }
    pp1[ancol] = ancol;

   
    
    mxFree (pp) ;
    mxFree (Lp) ;
    mxFree (Li) ;
    mxFree (Lx) ;
    mxFree (Up) ;
    mxFree (Ui) ;
    mxFree (Ux) ;
    
}
