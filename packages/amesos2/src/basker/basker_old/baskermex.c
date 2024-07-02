// @HEADER
// *****************************************************************************
//                   Basker: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the Basker contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <mex.h>
#include <matrix.h>
#include <stdlib.h>
#include <string.h>
#include "basker.h"

/* TODO: Need a better way to do this with matlab.*/
#define Int long

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
        mexErrMsgTxt (" Incorrect number of arguments to sproductmex \n") ;
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

    /* O(n) initialization */
    ws = (mwIndex *) mxCalloc ( (ancol)+(4*anrow), sizeof(mwIndex)) ;
    X = (double *) mxCalloc ( 2*anrow, sizeof(double)) ;
    pinv = (mwIndex *) mxCalloc ( ancol, sizeof(mwIndex)) ;

    Lp = (mwIndex *) mxCalloc ( ancol+1, sizeof(mwIndex)) ;
    Li = (mwIndex *) mxCalloc ( lnnz, sizeof(Int)) ;
    Lx = (double *) mxCalloc ( lnnz, sizeof(double)) ;

    Up = (mwIndex *) mxCalloc ( ancol+1, sizeof(mwIndex)) ;
    Ui = (mwIndex *) mxCalloc ( unnz, sizeof(Int)) ;
    Ux = (double *) mxCalloc ( unnz, sizeof(double)) ;

    PRINT(("Calling basker \n"));
    result = basker_basker_l(Ap, Ai, Ax, anrow, ancol, ws , X,
                Lp, &Li, &Lx, Up, &Ui, &Ux, &lnnz, &unnz, pinv) ;
    PRINT(("Back in mex function %d \n",Lp[ancol]));

    if (result)
    {
        mxFree (X) ;
        mxFree (ws) ;
        mxFree (Lp) ;
        mxFree (Li) ;
        mxFree (Lx) ;
        mxFree (Up) ;
        mxFree (Ui) ;
        mxFree (Ux) ;
        mexErrMsgTxt (" basker failed \n") ;
    }

    plhs[0] = mxCreateSparse (anrow, ancol, Lp[ancol], mxREAL) ;
    Lp1 = mxGetJc (plhs[0]) ;
    Li1 = mxGetIr (plhs[0]) ;
    Lx1 = mxGetPr (plhs[0]) ;

    plhs[1] = mxCreateSparse (anrow, ancol, Up[ancol], mxREAL) ;
    Up1 = mxGetJc (plhs[1]) ;
    Ui1 = mxGetIr (plhs[1]) ;
    Ux1 = mxGetPr (plhs[1]) ;

    plhs[2] = mxCreateSparse (ancol, ancol, ancol, mxREAL) ;
    pp = mxGetJc (plhs[2]) ;
    pi = mxGetIr (plhs[2]) ;
    px = mxGetPr (plhs[2]) ;

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

    for ( i = 0 ; i < ancol ; i++)
    {
      pp[i] = i ;
      pi[ pinv[i] ] = i ;
      px[i] = 1 ;
    }
    pp[anrow] = ancol ;

    /*for ( i = 0 ; i < ancol ; i++)
        printf("pinv[%d] = %d\n", i, pinv[i]);*/

    for ( i = 0 ; i < ancol ; i++)
    {
        for ( j = Lp1[i] ; j < Lp1[i+1] ; j++ )
        {
            Li1[j] = pinv[Li1[j]] ;
        }
    }


    mxFree (X) ;
    mxFree (ws) ;
    mxFree (pinv) ;
    mxFree (Lp) ;
    mxFree (Li) ;
    mxFree (Lx) ;
    mxFree (Up) ;
    mxFree (Ui) ;
    mxFree (Ux) ;
}
