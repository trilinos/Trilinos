// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_ORDER_MATCH_HPP
#define SHYLUBASKER_ORDER_MATCH_HPP


//Include statements
//#include "mc64ad.hpp"
/* Note, come back and rewrite these 
 * function yourself at somepoint to clean up
 * In the mean time to test how well this works, we will use the auto rewritten from SuperLU_dist with modifications
 */

#include "mwm2.hpp"
//#if defined(HAVE_AMESOS2_SUPERLUDIST) && !defined(BASKER_MC64)
//  #define BASKER_SUPERLUDIS_MC64
//#endif
#ifdef BASKER_MC64
#include "mc64ad.hpp"
#elif defined(BASKER_SUPERLUDIS_MC64)
// FIX this
extern  "C"
{
  void mc64id_dist(int *);

  void mc64ad_dist(int*, int*, int*, int*, int*,  double*,
                   int*, int*, int*, int*, int*,  double*,
                   int*, int*);
}
#endif


namespace BaskerNS
{
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::mwm(BASKER_MATRIX &M, INT_1DARRAY perm)
  {
    Int num = 0;
    mwm_order::mwm(M.nrow, M.nnz, 
		   &(M.col_ptr[0]), &(M.row_idx[0]),
		   &(M.val[0]), &(perm[0]), 
		   num);
    return num;
  }//end mwm


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::mc64(Int n_, Int nnz_, Int *colptr, Int *row_idx, Entry *val,
                                         Int _job, Int *_perm, Entry *_scale_row, Entry *_scale_col)
  {
    //Note using primative types to match fortran
  #ifdef BASKER_SUPERLUDIS_MC64
    typedef     int         int_t;
  #else
    typedef     Int         int_t;
  #endif
    typedef     double      entry_t;

    int_t liw, ldw, num;
    int_t *iw;
    #if defined(BASKER_MC64) || defined(BASKER_SUPERLUDIS_MC64)
    int_t icntl[10], info[10];
    #endif
    int_t job  = _job;
    entry_t *dw;
    entry_t *nzval_abs;
    int_t n = n_;
    int_t nnz = nnz_;

    nzval_abs = (entry_t*)malloc(nnz*sizeof(entry_t));

    liw = 5*n;
    if(job == 3) 
    { liw = 10*n + n; }
    iw = (int_t*) malloc(liw*sizeof(int_t));
    ldw = 3*n+nnz;
    dw = (entry_t*) malloc(ldw*sizeof(entry_t));

    //Convert to 1 formatting
    for(Int i = 0; i <= n; ++i)
      colptr[i] = colptr[i]+1;
    for(Int i = 0; i < nnz; ++i)
      row_idx[i] = row_idx[i]+1;
    
    //call init
    #ifdef BASKER_MC64
    mc64id_(icntl);
    #elif defined(BASKER_SUPERLUDIS_MC64)
    mc64id_dist(icntl);
    #endif

    for(Int i = 0; i < nnz; ++i)
    { nzval_abs[i] = abs(val[i]); }

    Int *perm;
    perm = (Int*) malloc(n*sizeof(Int));

    #ifdef BASKER_MC64
    mc64ad_(&job, &n, &nnz, colptr, row_idx, nzval_abs,
	    &num, perm, &liw, iw, &ldw, dw, icntl, info);
    #elif defined(BASKER_SUPERLUDIS_MC64)
    mc64ad_dist(&job, &n, &nnz, colptr, row_idx, nzval_abs,
	        &num, perm, &liw, iw, &ldw, dw, icntl, info);
    #endif

    //debug

    //convert indexing back
    for(Int i=0; i <= n; ++i)
    { colptr[i] = colptr[i]-1; }
    for(Int i=0; i < nnz; ++i)
    { row_idx[i] = row_idx[i]-1; }
    for(Int i=0; i < n; ++i)
    { _perm[i] = perm[i]-1; }
    for(Int i =0; i < n; ++i)
    { _scale_row[i] = exp (dw[i]); }
    for(Int i =0; i < n; ++i)
    { _scale_col[i] = exp (dw[i+n]); }

    //add job 5 special 

    free(iw);
    free(dw);
    free(nzval_abs);

    return 0;

  }//end mc64
  
  // With user-specified matrix
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::mc64(BASKER_MATRIX &M, Int _job, INT_1DARRAY _perm,
                                         ENTRY_1DARRAY _scale_row, ENTRY_1DARRAY _scale_col)
  {
    int info;
    info = mc64(M.nrows, M.nnz, &(M.colptr(0)), &(M.rowidx(0)), &(M.val(0)),
                _job, &(_perm(0)), &(_scale_row(0)), &(_scale_col(0)));
    return info;
  }


  //With default matirix
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::mc64(Int _job, INT_1DARRAY _perm, ENTRY_1DARRAY _scale_row, ENTRY_1DARRAY _scale_col)
  {
    //Note using primative types to match fortran
    #ifdef BASKER_SUPERLUDIS_MC64
    typedef     int         int_t;
    #else
    typedef     Int         int_t;
    #endif
    typedef     double      entry_t;

    #if defined(BASKER_MC64) || defined(BASKER_SUPERLUDIS_MC64)
    int_t num;
    int_t icntl[10], info[10];
    #endif
    int_t liw, ldw;
    int_t *iw;
    int_t job  = _job;
    entry_t *dw;
    entry_t *nzval_abs;
    int_t n = A.nrow;
    int_t nnz = A.nnz;

    nzval_abs = (entry_t*)malloc(A.nnz*sizeof(entry_t));

    liw = 5*n;
    if(job == 3) {
      liw = 10*A.nrow + A.nnz;
    }
    iw = (int_t*) malloc(liw*sizeof(int_t));
    ldw = 3*n+nnz;
    dw = (entry_t*) malloc(ldw*sizeof(entry_t));

    //Convert to 1 formatting
    for(Int i = 0; i <= A.nrow; ++i)
      A.col_ptr[i] = A.col_ptr[i]+1;
    for(Int i = 0; i< A.nnz; ++i)
      A.row_idx[i] = A.row_idx[i]+1;
    
    //call init
    #ifdef BASKER_MC64
    mc64id_(icntl);
    #elif defined(BASKER_SUPERLUDIS_MC64)
    mc64id_dist(icntl);
    #endif

    for(Int i = 0; i < A.nnz; ++i)
    { nzval_abs[i] = abs(A.val[i]); }
    
    Int *perm;
    perm = (Int*) malloc(A.nrow*sizeof(Int));
   
    #if defined(BASKER_MC64) || defined(BASKER_SUPERLUDIS_MC64)
    Int* colptr = &(A.col_ptr[0]);
    Int* rowidx = &(A.row_idx[0]);
    #ifdef BASKER_MC64
    mc64ad_(&job, &n, &nnz, colptr, rowidx, nzval_abs,
	    &num, perm,
            &liw, iw, &ldw, dw, icntl, info);
    #elif defined(BASKER_SUPERLUDIS_MC64)
    mc64ad_dist(&job, &n, &nnz, colptr, rowidx, nzval_abs,
	        &num, perm,
                &liw, iw, &ldw, dw, icntl, info);
    #endif
    #endif

    //debug

    //convert indexing back
    for(Int i =0; i <=A.nrow; ++i)
      A.col_ptr[i] = A.col_ptr[i]-1;
    for(Int i =0; i < A.nnz; ++i)
      A.row_idx[i] = A.row_idx[i]-1;
    for(Int i =0; i < A.nrow; ++i)
      _perm[i] = perm[i] -1;
    for(Int i =0; i < A.nrow; ++i)
      _scale_row[i] = exp (dw[i]);
    for(Int i =0; i < A.nrow; ++i)
      _scale_col[i] = exp (dw[i+A.nrow]);

    free(iw);
    free(dw);
    free(nzval_abs);

    return 0;

  }//end mc64
  
 
}//end namespace BaskerNS

#endif //end ifndef BASKER_ORDERING_MATCHING_HPP
