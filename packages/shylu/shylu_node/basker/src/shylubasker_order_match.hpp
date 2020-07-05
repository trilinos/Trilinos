#ifndef SHYLUBASKER_ORDER_MATCH_HPP
#define SHYLUBASKER_ORDER_MATCH_HPP


//Include statements
#if defined(HAVE_AMESOS2_SUPERLUDIST)
  #define BASKER_MC64
#endif
#ifdef BASKER_MC64
//#include "mc64ad.hpp"
/* Note, come back and rewrite these 
 * function yourself at somepoint to clean up
 * In the mean time to test how well this works, we will use the auto rewritten from SuperLU_dist with modifications
 */

extern  "C"
{
  void mc64id_dist(int *);

  void mc64ad_dist(int*, int*, int*, int*, int*,  double*,
                   int*, int*, int*, int*, int*,  double*,
                   int*, int*);
}
#endif

#include "mwm2.hpp"



namespace BaskerNS
{
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::mwm(BASKER_MATRIX &M, INT_1DARRAY perm)
  {
    //DEBUG put to get compiled
    
    Int num = 0;
    mwm_order::mwm(M.nrow, M.nnz, 
		   &(M.col_ptr[0]), &(M.row_idx[0]),
		   &(M.val[0]), &(perm[0]), 
		   num);

    return 0;
  }//end mwm


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::mc64(BASKER_MATRIX &M, Int _job, INT_1DARRAY _perm,
                                         ENTRY_1DARRAY _scale_row, ENTRY_1DARRAY _scale_col)
  {
    //Note using primative types to match fortran
  #ifdef BASKER_MC64
    typedef     int         int_t;
  #else
    typedef     Int         int_t;
  #endif
    typedef     double      entry_t;

    int_t i, liw, ldw, num;
    int_t *iw, icntl[10], info[10];
    int_t job  = _job;
    entry_t *dw;
    entry_t *nzval_abs;
    int_t n = M.nrow;
    int_t nnz = M.nnz;

    nzval_abs = (entry_t*)malloc(M.nnz*sizeof(entry_t));

    liw = 5*n;
    if(job == 3) 
    { liw = 10*M.nrow + M.nnz; }
    iw = (int_t*) malloc(liw*sizeof(int_t));
    ldw = 3*n+nnz;
    dw = (entry_t*) malloc(ldw*sizeof(entry_t));

    //Convert to 1 formatting
    for(Int i = 0; i <= M.nrow; ++i)
      M.col_ptr[i] = M.col_ptr[i]+1;
    for(Int i = 0; i < M.nnz; ++i)
      M.row_idx[i] = M.row_idx[i]+1;
    
    //call init
    #ifdef BASKER_MC64
    mc64id_dist(icntl);
    #endif

    for(Int i = 0; i < M.nnz; ++i)
    { nzval_abs[i] = abs(M.val[i]); }

    Int* colptr = &(M.col_ptr[0]);
    Int* rowidx = &(M.row_idx[0]);
    Entry* val  = &(M.val[0]);

    Int *perm;
    perm = (Int*) malloc(M.nrow*sizeof(Int));

    #ifdef BASKER_MC64
    mc64ad_dist(&job, &n, &nnz, colptr, rowidx, nzval_abs,
	        &num, perm, &liw, iw, &ldw, dw, icntl, info);
    #endif

    //debug

    //convert indexing back
    for(Int i=0; i <= M.nrow; ++i)
    { M.col_ptr[i] = M.col_ptr[i]-1; }
    for(Int i=0; i < M.nnz; ++i)
    { M.row_idx[i] = M.row_idx[i]-1; }
    for(Int i=0; i < M.nrow; ++i)
    { _perm[i] = perm[i] -1; }
    for(Int i =0; i < A.nrow; ++i)
    { _scale_row[i] = exp (dw[i]); }
    for(Int i =0; i < A.nrow; ++i)
    { _scale_col[i] = exp (dw[i+A.nrow]); }

    //add job 5 special 

    free(iw);
    free(dw);
    free(nzval_abs);

    return 0;

  }//end mc64
  

  //With default matirix
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::mc64(Int _job, INT_1DARRAY _perm, ENTRY_1DARRAY _scale_row, ENTRY_1DARRAY _scale_col)
  {
    //Note using primative types to match fortran
    #ifdef BASKER_MC64
    typedef     int         int_t;
    #else
    typedef     Int         int_t;
    #endif
    typedef     double      entry_t;

    #ifdef BASKER_MC64
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
    if(job == 3) 
      { liw = 10*A.nrow + A.nnz; }
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
    mc64id_dist(icntl);
    #endif

    for(Int i = 0; i < A.nnz; ++i)
    { nzval_abs[i] = abs(A.val[i]); }
    
    Int *perm;
    perm = (Int*) malloc(A.nrow*sizeof(Int));
   
    #ifdef BASKER_MC64
    Int* colptr = &(A.col_ptr[0]);
    Int* rowidx = &(A.row_idx[0]);
    mc64ad_dist(&job, &n, &nnz, colptr, rowidx, nzval_abs,
	        &num, perm,
                &liw, iw, &ldw, dw, icntl, info);
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
