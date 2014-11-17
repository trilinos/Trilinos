#ifndef BASKER_DECL_HPP
#define BASKER_DECL_HPP

#include "basker_types.hpp"

namespace Basker{

  template <class Int, class Entry>
  class Basker
  {
    
  public:
    Basker();
    Basker(Int nnzL, Int nnzU);
    ~Basker();
    int preorder(Int *row_perm, Int *col_perm);
    int factor(Int nrow, Int ncol , Int nnz, Int *col_ptr, Int *row_idx, Entry *val);
    int returnL(Int *dim, Int *nnz, Int **col_ptr, Int **row_idx, Entry **val);
    int returnU(Int *dim, Int *nnz, Int **col_ptr, Int **row_idx, Entry **val);
    int returnP(Int **p);
    int solve( Entry* b, Entry* x);
    int solveMultiple(Int nrhs, Entry *b, Entry *x);
    //int solve();
    
  private:
    int basker_dfs(
		   Int n,
		   Int j, 
		   Int *Li, 
		   Int *Lp, 
		   Int *color, 
		   Int *pattern, /* o/p */
		   Int *top,       /* o/p */ 
		   Int k,
		   Int *pinv,
		   Int *stack
		   );
    void free_factor();
    void free_perm_matrix();
    int low_tri_solve_csc(Int n, Int* col_ptr, Int *row_idx, Entry *val,  Entry *x,  Entry *b); 
    int up_tri_solve_csc(Int n, Int* col_ptr, Int *row_idx, Entry *val, Entry *x, Entry *b); 
    int permute_row(Int *p, basker_matrix<Int,Entry> *B);
    int permute_column(Int *p, basker_matrix<Int, Entry> *B);
    int sort_factors();
    basker_matrix<Int, Entry> *A;
    basker_matrix<Int, Entry> *L;
    basker_matrix<Int, Entry> *U;
    Int *in_perm;
    Int *pinv;
    bool been_fact;
    bool perm_flag;

  };

}/*End namespace*/
#endif
