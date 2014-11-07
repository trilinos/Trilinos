#ifndef BASKER_DEF_HPP
#define BASKER_DEF_HPP

#include "basker_decl.hpp"
#include "basker_scalartraits.hpp"
#include "basker.hpp"


//#include <assert.h>
#include <iostream>
#include <stdio.h>

using namespace std;

//#define BASKER_DEBUG 1
//#undef UDEBUG

namespace Basker{
   
  template <class Int, class Entry>
  Basker<Int, Entry>::Basker()
  {
    //A = new basker_matrix<Int,Entry>;
    A = (basker_matrix<Int,Entry> *) malloc(sizeof(basker_matrix<Int,Entry>));
    //L = new basker_matrix<Int,Entry>;
    L = (basker_matrix<Int,Entry> *) malloc(sizeof(basker_matrix<Int,Entry>));
    L->nnz = 0;
    //U = new basker_matrix<Int,Entry>;
    U = (basker_matrix<Int,Entry> *) malloc(sizeof(basker_matrix<Int,Entry>));
    U->nnz = 0;

    been_fact = false;
    perm_flag = false;
  }
 
 
  template <class Int, class Entry>
  Basker<Int, Entry>::Basker(Int nnzL, Int nnzU)
  {
    A = new basker_matrix<Int,Entry>;
    L = new basker_matrix<Int,Entry>;
    L->nnz = nnzL;
    U = new basker_matrix<Int,Entry>;
    U->nnz = nnzU;

    been_fact = false;
    perm_flag = false;
  }
 

  template <class Int, class Entry>
  Basker<Int, Entry>::~Basker()
  {
    //free factor
    if(been_fact)
      {
	free_factor();
	FREE(pinv);
      }
    if(perm_flag)
      {
	free_perm_matrix();
      }
    FREE(A);
    FREE(L);
    FREE(U);
    
  }


  template <class Int, class Entry>
  int Basker<Int,Entry>:: basker_dfs
  ( 
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
  )
  {
    Int have_elements = 1;
    Int i, t, i1, head ; 
    Int start, end, done, *store ;
   

    store = stack + n ;
    head = 0;
    stack[head] = j;
    bool has_elements = true;
     
    while(has_elements)
      {
	j = stack[head] ;
#ifdef BASKER_DEBUG
	std::cout << "DFS: " << j << "COLOR: " << color[j] << std::endl;
#endif
	t = pinv [j] ;
	if (color[j] == 0)
	  {
	    /* Seeing this column for first time */
	    color[j] = 1 ;
	    start = Lp[t] ;
	  }
	else
	  {
	    ASSERT (color[j] == 1) ; /* color cannot be 2 when we are here */
	    start = store[j];
	  }
	done = 1;

	if ( t != n )
	  {
	    end = Lp[t+1] ;
	    for ( i1 = start ; i1 < end ; i1++ )
	      {
		i = Li[i1] ;
		if ( color[i] == 0 )
		  {
		    stack[++head] = i;
		    store[j] = i1+1;
		    done = 0;
		    break;
		  }
	      }
	  }
	if (done)
	  {
	    pattern[--*top] = j ;
	    color[j] = 2 ;
	    if(head == 0)
	      {
		has_elements = false;
	      }
	    else
	      {
		head--;
	      }
	  }
      }
#ifdef BASKER_DEBUG
    std::cout << "Out of DFS: " << j << std::endl;
#endif
    return 0;
  } //End dfs



  template <class Int, class Entry>
  int Basker<Int,Entry>::factor(Int nrow, Int ncol , Int nnz, Int *col_ptr, Int *row_idx, Entry *val)
  {

    /*Initalize A basker matrix struc */
#ifdef  BASKER_DEBUG

    ASSERT(nrow > 0);
    ASSERT(ncol > 0);
    ASSERT(nnz > 0);

#endif

    A->nrow = nrow;
    A->ncol = ncol;
    A->nnz = nnz;
    A->col_ptr = col_ptr;
    A->row_idx = row_idx;
    A->val = val;
    /*End initalize A*/

    /*Creating space for L and U*/
    L->nrow = nrow;
    L->ncol = ncol;
    if(L->nnz == 0)
      {
	L->nnz = 2*A->nnz;
      }
    L->col_ptr = (Int *) CALLOC(ncol+1, sizeof(Int));
    L->row_idx = (Int *) CALLOC(L->nnz, sizeof(Int));
    L->val =     (Entry *) CALLOC(L->nnz, sizeof(Entry));

    U->nrow = nrow;
    U->ncol = ncol;
    if(U->nnz == 0)
      {
	U->nnz = 2*A->nnz;
      }
    U->col_ptr = (Int *) CALLOC(ncol+1, sizeof(Int));
    U->row_idx = (Int *) CALLOC(U->nnz, sizeof(Int));
    U->val =     (Entry *) CALLOC(U->nnz, sizeof(Entry));
    /*End creating space for L and U*/

    /*Creating working space*/
    Int *tptr;
    Entry *X;
    tptr = (Int *)   CALLOC( (ncol)+(4*nrow), sizeof(Int));
    X =    (Entry *) CALLOC(2*nrow, sizeof(Entry));
    pinv = (Int * )  CALLOC(ncol+1, sizeof(Int)); //Note extra pad
    /*End creating working space */
    
    /*Defining Variables Used*/
    Int i, j, k;
    Int *color, *pattern, *stack; // pointers into the work space
    Int top, top1, maxindex, t, j1, j2;
    Int lnnz, unnz, xnnz, lcnt, ucnt;
    Int cu_ltop, cu_utop;
    Int pp, p2, p;
    Int newsize;
    Entry pivot, value, xj;
    Entry absv, maxv;

    color = tptr;
    tptr += ncol;
    
    pattern = tptr;
    tptr += nrow;
    
    stack = tptr;
    tptr += 2*(nrow);

    
    cu_ltop = 0;
    cu_utop = 0;
    top =  ncol;
    top1 = ncol;
    lnnz = 0; //real found lnnz
    unnz = 0; //real found unnz

    for(k = 0 ; k < ncol; k++)
      {
	pinv[k] = ncol;      
      }

    /*For all columns in A .... */
    for (k = 0; k < ncol; k++)
      {

#ifdef BASKER_DEBUG
	cout << "k = " << k << endl;
#endif
	
	value = 0.0;
	pivot = 0.0;
	maxindex = ncol;
	j1 = 0;
	j2 = 0;
	lcnt = 0;
	ucnt = 0;

#ifdef BASKER_DEBUG
	ASSERT (top == ncol);
	
	for(i = 0; i < nrow; i++)
	  {
	    ASSERT(X[i] == 0);
	  }
	for(i = 0; i < ncol; i++)
	  {
	    ASSERT(color[i] ==0); 
	  }
#endif
	/* Reachability for every nonzero in Ak */
	for( i = col_ptr[k];  i < col_ptr[k+1]; i++)
	  {
	    j = row_idx[i];
	    X[j] = val[i];
	    
	    if(color[j] == 0)
	      {
		//do dfs
		basker_dfs(nrow, j, L->row_idx, L->col_ptr, color, pattern, &top, k, pinv, stack);

	      }

	  }//end reachable

	xnnz = ncol - top;
#ifdef BASKER_DEBUG
	cout << top << endl;
	cout << ncol << endl;
	cout << xnnz << endl;
	//ASSERT(xnnz <= nrow);
#endif
	/*Lx = b where x will be the column k in L and U*/
	top1 = top;
	for(pp = 0; pp < xnnz; pp++)
	  {
	    j = pattern[top1++];
	    color[j] = 0;
	    t = pinv[j];
	    
	    if(t!=ncol) //it has not been assigned
	      {
		xj = X[j];
		p2 = L->col_ptr[t+1];
		for(p = L->col_ptr[t]+1; p < p2; p++)
		  {
		    X[L->row_idx[p]] -= L->val[p] * xj;
		  }//over all rows
	      }

	  }

	/*get the pivot*/
	maxv = 0.0;
	for(i = top; i < nrow; i++)
	  {
	    j = pattern[i];
	    t = pinv[j];
	    value = X[j];
	    /*note may want to change this to traits*/
	    //absv = (value < 0.0 ? -value : value);
	    absv = BASKER_ScalarTraits<Entry>::abs(value);

	    if(t == ncol)
	      {
		lcnt++;
		if( absv > maxv)
		  {
		    maxv = absv;
		    pivot = value;
		    maxindex= j;
		  }
	      }
	  }
	ucnt = nrow - top -lcnt + 1;
	
	if(maxindex == ncol || pivot ==0)
	  {
	    cout << "Matrix is singular at index: " << maxindex << " pivot: " << pivot << endl;
	    return 1;
	  }

	pinv[maxindex] = k;
#ifdef BASKER_DEBUG
	if(maxindex != k )
	  {
	    cout << "Permuting pivot: " << k << " for row: " << maxindex << endl;
	  }
#endif	

	if(lnnz + lcnt >= L->nnz)
	  {
	    
	    newsize = L->nnz * 1.1 + 2*nrow + 1;
#ifdef BASKER_DEBUG
	    cout << "Out of memory -- Reallocating.  Old Size: " << L->nnz << " New Size: " << newsize << endl;
#endif
	    L->row_idx = (Int *) REALLOC(L->row_idx, newsize*sizeof(Int));
	    if(!(L->row_idx))
	      {
		cout << "WARNING: Cannot Realloc Memory" << endl;
		return 1;
	      }
	    L->val = (Entry *) REALLOC(L->val, newsize*sizeof(Entry));
	    if(!(L->val))
	      {
		cout << "WARNING: Cannot Realloc Memory" << endl;
		return 1;
	      }
	    L->nnz = newsize;
	    
	  }//realloc if L is out of memory
	
	if(unnz + ucnt >= U->nnz)
	  {
	    newsize = U->nnz*1.1 + 2*nrow + 1;
#ifdef BASKER_DEBUG
	    cout << "Out of memory -- Reallocating.  Old Size: " << L->nnz << " New Size: " << newsize << endl;
#endif
	    U->row_idx = (Int *) REALLOC(U->row_idx, newsize*sizeof(Int));
	    if(!(U->row_idx))
	      {
		cout << "WARNING: Cannot Realloc Memory" << endl;
		return 1;
	      }

	    U->val = (Entry *) REALLOC(U->val, newsize*sizeof(Entry));
	    if(!(U->val))
	      {
		cout << "WARNING: Cannot Realloc Memory" << endl;
		return 1;
	      }
	    U->nnz = newsize;
	  }//realloc if U is out of memory
	
	//L->col_ptr[lnnz] = maxindex;
	L->row_idx[lnnz] = maxindex;
	L->val[lnnz] = 1.0;
	lnnz++;

	for(i = top; i < nrow; i++)
	  {
	    j = pattern[i];
	    t = pinv[j];

	    /* check for numerical cancellations */
	    if(X[j] != 0)
	      {
		
		if(t != ncol)
		  {
		    if(unnz >= U->nnz)
		      {
			cout << "BASKER: Insufficent memroy for U" << endl;
			return 1;
		      }
		    U->row_idx[unnz] = pinv[j];
		    U->val[unnz] = X[j];
		    unnz++;
		    
		  }
		else if (t ==  ncol)
		  {
		    if(lnnz >= L->nnz)
		      {
			cout << "BASKER: Insufficent memroy for L" << endl;
			return 1;
		      }
		    
		    L->row_idx[lnnz]  = j;
		    //L->val[lnnz] = X[j]/pivot;
		    L->val[lnnz] = BASKER_ScalarTraits<Entry>::divide(X[j],pivot);
		    lnnz++;
		    
		  }

	      }
	    X[j] = 0;

	  }

	xnnz = 0;
	top = ncol;

	L->col_ptr[k] = cu_ltop;
	L->col_ptr[k+1] = lnnz;
	cu_ltop = lnnz;
	
	U->col_ptr[k] = cu_utop;
	U->col_ptr[k+1] = unnz;
	cu_utop = unnz;

      } //end for every column

#ifdef BASKER_DEBUG
    /*Print out found L and U*/
    for(k = 0; k < lnnz; k++)
      {
	printf("L[%d]=%g" , k , L->val[k]);	
      }
    cout << endl;
    for(k = 0; k < lnnz; k++)
      {
	printf("Li[%d]=%d", k, L->row_idx[k]);
      }
    cout << endl;
    for(k = 0; k < nrow; k++)
      {
	printf("p[%d]=%d", k, pinv[k]);
      }
    cout << endl;

#endif    
    /*  Repermute   */
    for( i = 0; i < ncol; i++)
      {
	for(k = L->col_ptr[i]; k < L->col_ptr[i+1]; k++)
	  {
	    L->row_idx[k] = pinv[L->row_idx[k]];
	  }
      }

#ifdef BASKER_DEBUG
    cout << "After Permuting" << endl;
    for(k = 0; k < lnnz; k++)
      {
	printf("Li[%d]=%d", k, L->row_idx[k]);
      }
    cout << endl;
#endif

    //FREE(X);
    //FREE(tptr);
    return 0;
  }//end factor

  template <class Int, class Entry>
  int Basker<Int, Entry>::returnL(Int *dim, Int *nnz, Int **col_ptr, Int **row_idx, Entry **val)
  {
    int i;
    *dim = L->nrow;
    *nnz = L->nnz;
    
    /*Does a bad copy*/
    
    *col_ptr = (Int *)   CALLOC(L->nrow+1, sizeof(Int));
    *row_idx = (Int *)   CALLOC(L->nnz, sizeof(Int));
    *val     = (Entry *) CALLOC(L->nnz, sizeof(Entry));
 
    
    for(i = 0; i < L->nrow+1; i++)
      {
	(*col_ptr)[i] = L->col_ptr[i];
      }
    
    for(i = 0; i < L->nnz; i++)
      {
	(*row_idx)[i] = L->row_idx[i];
	(*val)[i]     = L->val[i];
      }
    
  }

  template <class Int, class Entry>
  int Basker<Int, Entry>::returnU(Int *dim, Int *nnz, Int **col_ptr, Int **row_idx, Entry **val)
  {
    int i;
    *dim = U->nrow;
    *nnz = U->nnz;
    /*Does a bad copy*/
    *col_ptr = (Int *)   CALLOC(U->nrow+1, sizeof(Int));
    *row_idx = (Int *)   CALLOC(U->nnz, sizeof(Int));
    *val     = (Entry *) CALLOC(U->nnz, sizeof(Entry));
 
    for(i = 0; i < U->nrow+1; i++)
      {
	(*col_ptr)[i] = U->col_ptr[i];
      }
    for(i = 0; i < U->nnz; i++)
      {
	(*row_idx)[i] = U->row_idx[i];
	(*val)[i]     = U->val[i];
      }
  }

  template <class Int, class Entry>
  int Basker<Int, Entry>::returnP(Int** p)
  {
    Int i;
    *p = (Int *) CALLOC(A->nrow, sizeof(Int));
   
    for(i = 0; i < A->nrow; i++)
      {
	(*p)[pinv[i]] = i;  //Matlab perm-style
      }
  }
  
  template <class Int, class Entry>
  void Basker<Int, Entry>::free_factor()
  {
    //FREE L
    FREE(L->col_ptr);
    FREE(L->row_idx);
    FREE(L->val);
    

    //FREE U
    FREE(U->col_ptr);
    FREE(U->row_idx);
    FREE(U->val);
    
  }
  template <class Int, class Entry>
  void Basker<Int, Entry>::free_perm_matrix()
  {
    FREE(A->col_ptr);
    FREE(A->row_idx);
    FREE(A->val);
  }
  
  template <class Int, class Entry>
  int Basker<Int, Entry>::solve(Entry *b, Entry *x)
  {
    
    Int i;
    /*permute x to new ordering*/
    for(i = 0 ; i < A->ncol; i++) 
      {

		
	x[i] = b[pinv[i]];
	b[pinv[i]] = 0;
      }
    //ASSERT(L->nrow == L->ncol);
    low_tri_solve_csc(L->nrow, L->col_ptr, L->row_idx, L->val, b, x); 
    cout << " back solve " << endl;
    up_tri_solve_csc(U->nrow, U->col_ptr, U->row_idx, U->val, x, b);
    return 0;
  }
  
  template < class Int, class Entry>
  int Basker<Int, Entry>::low_tri_solve_csc( Int n, Int *col_ptr, Int *row_idx, Entry* val,  Entry *x, Entry *b)
  {
    Int i, j;
        /*for each column*/
    for(i = 0; i < n ; i++)
      {
	ASSERT(val[col_ptr[i]] != 0);
	
	x[i] = b[i]/val[col_ptr[i]]; //diag
	
	for(j = col_ptr[i]+1; j < (col_ptr[i+1]); j++) //update all rows
	  {
	    b[row_idx[j]] = b[row_idx[j]] - (val[j]*x[i]);
	  }
	
      }
    return 0;
  }

  template < class Int, class Entry>
  int Basker<Int, Entry>::up_tri_solve_csc( Int n, Int *col_ptr, Int *row_idx, Entry *val,  Entry *x, Entry *b)
  {
    Int i, j;
        /*for each column*/
    for(i = n; i > 0 ; i--)
      {
	int ii = i-1;
	ASSERT(val[col_ptr[i]-1] != 0);
			
	x[ii] = b[ii]/val[col_ptr[i]-1]; //diag
		
	for(j = (col_ptr[i]-2); j >= (col_ptr[ii]); j--)
	  {    
	    b[row_idx[j]] = b[row_idx[j]] - (val[j]*x[ii]);
	  }
      }
    return 0;
  }

  template <class Int, class Entry>
  int Basker<Int, Entry>::preorder(Int *row_perm, Int *col_perm)
  {
    
    basker_matrix <Int, Entry> *B;
    B = new basker_matrix<Int, Entry>;
    B->nrow = A->nrow;
    B->ncol = A->ncol;
    B->nnz = A->nnz;
    B->col_ptr = (Int *) CALLOC(A->ncol + 1, sizeof(Int));
    B->row_idx = (Int *) CALLOC(A->nnz, sizeof(Int));
    B->val     = (Entry *) CALLOC(A->val, sizeof(Int));

    int resultcol = permute_column(col_perm, B);
    int resultrow = permute_row(row_perm, B);

    /*Note: the csc matrices of A are the problem of the user
      therefore we will not free them*/
    A->col_ptr = B->col_ptr;
    A->row_idx = B->row_idx;
    A->val     = A->val;
    
    perm_flag = true; /*Now we will free A at the end*/

    return 0;
  }

  template <class Int, class Entry>
  int Basker <Int, Entry>::permute_column(Int *p, basker_matrix<Int,Entry> *B)
  {
    /*p(i) contains the destination of row i in the permuted matrix*/
    Int i,j, ii, jj;
    
    /*Determin column pointer of output matrix*/
    for(j=0; j < B->ncol; j++)
      {
	i = p[j];
	B->col_ptr[i+1] = A->col_ptr[j+1] - A->col_ptr[j];
      }
    /*get pointers from lengths*/
    B->col_ptr[0] = 0;
    for(j=0; j < B->ncol; j++)
      {
	B->col_ptr[j+1] = B->col_ptr[j+1] + B->col_ptr[j];
      }

    /*copy idxs*/
    Int k, ko;
    for(ii = 0 ; ii < B->ncol; ii++)
      {// old colum ii    new column p[ii]   k->pointer
	ko = B->col_ptr(p[ii]);
	for(k = A->col_ptr[ii]; k < A->col_ptr[ii+1]; k++)
	  {
	    B->row_index[ko] = A->row_index[k];
	    B->val[ko] = A->val[ko];
	    ko++;
	  }
      }
    return 0;
  }

  template <class Int, class Entry>
  int Basker <Int, Entry>::permute_row(Int *p,  basker_matrix<Int,Entry> *B)
  {
    Int k,i;
    for(k=0; k < A->nnz; k++)
      {
	B->row_idx[k] = p[A->row_idx[k]];
      }
    /*May want to sort values and index later*/
  }

}//end namespace
#endif  
  

