#ifndef BASKER_ORDER_AMD_HPP
#define BASKER_ORDER_AMD_HPP

//AMD is amazing on circuit problems in a way
//that ND can't

//This can be done to user the smaller domains

#ifdef HAVE_AMESOS
#include "amesos_amd.h"
#include "amesos_colamd.h"
#include "amesos_ccolamd.h"
#endif

namespace BaskerNS
{

   //===============================AMD===================
  
  template <class Int>
  BASKER_FINLINE
  int amesos_amd
  (
   Int n,
   Int *Ap, 
   Int *Ai, 
   Int *p, 
   double *Control,
   double *Info
   )
  {
    return -1;
  }//end amesos_amd()
    
  template <>
  BASKER_FINLINE
  int amesos_amd<>
  (
   int n,
   int *Ap,
   int *Ai,
   int *p,
   double *Control,
   double *Info
   )
  {
    amesos_amd_order(n,Ap,Ai,p,Control,Info);
    return 0;
  }//end amesos_amd<int>
 

  template <>
  BASKER_FINLINE
  int amesos_amd<>
  (
   long n,
   long  *Ap,
   long *Ai,
   long *p,
   double   *Control,
   double   *Info
   )
  {
    amesos_amd_l_order(n,Ap,Ai,p,Control,Info);
    return 0;
  }//end amesos_amd<long int>
  

    //==========================csymamd===================

  template <class Int>
  BASKER_FINLINE
  int my_amesos_csymamd
  (
   Int n, 
   Int *Ap,
   Int *Ai,
   Int *p, 
   Int *cmember
   )
  {
    return -1;
  }//end my_amesos_csymamd


  template <>
  BASKER_FINLINE
  int my_amesos_csymamd <>
  (
   int n, 
   int *Ap,
   int *Ai,
   int *p, 
   int *cmember
   )
  {
    
    double knobs[CCOLAMD_KNOBS];
    int    stats[CCOLAMD_STATS];

    //use default knob settings
    amesos_ccolamd_set_defaults(knobs);
    knobs[0] = 10;
    knobs[1] = 0;
    knobs[2] = 2;

    amesos_csymamd(n, Ai, Ap,  p, knobs, stats, 
		     &(calloc), &(free), 
		     cmember, 0);

    amesos_csymamd_report(stats);
    
    return 0;
  }


  template <>
  BASKER_FINLINE
  int my_amesos_csymamd <>
  (
   long n, 
   long *Ap,
   long *Ai,
   long *p, 
   long *cmember
   )
  {
    
    double knobs[CCOLAMD_KNOBS];
    long    stats[CCOLAMD_STATS];

    //use default knob settings
    amesos_ccolamd_l_set_defaults(knobs);
    knobs[0] = 10;
    knobs[1] = 0;
    knobs[2] = 2;

    amesos_csymamd_l(n, Ai, Ap,  p, knobs, stats, 
		     &(calloc), &(free), 
		     cmember, 0);

    amesos_csymamd_l_report(stats);
    
    return 0;
  }






  template <class Int, class Entry, class Exe_Space>
  BASKER_FINLINE
  void Basker<Int,Entry,Exe_Space>::amd_order
  (
   BASKER_MATRIX &M,
   INT_1DARRAY   p
   )
  {

    double amd_info[AMD_INFO];
    amesos_amd(M.ncol, &(M.col_ptr(0)), 
	       &(M.row_idx(0)), &(p(0)),
	       NULL, amd_info);


  }//end amd_order()
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_FINLINE
  void Basker<Int, Entry,Exe_Space>::csymamd_order
  (
   BASKER_MATRIX &M,
   INT_1DARRAY p,
   INT_1DARRAY cmember
   )
  {    

    amd_flag = BASKER_TRUE;

    //Debug,
    #ifdef BASKER_DEBUG_ORDER_AMD
    printf("cmember: \n");
    for(Int i = 0; i < M.ncol; ++i)
      {
	printf("(%d, %d), ", i, cmember(i));
      }
    printf("\n"); 
    #endif



    INT_1DARRAY temp_p;
    BASKER_ASSERT(M.ncol > 0, "AMD perm not long enough");
    MALLOC_INT_1DARRAY(temp_p, M.ncol+1);
    init_value(temp_p, M.ncol+1, (Int) 0);
    
    my_amesos_csymamd(M.ncol, &(M.col_ptr(0)), &(M.row_idx(0)),
		     &(temp_p(0)), &(cmember(0)));


    for(Int i = 0; i < M.ncol; ++i)
      {
	p(temp_p(i)) = i;
      }


  }//end csymamd()

 
  //======================COLAMD=======================

  template <class Int>
  BASKER_FINLINE
  int amesos_colamd
  (
   Int n_row, 
   Int n_col,
   Int Alen,
   Int *A,
   Int *p,
   double *knobs, 
   Int *stats
   )
  {
    return -1;
  }//end amesos_colamd()
  
 
  template < >
  BASKER_FINLINE
  int amesos_colamd<>
  (
   int n_row,
   int n_col, 
   int Alen,
   int *A,
   int *p,
   double *knobs,
   int *stats
   )
  {
    amesos_colamd(n_row,n_col,Alen,A,p,knobs,stats);
    return 0;
  }//end amesos_colamd<int>
  

  //template<class Entry, class Exe_Space>
  template <>
  BASKER_FINLINE
  int amesos_colamd<>
  (
   long n_row,
   long n_col,
   long Alen,
   long *A,
   long *p,
   double *knobs,
   long *stats
   )
  {
    amesos_colamd_l(n_row, n_col, Alen, A, p, knobs, stats);
    return 0;
  }
  


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::blk_amd(BASKER_MATRIX &M, INT_1DARRAY p)
  {
    
    //p == length(M)
    //Scan over all blks
    //Note, that this needs to be made parallel in the 
    //future (Future Josh will be ok with this, right?)

    //This is a horrible way to do this!!!!!
    //KLU does this very nice, but they also make all the little blks
    INT_1DARRAY temp_col;
    MALLOC_INT_1DARRAY(temp_col, M.ncol+1);
    INT_1DARRAY temp_row;
    MALLOC_INT_1DARRAY(temp_row, M.nnz);


    for(Int b = btf_tabs_offset; b < btf_nblks; b++)
      {
	Int blk_size = btf_tabs(b+1) - btf_tabs(b);
	if(blk_size < 3)
	  {
	    continue;
	  }
	
	INT_1DARRAY tempp;
	MALLOC_INT_1DARRAY(tempp, blk_size+1);
	
	
	//Fill in temp matrix
	Int nnz = 0;
	Int column = 1;
	temp_col(0) = 0;
	for(Int k = btf_tabs(b); k < btf_tabs(b+1); k++)
	  {
	    for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); i++)
	      {
		temp_row(nnz) = M.row_idx(i);
		nnz++;
	      }// end over all row_idx
	    temp_col(column) = nnz;
	    column++;
	  }//end over all columns k
	
       

      }//over all blk_tabs


  }//edn blk_amd()
      


}//end namespace BaskerNS

#endif
