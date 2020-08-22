#ifndef SHYLUBASKER_ORDER_AMD_HPP
#define SHYLUBASKER_ORDER_AMD_HPP

//AMD is amazing on circuit problems in a way
//that ND can't

//This can be done to user the smaller domains

#include "trilinos_amd.h"
#include "trilinos_colamd.h"
#include "trilinos_ccolamd.h"
//#if defined(HAVE_AMESOS2_SUPERLUDIST) && !defined(BASKER_MC64)
//  #define BASKER_SUPERLUDIS_MC64
//#endif

namespace BaskerNS
{
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
    double knobs[TRILINOS_CCOLAMD_KNOBS];
    int    stats[TRILINOS_CCOLAMD_STATS];

    //use default knob settings
    trilinos_ccolamd_set_defaults(knobs);
    knobs[0] = 10;
    knobs[1] = 0;
    knobs[2] = 2;

    trilinos_csymamd(n, Ai, Ap,  p, knobs, stats, 
        &(calloc), &(free), 
        cmember, 0);

    //trilinos_csymamd_report(stats);

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
    double knobs[TRILINOS_CCOLAMD_KNOBS];
    long    stats[TRILINOS_CCOLAMD_STATS];

    //use default knob settings
    trilinos_ccolamd_l_set_defaults(knobs);
    knobs[0] = 10;
    knobs[1] = 0;
    knobs[2] = 2;

    trilinos_csymamd_l(n, Ai, Ap,  p, knobs, stats, 
		     &(calloc), &(free), 
		     cmember, 0);

    //trilinos_csymamd_l_report(stats);
    
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

    double amd_info[TRILINOS_AMD_INFO];
    trilinos_amd(M.ncol, &(M.col_ptr(0)), 
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

    //If doing  iluk, we will not want this.
    //See amd blk notes
    if(Options.incomplete == BASKER_TRUE)
    {
      for(Int i = 0; i < M.ncol; i++)
      {
        p(i) = i;
      }
      return;
    }

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
  int trilinos_colamd
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
  }//end trilinos_colamd()
  
 
  template < >
  BASKER_FINLINE
  int trilinos_colamd<>
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
    trilinos_colamd(n_row,n_col,Alen,A,p,knobs,stats);
    return 0;
  }//end trilinos_colamd<int>
  

  //template<class Entry, class Exe_Space>
  template <>
  BASKER_FINLINE
  int trilinos_colamd<>
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
    trilinos_colamd_l(n_row, n_col, Alen, A, p, knobs, stats);
    return 0;
  }
  


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::blk_amd( BASKER_MATRIX &M, INT_1DARRAY p )
  {
    //p == length(M)
    //Scan over all blks
    //Note, that this needs to be made parallel in the 
    //future (Future Josh will be ok with this, right?)

    //This is a horrible way to do this!!!!!
    //KLU does this very nice, but they also make all the little blks
    INT_1DARRAY temp_col;
    MALLOC_INT_1DARRAY(temp_col, M.ncol+1);
    INT_1DARRAY   temp_row;
    MALLOC_INT_1DARRAY  (temp_row, M.nnz);

    for(Int b = btf_tabs_offset; b < btf_nblks; b++)
    {
      Int blk_size = btf_tabs(b+1) - btf_tabs(b);
      if(blk_size < 3)
      {
        for(Int ii = 0; ii < blk_size; ++ii)
        {
          //printf("set %d \n", btf_tabs(b)+ii-M.scol);
          p(ii+btf_tabs(b)) = btf_tabs(b)+ii-M.scol;
        }
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
          if(M.row_idx(i) < btf_tabs(b))
          { continue; }

          temp_row(nnz) = M.row_idx(i) - btf_tabs(b);
          nnz++;
        }// end over all row_idx

        temp_col(column) = nnz;
        column++;
      }//end over all columns k

      #ifdef BASKER_DEBUG_ORDER_AMD
      printf("col_ptr: ");
      for(Int i = 0 ; i < blk_size+1; i++)
      {
        printf("%d, ", temp_col(i));
      }
      printf("\n");
      printf("row_idx: ");
      for(Int i = 0; i < nnz; i++)
      {
        printf("%d, ", temp_row(i));
      }
      printf("\n");
      #endif

      BaskerSSWrapper<Int>::amd_order(blk_size, &(temp_col(0)), 
          &(temp_row(0)),&(tempp(0)));

      #ifdef BASKER_DEBUG_ORDER_AMD
      printf("blk: %d order: \n", b);
      for(Int ii = 0; ii < blk_size; ii++)
      {
        //printf("%d, ", tempp(ii));
        std::cout << ii << " " << tempp(ii) << std::endl;
      }
      #endif

      //Add to the bigger perm vector
      for(Int ii = 0; ii < blk_size; ii++)
      {
        p(tempp(ii)+btf_tabs(b)) = ii+btf_tabs(b);
      }

      FREE_INT_1DARRAY(tempp);

    }//over all blk_tabs

    #ifdef BASKER_DEBUG_AMD_ORDER
    printf("blk amd final order\n");
    for(Int ii = 0; ii < M.ncol; ii++)
    {
      printf("%d, ", p(ii));
    }
    printf("\n");
    #endif

    FREE_INT_1DARRAY(temp_col);
    FREE_INT_1DARRAY(temp_row);

  }//end blk_amd()
      

  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::btf_blk_mwm_amd
  (
   Int b_start,
   Int b_num,
   BASKER_MATRIX &M, 
   INT_1DARRAY p_mwm, 
   INT_1DARRAY p_amd, 
   INT_1DARRAY btf_nnz, 
   INT_1DARRAY btf_work
  )
  {
    // printf("=============BTF_BLK_AMD_CALLED========\n");
    Int b_end = b_start+b_num;
    if(Options.incomplete == BASKER_TRUE)
    {
      //We note that AMD on incomplete ILUK
      //Seems realy bad and leads to a zero on the diag
      //Therefore, we simply return the natural ordering
      for(Int i = 0 ; i < M.ncol; i++)
      {
        p_amd(i) = i;
      }

      //We will makeup work to be 1, 
      //Since BTF is not supported in our iluk
      for(Int b = b_start; b < b_end; b++)
      {
        btf_nnz(b) = 1;
        btf_work(b) = 1;
      }

      return;
    }
 
    //p == length(M)
    //Scan over all blks
    //Note, that this needs to be made parallel in the 
    //future (Future Josh will be ok with this, right?)

    //This is a horrible way to do this!!!!!
    //KLU does this very nice, but they also make all the little blks
    INT_1DARRAY temp_col;
    MALLOC_INT_1DARRAY(temp_col, M.ncol+1);
    INT_1DARRAY temp_row;
    ENTRY_1DARRAY temp_val;
    MALLOC_INT_1DARRAY(temp_row, M.nnz);
    MALLOC_ENTRY_1DARRAY(temp_val, M.nnz);
    //printf("Done with btf_blk_mwm_amd malloc \n");
    //printf("blks: %d \n" , btf_nblks);

    //const int blk_size_threshold = 1; // was 3
    const int blk_size_threshold = 3; // was 3
    bool flag = Options.verbose;
    Entry one = (Entry)1.0;
    //MALLOC_ENTRY_1DARRAY (scale_row_array, A.nrow);
    //MALLOC_ENTRY_1DARRAY (scale_col_array, A.nrow);
    if (flag) {
      std::cout << std::endl << " block mwm + amd (" << b_start << " : " << b_end-1 << ")" << std::endl;
    }
    for(Int b = b_start; b < b_end; b++)
    {
      Int blk_size = btf_tabs(b+1) - btf_tabs(b);

      if(blk_size < blk_size_threshold)
      {
        if (flag) {
          std::cout << " >> BLK_MWM_AMD::NO BLK MWM (blk=" << b << ", " << btf_tabs(b) << ":" << btf_tabs(b+1)-1 << ", too small) << " << std::endl;
        }
        for(Int ii = 0; ii < blk_size; ++ii)
        {
          //printf("set amd(%d) = mwm(%d) = %d (scol=%d)\n", ii+btf_tabs(b),ii+btf_tabs(b), btf_tabs(b)+ii,M.scol);
          //p_amd(ii+btf_tabs(b)) = btf_tabs(b)+ii-M.scol;
          //p_mwm(ii+btf_tabs(b)) = btf_tabs(b)+ii-M.scol;
          p_amd(ii+btf_tabs(b)) = btf_tabs(b)+ii;
          p_mwm(ii+btf_tabs(b)) = btf_tabs(b)+ii;
          scale_row_array(ii+btf_tabs(b)) = one;
          scale_col_array(ii+btf_tabs(b)) = one;
        }

        if (b < (Int)btf_work.extent(0)) {
          btf_work(b) = blk_size*blk_size*blk_size;
        }
        if (b < (Int)btf_nnz.extent(0)) {
          btf_nnz(b)  = (.5*(blk_size*blk_size) + blk_size);
        }
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
          if(M.row_idx(i) < btf_tabs(b))
          { continue; }

          temp_val(nnz) = M.val(i);
          temp_row(nnz) = M.row_idx(i) - btf_tabs(b);
          nnz++;
        }// end over all row_idx

        temp_col(column) = nnz;
        column++;
      }//end over all columns k

      #ifdef BASKER_DEBUG_ORDER_AMD
      printf("col_ptr: ");
      for(Int i = 0 ; i < blk_size+1; i++)
      {
        printf("%d, ", temp_col(i));
      }
      printf("\n");
      printf("row_idx: ");
      for(Int i = 0; i < nnz; i++)
      {
        printf("%d, ", temp_row(i));
      }
      printf("\n");
      #endif

      /*{
        std::cout << " m = [ " << std::endl;
        for(Int k = 0; k < blk_size; k++)
        {
          for(Int i = temp_col(k); i < temp_col(k+1); i++)
            printf("%d %d %e\n", temp_row(i), k, temp_val(i));
        }
        std::cout << " ]; " << std::endl;
      }*/
      if (Options.blk_matching == 0) {
        // no mwm, for debugging
        if (flag) {
          std::cout << " ** BLK_MWM_AMD::NO BLK MWM (blk=" << b << ", " << btf_tabs(b) << ":" << btf_tabs(b+1)-1 << ") ** " << std::endl;
          //flag = false;
        }
        for(Int ii = 0; ii < blk_size; ii++) {
          scale_row_array(btf_tabs(b)+ii) = one;
          scale_col_array(btf_tabs(b)+ii) = one;
        }
        for(Int ii = 0; ii < blk_size; ii++) tempp(ii) = ii;
      }
#if defined(BASKER_MC64) ||  defined(BASKER_SUPERLUDIS_MC64)
      else if (Options.blk_matching == 2) {
        if (flag) {
          std::cout << " ** BLK_MWM_AMD::MC64 MWM (blk=" << b << ", " << btf_tabs(b) << ":" << btf_tabs(b+1)-1 << ") ** " << std::endl;
          //flag = false;
        }
        Int job = 5; //2 is the default for SuperLU_DIST
        mc64(blk_size, nnz, &(temp_col(0)), &(temp_row(0)), &(temp_val(0)),
             job, &(tempp(0)), &(scale_row_array(btf_tabs(b))), &(scale_col_array(btf_tabs(b))));
        /*{
          std::cout << "p_mwm=[" << std::endl;
          for(Int ii = 0; ii < blk_size; ii++) {
            std::cout << btf_tabs(b)+ii << " " << tempp(ii) << " " << scale_row_array(ii) << " " << scale_col_array(ii) << std::endl;
          }
          std::cout << "];" << std::endl;
        }*/
      }
#endif
      else { //if (Options.blk_matching == 1)
        if (flag) {
          std::cout << " ** BLK_MWM_AMD::ShyLUBasker MWM (blk=" << b << ", " << btf_tabs(b) << ":" << btf_tabs(b+1)-1 << ") ** " << std::endl;
          //flag = false;
        }
        Int num = 0;
        mwm_order::mwm(blk_size, nnz,
                       &(temp_col(0)), &(temp_row(0)), &(temp_val(0)),
                       &(tempp(0)), num);
        for(Int ii = 0; ii < blk_size; ii++) {
          scale_row_array(btf_tabs(b)+ii) = one;
          scale_col_array(btf_tabs(b)+ii) = one;
        }
      }
      /*{
        std::cout << "m1=[" << std::endl;
        for(Int k = 0; k < blk_size; k++)
        {
          for(Int i = temp_col(k); i < temp_col(k+1); i++)
            printf("%d %d %e\n", temp_row(i), k, temp_val(i));
        }
        std::cout << "];" << std::endl;
      }*/

      // apply MWM to rows
      permute_row(nnz, &(temp_row(0)), &(tempp(0)));
      // sort for calling AMD
      sort_matrix(nnz, blk_size, &(temp_col(0)), &(temp_row(0)), &(temp_val(0)));

      //Add to the bigger perm vector
      for(Int ii = 0; ii < blk_size; ii++)
      {
        //printf( " mwm(%d) = %d + %d\n",ii+btf_tabs(b),tempp(ii),btf_tabs(b) );
        //p_mwm(tempp(ii)+btf_tabs(b)) = ii+btf_tabs(b);
        p_mwm(ii+btf_tabs(b)) = tempp(ii)+btf_tabs(b);
      }
      /*{
        for(Int k = 0; k < blk_size; k++)
        {
          printf( " %d\n",tempp(k) );
        }
        std::cout << "m2=[" << std::endl;
        for(Int k = 0; k < blk_size; k++)
        {
          for(Int i = temp_col(k); i < temp_col(k+1); i++)
            printf("%d %d %e\n", temp_row(i), k, temp_val(i));
        }
        std::cout << "];" << std::endl;
      }*/

      double l_nnz = 0;
      double lu_work = 0;
      BaskerSSWrapper<Int>::amd_order(blk_size, &(temp_col(0)), 
          &(temp_row(0)),&(tempp(0)), 
          l_nnz, lu_work);

      if (b < (Int)btf_nnz.extent(0)) {
        btf_nnz(b) = l_nnz;
      }
      if (b < (Int)btf_work.extent(0)) {
        btf_work(b) = lu_work;
      }

      #ifdef BASKER_DEBUG_ORDER_AMD
      printf("blk: %d order: \n", b);
      for(Int ii = 0; ii < blk_size; ii++)
      {
        printf("%d, ", tempp(ii));
      }
      #endif

      //Add to the bigger perm vector
      for(Int ii = 0; ii < blk_size; ii++)
      {
        //printf( " amd(%d) = %d + %d\n",tempp(ii)+btf_tabs(b),ii,btf_tabs(b) );
        p_amd(tempp(ii)+btf_tabs(b)) = ii+btf_tabs(b);
      }
      if (Options.verbose) {
        printf( " blk(%d: size=%d, rows=%d:%d)\n",b, btf_tabs(b+1)-btf_tabs(b), btf_tabs(b),btf_tabs(b+1)-1 );
      }
      /*std::cout << " p_amd = [ " << std::endl;
      for(Int ii = 0; ii < blk_size; ii++)
      {
        std::cout << tempp(ii) << std::endl;
      }
      std::cout << " ]; " << std::endl;*/

      FREE_INT_1DARRAY(tempp);

    }//over all blk_tabs

    #ifdef BASKER_DEBUG_AMD_ORDER
    printf("blk amd final order\n");
    for(Int ii = 0; ii < M.ncol; ii++)
    {
      printf("%d, ", p(ii));
    }
    printf("\n");
    #endif

    FREE_INT_1DARRAY(temp_col);
    FREE_INT_1DARRAY(temp_row);
    FREE_ENTRY_1DARRAY(temp_val);
  }//end blk_amd()

  // default: compute blk_mwm & blk_amd of all the blocks
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::btf_blk_mwm_amd
  (
   BASKER_MATRIX &M, 
   INT_1DARRAY p_mwm, 
   INT_1DARRAY p_amd, 
   INT_1DARRAY btf_nnz, 
   INT_1DARRAY btf_work
  )
  {
    // MWM + AMD is applied ** all ** the diagonal block (both A & C)
    // the diagonal blocks are split into A & C later
    Int b_start = 0; 
    Int b_num = btf_nblks;
    btf_blk_mwm_amd(b_start, b_num, M, p_mwm, p_amd, btf_nnz, btf_work);
  }

}//end namespace BaskerNS

#endif
