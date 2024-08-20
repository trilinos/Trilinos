// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
//#define BASKER_TIMER

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
          &(temp_row(0)),&(tempp(0)), Options.verbose);

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
  int Basker<Int,Entry,Exe_Space>::btf_blk_mwm_amd
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

      return BASKER_SUCCESS;
    }

    #ifdef BASKER_TIMER
    Kokkos::Timer timer_order;
    double mwm_time = 0.0;
    double amd_time = 0.0;
    #endif

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

    const int blk_size_threshold = 1; // was 3
    //const int blk_size_threshold = 3; // was 3
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
      if (Options.verbose) {
        std::cout << std::endl << " + blk(" << b << "): size=" << blk_size
                  << ", rows=" << btf_tabs(b) << ":" << btf_tabs(b+1)-1 << std::endl;
      }

      // -----------------------------------------------
      // compute MWM of b-th diagonal block
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
      Int frow = btf_tabs(b);
      Int nnz = 0;
      Int column = 1;
      temp_col(0) = 0;
      for(Int k = btf_tabs(b); k < btf_tabs(b+1); k++)
      {
        for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); i++)
        {
          if(M.row_idx(i) < frow)
          { continue; }

          temp_val(nnz) = M.val(i);
          temp_row(nnz) = M.row_idx(i) - frow;
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
      #ifdef BASKER_TIMER
      timer_order.reset();
      #endif
      if (Options.blk_matching == 0) {
        // no mwm
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
          std::cout << " ** BLK_MWM_AMD::ShyLUBasker MWM (blk=" << b << ", " << btf_tabs(b) << ":" << btf_tabs(b+1)-1
                    << ", nnz = " << nnz << ") ** " << std::endl;
          //flag = false;
        }
        Int num_match = 0;
        mwm_order::mwm(blk_size, nnz,
                       &(temp_col(0)), &(temp_row(0)), &(temp_val(0)),
                       &(tempp(0)), num_match);
        if(num_match < blk_size) {
          if(Options.verbose == BASKER_TRUE) {
            std::cout << " ++ Num of matches returned " << num_match
                      << " is less than blk_size(" << b << ") = " << blk_size 
                      << std::endl;
          }
          FREE_INT_1DARRAY(temp_col);
          FREE_INT_1DARRAY(temp_row);
          FREE_ENTRY_1DARRAY(temp_val);
          return BASKER_ERROR;
        }
        for(Int ii = 0; ii < blk_size; ii++) {
          scale_row_array(btf_tabs(b)+ii) = one;
          scale_col_array(btf_tabs(b)+ii) = one;
        }
      }
      #if 0
      printf( " > debug: set blk_mwm to identity\n" );
      for(Int ii = 0; ii < blk_size; ii++) {
        scale_row_array(btf_tabs(b)+ii) = one;
        scale_col_array(btf_tabs(b)+ii) = one;
      }
      for(Int ii = 0; ii < blk_size; ii++) tempp(ii) = ii;
      #endif
      #ifdef BASKER_TIMER
      mwm_time += timer_order.seconds();
      #endif

      // apply MWM to rows
      permute_row(nnz, &(temp_row(0)), &(tempp(0)));

      //Add to the bigger perm vector
      for(Int ii = 0; ii < blk_size; ii++)
      {
        //printf( " mwm(%d+%d) = %d + %d\n",btf_tabs(b),ii,tempp(ii),btf_tabs(b) );
        p_mwm(ii+btf_tabs(b)) = tempp(ii)+btf_tabs(b);
      }
      /*{
        //for(Int k = 0; k < blk_size; k++)
        //{
        //  printf( " %d\n",tempp(k) );
        //}
        std::cout << "m2=[" << std::endl;
        for(Int k = 0; k < blk_size; k++)
        {
          for(Int i = temp_col(k); i < temp_col(k+1); i++)
            printf("%d %d %e\n", temp_row(i), k, temp_val(i));
        }
        std::cout << "];" << std::endl;
      }*/

      // -----------------------------------------------
      // compute AMD of b-th diagonal block
      #ifdef BASKER_TIMER
      timer_order.reset();
      #endif
      // use nnz(A+A') as default estimate for nnz(L)...
      double l_nnz = nnz;
      double lu_work = nnz;
      if (!Options.amd_dom || Options.static_delayed_pivot != 0) {
        // skip AMD ordering
        if (flag) {
          std::cout << " >> + Basker::NO_AMD_ORDER" << std::endl;
        }
        for(Int ii = 0; ii < blk_size; ii++) {
          tempp(ii) = ii;
        }
      } else {
        #ifdef BASKER_SORT_MATRIX_FOR_AMD
        // sort for calling AMD
        if (flag) {
          std::cout << " >> + Basker::SORT_MATRIX" << std::endl;
        }
        sort_matrix(nnz, blk_size, &(temp_col(0)), &(temp_row(0)), &(temp_val(0)));
        #endif

        #if 0
        if (1) {
          // ND on diagonal block
          // form Basker Matrix storing the diagonal block
          BASKER_MATRIX B;
          B.set_shape(0, blk_size, 0, blk_size);
          B.col_ptr = temp_col;
          B.row_idx = temp_row;
          B.nnz = nnz;

          // compute A'+A
          BASKER_MATRIX BBT;
          AplusAT(B,BBT);

          // remove diagonal entries
          Int nnzi = 0;
          Int nnz2 = 0;
          for(Int i = 0; i < blk_size; i++) {
            for(Int k = nnzi; k <BBT.col_ptr(i+1); k++) {
              if(BBT.row_idx(k) != i) {
                BBT.row_idx(nnz2) = BBT.row_idx(k);
                nnz2++;
              }
            }
            nnzi = BBT.col_ptr(i+1);
            BBT.col_ptr(i+1) = nnz2;
          }

          Kokkos::Timer timer_metis;
          #if 1
          // allocate inv-permutation vectors
          INT_1DARRAY metis_iperm;
          MALLOC_INT_1DARRAY(metis_iperm, blk_size+1);

          // call METIS
          if (Options.verbose) {
            printf( " calling METIS_NodeND(n = %d, nnz = %d) \n",blk_size,nnz );
          }
          METIS_NodeND(&blk_size, &(BBT.col_ptr(0)), &(BBT.row_idx(0)), NULL, NULL, &(tempp(0)), &(metis_iperm(0)));

          FREE_INT_1DARRAY(metis_iperm); 
          #else
          SCOTCH_Strat stradat;
          SCOTCH_stratInit(&stradat);

          SCOTCH_Graph grafdat;
          Int *vwgts = NULL;
          if( SCOTCH_graphBuild(&grafdat, 0, blk_size, &(BBT.col_ptr(0)), &(BBT.col_ptr(1)), vwgts, NULL,
                                nnz2, &(BBT.row_idx(0)), NULL) !=0 ) {
            printf( "\n ** SCOTCH_graphBuild failed **\n\n" );
          }

          Int *permtab = (Int *)malloc((blk_size)  *sizeof(Int));
          Int *peritab = (Int *)malloc((blk_size)  *sizeof(Int));
          #if 1
          SCOTCH_Ordering ordedat;
          if (SCOTCH_graphOrderInit (&grafdat, &ordedat, permtab, peritab,
                                     NULL, NULL, NULL) != 0) {
            printf( "\n ** SCOTCH_graphOrderInit failed **\n\n" );
          }

          if (Options.verbose) {
            printf( " calling Scotch_graphOrderCompute(n = %d, nnz = %d) \n",blk_size,nnz );
          }
          if (SCOTCH_graphOrderCompute (&grafdat, &ordedat, &stradat) != 0) {
            printf( "\n ** SCOTCH_graphOrderCompute failed **\n\n" );
          }
          SCOTCH_graphOrderExit (&grafdat, &ordedat);
          #else
          Int num_levels = 0;
          double balrat = 0.2;
          Int flagval = SCOTCH_STRATLEAFSIMPLE | SCOTCH_STRATSEPASIMPLE;
          if (SCOTCH_stratGraphOrderBuild(&stradat, flagval, num_levels, balrat) != 0) {
            printf( "\n ** SCOTCH_stratGraphOrderBuild failed **\n\n" );
          }

          Int cblk = 0;
          Int *rangtab = (Int *)malloc((blk_size+1)*sizeof(Int));
          Int *treetab = (Int *)malloc((blk_size)  *sizeof(Int));
          if (Options.verbose) {
            printf( " calling Scotch_graphOrder(n = %d, nnz = %d) \n",blk_size,nnz );
          }
          if (SCOTCH_graphOrder(&grafdat &stradat, permtab, peritab,
                                &cblk, rangtab, treetab) != 0) {
            printf( "\n ** SCOTCH_graphOrder failed **\n\n" );
          }
          #endif

          SCOTCH_graphExit (&grafdat);
          SCOTCH_stratExit (&stradat);
          for(Int i = 0; i < blk_size; i++) tempp(i) = peritab[i];
          #endif
          //for(Int i = 0; i < blk_size; i++) tempp(i) = i;
          if (Options.verbose) {
            double time_metis = timer_metis.seconds();
            printf( " calling NodeND took %.2e seconds\n",time_metis );
          }

          /*printf( " M = [\n" );
          for(Int i = 0; i < blk_size; i++) {
            for(Int k = BBT.col_ptr(i); k < BBT.col_ptr(i+1); k++) printf( "%d %d\n",i,BBT.row_idx(k) );
          }
          printf( "];\n" );
          printf( " P = [\n" );
          for(Int i = 0; i < blk_size; i++) printf("%d\n",tempp(i) );
          printf( "];\n" );*/

          FREE(BBT);
        } else
        #endif
        {
          // AMD on diagonal block
          if (flag) {
            std::cout << " >> + Basker::AMD_ORDER + <<" << std::endl;
          }
          int amd_info =
          BaskerSSWrapper<Int>::amd_order(blk_size, &(temp_col(0)), &(temp_row(0)), 
                                          &(tempp(0)), l_nnz, lu_work, Options.verbose);
          if (amd_info == TRILINOS_AMD_OUT_OF_MEMORY || amd_info == TRILINOS_AMD_INVALID) {
            if(Options.verbose == BASKER_TRUE) {
              std::cout << " ++ amd_order returned " << amd_info
                        << " for blk(" << b << ") with lu_work = "
                        << lu_work << std::endl;
            }
            FREE_INT_1DARRAY(temp_col);
            FREE_INT_1DARRAY(temp_row);
            FREE_ENTRY_1DARRAY(temp_val);
            return BASKER_ERROR;
          } else if (Options.verbose == BASKER_TRUE) {
            std::cout << " ++ amd_order returned " << amd_info
                      << " for blk(" << b << ") with lu_work = " << lu_work
                      << " and l_nnz = " << l_nnz << " ++" << std::endl;
          }
        }
      }
      #ifdef BASKER_TIMER
      amd_time += timer_order.seconds();
      #endif
      #if 0
      printf( " >> debug: set blk_amd to identity <<\n" );
      for(Int ii = 0; ii < blk_size; ii++)
      {
        tempp(ii) = ii;
      }
      #endif
      if (b < (Int)btf_nnz.extent(0)) {
        btf_nnz(b) = l_nnz;
      }
      if (b < (Int)btf_work.extent(0)) {
        if (0 <= (Int)lu_work) {
          btf_work(b) = lu_work;
        } else {
          btf_work(b) = Teuchos::OrdinalTraits<Int>::max();
        }
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
        //printf( " amd(%d + %d = %d) = %d + %d = %d\n",tempp(ii),btf_tabs(b),tempp(ii)+btf_tabs(b), ii,btf_tabs(b),ii+btf_tabs(b) );
        p_amd(tempp(ii)+btf_tabs(b)) = ii+btf_tabs(b);
      }
      /*std::cout << " tempp_" << b << " = [ " << std::endl;
      for(Int ii = 0; ii < blk_size; ii++)
      {
        std::cout << tempp(ii) << std::endl;
      }
      std::cout << " ]; " << std::endl;
      for(Int ii = 0; ii < blk_size; ii++)
      {
        printf( " p_amd(%d + %d = %d) = %d + %d = %d\n",tempp(ii),btf_tabs(b),tempp(ii)+btf_tabs(b), ii,btf_tabs(b),ii+btf_tabs(b) );
      }*/

      FREE_INT_1DARRAY(tempp);

    }//over all blk_tabs
    #ifdef BASKER_TIMER
    std::cout << " ++ Basker order : MWM time: " << mwm_time << std::endl;
    std::cout << " ++ Basker order : AMD time: " << amd_time << std::endl;
    #endif

    #ifdef BASKER_DEBUG_AMD_ORDER
    printf("blk amd final order\n");
    for(Int ii = 0; ii < M.ncol; ii++)
    {
      printf(" > amdP[%d] = %d\n", ii,p_amd(ii));
    }
    printf("\n");
    #endif

    FREE_INT_1DARRAY(temp_col);
    FREE_INT_1DARRAY(temp_row);
    FREE_ENTRY_1DARRAY(temp_val);

    return BASKER_SUCCESS;
  }//end blk_amd()

  // default: compute blk_mwm & blk_amd of all the blocks
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::btf_blk_mwm_amd
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
    return btf_blk_mwm_amd(b_start, b_num, M, p_mwm, p_amd, btf_nnz, btf_work);
  }

}//end namespace BaskerNS

#endif
