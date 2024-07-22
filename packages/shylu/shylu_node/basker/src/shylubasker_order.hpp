// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_ORDER_HPP
#define SHYLUBASKER_ORDER_HPP

#include <Kokkos_Timer.hpp>

//Basker Includes
#include "shylubasker_types.hpp"
#include "shylubasker_util.hpp"

#include "shylubasker_order_match.hpp"
#include "shylubasker_order_scotch.hpp"
#include "shylubasker_order_btf.hpp"
#include "shylubasker_order_amd.hpp"

//#if defined(HAVE_AMESOS2_SUPERLUDIST) && !defined(BASKER_MC64)
//  #define BASKER_SUPERLUDIS_MC64
//#endif
//#undef BASKER_DEBUG_ORDER
//#define BASKER_TIMER


namespace BaskerNS
{
#ifdef BASKER_USE_QSORT
typedef struct basker_sort_pair_str {
  int perm;
  int rowid;
} basker_sort_pair;

static int basker_sort_matrix_col(const void *arg1, const void *arg2)
{
  const basker_sort_pair *val1 = (const basker_sort_pair *) arg1;
  const basker_sort_pair *val2 = (const basker_sort_pair *) arg2;
  return (val2->rowid < val1->rowid);
}
#endif

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::default_order()
  {
    match_ordering(1);
    partition(0);
    //camd?    
    return 0;
  }//end default order

  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::user_order
  (
   Int *perm,
   Int nblks,
   Int parts, 
   Int *row_tabs, Int *col_tabs,
   Int *tree_tabs
  )
  {
    //init tree structure
    init_tree(perm, nblks, parts, row_tabs, col_tabs, tree_tabs,0);

    //Find 2D structure
    matrix_to_views_2D(A);
    find_2D_convert(A);
    
    //Fill 2D structure
    #ifdef BASKER_KOKKOS
    kokkos_order_init_2D<Int,Entry,Exe_Space> iO(this);
    Kokkos::parallel_for(TeamPolicy(num_threads,1), iO);
    Kokkos::fence();
    #else
    //Comeback
    #endif
   
  }//end user_order


  //Need to completely reorganize

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::btf_order()
  {
    //1. Matching ordering on whole matrix
    //currently finds matching and permutes
    //found bottle-neck to work best with circuit problems
    sort_matrix(A);
    //printMTX("A_nonmatch.mtx", A);
    match_ordering(1);
    //printf("DEBUG1: done match\n");
    //for debuging
    sort_matrix(A);
    //printMTX("A_match.mtx", A);
   
    //2. BTF ordering on whole matrix
    // Gets estimate of work on all blocks
    //currently finds btf-hybrid and permutes
    //A -> [BTF_A, BTF_C; 0 , BTF B]

    MALLOC_INT_1DARRAY(btf_schedule, num_threads+1);
    init_value(btf_schedule, num_threads+1, 0);
    find_btf(A); 
   
    printf("Basker outer num_threads:%d, btf_tabs_offset:%d \n", num_threads, btf_tabs_offset);
    if(btf_tabs_offset != 0)
    {
      //  printf("A/B block stuff called\n");
      //3. ND on BTF_A
      //currently finds ND and permute BTF_A
      //Would like to change so finds permuation, 
      //and move into 2D-Structure
      //printMTX("A_BTF_FROM_A.mtx", BTF_A);
      sort_matrix(BTF_A);
      MALLOC_INT_1DARRAY(vals_order_scotch_array, BTF_A.nnz);
      scotch_partition(BTF_A);

      //need to do a row perm on BTF_B too
      if(btf_nblks > 1)
      {
        permute_row(BTF_B, part_tree.permtab);
      }
      //needed because  moving into 2D-Structure,
      //assumes sorted columns
      sort_matrix(BTF_A);
      if(btf_nblks > 1)
      {
        sort_matrix(BTF_B);
        sort_matrix(BTF_C);
      }
      //For debug
      //printMTX("A_BTF_PART_AFTER.mtx", BTF_A);

      //4. Init tree structure
      //This reduces the ND ordering into that fits,
      //thread counts
      init_tree_thread();

      //5. Permute BTF_A
      //Constrained symamd on A
      INT_1DARRAY cmember;
      MALLOC_INT_1DARRAY(cmember, BTF_A.ncol+1);
      init_value(cmember,BTF_A.ncol+1,(Int) 0);
      for(Int i = 0; i < tree.nblks; ++i)
      {
        for(Int j = tree.col_tabs(i); j < tree.col_tabs(i+1); ++j)
        {
          cmember(j) = i;
        }
      }
      MALLOC_INT_1DARRAY(order_csym_array, BTF_A.ncol);
      init_value(order_csym_array, BTF_A.ncol, (Int)0);
      csymamd_order(BTF_A, order_csym_array, cmember);

      permute_col(BTF_A, order_csym_array);
      sort_matrix(BTF_A);
      permute_row(BTF_A, order_csym_array);
      sort_matrix(BTF_A);
      //printMTX("A_BTF_AMD.mtx", BTF_A);

      if(btf_nblks > 1)
      {
        permute_row(BTF_B, order_csym_array);
        sort_matrix(BTF_B);
        //printMTX("B_BTF_AMD.mtx", BTF_B);
        sort_matrix(BTF_C);
        //printMTX("C_BTF_AMD.mtx", BTF_C);
      }

      //6. Move to 2D Structure
      //finds the shapes for both view and submatrices,
      //need to be changed over to just submatrices
      matrix_to_views_2D(BTF_A);

      //finds the starting point of A for submatrices
      find_2D_convert(BTF_A);

      //now we can fill submatrices 
      //(strictly-lower blocks in ALM, and upper blocks in AVM)
#ifdef BASKER_KOKKOS
      kokkos_order_init_2D<Int,Entry,Exe_Space> iO(this);
      Kokkos::parallel_for(TeamPolicy(num_threads,1), iO);
      Kokkos::fence();
#else
      //Comeback
#endif
      //printMTX("BTF_A.mtx", BTF_A); 

    }//if btf_tab_offset == 0
    
    if(btf_nblks > 1)
    {
      sort_matrix(BTF_C);
      //printMTX("C_TEST.mtx", BTF_C);
      //Permute C

      MALLOC_INT_1DARRAY(order_c_csym_array, BTF_C.ncol+1);
      init_value(order_c_csym_array, BTF_C.ncol+1,(Int) 0);

      printf("BEFORE \n");

      //csymamd_order(BTF_C, order_c_csym_array, cmember);

      blk_amd(BTF_C, order_c_csym_array);

      printf("After perm\n");

      permute_col(BTF_C, order_c_csym_array);
      sort_matrix(BTF_C);
      permute_row(BTF_C, order_c_csym_array);
      sort_matrix(BTF_C);

      if(BTF_E.ncol > 0)
      {
        permute_col(BTF_E, order_c_csym_array, BTF_A.ncol);
      }
      if(BTF_B.ncol > 0)
      {
        permute_col(BTF_B, order_c_csym_array);
        sort_matrix(BTF_B);
        //printMTX("BTF_B.mtx", BTF_B);
      }
      //printMTX("BTF_C.mtx", BTF_C);
    }
    return 0;
  }//end btf_order


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::btf_order2()
  {
    #ifdef BASKER_TIMER
    double order_time = 0.0;
    Kokkos::Timer timer_order;
    timer_order.reset();
    #endif

    //================================================================
    //1. Matching ordering on whole matrix
    //currently finds matching and permutes
    //found bottle-neck to work best with circuit problems
   
    //new for sfactor_copy2 replacement
    //this will store the composition of permutations and sorts due to btf with respect to the full val array
    btfd_nnz = 0;
    btfe_nnz = 0;
    btfa_nnz = 0;
    btfb_nnz = 0;
    btfc_nnz = 0;

    // NDE: New for Amesos2 CRS changes
    // compose the transpose+sort permutations with the running 'total' composition
    // Note: This is wasted work if the 'special case' is not hit
    MALLOC_INT_1DARRAY(vals_perm_composition, A.nnz);
    for( Int i = 0; i < A.nnz; ++i ){
      vals_perm_composition(i) = i;
    }
    permute_inv(vals_perm_composition, vals_crs_transpose, A.nnz);
    //A.print_matrix("a.dat");

    //--------------------------------------------------
    // 0: Basker, 1: trilinos, 2: MC64
    if(Options.verbose == BASKER_TRUE)
    {
        std::cout << " ++ calling match_ordering( " << Options.btf_matching << " )" << std::endl;
    }
    if (match_ordering(Options.btf_matching) != BASKER_SUCCESS) {
        if(Options.verbose == BASKER_TRUE) {
            std::cout << " ++ match_ordering failed ++ " << std::endl;
        }
        return BASKER_ERROR;
    }
    /*printf( " matchP=[\n" );
    for (int i = 0; i < A.ncol; i++) printf( "%d\n", order_match_array(i));
    printf( "];\n" );
    A.print_matrix("c.dat");*/


    #ifdef BASKER_TIMER
    order_time = timer_order.seconds();
    std::cout << " ++ Basker order : matching time: " << order_time << std::endl;
    timer_order.reset();
    #endif

    MALLOC_INT_1DARRAY(btf_schedule, num_threads+1);
    init_value(btf_schedule, num_threads+1, (Int)0);

    //================================================================
    //2. BTF ordering on whole matrix
    //   NOTE: Its also calls btf_blk_mwm_amd, 
    //         even if blk_mwm is enabled (for allocating workspace since amd provides stats for block LU)
    //currently finds btf-hybrid and permutes
    //A -> [BTF_A, BTF_C; 0 , BTF B]
    // where A is one "large" diagonal block for threaded factorization, and
    //       B contains "small" diagonabl blocks for sequential factorization
    //printf("outer num_threads:%d \n", num_threads);
    if (find_btf2(A) != BASKER_SUCCESS) {
        if(Options.verbose == BASKER_TRUE) {
            std::cout << " ++ find_btf2 failed ++ " << std::endl;
        }
        return BASKER_ERROR;
    }
    /*printf( " btfP=[\n" );
    for (int i = 0; i < A.ncol; i++) printf( "%d\n", order_btf_array(i));
    printf( "];\n" );*/

    #ifdef BASKER_TIMER
    order_time = timer_order.seconds();
    std::cout << " ++ Basker order : BTF time     : " << order_time << std::endl;
    timer_order.reset();
    #endif

    //================================================================
    // allocate mapping into D & E blocks
    if (btf_top_tabs_offset > 0) {
      btfd_nnz = BTF_D.nnz;
      btfe_nnz = BTF_E.nnz;

      if ( BTF_D.nnz > 0 ) {
        MALLOC_INT_1DARRAY(vals_order_ndbtfd_array, BTF_D.nnz);
        MALLOC_INT_1DARRAY(inv_vals_order_ndbtfd_array, BTF_D.nnz);
        for (Int i = 0; i < BTF_D.nnz; ++i) {
          vals_order_ndbtfd_array(i) = i;
          inv_vals_order_ndbtfd_array(i) = i;
        }
      }

      if ( BTF_E.nnz > 0 ) {
        MALLOC_INT_1DARRAY(vals_order_ndbtfe_array, BTF_E.nnz);
        MALLOC_INT_1DARRAY(inv_vals_order_ndbtfe_array, BTF_E.nnz);
        for (Int i = 0; i < BTF_E.nnz; ++i) {
          vals_order_ndbtfe_array(i) = i;
          inv_vals_order_ndbtfe_array(i) = i;
        }
      }
    }

    //================================================================
    //3. ND ordering of the big A block
    if (btf_tabs_offset != 0) // BTF_A exists and is not a btf_nblks > 1
    {
      // allocation perm vectors
      btfa_nnz = BTF_A.nnz;
      btfb_nnz = BTF_B.nnz;
      MALLOC_INT_1DARRAY(part_tree.permtab, BTF_A.ncol);
      MALLOC_INT_1DARRAY(order_csym_array,  BTF_A.ncol);
      MALLOC_INT_1DARRAY(order_csym_inv,    BTF_A.ncol);
      MALLOC_INT_1DARRAY(order_nd_inv,      BTF_A.ncol);
      if ( BTF_A.nnz > 0 ) {
        MALLOC_INT_1DARRAY(    vals_order_ndbtfa_array, BTF_A.nnz); //track nd perms
        MALLOC_INT_1DARRAY(inv_vals_order_ndbtfa_array, BTF_A.nnz);

        MALLOC_INT_1DARRAY(vals_order_scotch_array, BTF_A.nnz);
        MALLOC_INT_1DARRAY(vals_order_csym_array,   BTF_A.nnz);
      }
      if ( BTF_B.nnz > 0 ) {
        MALLOC_INT_1DARRAY(    vals_order_ndbtfb_array, BTF_B.nnz); 
        MALLOC_INT_1DARRAY(inv_vals_order_ndbtfb_array, BTF_B.nnz);
      }
      if ( BTF_C.nnz > 0 ) {
        MALLOC_INT_1DARRAY(    vals_order_ndbtfc_array, BTF_C.nnz); 
        MALLOC_INT_1DARRAY(inv_vals_order_ndbtfc_array, BTF_C.nnz);
      }

      // init perm to identity
      for (Int i = 0; i < BTF_A.nnz; ++i) {
        vals_order_ndbtfa_array(i) = i;
        inv_vals_order_ndbtfa_array(i) = i;

        vals_order_scotch_array(i) = i;
        vals_order_csym_array(i) = i;
      }
      for (Int i = 0; i < BTF_B.nnz; ++i) {
        vals_order_ndbtfb_array(i) = i;
        inv_vals_order_ndbtfb_array(i) = i;
      }
      for (Int i = 0; i < BTF_C.nnz; ++i) {
        vals_order_ndbtfc_array(i) = i;
        inv_vals_order_ndbtfc_array(i) = i;
      }
      for (Int i = 0; i < BTF_A.ncol; ++i) {
        part_tree.permtab(i) = i;
      }
      for (Int i = 0; i < BTF_A.ncol; ++i) {
        order_csym_array(i) = i;
      }

      // compute ND (if no MWM is requested)
      if (Options.blk_matching == 0 || Options.static_delayed_pivot != 0) {
        BASKER_BOOL keep_zeros = true;
        BASKER_BOOL compute_nd = true;
        BASKER_BOOL apply_nd   = true; //!(Options.static_delayed_pivot);
        int info_scotch = apply_scotch_partition(keep_zeros, compute_nd, apply_nd);
        if (info_scotch != BASKER_SUCCESS) {
          return info_scotch;
        }
      }
    } //if btf_tabs_offset != 0

    #ifdef BASKER_TIMER
    order_time = timer_order.seconds();
    std::cout << " ++ Basker order : sort(A) time : " << order_time << std::endl;
    timer_order.reset();
    #endif

    //================================================================
    // allocate mapping into C block
    btfc_nnz = BTF_C.nnz;
    if ( btf_tabs_offset == 0 && BTF_C.nnz > 0 ) {
      // NDE: May need to add permutation for this case...
      //new for sfactor_copy2 replacement
      MALLOC_INT_1DARRAY(vals_order_ndbtfc_array, BTF_C.nnz); //track nd perms; BTF_A must be declared here, else it does not exist
      MALLOC_INT_1DARRAY(inv_vals_order_ndbtfc_array, BTF_C.nnz);
      for (Int i = 0; i < BTF_C.nnz; ++i) {
        vals_order_ndbtfc_array(i) = i;
        inv_vals_order_ndbtfc_array(i) = i;
      }
      // NDE: already sorted above; this is redundant, unless btf_tabs_offset = 0
      //sort_matrix_store_valperms(BTF_C, vals_order_ndbtfc_array);
      //permute_inv(inv_vals_order_ndbtfc_array, vals_order_ndbtfc_array, BTF_C.nnz);
    }

    if(Options.verbose_matrix_out == BASKER_TRUE)
    {
      printMTX("C_Symbolic.mtx", BTF_C);
    }
    #ifdef BASKER_TIMER
    order_time = timer_order.seconds();
    std::cout << " ++ Basker order : sort(C) time : " << order_time << std::endl;
    #endif
    
    return BASKER_SUCCESS;
  }//end btf_order2

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::order_incomplete()
  {

    Kokkos::Timer timer_one;

    if(Options.matching == BASKER_TRUE)
    {
      match_ordering(1);
    }
    else
    {
      match_flag = BASKER_FALSE;
    }

    if(Options.btf == BASKER_TRUE)
    {
      //2. BTF ordering on whole matrix
      //currently finds btf-hybrid and permutes
      //A -> [BTF_A, BTF_C; 0 , BTF B]
      //printf("outter num_threads:%d \n", num_threads);
      MALLOC_INT_1DARRAY(btf_schedule, num_threads+1);
      init_value(btf_schedule, num_threads+1, (Int)0);
      find_btf2(A);
    }
    else
    {
      btf_flag = BASKER_FALSE;
    }

    std::cout << "Order Timer 1: "
      << timer_one.seconds()
      << std::endl;

    timer_one.reset();

    //ND Order (We should always do this for use)
    if(Options.btf == BASKER_TRUE && btf_tabs_offset != 0)
    {
      MALLOC_INT_1DARRAY(vals_order_scotch_array, BTF_A.nnz);
      scotch_partition(BTF_A);
      if(btf_nblks > 1)
      {
        permute_row(BTF_B, part_tree.permtab);
      }
      init_tree_thread();
      if(Options.amd_domains == BASKER_TRUE)
      {
        INT_1DARRAY cmember;
        MALLOC_INT_1DARRAY(cmember, BTF_A.ncol+1);
        init_value(cmember,BTF_A.ncol+1,(Int) 0);
        for(Int i = 0; i < tree.nblks; ++i)
        {
          for(Int j = tree.col_tabs(i); 
              j < tree.col_tabs(i+1); ++j)
          {
            cmember(j) = i;
          }
        }
        MALLOC_INT_1DARRAY(order_csym_array, BTF_A.ncol);
        init_value(order_csym_array, BTF_A.ncol, (Int)0);
        csymamd_order(BTF_A, order_csym_array, cmember, Options.verbose);

        permute_col(BTF_A, order_csym_array);
        sort_matrix(BTF_A);
        permute_row(BTF_A, order_csym_array);
        sort_matrix(BTF_A);

        if(btf_nblks > 1)
        {
          permute_row(BTF_B, order_csym_array);
          sort_matrix(BTF_B);
        }
      }//If amd order domains
      matrix_to_views_2D(BTF_A);
      find_2D_convert(BTF_A);
      kokkos_order_init_2D<Int,Entry,Exe_Space> iO(this);
      Kokkos::parallel_for(TeamPolicy(num_threads,1), iO);
      Kokkos::fence();
    }
    else
    {
      MALLOC_INT_1DARRAY(vals_order_scotch_array, A.nnz);
      scotch_partition(A);
      init_tree_thread();
      //Add domain order options
      matrix_to_views_2D(A);
      find_2D_convert(A);
      kokkos_order_init_2D<Int,Entry,Exe_Space> iO(this);
      Kokkos::parallel_for(TeamPolicy(num_threads,1), iO);
      Kokkos::fence();
    }

    std::cout << "Order Timer 2: "
      << timer_one.seconds()
      << std::endl;

    return;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry,Exe_Space>::partition(int option)
  {
    //Option does nothing right now
    if(Options.btf == BASKER_FALSE)
    {
      MALLOC_INT_1DARRAY(vals_order_scotch_array, A.nnz);
      scotch_partition(A);
    }
    else
    {
      MALLOC_INT_1DARRAY(vals_order_scotch_array, BTF_A.nnz);
      scotch_partition(BTF_A);
    }
    return 0;
  }//end partition()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry,Exe_Space>::match_ordering(int option)
  {
    /* ---- Tests --------

       INT_1DARRAY mperm;
       MALLOC_INT_1DARRAY(mperm, A.nrow);
       mc64(2,mperm);

       INT_1DARRAY mperm2;
       MALLOC_INT_1DARRAY(mperm2, A.nrow);
       mwm(A,mperm2);

       return 0;
    */

    //You would think a match order would help!
    //It does not!

    if (option < 0) {
      // no matching
      match_flag = BASKER_FALSE;
      if(Options.verbose == BASKER_TRUE) {
        std::cout << " ++ calling NO matching ++ " << std::endl;
      }
      if (Options.matrix_scaling != 0) {
        if(Options.verbose == BASKER_TRUE) {
          std::cout << " ++ allocate matrix sscaling ++ " << std::endl;
        }
        MALLOC_ENTRY_1DARRAY (scale_row_array, A.nrow);
        MALLOC_ENTRY_1DARRAY (scale_col_array, A.nrow);
      }
    } else {
      match_flag = BASKER_TRUE;
      //A.print_matrix("A.dat");

      Entry one(1.0);
      int num_match = min(A.nrow, A.ncol);
      MALLOC_INT_1DARRAY(order_match_array, A.nrow);
      MALLOC_ENTRY_1DARRAY (scale_row_array, A.nrow);
      MALLOC_ENTRY_1DARRAY (scale_col_array, A.nrow);
      for(Int i = 0; i < A.nrow; i++)
      {
        order_match_array(i) = i;
        scale_row_array(i) = one;
        scale_col_array(i) = one;
      }
      if(Options.incomplete == BASKER_FALSE)
      {
        if (option == 0) {
          if(Options.verbose == BASKER_TRUE) {
            std::cout << " ++ no BTF matching ++ " << std::endl;
          }
        }
        else if (option == 1) {
          if(Options.verbose == BASKER_TRUE) {
            std::cout << " ++ calling ShyLUBasker::MWM (" << A.nrow << " x " << A.ncol << ") ++ " << std::endl;
          }
          num_match = mwm(A, order_match_array);
        } 
        #if defined(BASKER_MC64) || defined(BASKER_SUPERLUDIS_MC64)
        else if (option == 3) {
          Int job = 5; //2 is the default for SuperLU_DIST
          // call mc64
          if(Options.verbose == BASKER_TRUE) {
            std::cout << " ++ calling MC64 ++ " << std::endl;
          }
          mc64(job, order_match_array, scale_row_array, scale_col_array);
          #if 0
          for(Int j = 0; j < A.ncol; j++) {
            for(Int k = A.col_ptr[j]; k < A.col_ptr[j+1]; k++) {
              A.val[k] *= (scale_col_array(j) * scale_row_array(A.row_idx[k]));
            }
          }
          #endif
        }
        #endif
        else { // if (option == 2 or default)
          double maxwork = 0.0;
          double work;
          INT_1DARRAY WORK;
          MALLOC_INT_1DARRAY(WORK, 5 * (A.ncol));
          if(Options.verbose == BASKER_TRUE) {
            std::cout << " ++ calling TRILINOS_BTF_MAXTRANS (" << A.nrow << " x " << A.ncol << ") ++ " << std::endl;
          }
          if (std::is_same<Int, int>::value) {
            int *col_ptr = reinterpret_cast <int*> (&(A.col_ptr(0)));
            int *row_idx = reinterpret_cast <int*> (&(A.row_idx(0)));
            int *order   = reinterpret_cast <int*> (&(order_match_array(0)));
            int *iwork   = reinterpret_cast <int*> (&(WORK(0)));
            num_match = trilinos_btf_maxtrans ((int)A.nrow, (int)A.ncol, col_ptr, row_idx, maxwork, &work, order, iwork);
          } else {
            long int *col_ptr = reinterpret_cast <long int*> (&(A.col_ptr(0)));
            long int *row_idx = reinterpret_cast <long int*> (&(A.row_idx(0)));
            long int *order   = reinterpret_cast <long int*> (&(order_match_array(0)));
            long int *iwork   = reinterpret_cast <long int*> (&(WORK(0)));
            num_match = trilinos_btf_l_maxtrans ((long int)A.nrow, (long int)A.ncol, col_ptr, row_idx, maxwork, &work, order, iwork);
          }
          #if 0
          printf( " >> debug: set order_match to identity <<\n" );
          for(Int i = 0; i < A.nrow; i++)
          {
            order_match_array(i) = i;
          }
          #endif
          FREE_INT_1DARRAY(WORK);
        }
        // apply matching to the rows of A
        permute_row(A, order_match_array);
      }
      else
      {
        for(Int i = 0; i < A.nrow; i++)
        {
          order_match_array(i) = i;
        }
      }
      //printf( " match_array\n" );
      //for(Int j = 0; j < A.ncol; j++) printf( " > %d\n",order_match_array(j) );
      //A.print_matrix("B.dat");
      if(num_match < min(A.nrow, A.ncol)) {
        if(Options.verbose == BASKER_TRUE) {
          std::cout << " ++ Num of matches returned " << num_match
                    << " is less than nrow = " << A.nrow << " or ncol = " << A.ncol
                    << std::endl;
        }
        return BASKER_ERROR;
      }

      //We want to test what the match ordering does if
      //have explicit zeros
#ifdef BASKER_DEBUG_ORDER
      FILE *fp;
      fp = fopen("match_order.txt", "w");
      printf("Matching Perm \n");
      for(Int i = 0; i < A.nrow; i++)
      {
        printf("%d, \n", order_match_array(i));
        fprintf(fp, "%d \n", order_match_array(i));
      }
      printf("\n");
      fclose(fp);
#endif
    }

    //May have to call row_idx sort
    return BASKER_SUCCESS;
  }//end match_ordering()


  //--------------------------------------------------------------
  template <class Int>
  struct kokkos_amd_order
  {
    Int nleaves;
    Int nblks;
    INT_1DARRAY col_tabs;
    INT_1DARRAY col_ptr;
    INT_1DARRAY row_idx;

    INT_1DARRAY tempp;
    INT_1DARRAY temp_col;
    INT_1DARRAY temp_row;

    INT_1DARRAY order_csym_array;
    BASKER_BOOL verbose;

    kokkos_amd_order(Int _nleaves,
                     Int _nblks,
                     INT_1DARRAY _col_tabs,
                     INT_1DARRAY _col_ptr,
                     INT_1DARRAY _row_idx,
                     //
                     INT_1DARRAY _tempp,
                     INT_1DARRAY _temp_col,
                     INT_1DARRAY _temp_row,
                     //
                     INT_1DARRAY _order_csym_array,
                     BASKER_BOOL _verbose) :
    nleaves (_nleaves),
    nblks   (_nblks),
    col_tabs(_col_tabs),
    col_ptr (_col_ptr),
    row_idx (_row_idx),
    tempp   (_tempp),
    temp_col(_temp_col),
    temp_row(_temp_row),
    order_csym_array(_order_csym_array),
    verbose(_verbose)
    {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int id) const {
      for (Int b = id; b < nblks; b += nleaves) {
        Int frow = col_tabs(b);
        Int erow = col_tabs(b+1);
        Int fnnz = col_ptr(frow);
        Int blk_size = erow - frow;
        //std::cout << " amd_order_functor ( " << id << " ): block " << b << " size = " << blk_size << std::endl;

        Int nnz = 0;
        temp_col(frow+b) = 0;
        //printf( " id=%d: frow=%d, erow=%d, fnnz=%d (nblks=%d, nleaves=%d)\n",id,frow,erow,fnnz,nblks,nleaves );
        //printf( "C=[\n");
        for(Int k = frow; k < erow; k++) {
          for(Int i = col_ptr(k); i < col_ptr(k+1); i++) {
            if(row_idx(i) >= frow && row_idx(i) < erow) {
              temp_row(fnnz + nnz) = row_idx(i) - frow;
              //printf( " %d %d %d\n",id,temp_row(fnnz + nnz),k );
              nnz++;
            }
          }
          temp_col(b+k+1) = nnz;
        }
        //printf( "];\n");
        #ifdef BASKER_SORT_MATRIX_FOR_AMD
        sort_matrix(nnz, blk_size, &(temp_col(frow+b)), &(temp_row(fnnz)), &(temp_val(fnnz)));
        #endif

        double l_nnz = 0;
        double lu_work = 0;
        BaskerSSWrapper<Int>::amd_order(blk_size, &(temp_col(frow+b)), &(temp_row(fnnz)),
                                        &(tempp(frow)), l_nnz, lu_work, verbose);
        for(Int k = 0; k < blk_size; k++)
        {
          order_csym_array(frow+tempp(frow + k)) = frow+k;
        }
      }
    }
  };
  //--------------------------------------------------------------

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::apply_scotch_partition(BASKER_BOOL keep_zeros, BASKER_BOOL compute_nd, BASKER_BOOL apply_nd)
  {
    #if 1//def BASKER_TIMER
    Kokkos::Timer scotch_timer;
    #endif

    //--------------------------------------------------------------
    //3. ND on BTF_A (compute ND and apply/permute BTF_A)
    //currently finds ND and permute BTF_A
    //Would like to change so finds permuation, 
    //and move into 2D-Structure
    //printf(" in apply_scotch_partition\n" );
    //AAT.print_matrix("T.dat");
    BASKER_MATRIX AAT;
    if(Options.symmetric == BASKER_TRUE) {
      AAT = BTF_A;
    } else {
      // compute AAT here
      AplusAT(BTF_A, AAT, keep_zeros);
    }
    //AAT.print_matrix("AAT_.dat");
    #ifdef BASKER_TIMER
    double apat_time = scotch_timer.seconds();
    std::cout << " ++ Basker apply_scotch : ++ AplusAT  : nnz = " << BTF_A.nnz << " -> " << AAT.nnz
              << " ( symmetric = " << Options.symmetric << ", keep_zeros = " << keep_zeros << " ) "
              << apat_time << " seconds" << std::endl;
    scotch_timer.reset();
    #endif

    // tree is prepped; permutation then applied to BTF_A
    /*printf(" Before permute:\n" );
    BTF_A.print_matrix("A1.dat");*/
    int info_scotch = 0;
    if (compute_nd) {
      if (Options.symmetric == BASKER_TRUE) {
        info_scotch = scotch_partition(BTF_A, apply_nd);
      } else {
        info_scotch = scotch_partition(BTF_A, AAT, apply_nd);
      }
    } else if (apply_nd) {
      // permtab is applied inside scotch_partition (also stored in vals_order_scotch_array, which is not needed) if compute_nd
      permute_row(BTF_A, part_tree.permtab);
      permute_col(BTF_A, part_tree.permtab);
    }
    //printf(" After permute:\n" );
    //for(Int j = 0; j < AAT.ncol; j++) printf( " %d %d %d\n",j,part_tree.permtab(j),part_tree.ipermtab(j) );
    //BTF_A.print_matrix("A2.dat");

    if (info_scotch != BASKER_SUCCESS || !apply_nd) {
      if(Options.verbose == BASKER_TRUE) {
        std::cout << " > scotch_partition returned info = " << info_scotch << " with apply_nd = " << apply_nd << std::endl;
      }
      return info_scotch;
    }
    if (Options.symmetric != BASKER_TRUE) {
      // apply ND to AAT (scotch_partition applies to BTF_A, only)
      // TODO: just call permute_row & permute_col?
      INT_1DARRAY col_ptr;
      INT_1DARRAY row_idx;
      MALLOC_INT_1DARRAY(col_ptr, AAT.ncol+1);
      MALLOC_INT_1DARRAY(row_idx, AAT.nnz);

      Int nnz = 0;
      col_ptr[0] = 0;
      for (Int j = 0; j < AAT.ncol; j++) {
        Int col = part_tree.ipermtab[j];
        for (Int k = AAT.col_ptr[col]; k < AAT.col_ptr[col+1]; k++) {
          row_idx[nnz] = part_tree.permtab[AAT.row_idx[k]];
          nnz++;
        }
        col_ptr[j+1] = nnz;
      }
      for (Int j = 0; j <= AAT.ncol; j++) {
        AAT.col_ptr[j] = col_ptr[j];
      }
      for (Int k = 0; k < nnz; k++) {
        AAT.row_idx[k] = row_idx[k];
      }
    }
    #ifdef BASKER_TIMER
    double scotch_time = scotch_timer.seconds();
    std::cout << " ++ Basker apply_scotch : scotch partition time     : " << scotch_time << std::endl;
    scotch_timer.reset();
    #endif
    /*printf(" After ND\n");
    printf(" p = [\n" );
    for(Int j = 0; j < BTF_A.ncol; j++) {
      printf("%d\n",part_tree.permtab(j) );
    }
    printf("];\n");
    BTF_A.print_matrix("pA.dat");
    AAT.print_matrix("AAT.dat");*/

    //need to do a col perm on BTF_E and a row perm on BTF_B too
    if (BTF_E.ncol > 0) {
      //permute_col(BTF_E, part_tree.permtab);
      permute_col_store_valperms(BTF_E, part_tree.permtab, inv_vals_order_ndbtfe_array);
      if (Options.blk_matching == 0) {
        permute_inv(vals_order_ndbtfe_array, inv_vals_order_ndbtfe_array, BTF_E.nnz);
      }
    }
    if(btf_nblks > 1)
    {
      permute_row(BTF_B, part_tree.permtab);
    }
    //BTF_A.print_matrix("pT.dat");

    //--------------------------------------------------------------
    //4. Init tree structure
    //This reduces the ND ordering into that fits,
    //thread counts
    init_tree_thread();
    #ifdef BASKER_TIMER
    Kokkos::Timer scotch_amd_timer;
    double tree_time = scotch_timer.seconds();
    std::cout << " ++ Basker apply_scotch : init tree time            : " << tree_time << std::endl;
    scotch_timer.reset();
    double csymamd_time = 0.0;
    #endif


    //--------------------------------------------------------------
    //5. Constrained symamd on A
    //Init for Constrained symamd on A
    bool run_nd_on_leaves  = Options.run_nd_on_leaves;
    bool run_amd_on_leaves = Options.run_amd_on_leaves;
    //if (Options.amd_dom && Options.static_delayed_pivot == 0)
    if (Options.amd_dom && (!run_nd_on_leaves && !run_amd_on_leaves)) // still compute csymamd in symbolic for delayed pivot (to reduce cost of setup in numeric)
    {
      if (Options.symmetric != BASKER_TRUE) { // TODO: replace with parameter, e.g., use_csymamd
        // flag for permute_composition_for_solve
        amd_flag = BASKER_TRUE;

        Int nblks = tree.nblks;
        INT_1DARRAY tempp;
        INT_1DARRAY temp_col;
        INT_1DARRAY temp_row;
        MALLOC_INT_1DARRAY  (tempp,    AAT.ncol);
        MALLOC_INT_1DARRAY  (temp_col, AAT.ncol+nblks);
        MALLOC_INT_1DARRAY  (temp_row, AAT.nnz);

        #ifdef BASKER_TIMER
        scotch_amd_timer.reset();
        #endif
        #if 0
        Int nleaves = 1;
        printf( " * debug nleaves = 1 *\n" );
        #else
        Int nleaves = num_threads;
        #endif
        kokkos_amd_order<Int> amd_functor(nleaves, nblks, tree.col_tabs, AAT.col_ptr, AAT.row_idx,
                                          tempp, temp_col, temp_row, order_csym_array, Options.verbose);
        Kokkos::parallel_for("BLK_AMD on A", Kokkos::RangePolicy<Exe_Space>(0, nleaves), amd_functor);
        Kokkos::fence();
        if(Options.verbose == BASKER_TRUE) {
          std::cout << " ++ Basker AMD_functor on A ++ " << std::endl << std::endl;
        }
        #ifdef BASKER_TIMER
        csymamd_time += scotch_amd_timer.seconds();
        std::cout << " ++ Basker apply_scotch : constrained symm amd time : " << csymamd_time
                  << " using amd" << std::endl;
        scotch_timer.reset();
        #endif
      } else
      {
        INT_1DARRAY cmember;
        MALLOC_INT_1DARRAY(cmember, BTF_A.ncol+1);
        for(Int i = 0; i < tree.nblks; ++i) {
          for(Int j = tree.col_tabs(i); j < tree.col_tabs(i+1); ++j) {
            cmember(j) = i;
          }
        }
        //--------------------------------------------------------------
        // call Constrained symamd on A
        //printf("============CALLING CSYMAMD============\n");
        init_value(order_csym_array, BTF_A.ncol, (Int)0);

        #ifdef BASKER_TIMER
        scotch_amd_timer.reset();
        #endif
        csymamd_order(BTF_A, order_csym_array, cmember);
        #ifdef BASKER_TIMER
        csymamd_time = scotch_amd_timer.seconds();
        #endif

        if(Options.verbose == BASKER_TRUE) {
          std::cout << " ++ Basker CSYMAMD_functor on A ++ " << std::endl << std::endl;
        }
        #ifdef BASKER_TIMER
        double csymamd_tot_time = scotch_timer.seconds();
        std::cout << " ++ Basker apply_scotch : constrained symm amd time : " << csymamd_tot_time
                  << " (" << csymamd_time << ") using csymamd" << std::endl;
        scotch_timer.reset();
        #endif
      }
      #if 0 // reset to I for debug
      printf( " >> debug: set order_csym_array to identity <<\n" );
      for (Int i = 0; i < BTF_A.ncol; ++i) {
        order_csym_array(i) = i;
      }
      #endif
    } else {
      if(Options.verbose == BASKER_TRUE) {
        std::cout << " ++ Basker no AMD on A ++ " << std::endl << std::endl;
      }
    }
    //sort_matrix(BTF_A); // unnecessary?

    //new for sfactor_copy2 replacement
    permute_col_store_valperms(BTF_A, order_csym_array, vals_order_csym_array); //NDE: Track movement of vals (lin_ind of row,col) here
    if (Options.blk_matching == 0) {
      permute_inv(vals_order_ndbtfa_array, vals_order_csym_array, BTF_A.nnz);     //must permute the array holding the perms
    }
    //sort_matrix_store_valperms(BTF_A, vals_order_ndbtfa_array);
    permute_row(BTF_A, order_csym_array);

    //printf("pp=[");
    //for(Int j = 0; j < BTF_A.ncol; j++) printf("%d\n",order_csym_array(j));
    //printf("]\n");
    //BTF_A.print_matrix("ppT.dat");

    // ----------------------------------------------
    // sort matrix(BTF_A);
    #ifdef BASKER_TIMER
    Kokkos::Timer scotch_sort_timer;
    #endif
    //for (Int j = 0; j < BTF_A.nnz; j++) printf( " - vals_ndbtfa(%d) = %d\n",j,vals_order_ndbtfa_array(j) );
#if defined(SHYLU_BASKER_SORT_BLOCK_A)
    sort_matrix_store_valperms(BTF_A, vals_order_ndbtfa_array); //new for replacement
#else
    ndsort_matrix_store_valperms(BTF_A, vals_order_ndbtfa_array);
#endif
    #ifdef BASKER_TIMER
    double sortA_time = scotch_timer.seconds();
    double sortA_time1 = scotch_sort_timer.seconds();
    std::cout << " ++ Basker apply_scotch : sort(A) time              : " << sortA_time
              << " ( " << sortA_time1 << " )" << std::endl;
    scotch_timer.reset();
    #endif

    if (BTF_E.ncol > 0) {
      //permute_col(BTF_E, order_csym_array);
      permute_col_store_valperms(BTF_E, order_csym_array, inv_vals_order_ndbtfe_array);
      if (Options.blk_matching == 0) {
        permute_inv(vals_order_ndbtfe_array, inv_vals_order_ndbtfe_array, BTF_E.nnz);
      }
    }
    if(btf_nblks > 1) {
      permute_row(BTF_B, order_csym_array);
      //new for sfactor_copy2 replacement
      if ( BTF_B.nnz > 0 ) {
        #ifdef BASKER_TIMER
        double sortB_time = scotch_timer.seconds();
        std::cout << " ++ Basker apply_scotch : sort(B) time              : " << sortB_time << std::endl;
        scotch_timer.reset();
        #endif
      }
      if ( BTF_C.nnz > 0 ) {
        #ifdef BASKER_TIMER
        double sortC_time = scotch_timer.seconds();
        std::cout << " ++ Basker apply_scotch : sort(C) time              : " << sortC_time << std::endl;
        scotch_timer.reset();
        #endif
      }
    }

    //new for sfactor_copy2 replacement
    // Steps below necessary for later use in Factor(...) to follow permutation convention used up to this point
    if (BTF_D.nnz > 0) {
      permute_inv(inv_vals_order_ndbtfd_array, vals_order_ndbtfd_array, BTF_D.nnz);
    }
    if (BTF_E.nnz > 0) {
      for (int i = 0; i < BTF_E.nnz; i++) inv_vals_order_ndbtfe_array(i) = i;
      permute_inv(inv_vals_order_ndbtfe_array, vals_order_ndbtfe_array, BTF_E.nnz);
    }
    if (BTF_A.nnz > 0) {
      permute_inv(inv_vals_order_ndbtfa_array, vals_order_ndbtfa_array, BTF_A.nnz);
    }
    if (BTF_B.nnz > 0) { //this may not be the best way to test...
      permute_inv(inv_vals_order_ndbtfb_array, vals_order_ndbtfb_array, BTF_B.nnz);
    }
    if (BTF_C.nnz > 0) { //this may not be the best way to test...
      permute_inv(inv_vals_order_ndbtfc_array, vals_order_ndbtfc_array, BTF_C.nnz);
    }
    //printMTX("BTF_A.mtx", BTF_A);
    #ifdef BASKER_TIMER
    double perm_time = scotch_timer.seconds();
    std::cout << " ++ Basker apply_scotch : perm matrix time          : " << perm_time << std::endl;
    scotch_timer.reset();
    #endif

    //--------------------------------------------------------------
    //7. Move to 2D Structure
    //finds the shapes for both view and submatrices,
    //need to be changed over to just submatrices
    matrix_to_views_2D(BTF_A);
    //finds the starting point of A for submatrices
    find_2D_convert(BTF_A);
    //now we can fill submatrices
    #ifdef BASKER_KOKKOS
    kokkos_order_init_2D<Int,Entry,Exe_Space> iO(this);
    Kokkos::parallel_for(TeamPolicy(num_threads,1), iO);
    Kokkos::fence();
    #else
    //Comeback
    #endif
    #ifdef BASKER_TIMER
    double init_2d_time = scotch_timer.seconds();
    std::cout << " ++ Basker apply_scotch : init 2D time              : " << init_2d_time << std::endl;
    scotch_timer.reset();
    #endif

    #if 1
    {
      // initialize threading  
      basker_barrier.init(num_threads, 16, tree.nlvls );

      #ifdef BASKER_TIMER
      double barrier_init_time = scotch_timer.seconds();
      std::cout << " ++ Basker apply_scotch : barrier init time         : " << barrier_init_time << std::endl;
      #endif
    }
    #endif

    return 0;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::scotch_partition(BASKER_MATRIX &M, BASKER_BOOL apply_nd)
  {
    int info_scotch = 0;
    if(Options.symmetric == BASKER_TRUE) {
      info_scotch = scotch_partition(M, M, apply_nd);
    } else {
      BASKER_MATRIX MMT;
      AplusAT(M,MMT);
      info_scotch = scotch_partition(M, MMT, apply_nd);
      FREE(MMT);
    }
    return info_scotch;
  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::scotch_partition(BASKER_MATRIX &M, BASKER_MATRIX &MMT, BASKER_BOOL apply_nd)
  {
    nd_flag = BASKER_FALSE;

    int info_scotch = part_scotch(MMT, part_tree);
    if (info_scotch != BASKER_SUCCESS || !apply_nd) {
      if(Options.verbose == BASKER_TRUE) {
        std::cout << " > scotch_partition returned with info = " << info_scotch << " and apply_nd = " << apply_nd << std::endl;
      }
      return info_scotch;
    }

    nd_flag = BASKER_TRUE;
    //permute
    permute_row(M, part_tree.permtab);
    //for(Int i = 0; i < M.nrow; i++) printf("%d %d\n",i,part_tree.permtab(i));

    // new sfactor_copy2 replacement changes
    //permute_col(M, part_tree.permtab); //old, try the new below

    for (Int i = 0; i < M.nnz; ++i) {
      vals_order_scotch_array(i) = i; //init
    }
    permute_col_store_valperms(M, part_tree.permtab, vals_order_scotch_array); //NDE: Track movement of vals (lin_ind of row,col) here
    if (Options.blk_matching == 0 || Options.static_delayed_pivot != 0) {
      permute_inv(vals_order_ndbtfa_array, vals_order_scotch_array, M.nnz); //must permute the array holding the perms
    }

    //May need to sort row_idx
    return BASKER_SUCCESS; 
  }//end scotch_partition()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::permute_inv
  (
   INT_1DARRAY vec,
   INT_1DARRAY p,
   Int n
  )
  {
    if (n > 0) {
      INT_1DARRAY temp;
      MALLOC_INT_1DARRAY(temp,n);

      for(Int i = 0; i < n; ++i)
      {
        temp(p(i)) = vec(i);
      }
      for(Int i = 0; i < n; ++i)
      {
        vec(i) = temp(i);
      }

      FREE_INT_1DARRAY(temp);
    }

    return BASKER_SUCCESS;
  }//end permute_inv (int,int)
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_inv
  (
   ENTRY_1DARRAY vec,
   INT_1DARRAY   p, 
   Int n
  )
  {
    ENTRY_1DARRAY temp;
    MALLOC_ENTRY_1DARRAY(temp, n);
    init_value(temp, n, (Entry) 0.0);

    //Permute
    for(Int i = 0; i < n; i++)
    {
      temp(p(i)) = vec(i);
    }
    //Copy back
    for(Int i = 0; i < n; i++)
    {
      vec(i) = temp(i);
    }

    FREE_ENTRY_1DARRAY(temp);

    return 0;
  }


  //JDB: need to comeback
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_inv
  (
   ENTRY_1DARRAY vec,
   INT_1DARRAY   p, 
   Int n, //size(vec) //size(vec) > size(p)
   Int m, //size(p)
   Int start
  )
  {
    if(m > n)
    {
      printf("ERROR Permute inv \n");
      return BASKER_ERROR;
    }

    ENTRY_1DARRAY temp;
    MALLOC_ENTRY_1DARRAY(temp, n);
    init_value(temp, n, (Entry) 0.0);

    //Permute
    for(Int i = 0; i < n; i++)
    {
      //temp(i) = vec(p(i));
      temp(p(i)+start) = vec(i+start);
    }
    //Copy back
    for(Int i = 0; i < n; i++)
    {
      vec(i+start) = temp(i+start);
    }

    FREE_ENTRY_1DARRAY(temp);

    return 0;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute
  (
   ENTRY_1DARRAY vec,
   INT_1DARRAY   p, 
   Int n
  )
  {
    ENTRY_1DARRAY temp;
    MALLOC_ENTRY_1DARRAY(temp, n);
    init_value(temp, n, (Entry) 0.0);

    //Permute
    for(Int i = 0; i < n; i++)
    {
      temp(i) = vec(p(i));
    }
    //Copy back
    for(Int i = 0; i < n; i++)
    {
      vec(i) = temp(i);
    }

    FREE_ENTRY_1DARRAY(temp);

    return 0;
  }


  //JDB:: come back
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute
  (
   ENTRY_1DARRAY vec,
   INT_1DARRAY   p, 
   Int n,  //n = size(vec) //n > m 
   Int m,  //m = size(p) 
   Int start
  )
  {

    if(m > n)
    {
      printf("ERROR permtue \n");
      return BASKER_ERROR;
    }

    ENTRY_1DARRAY temp;
    MALLOC_ENTRY_1DARRAY(temp, n);
    init_value(temp, n, (Entry) 0.0);

    //Permute
    for(Int i = 0; i < n; i++)
    {
      temp(i) = vec(p(i));
    }
    //Copy back
    for(Int i = 0; i < n; i++)
    {
      vec(i) = temp(i);
    }

    FREE_ENTRY_1DARRAY(temp);

    return 0;
  }
  

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute
  (
   BASKER_MATRIX &M,
   INT_1DARRAY row,
   INT_1DARRAY col
   )
  {
    permute_col(M,col);
    permute_row(M,row);
    return 0;
  }//end permute(int, int)


  // NDE - added functions for solve performance improvements
  // This function combines two steps: inverse permutation of the rhs and copying to local views
  // The structure of the solve codes expect that x (at function call and at end as lhs) 
  // acts as the rhs during the solve routines, and that y (at function call the rhs) 
  // stores the lhs (until the final ND block solves in any case) and is init to 0 
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_inv_and_init_for_solve
  (
   Entry* y,
   ENTRY_1DARRAY &xcon,
   ENTRY_1DARRAY &ycon,
   INT_1DARRAY  &p, 
   Int n
  )
  {
    //Permute
    for(Int i = 0; i < n; i++) {
      xcon(p(i))  = y[i];
      ycon(i)     = (Entry) 0.0;
    }
    return 0;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_and_finalcopy_after_solve
  (
   Entry* x,
   ENTRY_1DARRAY &xconv,
   ENTRY_1DARRAY &yconv,
   INT_1DARRAY  &p,
   Int n
  )
  {
    /*
    // Steps from original code - this works with ND
    for(Int i = btf_tabs(btf_tabs_offset); i < n; i++) // final step from serial_btf_solve
    {
      xconv(i) = yconv(i);
    }

    for(Int i = 0; i < n; i++) //perm xconv back to original ordering and copy back to raw lhs pointer
    { x[i] = xconv(p(i)); }
    */

    const Int poffset = btf_tabs(btf_tabs_offset);
    for(Int i = 0; i < n; i++) //perm xconv back to original ordering and copy back to raw lhs pointer
    { 
      Int permi = p(i);
      if ( permi < poffset )
      {
      // ND blocks
        //x[i] = xconv(p(i)); 
        x[i] = xconv(permi);
      } 
      else {
      // btf blocks
        //x[i] = yconv(p(i)); 
        x[i] = yconv(permi); 
      }
    }

    return 0;
  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_and_init_for_solve
  (
   Entry* y,
   ENTRY_1DARRAY &xcon,
   ENTRY_1DARRAY &ycon,
   INT_1DARRAY  &p, 
   Int n
  )
  {
    //Permute
    for(Int i = 0; i < n; i++) {
      xcon(i)  = y[p(i)];
      ycon(i)  = (Entry) 0.0;
    }
    return 0;
  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_inv_and_finalcopy_after_solve
  (
   Entry* x,
   ENTRY_1DARRAY &xconv,
   ENTRY_1DARRAY &yconv,
   INT_1DARRAY  &p,
   Int n
  )
  {

    const Int poffset = btf_tabs(btf_tabs_offset);
    // pre-offset indices of xconv solution in ND block partition
    // >= poffset indices of yconv solution in small BTF block partition
    for(Int i = 0; i < n; i++) //perm xconv back to original ordering and copy back to raw lhs pointer
    { 
      Int permi = p(i);
      if ( i < poffset )
      {
      // ND blocks
        x[permi] = xconv(i); 
      } 
      else {
      // small btf blocks
        x[permi] = yconv(i); 
      }
    }

    return 0;
  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute
  (
   Entry* vec,
   INT_1DARRAY   p,
   Int n
  )
  {
    ENTRY_1DARRAY temp;
    MALLOC_ENTRY_1DARRAY(temp, n);
    init_value(temp, n, (Entry) 0.0);

    //Permute
    for(Int i = 0; i < n; i++)
    {
      temp(i) = vec[p(i)];
    }
    //Copy back
    for(Int i = 0; i < n; i++)
    {
      vec[i] = temp(i);
    }
    FREE_ENTRY_1DARRAY(temp);
    return 0;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute
  (
   INT_1DARRAY vec, 
   INT_1DARRAY   p, 
   Int n
  )
  {
    INT_1DARRAY temp;
    MALLOC_INT_1DARRAY(temp, n);
    init_value(temp, n, (Int) 0);

    //Permute
    for(Int i = 0; i < n; i++)
    {
      temp(i) = vec(p(i));
    }
    //Copy back
    for(Int i = 0; i < n; i++)
    {
      vec(i) = temp(i);
    }
    FREE_INT_1DARRAY(temp);
    return 0;
  }


  // NDE: permute_with_workspace
  //      this routine uses a pre-allocated array 
  //      to serve as the 'temp' array during permutation
  //      the workspace size is gn; may want to force that here, 
  //      or init value < gn - n (input) to 0
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_with_workspace
  (
   INT_1DARRAY & vec,
   INT_1DARRAY & p,
   Int n,
   Int istart,
   Int offset
  )
  {
    //Permute
    for(Int i = 0; i < n; i++)
    {
      perm_comp_iworkspace_array(i) = vec(istart+p(i+offset));
    }
    //Copy back
    for(Int i = 0; i < n; i++)
    {
      vec(istart+i) = perm_comp_iworkspace_array(i);
    }
    return BASKER_SUCCESS;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::permute_inv_with_workspace
  (
   INT_1DARRAY & vec,
   INT_1DARRAY & p,
   Int n,
   Int istart,
   Int offset
  )
  {
    //Permute
    for(Int i = 0; i < n; ++i)
    {
      perm_comp_iworkspace_array(p(i+offset)) = vec(istart+i);
    }
    //Copy back
    for(Int i = 0; i < n; ++i)
    {
      vec(istart+i) = perm_comp_iworkspace_array(i);
    }

    return BASKER_SUCCESS;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_array_with_workspace
  (
   Entry* vec,
   INT_1DARRAY & p,
   Int n
  )
  {
    //Permute
    for(Int i = 0; i < n; i++)
    {
      perm_comp_fworkspace_array(i) = vec[p(i)];
    }
    //Copy back
    for(Int i = 0; i < n; i++)
    {
      vec[i] = perm_comp_fworkspace_array(i);
    }
    return BASKER_SUCCESS;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_with_workspace
  (
   ENTRY_1DARRAY & vec,
   INT_1DARRAY & p,
   Int n
  )
  {
    //Permute
    for(Int i = 0; i < n; i++)
    {
      perm_comp_fworkspace_array(i) = vec(p(i));
    }
    //Copy back
    for(Int i = 0; i < n; i++)
    {
      vec(i) = perm_comp_fworkspace_array(i);
    }
    return BASKER_SUCCESS;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::permute_inv_array_with_workspace
  (
   Entry* vec,
   INT_1DARRAY & p,
   Int n
  )
  {
    //Permute
    for(Int i = 0; i < n; ++i)
    {
      perm_comp_fworkspace_array(p(i)) = vec[i];
    }
    //Copy back
    for(Int i = 0; i < n; ++i)
    {
      vec[i] = perm_comp_fworkspace_array(i);
    }

    return BASKER_SUCCESS;
  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::permute_inv_with_workspace
  (
   ENTRY_1DARRAY & vec,
   INT_1DARRAY & p,
   Int n
  )
  {
    //Permute
    for(Int i = 0; i < n; ++i)
    {
      perm_comp_fworkspace_array(p(i)) = vec(i);
    }
    //Copy back
    for(Int i = 0; i < n; ++i)
    {
      vec(i) = perm_comp_fworkspace_array(i);
    }

    return BASKER_SUCCESS;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::permute_composition_for_solve( Int ncols )
  {
    Int gn_ = ncols;
    // This macro is to test a minor improvement to this routine
    #define REDUCE_PERM_COMP_ARRAY 1

    for(Int i = 0; i < gn_; i++) // Move this to a different init function
    {
      perm_inv_comp_array(i) = i; 
      #ifndef REDUCE_PERM_COMP_ARRAY
      perm_comp_array(i) = i;
      #endif
    }

    // perm_inv_comp_array calculation
    //
    // This replaces the inverse perms in solve_interface
    // i.e. (...o p2^-1 o p1^-1)(v) = (p1 o p2 o ... )^-1(v)
    // BUT only p1 o p2 o ... is computed here and stored in perm_inv_comp_array 
    // In solve_interface, the inverse of this composition will be applied to v
    //
    // Note: The order had to be reversed as well as the type of permutation call
    //
    // p6 inv
    //      permute(perm_inv_comp_array, gperm, gn_);

    // NDE: This is moved out and initialized once - reused in 4 possible permutations
    using range_type = Kokkos::pair<int, int>;
    Int offset_amd = gn_;
    Int offset_nd  = gn_+gn_;
    auto iwork_amd = Kokkos::subview(perm_comp_iworkspace_array,
                                     range_type(offset_amd, offset_amd+gn_));
    auto iwork_nd  = Kokkos::subview(perm_comp_iworkspace_array,
                                     range_type(offset_nd, offset_nd+gn_));

    MALLOC_INT_1DARRAY(numeric_row_iperm_array, gn_);
    MALLOC_INT_1DARRAY(numeric_col_iperm_array, gn_);

    MALLOC_INT_1DARRAY(symbolic_row_iperm_array, gn_);
    MALLOC_INT_1DARRAY(symbolic_col_iperm_array, gn_);
    for(Int i = 0; i < gn_; ++i)
    {
      symbolic_col_iperm_array(i) = i;
      symbolic_row_iperm_array(i) = i;
    }

    // p5 inv
    if(amd_flag == BASKER_TRUE)
    {
      //printVec("amd.txt",order_csym_array, gn_);
      Int scol_top = btf_tabs[btf_top_tabs_offset];
      for(Int i = 0; i < scol_top; ++i) {
        iwork_amd(i) = i;
      }
      for(Int i = 0; i < BTF_A.ncol; ++i) {
        iwork_amd(scol_top + i) = scol_top + order_csym_array(i);
      }
      for(Int i = scol_top+BTF_A.ncol; i < gn_; ++i) {
        iwork_amd(i) = i;
      }
      //for(Int i = 0; i < gn_; ++i) {
      //  printf( " + perm[%d] = %d\n",i,iwork_amd(i));
      //}
      //if(Options.static_delayed_pivot == 0)
      {
        if(Options.verbose == BASKER_TRUE) {
          printf(" > CSYM AMD order in A(scol=%d, ncol=%d) and gn=%d\n", (int)BTF_A.scol, (int)BTF_A.ncol, (int)gn_);
        }
        permute_with_workspace(perm_inv_comp_array, perm_comp_iworkspace_array, gn_, 0, offset_amd);
        //for(Int i = 0; i < gn_; i++) printf( " %d, %d\n",i,perm_inv_comp_array(i) );
        //printf("\n");
      }
    }
    // p4 inv
    if(nd_flag == BASKER_TRUE)
    {
      Int scol_top = btf_tabs[btf_top_tabs_offset];
      //printVec("nd.txt", part_tree.permtab, gn_);
      for(Int i = 0; i < scol_top; ++i) {
        iwork_nd(i) = i;
      }
      for(Int i = 0; i < BTF_A.ncol; ++i) {
        iwork_nd(scol_top + i) = scol_top + part_tree.permtab(i);
      }
      for(Int i = scol_top+BTF_A.ncol; i < gn_; ++i) {
        iwork_nd(i) = i;
      }
      //for(Int i = 0; i < gn_; ++i) {
      //  printf( " + perm[%d] = %d\n",i,iwork_nd(i));
      //}
      //if(Options.static_delayed_pivot == 0)
      {
        if(Options.verbose == BASKER_TRUE) {
          printf(" > ND order in \n");
        }
        permute_with_workspace(perm_inv_comp_array, perm_comp_iworkspace_array, gn_, 0, offset_nd);
        //for(Int i = 0; i < gn_; i++) printf( " %d, %d\n",i,perm_inv_comp_array(i) );
        //printf("\n");
      }
    }
    // p2, p3 inv
    if(btf_flag == BASKER_TRUE)
    {
      // even if blk_matching is enabled, blk_amd & blk_mwm are still applied 
      //p3
      if(Options.verbose == BASKER_TRUE) {
        printf(" > blk_amd order in\n");
        //for(Int i = 0; i < gn_; ++i) {
        //  printf( " perm[%d] = %d\n",i,order_blk_amd_array(i));
        //}
      }
      permute_with_workspace(perm_inv_comp_array, order_blk_amd_array, gn_);
      //for(Int i = 0; i < gn_; i++) printf( " %d, %d\n",i,perm_inv_comp_array(i) );
      //printf("\n");

      if(Options.verbose == BASKER_TRUE) {
        printf(" > blk_mwm order in\n");
        //for(Int i = 0; i < gn_; ++i) {
        //  printf( " perm[%d] = %d\n",i,order_blk_mwm_array(i));
        //}
      }
      permute_with_workspace(perm_inv_comp_array, order_blk_mwm_array, gn_);
      //for(Int i = 0; i < gn_; i++) printf( " %d, %d\n",i,perm_inv_comp_array(i) );
      //printf("\n");

      #if 0 // should go to numeric?
      if(Options.static_delayed_pivot != 0)
      {
        if(amd_flag == BASKER_TRUE) {
          if(Options.verbose == BASKER_TRUE) {
            printf(" >> CSYM AMD order in A(scol=%d, ncol=%d) and gn=%d, with delayed\n",BTF_A.scol,BTF_A.ncol,gn_);
          }
          for(Int i = 0; i < gn_; ++i) {
            printf( " + perm[%d + %d] = %d\n",offset_amd,i,perm_comp_iworkspace_array(offset_amd + i));
          }
          permute_with_workspace(perm_inv_comp_array, perm_comp_iworkspace_array, gn_, offset_amd);
        }
        // p4 inv
        if(nd_flag == BASKER_TRUE) {
          if(Options.verbose == BASKER_TRUE) {
            printf(" >> ND order in, with delayed\n");
          }
          for(Int i = 0; i < gn_; ++i) {
            printf( " + perm[%d + %d] = %d\n",offset_nd,i,perm_comp_iworkspace_array(offset_nd + i));
          }
          permute_with_workspace(perm_inv_comp_array, perm_comp_iworkspace_array, gn_, offset_nd);
        }
      }
      #endif
      //for (int i=0; i<(int)perm_inv_comp_array.extent(0); i++) printf( "%d %d\n",i,perm_inv_comp_array(i) );
      //printVec("btf.txt", order_btf_array, gn_);
      //p2
      if(Options.verbose == BASKER_TRUE) {
        printf(" > btf order in\n");
        //for(Int i = 0; i < gn_; ++i) {
        //  printf( " perm[%d] = %d\n",i,order_btf_array(i));
        //}
      }
      permute_with_workspace(perm_inv_comp_array, order_btf_array, gn_);
      //for(Int i = 0; i < gn_; i++) printf( " %d, %d\n",i,perm_inv_comp_array(i) );
      //printf("\n");

      // compose perm for symbolic phase
      permute_with_workspace(symbolic_row_iperm_array, order_btf_array, gn_);
      permute_with_workspace(symbolic_col_iperm_array, order_btf_array, gn_);
    }

    // p1 inv
    if(match_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE) {
        printf(" > match order in\n");
        //for(Int i = 0; i < gn_; ++i) {
        //  printf( " perm[%d] = %d\n",i,order_match_array(i));
        //}
      }
      permute_with_workspace(perm_inv_comp_array, order_match_array, gn_);
      //for(Int i = 0; i < gn_; i++) printf( " %d, %d\n",i,perm_inv_comp_array(i) );
      //printf("\n");

      // compose perm for symbolic phase
      permute_with_workspace(symbolic_row_iperm_array, order_match_array, gn_);
    }
    MALLOC_INT_1DARRAY(symbolic_row_perm_array, gn_);
    MALLOC_INT_1DARRAY(symbolic_col_perm_array, gn_);
    for(Int i = 0; i < gn_; ++i)
    {
      symbolic_col_perm_array(i) = i;
      symbolic_row_perm_array(i) = i;
    }
    permute_inv(symbolic_row_perm_array, symbolic_row_iperm_array, gn_);
    permute_inv(symbolic_col_perm_array, symbolic_col_iperm_array, gn_);

    // ===============================================================================
    // perm_comp_array calculation
    // Note: don't need to inverse a row only perm
    // q1
    #if REDUCE_PERM_COMP_ARRAY
    for(Int i = 0; i < gn_; i++) // NDE: Move this to a different init function
    {
      perm_comp_array(i) = perm_inv_comp_array(i);
      //if (i != perm_comp_array(i)) printf( " perm_comp_array(%d) = %d\n",i,perm_comp_array(i) );
    }
    if(match_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE)
      {printf(" x match order out\n");}
      //printVec("match.txt", order_match_array, gn_);
      permute_inv_with_workspace(perm_comp_array, order_match_array, gn_);
      //for(Int i = 0; i < gn_; i++) if (i != perm_comp_array(i)) printf( " perm_comp_array(%d) = %d\n",i,perm_comp_array(i) );
    }
    if(btf_flag == BASKER_TRUE)
    {
      // invert btf first
      permute_inv_with_workspace(perm_comp_array, order_btf_array, gn_);
      #if 0 // should go to numeric?
      if (Options.static_delayed_pivot != 0) {
        if (nd_flag) {
          if(Options.verbose == BASKER_TRUE)
          { printf(" xx ND order out, delayed\n");}
          permute_inv_with_workspace(perm_comp_array, perm_comp_iworkspace_array, gn_, offset_nd);
        }
        if (nd_flag) {
          if(Options.verbose == BASKER_TRUE)
          { printf(" xx CSYM AMD order out, delayed\n");}
          permute_inv_with_workspace(perm_comp_array, perm_comp_iworkspace_array, gn_, offset_amd);
        }
      }
      #endif
      // now invert mwm
      if(Options.verbose == BASKER_TRUE)
      { printf(" x blk_mwm order out\n");}
      permute_inv_with_workspace(perm_comp_array, order_blk_mwm_array, gn_);

      // finally reapply btf
      #if 0 // should go to numeric?
      if (Options.static_delayed_pivot != 0) {
        if (nd_flag) {
          permute_with_workspace(perm_comp_array, perm_comp_iworkspace_array, gn_, offset_amd);
        }
        if (nd_flag) {
          permute_with_workspace(perm_comp_array, perm_comp_iworkspace_array, gn_, offset_nd);
        }
      }
      #endif
      permute_with_workspace(perm_comp_array, order_btf_array, gn_);
      //for(Int i = 0; i < gn_; i++) if (i != perm_comp_array(i)) printf( " perm_comp_array(%d) = %d\n",i,perm_comp_array(i) );
    }
    #else
    if(amd_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE)
      {printf("CSYM AMD order out \n");}
      //printVec(order_csym_array, gn_);
      for(Int i = 0; i < BTF_A.ncol; ++i)
      {
        perm_comp_iworkspace_array(i) = order_csym_array(i);
      }
      permute_with_workspace(perm_comp_array,perm_comp_iworkspace_array, gn_);
    }
    // q2
    if(nd_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE)
      {printf("ND order out \n");}
      //printVec(part_tree.permtab, gn_);
      for(Int i = 0; i < BTF_A.ncol; ++i)
      {
        perm_comp_iworkspace_array(i) = part_tree.permtab(i);
      }
      permute_with_workspace(perm_comp_array,perm_comp_iworkspace_array, gn_);
    }
    // q3, q4
    if(btf_flag == BASKER_TRUE)
    {
      //printVec(order_btf_array, gn_);
      // q3
      if(Options.verbose == BASKER_TRUE)
      {printf("blk_amd order out\n");}
      permute_with_workspace(perm_comp_array,order_blk_amd_array, gn_);
      // q4
      if(Options.verbose == BASKER_TRUE)
      {printf("btf order out\n");}
      permute_with_workspace(perm_comp_array,order_btf_array, gn_);
    }
    #endif
    #undef REDUCE_PERM_COMP_ARRAY
  }
  //NDE end


   /*GOBACK and make this good*/
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_col_store_valperms
  (
   BASKER_MATRIX &M,
   INT_1DARRAY &col,
   INT_1DARRAY &vals_order_perm
  )
  {
    if((M.ncol == 0)||(M.nnz == 0))
    { return 0; }

    Int n = M.ncol;
    Int num_col = (Int)(col.extent(0));
    if (n > num_col) {
      // ** we'll just permute the first col.extent(0) or M.ncol columns **
      // (e.g., E may have more columns than perm)
      if(Options.verbose == BASKER_TRUE) {
        printf( " > note: permute_col only %d columns since perm has only %d columns (%d:%d))\n", (int)num_col, (int)num_col, 0, (int)n);
      }
      n = num_col;
    }
    Int nnz = M.col_ptr(n);

    INT_1DARRAY temp_p;
    INT_1DARRAY temp_i;
    ENTRY_1DARRAY temp_v;

    MALLOC_INT_1DARRAY(temp_p, n+1);
    MALLOC_INT_1DARRAY(temp_i, nnz);
    MALLOC_ENTRY_1DARRAY(temp_v, nnz);

    //Determine column ptr of output matrix
    Kokkos::parallel_for(
      "permute_col", n,
      KOKKOS_LAMBDA(const int j) {
        Int i = col (j);
        temp_p (i+1) = M.col_ptr (j+1) - M.col_ptr (j);
      });
    Kokkos::fence();

    //Get ptrs from lengths
    temp_p (0) = 0;
    for(Int j = 0; j < n; j++)
    {
      temp_p (j+1) = temp_p (j+1) + temp_p (j);
    }

    //copy idxs
    Kokkos::parallel_for(
      "permute_col", n,
      KOKKOS_LAMBDA(const int ii) {
        Int ko = temp_p (col (ii) );
        for(Int k = M.col_ptr (ii); k < M.col_ptr (ii+1); k++)
        {
          temp_i (ko) = M.row_idx (k); // preserves order of indices in row_idx (they may not be ordered, but won't be reshuffled)
          temp_v (ko) = M.val (k);     // and corresponding order with vals

          vals_order_perm(k) = ko; //this works for first perm or individual; how to best compose subsequent perms? track separately and compose after?
          ko++;
        }
      });
    Kokkos::fence();
    // copy the indexes for the remaining columns
    Kokkos::parallel_for(
      "permute_col", Kokkos::RangePolicy<Exe_Space> (n, M.ncol),
      KOKKOS_LAMBDA(const int ii) {
        for(Int k = M.col_ptr (ii); k < M.col_ptr (ii+1); k++)
        {
          vals_order_perm(k) = k;
        }
      });

    //copy back into A
    Kokkos::parallel_for(
      "permute_col", n+1,
      KOKKOS_LAMBDA(const int ii) {
        M.col_ptr (ii) = temp_p (ii);
      });
    Kokkos::parallel_for(
      "permute_col", nnz,
      KOKKOS_LAMBDA(const int ii) {
        M.row_idx (ii) = temp_i (ii);
        M.val (ii) = temp_v (ii);
      });
    Kokkos::fence();

    FREE_INT_1DARRAY(temp_p);
    FREE_INT_1DARRAY(temp_i);
    FREE_ENTRY_1DARRAY(temp_v);

    return 0;
  }//end permute_col(int) 

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_col
  (
   BASKER_MATRIX &M,
   INT_1DARRAY col,
   Int frow
  )
  {
    if((M.ncol == 0)||(M.nnz == 0))
    { return 0; }

    Int n = M.ncol;
    Int num_col = (Int)(col.extent(0));
    if (n > num_col) {
      // ** we'll just permute the first col.extent(0) or M.ncol columns **
      // (e.g., E may have more columns than perm)
      if(Options.verbose == BASKER_TRUE) {
        printf( " > note: permute_col only %d columns since perm has only %d columns (%d:%d))\n", (int)num_col, (int)num_col, (int)frow, (int)(frow+n-1));
      }
      n = num_col;
    }
    Int nnz = M.col_ptr(frow+n) - M.col_ptr(frow);
    //printf( " ptr(%d) = %d, ptr(%d) = %d -> nnz = %d (M.n = %d, M.nnz = %d)\n",frow,M.col_ptr(frow),frow+n,M.col_ptr(frow+n),nnz,M.ncol,M.nnz );
    if (nnz == 0) {
      return 0;
    }

    INT_1DARRAY temp_p;
    INT_1DARRAY temp_i;
    ENTRY_1DARRAY temp_v;

    MALLOC_INT_1DARRAY(temp_p, n+1);
    MALLOC_INT_1DARRAY(temp_i, nnz);
    MALLOC_ENTRY_1DARRAY(temp_v, nnz);

    //init_value(temp_p, n+1, (Int)0);
    //init_value(temp_i, nnz, (Int)0);
    //init_value(temp_v, nnz, (Entry)0.0);

    //Determine column ptr of output matrix
    for(Int j = frow; j < frow+n; j++)
    {
      Int i = col (j-frow);
      temp_p (i+1) = M.col_ptr (j+1) - M.col_ptr (j);
    }
    //Get ptrs from lengths
    temp_p (0) = 0;
    for(Int j = 0; j < n; j++)
    {
      temp_p (j+1) = temp_p (j+1) + temp_p (j);
    }

    //copy idxs
    Kokkos::parallel_for(
      "permute_col", Kokkos::RangePolicy<Exe_Space> (frow, frow+n),
      KOKKOS_LAMBDA(const int ii) {
        Int ko = temp_p (col (ii-frow));
        for(Int k = M.col_ptr (ii); k < M.col_ptr (ii+1); k++)
        {
          temp_i (ko) = M.row_idx (k); // preserves order of indices in row_idx (they may not be ordered, but won't be reshuffled)
          temp_v (ko) = M.val (k);     // and corresponding order with vals
          ko++;
        }
      });
    Kokkos::fence();

    //copy back into A
    Kokkos::parallel_for(
      "permute_col", Kokkos::RangePolicy<Exe_Space> (frow, frow+n),
      KOKKOS_LAMBDA(const int ii) {
        M.col_ptr (ii+1) = M.col_ptr (frow) + temp_p (ii-frow+1);
      });
    Kokkos::fence();

    Kokkos::parallel_for(
      "permute_col", nnz,
      KOKKOS_LAMBDA(const int ii) {
        M.row_idx (M.col_ptr (frow) + ii) = temp_i (ii);
        M.val (M.col_ptr (frow) + ii) = temp_v (ii);
      });
    Kokkos::fence();

    FREE_INT_1DARRAY(temp_p);
    FREE_INT_1DARRAY(temp_i);
    FREE_ENTRY_1DARRAY(temp_v);

    return 0;
  }//end permute_col(int) 
  

  /*GOBACK and make this good*/
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::permute_row
  (
   Int nnz,
   Int *row_idx,
   Int *row
  )
  {
    if(nnz == 0)
    { return 0; }

    //permute
    Kokkos::parallel_for(
      "permute_row", nnz,
      KOKKOS_LAMBDA(const int &k) {
        row_idx[k] = row[row_idx[k]];
      });
    Kokkos::fence();

    return 0;
  }//end permute_row(matrix,int)

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::permute_row
  (
   BASKER_MATRIX &M,
   INT_1DARRAY row
  )
  {
    permute_row(M.nnz, &(M.row_idx(0)), &(row(0)));
    return 0;
  }//end permute_row(matrix,int)

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::sort_matrix_store_valperms( BASKER_MATRIX &M, INT_1DARRAY &order_vals_perms )
  {
    if(M.nnz == 0)
    { return 0; }

    #define USE_WORKSPACE_FOR_SORT
    #if defined(USE_WORKSPACE_FOR_SORT)
    INT_1DARRAY perm;
    MALLOC_INT_1DARRAY (perm, M.nrow);

    ENTRY_1DARRAY dwork;
    MALLOC_ENTRY_1DARRAY (dwork, M.nrow);
    INT_1DARRAY iwork;
    MALLOC_INT_1DARRAY (iwork, M.nrow);
    #endif

    #ifdef BASKER_TIMER
    std::cout << " ~~ sort_matrix_store_valperms(n = " << M.ncol << " nnz = " << M.nnz << ")" << std::endl;
    #endif
    //std::cout << " M = [" << std::endl;
    #ifdef BASKER_USE_QSORT
    Int max_nnz_row = 0;
    for(Int k = 0; k < M.ncol; k++) {
      Int nnz_row = M.col_ptr[k+1] - M.col_ptr[k];
      if (max_nnz_row < nnz_row) max_nnz_row = nnz_row;
    }
    basker_sort_pair *sort_pairs = (basker_sort_pair*)malloc(max_nnz_row * sizeof(basker_sort_pair));
    #endif
    for(Int k = 0; k < M.ncol; k++)
    {
      #ifdef BASKER_TIMER
      Kokkos::Timer timer_order;
      timer_order.reset();
      #endif
      Int start_row = M.col_ptr[k];
      #if defined(USE_WORKSPACE_FOR_SORT)
      perm(0) = start_row;
      #endif
      #ifdef BASKER_USE_QSORT
      {
        Int nnz = M.col_ptr[k+1] - M.col_ptr[k];
        for (Int i = 0; i < nnz; i++) {
          sort_pairs[i].perm = start_row + i;
          sort_pairs[i].rowid = M.row_idx[start_row + i];
        }
        qsort(sort_pairs, (size_t) nnz, sizeof(basker_sort_pair),
              &basker_sort_matrix_col);
        for (Int i = 0; i < nnz; i++) {
          perm(i) = sort_pairs[i].perm;
        }
      }
      #else
      for(Int i = M.col_ptr[k]+1; i < M.col_ptr[k+1]; i++)
      {
        #if 1
          // binary search
          Int jj_start = start_row;
          Int jj_end   = i;
          Int id = M.row_idx[i];
          while(jj_end > jj_start)
          {
            Int jj_mid = jj_start + (jj_end-jj_start)/2;
            #if defined(USE_WORKSPACE_FOR_SORT)
            Int id_mid = M.row_idx[perm(jj_mid-start_row)];
            #else
            Int id_mid = M.row_idx[jj_mid];
            #endif
            if (id_mid == id) {
                break;
            } else if (id_mid > id) {
                jj_end = jj_mid-1;
            } else {
                jj_start = jj_mid+1;
            }
          }
          Int jj = jj_start;
          #if defined(USE_WORKSPACE_FOR_SORT)
          if (jj < i && M.row_idx[perm(jj-start_row)] < id) 
          #else
          if (jj < i && M.row_idx[jj] < id) 
          #endif
          {
            jj++;
          }
        #else
          //just use insertion sort
          Int jj = i;
          while((jj > start_row) && (M.row_idx[jj-1] > M.row_idx[i]))
          {
            jj = jj-1;
          }   //end while jj
        #endif
        // insert i at jj
        // NOTE: the complexity of "shift" depends on the inputs 
        //       e.g., if already sorted, then no shift
        #if !defined(USE_WORKSPACE_FOR_SORT)
        Int   t_row_idx = M.row_idx[i];
        Entry t_row_val = M.val[i];
        Int   tmp_index = order_vals_perms(i);
        #endif
        // > shift to make space at jj
        for (Int j = i; j > jj; j--) {
            #if defined(USE_WORKSPACE_FOR_SORT)
            perm(j-start_row) = perm(j-start_row-1);
            #else
            M.row_idx[j] = M.row_idx[j-1];
            M.val[j]     = M.val[j-1];
            order_vals_perms(j) = order_vals_perms(j-1);
            #endif
        }
        // > insert i at jj
        #if defined(USE_WORKSPACE_FOR_SORT)
        perm(jj-start_row) = i;
        #else
        M.row_idx[jj] = t_row_idx;
        M.val[jj]     = t_row_val;
        order_vals_perms(jj) = tmp_index;
        #endif
      } //end for i
      #endif

      #if defined(USE_WORKSPACE_FOR_SORT)
      for(Int i = 0; i < M.col_ptr[k+1]-start_row; i++)
      {
          //printf( " perm[%d] = %d, row_idx[%d] = %d\n",(int)i,(int)perm(i),(int)perm(i),(int)M.row_idx[perm(i)] );
          iwork[i] = order_vals_perms(perm(i));
          dwork[i] = M.val[perm(i)];
          perm(i)  = M.row_idx[perm(i)];
      }
      //printf( "\n" );
      for(Int i = 0; i < M.col_ptr[k+1]-start_row; i++)
      {
          Int id = start_row+i;
          order_vals_perms(id) = iwork[i];
          M.val[id]     = dwork[i];
          M.row_idx[id] = perm(i);
      }
      #endif
      //for(Int i = M.col_ptr[k]; i < M.col_ptr[k+1]; i++) {
      //    std::cout << k << " " << M.row_idx[i] << " " << M.val[i] << std::endl;
      //}
      #ifdef BASKER_TIMER
      double order_time = timer_order.seconds();
      if (order_time > 1) {
        std::cout << "    > Basker sort time col=" << k << " : " << order_time << " (nnz = " << M.col_ptr[k+1]-M.col_ptr[k] << ")" <<std::endl;
      }
      #endif
    }//end over all columns k
    #ifdef BASKER_USE_QSORT
    free(sort_pairs);
    #endif

    return 0;
  }//end sort_matrix()


  //--------------------------------------------------------------
  template <class Int, class Entry>
  struct kokkos_nd_sorter
  {
    Int nids;
    Int nrow;
    Int nblks;
    INT_1DARRAY   nd_map;
    INT_1DARRAY   order_vals_perms;

    INT_1DARRAY   col_ptr;
    INT_1DARRAY   row_idx;
    ENTRY_1DARRAY val;

    INT_1DARRAY   nd_ptr_global;
    INT_1DARRAY   perm_global;

    INT_1DARRAY   iwork_global;
    ENTRY_1DARRAY dwork_global;

    kokkos_nd_sorter(Int _nids,
                     Int _nrow,
                     Int _nblks,
                     //
                     INT_1DARRAY   _nd_map,
                     INT_1DARRAY   _order_vals_perms,
                     INT_1DARRAY   _col_ptr,
                     INT_1DARRAY   _row_idx,
                     ENTRY_1DARRAY _val,
                     //
                     INT_1DARRAY   _nd_ptr_global,
                     INT_1DARRAY   _perm_global,
                     INT_1DARRAY   _iwork_global,
                     ENTRY_1DARRAY _dwork_global) :
    nids(_nids),
    nrow(_nrow),
    nblks(_nblks),
    nd_map(_nd_map),
    order_vals_perms(_order_vals_perms),
    col_ptr(_col_ptr),
    row_idx(_row_idx),
    val(_val),
    nd_ptr_global(_nd_ptr_global),
    perm_global(_perm_global),
    iwork_global(_iwork_global),
    dwork_global(_dwork_global)
    {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int id) const {
      //#define BASKER_TIMER_AMD_FUNCTOR
      #if defined(BASKER_TIMER_AMD_FUNCTOR)
      Kokkos::Timer timer1;
      Kokkos::Timer timer_id;
      double time1 = 0.0;
      double time2 = 0.0;
      double time3 = 0.0;
      #endif

      Int offset = id*(1+nblks);
      using range_type = Kokkos::pair<int, int>;
      auto nd_ptr = Kokkos::subview(nd_ptr_global,
                                    range_type(offset, offset+nblks+1));
      offset = id*nrow;
      auto perm  = Kokkos::subview(perm_global,
                                   range_type(offset, offset+nrow));
      auto dwork = Kokkos::subview(dwork_global,
                                   range_type(offset, offset+nrow));
      auto iwork = Kokkos::subview(iwork_global,
                                   range_type(offset, offset+nrow));

      Int mloc = (nrow + nids - 1)/nids;
      Int start_k = id * mloc;
      Int end_k = start_k + mloc;
      if (end_k > nrow) {
        end_k = nrow;
      }
      #if defined(BASKER_TIMER_AMD_FUNCTOR)
      printf( " > sort subview time (%d) : %e seconds\n", id, timer_id.seconds() );
      timer_id.reset();
      printf( " %d: sort(k = %d : %d)\n",id,start_k,end_k-1 );
      #endif
      for(Int k = start_k; k < end_k; k++)
      {
        #if defined(BASKER_TIMER_AMD_FUNCTOR)
        timer1.reset();
        #endif
        // count nnz in each ND interior/separator
        for (Int i = 0; i <= nblks; i++) {
          nd_ptr(i) = 0;
        }
        // number of nz in the upper-part of diagonal block (for find_2D_convert, upper part needs to come before lower)
        Int num_upper = 0;
        Int col_id = nd_map(k);
        for(Int i = col_ptr[k]; i < col_ptr[k+1]; i++)
        {
          Int row_id = nd_map(row_idx[i]);
          nd_ptr(row_id+1) ++;
          if (row_id == col_id && row_idx[i] <= k) {
            // count nz in upper part of diagonal block
            num_upper ++;
          }
        }
        for (Int i = 0; i < nblks; i++) {
          nd_ptr(i+1) += nd_ptr(i);
        }
        #if defined(BASKER_TIMER_AMD_FUNCTOR)
        time1 += timer1.seconds();
        timer1.reset();
        #endif
        // sort into workspace
        Int diag_lower_ptr = nd_ptr(col_id) + num_upper;
        for(Int i = col_ptr[k]; i < col_ptr[k+1]; i++)
        {
          Int row_id = nd_map(row_idx[i]);
          Int idx = nd_ptr(row_id);
          if (row_id == col_id && row_idx[i] > k) {
            // shift for nz in lower part of diagonal block
            idx = diag_lower_ptr;
            diag_lower_ptr ++;
          } else {
            nd_ptr(row_id) ++;
          }
          perm[idx]   = order_vals_perms(i);
          dwork[idx]  = val[i];
          iwork[idx]  = row_idx[i];
        }
        #if defined(BASKER_TIMER_AMD_FUNCTOR)
        time2 += timer1.seconds();
        timer1.reset();
        #endif
        // copy from workspace
        Int start_row = col_ptr[k];
        for(Int i = 0; i < col_ptr[k+1] - col_ptr[k]; i++)
        {
          order_vals_perms(start_row + i) = perm[i];
          val[start_row + i]     = dwork[i];
          row_idx[start_row + i] = iwork[i];
        }
        #if defined(BASKER_TIMER_AMD_FUNCTOR)
        time3 += timer1.seconds();
        #endif
      }
      #if defined(BASKER_TIMER_AMD_FUNCTOR)
      double done_time = timer_id.seconds();
      printf( " > sort-functor done time (%d) : %e + %e + %e -> %e seconds\n", id, time1,time2,time3,done_time );
      #endif
    }
  };
  //--------------------------------------------------------------

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::ndsort_matrix_store_valperms( BASKER_MATRIX &M, INT_1DARRAY &order_vals_perms, BASKER_BOOL track_perm )
  {
    //#define BASKER_TIMER_AMD
    #ifdef BASKER_TIMER_AMD
    Kokkos::Timer timer_order;
    #endif
    Int nblks = tree.nblks;
    INT_1DARRAY nd_map;
    MALLOC_INT_1DARRAY (nd_map, M.nrow);

    //printf( " nrow = %d, nblks = %d\n",M.nrow,nblks );
    //for(Int k = 0; k <= nblks; k++) printf( "tree[%d] = %d\n",k,tree.row_tabs[k] );

    Int nids = num_threads;
    Kokkos::parallel_for(
      "ndsort_matrix_store_valperms", nids,
      KOKKOS_LAMBDA(const int id) {
        for (Int k = id; k < nblks; k += nids) {
          for (Int i = tree.row_tabs[k]; i < tree.row_tabs[k+1]; i++) {
            nd_map(i) = k;
          }
        }
      });
    Kokkos::fence();

    #ifdef BASKER_TIMER_AMD
    std::cout << std::endl << " > sort malloc time : " << timer_order.seconds() << std::endl;
    timer_order.reset();
    #endif
#if 1
    // sequential
    INT_1DARRAY nd_ptr;
    MALLOC_INT_1DARRAY (nd_ptr, nblks+1);

    INT_1DARRAY   perm;
    INT_1DARRAY   iwork;
    ENTRY_1DARRAY dwork;
    MALLOC_ENTRY_1DARRAY (dwork, M.nrow);
    MALLOC_INT_1DARRAY   (iwork, M.nrow);
    MALLOC_INT_1DARRAY   (perm,  M.nrow);
    #ifdef BASKER_TIMER_AMD_FUNCTOR
    printf( " > sort subview time (%d) : %e seconds\n", id, timer_id.seconds() );
    timer_id.reset();
    #endif
    for(Int k = 0; k < M.ncol; k ++)
    {
      #ifdef BASKER_TIMER_AMD_FUNCTOR
      timer1.reset();
      #endif
      // count nnz in each ND interior/separator
      for (Int i = 0; i <= nblks; i++) {
        nd_ptr(i) = 0;
      }
      // number of nz in the upper-part of diagonal block (for find_2D_convert, upper part needs to come before lower)
      Int num_upper = 0;
      Int col_id = nd_map(k);
      for(Int i = M.col_ptr[k]; i < M.col_ptr[k+1]; i++)
      {
        Int row_id = nd_map(M.row_idx[i]);
        nd_ptr(row_id+1) ++;
        if (row_id == col_id && M.row_idx[i] <= k) {
          // count nz in upper part of diagonal block
          num_upper ++;
        }
      }
      for (Int i = 0; i < nblks; i++) {
        nd_ptr(i+1) += nd_ptr(i);
      }
      #ifdef BASKER_TIMER_AMD_FUNCTOR
      time1 += timer1.seconds();
      timer1.reset();
      #endif
      // sort into workspace
      Int diag_lower_ptr = nd_ptr(col_id) + num_upper;
      for(Int i = M.col_ptr[k]; i < M.col_ptr[k+1]; i++)
      {
        Int row_id = nd_map(M.row_idx[i]);
        Int idx = nd_ptr(row_id);
        if (row_id == col_id && M.row_idx[i] > k) {
          // shift for nz in lower part of diagonal block
          idx = diag_lower_ptr;
          diag_lower_ptr ++;
        } else {
          nd_ptr(row_id) ++;
        }
        if (track_perm) {
          perm[idx]  = order_vals_perms(i);
        }
        dwork[idx]  = M.val[i];
        iwork[idx]  = M.row_idx[i];
      }
      #ifdef BASKER_TIMER_AMD_FUNCTOR
      time2 += timer1.seconds();
      timer1.reset();
      #endif
      // copy from workspace
      Int start_row = M.col_ptr[k];
      Int num_rows = M.col_ptr[k+1] - M.col_ptr[k];
      #if 0 // seems not enough work within each column to parallelize
      Kokkos::parallel_for(
        "permute_col", num_rows,
        KOKKOS_LAMBDA(const int i) {
          order_vals_perms(start_row + i) = perm[i];
          M.val[start_row + i]     = dwork[i];
          M.row_idx[start_row + i] = iwork[i];
        });
      Kokkos::fence();
      #else
      for(Int i = 0; i < num_rows; i++)
      {
        if (track_perm) {
          order_vals_perms(start_row + i) = perm[i];
        }
        M.val[start_row + i]     = dwork[i];
        M.row_idx[start_row + i] = iwork[i];
      }
      #endif
      #ifdef BASKER_TIMER_AMD_FUNCTOR
      time3 += timer1.seconds();
      #endif
    }
    #ifdef BASKER_TIMER_AMD_FUNCTOR
    printf( " > sort done time (%d) : %e + %e + %e -> %e seconds\n", id, time1,time2,time3,timer_id.seconds() );
    #endif
    FREE_INT_1DARRAY (nd_ptr);

    FREE_ENTRY_1DARRAY (dwork);
    FREE_INT_1DARRAY   (iwork);
    FREE_INT_1DARRAY   (perm);
#else
    INT_1DARRAY nd_ptr_global;
    MALLOC_INT_1DARRAY (nd_ptr_global, nids*(nblks+1));

    INT_1DARRAY perm_global;
    INT_1DARRAY iwork_global;
    ENTRY_1DARRAY dwork_global;
    MALLOC_ENTRY_1DARRAY (dwork_global, nids*(M.nrow));
    MALLOC_INT_1DARRAY   (iwork_global, nids*(M.nrow));
    MALLOC_INT_1DARRAY   (perm_global,  nids*(M.nrow));

    #if 0
    kokkos_nd_sorter<Int, Entry> sorter_functor(nids, M.nrow, nblks, nd_map, order_vals_perms, M.col_ptr, M.row_idx, M.val,
                                                nd_ptr_global, perm_global, iwork_global, dwork_global);
    Kokkos::parallel_for("ND SORTER on A", Kokkos::RangePolicy<Exe_Space>(0, nids), sorter_functor);
    #else
    using range_type = Kokkos::pair<int, int>;
    Kokkos::parallel_for(
      "ndsort_matrix_store_valperms", nids,
      KOKKOS_LAMBDA(const int id) {

        //#define BASKER_TIMER_AMD_FUNCTOR
        #ifdef BASKER_TIMER_AMD_FUNCTOR
        Kokkos::Timer timer1;
        Kokkos::Timer timer_id;
        double time1 = 0.0;
        double time2 = 0.0;
        double time3 = 0.0;
        #endif
        Int offset = id*(1+nblks);
        auto nd_ptr = Kokkos::subview(nd_ptr_global,
                                      range_type(offset, offset+nblks+1));
        offset = id*(M.nrow);
        auto perm  = Kokkos::subview(perm_global,
                                     range_type(offset, offset+M.nrow));
        auto dwork = Kokkos::subview(dwork_global,
                                     range_type(offset, offset+M.nrow));
        auto iwork = Kokkos::subview(iwork_global,
                                     range_type(offset, offset+M.nrow));

        Int mloc = (M.ncol + nids - 1)/nids;
        Int start_k = id * mloc;
        Int end_k = start_k + mloc;
        if (end_k > M.ncol) {
          end_k = M.ncol;
        }
        #ifdef BASKER_TIMER_AMD_FUNCTOR
        printf( " > sort subview time (%d) : %e seconds\n", id, timer_id.seconds() );
        timer_id.reset();
        printf( " %d: sort(k = %d : %d)\n",id,start_k,end_k-1 );
        #endif
        for(Int k = start_k; k < end_k; k++)
        {
          #ifdef BASKER_TIMER_AMD_FUNCTOR
          timer1.reset();
          #endif
          // count nnz in each ND interior/separator
          for (Int i = 0; i <= nblks; i++) {
            nd_ptr(i) = 0;
          }
          // number of nz in the upper-part of diagonal block (for find_2D_convert, upper part needs to come before lower)
          Int num_upper = 0;
          Int col_id = nd_map(k);
          for(Int i = M.col_ptr[k]; i < M.col_ptr[k+1]; i++)
          {
            Int row_id = nd_map(M.row_idx[i]);
            nd_ptr(row_id+1) ++;
            if (row_id == col_id && M.row_idx[i] <= k) {
              // count nz in upper part of diagonal block
              num_upper ++;
            }
          }
          for (Int i = 0; i < nblks; i++) {
            nd_ptr(i+1) += nd_ptr(i);
          }
          #ifdef BASKER_TIMER_AMD_FUNCTOR
          time1 += timer1.seconds();
          timer1.reset();
          #endif
          // sort into workspace
          Int diag_lower_ptr = nd_ptr(col_id) + num_upper;
          for(Int i = M.col_ptr[k]; i < M.col_ptr[k+1]; i++)
          {
            Int row_id = nd_map(M.row_idx[i]);
            Int idx = nd_ptr(row_id);
            if (row_id == col_id && M.row_idx[i] > k) {
              // shift for nz in lower part of diagonal block
              idx = diag_lower_ptr;
              diag_lower_ptr ++;
            } else {
              nd_ptr(row_id) ++;
            }
            perm[idx]   = order_vals_perms(i);
            dwork[idx]  = M.val[i];
            iwork[idx]  = M.row_idx[i];
          }
          #ifdef BASKER_TIMER_AMD_FUNCTOR
          time2 += timer1.seconds();
          timer1.reset();
          #endif
          // copy from workspace
          Int start_row = M.col_ptr[k];
          for(Int i = 0; i < M.col_ptr[k+1] - M.col_ptr[k]; i++)
          {
            order_vals_perms(start_row + i) = perm[i];
            M.val[start_row + i]     = dwork[i];
            M.row_idx[start_row + i] = iwork[i];
          }
          #ifdef BASKER_TIMER_AMD_FUNCTOR
          time3 += timer1.seconds();
          #endif
        }
        #ifdef BASKER_TIMER_AMD_FUNCTOR
        double done_time = timer_id.seconds();
        printf( " > sort done time (%d) : %e + %e + %e -> %e seconds\n", id, time1,time2,time3,done_time );
        #endif
      });
    #endif
    Kokkos::fence();
    FREE_INT_1DARRAY (nd_ptr_global);

    FREE_ENTRY_1DARRAY (dwork_global);
    FREE_INT_1DARRAY   (iwork_global);
    FREE_INT_1DARRAY   (perm_global);
#endif
    #ifdef BASKER_TIMER_AMD
    std::cout << " > sort time : " << timer_order.seconds() << std::endl;
    #endif

    FREE_INT_1DARRAY (nd_map);

    return 0;
  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::sort_matrix( Int nnz, Int ncol, Int *col_ptr, Int *row_idx, Entry *val )
  {
    if(nnz == 0)
    { return 0; }

    //just use insertion sort
    for(Int k = 0; k < ncol; k++)
    {
      Int start_row = col_ptr[k];
      for(Int i = col_ptr[k]+1; i < col_ptr[k+1]; i++)
      {
        Int jj = i;
        while((jj > start_row) && (row_idx[jj-1] > row_idx[jj]))
        {
          //swap
          Int   t_row_idx = row_idx[jj-1];
          Entry t_row_val = val[jj-1];

          row_idx[jj-1] = row_idx[jj];
          val[jj-1]     = val[jj];

          row_idx[jj] = t_row_idx;
          val[jj]     = t_row_val;

          jj = jj-1;
        } //end while jj
      } //end for i
    }//end over all columns k

    return 0;
  }//end sort_matrix()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::sort_matrix( BASKER_MATRIX &M )
  {
    int info = 0;
    info = sort_matrix( M.nnz, M.ncol,&(M.col_ptr(0)), &(M.row_idx(0)), &(M.val(0)) );
    return info;
  }

}//end namespace basker

#undef BASKER_TIMER
#endif //end ifndef basker_order_hpp
