#ifndef SHYLUBASKER_ORDER_HPP
#define SHYLUBASKER_ORDER_HPP

#include <impl/Kokkos_Timer.hpp>

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

    printf("Basker outer num_threads:%d \n", num_threads);
    MALLOC_INT_1DARRAY(btf_schedule, num_threads+1);
    init_value(btf_schedule, num_threads+1, 0);
    find_btf(A); 
   
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
      MALLOC_INT_1DARRAY(order_csym_array, BTF_A.ncol+1);
      init_value(order_csym_array, BTF_A.ncol+1,(Int) 0);

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

      if(btf_tabs_offset != 0)
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

    /*printf( " a=[\n" );
    for(Int j = 0; j < A.ncol; j++) {
      for(Int k = A.col_ptr[j]; k < A.col_ptr[j+1]; k++) {
        printf( "%d %d %e\n", A.row_idx[k], j, A.val[k]);
      }
    }
    printf( "];\n" );*/

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
    printf( " c=[\n" );
    for(Int j = 0; j < A.ncol; j++) {
      for(Int k = A.col_ptr[j]; k < A.col_ptr[j+1]; k++) {
        printf( "%d %d %e\n", A.row_idx[k], j, A.val[k]);
      }
    }
    printf( "];\n" );*/


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
    find_btf2(A);

    /*printf( " btfP=[\n" );
    for (int i = 0; i < A.ncol; i++) printf( "%d\n", order_btf_array(i));
    printf( "];\n" );*/

    #ifdef BASKER_TIMER
    order_time = timer_order.seconds();
    std::cout << " ++ Basker order : BTF time     : " << order_time << std::endl;
    timer_order.reset();
    #endif

    //================================================================
    //3. ND ordering of the big A block
    if (btf_tabs_offset != 0) // BTF_A exists and is not a btf_nblks > 1
    {
      // allocation perm vectors
      btfa_nnz = BTF_A.nnz;
      btfb_nnz = BTF_B.nnz;
      MALLOC_INT_1DARRAY(part_tree.permtab, BTF_A.ncol+1);
      MALLOC_INT_1DARRAY(order_csym_array,  BTF_A.ncol+1);
      MALLOC_INT_1DARRAY(order_csym_inv,    BTF_A.ncol+1);
      MALLOC_INT_1DARRAY(order_nd_inv,      BTF_A.ncol+1);
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
      for (Int i = 0; i <= BTF_A.ncol; ++i) {
        part_tree.permtab(i) = i;
      }
      for (Int i = 0; i <= BTF_A.ncol; ++i) {
        order_csym_array(i) = i;
      }

      // compute ND (if no MWM is requested)
      if (Options.blk_matching == 0) {
        apply_scotch_partition();
      }
    } //if btf_tabs_offset != 0

    #ifdef BASKER_TIMER
    order_time = timer_order.seconds();
    std::cout << " ++ Basker order : sort(A) time : " << order_time << std::endl;
    timer_order.reset();
    #endif

    //================================================================
    if(btf_nblks > 1) //else only BTF_A exists, A is assigned directly to it...
    {
      if ( btf_tabs_offset == 0 && BTF_C.nnz > 0 ) {
        // NDE: May need to add permutation for this case...
        //new for sfactor_copy2 replacement
        btfc_nnz = BTF_C.nnz;
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

    Kokkos::Impl::Timer timer_one;

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
    if((Options.btf == BASKER_TRUE) &&
        (btf_tabs_offset != 0))
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
        MALLOC_INT_1DARRAY(order_csym_array, BTF_A.ncol+1);
        init_value(order_csym_array, BTF_A.ncol+1,(Int) 0);
        csymamd_order(BTF_A, order_csym_array, cmember);

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
      if(Options.verbose == BASKER_TRUE) {
        std::cout << " ++ calling NO matching ++ " << std::endl;
      }
      match_flag = BASKER_FALSE;
    } else {
      match_flag = BASKER_TRUE;
      /*printf( " A = [\n" );
      for(Int j = 0; j < A.ncol; j++) {
        for(Int k = A.col_ptr[j]; k < A.col_ptr[j+1]; k++) {
          printf( " %d %d %e\n", A.row_idx[k],j,A.val[k]);
        }
      }
      printf("];\n");*/

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
          num_match = trilinos_btf_maxtrans (A.nrow, A.ncol, &(A.col_ptr(0)), &(A.row_idx(0)), maxwork,
                                             &work, &(order_match_array(0)), &(WORK(0)));
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
      /*for(Int j = 0; j < A.ncol; j++) printf( " > %d\n",order_match_array(j) );
      printf( " B = [\n" );
      for(Int j = 0; j < A.ncol; j++) {
        for(Int k = A.col_ptr[j]; k < A.col_ptr[j+1]; k++) {
          printf( " %d %d %e\n", A.row_idx[k],j,A.val[k]);
        }
      }
      printf("];\n");*/
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

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::apply_scotch_partition()
  {
    //--------------------------------------------------------------
    //3. ND on BTF_A
    //currently finds ND and permute BTF_A
    //Would like to change so finds permuation, 
    //and move into 2D-Structure
#if 0 // FIX: Do we need to sort?
    sort_matrix_store_valperms(BTF_A, vals_order_ndbtfa_array);
#endif

    /*printf(" in apply_scotch_partition\n" );
    printf(" T = [\n" );
    for(Int j = 0; j < BTF_A.ncol; j++) {
      for(Int k = BTF_A.col_ptr[j]; k < BTF_A.col_ptr[j+1]; k++) {
        printf("%d %d %.16e\n", BTF_A.row_idx[k], j, BTF_A.val[k]);
      }
    }
    printf("];\n");*/

    scotch_partition(BTF_A); // tree is prepped; permutation then applied to BTF_A
    /*printf(" After ND\n");
    printf(" p = [\n" );
    for(Int j = 0; j < BTF_A.ncol; j++) {
      printf("%d\n",part_tree.permtab(j) );
    }
    printf("];\n");
    printf(" ppT = [\n" );
    for(Int j = 0; j < BTF_A.ncol; j++) {
      for(Int k = BTF_A.col_ptr[j]; k < BTF_A.col_ptr[j+1]; k++) {
        printf("%d %d %.16e\n", BTF_A.row_idx[k], j, BTF_A.val[k]);
      }
    }
    printf("];\n");*/

    //need to do a row perm on BTF_B too
    if(btf_nblks > 1)
    {
      permute_row(BTF_B, part_tree.permtab);
    }
    //needed because  moving into 2D-Structure,
    //assumes sorted columns
    //sort_matrix(BTF_A);
#if 0 // FIX: do we need to sort (eg for cAMD)?
    sort_matrix_store_valperms(BTF_A, vals_order_ndbtfa_array); // BTF_A( perm(i) ) <- BTF_A( i )
    if(btf_nblks > 1)
    {
      //new for sfactor_copy2 replacement
      if ( BTF_B.nnz > 0 ) {
        sort_matrix_store_valperms(BTF_B, vals_order_ndbtfb_array);
      }
      if ( BTF_C.nnz > 0 ) {
        sort_matrix_store_valperms(BTF_C, vals_order_ndbtfc_array);
      }
    }
#endif
    //For debug
    //printMTX("A_BTF_PART_AFTER.mtx", BTF_A);

    //--------------------------------------------------------------
    //4. Init tree structure
    //This reduces the ND ordering into that fits,
    //thread counts
    init_tree_thread();

    //--------------------------------------------------------------
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
    init_value(order_csym_array, BTF_A.ncol+1,(Int) 0);
    //printf("============CALLING CSYMAMD============\n");
    csymamd_order(BTF_A, order_csym_array, cmember);
    #if 0 // reset to I for debug
    for (Int i = 0; i < BTF_A.ncol; ++i) {
      order_csym_array(i) = i;
    }
    #endif
    //permute_col(BTF_A, order_csym_array);
    //sort_matrix(BTF_A); // unnecessary?

    //new for sfactor_copy2 replacement
    permute_col_store_valperms(BTF_A, order_csym_array, vals_order_csym_array); //NDE: Track movement of vals (lin_ind of row,col) here
    if (Options.blk_matching == 0) {
      permute_inv(vals_order_ndbtfa_array, vals_order_csym_array, BTF_A.nnz);     //must permute the array holding the perms
    }
//    sort_matrix_store_valperms(BTF_A, vals_order_ndbtfa_array);

    permute_row(BTF_A, order_csym_array);

    /*printf(" After cAMD\n");
    printf(" p = [\n" );
    for(Int j = 0; j < BTF_A.ncol; j++) {
      printf("%d\n",order_csym_array(j) );
    }
    printf("];\n");
    printf(" ppT = [\n" );
    for(Int j = 0; j < BTF_A.ncol; j++) {
      for(Int k = BTF_A.col_ptr[j]; k < BTF_A.col_ptr[j+1]; k++) {
        printf("%d %d %.16e\n", BTF_A.row_idx[k], j, BTF_A.val[k]);
      }
    }
    printf("];\n");*/

    // ----------------------------------------------
    // sort matrix(BTF_A);
    //for (Int j = 0; j < BTF_A.nnz; j++) printf( " - vals_ndbtfa(%d) = %d\n",j,vals_order_ndbtfa_array(j) );
    sort_matrix_store_valperms(BTF_A, vals_order_ndbtfa_array); //new for replacement
    //printMTX("A_BTF_AMD1.mtx", BTF_A);
    /*for (Int j = 0; j < BTF_A.nnz; j++) printf( " + vals_ndbtfa(%d) = %d\n",j,vals_order_ndbtfa_array(j) );
    printf(" ppT = [\n" );
    for(Int j = 0; j < BTF_A.ncol; j++) {
      for(Int k = BTF_A.col_ptr[j]; k < BTF_A.col_ptr[j+1]; k++) {
        printf("%d %d %.16e\n", BTF_A.row_idx[k], j, BTF_A.val[k]);
      }
    }
    printf("];\n");*/

    if(btf_nblks > 1)
    {
      permute_row(BTF_B, order_csym_array);
      //new for sfactor_copy2 replacement
      if ( BTF_B.nnz > 0 ) {
        sort_matrix_store_valperms(BTF_B, vals_order_ndbtfb_array);
      }
      if ( BTF_C.nnz > 0 ) {
        sort_matrix_store_valperms(BTF_C, vals_order_ndbtfc_array);
      }
    }

    //new for sfactor_copy2 replacement
    // Steps below necessary for later use in Factor(...) to follow permutation convention used up to this point
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

    //--------------------------------------------------------------
    //6. Move to 2D Structure
    //finds the shapes for both view and submatrices,
    //need to be changed over to just submatrices
    matrix_to_views_2D(BTF_A);
    //finds the starting point of A for submatrices
    find_2D_convert(BTF_A);
    //now we can fill submatrices
    //printf("AFTER CONVERT\n");
    #ifdef BASKER_KOKKOS
    kokkos_order_init_2D<Int,Entry,Exe_Space> iO(this);
    Kokkos::parallel_for(TeamPolicy(num_threads,1), iO);
    Kokkos::fence();
    #else
    //Comeback
    #endif

    #if 1
    {
      // initialize threading  
      #ifdef BASKER_TIMER
      Kokkos::Timer timer_init;
      double barrier_init_time = 0.0;
      #endif

      basker_barrier.init(num_threads, 16, tree.nlvls );

      #ifdef BASKER_TIMER
      barrier_init_time += timer_init.seconds();
      std::cout << " ++ Basker order : sort(C) time : " << order_time << std::endl;
      #endif
    }
    #endif

    return 0;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::scotch_partition( BASKER_MATRIX &M )
  { 
    nd_flag = BASKER_TRUE;

    if(Options.symmetric == BASKER_TRUE)
    {
      //printf("Scotch Symmetric\n");
      part_scotch(M, part_tree);
    }
    else
    {
      //printf("Scotch Nonsymmetrix\n");
      BASKER_MATRIX MMT;
      AplusAT(M,MMT);
      //printMTX("AAT.mtx", MMT);
      part_scotch(MMT, part_tree);
      FREE(MMT);
    }

    nd_flag = BASKER_TRUE;
    //permute
    permute_row(M, part_tree.permtab);
    // new sfactor_copy2 replacement changes
    //permute_col(M, part_tree.permtab); //old, try the new below

    for (Int i = 0; i < M.nnz; ++i) {
      vals_order_scotch_array(i) = i; //init
    }
    permute_col_store_valperms(M, part_tree.permtab, vals_order_scotch_array); //NDE: Track movement of vals (lin_ind of row,col) here
    permute_inv(vals_order_ndbtfa_array, vals_order_scotch_array, M.nnz); //must permute the array holding the perms

    //May need to sort row_idx
    return 0; 
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
    INT_1DARRAY temp;
    MALLOC_INT_1DARRAY(temp,n);
    init_value(temp, n, (Int) 0);

    for(Int i = 0; i < n; ++i)
    {
      temp(p(i)) = vec(i);
    }
    for(Int i = 0; i < n; ++i)
    {
      vec(i) = temp(i);
    }

    FREE_INT_1DARRAY(temp);

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
   Int n
  )
  {
    //Permute
    for(Int i = 0; i < n; i++)
    {
      perm_comp_iworkspace_array(i) = vec(p(i));
    }
    //Copy back
    for(Int i = 0; i < n; i++)
    {
      vec(i) = perm_comp_iworkspace_array(i);
    }
    return BASKER_SUCCESS;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::permute_inv_with_workspace
  (
   INT_1DARRAY & vec,
   INT_1DARRAY & p,
   Int n
  )
  {
    //Permute
    for(Int i = 0; i < n; ++i)
    {
      perm_comp_iworkspace_array(p(i)) = vec(i);
    }
    //Copy back
    for(Int i = 0; i < n; ++i)
    {
      vec(i) = perm_comp_iworkspace_array(i);
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
    Int gn = ncols;
    // This macro is to test a minor improvement to this routine
    #define REDUCE_PERM_COMP_ARRAY 1

    for(Int i = 0; i < gn; i++) // Move this to a different init function
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
    //      permute(perm_inv_comp_array, gperm, gn);

    // NDE: This is moved out and initialized once - reused in 4 possible permutations
    for(Int i = BTF_A.ncol; i < gn; ++i)
    {
      perm_comp_iworkspace_array(i) = i;
    }

    MALLOC_INT_1DARRAY(numeric_row_iperm_array, gn);
    MALLOC_INT_1DARRAY(numeric_col_iperm_array, gn);

    MALLOC_INT_1DARRAY(symbolic_row_iperm_array, gn);
    MALLOC_INT_1DARRAY(symbolic_col_iperm_array, gn);
    for(Int i = 0; i < gn; ++i)
    {
      symbolic_col_iperm_array(i) = i;
      symbolic_row_iperm_array(i) = i;
    }

    // p5 inv
    if(amd_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE)
      {printf("AMD order in \n");}
      //printVec("amd.txt",order_csym_array, gn);
      for(Int i = 0; i < BTF_A.ncol; ++i)
      {
        perm_comp_iworkspace_array(i) = order_csym_array(i);
      }
      permute_with_workspace(perm_inv_comp_array, perm_comp_iworkspace_array, gn);
    }
    // p4 inv
    if(nd_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE)
      {printf("ND order in \n");}
      //printVec("nd.txt", part_tree.permtab, gn);
      for(Int i = 0; i < BTF_A.ncol; ++i)
      {
        perm_comp_iworkspace_array(i) = part_tree.permtab(i);
      }
      permute_with_workspace(perm_inv_comp_array, perm_comp_iworkspace_array, gn);
    }
    // p2, p3 inv
    if(btf_flag == BASKER_TRUE)
    {
      //printVec("btf_amd.txt", order_c_csym_array, gn);
      {
        // even if blk_matching is enabled, blk_amd & blk_mwm are still applied 
        //p3
        if(Options.verbose == BASKER_TRUE)
        { printf(" > blk_amd order in\n");}
        permute_with_workspace(perm_inv_comp_array, order_blk_amd_array, gn);

        if(Options.verbose == BASKER_TRUE)
        { printf(" > blk_mwm order in\n");}
        permute_with_workspace(perm_inv_comp_array, order_blk_mwm_array, gn);
      }
      //for (int i=0; i<(int)perm_inv_comp_array.extent(0); i++) printf( "%d %d\n",i,perm_inv_comp_array(i) );
      //printVec("btf.txt", order_btf_array, gn);
      //p2
      if(Options.verbose == BASKER_TRUE)
      { printf("btf order in\n"); }
      permute_with_workspace(perm_inv_comp_array, order_btf_array, gn);
      //for (int i=0; i<(int)perm_inv_comp_array.extent(0); i++) printf( "%d %d\n",i,perm_inv_comp_array(i) );

      // compose perm for symbolic phase
      permute_with_workspace(symbolic_row_iperm_array, order_btf_array, gn);
      permute_with_workspace(symbolic_col_iperm_array, order_btf_array, gn);
    }
    // p1 inv
    if(match_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE)
      {printf("match order in\n");}
      //printVec("match.txt", order_match_array, gn);
      permute_with_workspace(perm_inv_comp_array, order_match_array, gn);

      // compose perm for symbolic phase
      permute_with_workspace(symbolic_row_iperm_array, order_match_array, gn);
    }
    MALLOC_INT_1DARRAY(symbolic_row_perm_array, gn);
    MALLOC_INT_1DARRAY(symbolic_col_perm_array, gn);
    for(Int i = 0; i < gn; ++i)
    {
      symbolic_col_perm_array(i) = i;
      symbolic_row_perm_array(i) = i;
    }
    permute_inv(symbolic_row_perm_array, symbolic_row_iperm_array, gn);
    permute_inv(symbolic_col_perm_array, symbolic_col_iperm_array, gn);

    // ===============================================================================
    // perm_comp_array calculation
    // Note: don't need to inverse a row only perm
    // q1
    #if REDUCE_PERM_COMP_ARRAY
    for(Int i = 0; i < gn; i++) // NDE: Move this to a different init function
    {
      perm_comp_array(i) = perm_inv_comp_array(i);
    }
    if(match_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE)
      {printf("match order out\n");}
      //printVec("match.txt", order_match_array, gn);
      permute_inv_with_workspace(perm_comp_array, order_match_array, gn);
    }
    if(btf_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE)
      { printf(" > blk_mwm order out\n");}
      // invert btf first
      permute_inv_with_workspace(perm_comp_array, order_btf_array, gn);
      // now invert mwm
      permute_inv_with_workspace(perm_comp_array, order_blk_mwm_array, gn);
      // finally reapply btf
      permute_with_workspace(perm_comp_array, order_btf_array, gn);
    }
    #else
    if(amd_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE)
      {printf("AMD order out \n");}
      //printVec(order_csym_array, gn);
      for(Int i = 0; i < BTF_A.ncol; ++i)
      {
        perm_comp_iworkspace_array(i) = order_csym_array(i);
      }
      permute_with_workspace(perm_comp_array,perm_comp_iworkspace_array, gn);
    }
    // q2
    if(nd_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE)
      {printf("ND order out \n");}
      //printVec(part_tree.permtab, gn);
      for(Int i = 0; i < BTF_A.ncol; ++i)
      {
        perm_comp_iworkspace_array(i) = part_tree.permtab(i);
      }
      permute_with_workspace(perm_comp_array,perm_comp_iworkspace_array, gn);
    }
    // q3, q4
    if(btf_flag == BASKER_TRUE)
    {
      //printVec(order_btf_array, gn);
      // q3
      if(Options.verbose == BASKER_TRUE)
      {printf("blk_amd order out\n");}
      permute_with_workspace(perm_comp_array,order_blk_amd_array, gn);
      // q4
      if(Options.verbose == BASKER_TRUE)
      {printf("btf order out\n");}
      permute_with_workspace(perm_comp_array,order_btf_array, gn);
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
    Int nnz = M.nnz;

    INT_1DARRAY temp_p;
    MALLOC_INT_1DARRAY(temp_p, n+1);
    init_value(temp_p, n+1, (Int)0);
    INT_1DARRAY temp_i;
    MALLOC_INT_1DARRAY(temp_i, nnz);
    init_value(temp_i, nnz, (Int)0);
    ENTRY_1DARRAY temp_v;
    MALLOC_ENTRY_1DARRAY(temp_v, nnz);
    init_value(temp_v, nnz, (Entry)0.0);

    //Determine column ptr of output matrix
    for(Int j = 0; j < n; j++)
    {
      Int i = col (j);
      temp_p (i+1) = M.col_ptr (j+1) - M.col_ptr (j);
    }
    //Get ptrs from lengths
    temp_p (0) = 0;

    for(Int j = 0; j < n; j++)
    {
      temp_p (j+1) = temp_p (j+1) + temp_p (j);
    }
    //copy idxs
    for(Int ii = 0; ii < n; ii++)
    {
      Int ko = temp_p (col (ii) );
      for(Int k = M.col_ptr (ii); k < M.col_ptr (ii+1); k++)
      {
        temp_i (ko) = M.row_idx (k); // preserves order of indices in row_idx (they may not be ordered, but won't be reshuffled)
        temp_v (ko) = M.val (k);     // and corresponding order with vals

        vals_order_perm(k) = ko; //this works for first perm or individual; how to best compose subsequent perms? track separately and compose after?

        ko++;
      }
    }

    //copy back into A
    for(Int ii=0; ii < n+1; ii++)
    {
      M.col_ptr (ii) = temp_p (ii);
    }
    for(Int ii=0; ii < nnz; ii++)
    {
      M.row_idx (ii) = temp_i (ii);
      M.val (ii) = temp_v (ii);
    }
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
   INT_1DARRAY col
  )
  {
    if((M.ncol == 0)||(M.nnz == 0))
    { return 0; }

    Int n = M.ncol;
    Int nnz = M.nnz;

    INT_1DARRAY temp_p;
    MALLOC_INT_1DARRAY(temp_p, n+1);
    init_value(temp_p, n+1, (Int)0);
    INT_1DARRAY temp_i;
    MALLOC_INT_1DARRAY(temp_i, nnz);
    init_value(temp_i, nnz, (Int)0);
    ENTRY_1DARRAY temp_v;
    MALLOC_ENTRY_1DARRAY(temp_v, nnz);
    init_value(temp_v, nnz, (Entry)0.0);

    //Determine column ptr of output matrix
    for(Int j = 0; j < n; j++)
    {
      Int i = col (j);
      temp_p (i+1) = M.col_ptr (j+1) - M.col_ptr (j);
    }
    //Get ptrs from lengths
    temp_p (0) = 0;

    for(Int j = 0; j < n; j++)
    {
      temp_p (j+1) = temp_p (j+1) + temp_p (j);
    }
    //copy idxs
    for(Int ii = 0; ii < n; ii++)
    {
      Int ko = temp_p (col (ii));
      for(Int k = M.col_ptr (ii); k < M.col_ptr (ii+1); k++)
      {
        temp_i (ko) = M.row_idx (k); // preserves order of indices in row_idx (they may not be ordered, but won't be reshuffled)
        temp_v (ko) = M.val (k);     // and corresponding order with vals
        ko++;
      }
    }
    //copy back into A
    for(Int ii=0; ii < n+1; ii++)
    {
      M.col_ptr (ii) = temp_p (ii);
    }
    for(Int ii=0; ii < nnz; ii++)
    {
      M.row_idx (ii) = temp_i (ii);
      M.val (ii) = temp_v (ii);
    }
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

    // TODO: no need for the workspace?
    INT_1DARRAY temp_i;
    MALLOC_INT_1DARRAY(temp_i, nnz);
    init_value(temp_i, nnz, (Int)0);

    //permute
    for(Int k = 0; k < nnz; k++)
    {
      temp_i[k] = row[row_idx[k]];
    }
    //Copy back
    for(Int k = 0; k < nnz; k++)
    {
      row_idx[k] = temp_i[k];
    }

    FREE_INT_1DARRAY(temp_i);

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
#if 1
    permute_row(M.nnz, &(M.row_idx(0)), &(row(0)));
#else
    if(M.nnz == 0)
    { return 0; }

    Int nnz = M.nnz;
    INT_1DARRAY temp_i;
    MALLOC_INT_1DARRAY(temp_i, nnz);
    init_value(temp_i, nnz, (Int)0);

    //permute
    for(Int k = 0; k < nnz; k++)
    {
      temp_i[k] = row[M.row_idx[k]]; // where does row_id=M.row_idx[k] move to? moves to row[row_id] looks like inverse perm...
    }
    //Copy back
    for(Int k = 0; k < nnz; k++)
    {
      M.row_idx[k] = temp_i[k];
    }

    FREE_INT_1DARRAY(temp_i);
#endif
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
      #if defined(USE_WORKSPACE_FOR_SORT)
      for(Int i = 0; i < M.col_ptr[k+1]-start_row; i++)
      {
          iwork[i] = order_vals_perms(perm(i));
          dwork[i] = M.val[perm(i)];
          perm(i)  = M.row_idx[perm(i)];
      }
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
        std::cout << "    > Basker sort time     : " << order_time << " (nnz = " << M.col_ptr[k+1]-M.col_ptr[k] << ")" <<std::endl;
      }
      #endif
    }//end over all columns k
    //std::cout << " ];" << std::endl;

    return 0;
  }//end sort_matrix()


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
#endif //end ifndef basker_order_hpp
