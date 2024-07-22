// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_DEF_HPP
#define SHYLUBASKER_DEF_HPP

/*Basker Includes*/
#include "shylubasker_decl.hpp"
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_def.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"
#include "shylubasker_tree.hpp"
#include "shylubasker_sfactor.hpp"
#include "shylubasker_sfactor_inc.hpp"
#include "shylubasker_nfactor.hpp"
#include "shylubasker_nfactor_inc.hpp"
#include "shylubasker_solve_rhs.hpp"
#include "shylubasker_util.hpp"
#include "shylubasker_stats.hpp"
#include "shylubasker_order.hpp"

#include "shylubasker_solve_rhs_tr.hpp"

/*Kokkos Includes*/
#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#else
#include <omp.h>
#endif

/*System Includes*/
#include <iostream>

//#if defined(HAVE_AMESOS2_SUPERLUDIST) && !defined(BASKER_MC64)
//  #define BASKER_SUPERLUDIS_MC64
//#endif
//#define BASKER_TIMER 


namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Basker<Int, Entry, Exe_Space>::Basker()
  {   
    //Presetup flags
    matrix_flag       = BASKER_FALSE;
    order_flag        = BASKER_FALSE;
    tree_flag         = BASKER_FALSE;
    symb_flag         = BASKER_FALSE;
    factor_flag       = BASKER_FALSE;
    workspace_flag    = BASKER_FALSE;
    rhs_flag          = BASKER_FALSE;
    solve_flag        = BASKER_FALSE;
    nd_flag           = BASKER_FALSE;
    amd_flag          = BASKER_FALSE;
    same_pattern_flag = BASKER_FALSE;

    // for delayed pivot
    part_tree_saved   = BASKER_FALSE;

    //Default number of threads
    num_threads = 1;
    global_nnz  = 0;
    gn = 0;

    btf_total_work = 0;

    btf_tabs_offset = 0;
    btf_nblks = 0;
    btf_top_tabs_offset = 0;
    btf_top_nblks = 0;
  }//end Basker()
  

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Basker<Int ,Entry, Exe_Space>::~Basker()
  {}//end ~Basker()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::Finalize()
  {
    //finalize all matrices
    A.Finalize();
    At.Finalize(); //??? is At even used
    BTF_A.Finalize();
    BTF_C.Finalize();
    BTF_B.Finalize();
    BTF_D.Finalize();
    BTF_E.Finalize();
   
    //finalize array of 2d matrics
    FREE_MATRIX_2DARRAY(AVM, tree.nblks);
    FREE_MATRIX_2DARRAY(ALM, tree.nblks);
    
    FREE_MATRIX_2DARRAY(LL, tree.nblks);
    FREE_MATRIX_2DARRAY(LU, tree.nblks);
   
    FREE_INT_1DARRAY(LL_size);
    FREE_INT_1DARRAY(LU_size);
    
    //BTF structure
    FREE_INT_1DARRAY(btf_tabs);
    FREE_INT_1DARRAY(btf_blk_work);
    FREE_INT_1DARRAY(btf_blk_nnz);
    FREE_MATRIX_1DARRAY(LBTF);
    FREE_MATRIX_1DARRAY(UBTF);
       
    //Thread Array
    FREE_THREAD_1DARRAY(thread_array);
    basker_barrier.Finalize();
       
    //S (Check on this)
    FREE_INT_2DARRAY(S, tree.nblks);
    
    //Permuations
    FREE_INT_1DARRAY(gperm);
    FREE_INT_1DARRAY(gpermi);
    if(match_flag == BASKER_TRUE)
    {
      FREE_INT_1DARRAY(order_match_array);
      match_flag = BASKER_FALSE;
    }
    if(btf_flag == BASKER_TRUE)
    {
    //NDE: sfactor_copy2 replacement
      FREE_INT_1DARRAY(vals_order_btf_array);
      FREE_INT_1DARRAY(vals_order_blk_amd_array);
      FREE_INT_1DARRAY_PAIRS(vals_block_map_perm_pair); //this will map perm(val) indices to BTF_A, BTF_B, and BTF_C 

      FREE_INT_1DARRAY(vals_order_ndbtfd_array);
      FREE_INT_1DARRAY(vals_order_ndbtfe_array);
      FREE_INT_1DARRAY(vals_order_ndbtfa_array); //track nd perms; BTF_A must be declared here, else it does not exist
      FREE_INT_1DARRAY(vals_order_ndbtfb_array); //track nd perms; BTF_A must be declared here, else it does not exist
      FREE_INT_1DARRAY(vals_order_ndbtfc_array); //track nd perms; BTF_A must be declared here, else it does not exist

      FREE_INT_1DARRAY(inv_vals_order_ndbtfd_array);
      FREE_INT_1DARRAY(inv_vals_order_ndbtfe_array);
      FREE_INT_1DARRAY(inv_vals_order_ndbtfa_array);
      FREE_INT_1DARRAY(inv_vals_order_ndbtfb_array);
      FREE_INT_1DARRAY(inv_vals_order_ndbtfc_array);

      FREE_INT_1DARRAY(order_btf_array);
      FREE_INT_1DARRAY(order_blk_amd_array); // previous memory leak due to this array, likely
      btf_flag = BASKER_FALSE;
    }
    if(nd_flag == BASKER_TRUE)
    {
      FREE_INT_1DARRAY(vals_order_scotch_array);

      FREE_INT_1DARRAY(order_scotch_array);
      nd_flag = BASKER_FALSE;
    }
    if(amd_flag == BASKER_TRUE)
    {
      FREE_INT_1DARRAY(vals_order_csym_array);

      FREE_INT_1DARRAY(order_csym_array);
      amd_flag = BASKER_FALSE;
    }

    //NDE: sfactor_copy2 replacement
    FREE_INT_1DARRAY(vals_perm_composition);
    FREE_INT_1DARRAY(vals_crs_transpose);
    //FREE_ENTRY_1DARRAY(input_vals_unordered);

    //NDE: Free workspace and permutation arrays
    FREE_INT_1DARRAY(perm_comp_array);
    FREE_INT_1DARRAY(perm_inv_comp_array);
    FREE_INT_1DARRAY(perm_comp_iworkspace_array);
    FREE_ENTRY_1DARRAY(perm_comp_fworkspace_array);

    FREE_ENTRY_1DARRAY(x_view_ptr_copy);
    FREE_ENTRY_1DARRAY(y_view_ptr_copy);

    FREE_ENTRY_1DARRAY(x_view_ptr_scale);
    FREE_ENTRY_1DARRAY(y_view_ptr_scale);

    //Structures
    part_tree.Finalize();
    tree.Finalize();
    stree.Finalize();
    stats.Finalize();
  }//end Finalize()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::InitMatrix(string filename)
  { 
    //Note: jdb comeback to add trans option
    readMTX(filename, A);
    A.srow = 0;
    A.scol = 0;
    matrix_flag = true;
    return 0;
  }//end InitMatrix (file)

  template <class Int, class Entry, class Exe_Space >
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::InitMatrix
  (
   Int nrow, 
   Int ncol, 
   Int nnz, 
   Int *col_ptr,
   Int *row_idx, 
   Entry *val
  )
  {
    //Note: jdb comeback to add trans option
    A.init_matrix("Original Matrix", nrow, ncol, nnz, col_ptr, row_idx, val);
    A.scol = 0;
    A.srow = 0;
    sort_matrix(A);
    matrix_flag = true;
    return 0;
  }//end InitMatrix (int, int , int, int *, int *, entry *)


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Order(Int option)
  {
    //Option = 0, FAIL NATURAL WITHOUT BOX
    //Option = 1, BASKER Standard
    //Option = 2, BTF BASKER

    if(option == 1)
    {	
      default_order();
    }
    else if(option == 2)
    {
      btf_order();
    }
    else
    {
      printf("\n\n ERROR---No Order Selected \n\n");
      return -1;
    }

    basker_barrier.init(num_threads, 16, tree.nlvls);

    order_flag = true;
    return 0;
  }//end Order()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::InitOrder(Int option)
  {
    tree_flag = true;
    return 0;
  }//end InitOrder


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::InitOrder
  (
   Int *perm, 
   Int nblks, 
   Int parts, 
   Int *row_tabs, Int *col_tabs,
   Int *tree_tabs
  )
  {
    /*------------OLD
    init_tree(perm, nblks, parts, row_tabs, col_tabs, tree_tabs, 0);
    #ifdef BASKER_2DL
    matrix_to_views_2D(A);
    find_2D_convert(A);
    #else
    matrix_to_views(A,AV);
    #endif
    ----------*/

    user_order(perm,nblks,parts,row_tabs,col_tabs, tree_tabs);

    basker_barrier.init(num_threads, 16, tree.nlvls );

    //printf("done with init order\n");

    tree_flag = true;
    return 0;
  }//end InitOrder


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Symbolic( Int option )
  {
    #ifdef BASKER_TIMER 
    Kokkos::Timer timer;
    double time = 0.0;
    #endif

    //symmetric_sfactor();
    sfactor(); // NDE: This is the 'old' routine or alternative? When is this case used?

    if(option == 0)
    {}
    else if(option == 1)
    {}

    #ifdef BASKER_TIMER
    time = timer.seconds();
    stats.time_sfactor += time;
    std::cout << "Basker Symbolic total time: " << time << std::endl;
    #endif

    // NDE store matrix dims here
    sym_gn = A.ncol;
    sym_gm = A.nrow;
    MALLOC_ENTRY_1DARRAY(x_view_ptr_copy, sym_gn); //used in basker_solve_rhs - move alloc
    MALLOC_ENTRY_1DARRAY(y_view_ptr_copy, sym_gm);

    MALLOC_ENTRY_1DARRAY(x_view_ptr_scale, sym_gn); //used in basker_solve_rhs - move alloc
    MALLOC_ENTRY_1DARRAY(y_view_ptr_scale, sym_gm);

    MALLOC_INT_1DARRAY(perm_inv_comp_array , sym_gm); //y
    MALLOC_INT_1DARRAY(perm_comp_array, sym_gn); //x

    Int lwork = 3 * sym_gn;
    if (btf_tabs_offset != 0) {
      lwork = (BTF_A.nnz > lwork ? BTF_A.nnz : lwork);
      lwork = (BTF_B.nnz > lwork ? BTF_B.nnz : lwork);
      lwork = (BTF_C.nnz > lwork ? BTF_C.nnz : lwork);
    }
    MALLOC_INT_1DARRAY(perm_comp_iworkspace_array, lwork); 
    MALLOC_ENTRY_1DARRAY(perm_comp_fworkspace_array, sym_gn);
    permute_composition_for_solve(sym_gn);

    return 0;
  }//end Symbolic


  //This is the interface for Amesos2
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Symbolic
  (
   Int nrow,
   Int ncol,
   Int nnz, 
   Int *col_ptr, 
   Int *row_idx, 
   Entry *val,
   bool crs_transpose_needed_
  )
  {
    #ifdef BASKER_TIMER 
    std::ios::fmtflags old_settings = cout.flags();
    int old_precision = std::cout.precision();
    std::cout.setf(ios::fixed, ios::floatfield);
    std::cout.precision(8);
    double time = 0.0;
    Kokkos::Timer timer;
    #endif

    if(Options.verbose == BASKER_TRUE)
    {
      std::cout << "\n == Basker Symbolic ==" << std::endl;
      std::cout << "Matrix dims: " << nrow << " x " << ncol << ", nnz = " << nnz << std::endl;
    }

    //Init Matrix A.
    if(matrix_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE) {
        printf("BASKER: YOU CANNOT RERUN SYMBOLIC\n");
      }
      return BASKER_ERROR;
    }
    else
    {
      #ifdef BASKER_TIMER
      double init_time = 0.0;
      Kokkos::Timer timer_init;
      #endif

      // NDE: New for Amesos2 CRS mods
      MALLOC_INT_1DARRAY(vals_crs_transpose,nnz);
      for( Int i = 0; i < nnz; ++i ){
        vals_crs_transpose(i) = i; // init vals permutation due to transpose
      }

      A.init_matrix("Original Matrix", nrow, ncol, nnz, col_ptr, row_idx, val);
      A.scol = 0;
      A.srow = 0;

      if(Options.transpose == BASKER_FALSE)
      {
        // NDE: New for Amesos2 CRS mods
        if ( crs_transpose_needed_ ) {
          matrix_transpose(0, nrow, 0, ncol, nnz, col_ptr, row_idx, val, A, vals_crs_transpose);
        }
      }
      else
      {
        if ( !crs_transpose_needed_ ) {
          matrix_transpose(0, nrow, 0, ncol, nnz, col_ptr, row_idx, val, A, vals_crs_transpose);
        }
      }
      sort_matrix(A);
      this->crs_transpose_needed = crs_transpose_needed_;

      if(Options.verbose == BASKER_TRUE)
      {
        printf("Basker Matrix Loaded \n");
      }

      if(Options.verbose_matrix_out == BASKER_TRUE)
      {
        printMTX("A_Symbolic.mtx", A);
      }

      matrix_flag = BASKER_TRUE;

      #ifdef BASKER_TIMER
      init_time += timer_init.seconds();
      std::cout << std::endl << "Basker Symbolic matrix init time: " << init_time << std::endl;
      #endif
    }
    /*FILE *fp = fopen("symbolic.dat","w");
    printf( " start symbolic with\n original A = [\n" );
    for(Int j = 0; j < A.ncol; j++) {
      for(Int k = col_ptr[j]; k < col_ptr[j+1]; k++) {
        fprintf(fp, " %d %d %e\n", row_idx[k],j,val[k]);
      }
    }
    fclose(fp);
    printf("];\n");*/

    /*printf( " start symbolic with\n A = [\n" );
    for(Int j = 0; j < A.ncol; j++) {
      for(Int k = A.col_ptr[j]; k < A.col_ptr[j+1]; k++) {
        printf( " %d %d %e\n", A.row_idx[k],j,A.val[k]);
      }
    }
    printf("];\n");*/

    //Init Ordering
    //Always will do btf_ordering
    //This should also call create tree
    if(order_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE) {
        printf("BASKER: YOU CANNOT RERUN ORDER\n");
      }
      return BASKER_ERROR;
    }
    else
    {
      #ifdef BASKER_TIMER
      double order_time = 0.0;
      double barrier_init_time = 0.0;
      Kokkos::Timer timer_order;
      #endif
      //-------------------------------------------------
      //Find BTF ordering
      if(btf_order2() != BASKER_SUCCESS)
      {
        return BASKER_ERROR;
      }

      if(Options.verbose == BASKER_TRUE)
      {
        printf("Basker Ordering Found \n");
      }

      /*if((Options.btf == BASKER_TRUE) && (btf_tabs_offset != 0))
      {
        #ifdef BASKER_TIMER
        Kokkos::Timer timer_init;
        #endif

        basker_barrier.init(num_threads, 16, tree.nlvls );

        #ifdef BASKER_TIMER
        barrier_init_time += timer_init.seconds();
        #endif
      }*/
      order_flag = BASKER_TRUE;

      if(Options.verbose == BASKER_TRUE)
      {
        printf("Basker P2P Thread Barriers Init\n");
      }

      #ifdef BASKER_TIMER
      order_time += timer_order.seconds();
      std::cout << "Basker Symbolic order arrays time: "    << order_time << std::endl;
      std::cout << " > Basker Symbolic init barrier time: " << barrier_init_time << std::endl;
      #endif
    }

    if(symb_flag == BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE) {
        printf("BASKER: YOU CANNOT RERUN SFACTOR\n");
      }
      return BASKER_ERROR;
    }
    else
    {
      #ifdef BASKER_TIMER
      double sfactor_time = 0.0;
      Kokkos::Timer timer_sfactor;
      #endif
      //-------------------------------------------------
      //Do symbolic factorization
      if(Options.incomplete == BASKER_FALSE)
      {
        sfactor();
      }
      else
      {
        sfactor_inc();
      }
      #ifdef BASKER_TIMER
      sfactor_time += timer_sfactor.seconds();
      std::cout << "Basker Symbolic sfactor time: " << sfactor_time << std::endl;
      #endif

      if(Options.verbose == BASKER_TRUE)
      {
        printf("Basker Nonzero Counts Found \n");
      }
      symb_flag = BASKER_TRUE;
    }


    if(Options.verbose == BASKER_TRUE)
    {
      printf(" == Basker Symbolic Done ==\n\n");
    }

    #ifdef BASKER_TIMER
    time = timer.seconds();
    stats.time_sfactor += time;
    std::cout << "Basker Symbolic total time: " << time
              << std::endl << std::endl;
    std::cout.precision(old_precision);
    std::cout.flags(old_settings);
    #endif

    // NDE store matrix dims here for comparison in Factor
    sym_gn = A.ncol;
    sym_gm = A.nrow;
    MALLOC_ENTRY_1DARRAY(x_view_ptr_copy, sym_gn); //used in basker_solve_rhs
    MALLOC_ENTRY_1DARRAY(y_view_ptr_copy, sym_gm);

    MALLOC_ENTRY_1DARRAY(x_view_ptr_scale, sym_gn); //used in basker_solve_rhs
    MALLOC_ENTRY_1DARRAY(y_view_ptr_scale, sym_gm);

    MALLOC_INT_1DARRAY(perm_inv_comp_array, sym_gm); //y
    MALLOC_INT_1DARRAY(perm_comp_array, sym_gn); //x

    Int lwork = 3 * sym_gn;
    if (btf_tabs_offset != 0) {
      lwork = (BTF_A.nnz > lwork ? BTF_A.nnz : lwork);
      lwork = (BTF_B.nnz > lwork ? BTF_B.nnz : lwork);
      lwork = (BTF_C.nnz > lwork ? BTF_C.nnz : lwork);
    }
    MALLOC_INT_1DARRAY(perm_comp_iworkspace_array, lwork);
    MALLOC_ENTRY_1DARRAY(perm_comp_fworkspace_array, sym_gn);
    permute_composition_for_solve(sym_gn);

    /*printf( " end symbolic with\n original A = [\n" );
    for(Int j = 0; j < ncol; j++) {
      for(Int k = col_ptr[j]; k < col_ptr[j+1]; k++) {
        printf( " %d %d %e\n", row_idx[k],j,val[k]);
      }
    }
    printf("];\n");*/

    return 0;
  } //end Symbolic()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Factor(Int option)
  {
    #ifdef BASKER_TIMER
    Kokkos::Timer timer;
    double time = 0.0;
    #endif

    if(symb_flag != BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE) {
        printf("BASKER: YOU NEED TO RUN SYMBOLIC BEFORE FACTOR\n");
      }
      return BASKER_ERROR;
    }

    //Reset error codes
    reset_error();

    //Do numerical factorization
    factor_notoken(option);

    #ifdef BASKER_TIMER
    time += timer.seconds();
    stats.time_nfactor += time;
    std::cout << "Basker factor_notoken time: " << time << std::endl;
    std::cout << "Basker Factor total   time: " << time << std::endl;
    #endif

    factor_flag = BASKER_TRUE;

    return 0;
  }//end Factor()


  //This is the interface for Amesos2
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::Factor
  (
   Int nrow, 
   Int ncol,
   Int nnz, 
   Int *col_ptr, 
   Int *row_idx, 
   Entry *val
  ) 
  {
    //Reset error codes
    int err = 0; //init for return value from sfactor_copy2, factor_notoken etc.

    #ifdef BASKER_TIMER
    std::ios::fmtflags old_settings = cout.flags();
    int old_precision = std::cout.precision();
    std::cout.setf(ios::fixed, ios::floatfield);
    std::cout.precision(8);
    double time = 0.0;
    Kokkos::Timer timer;
    #endif
    if(Options.verbose == BASKER_TRUE)
    {
      printf("\n == Basker Factor ==\n\n");
    }

    if(symb_flag != BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE) {
        printf("BASKER: YOU NEED TO RUN SYMBOLIC BEFORE FACTOR\n");
      }
      return BASKER_ERROR;
    }

    //sfactor_copy2 stuff
    // This part is stored in case a matrix_transpose will be needed (if input is passed in as CRS)
    // WARNING! This may require a sort step first if CRS matrix input - make sure Symbolic matches input pattern ie CCS vs CRS
    //BASKER_MATRIX M;
    //M.init_matrix("Original Matrix", nrow, ncol, nnz, col_ptr, row_idx, val);
    //sort_matrix(M);
    //MALLOC_ENTRY_1DARRAY(input_vals_unordered,nnz);
    //for( Int i = 0; i < nnz; ++i ) { //sfactor_copy2 replacement: setup for testing success of perms
    //  input_vals_unordered(i) = M.val(i);
    //}

    // new vars needed
    //inv_vals_order_ndbtfa_array, inv_vals_order_ndbtfb_array, inv_vals_order_ndbtfc_array
    //vals_block_map_perm_pair, vals_perm_composition


    // Summary:
    // When Symbolic is called, the compressed matrix pointer data is copied into Basker. 
    // If this is done through Amesos2 with a single process, the CRS format is not converted to CCS; 
    // that occurs by taking the transpose of the input CRS structure matrix, then performing symbolic factorization.
    // The results of the Symbolic step are stored as permutation vectors. 
    // When Factor is called, the input pointers will match those that were passed to Symbolic (i.e. if done through
    // Amesos2 on a single process, they will be in CRS format). The transpose operation, if applied to the pointers, 
    // is encoded within the perm_composition arrays. 
    // Factor starts by copying the input pointer data into the local pointers. This is done by applying the permutations
    // stored during Symbolic phase to the input pointers and storing the results in the local pointers. 
    // This permute+copy operation stores the values in the A.val array as well as the blocks BTF_*.val; colptr and rowidx 
    // are assumed to be unchtetBasis.getValues(dbasisAtLattice, lattice, OPERATOR_D3);anged from the time Symbolic was called, and so do not need to be copied or recomputed.
    // After this, sfactor_copy2 is called with copies the results from BTF_A to the 2d ND blocks
    //
    // This assumes the ColPtr and RowIdx values SHOULD be identical if reusing the Symbolic structure, but if there
    // is a change, for example reusing the Symbolic structure but a change in the values of the matrix entries and if 
    // some diagonal entry ends up becoming a zero, then Symbolic should be rerun
    //
/*
    // Look for zeros on the diagonal - an Option / parameter list option should be used to enable this
    // col_ptr has ncol+1 entries - final entry = nnz
    // row_idx and val have nnz entries
    // Scan through the matrix for 0's on the diagonal
    if ( nrow == ncol ) // square matrices, CRS vs CCS is irrelevant (i.e. transpose does not need to be applied first) - but needs to be known for Symbolic...
    {
      std::cout << "  Check diagonal entries..." << std::endl;
      Int total_diag_found = 0;
      for ( Int colid = 0; colid < ncol; ++colid ) {
        Int rowids_start = col_ptr[colid];
        Int rowids_end = col_ptr[colid+1];
        Int num_row_entries = rowids_end - rowids_start;
        Int counter = 0;
        Int counter_diag_found = 0;
        for ( Int j = rowids_start; j < rowids_end; ++j ) {
          Int rowid = row_idx[j];
          if ( rowid != colid ) {
            ++counter; // track how many are non-diag entries
          }
          if ( rowid == colid ) {
            ++counter_diag_found;
            ++total_diag_found;
            //check non-zero - this assumes there IS an entry for every diagonal
            if ( val[j] == 0 )
              std::cout << "ZERO ON DIAG!!!" << std::endl;
          }
        }
        if ( counter == num_row_entries ) { // if this happens, a diag was never found i.e. assumed 0
          std::cout << " colid = " << colid << "  counter = " << counter << "  num_row_entries = " << num_row_entries << std::endl;
          std::cout << " ZERO DIAG!" << std::endl;
        }
//        std::cout << "colid: " << colid << "   diag_found: " << counter_diag_found << std::endl;
      }
      if ( total_diag_found != ncol ) {
        std::cout << "MISSING DIAG ENTRY!  Found: " << total_diag_found << "  expected: " << ncol << std::endl;
      }
    }
*/

    /*{
      FILE *fp;
      if(Options.transpose) {
        fp = fopen("numeric_t.dat","w");
      } else {
        fp = fopen("numeric.dat","w");
      }
      printf( " start numeric with\n original A = [\n" );
      for(Int j = 0; j < ncol; j++) {
        for(Int k = col_ptr[j]; k < col_ptr[j+1]; k++) {
          //fprintf(fp, " %d %d %e\n", row_idx[k],j,val[k]);
          printf(" %d %d %e\n", row_idx[k],j,val[k]);
        }
      }
      fclose(fp);
      printf("];\n");
    }*/
    //BTF_A.print_matrix("AA.dat");
    //BTF_B.print_matrix("BB.dat");
    //BTF_C.print_matrix("CC.dat");
    //BTF_D.print_matrix("DD.dat");
    //BTF_E.print_matrix("EE.dat");

    using range_type = Kokkos::pair<int, int>;
    if (Options.static_delayed_pivot != 0 || Options.blk_matching != 0) {
        Kokkos::Timer resetperm_timer;

        //Int b_first = btf_tabs_offset;
        Int nfirst = btf_tabs(btf_tabs_offset);
        // to revert BLK_AMD ordering
        for (Int k = 0; k < (Int)ncol; k++) order_blk_amd_inv(k) = k;
        permute_inv(order_blk_amd_inv, order_blk_amd_array, ncol); 

        // to revert BLK_MWM ordering
        for (Int k = 0; k < (Int)ncol; k++) order_blk_mwm_inv(k) = k;
        permute_inv(order_blk_mwm_inv, order_blk_mwm_array, ncol);

        // --------------------------------------------
        // reset the large A block
        if (btf_tabs_offset != 0) {
          Int a_first = btf_top_tabs_offset;
          Int a_nfirst = btf_tabs(a_first);
          Int a_nlast  = a_nfirst + BTF_A.ncol;

          if (Options.blk_matching != 0 && Options.static_delayed_pivot == 0) {
            // revert the sort matrix for ND
            permute_with_workspace(BTF_A.row_idx, inv_vals_order_ndbtfa_array, BTF_A.nnz);
            permute_with_workspace(BTF_B.row_idx, inv_vals_order_ndbtfb_array, BTF_B.nnz);
            permute_with_workspace(BTF_C.row_idx, inv_vals_order_ndbtfc_array, BTF_C.nnz);

            // --------------------------------------------
            // reinitialize ordering vector to track matrix sort
            for (Int i = 0; i < BTF_A.nnz; ++i) {
              vals_order_ndbtfa_array(i) = i;
              inv_vals_order_ndbtfa_array(i) = i;
            }
            for (Int i = 0; i < BTF_B.nnz; ++i) {
              vals_order_ndbtfb_array(i) = i;
              inv_vals_order_ndbtfb_array(i) = i;
            }
            for (Int i = 0; i < BTF_C.nnz; ++i) {
              vals_order_ndbtfc_array(i) = i;
              inv_vals_order_ndbtfc_array(i) = i;
            }
            #if 1
            // revert camd & nd ordering
            // > revert camd
            if (order_csym_array.extent(0) > 0) {
              for (Int k = 0; k < (Int)BTF_A.ncol; k++) order_csym_inv(k) = k;
              permute_inv(order_csym_inv, order_csym_array, BTF_A.ncol);

              // Revert camd ordering to A
              permute_row(BTF_A, order_csym_inv);
              permute_col(BTF_A, order_csym_inv);
              if (BTF_E.ncol > 0) {
                // Revert CAMD perm to cols of E
                permute_col(BTF_E, order_csym_inv);
              }
              if (btf_tabs_offset < btf_nblks) {
                // Revert ND perm to rows of B
                permute_row(BTF_B, order_csym_inv);
              }
            }

            // > revert nd
            if (part_tree.permtab.extent(0) > 0) {
              for (Int k = 0; k < (Int)BTF_A.ncol; k++) order_nd_inv(k) = k;
              permute_inv(order_nd_inv, part_tree.permtab, BTF_A.ncol); 
              //for (Int k = 0; k < (Int)BTF_A.ncol; k++) printf( " > nd(%d)=%d\n",k,part_tree.permtab(k) );

              // Revert ND ordering to A
              permute_row(BTF_A, order_nd_inv);
              permute_col(BTF_A, order_nd_inv);
              if (BTF_E.ncol > 0) {
                // Revert ND perm to cols of E
                permute_col(BTF_E, order_nd_inv);
              }
              if (btf_tabs_offset < btf_nblks) {
                // Revert ND perm to rows of B
                permute_row(BTF_B, order_nd_inv);
              }
            }
            #endif
          }
          if (Options.static_delayed_pivot != 0 && part_tree_saved == BASKER_TRUE) {
            // revert sort after ND fix
            permute_inv_with_workspace(BTF_A.row_idx, vals_order_ndbtfa_array, BTF_A.nnz);

            // reload original ND
            for (Int j = 0; j < (Int)BTF_A.ncol; j++) {
              part_tree.permtab[j] = part_tree_orig.permtab[j];
              part_tree.ipermtab[j] = part_tree_orig.ipermtab[j];
            }
            for (Int k = 0; k <= part_tree.nblks; k++) {
              part_tree.row_tabs(k) = part_tree_orig.row_tabs(k);
              part_tree.col_tabs(k) = part_tree_orig.col_tabs(k);

              tree.row_tabs(k) = part_tree_orig.row_tabs(k);
              tree.col_tabs(k) = part_tree_orig.col_tabs(k);
            }
          }

          // revert AMD perm to rows of A
          auto order_nd_amd = Kokkos::subview(order_blk_amd_inv,
                                              range_type(a_nfirst, a_nlast));
          for (Int i = 0; i < (Int)BTF_A.ncol; i++) {
            order_nd_amd(i) -= a_nfirst;
          }
          permute_row(BTF_A, order_nd_amd);
          permute_col(BTF_A, order_nd_amd);
          if (BTF_E.ncol > 0) {
            // Revert ND perm to cols of E
            permute_col(BTF_E, order_nd_amd);
          }
          if (btf_tabs_offset < btf_nblks) {
            // revert AMD perm to rows and cols
            permute_row(BTF_B, order_nd_amd);
          }
          for (Int i = 0; i < (Int)BTF_A.ncol; i++) {
            order_nd_amd(i) += a_nfirst;
          }

          if (Options.static_delayed_pivot != 0 || Options.blk_matching != 0) {
            // revert MWM perm to rows of A
            auto order_nd_mwm = Kokkos::subview(order_blk_mwm_inv, 
                                                range_type(a_nfirst, a_nlast));
            for (Int i = 0; i < (Int)BTF_A.ncol; i++) {
              order_nd_mwm(i) -= a_nfirst;
            }
            permute_row(BTF_A, order_nd_mwm);
            if (btf_tabs_offset < btf_nblks) {
              permute_row(BTF_B, order_nd_mwm);
            }
            for (Int i = 0; i < (Int)BTF_A.ncol; i++) {
              order_nd_mwm(i) += a_nfirst;
            }
          }
        }

        // --------------------------------------------
        // reset the small C blocks
        if (btf_nblks > btf_tabs_offset) {
          // revert BLK_AMD ordering
          auto order_blk_amd_c = Kokkos::subview(order_blk_amd_inv,
                                                 range_type(nfirst, ncol));
          for (Int i = 0; i < (Int)BTF_C.ncol; i++) {
            order_blk_amd_c(i) -= nfirst;
          }
          permute_row(BTF_C, order_blk_amd_c);
          permute_col(BTF_C, order_blk_amd_c);
          if (BTF_E.ncol > 0) {
            // Apply AMD perm to cols
            permute_col(BTF_E, order_blk_amd_c, BTF_A.ncol);
          }
          if (BTF_B.ncol > 0) {
            // Apply AMD perm to cols
            permute_col(BTF_B, order_blk_amd_c);
          }
          for (Int i = 0; i < (Int)BTF_C.ncol; i++) {
            order_blk_amd_c(i) += nfirst;
          }

          if (Options.blk_matching != 0) {
            // revert BLK_MWM ordering
            auto order_blk_mwm_c = Kokkos::subview(order_blk_mwm_inv, 
                                                   range_type (nfirst, ncol));
            for (Int i = 0; i < (Int)BTF_C.ncol; i++) {
              order_blk_mwm_c(i) -= nfirst;
            }
            permute_row(BTF_C, order_blk_mwm_c);
            for (Int i = 0; i < (Int)BTF_C.ncol; i++) {
              order_blk_mwm_c(i) += nfirst;
            }
          }
        }

        #define SHYLU_BASKER_AMD_ON_D
        #ifdef SHYLU_BASKER_AMD_ON_D
        // --------------------------------------------
        // reset the small D blocks
        if (btf_top_tabs_offset  > 0) {
          Int d_last = btf_top_tabs_offset;
          Int ncol = btf_tabs(d_last);

          // revert BLK_AMD ordering
          auto order_blk_amd_d = Kokkos::subview(order_blk_amd_inv,
                                                 range_type(0, ncol));
          permute_row(BTF_D, order_blk_amd_d);
          permute_col(BTF_D, order_blk_amd_d);
          if (BTF_E.ncol > 0) {
            // Apply AMD perm to cols
            permute_row(BTF_E, order_blk_amd_d);
          }

          if (Options.blk_matching != 0) {
            // revert BLK_MWM ordering
            auto order_blk_mwm_d = Kokkos::subview(order_blk_mwm_inv, 
                                                   range_type (0, ncol));
            permute_row(BTF_D, order_blk_mwm_d);
            if (BTF_E.ncol > 0) {
              // Apply MWM perm to cols
              permute_row(BTF_E, order_blk_mwm_d);
            }
          }
        }
        #endif
        if(Options.verbose == BASKER_TRUE) {
          std::cout<< "Basker Factor: Time to revert all the numerical perm: " << resetperm_timer.seconds() << std::endl;
        }
    } // end of if(Options.blk_matching != 0)

    Kokkos::Timer copyperm_timer;
    //printf( " A.nnz= %d vs (%d, %d) nblks=%d, btfa_nnz=%d, btfb_nnz=%d, btfc_nnz=%d\n",(int)nnz, (int)A.nnz,(int)A.val.extent(0),
    //        btf_nblks,btfa_nnz,btfb_nnz,btfc_nnz );
    if (btf_nblks == 0) {
      std::cout << "Basker Factor error: Case for btf_nbkls = 0 is not implemented" << std::endl;
        //A.val(i) = val[ i ]; // may need to apply matching or nd order permutation...
      return BASKER_ERROR;
    } else {
    #ifdef KOKKOS_ENABLE_OPENMP
    #pragma omp parallel for
    #endif
      for( Int i = 0; i < nnz; ++i ) {
        A.val(i) = val[ vals_perm_composition(i) ];
        if ( btfd_nnz != 0 ) {
          if ( vals_block_map_perm_pair(i).first == -1 ) { //in BTF_D
            //Int k = inv_vals_order_ndbtfd_array( vals_block_map_perm_pair(i).second );
            //printf( " BTF_D.val(inv(%d -> %d) = %d) = val[perm(%d) = %d] = %e (%d)\n",i,vals_block_map_perm_pair(i).second, k,
            //        i,vals_perm_composition(i), val[ vals_perm_composition(i) ], BTF_D.row_idx(k));
            BTF_D.val( inv_vals_order_ndbtfd_array( vals_block_map_perm_pair(i).second ) ) = A.val(i);
          }
        }
        if ( btfe_nnz != 0 ) {
          if ( vals_block_map_perm_pair(i).first == -2 ) { //in BTF_E
            //Int k = inv_vals_order_ndbtfe_array( vals_block_map_perm_pair(i).second );
            //printf( " BTF_E.val(inv(%d -> %d) = %d) = val[perm(%d) = %d] = %e (%d)\n",i,vals_block_map_perm_pair(i).second, k,
            //        i,vals_perm_composition(i), val[ vals_perm_composition(i) ], BTF_E.row_idx(k));
            BTF_E.val( inv_vals_order_ndbtfe_array( vals_block_map_perm_pair(i).second ) ) = A.val(i);
          }
        }

        //printf( " A.val(%d) = %e (first = %d, btf_nnz=%d,%d,%d)\n",i,val[ vals_perm_composition(i) ], vals_block_map_perm_pair(i).first,btfa_nnz,btfb_nnz,btfc_nnz );
        if ( btfa_nnz != 0 ) { //is this unnecessary? yes, but shouldn't the first label account for this anyway?
          if ( vals_block_map_perm_pair(i).first == 0 ) { //in BTF_A
            //Int k = inv_vals_order_ndbtfa_array( vals_block_map_perm_pair(i).second );
            //printf( " BTF_A.val(inv(%d -> %d) = %d) = val[perm(%d) = %d] = %e (%d)\n",i,vals_block_map_perm_pair(i).second, k,
            //        i,vals_perm_composition(i), val[ vals_perm_composition(i) ], BTF_A.row_idx(k));
            BTF_A.val( inv_vals_order_ndbtfa_array( vals_block_map_perm_pair(i).second ) ) = A.val(i);
          }
        }
        if ( btfb_nnz != 0 ) {
          if ( vals_block_map_perm_pair(i).first == 1 ) { //in BTF_B
            //printf( " BTF_B.val(%d) = %e\n",inv_vals_order_ndbtfb_array( vals_block_map_perm_pair(i).second ),val[ vals_perm_composition(i) ] );
            BTF_B.val( inv_vals_order_ndbtfb_array( vals_block_map_perm_pair(i).second ) ) = A.val(i);
          }
        }
        // if ( BTF_C.nnz != 0 ) // this causes compiler error, and nnz blocks values are different with this command with small blocks matrix (not nd) - why?
        if ( btfc_nnz != 0 ) {
          if ( vals_block_map_perm_pair(i).first == 2 ) { //in BTF_C
            //printf( " BTF_C.val(%d) = %e\n",inv_vals_order_ndbtfc_array( vals_block_map_perm_pair(i).second ), val[ vals_perm_composition(i) ] );
            BTF_C.val( inv_vals_order_ndbtfc_array( vals_block_map_perm_pair(i).second ) ) = A.val(i);
          }
        }
      } //end for
    } //end if
    if(Options.verbose == BASKER_TRUE) {
      std::cout << "Basker Factor: Time to permute and copy from input vals to new vals and blocks: "
                << copyperm_timer.seconds() << std::endl;
    }
    //end sfactor_copy2 replacement stuff


    /*printf(" K = [\n" );
    for(Int j = 0; j < A.ncol; j++) {
      for(Int k = col_ptr[j]; k < col_ptr[j+1]; k++) {
        printf("%d %d %e\n", row_idx[k], j, val[k]);
      }
    }
    printf("];\n");*/

    //BTF_A.print_matrix("A_.dat");
    //BTF_B.print_matrix("B_.dat");
    //BTF_C.print_matrix("C_.dat");
    //BTF_D.print_matrix("D_.dat");
    //BTF_E.print_matrix("E_.dat");


    if ((Options.replace_zero_pivot || Options.replace_tiny_pivot) && BTF_A.nnz > 0) {
      // to compute one-norm of global A and BTF_A
      Kokkos::Timer normA_timer;
      using STS = Teuchos::ScalarTraits<Entry>;
      using Mag = typename STS::magnitudeType;
      const Entry zero (0.0);
      A.anorm = abs(zero);
      for (Int j = 0; j < (Int)A.ncol; j++) {
        Mag anorm_j = abs(zero);
        for (Int k = A.col_ptr(j); k < A.col_ptr(j+1); k++) {
          anorm_j += abs(A.val(k));
        }
        if (anorm_j > A.anorm) {
          A.anorm = anorm_j;
        }
      }
      A.gnorm = A.anorm;

      for (Int j = 0; j < (Int)BTF_A.ncol; j++) {
        Mag anorm_j = abs(zero);
        for (Int k = BTF_A.col_ptr(j); k < BTF_A.col_ptr(j+1); k++) {
          anorm_j += abs(BTF_A.val(k));
        }
        if (anorm_j > BTF_A.anorm) {
          BTF_A.anorm = anorm_j;
        }
      }
      BTF_A.gnorm = A.anorm;
      if(Options.verbose == BASKER_TRUE) {
         cout<< " Basker Factor: Time to compute " 
             << " norm(A) = "     << BTF_A.gnorm << " with n = " << A.ncol << ", and "
             << " norm(BTF_A) = " << BTF_A.anorm << " with n = " << BTF_A.ncol
             << " : " << normA_timer.seconds() << std::endl;
      }
    }

    // ==================================================================================================
    // compute blk_mwm & blk_amd
    if(Options.verbose == BASKER_TRUE) {
      std::cout << " --- Factor::blk_matching = " << Options.blk_matching << " ---" << std::endl;
      std::cout << "     btf_tabs_offset = " << btf_tabs_offset << " btf_nblks = " << btf_nblks << std::endl;
    }
    if (Options.static_delayed_pivot != 0 || Options.blk_matching != 0) {
      // reinitialize the whole perm and scale (for now)
      Entry one (1.0);
      for (Int i = 0; i < (Int)A.nrow; i++) {
        order_blk_mwm_array(i) = i;
        order_blk_amd_array(i) = i;
        scale_row_array(i) = abs(one);
        scale_col_array(i) = abs(one);

        numeric_row_iperm_array(i) = i;
        numeric_col_iperm_array(i) = i;
      }

      #ifdef SHYLU_BASKER_AMD_ON_D
      if (btf_top_tabs_offset > 0) {
        Kokkos::Timer mwm_amd_perm_timer;
        Int d_last = btf_top_tabs_offset;
        Int ncol = btf_tabs(d_last);

        auto order_blk_mwm_d = Kokkos::subview(order_blk_mwm_array,
                                               range_type(0, ncol));
        auto order_blk_amd_d = Kokkos::subview(order_blk_amd_array,
                                               range_type(0, ncol));
        if(Options.verbose == BASKER_TRUE) {
          std::cout << " calling MWM on D" << std::endl;
          std::cout << " btf_blk_mwm: btf_top_tabs(" << btf_top_tabs_offset << ")=" << ncol << std::endl;
        }

        // ----------------------------------------------------------------------------------------------
        // recompute MWM and AMD on each block of D
        INT_1DARRAY blk_nnz;
        INT_1DARRAY blk_work;
        btf_blk_mwm_amd(0, d_last, BTF_D,
                        order_blk_mwm_d, order_blk_amd_d,
                        blk_nnz, blk_work);
        #if 0 //debug
        printf( " >> debug: set order_blk_mwm/amd on D to identity <<\n" );
        for (Int i = 0; i < (Int)BTF_D.nrow; i++) {
          order_blk_mwm_d(i) = i;
          order_blk_amd_d(i) = i;
        }
        #endif

        // Apply MWM perm to rows
        permute_row(BTF_D, order_blk_mwm_d);

        // Apply AMD perm to rows and cols
        // NOTE: no need to udpate vals_order_blk_amd_array (it will read in A without permutations, and compute perm here)
        permute_row(BTF_D, order_blk_amd_d);
        permute_col(BTF_D, order_blk_amd_d);
        if (BTF_E.ncol > 0) {
          // Apply AMD perm to cols
          permute_row(BTF_E, order_blk_amd_d);
        }

        if(Options.verbose == BASKER_TRUE) {
          std::cout << std::endl << std::endl;
          std::cout << "Basker Factor: Time to compute and apply MWM+AMD on diagonal blocks: " << mwm_amd_perm_timer.seconds() << std::endl;
        }
      }
      #endif

      if (btf_tabs_offset != 0) {
        Kokkos::Timer nd_perm_timer;

        // ----------------------------------------------------------------------------------------------
        // compute MWM on a big block A
        if(Options.verbose == BASKER_TRUE) {
          std::cout << " calling MWM on A(n=" << BTF_A.ncol << ", nnz=" << BTF_A.nnz 
                    << ", btf_tabs_offset=" << btf_tabs_offset << ", " << btf_tabs(0) << " : " << btf_tabs(btf_tabs_offset)-1
                    << " )" << std::endl;
        }

        // ----------------------------------------------------------------------------------------------
        // compute MWM + AMD ordering
        Kokkos::Timer nd_mwm_amd_timer;
        Int a_first = btf_top_tabs_offset;
        Int nfirst = btf_tabs(a_first);
        Int nlast  = nfirst + BTF_A.ncol;
        auto order_nd_mwm = Kokkos::subview(order_blk_mwm_array,
                                            range_type(nfirst, nlast));
        auto order_nd_amd = Kokkos::subview(order_blk_amd_array,
                                            range_type(nfirst, nlast));
        #if 1
        // run MWM + AMD on A as one block
        int num_blks_A = 1;
        int save_tab = btf_tabs(1);
        //btf_tabs(1) = btf_tabs(btf_tabs_offset);
        btf_tabs(1) = BTF_A.ncol;
        #else
        // run MWM + AMD on each smaller blocks inside A
        int num_blks_A = btf_tabs_offset;
        int save_tab = btf_tabs(1);
        #endif

        // save AMD option
        BASKER_BOOL amd_dom = Options.amd_dom;
        Options.amd_dom = false; // since we do ND ?
        // save MWM option
        Int blk_matching_A = Options.blk_matching;
        if (Options.static_delayed_pivot != 0) {
          if(Options.verbose == BASKER_TRUE) {
            std::cout << "  -> switching MWM option from " << Options.blk_matching
                      << " to " << Options.static_delayed_pivot
                      << " for delayed pivots " << std::endl;
            if (!Options.amd_dom) {
              std::cout << "  -> turning off AMD option for ND block" << std::endl;
            }
          }
          Options.blk_matching = Options.static_delayed_pivot;
        }
        INT_1DARRAY blk_nnz;
        INT_1DARRAY blk_work;
        btf_blk_mwm_amd(0, num_blks_A, BTF_A,
                        order_nd_mwm, order_nd_amd,
                        blk_nnz, blk_work);
        if (Options.static_delayed_pivot != 0) {
          // revert MWM option
          Options.blk_matching = blk_matching_A;
        }
        #if 0 //debug:
        // reset for debug
        printf( " >> debug: set order_nd_mwm/amd on BFT_A to identity <<\n" );
        for (Int i = 0; i < (Int)BTF_A.nrow; i++) {
          order_nd_mwm(i) = i;
          order_nd_amd(i) = i;
        }
        #endif
        //printf( "order_nd_mwm\n");
        //for (Int j = 0; j < BTF_A.nrow; j++) printf( " %d %d\n",j, order_nd_mwm(j));

        if (Options.static_delayed_pivot != 0) {
          // ----------------------------------------------------------------------------------------------
          // skip applying AMD: btf_blk_mwm_amd does not do AMD with delayed pivots
          //for (Int i = 0; i < (Int)BTF_A.nrow; i++) {
          //  order_nd_amd(i) = i;
          //}
          if(Options.verbose == BASKER_TRUE) {
            std::cout<< " > Basker Factor: Time to compute MWM on a big block A: " << nd_mwm_amd_timer.seconds() << std::endl;
          }
        } else {
          // ---------------------------------------------------------------------------------------------
          // apply MWM on a big block A
          //for (int i=0; i<BTF_A.ncol; i++) printf( " - mxm_a(%d)=%d\n",i,order_nd_mwm(i));
          permute_row(BTF_A, order_nd_mwm);
          if (BTF_B.nrow > 0) {
            // Apply AMD perm to rows of B
            permute_row(BTF_B, order_nd_mwm);
          }

          // ----------------------------------------------------------------------------------------------
          // Apply AMD perm to rows and cols
          // NOTE: no need to udpate vals_order_blk_amd_array (it will read in A without permutations, and compute perm here)
          #if 0
          // skip applying AMD since we do ND (might helpful, e.g., one thread)
          for (Int i = 0; i < (Int)BTF_A.nrow; i++) {
            order_nd_amd(i) = i;
          }
          #endif
          if (Options.amd_dom) {
            permute_row(BTF_A, order_nd_amd);
            permute_col(BTF_A, order_nd_amd);
            if (btf_tabs_offset < btf_nblks) {
              // Apply AMD perm to rows and cols
              permute_row(BTF_B, order_nd_amd);
            }
            if (BTF_E.ncol > 0) {
              // Apply AMD perm to cols
              permute_col(BTF_E, order_nd_amd);
            }
          }
          // ----------------------------------------------------------------------------------------------
          // shift
          if(nfirst > 0) {
            for (Int i = 0; i < (Int)BTF_A.nrow; i++) {
              order_nd_mwm(i) += nfirst;
              order_nd_amd(i) += nfirst;
              //printf( " nd_mwm(%d) = %d, nd_amd(%d) = %d\n",i,order_nd_mwm(i), i,order_nd_amd(i) );
            }
          }
          if(Options.verbose == BASKER_TRUE) {
            std::cout<< " > Basker Factor: Time to compute and apply MWM+AMD on a big block A: " << nd_mwm_amd_timer.seconds() << std::endl;
          }
        }
        // reset AMD option
        Options.amd_dom = amd_dom;
        // reset tabs
        btf_tabs(1) = save_tab;

        int info_scotch = 0;
        BASKER_BOOL keep_zeros = true;
        Kokkos::Timer nd_nd_timer;
        if (Options.static_delayed_pivot != 0) {
          nd_mwm_amd_timer.reset();

          // instead of pivoting from higher separator in ND tree, we push the row to the separator
          // > also update sg.rangtab
          Int nblks = part_tree.nblks;
          INT_1DARRAY nd_map;
          INT_1DARRAY nd_sizes;
          MALLOC_INT_1DARRAY (nd_map,   BTF_A.nrow);
          MALLOC_INT_1DARRAY (nd_sizes, nblks+1);

          // original domain sizes
          //for (Int k = 0; k <= nblks; k++) {
          //  printf( " + row_tabs(%d) = %d\n",k,part_tree.row_tabs[k] );
          //}
          //for (Int j = 0; j < BTF_A.nrow; j++) printf( "permtab(%d) = %d, ipermtab(%d) = %d\n",j,part_tree.permtab[j], j,part_tree.ipermtab[j]);
          //printf( "\n" );
          nd_sizes(0) = 0;
          Kokkos::fence();
          Kokkos::parallel_for(
            "ndsort_matrix_store_valperms", num_threads,
            KOKKOS_LAMBDA(const int id) {
              for (Int k = id; k < nblks; k += num_threads) {
                for (Int i = part_tree.row_tabs[k]; i < part_tree.row_tabs[k+1]; i++) {
                  nd_map(i) = k;
                }
                nd_sizes(k+1) = part_tree.row_tabs[k+1] - part_tree.row_tabs[k];
              }
            });
          Kokkos::fence();
          if(Options.verbose == BASKER_TRUE) {
            printf( " > Initial ND block sizes:\n" );
            for (Int k = 0; k <= nblks; k++) {
              printf( " + nd_sizes(%d) = %d\n", (int)k, (int)nd_sizes(k) );
            }
            printf( "\n" );
          }

          // compute perm invert
          INT_1DARRAY iorder_nd_mwm;
          MALLOC_INT_1DARRAY (iorder_nd_mwm, BTF_A.nrow);
          for (Int id = 0; id < (Int)BTF_A.ncol; id++) {
            iorder_nd_mwm(order_nd_mwm(id)) = id;
          }
          //for (Int j = 0; j < BTF_A.nrow; j++) printf( " order_mwm(%d) = %d, iorder_mwm(%d) = %d\n",j, order_nd_mwm(j), j,iorder_nd_mwm(j));

          // look for invalid MWM swaps
          for (Int id = 0; id < (Int)BTF_A.ncol; id++) {
            Int j = id; //part_tree.permtab[id];
            Int k = iorder_nd_mwm(j);

            //printf( " j = %d (dom=%d), k = %d (dom=%d)\n",j,nd_map(j),k,nd_map(k) );
            // moving from dom1 (lower in ND tree) to dom2 (higher in ND tree)
            Int dom1 = nd_map(j);
            Int dom2 = nd_map(k);
            if (dom2 > dom1) 
            {
              nd_sizes(dom1+1) --;
              nd_sizes(dom2+1) ++;
            }
          }
          if(Options.verbose == BASKER_TRUE) {
            printf( " > ND block sizes after MWM and fixes:\n" );
            for (Int k = 0; k <= nblks; k++) {
              printf( " - nd_sizes(%d) = %d\n", (int)k, (int)nd_sizes(k) );
            }
            printf( "\n" );
          }

          // convert to offset
          for (Int k = 0; k < nblks; k++) {
            nd_sizes(k+1) += nd_sizes(k);
          }

          INT_1DARRAY order_nd_mwm2;
          MALLOC_INT_1DARRAY (order_nd_mwm2, BTF_A.nrow);
          for (Int id = 0; id < (Int)BTF_A.ncol; id++) {
            Int j = id; //part_tree.permtab[id];
            Int k = iorder_nd_mwm(j);

            // moving into dom at higher level in ND tree
            Int dom = (nd_map(j) > nd_map(k) ? nd_map(j) : nd_map(k));
            order_nd_mwm2(nd_sizes(dom)) = j;
            nd_sizes(dom) ++;
          }

          // save original ND
          if (part_tree_saved == BASKER_FALSE) {
            MALLOC_INT_1DARRAY(part_tree_orig.permtab,  BTF_A.ncol);
            MALLOC_INT_1DARRAY(part_tree_orig.ipermtab, BTF_A.ncol);
            for (Int j = 0; j < (Int)BTF_A.ncol; j++) {
              part_tree_orig.permtab[j] = part_tree.permtab[j];
              part_tree_orig.ipermtab[j] = part_tree.ipermtab[j];
            }

            MALLOC_INT_1DARRAY(part_tree_orig.row_tabs, 1+nblks);
            MALLOC_INT_1DARRAY(part_tree_orig.col_tabs, 1+nblks);
            for (Int k = 0; k <= nblks; k++) {
              part_tree_orig.row_tabs(k) = part_tree.row_tabs(k);
              part_tree_orig.col_tabs(k) = part_tree.col_tabs(k);
            }
            part_tree_saved = BASKER_TRUE;
          }

          // update permtab
          INT_1DARRAY order_nd_permtab;
          MALLOC_INT_1DARRAY (order_nd_permtab, BTF_A.nrow);
          for (Int id = 0; id < (Int)BTF_A.ncol; id++) {
            order_nd_permtab[id] = part_tree.permtab(order_nd_mwm2(id));
          }
          for (Int id = 0; id < (Int)BTF_A.ncol; id++) {
            part_tree.permtab[id] = order_nd_permtab(id);
          }

          // update ipermtab
          for (Int j = 0; j < (Int)BTF_A.ncol; j++) {
            part_tree.ipermtab[part_tree.permtab[j]] = j;
          }

          // update rangtab (row_tabs/col_tabs)
          for (Int k = nblks; k > 0; k--) {
            part_tree.row_tabs(k) = nd_sizes(k-1);
            part_tree.col_tabs(k) = nd_sizes(k-1);

            tree.row_tabs(k) = nd_sizes(k-1);
            tree.col_tabs(k) = nd_sizes(k-1);
          }
          part_tree.row_tabs(0) = 0;
          part_tree.col_tabs(0) = 0;
          //printf( "order_nd_mwm\n");
          //for (Int j = 0; j < BTF_A.nrow; j++) printf( " mwm (%d) = %d\n",j, order_nd_mwm(j));
          //for (Int j = 0; j < BTF_A.nrow; j++) printf( " mwm2(%d) = %d\n",j, order_nd_mwm2(j));
          //for (Int j = 0; j < BTF_A.nrow; j++) printf( "permtab(%d) = %d, ipermtab(%d) = %d\n",j,part_tree.permtab[j], j,part_tree.ipermtab[j]);
          //for (Int k = 0; k <= nblks; k ++) {
          //  printf( " > row_tabs(%d) = %d, col_tabs(%d) = %d\n",k,part_tree.row_tabs(k), k,part_tree.col_tabs(k) );
          //}
          //for (Int k = 0; k < nblks; k ++) {
          //  printf( " row_tabs[%d] -> %d\n",k,part_tree.row_tabs[k+1]-part_tree.row_tabs[k]);
          //}

          // ---------------------------------------------------------------------------------------------
          // apply MWM on a big block A (may break ND)
          permute_row(BTF_A, order_nd_mwm);
          if (BTF_B.nrow > 0) {
            // Apply AMD perm to rows of B
            permute_row(BTF_B, order_nd_mwm);
          }
          //for (int i=0; i<BTF_A.ncol; i++) printf( " - mxm_a(%d)=%d\n",i,order_nd_mwm(i));

          // ----------------------------------------------------------------------------------------------
          // apply nd-fix 
          INT_1DARRAY order_nd_amd2;
          MALLOC_INT_1DARRAY (order_nd_amd2, BTF_A.nrow);
          for (Int j = 0; j < (Int)BTF_A.ncol; j++) {
            order_nd_amd2(order_nd_mwm2(j)) = order_nd_amd(j);
          }

          // Apply fix to A
          permute_row(BTF_A, order_nd_amd2);
          permute_col(BTF_A, order_nd_amd2);
          if (btf_tabs_offset < btf_nblks) {
            // Apply fix to rows of B
            permute_row(BTF_B, order_nd_amd2);
          }
          if (BTF_E.ncol > 0) {
            // Apply fix to cols of E
            permute_col(BTF_E, order_nd_amd2);
          }
          if(Options.verbose == BASKER_TRUE) {
            std::cout<< "   + Basker Factor: Time to fix ND on a big block A: " << nd_mwm_amd_timer.seconds() << std::endl;
            nd_mwm_amd_timer.reset();
          }

          #if 0 // commented out (rely on csyamd in symbolic)
          // compute csymamd on fixed ND
          if (1) {
            BASKER_MATRIX AAT;
            AplusAT(BTF_A, AAT);

            Int nblks = tree.nblks;
            INT_1DARRAY tempp;
            INT_1DARRAY temp_col;
            INT_1DARRAY temp_row;
            MALLOC_INT_1DARRAY  (tempp,    AAT.ncol);
            MALLOC_INT_1DARRAY  (temp_col, AAT.ncol+nblks);
            MALLOC_INT_1DARRAY  (temp_row, AAT.nnz);
            if(Options.verbose == BASKER_TRUE) {
              std::cout<< "   + Basker Factor: Time to compute A+A^T: " << nd_mwm_amd_timer.seconds() << std::endl;
              nd_mwm_amd_timer.reset();
            }

            #if 1
            Int nleaves = num_threads;
            kokkos_amd_order<Int> amd_functor(nleaves, nblks, tree.col_tabs, AAT.col_ptr, AAT.row_idx,
                                              tempp, temp_col, temp_row, order_nd_amd, Options.verbose);
            Kokkos::parallel_for("BLK_AMD on A", Kokkos::RangePolicy<Exe_Space>(0, nleaves), amd_functor);
            Kokkos::fence();
            #else
            for(Int b = 0; b < tree.nblks; ++b) {
              Int frow = tree.col_tabs(b);
              Int erow = tree.col_tabs(b+1);
              Int fnnz = AAT.col_ptr(frow);

              Int nnz = 0;
              temp_col(frow+b) = 0;
              for(Int k = frow; k < erow; k++) {
                for(Int i = AAT.col_ptr(k); i < AAT.col_ptr(k+1); i++) {
                  if(AAT.row_idx(i) >= frow && AAT.row_idx(i) < erow) {
                    temp_row(fnnz + nnz) = AAT.row_idx(i) - frow;
                    nnz++;
                  }
                }
                temp_col(b+k+1) = nnz;
              }
              Int blk_size = erow - frow;
              double l_nnz = 0;
              double lu_work = 0;
              BaskerSSWrapper<Int>::amd_order(blk_size, &(temp_col(frow+b)), &(temp_row(fnnz)),
                                              &(tempp(frow)), l_nnz, lu_work, Options.verbose);
              for(Int k = 0; k < blk_size; k++)
              {
                order_nd_amd(frow+tempp(frow + k)) = frow+k;
              }
            }
            #endif
          } else {
            INT_1DARRAY cmember;
            MALLOC_INT_1DARRAY(cmember, BTF_A.ncol+1);
            for(Int i = 0; i < tree.nblks; ++i) {
              for(Int j = tree.col_tabs(i); j < tree.col_tabs(i+1); ++j) {
                cmember(j) = i;
              }
            }
            init_value(order_nd_amd, BTF_A.ncol, (Int)0);
            csymamd_order(BTF_A, order_nd_amd, cmember);
          }
          //for (int i=0; i<BTF_A.ncol; i++) printf( " - csymamd_a(%d)=%d\n",i,order_nd_amd(i));
          if(Options.verbose == BASKER_TRUE) {
            std::cout<< "   + Basker Factor: Time to compute AMD on a big block A: " << nd_mwm_amd_timer.seconds() << std::endl;
            nd_mwm_amd_timer.reset();
          }

          // apply csymamd
          permute_row(BTF_A, order_nd_amd);
          permute_col(BTF_A, order_nd_amd);
          if (btf_tabs_offset < btf_nblks) {
            // Apply AMD perm to rows
            permute_row(BTF_B, order_nd_amd);
          }
          if (BTF_E.ncol > 0) {
            // Apply AMD perm to cols
            permute_col(BTF_E, order_nd_amd);
          }

          // combine csymamd into fix
          INT_1DARRAY order_nd_amd3;
          MALLOC_INT_1DARRAY (order_nd_amd3, BTF_A.nrow);
          for (Int j = 0; j < (Int)BTF_A.ncol; j++) {
            order_nd_amd3(j) = order_nd_amd(order_nd_amd2(j));
          }
          for (Int j = 0; j < (Int)BTF_A.ncol; j++) {
            order_nd_amd(j) = order_nd_amd3(j);
          }
          if(Options.verbose == BASKER_TRUE) {
            std::cout<< "   + Basker Factor: Time to apply AMD on a big block A: " << nd_mwm_amd_timer.seconds() << std::endl;
            nd_mwm_amd_timer.reset();
          }
          #else
          for (Int j = 0; j < (Int)BTF_A.ncol; j++) {
            order_nd_amd(j) = order_nd_amd2(j);
          }
          #endif
          //for (int i=0; i<BTF_A.ncol; i++) printf( " - amd_a(%d)=%d\n",i,order_nd_amd(i));

          // ----------------------------------------------------------------------------------------------
          // shift
          if(nfirst > 0) {
            #if 0
            for (Int i = 0; i < (Int)BTF_A.nrow; i++) {
              order_nd_mwm(i) += nfirst;
              order_nd_amd(i) += nfirst;
            }
            #else
            Kokkos::parallel_for(
              "reset ndbtfa", BTF_A.nrow,
              KOKKOS_LAMBDA(const int i) {
                order_nd_mwm(i) += nfirst;
                order_nd_amd(i) += nfirst;
              });
            #endif
          }

	  // ----------------------------------------------------------------------------------------------
          // 'sort' rows of BTF_A into ND structure
          #if 0
          for (Int i = 0; i < BTF_A.nnz; ++i) {
            vals_order_ndbtfa_array(i) = i;
          }
          #else
          Kokkos::parallel_for(
            "reset ndbtfa", BTF_A.nnz,
            KOKKOS_LAMBDA(const int i) {
              vals_order_ndbtfa_array(i) = i;
            });
          Kokkos::fence();
          #endif
          BASKER_BOOL track_perm = BASKER_TRUE;
          ndsort_matrix_store_valperms(BTF_A, vals_order_ndbtfa_array, track_perm);
          if(Options.verbose == BASKER_TRUE) {
            std::cout<< "   + Basker Factor: Time to sort into ND blocks on a big block A: " << nd_mwm_amd_timer.seconds() << std::endl;
            nd_mwm_amd_timer.reset();
          }

          // ----------------------------------------------------------------------------------------------
          // Allocate & Initialize blocks
          #ifdef BASKER_KOKKOS
          matrix_to_views_2D(BTF_A);
          #endif
          if(Options.verbose == BASKER_TRUE) {
            std::cout<< "   + Basker Factor: Time to convert a big block A into views: " << nd_mwm_amd_timer.seconds() << std::endl;
          }
        } else {
          // ----------------------------------------------------------------------------------------------
          // compute & apply ND on a big block A
          info_scotch = apply_scotch_partition(keep_zeros);
        }
        if (info_scotch != BASKER_SUCCESS) {
          return info_scotch;
        }
        if(Options.verbose == BASKER_TRUE) {
          std::cout<< " > Basker Factor: Time to compute and apply ND on a big block A: " << nd_nd_timer.seconds() << std::endl;
        }

        nd_nd_timer.reset();
        // ----------------------------------------------------------------------------------------------
        // do symbolic & workspace allocation
        Kokkos::Timer nd_symbolic_timer;
        symmetric_sfactor();
        if(Options.verbose == BASKER_TRUE) {
          std::cout<< " > Basker Factor: Time for symbolic after ND on a big block A: " << nd_symbolic_timer.seconds() << std::endl;
        }

        Kokkos::Timer nd_last_dense_timer;
        bool flag = true;
        btf_last_dense(flag);
        if(Options.verbose == BASKER_TRUE) {
          std::cout<< " > Basker Factor: Time for last-dense after ND on a big block A: " << nd_last_dense_timer.seconds() << std::endl;
        }


        #ifdef BASKER_KOKKOS
        // ----------------------------------------------------------------------------------------------
        // Allocate & Initialize blocks
        kokkos_sfactor_init_factor<Int,Entry,Exe_Space>
          iF(this);
        Kokkos::parallel_for(TeamPolicy(num_threads,1), iF);
        Kokkos::fence();

        /*kokkos_sfactor_init_workspace<Int,Entry,Exe_Space>
          iWS(flag, this);
        Kokkos::parallel_for(TeamPolicy(num_threads,1), iWS);
        Kokkos::fence();*/
        #endif
        if(Options.verbose == BASKER_TRUE) {
          std::cout<< " > Basker Factor: Time for init factors after ND on a big block A: " << nd_nd_timer.seconds() << std::endl;
        }

        if(Options.verbose == BASKER_TRUE) {
          std::cout<< "Basker Factor: Time to compute ND and setup: " << nd_perm_timer.seconds() << std::endl << std::endl;
        }
      } // end of ND & setups

      if ((Options.static_delayed_pivot != 0 || Options.blk_matching != 0) && btf_nblks > btf_tabs_offset) {
        Kokkos::Timer mwm_amd_perm_timer;

        // save MWM option
        Int blk_matching_A = Options.blk_matching;
        if (Options.static_delayed_pivot != 0) {
          if(Options.verbose == BASKER_TRUE) {
            std::cout << "  -> switching MWM option from " << Options.blk_matching
                      << " to " << Options.static_delayed_pivot
                      << " for delayed pivots " << std::endl;
            if (!Options.amd_dom) {
              std::cout << "  -> turning off AMD option for ND block" << std::endl;
            }
            std::cout << std::endl;
          }
          Options.blk_matching = Options.static_delayed_pivot;
        }
        Int b_first = btf_tabs_offset;
        Int nfirst = btf_tabs(b_first);
        auto order_blk_mwm_c = Kokkos::subview(order_blk_mwm_array,
                                               range_type(nfirst, ncol));
        auto order_blk_amd_c = Kokkos::subview(order_blk_amd_array,
                                               range_type(nfirst, ncol));
        if(Options.verbose == BASKER_TRUE) {
          std::cout << " calling MWM on C" << std::endl;
          std::cout << " btf_blk_mwm: btf_tabs(" << btf_tabs_offset << ")=" << nfirst << std::endl;
        }

        // ----------------------------------------------------------------------------------------------
        // recompute MWM and AMD on each block of C
        INT_1DARRAY blk_nnz;
        INT_1DARRAY blk_work;
        for (Int i = b_first; i <= btf_nblks; i++) {
          //printf( " > btf_tabs(%d) = %d -> %d\n",i,btf_tabs(i),btf_tabs(i)-nfirst );
          btf_tabs(i) -= nfirst;
        }
        btf_blk_mwm_amd(b_first, btf_nblks-b_first, BTF_C,
                        order_blk_mwm_c, order_blk_amd_c,
                        blk_nnz, blk_work);
        if (Options.static_delayed_pivot != 0) {
          // revert MWM option
          Options.blk_matching = blk_matching_A;
        }
        #if 0 //debug
        printf( " >> debug: set order_blk_mwm/amd on BTF_C to identity <<\n" );
        for (Int i = 0; i < (Int)BTF_C.nrow; i++) {
          order_blk_mwm_c(i) = i;
          order_blk_amd_c(i) = i;
        }
        #endif
        //for (Int i = 0; i < (Int)BTF_C.nrow; i++) printf( " x mwm_c(%d)=%d amd_c(%d)=%d\n",i,order_blk_mwm_c(i), i,order_blk_amd_c(i));

        // Apply MWM perm to rows
        permute_row(BTF_C, order_blk_mwm_c);

        // Apply AMD perm to rows and cols
        // NOTE: no need to udpate vals_order_blk_amd_array (it will read in A without permutations, and compute perm here)
        permute_row(BTF_C, order_blk_amd_c);
        permute_col(BTF_C, order_blk_amd_c);
        if (BTF_E.ncol > 0) {
          // Apply AMD perm to cols
          permute_col(BTF_E, order_blk_amd_c, BTF_A.ncol);
        }
        if (btf_tabs_offset > 0) {
          // Apply AMD perm to cols
          permute_col(BTF_B, order_blk_amd_c);
        }

        // shift back perm to global index
        for (Int i = b_first; i <= btf_nblks; i++) {
          btf_tabs(i) += nfirst;
        }
        for (Int i = 0; i < (Int)BTF_C.ncol; i++) {
          order_blk_mwm_c(i) += nfirst;
          order_blk_amd_c(i) += nfirst;
        }

        if(Options.verbose == BASKER_TRUE) {
          std::cout << std::endl << std::endl;
          std::cout << "Basker Factor: Time to compute and apply MWM+AMD on diagonal blocks: " << mwm_amd_perm_timer.seconds() << std::endl;
        }
      }
      //for (Int i = 0; i < (Int)ncol; i++) {
      //  printf( " - mwm(%d)=%d amd(%d)=%d\n",i,order_blk_mwm_array(i), i,order_blk_amd_array(i));
      //}
      //for (Int i = 0; i < (Int)BTF_A.ncol; i++) {
      //  printf( " + nd(%d)=%d camd(%d)=%d\n",i,part_tree.permtab(i), i,order_csym_array(i));
      //}

      // compose row permute
      //printf( " > ND and CSYMAMD in numeric setup\n" );
      //for (int i=0; i<gn; i++) printf( " > numeric_row_iperm[%d] = %d\n",i,numeric_row_iperm_array(i) );
      //printf( "  + permtab, csym_array\n" );
      //for (int i=0; i<BTF_A.ncol; i++) printf( " - %d %d %d\n",i, part_tree.permtab(i),   order_csym_array(i) );
      //printf( "  + blk_mwm, blk_amd\n" );
      //for (int i=0; i<ncol;       i++) printf( " + %d %d %d\n",i, order_blk_mwm_array(i), order_blk_amd_array(i) );
      permute_inv_with_workspace(numeric_row_iperm_array, order_blk_mwm_array, ncol);
      permute_inv_with_workspace(numeric_row_iperm_array, order_blk_amd_array, ncol);
      if (BTF_A.ncol > 0 && Options.static_delayed_pivot == 0) {
        Int scol_top = btf_tabs[btf_top_tabs_offset]; 

        //printf( " A.scol = %d, scol_top = %d (btf_tabs = %d, btf_top_tabs = %d)\n",BTF_A.scol, scol_top, btf_tabs_offset, btf_top_tabs_offset );
        //for (Int i = 0; i < btf_top_tabs_offset; i++) printf( " x %d: %d\n",i,btf_tabs[i+1]-btf_tabs[i] );
        //for (Int i = btf_top_tabs_offset; i < btf_tabs_offset; i++) printf( " + %d: %d\n",i,btf_tabs[i+1]-btf_tabs[i] );
        //for (Int i = btf_tabs_offset; i < btf_nblks; i++) printf( " - %d: %d\n",i,btf_tabs[i+1]-btf_tabs[i] );

        permute_inv_with_workspace(numeric_row_iperm_array, part_tree.permtab, BTF_A.ncol, scol_top);
        permute_inv_with_workspace(numeric_row_iperm_array,  order_csym_array, BTF_A.ncol, scol_top);
      }
      // (int i=0; i<gn; i++) printf( " < numeric_row_iperm[%d] = %d\n",i,numeric_row_iperm_array(i) );
      // compose col permute
      if (BTF_A.ncol > 0 && Options.static_delayed_pivot == 0) {
        Int scol_top = btf_tabs[btf_top_tabs_offset]; 
        permute_with_workspace(numeric_col_iperm_array,  order_csym_array, BTF_A.ncol, scol_top);
        permute_with_workspace(numeric_col_iperm_array, part_tree.permtab, BTF_A.ncol, scol_top);
      }
      permute_with_workspace(numeric_col_iperm_array, order_blk_amd_array, ncol);
      //for (int i = 0; i < ncol; i++) printf( " %d: %d %d\n",i, numeric_row_iperm_array(i),numeric_col_iperm_array(i) );
    }
    else {
      if(Options.verbose == BASKER_TRUE) {
        std::cout << " No MWM on diagonal blocks " << std::endl;
      }
      if (Options.matrix_scaling != 0)
      {
        using STS = Teuchos::ScalarTraits<Entry>;
        using Real = typename STS::magnitudeType;
        const Real one (1.0);
        const Real eps = Teuchos::ScalarTraits<Entry>::eps ();
        const Real normA = BTF_A.anorm;

        Int scol_top = btf_tabs[btf_top_tabs_offset];

        ////////////////////
        // row and column scale A
        if (Options.matrix_scaling == 1) {
          // find diagonal scaling factors 
          for(Int j = 0; j < BTF_A.ncol; j++) {
            Int col = scol_top+j;
            scale_row_array[col] = one;
            scale_col_array[col] = one;
            for(Int k = BTF_A.col_ptr[j]; k < BTF_A.col_ptr[j+1]; k++) {
              if (BTF_A.row_idx[k] == j) {
                Entry ajj = abs(BTF_A.val[k]);
                if (abs(ajj) <= eps*normA) {
                  ajj = eps*normA;
                }
                scale_row_array[col] = one / sqrt(ajj);
                scale_col_array[col] = one / sqrt(ajj);
              }
            }
          }

          // scale matrix with symmetric diagonal scale
          for(Int j = 0; j < BTF_A.ncol; j++) {
            Int col = scol_top+j;
            for(Int k = BTF_A.col_ptr[j]; k < BTF_A.col_ptr[j+1]; k++) {
              Int row = scol_top+BTF_A.row_idx[k];
              BTF_A.val[k] *= scale_col_array[col];
              BTF_A.val[k] *= scale_row_array[row];
            }
          }
        } else {
          // find max non-zero entry on row
          for(Int j = 0; j < BTF_A.ncol; j++) {
            for(Int k = BTF_A.col_ptr[j]; k < BTF_A.col_ptr[j+1]; k++) {
              Int row = scol_top+BTF_A.row_idx[k];
              if (abs(scale_row_array[row]) < abs(BTF_A.val[k])) {
                scale_row_array[row] = abs(BTF_A.val[k]);
              }
            }
          }
          for(Int i = 0; i < BTF_A.ncol; i++) {
            Int row = scol_top+i;
            if (abs(scale_row_array[row]) <= eps*normA) {
              scale_row_array[row] = eps*normA;
            }
            scale_row_array[row] = one / scale_row_array[row];
          }

          // scale matrix with row scale
          for(Int j = 0; j < BTF_A.ncol; j++) {
            for(Int k = BTF_A.col_ptr[j]; k < BTF_A.col_ptr[j+1]; k++) {
              Int row = scol_top+BTF_A.row_idx[k];
              BTF_A.val[k] *= scale_row_array[row];
            }
          }

          // find max non-zero entry on col, scale matrix
          for(Int j = 0; j < BTF_A.ncol; j++) {
            // find max entry
            Int col = scol_top+j;
            for(Int k = BTF_A.col_ptr[j]; k < BTF_A.col_ptr[j+1]; k++) {
              if (abs(scale_col_array[col]) < abs(BTF_A.val[k])) {
                scale_col_array[col] = abs(BTF_A.val[k]);
              }
            }
            if (abs(scale_col_array[col]) <= eps*normA) {
              scale_col_array[col] = eps*normA;
            }
            // scale
            scale_col_array[col] = one / scale_col_array[col];
            for(Int k = BTF_A.col_ptr[j]; k < BTF_A.col_ptr[j+1]; k++) {
              BTF_A.val[k] *= scale_col_array[col];
            }
          }
        }

        ////////////////////
        // row-scale B
        for(Int j = 0; j < BTF_B.ncol; j++) {
          for(Int k = BTF_B.col_ptr[j]; k < BTF_B.col_ptr[j+1]; k++) {
            Int row = scol_top+BTF_B.row_idx[k];
            BTF_B.val[k] *= scale_row_array[row];
          }
        }
        ////////////////////
        // col-scale E
        for(Int j = 0; j < BTF_E.ncol; j++) {
          Int col = scol_top+j;
          for(Int k = BTF_E.col_ptr[j]; k < BTF_E.col_ptr[j+1]; k++) {
            BTF_E.val[k] *= scale_col_array[col];
          }
        }
      }
    }
    reset_error();

    #ifdef BASKER_DEBUG_DIAG
    // Look for zeros on the diagonal - an Option / parameter list option should be used to enable this
    // col_ptr has ncol+1 entries - final entry = nnz
    // row_idx and val have nnz entries
    // Scan through the matrix for 0's on the diagonal
    if ( nrow == ncol ) // square matrices, CRS vs CCS is irrelevant (i.e. transpose does not need to be applied first) - but needs to be known for Symbolic...
    {
      std::cout << "Basker Factor(...) debugging:  Check diagonal entries for zeros" << std::endl;
      Int total_diag_found = 0;
      Int total_zero_diag_found = 0;
      double min = 100000000;

      for ( Int colid = 0; colid < ncol; ++colid ) {
        Int rowids_start = A.col_ptr[colid];
        Int rowids_end = A.col_ptr[colid+1];
        Int num_row_entries = rowids_end - rowids_start;
        Int counter_off_diag = 0;

        for ( Int j = rowids_start; j < rowids_end; ++j ) {
          Int rowid = A.row_idx[j];
          if ( rowid != colid ) {
            ++counter_off_diag; // track how many are off-diagonal entries
          }
          if ( rowid == colid ) {
            ++total_diag_found; // diagonal entry is present; next check if its value is 0

            //check non-zero - this occurs when there IS an entry for diagonal
            if ( A.val[j] == static_cast<Entry>(0) ) {
              ++total_zero_diag_found;
              std::cout << "ZERO ON DIAG!!!" << std::endl;
            }
            if ( std::abs(A.val[j]) < min ) {
              min = std::abs(A.val[j]); // track minimum values (absolute value taken to see entries close to 0)
            }
          }
        } // end scanning through row entries in column colid

        if ( counter_off_diag + 1 != num_row_entries ) { // if this happens, a diag entry was never found for this column - assumed 0
          std::cout << "ZERO DIAG Found!" << std::endl;
          std::cout << " colid = " << colid << "  counter_off_diag = " << counter_off_diag << "  num_row_entries = " << num_row_entries << std::endl;
        }
      }
      if ( total_diag_found == ncol ) {
        std::cout << "NO MISSING DIAG ENTRY!  Found: " << total_diag_found << "  expected: " << ncol << std::endl;
      }
      else {
        std::cout << "MISSING DIAG ENTRY!  Non-zero entries found: " << total_diag_found << "  expected: " << ncol << std::endl;
      }
      std::cout << "Total zero diags found: " << total_zero_diag_found << std::endl;
      std::cout << "Min diag entry = " << min << std::endl;
    }
    #endif

    // sfactor_copy2 is now only responsible for the copy from BTF_A to 2D blocks
    Kokkos::Timer timer_sfactorcopy;
    double sfactorcopy_time = 0.0;
    if (btf_tabs_offset != 0) {
      bool flag = true;
      #ifdef BASKER_KOKKOS
      Kokkos::Timer nd_setup1_timer;
      /*kokkos_sfactor_init_factor<Int,Entry,Exe_Space>
        iF(this);
      Kokkos::parallel_for(TeamPolicy(num_threads,1), iF);
      Kokkos::fence();
      if(Options.verbose == BASKER_TRUE) {
        std::cout<< " > Basker Factor: Time for init factors after ND on a big block A: " << nd_setup1_timer.seconds() << std::endl;
      }*/

      Kokkos::Timer nd_setup2_timer;
      kokkos_sfactor_init_workspace<Int,Entry,Exe_Space>
        iWS(flag, this);
      Kokkos::parallel_for(TeamPolicy(num_threads,1), iWS);
      Kokkos::fence();
      if(Options.verbose == BASKER_TRUE) {
        std::cout<< " > Basker Factor: Time for workspace allocation after ND on a big block A: " << nd_setup2_timer.seconds() << std::endl;
      }
      #endif
    }
    bool copy_BTFA = (Options.blk_matching == 0 || Options.static_delayed_pivot != 0);
    bool alloc_BTFA = (Options.static_delayed_pivot != 0);
    err = sfactor_copy2(alloc_BTFA, copy_BTFA);

    if(Options.verbose == BASKER_TRUE) {
      sfactorcopy_time += timer_sfactorcopy.seconds();
      std::cout << "Basker Factor sfactor_copy2 time: " << sfactorcopy_time << std::endl;
      std::cout << " >> error = " << err << std::endl;
    }
    if(err == BASKER_ERROR)
    { return BASKER_ERROR; }

    Kokkos::Timer timer_factornotoken;
    double fnotoken_time = 0.0;
    if(Options.incomplete == BASKER_FALSE)    
    {
      err = factor_notoken(0);
    }
    else
    {
      err = factor_inc_lvl(0);
    }


    if(Options.verbose == BASKER_TRUE) {
      fnotoken_time += timer_factornotoken.seconds();
      std::cout << "Basker factor_notoken total time: " << fnotoken_time << std::endl;
    }

    if(err == BASKER_ERROR)
    {
      printf("ShyLUBasker factor_notoken/inc_lvl error returned (err=%d)\n",err);
      return BASKER_ERROR; 
    }

    if ( sym_gn != gn || sym_gm != gm ) {
      printf( "ShyLUBasker Factor error: Matrix dims at Symbolic and Factor stages do not agree (sym_gm=%d, gn=%d, gm=%d)",(int)sym_gn,(int)gn,(int)gm);
      printf( " - Symbolic reordered structure will not apply.\n");
      //exit(EXIT_FAILURE);
      return BASKER_ERROR; 
    }

    if(Options.verbose == BASKER_TRUE)
    { printf(" == Basker Factor Done ==\n\n"); }

    //DEBUG_PRINT();

    #ifdef BASKER_TIMER
    time += timer.seconds();
    std::cout << "Basker Factor total time: " << time << std::endl;
    std::cout.precision(old_precision);
    std::cout.flags(old_settings);
    #endif

    factor_flag = BASKER_TRUE;

    return 0;
  }//end Factor()


  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::Factor_Inc(Int _Options)
  {
    factor_inc_lvl(_Options);

    return 0;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Solve(Entry *b, Entry *x, bool transpose)
  {
    #ifdef BASKER_TIMER 
    Kokkos::Timer timer;
    double time = 0.0;
    #endif

    if(Options.verbose == BASKER_TRUE)
    {
      printf("Basker Solve Called (%s)\n",(transpose ? " transpose" : "non-transpose"));
    }

    if(factor_flag != BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE) {
        printf("BASKER: YOU NEED TO RUN FACTOR  BEFORE SOLVE\n");
      }
      return BASKER_ERROR;
    }

    if (transpose == false)
      solve_interface(x,b);
    else
      solve_interfacetr(x,b);

    if(Options.verbose == BASKER_TRUE)
    {
      printf("Basker Solve Done \n");
    }

    #ifdef BASKER_TIMER
    time += timer.seconds();
    std::cout << "Basker Solve total time: " << time << std::endl;
    #endif

    return 0;
  }//Solve(Entry *, Entry *);


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::Solve(Int _nrhs, Entry *b, Entry *x, bool transpose)
  {
    #ifdef BASKER_TIMER 
    Kokkos::Timer timer;
    double time = 0.0;
    #endif

    if(Options.verbose == BASKER_TRUE)
    {
      printf("Basker MultiSolve Called with %d RHSs (%s)\n",(int)_nrhs,(transpose ? "transpose" : "non-transpose"));
    }

    if(factor_flag != BASKER_TRUE)
    {
      if(Options.verbose == BASKER_TRUE) {
        printf("BASKER: YOU NEED TO RUN FACTOR  BEFORE SOLVE\n");
      }
      return BASKER_ERROR;
    }

    if (transpose == false)
      solve_interface(_nrhs,x,b);
    else
      solve_interfacetr(_nrhs,x,b);

    if(Options.verbose == BASKER_TRUE)
    {
      printf("Basker Multisolve Done \n");
    }
    
    #ifdef BASKER_TIMER
    time += timer.seconds();
    std::cout << "Basker Solve total time: " << time << std::endl;
    #endif

    return 0;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Solve(ENTRY_1DARRAY b, ENTRY_1DARRAY x, bool transpose)
  {
    printf("Basker: This solve call not implemented\n");
    return -1;
  }//Solve(ENTRY_1D, ENTRY_1D);


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Solve(Int _nrhs, Entry *b, Entry *x, Int option, bool transpose)
  {    
    int err = 0;
    printf("Basker: This solve call not implemented\n");
    if(solve_flag == false) //never solved before
    {
      //err = malloc_init_solve(_nrhs, x, b);
    }
    if(solve_flag == true) //fix data
    {
      //Come back to add options for this case
      return -1;
    }

    //err = solve(sol,rhs);

    return err;
  }//end Solve()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::SetThreads(Int nthreads)
  {
    //Need to test if power of nparts
    //TODO: hard-coded to be two. It is also hard-coded in shylubasker_structs.hpp
    double nparts = 2.0;
    if (pow(nparts, log((double)nthreads)/log(nparts)) != nthreads)
    {
      BASKER_ASSERT(0==1, "Basker SetThreads Assert: Number of thread error - not a power of 2");
      //Set default 1
      num_threads = 1;
      return BASKER_ERROR;
    }

    //Next test if Kokkos has that many threads!
    //This is a common mistake in mpi-based apps
    #ifdef KOKKOS_ENABLE_OPENMP
    int check_value = Kokkos::OpenMP::impl_max_hardware_threads();
    if(nthreads > check_value)
    {
      BASKER_ASSERT(0==1, "Basker SetThreads Assert: Number of thread not available");
      num_threads =  1;
      return BASKER_ERROR;
    }
    #else
    nthreads = 1;
    #endif

    num_threads = nthreads;
    return BASKER_SUCCESS;
  }//end SetThreads()


  //Return nnz of L
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::GetLnnz(Int &Lnnz)
  {
    (Lnnz) = get_Lnnz();
    if(Lnnz == 0)
    { return BASKER_ERROR; }
    else
    { return BASKER_SUCCESS; }
  }//end GetLnnz();


  //Return nnz of U
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::GetUnnz(Int &Unnz)
  {
    (Unnz) = get_Unnz();
    if(Unnz == 0)
    { return BASKER_ERROR; }
    else
    { return BASKER_SUCCESS; }
  }//end GetUnnz()


  //Returns assembled L
  template<class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::GetL
  (
   Int &n, 
   Int &nnz, 
   Int **col_ptr, 
   Int **row_idx, 
   Entry **val
  )
  {
    get_L(n,nnz,col_ptr, row_idx, val);
    
    return BASKER_SUCCESS;
  }//end GetL()
  

  //returns assembles U
  template<class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::GetU
  (
   Int &n, 
   Int &nnz, 
   Int **col_ptr, 
   Int **row_idx, 
   Entry **val
  )
  {
    get_U(n, nnz, col_ptr, row_idx, val);
    return BASKER_SUCCESS;
  }//end GetU()


  //returns global P
  template<class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::GetPerm(Int **lp, Int **rp)
  {
    INT_1DARRAY lp_array;
    MALLOC_INT_1DARRAY(lp_array, gn);
    INT_1DARRAY rp_array;
    MALLOC_INT_1DARRAY(rp_array, gn);

    get_total_perm(lp_array, rp_array);

    (*lp) = new Int[gn];
    (*rp) = new Int[gn];

    for(Int i = 0; i < gn; ++i)
    {
      (*lp)[i] = lp_array(i);
      (*rp)[i] = rp_array(i);
    }

    FREE_INT_1DARRAY(lp_array);
    FREE_INT_1DARRAY(rp_array);

    return BASKER_SUCCESS;
  }//end GetPerm()


  //Timer Information function
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int, Entry, Exe_Space>::PrintTime()
  {
    // stats.print_time();

    /*
    print_local_time_stats();

    std::cout << std::endl 
              << "---------------TIME-------------------"<< std::endl
              << "Tree:    " << stats.tree_time    << std::endl
              << "SFactor: " << stats.sfactor_time << std::endl
              << "Nfactor: " << stats.nfactor_time << std::endl
              << "LSolve:  " << stats.lower_solve_time << std::endl
              << "USolve:  " << stats.upper_solve_time << std::endl
              << "-------------END TIME------------------"
              << std::endl << std::endl;

    stats.sfactor_time = 0;
    stats.nfactor_time = 0;
    stats.lower_solve_time = 0;
    stats.upper_solve_time = 0;
    */
  }


  //Debug tester function
  template<class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::DEBUG_PRINT()
  {
    //print_sep_bal();

    #ifdef BASKER_2DL
    printL2D();
    printLMTX();
    #else
    //printL();
    #endif
    std::cout << "L printed " << std::endl;
    printU();
    printUMTX();
    std::cout << "U printed" << std::endl;
    //printRHS();
    std::cout << "RHS printed" << std::endl;
    //printSOL();
    std::cout << "SOL printed" << std::endl;
    //printTree();
    std::cout << "Tree printed" << std::endl;

    //Print out vectors
    if(match_flag == BASKER_TRUE)
    {
      printVec("match.csc", order_match_array,
          order_match_array.extent(0));
    }
    if(btf_flag == BASKER_TRUE)
    {
      printVec("btf.csc", order_btf_array,
          order_btf_array.extent(0));
      printVec("amdblk.csc", order_blk_amd_array,
          order_blk_amd_array.extent(0));
    }
    if(btf_tabs_offset != 0)
    {
      printVec("ND.csc", part_tree.permtab, 
          part_tree.permtab.extent(0));
    }
    if(amd_flag == BASKER_TRUE)
    {
      printVec("camd.csc", order_csym_array,
          order_csym_array.extent(0));
    }
  }//end DEBUG_PRINT()

  template<class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::PRINT_C()
  {
    BASKER_MATRIX  &M = BTF_C;
    printf( " > B = [\n" );
    for(Int k = 0; k < M.ncol; ++k) {
      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i) {
        printf( "%d %d %e\n", M.row_idx(i), k, M.val(i));
      }
    }
    printf( " ];\n" );
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Info()
  {
    std::cout << "---------BASKER <2D>---------" 
              << "---------V.  0.0.3 ------------- "
              << "Written by Joshua Dennis Booth"
              << "jdbooth@sandia.gov"
              << "Sandia National Labs"
              << "---------------------------------"
              << std::endl;
    return 0;
  }//end info

}//End namespace

#undef BASKER_TIMER
#endif //End ifndef
