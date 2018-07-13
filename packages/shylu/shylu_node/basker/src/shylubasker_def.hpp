#ifndef SHYLUBASKER_DEF_HPP
#define SHYLUBASKER_DEF_HPP

//#define BASKER_TIME

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

/*Kokkos Includes*/
#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#else
#include <omp.h>
#endif

/*System Includes*/
#include <iostream>

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

    //Default number of threads
    num_threads = 1;
    global_nnz  = 0;
    gn = 0;

    btf_total_work = 0;
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
    FREE_MATRIX_VIEW_2DARRAY(AV, tree.nblks);
    FREE_MATRIX_VIEW_2DARRAY(AL, tree.nblks);
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

      FREE_INT_1DARRAY(vals_order_ndbtfa_array); //track nd perms; BTF_A must be declared here, else it does not exist
      FREE_INT_1DARRAY(vals_order_ndbtfb_array); //track nd perms; BTF_A must be declared here, else it does not exist
      FREE_INT_1DARRAY(vals_order_ndbtfc_array); //track nd perms; BTF_A must be declared here, else it does not exist
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
    MALLOC_INT_1DARRAY(perm_inv_comp_array , sym_gm); //y
    MALLOC_INT_1DARRAY(perm_comp_array, sym_gn); //x

    MALLOC_INT_1DARRAY(perm_comp_iworkspace_array, sym_gn); 
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
   bool crs_transpose_needed
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
      std::cout << "Basker Symbolic" << std::endl;
      std::cout << "Matrix dims: " << nrow << " " << ncol << " " << nnz << std::endl;
    }
    //Init Matrix A.
    if(matrix_flag == BASKER_TRUE)
    {
      printf("BASKER: YOU CANNOT RERUN SYMBOLIC\n");
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

      if(Options.transpose == BASKER_FALSE)
      {
        A.init_matrix("Original Matrix", nrow, ncol, nnz, col_ptr, row_idx, val);
        A.scol = 0;
        A.srow = 0;

        // NDE: New for Amesos2 CRS mods
        if ( crs_transpose_needed ) {
          matrix_transpose(0, nrow, 0, ncol, nnz, col_ptr, row_idx, val, A, vals_crs_transpose);
        }
      }
      else
      {
        //Will transpose and put in A using little extra
        //if ( crs_transpose_needed ) {
        //  matrix_transpose(0, nrow, 0, ncol, nnz, col_ptr, row_idx, val, A, vals_crs_transpose);
        //}
        matrix_transpose(0, nrow, 0, ncol, nnz, col_ptr, row_idx, val, A);
        // NDE: How should transpose be handled special case (when CRS format comes in)? null op, i.e. transpose of transpose yields original input?
      }
      sort_matrix(A);

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
      std::cout << "Basker Symbolic matrix init time: " << init_time << std::endl;
      #endif
    }

    //Init Ordering
    //Always will do btf_ordering
    //This should also call create tree
    if(order_flag == BASKER_TRUE)
    {
      printf("BASKER: YOU CANNOT RERUN ORDER\n");
      return BASKER_ERROR;
    }
    else
    {
      //btf_order();
      #ifdef BASKER_TIMER
      double order_time = 0.0;
      Kokkos::Timer timer_order;
      #endif
      /*
         if(Options.incomplete == BASKER_TRUE)
         {
           order_incomplete();
         }
         else
         {
           btf_order2();
         }
      */

      btf_order2();

      if(Options.verbose == BASKER_TRUE)
      {
        printf("Basker Ordering Found \n");
      }

      if((Options.btf == BASKER_TRUE) && (btf_tabs_offset != 0))
      {
        basker_barrier.init(num_threads, 16, tree.nlvls );
      }
      order_flag = BASKER_TRUE;

      if(Options.verbose == BASKER_TRUE)
      {
        printf("Basker P2P Thread Barriers Init\n");
      }

      #ifdef BASKER_TIMER
      order_time += timer_order.seconds();
      std::cout << "Basker Symbolic order arrays time: " << order_time << std::endl;
      #endif
    }

    if(symb_flag == BASKER_TRUE)
    {
      printf("BASKER: YOU CANNOT RERUN SFACTOR\n");
      return BASKER_ERROR;
    }
    else
    {
      #ifdef BASKER_TIMER
      double sfactor_time = 0.0;
      Kokkos::Timer timer_sfactor;
      #endif
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
      printf("Basker Symbolic Done \n");
    }

    #ifdef BASKER_TIMER
    time = timer.seconds();
    stats.time_sfactor += time;
    std::cout << "Basker Symbolic total time: " << time << std::endl;
    std::cout.precision(old_precision);
    std::cout.flags(old_settings);
    #endif

    // NDE store matrix dims here for comparison in Factor
    sym_gn = A.ncol;
    sym_gm = A.nrow;
    MALLOC_ENTRY_1DARRAY(x_view_ptr_copy, sym_gn); //used in basker_solve_rhs
    MALLOC_ENTRY_1DARRAY(y_view_ptr_copy, sym_gm);
    MALLOC_INT_1DARRAY(perm_inv_comp_array , sym_gm); //y
    MALLOC_INT_1DARRAY(perm_comp_array, sym_gn); //x

    MALLOC_INT_1DARRAY(perm_comp_iworkspace_array, sym_gn); 
    MALLOC_ENTRY_1DARRAY(perm_comp_fworkspace_array, sym_gn);
    permute_composition_for_solve(sym_gn);

    return 0;
  }//end Symbolic()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Factor(Int option)
  {
    #ifdef BASKER_TIMER
    Kokkos::Timer timer;
    double time = 0.0;
    #endif

    factor_notoken(option);

    #ifdef BASKER_TIMER
    time += timer.seconds();
    stats.time_nfactor += time;
    std::cout << "Basker factor_notoken time: " << time << std::endl;
    timer.reset();
    #endif

    #ifdef BASKER_TIMER
    time += timer.seconds();
    std::cout << "Basker Factor total time: " << time << std::endl;
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
    #ifdef BASKER_TIMER
    std::ios::fmtflags old_settings = cout.flags();
    int old_precision = std::cout.precision();
    std::cout.setf(ios::fixed, ios::floatfield);
    std::cout.precision(8);
    double time = 0.0;
    Kokkos::Timer timer;
    #endif

    int err = 0; //init for return value from sfactor_copy2, factor_notoken etc.

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

    #ifdef BASKER_TIMER
    Kokkos::Timer copyperm_timer;
    #endif

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

    if ( btf_nblks > 1 ) { //non-single block case
    #ifdef KOKKOS_ENABLE_OPENMP
    #pragma omp parallel for
    #endif
      for( Int i = 0; i < nnz; ++i ) {
        A.val(i) = val[ vals_perm_composition(i) ];
        if ( btfa_nnz != 0 ) { //is this unnecessary? yes, but shouldn't the first label account for this anyway?
          if ( vals_block_map_perm_pair(i).first == 0 ) { //in BTF_A
            BTF_A.val( inv_vals_order_ndbtfa_array( vals_block_map_perm_pair(i).second ) ) = val[ vals_perm_composition(i) ];
          }
        }
        if ( btfb_nnz != 0 ) {
          if ( vals_block_map_perm_pair(i).first == 1 ) { //in BTF_B
            BTF_B.val( inv_vals_order_ndbtfb_array( vals_block_map_perm_pair(i).second ) ) = val[ vals_perm_composition(i) ];
          }
        }
        //      if ( BTF_C.nnz != 0 ) // this causes compiler error, and nnz blocks values are different with this command with small blocks matrix (not nd) - why?
        if ( btfc_nnz != 0 ) {
          if ( vals_block_map_perm_pair(i).first == 2 ) { //in BTF_C
            BTF_C.val( inv_vals_order_ndbtfc_array( vals_block_map_perm_pair(i).second ) ) = val[ vals_perm_composition(i) ];
          }
        }
      } //end for
    } //end if
    else if ( btf_nblks == 1 )
    {
    #ifdef KOKKOS_ENABLE_OPENMP
    #pragma omp parallel for
    #endif
      for( Int i = 0; i < nnz; ++i ) {
//        A.val(i) = val[ i ]; //this along with BTF_A = A (without permuting)  works with the SolverFactory test matrix... - maybe the btf ordering and nd ordering turned out to be inverses for that matrix...
        BTF_A.val( inv_vals_order_ndbtfa_array(i) ) = val[ vals_perm_composition(i) ]; // BTF_A = A assigned during Symbolic (break_into_parts2) for this case; thus, identity map between A.val indices and BTF_A.val indices
        // for btf_nblks = 1 case, btf_tabs_offset set to 1 and possible to still apply nested dissection on BTF_A (i.e. A itself)
      }
//      BTF_A = A; //unnecessary - this equality was set during break_into_parts2, they point to the same data; for safety, should this simply be copied instead (i.e. deep copy the data)?
    } //end single block case
    else {
      std::cout << "Basker Factor error: Case for btf_nbkls = 0 is not implemented" << std::endl;
        //A.val(i) = val[ i ]; // may need to apply matching or nd order permutation...
      return BASKER_ERROR;
    }

    #ifdef BASKER_TIMER
    std::cout<< "Basker Factor: Time to permute and copy from input vals to new vals and blocks: " << copyperm_timer.seconds() << std::endl;
    #endif
    //end sfactor_copy2 replacement stuff


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


    #ifdef BASKER_TIMER
    Kokkos::Timer timer_sfactorcopy;
    double sfactorcopy_time = 0.0;
    #endif

    // sfactor_copy2 is now only responsible for the copy from BTF_A to 2D blocks
    err = sfactor_copy2();

    #ifdef BASKER_TIMER
    sfactorcopy_time += timer_sfactorcopy.seconds();
    std::cout << "Basker Factor sfactor_copy2 time: " << sfactorcopy_time << std::endl;
    #endif
    if(err == BASKER_ERROR)
    { return BASKER_ERROR; }

    #ifdef BASKER_TIMER
    Kokkos::Timer timer_factornotoken;
    double fnotoken_time = 0.0;
    #endif

    if(Options.incomplete == BASKER_FALSE)    
    {
      err = factor_notoken(0);
    }
    else
    {
      err = factor_inc_lvl(0);
    }

    #ifdef BASKER_TIMER
    fnotoken_time += timer_factornotoken.seconds();
    std::cout << "Basker factor_notoken total time: " << fnotoken_time << std::endl;
    #endif

    if(err == BASKER_ERROR)
    { 
      printf("ShyLUBasker factor_notoken/inc_lvl error returned\n");
      return BASKER_ERROR; 
    }

    if ( sym_gn != gn || sym_gm != gm ) {
      printf( "ShyLUBasker Factor error: Matrix dims at Symbolic and Factor stages do not agree - Symbolic reordered structure will not apply.\n");
      //exit(EXIT_FAILURE);
      return BASKER_ERROR; 
    }

    if(Options.verbose == BASKER_TRUE)
    { printf("Basker Factor Done \n"); }

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
  int Basker<Int,Entry,Exe_Space>::Factor_Inc(Int Options)
  {
    factor_inc_lvl(Options);

    return 0;
  }


  //Interface for solve.... only doing parallel solve right now.
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::SolveTest()
  {
    test_solve();
    return 0;
  }//end SolveTest


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Solve(Entry *b, Entry *x)
  {
    #ifdef BASKER_TIMER 
    Kokkos::Timer timer;
    double time = 0.0;
    #endif

    if(Options.verbose == BASKER_TRUE)
    {
      printf("Basker Solve Called \n");
    }

    solve_interface(x,b);

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
  int Basker<Int,Entry,Exe_Space>::Solve(Int nrhs, Entry *b, Entry *x)
  {
    #ifdef BASKER_TIMER 
    Kokkos::Timer timer;
    double time = 0.0;
    #endif

    if(Options.verbose == BASKER_TRUE)
    {
      printf("Basker MultiSolve Called \n");
    }

    solve_interface(nrhs,x,b);

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
  int Basker<Int, Entry, Exe_Space>::Solve(ENTRY_1DARRAY b, ENTRY_1DARRAY x)
  {
    printf("Basker: This solve call not implemented\n");
    return -1;
  }//Solve(ENTRY_1D, ENTRY_1D);


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Solve(Int nrhs, Entry *b, Entry *x, Int option)
  {    
    int err = 0;
    printf("Basker: This solve call not implemented\n");
    if(solve_flag == false) //never solved before
    {
      //err = malloc_init_solve(nrhs, x, b);
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
    //Need to test if power of 2.
    if((nthreads != 1) && (nthreads%2 != 0))
    {
      BASKER_ASSERT(0==1, "Basker SetThreads Assert: Number of thread error - not a multiple of 2");
      //Set default 1
      num_threads = 1;
      return BASKER_ERROR;
    }

    //Next test if Kokkos has that many threads!
    //This is a common mistake in mpi-based apps
    #ifdef KOKKOS_ENABLE_OPENMP
    int check_value = Kokkos::OpenMP::max_hardware_threads();
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

#endif //End ifndef
