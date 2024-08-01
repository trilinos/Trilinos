// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_ORDER_BTF_HPP
#define SHYLUBASKER_ORDER_BTF_HPP

#include "shylubasker_types.hpp"
#include "shylubasker_sswrapper.hpp"

#include "trilinos_btf_decl.h"

//#define BASKER_DEBUG_ORDER_BTF
//#define BASKER_TIMER

namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::find_btf_schedule
  (
   BASKER_MATRIX &M,
   Int           nblks,
   INT_1DARRAY  _btf_tabs
  )
  {
    //Find total work estimate
    Int total_work_estimate = 0;
    for(Int b = btf_tabs_offset; b < nblks; b++)
    {
     	total_work_estimate += btf_blk_work(b);
    }
   
    //Int break_size    = ceil((double)total_work_estimate*(
    //			      ((double)1/num_threads)));

    Int break_size    = ceil(  (double)total_work_estimate*(
			                        ((double)1/num_threads) + 
                              ((double)BASKER_BTF_IMBALANCE)) );

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Total schedul size: %ld \n", (long)total_work_estimate);
    printf("Break schedule size: %ld \n", (long)break_size);
    printf("Total num blks: %ld \n", (long)btf_nblks);
    #endif
    
    Int t_size      = 0;
    Int t_loc       = 0;
    btf_schedule(0) = btf_tabs_offset;
    //BASKER_BOOL  move_fwd = BASKER_TRUE; //NU

    for(Int b = btf_tabs_offset; b < btf_nblks; b++)
    {
      Int blk_work = btf_blk_work(b);
      t_size += blk_work;

      #ifdef BASKER_DEBUG_ORDER_BTF
      printf("t: %d blk: %d work: %d twork: %d \n", 
             t_loc,b, blk_work, t_size);
      #endif

      if(((t_size > break_size) && (t_loc < num_threads-1)) || (b == btf_nblks-1))
      {
        #ifdef BASKER_DEBUG_ORDER_BTF
        printf("New Schedule BLK: thread: %d size: %d blk: %d %d \n", t_loc, t_size, btf_schedule(t_loc), b); 
        #endif

        t_loc++;
        btf_schedule(t_loc) = b+1;
        t_size = blk_work;
      }
    }

  }//end find_btf_schedule()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::find_btf( BASKER_MATRIX &M )
  {
    Int nblks = 0;

    strong_component(M, nblks, order_btf_array, btf_tabs);

    btf_flag = BASKER_TRUE;

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("BTF nblks returned: %d \n", nblks);
    BASKER_ASSERT(nblks>1, "NOT ENOUGH BTF BLOCKS");
    #endif

    #ifdef BASKER_DEBUG_ORDER_BTF
    if(nblks<2)
    {
      printf("BTF did not find enough blks\n");
    }
    #endif

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("\n\nBTF tabs: \n");
    for(Int i=0; i < nblks+1; i++)
    {
      printf("%d, ", btf_tabs(i));
    }
    printf("\n");
    #endif

    permute_col(M, order_btf_array);
    permute_row(M, order_btf_array);

    break_into_parts(M, nblks, btf_tabs);

    btf_nblks = nblks;

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("------------BTF CUT: %d --------------\n", 
	  btf_tabs(btf_tabs_offset));
    #endif

    return 0;
  }//end find BTF


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::find_btf2
  (
   BASKER_MATRIX &M
  )
  {
    #ifdef BASKER_TIMER
    double order_time = 0.0;
    Kokkos::Timer timer_order;
    timer_order.reset();
    #endif
    Int nblks = 0;

    //================================================================
    //compute BTF
    strong_component(M, nblks, order_btf_array, btf_tabs);
    #ifdef BASKER_TIMER
    order_time = timer_order.seconds();
    std::cout << " >>> Basker order : strong comp time  : " << order_time << std::endl;
    #endif

    btf_nblks = nblks;
    btf_flag = BASKER_TRUE;

    if(Options.verbose == BASKER_TRUE)
    {
      printf("Basker: BTF nblks returned: %ld \n", (long)nblks);
    }
    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Basker: BTF nblks returned: %d \n", nblks);
    //BASKER_ASSERT(nblks>1, "NOT ENOUGH BTF BLOCKS");
    if(nblks<2)
    { printf("Basker: BTF did not find enough blks during strong_component\n"); }
    #endif


    if (Options.verbose == BASKER_TRUE)
    {
      printf("Basker: num_threads: %d \n", (int)num_threads);
      //printf("Basker BTF tabs: \n");
      //for(Int i=0; i < nblks+1; i++)
      //{
      //  printf(" btf_tabs[%d] = %ld, ", (int)i, (long)btf_tabs(i));
      //}
      //printf("\n");
    }//if verbose

    /*printf(" A = [\n" );
    for(Int j = 0; j < M.ncol; j++) {
      for(Int k = M.col_ptr[j]; k < M.col_ptr[j+1]; k++) {
        printf("%d %d %.16e\n", M.row_idx[k], j, M.val[k]);
      }
    }
    printf("];\n");*/

    MALLOC_INT_1DARRAY(vals_order_btf_array, M.nnz);
    //permute_col(M, order_btf_array);                                    //NDE: col-order M
    permute_col_store_valperms(M, order_btf_array, vals_order_btf_array); //NDE: col-order M & Track movement
    permute_row(M, order_btf_array);

    permute_inv(vals_perm_composition, vals_order_btf_array, M.nnz);


    //================================================================
    //compute (optionally MWM and) AMD on each block
    //NOTE: MWM is computed though both MWM and AMD will be recomputed in numerics
    //      AMD provides stats about LU
    MALLOC_INT_1DARRAY(order_blk_amd_array, M.ncol);
    init_value(order_blk_amd_array, M.ncol, (Int)0);
    MALLOC_INT_1DARRAY(order_blk_mwm_array, M.ncol);
    init_value(order_blk_mwm_array, M.ncol, (Int)0);

    MALLOC_INT_1DARRAY(order_blk_amd_inv, M.ncol);
    MALLOC_INT_1DARRAY(order_blk_mwm_inv, M.ncol);

    MALLOC_INT_1DARRAY(btf_blk_nnz, nblks+1);
    init_value(btf_blk_nnz, nblks+1, (Int)0);

    MALLOC_INT_1DARRAY(btf_blk_work, nblks+1);
    init_value(btf_blk_work, nblks+1, (Int)0);

    //=====================================================================
    //Find MWM + AMD blk ordering, also compute nnz and get workspace size
    //NOTE: ordering is computed for each of **ALL** the diagonal blocks
    //      (i.e., both A & C) since they are split to A & C after the 
    //      ordering is computed
    if(Options.verbose == BASKER_TRUE)
    {
      printf("Basker: block MWM+AMD(blk_matching = %d) \n", (int)Options.blk_matching);
    }
    #ifdef BASKER_TIMER
    timer_order.reset();
    #endif
    int blk_mwm_info = btf_blk_mwm_amd(M, order_blk_mwm_array, order_blk_amd_array, btf_blk_nnz, btf_blk_work);
    if (blk_mwm_info != BASKER_SUCCESS) {
      if(Options.verbose == BASKER_TRUE)
      {
        printf("Basker: error btf_blk_mwm_amd returned %d \n", (int)blk_mwm_info);
      }
      return blk_mwm_info;
    }
    #if 0 //debug
    printf( " >> debug: set order_blk_mwm/amd on global matrix to identity <<\n" );
    for (Int i = 0; i < (Int)M.nrow; i++) {
      order_blk_mwm_array(i) = i;
      order_blk_amd_array(i) = i;
    }
    #endif
    #ifdef BASKER_TIMER
    order_time = timer_order.seconds();
    std::cout << " >>> Basker order : Block MWM+AMD time: " << order_time << std::endl;
    #endif

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Basker blk_perm:\n");
    for(Int i = 0; i < M.ncol; i++)
    {
      printf("(%d,%d) ", i, order_blk_amd_array(i));
    }
    printf("\n");
    printf("id/blk_size/blk_nnz/work: \n");
    for(Int i = 0; i < nblks; i++)
    {
      printf("(%d, %d, %d, %d) ", i,
          btf_tabs(i+1)-btf_tabs(i), 
          btf_blk_nnz(i), btf_blk_work(i));
    }
    printf("\n");
    #endif

    //printMTX("A_BEFORE.mtx", M);
    //printVec("AMD.txt", order_blk_amd_array, M.ncol);

    #ifdef BASKER_TIMER
    timer_order.reset();
    #endif
    MALLOC_INT_1DARRAY(vals_order_blk_amd_array, M.nnz);
    if (Options.blk_matching == 0) // no blk_matching (TODO: should we add cardinality-matrching on each block?)
    {
      if (Options.verbose == BASKER_TRUE)
      {
        printf("Basker find BTF: apply BLK AMD\n");
      }
      //M.print_matrix("B.dat");

      // > apply AMD to cols & rows
      permute_col_store_valperms(M, order_blk_amd_array, vals_order_blk_amd_array); //NDE: col-order M & Track movement
      permute_row(M, order_blk_amd_array);
      //M.print_matrix("T.dat");

      #ifdef BASKER_TIMER
      order_time = timer_order.seconds();
      std::cout << " >>> Basker order : val-perm time     : " << order_time << std::endl;
      timer_order.reset();
      #endif

      permute_inv(vals_perm_composition, vals_order_blk_amd_array, M.nnz);
      #ifdef BASKER_TIMER
      order_time = timer_order.seconds();
      std::cout << " >>> Basker order : invert perm time  : " << order_time << std::endl;
      #endif
    } else {
      // reset matrix order and scale since they will be computed during numerical factorization
      if (Options.verbose == BASKER_TRUE)
      {
        printf("Basker find BTF: skip applying BLK AMD\n");
      }
      Entry one (1.0);
      for (Int i = 0; i < (Int)M.nrow; i++) {
        order_blk_mwm_array(i) = i;
        order_blk_amd_array(i) = i;
        scale_row_array(i) = one;
        scale_col_array(i) = one;
      }
    }

    //================================================================
    //Split the matrix M into blocks
    //NDE at this point, vals_perm_composition stores the permutation of the vals array; will be needed during break_into_parts
    break_into_parts2(M, nblks, btf_tabs);
    #ifdef BASKER_TIMER
    order_time = timer_order.seconds();
    std::cout << " >>> Basker order : partition time    : " << order_time << std::endl;
    timer_order.reset();
    #endif
    //M.print_matrix("M.dat");
    //BTF_A.print_matrix("A.dat");
    //BTF_D.print_matrix("D.dat");
    //BTF_E.print_matrix("E.dat");
    //BTF_C.print_matrix("C.dat");
    //BTF_B.print_matrix("B.dat");


    //================================================================
    //find schedule
    find_btf_schedule(M, nblks, btf_tabs);
    #ifdef BASKER_TIMER
    order_time = timer_order.seconds();
    std::cout << " >>> Basker order : schedule time     : " << order_time << std::endl;
    #endif

    if (Options.verbose == BASKER_TRUE)
    {
      printf("Basker BTF Cut: %ld \n", (long)btf_tabs(btf_tabs_offset));
    }

#ifdef BASKER_DEBUG_ORDER_BTF
    printf("------------BTF CUT: %d --------------\n",  btf_tabs(btf_tabs_offset));
#endif

    return BASKER_SUCCESS;
  }//end find_btf2


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry,Exe_Space>::break_into_parts
  (
   BASKER_MATRIX &M,
   Int           nblks,
   INT_1DARRAY  _btf_tabs
  )
  {
    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("break_into_parts called \n");
    printf("nblks: %d \n", nblks);
    #endif

    Options.btf = BASKER_TRUE;

    //Alg.  
    // A -> [BTF_A  BTF_B] 
    //      [0      BTF_C]
    //1. Run backward through the btf_tabs to find size C
    //2. Form A,B,C based on size in 1.

    //Step 1.
    Int t_size            = 0;
    Int scol              = M.ncol;
    Int blk_idx           = nblks;
    BASKER_BOOL  move_fwd = BASKER_TRUE;
    while(move_fwd==BASKER_TRUE)
    {

      Int blk_size = _btf_tabs(blk_idx) - _btf_tabs(blk_idx-1);

    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("move_fwd loop \n");
      BASKER_ASSERT(blk_idx>=0, "btf blk idx off");
      BASKER_ASSERT(blk_size>0, "btf blk size wrong");
      printf("blk_idx: %d blk_size: %d \n", 
          blk_idx, blk_size);
      std::cout << blk_size << std::endl;
    #endif


      if((blk_size < Options.btf_large) &&
          ((((double)t_size+blk_size)/(double)M.ncol) < Options.btf_max_percent))
      {
    #ifdef BASKER_DEBUG_ORDER_BTF
        printf("first choice \n");
        printf("blksize test: %d %d %d \n",
            blk_size, Options.btf_large, 
            BASKER_BTF_LARGE);
        printf("blkpercent test: %f %f %f \n", 
            ((double)t_size+blk_size)/(double)M.ncol, 
            Options.btf_max_percent, 
            (double) BASKER_BTF_MAX_PERCENT);
    #endif

        t_size  = t_size+blk_size;
        blk_idx = blk_idx-1;
        scol    = _btf_tabs[blk_idx];
      }
      else
      {
        //printf("second choice \n");
        //#ifdef BASKER_DEBUG_ORDER_BTF
        printf("Cut: blk_size: %ld percent: %lf \n",
            (long)blk_size, ((double)t_size+blk_size)/(double)M.ncol);

        if((((double)t_size+blk_size)/(double)M.ncol) == 1.0)
        {
          blk_idx = 0;
          t_size = t_size + blk_size;
          scol   = _btf_tabs[blk_idx];
        }

        //#endif
        move_fwd = BASKER_FALSE;
      }
    }//end while(move_fwd)

  #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Done finding BTF cut.  Cut size: %d scol: %d \n",
        t_size, scol);
    //BASKER_ASSERT(t_size > 0, "BTF CUT SIZE NOT BIG ENOUGH\n");
    BASKER_ASSERT((scol >= 0) && (scol < M.ncol), "SCOL\n");
  #endif

    //Comeback and change
    btf_tabs_offset = blk_idx;

    //Step 2. Move into Blocks
    if(btf_tabs_offset != 0)
    {
      //--Move A into BTF_A;
      BTF_A.set_shape(0, scol, 0, scol);
      BTF_A.nnz = M.col_ptr(scol);

    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("Init BTF_A. ncol: %d nnz: %d \n",
          scol, BTF_A.nnz);
    #endif

      if(BTF_A.v_fill == BASKER_FALSE)
      {
        BASKER_ASSERT(BTF_A.ncol >= 0, "BTF_A, col_ptr");
        MALLOC_INT_1DARRAY(BTF_A.col_ptr, BTF_A.ncol+1);
        BASKER_ASSERT(BTF_A.nnz > 0, "BTF_A, nnz");
        MALLOC_INT_1DARRAY(BTF_A.row_idx, BTF_A.nnz);
        MALLOC_ENTRY_1DARRAY(BTF_A.val, BTF_A.nnz);
        BTF_A.fill();
      }

      Int annz = 0;
      for(Int k = 0; k < scol; ++k)
      {
      #ifdef BASKER_DEBUG_ORDER_BTF
        printf("copy column: %d into A_BTF, [%d %d] \n", 
            k, M.col_ptr(k), M.col_ptr(k+1));
      #endif

        for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
        {
          //printf("annz: %d i: %d \n", annz, i);
          BTF_A.row_idx(annz) = M.row_idx(i);
          BTF_A.val(annz)     = M.val(i);
          annz++;
        }

        BTF_A.col_ptr(k+1) = annz;
      }

    }//no A

    //Fill in B and C at the same time
    INT_1DARRAY cws;
    BASKER_ASSERT((M.ncol-scol+1) > 0, "BTF_SIZE MALLOC");
    MALLOC_INT_1DARRAY(cws, M.ncol-scol+1);
    init_value(cws, M.ncol-scol+1, (Int)M.ncol);
    BTF_B.set_shape(0 , scol,
        scol, M.ncol-scol);
    BTF_C.set_shape(scol, M.ncol-scol,
        scol, M.ncol-scol);

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Set Shape BTF_B: %d %d %d %d \n",
        BTF_B.srow, BTF_B.nrow,
        BTF_B.scol, BTF_B.ncol);
    printf("Set Shape BTF_C: %d %d %d %d \n",
        BTF_C.srow, BTF_C.nrow,
        BTF_C.scol, BTF_C.nrow);
    #endif

    //Scan and find nnz
    //We can do this much better!!!!
    Int bnnz = 0;
    Int cnnz = 0;
    for(Int k = scol; k < M.ncol; ++k)
    {
    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("Scanning nnz, k: %d \n", k);
    #endif

      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
      {
        if(M.row_idx(i) < scol)
        {
        #ifdef BASKER_DEBUG_ORDER_BTF
          printf("Adding nnz to Upper, %d %d \n",
              scol, M.row_idx(i));
        #endif
          bnnz++;
        }
        else
        {
        #ifdef BASKER_DEBUG_ORDER_BTF
          printf("Adding nnz to Lower, %d %d \n",
              scol, M.row_idx(i));
        #endif
          cnnz++;
        }
      }//over all nnz in k
    }//over all k

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("BTF_B nnz: %d \n", bnnz);
    printf("BTF_C nnz: %d \n", cnnz);
    #endif

    BTF_B.nnz = bnnz;
    BTF_C.nnz = cnnz;

    //Malloc need space
    if((BTF_B.v_fill == BASKER_FALSE) &&
        (BTF_B.nnz > 0))
    {
      BASKER_ASSERT(BTF_B.ncol >= 0, "BTF_B ncol");
      MALLOC_INT_1DARRAY(BTF_B.col_ptr, BTF_B.ncol+1);
      BASKER_ASSERT(BTF_B.nnz > 0, "BTF_B.nnz");
      MALLOC_INT_1DARRAY(BTF_B.row_idx, BTF_B.nnz);
      MALLOC_ENTRY_1DARRAY(BTF_B.val, BTF_B.nnz);
      BTF_B.fill();
    }
    if(BTF_C.v_fill == BASKER_FALSE)
    {
      BASKER_ASSERT(BTF_C.ncol >= 0, "BTF_C.ncol");
      MALLOC_INT_1DARRAY(BTF_C.col_ptr, BTF_C.ncol+1);
      BASKER_ASSERT(BTF_C.nnz > 0, "BTF_C.nnz");
      MALLOC_INT_1DARRAY(BTF_C.row_idx, BTF_C.nnz);
      MALLOC_ENTRY_1DARRAY(BTF_C.val, BTF_C.nnz);
      BTF_C.fill();
    }

    //scan again (Very bad!!!)
    bnnz = 0;
    cnnz = 0;
    for(Int k = scol; k < M.ncol; ++k)
    {
    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("Scanning nnz, k: %d \n", k);
    #endif

      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
      {
        if(M.row_idx(i) < scol)
        {
        #ifdef BASKER_DEBUG_ORDER_BTF
          printf("Adding nnz to Upper, %d %d \n",
              scol, M.row_idx[i]);
        #endif

          BASKER_ASSERT(BTF_B.nnz > 0, "BTF B uninit");
          //BTF_B.row_idx[bnnz] = M.row_idx[i];
          //Note: do not offset because B srow = 0
          BTF_B.row_idx(bnnz) = M.row_idx(i);
          BTF_B.val(bnnz)     = M.val(i);
          bnnz++;
        }
        else
        {
        #ifdef BASKER_DEBUG_ORDER_BTF
          printf("Adding nnz Lower,k: %d  %d %d %f \n",
              k, scol, M.row_idx[i], 
              M.val(i));
        #endif
          //BTF_C.row_idx[cnnz] = M.row_idx[i];
          BTF_C.row_idx(cnnz) = M.row_idx(i)-scol;
          BTF_C.val(cnnz)     = M.val(i);
          cnnz++;
        }
      }//over all nnz in k
      if(BTF_B.nnz > 0)
      {
        BTF_B.col_ptr(k-scol+1) = bnnz;
      }
      BTF_C.col_ptr(k-scol+1) = cnnz;
    }//over all k

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("After BTF_B nnz: %d \n", bnnz);
    printf("After BTF_C nnz: %d \n", cnnz);
    #endif

    return 0;
  }//end break_into_parts


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry,Exe_Space>::break_into_parts2
  (
   BASKER_MATRIX &M,
   Int           nblks,
   INT_1DARRAY  _btf_tabs
  )
  {
  #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Basker: break_into_parts2 called \n");
    printf("nblks: %d \n", nblks);
  #endif

    Options.btf = BASKER_TRUE;

    //Alg. 
    // Old: 
    // A -> [BTF_A  BTF_B] 
    //      [0      BTF_C]
    // New: 
    // A -> [BTF_D    BTF_E     ]  (BTF_E stores interface to both BTF_A and BTF_C)
    //      [0      BTF_A  BTF_B]  (BTF_A has just one block)
    //      [0      0      BTF_C]
    //1. Run backward through the btf_tabs to find size of C
    //2. Form A,B,C ,and D & E,  based on step 1.

    //Short circuit, 
    //If nblks  == 1, than only BTF_A exists
    // NDE: In this case, vals_block_map_perm_pair is not allocated nor used - A is assigned to BTF_A directly
    if(nblks == 1)
    {
    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("Basker: break_into_parts2 - short circuit for single block case\n");
    #endif
      #if !defined (HAVE_SHYLU_NODEBASKER_METIS) & !defined(HAVE_SHYLU_NODEBASKER_SCOTCH)
      if (Options.run_nd_on_leaves == BASKER_TRUE) {
        if(Options.verbose == BASKER_TRUE) {
          printf("Basker: turning off ND-on-leaves option since no METIS nor SCOTCH (hence sequential)\n");
        }
        Options.run_nd_on_leaves = BASKER_FALSE;
      }
      #endif
      if (Options.replace_zero_pivot == BASKER_TRUE) {
        if(Options.verbose == BASKER_TRUE) {
          printf("Basker: turning off replace-zero-pivot option since one block (to identify singular matrix)\n");
        }
        Options.replace_zero_pivot = BASKER_FALSE;
      }
    }

    //Short circuit for incomplete
    if(Options.incomplete == BASKER_TRUE)
    {
    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("Basker: break_into_parts2 - short ciruit incomplete\n");
    #endif
      BTF_A = A;
      btf_nblks = 1;
      btf_tabs_offset = 1;
      return 0;
    }

    Int scol_top          = 0;      // starting column of the BTF_A bloc
    Int scol              = M.ncol; // starting column of the BTF_C blocks (end of BTF_A block)
    Int blk_idx           = nblks;  // start at lower right corner block, move left and up the diagonal
    #if !defined (HAVE_SHYLU_NODEBASKER_METIS) & !defined(HAVE_SHYLU_NODEBASKER_SCOTCH)
    // use Metis on a large block even with one thread
    if (num_threads == 1) {
      // Short circuit for single thread = no big block A
      scol = 0;
      blk_idx = 0;
      if(Options.verbose == BASKER_TRUE) {
        printf("Basker: short-circuit for one thread\n");
      }
    } else
    #endif
    {
      //Step 1.
      //Find total work estimate
      double total_work_estimate = 0;
      for(Int b = 0; b < nblks; b++) //nblks is input; determined during btf ordering - total BTF_A blocks AND BTF_C blocks
      {
        Int blk_size = _btf_tabs(b+1) - _btf_tabs(b);
        Int wrk_size = btf_blk_work(b); //determined prior, during btf ordering
        total_work_estimate += (blk_size > wrk_size ? blk_size : wrk_size);
      }
      //Set a class variable to use later
      btf_total_work = total_work_estimate;
      //printf("num_threads: %d epsilon: %f \n",
      //	   num_threads, 
      //	   ((double)1/num_threads) +
      //	   ((double)BASKER_BTF_IMBALANCE));
      #if 0 // forcing to have the big A bloock for debug
      double break_work_size = 0.0;
      //double break_block_size = 0.0;
      double break_block_size = 0.0;
      printf( " > debug: break_size = %f, %f\n",break_work_size,break_block_size );
      #else
      // A block if it is larger than work esitimate assigned to one thread
      double break_fact = 0.7;
      double break_work_size = ceil(total_work_estimate*(break_fact * ((double)1.0/num_threads) + ((double)BASKER_BTF_IMBALANCE)));
      double break_block_size = 20 * num_threads; //0;
      #endif
      if(Options.verbose == BASKER_TRUE) {
        printf("Basker: Break size for workspace and size: %d and %d with %d threads (total work estimate = %f)\n",
                (int)break_work_size, (int)break_block_size, (int)num_threads, total_work_estimate);
      }

      scol_top          = 0;      // starting column of the BTF_A bloc
      scol              = M.ncol; // starting column of the BTF_C blocks (end of BTF_A block)
      blk_idx           = nblks;  // start at lower right corner block, move left and up the diagonal
                                  // note: blk_idx index acts as block id + 1; btf_tabs(nblks) is likely the end column number

      Int t_size            = 0;      //total size of cols from 'small' blocks in BTF_C: matrix ncols - t_size = BTF_A ncols
      BASKER_BOOL  move_fwd = BASKER_TRUE;
      while(move_fwd==BASKER_TRUE)
      {
        Int blk_work = btf_blk_work(blk_idx-1);
        // subtract the bounding column ids to determine size of the (square) block
        Int blk_size = _btf_tabs(blk_idx) - _btf_tabs(blk_idx-1);

        #ifdef BASKER_DEBUG_ORDER_BTF
        printf(" \n move_fwd loop \n");
        BASKER_ASSERT(blk_idx>=0, "btf blk idx off");
        BASKER_ASSERT(blk_work>=0, "btk_work wrong");
        BASKER_ASSERT(blk_size>0, "btf blk size wrong");
        printf("blk_idx: %d blk_work: %d break_size: %d \n",
            blk_idx, blk_work, break_size);
        #endif

        //Should be end
        //if(((blk_work < break_size) ||
        //  (blk_size < BASKER_BTF_SMALL)) &&
        //  (blk_idx > 1))

        //Continue to be in btf
        if( (Options.use_sequential_diag_facto || (blk_work <= break_work_size || blk_size < break_block_size)) && (blk_idx > 1) )
        {
          #ifdef BASKER_DEBUG_ORDER_BTF
           printf("Basker(blk_idx=%d, blk_size=%d, blk_work=%d, break_size=%d): continue with fine structure btf blocks\n",
                   (int)blk_idx,(int)blk_size,(int)blk_work,(int)break_size);
          #endif
          t_size  = t_size+blk_size;
          blk_idx = blk_idx-1;
          scol    = _btf_tabs[blk_idx];
        }
        //break due to size i.e. entered non-trivial large BTF_A block
        else if( blk_work > break_work_size && blk_size >= break_block_size)
        {
          if(Options.verbose == BASKER_TRUE) {
            printf("Basker: blk=%d break due to size (work: %d > %d, size: %d > %d)\n",(int)blk_idx-1, (int)blk_work,(int)break_work_size, (int)blk_size,(int)break_block_size);
          }
          move_fwd = BASKER_FALSE;
        }
        //break due to end i.e. no 'large' BTF_A block for ND; only fine BTF structure
        else if(blk_idx == 1)
        {
          #ifdef BASKER_DEBUG_ORDER_BTF
           printf("Basker: last block\n");
          #endif
          //printf("break last blk\n");
          blk_idx = 0;
          t_size = t_size + blk_size;
          scol = _btf_tabs[blk_idx];	
          move_fwd = BASKER_FALSE;
        }
        //should not be called
        else
        {
          if(Options.verbose == BASKER_TRUE) {
            printf("Basker: invalid blk=%d (work: %d vs %d, size: %d vs %d)\n",(int)blk_idx-1, (int)blk_work,(int)break_work_size, (int)blk_size,(int)break_block_size);
          }
          BASKER_ASSERT(1==0, "btf order break");
          move_fwd = BASKER_FALSE;
        }
      }//end while(move_fwd)
    }

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Basker: Done finding BTF2 cut.  Cut size: %d scol: %d \n", t_size, scol);
    printf("Basker: Done finding BTF2 cut. blk_idx: %d \n", blk_idx);
    //BASKER_ASSERT(t_size > 0, "BTF CUT SIZE NOT BIG ENOUGH\n");
    BASKER_ASSERT((scol >= 0) && (scol <= M.ncol), "SCOL\n");
    #endif

    //Comeback and change
    // btf_tabs_offset is offset to id of first block after diag blocks in BTF_A
    // (i.e., the first block in BTF_C)
    btf_tabs_offset = blk_idx; // 

    //Step 2. Move into Blocks 
    MALLOC_INT_1DARRAY_PAIRS(vals_block_map_perm_pair, M.nnz); //this will store and map A.val indices to val indices of BTF_A, BTF_B, and BTF_C 
    if(btf_tabs_offset != 0) // if the id of the first block in BTF_C != 0
    {                        // -> we have BTF_A, and maybe blocks in BTF_D
#if defined(BASKER_SPLIT_A)
      #if 0
      // debuging with A merging all the top blocks
      printf( " >> debug: remove D/E blocks <<\n" );
      btf_top_tabs_offset = 0;
      btf_top_nblks = 0;

      scol_top = 0;
      #else
      btf_top_tabs_offset = btf_tabs_offset-1;  // starting block ID of BTF_A block (for now, there is just one big block in A)
      btf_top_nblks = btf_top_tabs_offset;      // number of blocks in BTF_D

      scol_top = _btf_tabs[btf_top_tabs_offset]; // the first column index of A
      #endif

      // extract BTF_D and BTF_E
      Int dnnz = M.col_ptr(scol_top);
      Int ennz = 0;
      for(Int k = scol_top; k < M.ncol; ++k) {
        for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i) {
          if(M.row_idx(i) < scol_top) {
            ennz ++;
          }
        }
      }
      BTF_D.set_shape(0, scol_top, 0,        scol_top);
      BTF_E.set_shape(0, scol_top, scol_top, M.ncol-scol_top);
      BTF_D.nnz = dnnz;
      BTF_E.nnz = ennz;
      if(BTF_D.v_fill == BASKER_FALSE)
      {
        MALLOC_INT_1DARRAY(BTF_D.col_ptr, BTF_D.ncol+1);
        if (BTF_D.nnz > 0) {
          MALLOC_INT_1DARRAY(BTF_D.row_idx, BTF_D.nnz);
          MALLOC_ENTRY_1DARRAY(BTF_D.val, BTF_D.nnz);
        }
        BTF_D.fill();
      }
      if(BTF_E.v_fill == BASKER_FALSE)
      {
        MALLOC_INT_1DARRAY(BTF_E.col_ptr, BTF_E.ncol+1);
        if (BTF_E.nnz > 0) {
          MALLOC_INT_1DARRAY(BTF_E.row_idx, BTF_E.nnz);
          MALLOC_ENTRY_1DARRAY(BTF_E.val, BTF_E.nnz);
        }
        BTF_E.fill();
      }
      if (BTF_D.nnz > 0) {
        dnnz = 0;
        for(Int k = 0; k < scol_top; ++k) {
          for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i) {
            BTF_D.row_idx(dnnz) = M.row_idx(i);
            BTF_D.val(dnnz)     = M.val(i);
            vals_block_map_perm_pair(i) = std::pair<Int,Int>(-1, dnnz);

            dnnz ++;
          }
          BTF_D.col_ptr(k+1) = dnnz;
        }
      }
      ennz = 0;
      if (BTF_E.nnz > 0) {
        for(Int k = scol_top; k < M.ncol; ++k) {
          for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i) {
            if(M.row_idx(i) < scol_top) {
              BTF_E.row_idx(ennz) = M.row_idx(i);
              BTF_E.val(ennz)     = M.val(i);
              vals_block_map_perm_pair(i) = std::pair<Int,Int>(-2, ennz);

              ennz++;
            }
          }
          BTF_E.col_ptr(k-scol_top+1) = ennz;
        }
      }
#endif

      //--Move A into BTF_A;
      Int annz = 0;
      if (scol_top > 0) {
        for(Int k = scol_top; k < scol; ++k) {
          for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i) {
            if (M.row_idx(i) >= scol_top) {
              annz ++;
            }
          }
        }
      } else {
        annz = M.col_ptr(scol);
      }
      //BTF_A.set_shape(scol_top, scol-scol_top, scol_top, scol-scol_top);
      BTF_A.set_shape(0, scol-scol_top, 0, scol-scol_top);
      BTF_A.nnz = annz;
      if(Options.verbose == BASKER_TRUE) {
        printf( " + scol = %d, scol_top = %d, btf_tabs = %d, btf_top_tabs = %d, annz = %d\n",
                (int)scol, (int)scol_top, (int)btf_tabs_offset, (int)btf_top_tabs_offset, (int)annz );
        printf( " +  BTF_A.srow = %d, BTF_A.scol = %d, BTF_A.nrow = %d, BTF_A.ncol = %d\n",
                (int)BTF_A.srow, (int)BTF_A.scol, (int)BTF_A.nrow, (int)BTF_A.ncol );
      }

    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("Basker Init BTF_A. ncol: %d nnz: %d \n", scol, BTF_A.nnz);
    #endif

      if(BTF_A.v_fill == BASKER_FALSE)
      {
        BASKER_ASSERT(BTF_A.ncol >= 0, "BTF_A, col_ptr");
        MALLOC_INT_1DARRAY(BTF_A.col_ptr, BTF_A.ncol+1);
        if(BTF_A.nnz > 0) {
          BASKER_ASSERT(BTF_A.nnz > 0, "BTF_A, nnz");
          MALLOC_INT_1DARRAY(BTF_A.row_idx, BTF_A.nnz);
          MALLOC_ENTRY_1DARRAY(BTF_A.val, BTF_A.nnz);
        }
        BTF_A.fill();
      }

      if(BTF_A.nnz > 0) {
        annz = 0;
        for(Int k = scol_top; k < scol; ++k)
        {
        #ifdef BASKER_DEBUG_ORDER_BTF
          printf("Basker copy column: %d into A_BTF, [%d %d] \n", 
              k, M.col_ptr(k), M.col_ptr(k+1));
        #endif

          for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
          {
            if (M.row_idx(i) >= scol_top) {
              BTF_A.row_idx(annz) = M.row_idx(i) - scol_top;
              BTF_A.val(annz)     = M.val(i);  //NDE: Track movement of vals (lin_ind of row,col) here

              vals_block_map_perm_pair(i) = std::pair<Int,Int>(0,annz);

              annz++;
            }
          }
          BTF_A.col_ptr(k-scol_top+1) = annz;
        }
      }
    }//no A

    if(Options.verbose == BASKER_TRUE) {
      printf( "\n > btf_tabs_offset = %d, btf_top_tabs_offset = %d\n", (int)btf_tabs_offset, (int)btf_top_tabs_offset );
      for (blk_idx = 0; blk_idx < btf_top_tabs_offset; blk_idx++) printf( " x %d: %d (%d)\n", (int)blk_idx, (int)(btf_tabs[blk_idx+1]-btf_tabs[blk_idx]),(int)btf_blk_work(blk_idx) );
      for (blk_idx = btf_top_tabs_offset; blk_idx < btf_tabs_offset; blk_idx++) printf( " + %d: %d (%d)\n", (int)blk_idx, (int)(_btf_tabs[blk_idx+1]-_btf_tabs[blk_idx]),(int)btf_blk_work(blk_idx) );
      for (blk_idx = btf_tabs_offset; blk_idx < nblks; blk_idx++) printf( " - %d: %d (%d)\n", (int)blk_idx, (int)(_btf_tabs[blk_idx+1]-_btf_tabs[blk_idx]),(int)btf_blk_work(blk_idx) );
      printf( "\n" );
    }

    //Fill in B and C at the same time
    INT_1DARRAY cws;
    BASKER_ASSERT((M.ncol-scol+1) > 0, "BTF_SIZE MALLOC");
    MALLOC_INT_1DARRAY(cws, M.ncol-scol+1);
    init_value(cws, M.ncol-scol+1, (Int)M.ncol);
    BTF_B.set_shape(scol_top,    scol-scol_top,  scol, M.ncol-scol);
    BTF_C.set_shape(scol,        M.ncol-scol,    scol, M.ncol-scol);

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf( " scol_top = %d, scol = %d, ncol = %d\n",scol_top, scol, M.ncol );
    printf("Basker Set Shape BTF_B: %d:%d, %d:%d \n",
        BTF_B.srow, BTF_B.nrow,
        BTF_B.scol, BTF_B.ncol);
    printf("Basker Set Shape BTF_C: %d:%d, %d:%d \n",
        BTF_C.srow, BTF_C.nrow,
        BTF_C.scol, BTF_C.nrow);
    #endif

    //Added check - this is allowed, single block case where BTF_A <- A; if this occurs simply skip out of this routine
    if((BTF_C.nrow == 0)||(BTF_C.ncol == 0))
    {
      //printf("ShyLUBasker: either BTF_C number of rows or columns is 0\n");
      return 0;
    }

    //Scan and find nnz
    //We can do this much better!!!!
    Int bnnz = 0;
    Int cnnz = 0;
    for(Int k = scol; k < M.ncol; ++k)
    {
    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("Basker Scanning nnz, k: %d \n", k);
    #endif

      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
      {
        if(M.row_idx(i) >= scol_top) {
          if(M.row_idx(i) < scol)
          {
          #ifdef BASKER_DEBUG_ORDER_BTF
            printf("Basker Adding nnz to Upper, %d %d \n", scol, M.row_idx(i));
          #endif
            bnnz++;
          }
          else
          {
          #ifdef BASKER_DEBUG_ORDER_BTF
            printf("Basker Adding nnz to Lower, %d %d \n", scol, M.row_idx(i));
          #endif
            cnnz++;
          }
        }
      }//over all nnz in k
    }//over all k

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("BTF_B nnz: %d \n", bnnz);
    printf("BTF_C nnz: %d \n", cnnz);
    #endif

    BTF_B.nnz = bnnz;
    BTF_C.nnz = cnnz;

    //Malloc need space
    if(BTF_B.v_fill == BASKER_FALSE)
    {
      BASKER_ASSERT(BTF_B.ncol >= 0, "BTF_B ncol");
      MALLOC_INT_1DARRAY(BTF_B.col_ptr, BTF_B.ncol+1);
      if(BTF_B.nnz > 0) {
        BASKER_ASSERT(BTF_B.nnz > 0, "BTF_B.nnz");
        MALLOC_INT_1DARRAY(BTF_B.row_idx, BTF_B.nnz);
        MALLOC_ENTRY_1DARRAY(BTF_B.val, BTF_B.nnz);
      }
      BTF_B.fill();
    }
    if(BTF_C.v_fill == BASKER_FALSE)
    {
      BASKER_ASSERT(BTF_C.ncol >= 0, "BTF_C.ncol");
      MALLOC_INT_1DARRAY(BTF_C.col_ptr, BTF_C.ncol+1);
      if(BTF_C.nnz > 0) {
        BASKER_ASSERT(BTF_C.nnz > 0, "BTF_C.nnz");
        MALLOC_INT_1DARRAY(BTF_C.row_idx, BTF_C.nnz);
        MALLOC_ENTRY_1DARRAY(BTF_C.val, BTF_C.nnz);
      }
      BTF_C.fill();
    }

    //scan again (Very bad!!!)
    bnnz = 0;
    cnnz = 0;
    for(Int k = scol; k < M.ncol; ++k)
    {
    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("Scanning nnz, k: %d \n", k);
    #endif

      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
      {
        if(M.row_idx(i) >= scol_top) {
          if(M.row_idx(i) < scol)
          {
          #ifdef BASKER_DEBUG_ORDER_BTF
            printf("Adding nnz to Upper, %d %d \n", scol, M.row_idx[i]);
          #endif

            BASKER_ASSERT(BTF_B.nnz > 0, "BTF B uninit");
            //BTF_B.row_idx[bnnz] = M.row_idx[i];
            //Note: do not offset because B srow = 0
            BTF_B.row_idx(bnnz) = M.row_idx(i) - scol_top;
            BTF_B.val(bnnz)     = M.val(i);

            vals_block_map_perm_pair(i) = std::pair<Int,Int>(1,bnnz);

            bnnz++;
          }
          else
          {
          #ifdef BASKER_DEBUG_ORDER_BTF
            printf("Adding nnz Lower,k: %d  %d %d %f \n",
                k, scol, M.row_idx[i], 
                M.val(i));
          #endif
            BTF_C.row_idx(cnnz) = M.row_idx(i)-scol;
            BTF_C.val(cnnz)     = M.val(i);

            vals_block_map_perm_pair(i) = std::pair<Int,Int>(2,cnnz);

            cnnz++;
          }
        }
      }//over all nnz in k
      //if(BTF_B.ncol > 0)
      {
        BTF_B.col_ptr(k-scol+1) = bnnz;
      }
      //if(BTF_C.ncol > 0)
      {
        BTF_C.col_ptr(k-scol+1) = cnnz;
      }
    }//over all k

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("After BTF_B nnz: %d \n", bnnz);
    printf("After BTF_C nnz: %d \n", cnnz);
    #endif

    return 0;
  }//end break_into_parts2 (based on imbalance)


  template <class Int,class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::strong_component
  (
   BASKER_MATRIX &M,
   Int           &nblks,
   INT_1DARRAY   &perm,
   INT_1DARRAY   &CC
  )
  {
    //printf("===Basker this strong comp called====");

    typedef long int   l_int;
    
    INT_1DARRAY perm_in;
    MALLOC_INT_1DARRAY(perm_in, M.ncol);
    MALLOC_INT_1DARRAY(perm, M.ncol);
    //JDB:Note, this needs to be changed just fixed for int/long
    MALLOC_INT_1DARRAY(CC, M.ncol+1);

    for(l_int i = 0; i < M.ncol; i++)
    {
      perm_in(i) = i;
    }

    if(Options.incomplete == BASKER_TRUE)
    {
      for(Int i = 0; i < M.ncol; i++)
      {
        perm(i) = i;
      }
      nblks = 1;
      CC(0) = 0;
      CC(1) = M.ncol;
      return 0;
    }

    BaskerSSWrapper<Int>::my_strong_component(M.ncol,
			&(M.col_ptr(0)),
			&(M.row_idx(0)),
			nblks,
			&(perm(0)),
			&(perm_in(0)), 
			&(CC(0)));

    if (Options.min_block_size > 0 && nblks > 1) {
      for (Int blk = 0; blk < nblks-1; blk++) {
        Int blk_size = CC(blk+1) - CC(blk);
        if (blk_size < Options.min_block_size) {
          if (blk == 0) {
            // merge this block to the next block
            for (Int I = blk+1; I < nblks; I++) {
              CC(I) = CC(I+1);
            }
            // the next block is blk-th block
            blk --;   // check this block again
            nblks --; // reducd # of blocks
          } else if (blk == nblks-1) {
            // merge this block to the prev block
            for (Int I = blk; I < nblks; I++) {
              CC(I) = CC(I+1);
            }
            // the prev block is (blk-1)th block,
            // the next block is blk-th block
            blk --;   // check this block again
            nblks --; // reducd # of blocks
            break;    // done
          } else {
            Int prev_size = CC(blk)-CC(blk-1);
            Int next_size = CC(blk+1)-CC(blk);
            if (prev_size < next_size) {
              // merge this block to the prev block
              for (Int I = blk; I <= nblks; I++) {
                CC(I) = CC(I-1);
              }
            } else {
              // merge this block to the next block
              for (Int I = blk+1; I < nblks; I++) {
                CC(I) = CC(I+1);
              }
            }
            blk --;   // check this block again
            nblks --; // reducd # of block
          }
        }
      }
    }

    #ifdef BASKER_DEBUG_ORDER_BTF
    FILE *fp;
    fp = fopen("btf.txt", "w");
    for(Int i = 0; i < M.ncol; i++)
    {
      fprintf(fp, "%ld \n", (long)perm(i));
    }
    fclose(fp);
    #endif
    
    return 0;
  }//end strong_component <long int>

}//end namespace BaskerNS
#undef BASKER_TIMER
#endif
