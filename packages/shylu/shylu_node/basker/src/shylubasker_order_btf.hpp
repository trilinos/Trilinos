#ifndef SHYLUBASKER_ORDER_BTF_HPP
#define SHYLUBASKER_ORDER_BTF_HPP

#include "shylubasker_types.hpp"
#include "shylubasker_sswrapper.hpp"

#include "trilinos_btf_decl.h"

//#define BASKER_DEBUG_ORDER_BTF

namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::find_btf_schedule
  (
   BASKER_MATRIX &M,
   Int           nblks,
   INT_1DARRAY   btf_tabs
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
    
    Int t_size            = 0;
    Int t_loc             = 0;
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

    strong_component(M,nblks,order_btf_array,btf_tabs);

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
    Int nblks = 0;

    strong_component(M,nblks,order_btf_array,btf_tabs);

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


    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Basker: num_threads: %d \n", num_threads);
    printf("\n\nBasker: BTF tabs: \n");
    for(Int i=0; i < nblks+1; i++)
    {
      printf("%d, ", btf_tabs(i));
    }
    printf("\n");
    #endif

    if (Options.verbose == BASKER_TRUE)
    {
      //Print first 10 and last 3 
      printf("Basker BTF tabs (first 10): \n");
      if (nblks < 10)
      {
        for(Int i=0; i < nblks+1; i++)
        {
          printf("%ld, ", (long)btf_tabs(i));
        }
      }
      else
      {
        printf("%ld, %ld, %ld, ...., ",
            (long)btf_tabs(0), (long)btf_tabs(1), (long)btf_tabs(2));
        printf("%ld, %ld, %ld",
            (long)btf_tabs(nblks-3), 
            (long)btf_tabs(nblks-2),
            (long)btf_tabs(nblks-1));
      }
      printf("\n");
    }//if verbose

    MALLOC_INT_1DARRAY(vals_order_btf_array, M.nnz);
    //permute_col(M, order_btf_array); //NDE: Track movement of vals (lin_ind of row,col) here
    permute_col_store_valperms(M, order_btf_array, vals_order_btf_array); //NDE: Track movement of vals (lin_ind of row,col) here
    permute_row(M, order_btf_array);

    permute_inv(vals_perm_composition, vals_order_btf_array, M.nnz);

    MALLOC_INT_1DARRAY(order_blk_amd_array, M.ncol);
    init_value(order_blk_amd_array, M.ncol, (Int)0);
    MALLOC_INT_1DARRAY(btf_blk_nnz, nblks+1);
    init_value(btf_blk_nnz, nblks+1, (Int) 0);
    MALLOC_INT_1DARRAY(btf_blk_work, nblks+1);
    init_value(btf_blk_work, nblks+1, (Int) 0);

    //Find AMD blk ordering, get nnz, and get work
    btf_blk_amd( M, order_blk_amd_array, btf_blk_nnz, btf_blk_work);

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

    MALLOC_INT_1DARRAY(vals_order_blk_amd_array, M.nnz);
    //permute_col(M, order_blk_amd_array); //NDE: Track movement of vals (lin_ind of row,col) here
    permute_col_store_valperms(M, order_blk_amd_array, vals_order_blk_amd_array); //NDE: Track movement of vals (lin_ind of row,col) here
    permute_row(M, order_blk_amd_array);

    permute_inv(vals_perm_composition, vals_order_blk_amd_array, M.nnz);

    // retry with original vals ordering
    sort_matrix_store_valperms(M, vals_perm_composition);

    //changed col to row, error.
    //print to see issue
    //printMTX("A_TOTAL.mtx", M);

    //NDE at this point, vals_perm_composition stores the permutation of the vals array; will be needed during break_into_parts
    break_into_parts2(M, nblks, btf_tabs);

    //find schedule
    find_btf_schedule(M, nblks, btf_tabs);


    if (Options.verbose == BASKER_TRUE)
    {
      printf("Basker BTF Cut: %ld \n", (long)btf_tabs(btf_tabs_offset));
    }

#ifdef BASKER_DEBUG_ORDER_BTF
    printf("------------BTF CUT: %d --------------\n",  btf_tabs(btf_tabs_offset));
#endif

    return 0;
  }//end find_btf2


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry,Exe_Space>::break_into_parts
  (
   BASKER_MATRIX &M,
   Int           nblks,
   INT_1DARRAY   btf_tabs
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

    //Short circuit, 
    //If nblks  == 1, than only BTF_A exists
    if(nblks == 1)
    {

    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("Short Circuit part_call \n");
    #endif
      BTF_A = A;
      //Options.btf = BASKER_FALSE;
      btf_tabs_offset = 1;
      return 0;
    }

    //Step 1.
    Int t_size            = 0;
    Int scol              = M.ncol;
    Int blk_idx           = nblks;
    BASKER_BOOL  move_fwd = BASKER_TRUE;
    while(move_fwd==BASKER_TRUE)
    {

      Int blk_size = btf_tabs(blk_idx)-
        btf_tabs(blk_idx-1);

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

        t_size = t_size+blk_size;
        blk_idx = blk_idx-1;
        scol   = btf_tabs[blk_idx];
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
          scol = btf_tabs[blk_idx];
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
   INT_1DARRAY   btf_tabs
  )
  {
  #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Basker: break_into_parts2 called \n");
    printf("nblks: %d \n", nblks);
  #endif

    Options.btf = BASKER_TRUE;

    //Alg.  
    // A -> [BTF_A  BTF_B] 
    //      [0      BTF_C]
    //1. Run backward through the btf_tabs to find size C
    //2. Form A,B,C based on size in 1.

    //Short circuit, 
    //If nblks  == 1, than only BTF_A exists
    // NDE: In this case, vals_block_map_perm_pair is not allocated nor used - A is assigned to BTF_A directly
    if(nblks == 1)
    {
    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("Basker: break_into_parts2 - short circuit for single block case\n");
    #endif
      BTF_A = A;
      //Options.btf = BASKER_FALSE; // NDE: how is this handled???
      btf_tabs_offset = 1;
      return 0;
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

    //Step 1.
    //Find total work estimate
    Int total_work_estimate = 0;
    for(Int b = 0; b < nblks; b++) //nblks is input; determined during btf ordering - total BTF_A blocks AND BTF_C blocks
    {
      total_work_estimate += btf_blk_work(b); //determined prior, during btf ordering
    }
    //Set a class variable to use later
    btf_total_work = total_work_estimate;
    //printf("num_threads: %d epsilon: %f \n",
    //	   num_threads, 
    //	   ((double)1/num_threads) +
    //	   ((double)BASKER_BTF_IMBALANCE));
    Int break_size    = ceil((double)total_work_estimate*(
          ((double)1/num_threads) + 
          ((double)BASKER_BTF_IMBALANCE)));

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Basker: Break size: %d \n", break_size);
    #endif

    Int t_size            = 0;      //total size of cols from 'small' blocks in BTF_C: matrix ncols - t_size = BTF_A ncols
    Int scol              = M.ncol; // starting column of the BTF_C blocks
    Int blk_idx           = nblks;  // start at lower right corner block, move left and up the diagonal
                                    // note: blk_idx index acts as block id + 1; btf_tabs(nblks) is likely the end column number
    BASKER_BOOL  move_fwd = BASKER_TRUE;

    while(move_fwd==BASKER_TRUE)
    {
      Int blk_work = btf_blk_work(blk_idx-1);
      Int blk_size  = btf_tabs(blk_idx) - 
        btf_tabs(blk_idx-1);       // subtract the bounding column ids to determine size of the (square) block

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
      if( (blk_work < break_size) && (blk_idx > 1) )
      {
      #ifdef BASKER_DEBUG_ORDER_BTF
        printf("Basker: continue with fine structure btf blocks\n");
      #endif

        t_size = t_size+blk_size;
        blk_idx = blk_idx-1;
        scol   = btf_tabs[blk_idx];
      }
      //break due to size i.e. entered non-trivial large BTF_A block
      else if( blk_work >= break_size )
      {
      #ifdef BASKER_DEBUG_ORDER_BTF
        printf("Basker: break due to size\n");
      #endif
        move_fwd = BASKER_FALSE;
      }
      //break due to end i.e. no 'large' BTF_A block for ND; only fine BTF structure
      else if(blk_idx == 1)
      {
        //printf("break last blk\n");
        blk_idx = 0;
        t_size = t_size + blk_size;
        scol = btf_tabs[blk_idx];	
        move_fwd = BASKER_FALSE;
      }
      //should not be called
      else
      {
        BASKER_ASSERT(1==0, "btf order break");
        move_fwd = BASKER_FALSE;
      }
    }//end while(move_fwd)

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Basker: Done finding BTF2 cut.  Cut size: %d scol: %d \n", t_size, scol);
    printf("Basker: Done finding BTF2 cut. blk_idx: %d \n", blk_idx);
    //BASKER_ASSERT(t_size > 0, "BTF CUT SIZE NOT BIG ENOUGH\n");
    BASKER_ASSERT((scol >= 0) && (scol <= M.ncol), "SCOL\n");
    #endif

    //Comeback and change
    // btf_tabs_offset is offset to id of first block after diag blocks in BTF_A 
    btf_tabs_offset = blk_idx;

    //Step 2. Move into Blocks 
    MALLOC_INT_1DARRAY_PAIRS(vals_block_map_perm_pair, M.nnz); //this will store and map A.val indices to val indices of BTF_A, BTF_B, and BTF_C 
    // Begin with BTF_A
    if(btf_tabs_offset != 0)
    {
      //--Move A into BTF_A;
      BTF_A.set_shape(0, scol, 0, scol);
      BTF_A.nnz = M.col_ptr(scol);

    #ifdef BASKER_DEBUG_ORDER_BTF
      printf("Basker Init BTF_A. ncol: %d nnz: %d \n", scol, BTF_A.nnz);
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
        printf("Basker copy column: %d into A_BTF, [%d %d] \n", 
            k, M.col_ptr(k), M.col_ptr(k+1));
      #endif

        for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
        {
          BTF_A.row_idx(annz) = M.row_idx(i);
          BTF_A.val(annz)     = M.val(i);  //NDE: Track movement of vals (lin_ind of row,col) here

          vals_block_map_perm_pair(i) = std::pair<Int,Int>(0,annz);

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
    BTF_B.set_shape(0, scol, scol, M.ncol-scol);
    BTF_C.set_shape(scol, M.ncol-scol, scol, M.ncol-scol);

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Baske Set Shape BTF_B: %d %d %d %d \n",
        BTF_B.srow, BTF_B.nrow,
        BTF_B.scol, BTF_B.ncol);
    printf("Basker Set Shape BTF_C: %d %d %d %d \n",
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
      }//over all nnz in k
    }//over all k

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("BTF_B nnz: %d \n", bnnz);
    printf("BTF_C nnz: %d \n", cnnz);
    #endif

    BTF_B.nnz = bnnz;
    BTF_C.nnz = cnnz;


    //Malloc need space
    if( (BTF_B.v_fill == BASKER_FALSE) && (BTF_B.nnz > 0) )
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
          printf("Adding nnz to Upper, %d %d \n", scol, M.row_idx[i]);
        #endif

          BASKER_ASSERT(BTF_B.nnz > 0, "BTF B uninit");
          //BTF_B.row_idx[bnnz] = M.row_idx[i];
          //Note: do not offset because B srow = 0
          BTF_B.row_idx(bnnz) = M.row_idx(i);
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
          //BTF_C.row_idx[cnnz] = M.row_idx[i];
          BTF_C.row_idx(cnnz) = M.row_idx(i)-scol;
          BTF_C.val(cnnz)     = M.val(i);

          vals_block_map_perm_pair(i) = std::pair<Int,Int>(2,cnnz);

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
#endif
