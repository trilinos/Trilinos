#ifndef BASKER_ORDER_BTF_HPP
#define BASKER_ORDER_BTF_HPP

#include "basker_types.hpp"
#include "basker_sswrapper.hpp"

//Depends on SuiteSparse in Amesos
#ifdef HAVE_AMESOS
#include "amesos_btf_decl.h"
#endif

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

    Int break_size    = ceil((double)total_work_estimate*(
			   ((double)1/num_threads) + 
	      		     ((double)BASKER_BTF_IMBALANCE)));


    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Total schedul size: %d \n", total_work_estimate);
    printf("Break schedule size: %d \n", break_size);
    printf("Total num blks: %d \n", btf_nblks);
    #endif
    
    Int t_size            = 0;
    Int t_loc             = 0;
    btf_schedule(0) = btf_tabs_offset;
    //BASKER_BOOL  move_fwd = BASKER_TRUE; //NU

    for(Int b = btf_tabs_offset; b < btf_nblks;
	b++)
      {
	Int blk_work = btf_blk_work(b);
	t_size += blk_work;
	
	#ifdef BASKER_DEBUG_ORDER_BTF
	printf("t: %d blk: %d work: %d twork: %d \n", 
	       t_loc,b, blk_work, t_size);
	#endif

	if(((t_size > break_size) && (t_loc < num_threads-1)) ||
	   (b == btf_nblks-1))
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
  int Basker<Int,Entry, Exe_Space>::find_btf(BASKER_MATRIX &M)
  {
    Int          nblks = 0;

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
    /*
    printf("\nBTF perm: \n");
    for(Int i=0; i <M.nrow; i++)
      {
	printf("%d, ", order_btf_array(i));
	//printf("%d, ", btf_perm(i));
      }
    */
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
    Int          nblks = 0;


    strong_component(M,nblks,order_btf_array,btf_tabs);

    btf_nblks = nblks;

    btf_flag = BASKER_TRUE;

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("BTF nblks returned: %d \n", nblks);
    //BASKER_ASSERT(nblks>1, "NOT ENOUGH BTF BLOCKS");
    #endif

    #ifdef BASKER_DEBUG_ORDER_BTF
    if(nblks<2)
      {
	printf("BTF did not find enough blks\n");
      }
    #endif


    #ifdef BASKER_DEBUG_ORDER_BTF
    /*
    printf("\nBTF perm: \n");
    for(Int i=0; i <M.nrow; i++)
      {
	printf("%d, ", order_btf_array(i));
	//printf("%d, ", btf_perm(i));
      }
    */
    printf("num_threads: %d \n", num_threads);
    printf("\n\nBTF tabs: \n");
    for(Int i=0; i < nblks+1; i++)
      {
	printf("%d, ", btf_tabs(i));
      }
    printf("\n");
    #endif

    permute_col(M, order_btf_array);
    permute_row(M, order_btf_array);

    MALLOC_INT_1DARRAY(order_blk_amd_array, M.ncol);
    init_value(order_blk_amd_array, M.ncol, (Int)0);
    MALLOC_INT_1DARRAY(btf_blk_nnz, nblks+1);
    init_value(btf_blk_nnz, nblks+1, (Int) 0);
    MALLOC_INT_1DARRAY(btf_blk_work, nblks+1);
    init_value(btf_blk_work, nblks+1, (Int) 0);


    //Find AMD blk ordering, get nnz, and get work
    btf_blk_amd( M, order_blk_amd_array,
		 btf_blk_nnz, btf_blk_work);


    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("blk_perm:\n");
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
    

    permute_col(M, order_blk_amd_array);
    permute_row(M, order_blk_amd_array);
    sort_matrix(M);

    //changed col to row, error.
    //print to see issue
    //printMTX("A_TOTAL.mtx", M);
    
       
    break_into_parts2(M, nblks, btf_tabs);

    //find schedule
    find_btf_schedule(M, nblks, btf_tabs);


    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("------------BTF CUT: %d --------------\n", 
	   btf_tabs(btf_tabs_offset));
    #endif

    return 0;
  }//end find BTF(nnz)

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
	    printf("Cut: blk_size: %d percent: %f \n",
		   blk_size, ((double)t_size+blk_size)/(double)M.ncol);
	    
	    if((((double)t_size+blk_size)/(double)M.ncol)
	       == 1.0)
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


    //2. Move into Blocks

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
    BTF_B.set_shape(0   , scol,
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
      //if(BTF_B.v_fill == BASKER_FALSE)
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

    //printf("\n\n");
    //printf("DEBUG\n");
    //BTF_C.print();
    //printf("\n\n");

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
    printf("break_into_parts2 called \n");
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

    //Short circuit for incomplete
    if(Options.incomplete == BASKER_TRUE)
      {
	//printf("Short Ciruit incomplete\n");
	BTF_A = A;
	btf_nblks = 1;
	btf_tabs_offset = 1;
	return 0;
      }

    //Step 1.
    //Find total work estimate
    Int total_work_estimate = 0;
    for(Int b = 0; b < nblks; b++)
      {
	total_work_estimate += btf_blk_work(b);
      }
    //Set a class variable to use later
    btf_total_work = total_work_estimate;
    //printf("Total work estimate: %d \n",
    //	   total_work_estimate);
    //printf("num_threads: %d epsilon: %f \n",
    //	   num_threads, 
    //	   ((double)1/num_threads) +
    //	   ((double)BASKER_BTF_IMBALANCE));
    Int break_size    = ceil((double)total_work_estimate*(
		((double)1/num_threads) + 
		((double)BASKER_BTF_IMBALANCE)));

    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Break size: %d \n", break_size);
    #endif
    
    Int t_size            = 0;
    Int scol              = M.ncol;
    Int blk_idx           = nblks;
    BASKER_BOOL  move_fwd = BASKER_TRUE;

    while(move_fwd==BASKER_TRUE)
      {

	//printf("------TEST blk_idx: %d \n",
	// blk_idx);
	Int blk_work = btf_blk_work(blk_idx-1);
	Int blk_size  = btf_tabs(blk_idx) - 
	                        btf_tabs(blk_idx-1);


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
	if(((blk_work < break_size) &&
	    (blk_idx > 1)))
	  {
	    #ifdef BASKER_DEBUG_ORDER_BTF
	    printf("first choice \n");
	    #endif

	    t_size = t_size+blk_size;
	    blk_idx = blk_idx-1;
	    scol   = btf_tabs[blk_idx];
	   
	  }
	//break due to size
	else if(blk_work >= break_size)
	  {
	    #ifdef BASKER_DEBUG_ORDER_BTF
	    printf("break due to size\n");
	    #endif
	    move_fwd = BASKER_FALSE;
	  }
	//break due to end
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
    printf("Done finding BTF2 cut.  Cut size: %d scol: %d \n",
	   t_size, scol);
    printf("Done finding BTF2 cut. blk_idx: %d \n", 
	   blk_idx);
    //BASKER_ASSERT(t_size > 0, "BTF CUT SIZE NOT BIG ENOUGH\n");
    
    BASKER_ASSERT((scol >= 0) && (scol <= M.ncol), "SCOL\n");
    #endif
    
    //Comeback and change
    btf_tabs_offset = blk_idx;


    //2. Move into Blocks

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
    BTF_B.set_shape(0   , scol,
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

    //Added check
    if((BTF_C.nrow == 0)||(BTF_C.ncol == 0))
      {
	return 0;
      }
    
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
      //if(BTF_B.v_fill == BASKER_FALSE)
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

    //printf("\n\n");
    //printf("DEBUG\n");
    //BTF_C.print();
    //printf("\n\n");
    
    return 0;

  }//end break_into_parts2 (based on imbalance)

  
 

#ifdef HAVE_AMESOS
  
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


    //printf("===this strong comp called====");

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

    //printf("SC one \n");
    //my_strong_component(M,nblks,perm,perm_in, CC);
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
        fprintf(fp, "%d \n", perm(i));
      }
    fclose(fp);
    #endif

    
    //printf("FOUND NBLKS: %d \n", nblks);

    return 0;

  }//end strong_component <long int>


#endif // End HAVE_AMESOS

}//end namespace BaskerNS

#endif //end BASKER_ORDER_BTF_HPP
