#ifndef BASKER_ORDER_BTF_HPP
#define BASKER_ORDER_BTF_HPP

#include "basker_types.hpp"
#include "basker_sswrapper.hpp"
#include <assert.h>

//Depends on SuiteSparse in Amesos
#ifdef HAVE_AMESOS
#include "amesos_btf_decl.h"
#endif

//#define BASKER_DEBUG_ORDER_BTF

namespace BaskerNS
{

  

  /*
  template <class Int>
  BASKER_INLINE
  int my_strong_component
   (
    //BASKER_MATRIX &M,
    //Int           &nblks,
    //INT_1DARRAY   &perm, 
    //INT_1DARRAY   &perm_in,
    //INT_1DARRAY   &CC
   Int           &n,
   Int           *col_ptr,
   Int           *row_idx,
   Int           &nblks,
   Int           *perm,
   Int           *perm_in,
   Int           *CC
   )
  {
    return -1;
  }//end strong_component


  template <>
  BASKER_INLINE
  int my_strong_component <int>
  (
   //BASKER_MATRIX &M,
   int           &n,
   int           *col_ptr,
   int           *row_idx,
   int           &nblks,
   int           *perm,
   int           *perm_in,
   int           *CC
   //INT_1DARRAY   &perm, 
   //INT_1DARRAY   &perm_in,
   //INT_1DARRAY   &CC
   )
  {
    //typedef long int  l_Int;
    typedef int l_Int;

    //l_Int p[M.nrow]; //output row_per
    l_Int p[n]; //output row_per
    //l_Int r[M.nrow+1]; //comp_tabs
    l_Int r[n+1]; //comp_tabs
    //We will want to add option to use q in the future
 
    l_Int work[n*4];
    //l_Int work[M.nrow*4];

    //printf("before amesos call \n");
   
    //nblks = amesos_btf_strongcomp(M.ncol,&(M.col_ptr[0]),
    //&(M.row_idx[0]), 
//				    &(perm_in[0]), p, r, work);
  //
  
    nblks = amesos_btf_strongcomp(n, col_ptr,
				    row_idx, 
				    perm_in, p, r, work);
    //printf("after amesos call \n");



    #ifdef BASKER_DEBUG_ORDER_BTF
    
    printf("\nBTF perm: \n");
  //for(Int i=0; i <M.nrow; i++)
    for(Int i=0; i < n; i++)
      {
	printf("%d, ", p[i]);
      }
    
    printf("\n\nBTF tabs: <right> \n");
    for(Int i=0; i < nblks+1; i++)
      {
	printf("%d, ", r[i]);
      }
    printf("\n");
    #endif

  //BASKER_ASSERT(M.nrow > 0, "M.nrow btf");
  BASKER_ASSERT(n > 0, "M.nrow btf");
  //MALLOC_INT_1DARRAY(perm,M.nrow);
  // MALLOC_INT_1DARRAY(perm, n);
  for(l_Int i = 0; i < n; i++)
    {
      perm[p[i]] = i;
    }
    BASKER_ASSERT((nblks+1) > 0, "nblks+1 btf");
    //MALLOC_INT_1DARRAY(CC, nblks+1);
    for(l_Int i = 0; i < nblks+1; i++)
      {
	CC[i] = r[i];
      }

    return 0;
  }

  template <>
  int my_strong_component <long>
  (
   //BASKER_MATRIX &M,
   //long           &nblks,
   //INT_1DARRAY   &perm, 
   //INT_1DARRAY   &perm_in,
   //INT_1DARRAY   &CC
   
   long           &n,
   long           *col_ptr,
   long          *row_idx,
   long           &nblks,
   long           *perm,
   long           *perm_in,
   long          *CC
   )
  {
    typedef long  l_Int;

    l_Int p[n]; //output row_per
    l_Int r[n+1]; //comp_tabs
    //We will want to add option to use q in the future
 
    l_Int work[n*4];

    //printf("before amesos call \n");
    
    //nblks = amesos_btf_l_strongcomp(M.ncol,&(M.col_ptr[0]),
//				    &(M.row_idx[0]), 
//				    &(perm_in[0]), p, r, work);
    
     nblks = amesos_btf_l_strongcomp(n,
				    col_ptr,
				     row_idx, 
				    perm_in, p, r, work);
    //printf("after amesos call \n");



    #ifdef BASKER_DEBUG_ORDER_BTF
    
    printf("\nBTF perm: \n");
    for(Int i=0; i <n; i++)
      {
	printf("%d, ", p[i]);
      }
    
    printf("\n\nBTF tabs: <right> \n");
    for(l_Int i=0; i < nblks+1; i++)
      {
	printf("%d, ", r[i]);
      }
    printf("\n");
    #endif

    BASKER_ASSERT(n > 0, "M.nrow btf");
    //MALLOC_INT_1DARRAY(perm,M.nrow);
    for(l_Int i = 0; i < n; i++)
      {
	perm[p[i]] = i;
      }
    BASKER_ASSERT((nblks+1) > 0, "nblks+1 btf");
    //MALLOC_INT_1DARRAY(CC, nblks+1);
    for(l_Int i = 0; i < nblks+1; i++)
      {
	CC[i] = r[i];
      }

    return 0;
  }//strong_component<long int, Entry, Exe_Space>

  */

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::find_btf(BASKER_MATRIX &M)
  {
    Int          nblks = 0;
    //INT_1DARRAY  btf_perm = order_btf_array;
    //INT_1DARRAY  btf_tabs;

    //printf("before strong comp \n");
    strong_component(M,nblks,order_btf_array,btf_tabs);
    //strong_component(M,nblks,btf_perm,btf_tabs);
    btf_flag = true;
    //printf("after strong comp \n");

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


    //printf("before perm_col\n");
    permute_col(M, order_btf_array);
    //permute_col(M, btf_perm);
    //printf("after perm_col\n");
    //printf("before perm_row \n");
    permute_row(M, order_btf_array);
    //permute_row(M, btf_perm);
    //printf("after perm_row\n");


    //printMTX("A_BTF_OR.mtx", M);

    //printf("before break \n");
    break_into_parts(M, nblks, btf_tabs);
    //printf("after break \n");


    btf_nblks = nblks;


    //printMTX("C_BTF.mtx", BTF_C, BASKER_TRUE);
    //printMTX("B_BTF_BEFORE_PERM.mtx", BTF_B);

    //#ifdef BASKER_DEBUG_ORDER_BTF
    printf("------------BTF CUT: %d --------------\n", 
	   btf_tabs(btf_tabs_offset));
    //#endif

    return 0;
  }//end find BTF

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
	
        //#ifdef BASKER_DEBUG_ORDER_BTF
	printf("Short Circuit part_call \n");
	//#endif
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

	Int blk_size = btf_tabs(blk_idx)-btf_tabs(blk_idx-1);


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

    //#ifdef BASKER_DEBUG_ORDER_BTF
    printf("Done finding BTF cut.  Cut size: %d scol: %d \n",
	   t_size, scol);
    //BASKER_ASSERT(t_size > 0, "BTF CUT SIZE NOT BIG ENOUGH\n");
    BASKER_ASSERT((scol >= 0) && (scol < M.ncol), "SCOL\n");
    //#endif
    
    //Comeback and change
    btf_tabs_offset = blk_idx;


    //2. Move into Blocks

    if(btf_tabs_offset != 0)
      {
    //--Move A into BTF_A;
    BTF_A.set_shape(0, scol, 0, scol);
    BTF_A.nnz = M.col_ptr(scol);
    
    //#ifdef BASKER_DEBUG_ORDER_BTF
    printf("Init BTF_A. ncol: %d nnz: %d \n",
	   scol, BTF_A.nnz);
    //#endif

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



  //JDB: break_into_part2 needs cleaned up to fit new interface

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
    printf("break_into_parts called \n");
    #endif
    
    //Alg.  
    //      [BTF_D  BTF_E]
    // A -> [BTF_A  BTF_B] 
    //      [0      BTF_C]
    //1. Run backward through the btf_tabs to find size C
    //2. Form A,B,C based on size in 1.


    //Step 1.
    Int t_size            = 0;
    Int scol              = M.ncol;
    Int blk_idx           = nblks;
    BASKER_BOOL  move_fwd = BASKER_TRUE;

    //Bottom right to upper left
    while(move_fwd==BASKER_TRUE)
      {

	Int blk_size = btf_tabs(blk_idx)-btf_tabs(blk_idx-1);


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
	    printf("Cut: blk_size: %d percent: %f \n",
		   blk_size, ((double)t_size+blk_size)/(double)M.ncol);
	    move_fwd = BASKER_FALSE;
	  }
      }//end while(move_fwd)

    //Upper left to bottom right
    Int top_t_size        = 0;
    Int top_scol          = M.scol;
    Int top_blk_idx       = 0;
    move_fwd = BASKER_TRUE;

    //Bottom right to upper left
    while(move_fwd==BASKER_TRUE)
      {

	Int blk_size = btf_tabs(top_blk_idx+1)-
	  btf_tabs(top_blk_idx);


	#ifdef BASKER_DEBUG_ORDER_BTF
	printf("move_back loop \n");
	BASKER_ASSERT(blk_idx>=0, "btf blk idx off");
	BASKER_ASSERT(blk_size>0, "btf blk size wrong");
	printf("blk_idx: %d blk_size: %d \n", 
	       blk_idx, blk_size);
	std::cout << blk_size << std::endl;
	#endif

	if((blk_size < Options.btf_large) &&
	   ((((double)top_t_size+blk_size)/(double)M.ncol) < Options.btf_max_percent))
	  {
	    #ifdef BASKER_DEBUG_ORDER_BTF
	    printf("first choice \n");
	    printf("blksize test: %d %d %d \n",
		   blk_size, Options.btf_large, 
		   BASKER_BTF_LARGE);
	    printf("blkpercent test: %f %f %f \n", 
		   ((double)top_t_size+blk_size)/(double)M.ncol, 
		   Options.btf_max_percent, 
		   (double) BASKER_BTF_MAX_PERCENT);
	    #endif

	    top_t_size  = top_t_size+blk_size;
	    top_blk_idx = top_blk_idx-1;
	    top_scol    = btf_tabs[top_blk_idx];
	  }
	else
	  {
	    //printf("second choice \n");
	    printf("Cut: blk_size: %d percent: %f \n",
		   blk_size, ((double)top_t_size+blk_size)/(double)M.ncol);
	    move_fwd = BASKER_FALSE;
	  }
      }//end while(move_fwd)



    //#ifdef BASKER_DEBUG_ORDER_BTF
    printf("Done finding Lower BTF cut.  Cut size: %d scol: %d \n",
	   t_size, scol);
    printf("Done finding Upper BTF cut.  Cut size: %d scol: %d \n", top_t_size, top_scol);
    BASKER_ASSERT(t_size > 0, "BTF CUT SIZE NOT BIG ENOUGH\n");
    BASKER_ASSERT((scol > 0) && (scol < M.ncol), "SCOL\n");
    //#endif
    
    //Comeback and change
    btf_tabs_offset = blk_idx;
    btf_top_tabs_offset = top_blk_idx;


    //2. Move into Blocks

    //--Move A into BTF_A;
    BTF_A.set_shape(top_scol+1, scol, top_scol+1, scol);
    //BTF_A.set_shape(0, scol, 0, scol);
    BTF_A.nnz = M.col_ptr(scol) - M.col_ptr(top_scol);
    //BTF_A.nnz = M.col_ptr(scol);
    
    #ifdef BASKER_DEBUG_ORDER_BTF
    printf("Init BTF_A. ncol: %d nnz: %d \n",
	   scol, BTF_A.nnz);
    #endif
    MALLOC_INT_1DARRAY(BTF_A.col_ptr, BTF_A.ncol+1);
    MALLOC_INT_1DARRAY(BTF_A.row_idx, BTF_A.nnz);
    MALLOC_ENTRY_1DARRAY(BTF_A.val, BTF_A.nnz);
    BTF_A.fill();
    
    Int annz = 0;
    for(Int k = top_scol; k < scol; ++k)
    //for(Int k = 0; k < scol; ++k)
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

    
    ///=========================================//
    //in move
 
    ///===========================================///


    //Fill in B and C at the same time
    INT_1DARRAY cws;
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
		       scol, M.row_idx[i]);
		#endif
		bnnz++;
	      }
	    else
	      {
		#ifdef BASKER_DEBUG_ORDER_BTF
		printf("Adding nnz to Lower, %d %d \n",
		       scol, M.row_idx[i]);
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
    MALLOC_INT_1DARRAY(BTF_B.col_ptr, BTF_B.ncol+1);
    MALLOC_INT_1DARRAY(BTF_B.row_idx, BTF_B.nnz);
    MALLOC_ENTRY_1DARRAY(BTF_B.val, BTF_B.nnz);
    BTF_B.fill();
    MALLOC_INT_1DARRAY(BTF_C.col_ptr, BTF_C.ncol+1);
    MALLOC_INT_1DARRAY(BTF_C.row_idx, BTF_C.nnz);
    MALLOC_ENTRY_1DARRAY(BTF_C.val, BTF_C.nnz);
    BTF_C.fill();

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
		       M.val[i]);
		#endif
		//BTF_C.row_idx[cnnz] = M.row_idx[i];
		BTF_C.row_idx(cnnz) = M.row_idx(i)-scol;
		BTF_C.val(cnnz)     = M.val(i);
		cnnz++;
	      }
	  }//over all nnz in k
	BTF_B.col_ptr(k-scol+1) = bnnz;
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


  /*
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::strong_component
  (
   BASKER_MATRIX &M,
   Int           &nblks,
   INT_1DARRAY   perm_in,
   INT_1DARRAY   perm,
   INT_1DARRAY   CC
   )
  {
    printf("General Form of strong_component called \n");
    return -1;
  }//strong_component<Int,Entry,Exe_Space>
  */

  /*
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::strong_component
  (
   BASKER_MATRIX &M,
   Int           &nblks,
   INT_1DARRAY   perm, 
   INT_1DARRAY   CC
   )
  {
    printf("General Form of strong_component called \n");
    return -1;
  }//string_component<Int,Entry,Exe_Space>
  */

#ifdef HAVE_AMESOS
  
  /*
  template <class Int,class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::strong_component
  (
   BASKER_MATRIX &M,
   int           &nblks,
   INT_1DARRAY   perm,
`   INT_1DARRAY   CC
   )
  {
    typedef int    l_int;
    
    INT_1DARRAY perm_in;
    MALLOC_INT_1DARRAY(perm_in, A.ncol);
    for(l_int i = 0; i < A.ncol; i++)
      {
	perm_in[i] = i;
      }
    strong_compnent(M,nblks,perm_in, perm, CC);
   
    return 0;

  }//end strong_component<int>
  */

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
    //printf("SC one \n");
    //my_strong_component(M,nblks,perm,perm_in, CC);
    BaskerSSWrapper<Int>::my_strong_component(M.ncol,
			&(M.col_ptr(0)),
			&(M.row_idx(0)),
			nblks,
			&(perm(0)),
			&(perm_in(0)), 
			&(CC(0)));

    FILE *fp;
    fp = fopen("btf.txt", "w");
    for(Int i = 0; i < M.ncol; i++)
      {
        fprintf(fp, "%d \n", perm(i));
      }
    fclose(fp);

    
    printf("FOUND NBLKS: %d \n", nblks);

    return 0;

  }//end strong_component <long int>


  /*
  template <class Entry, class Exe_Space>
  BASKER_INLINE
  Int Basker<int, Entry, Exe_Space>::strong_component
  (
   BASKER_MATRIX &M,
   int           &nblks,
   INT_1DARRAY   perm,
   INT_1DARRAY   perm_in,
   INT_1DARRAY   CC
   )
  {  
    int p[M.nrow]; //output row_per
    int r[M.nrow+1]; //comp_tabs
    //We will want to add option to use q in the future
 
    int work[M.nrow*4];
    
    nblks = amesos_btf_strongcomp(M.ncol,&(M.col_ptr[0]),
			  &(M.row_idx[0]), 
			  q, p, r, work);


    MALLOC_INT_1DARRAY(perm,M.nrow);
    for(int i = 0; i < M.nrow; i++)
      {
	perm[i] = p[i];
      }
    MALLOC_INT_1DARRAY(CC, nblks+1);
    for(int i = 0; i < nblks+1; i++)
      {
	CC[i] = r[i];
      }
    return -1;
  }//strong_component<int, Entry, Exe_Space>
  */

  /*
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::strong_component
  (
   BASKER_MATRIX &M,
   Int           &nblks,
   INT_1DARRAY   &perm, 
   INT_1DARRAY   &perm_in,
   INT_1DARRAY   &CC
   )
  {
    typedef long int  l_Int;

    l_Int p[M.nrow]; //output row_per
    l_Int r[M.nrow+1]; //comp_tabs
    //We will want to add option to use q in the future
 
    l_Int work[M.nrow*4];

    //printf("before amesos call \n");
    nblks = amesos_btf_l_strongcomp(M.ncol,&(M.col_ptr[0]),
				    &(M.row_idx[0]), 
				    &(perm_in[0]), p, r, work);
    //printf("after amesos call \n");



    #ifdef BASKER_DEBUG_ORDER_BTF
    
    printf("\nBTF perm: \n");
    for(Int i=0; i <M.nrow; i++)
      {
	printf("%d, ", p[i]);
      }
    
    printf("\n\nBTF tabs: <right> \n");
    for(Int i=0; i < nblks+1; i++)
      {
	printf("%d, ", r[i]);
      }
    printf("\n");
    #endif

    BASKER_ASSERT(M.nrow > 0, "M.nrow btf");
    MALLOC_INT_1DARRAY(perm,M.nrow);
    for(Int i = 0; i < M.nrow; i++)
      {
	perm(p[i]) = i;
      }
    BASKER_ASSERT((nblks+1) > 0, "nblks+1 btf");
    MALLOC_INT_1DARRAY(CC, nblks+1);
    for(Int i = 0; i < nblks+1; i++)
      {
	CC(i) = r[i];
      }

    return 0;
  }//strong_component<long int, Entry, Exe_Space>
  */
#endif // End HAVE_AMESOS

}//end namespace BaskerNS

#endif //end BASKER_ORDER_BTF_HPP
