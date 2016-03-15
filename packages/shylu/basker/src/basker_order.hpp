#ifndef BASKER_ORDER_HPP
#define BASKER_ORDER_HPP

#include <impl/Kokkos_Timer.hpp>

//Basker Includes
#include "basker_types.hpp"
#include "basker_util.hpp"

#include "basker_order_match.hpp"
#include "basker_order_scotch.hpp"
#include "basker_order_btf.hpp"
#include "basker_order_amd.hpp"

//#undef BASKER_DEBUG_ORDER

namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::default_order()
  {
    match_ordering(0);
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
    init_tree(perm, nblks, parts, row_tabs, 
	      col_tabs, tree_tabs,0);

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
    match_ordering(0);
    //printf("DEBUG1: done match\n");
    //for debuging
    sort_matrix(A);
    //printMTX("A_match.mtx", A);
   
    //2. BTF ordering on whole matrix
    // Gets estimate of work on all blocks
    //currently finds btf-hybrid and permutes
    //A -> [BTF_A, BTF_C; 0 , BTF B]


    printf("outter num_threads:%d \n", num_threads);
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
	//INT_1DARRAY csymamd_perm = order_csym_array;
	MALLOC_INT_1DARRAY(order_csym_array, BTF_A.ncol+1);
	//MALLOC_INT_1DARRAY(csymamd_perm, BTF_A.ncol+1);
	init_value(order_csym_array, BTF_A.ncol+1,(Int) 0);
	//init_value(csymamd_perm, BTF_A.ncol+1,(Int) 0);
	
	csymamd_order(BTF_A, order_csym_array, cmember);
	//csymamd_order(BTF_A, csymamd_perm, cmember);
	
	//permute(BTF_A, csymamd_perm, csymamd_perm);
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
  
    
    //printf("Done with ordering\n");
    
    return 0;
  }//end btf_order

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::btf_order2()
  {
    
    //Kokkos::Impl::Timer timer_one;
    
    //printf("-----------BTF2-------------\n");

    //1. Matching ordering on whole matrix
    //currently finds matching and permutes
    //found bottle-neck to work best with circuit problems
    sort_matrix(A); //May want to remove ? (need)
    match_ordering(0);
    sort_matrix(A); //May want to remove? (need)

   
    //2. BTF ordering on whole matrix
    //currently finds btf-hybrid and permutes
    //A -> [BTF_A, BTF_C; 0 , BTF B]
    //printf("outter num_threads:%d \n", num_threads);
    MALLOC_INT_1DARRAY(btf_schedule, num_threads+1);
    init_value(btf_schedule, num_threads+1, (Int)0);
  

    find_btf2(A); 

    //TEst 
    //sort_matrix(BTF_A);
    //printMTX("A_TEST.mtx", BTF_A);
    //sort_matrix(BTF_C);

    //btf_tabs_offset = 0;
    //btf_tabs
    
    //std::cout << "Order Timer 1: "
    //	      << timer_one.seconds()
    //	      << std::endl;

    if((btf_tabs_offset != 0) )
      {
	//printf("BTF_OFFSET CALLED\n");

        //  printf("A/B block stuff called\n");
	//3. ND on BTF_A
	//currently finds ND and permute BTF_A
	//Would like to change so finds permuation, 
	//and move into 2D-Structure
	//printMTX("A_BTF_FROM_A.mtx", BTF_A);
	sort_matrix(BTF_A);

	//printMTX("BTF_A_TEST.mtx", BTF_A);

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
	//INT_1DARRAY csymamd_perm = order_csym_array;
	MALLOC_INT_1DARRAY(order_csym_array, BTF_A.ncol+1);
	//MALLOC_INT_1DARRAY(csymamd_perm, BTF_A.ncol+1);
	init_value(order_csym_array, BTF_A.ncol+1,(Int) 0);
	//init_value(csymamd_perm, BTF_A.ncol+1,(Int) 0);
	//printf("============CALLING CSYMAMD============\n");
	csymamd_order(BTF_A, order_csym_array, cmember);
	//csymamd_order(BTF_A, csymamd_perm, cmember);
	
	//permute(BTF_A, csymamd_perm, csymamd_perm);
	permute_col(BTF_A, order_csym_array);
	sort_matrix(BTF_A);
	permute_row(BTF_A, order_csym_array);
	sort_matrix(BTF_A);
	

	//printMTX("A_BTF_AMD1.mtx", BTF_A);
	
	
	if(btf_nblks > 1)
	  {
	    permute_row(BTF_B, order_csym_array);
	    sort_matrix(BTF_B);
	    //printMTX("B_BTF_AMD.mtx", BTF_B);
	    sort_matrix(BTF_C);
	    //printMTX("C_BTF_AMD.mtx", BTF_C);
	  }
    

	//printMTX("BTF_A.mtx", BTF_A);
    
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
	
	/*
	printf("===PRINT==\n");
	printMTX("BTF_A.mtx", BTF_A); 
	*/

      }//if btf_tab_offset == 0

    
    if(btf_nblks > 1)
      {
	sort_matrix(BTF_C);
	//printMTX("C_TEST.mtx", BTF_C);
	//Permute C


	/*
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
	    printMTX("BTF_B.mtx", BTF_B);
	  }

	printMTX("BTF_C.mtx", BTF_C);
	*/
      }
	
  

    //printf("Done with ordering\n");
    
    return 0;
  }//end btf_order2

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::order_incomplete()
  {
     
    Kokkos::Impl::Timer timer_one;
    
    if(Options.matching == BASKER_TRUE)
      {
	match_ordering(0);
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
	scotch_partition(A);
      }
    else
      {
	scotch_partition(BTF_A);
      }
    return 0;
  }//end partition()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry,Exe_Space>::match_ordering(int option)
  {


    //printf("match_order called\n");

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

    //Int job = 2; //5 is the default for SuperLU_DIST
    //INT_1DARRAY mperm = order_match_array;
    MALLOC_INT_1DARRAY(order_match_array, A.nrow);
    //MALLOC_INT_1DARRAY(mperm, A.nrow);
    if(Options.incomplete == BASKER_FALSE)
      {
	mwm(A,order_match_array);
      }
    else
      {
	for(Int i = 0; i < A.nrow; i++)
	  {
	    order_match_array(i) = i;
	  }
      }
    //mc64(job,order_match_array);
    //mc64(job,mperm);

    match_flag = BASKER_TRUE;
    
    #ifdef BASKER_DEBUG_ORDER
    printf("Matching Perm \n");
    for(Int i = 0; i < A.nrow; i++)
      {
	printf("%d, \n", order_match_array(i));
	//printf("%d, \n", mperm[i]);
      }
    printf("\n");
    #endif

    //We want to test what the match ordering does if
    //have explicit zeros
    #ifdef BASKER_DEBUG_ORDER
    FILE *fp;
    fp = fopen("match_order.txt", "w");
    for(Int i = 0; i < A.nrow; i++)
      {
	fprintf(fp, "%d \n", order_match_array(i));
      }
    fclose(fp);
    #endif

    
    permute_row(A,order_match_array);
    //permute_row(A,mperm);
    //May have to call row_idx sort
    return 0;
    
  }//end match_ordering()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::scotch_partition
  (BASKER_MATRIX &M)
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
    //permute_col(M, part_tree.permtab);
    ///permute_row(M, part_tree.permtab);
    permute_row(M, part_tree.permtab);
    permute_col(M, part_tree.permtab);

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
	//temp(i) = vec(p(i));
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
	//temp(p(i)) = vec(i);
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
   Int n, //n = size(vec) //n > m 
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
	//temp(p(i)) = vec(i);
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


   /*GOBACK and make this good*/
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::permute_col
  (
   BASKER_MATRIX &M,
   INT_1DARRAY col
   )
  {
    if((M.ncol == 0)||(M.nnz == 0))
      return 0;

    Int n = M.ncol;
    Int nnz = M.nnz;
    //printf("Using n: %d nnz: %d \n", n, nnz);
    INT_1DARRAY temp_p;
    MALLOC_INT_1DARRAY(temp_p, n+1);
    init_value(temp_p, n+1, (Int)0);
    INT_1DARRAY temp_i;
    MALLOC_INT_1DARRAY(temp_i, nnz);
    init_value(temp_i, nnz, (Int)0);
    ENTRY_1DARRAY temp_v;
    MALLOC_ENTRY_1DARRAY(temp_v, nnz);
    init_value(temp_v, nnz, (Entry)0.0);
    //printf("done with init \n");
   
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
            temp_i (ko) = M.row_idx (k);
            temp_v (ko) = M.val (k);
            ko++;
          }
      }
    
    //copy back int A
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
  int Basker<Int,Entry, Exe_Space>::permute_row(
				      BASKER_MATRIX &M,
				       INT_1DARRAY row)
  {

    if(M.nnz == 0)
      {
	return 0;
      }
    Int nnz = M.nnz;
    INT_1DARRAY temp_i;
    MALLOC_INT_1DARRAY(temp_i, nnz);
    init_value(temp_i, nnz, (Int)0);

    //permute
    for(Int k = 0; k < nnz; k++)
      {
        temp_i[k] = row[M.row_idx[k]];
      }
    //Copy back
    for(Int k = 0; k < nnz; k++)
      {
        M.row_idx[k] = temp_i[k];
      }
    FREE_INT_1DARRAY(temp_i);
    return 0;
  }//end permute_row(matrix,int)


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::sort_matrix(BASKER_MATRIX &M)
  {
    
    if(M.nnz == 0)
      return 0;
    //printf("COMEBACK and add sort_matrix\n");
    //just use insertion sort
    for(Int k = 0; k < M.ncol; k++)
      {
	
	Int start_row = M.col_ptr[k];
	for(Int i = M.col_ptr[k]+1; i < M.col_ptr[k+1]; i++)
	  {
	    Int jj = i;
	    while((jj > start_row) && 
		  (M.row_idx[jj-1] > M.row_idx[jj]))
	      {

		/*
		printf("k: %d start: %d j: %d \n",
			k, start_row, jj);
		printf("index: %d %d val: %f %f \n",
		       M.row_idx[jj], M.row_idx[jj-1], 
		       M.val[jj], M.val[jj-1]);
		*/

		//swap
		Int   t_row_idx = M.row_idx[jj-1];
		Entry t_row_val = M.val[jj-1];
		
		M.row_idx[jj-1] = M.row_idx[jj];
		M.val[jj-1]     = M.val[jj];

		M.row_idx[jj] = t_row_idx;
		M.val[jj]     = t_row_val;
		
		jj = jj-1;

	      }

	  }

      }//end over all columns



    #ifdef BASKER_DEBUG_ORDER
    //printf("\n\n Matrix Sorted \n \n");
    //M.print();
    //printf("\n\n");
    #endif

    return 0;
  }//end sort_matrix()

}//end namespace basker
#endif //end ifndef basker_order_hpp
