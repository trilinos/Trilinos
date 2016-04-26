#ifndef BASKER_DEF_HPP
#define BASKER_DEF_HPP

//#define BASKER_TIME

/*Basker Includes*/
#include "shylubasker_decl.hpp"
#include "basker_matrix_decl.hpp"
#include "basker_matrix_def.hpp"
#include "basker_matrix_view_decl.hpp"
#include "basker_matrix_view_def.hpp"
#include "basker_tree.hpp"
#include "basker_sfactor.hpp"
#include "basker_sfactor_inc.hpp"
#include "basker_nfactor.hpp"
#include "basker_nfactor_inc.hpp"
#include "basker_solve_rhs.hpp"
#include "basker_util.hpp"
#include "basker_stats.hpp"
#include "basker_order.hpp"

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

    btf_total_work = 0;

  }//end Basker()
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Basker<Int ,Entry, Exe_Space>::~Basker()
  {
    
  }//end ~Basker()

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
        FREE_INT_1DARRAY(order_btf_array);
        btf_flag = BASKER_FALSE;
      }
    if(nd_flag == BASKER_TRUE)
      {
        FREE_INT_1DARRAY(order_scotch_array);
        nd_flag == BASKER_FALSE;
      }
    if(amd_flag == BASKER_TRUE)
      {
        FREE_INT_1DARRAY(order_csym_array);
        amd_flag == BASKER_FALSE;
      }
    
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
  int Basker<Int,Entry, Exe_Space>::InitMatrix(Int nrow, Int ncol, Int nnz, 
					       Int *col_ptr,
					       Int *row_idx, Entry *val)
  {
    //Note: jdb comeback to add trans option
    A.init_matrix("Original Matrix",
		  nrow, ncol, nnz, col_ptr, row_idx, val);
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
	//printf("btf_order called \n");
	btf_order();
	//printf("btf_order returned \n");
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
  int Basker<Int, Entry, Exe_Space>::Symbolic(Int option)
  {


    printf("calling symbolic \n");
    #ifdef BASKER_KOKKOS_TIME 
    Kokkos::Impl::Timer timer;
    #endif

    //symmetric_sfactor();
    sfactor();


    if(option == 0)
      {
        
      }
    else if(option == 1)
      {
	
        
      }

    #ifdef BASKER_KOKKOS_TIME
    stats.time_sfactor += timer.seconds();
    #endif

    return 0;
  }//end Symbolic


  //This is the interface for Amesos2
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Symbolic(Int nrow, Int ncol,
	         Int nnz, Int *col_ptr, Int *row_idx, Entry *val)
  {

    // printf("befor symbolic\n");
    if(Options.verbose == BASKER_TRUE)
      {
	printf("Basker Symbolic\n");
	printf("Matrix: %d %d %d \n",
	       nrow, ncol, nnz);
      }
    //Init Matrix A.
    if(matrix_flag == BASKER_TRUE)
      {
	printf("YOU CANNOT RERUN SYMBOLIC\n");
	return BASKER_ERROR;
      }
    else
      {

	//Kokkos::Impl::Timer timer_move;

	if(Options.transpose == BASKER_FALSE)
	  {
	    //printf("=======NO TRANS=====\n");
	    A.init_matrix("Original Matrix",
		  nrow, ncol, nnz, col_ptr, row_idx, val);
	    A.scol = 0;
	    A.srow = 0;
	    
	  }
	else
	  {
	    //printf("======TRANS=====\n");
	    //Will transpose and put in A using little extra
	    matrix_transpose(0, nrow,
			     0, ncol,
			     nnz,
			     col_ptr,
			     row_idx,
			     val,
			     A);
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

	//std::cout << "Transpose A: " << timer_move.seconds()
	//	  << std::endl;
      }

    //Init Ordering
    //Always will do btf_ordering
    //This should also call create tree
    if(order_flag == BASKER_TRUE)
      {
	printf("YOU CANNOT RERUN ORDER\n");
	return BASKER_ERROR;
      }
    else
      {
	//printf("btf_order called \n");
	//btf_order();
	
	Kokkos::Impl::Timer timer_order;
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

	//if(btf_tabs_offset != 0)
	if((Options.btf == BASKER_TRUE) &&
	   (btf_tabs_offset != 0))
	  {
	    basker_barrier.init(num_threads, 16, tree.nlvls );
	  }
	order_flag = BASKER_TRUE;

	if(Options.verbose == BASKER_TRUE)
	  {
	    printf("Basker P2P Thread Barriers Init\n");
	  }

	//std::cout << "Time Order/Init arrays " 
	//	  << timer_order.seconds()
	//	  << std::endl;

	//printf("btf_order done \n");
      }



    //printf("\n\n+++++++++++++++BREAKER BREAKER++++++++\n\n");

    

    if(symb_flag == BASKER_TRUE)
      {
        printf("YOU CANNOT RERUN SFACTOR\n");
	return BASKER_ERROR;
      }
    else
      {
	
       
	if(Options.incomplete == BASKER_FALSE)
	  {
	    sfactor();
	  }
	else
	  {
	    sfactor_inc();
	  }

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
    
    //printf("\nTEST ALM\n");
    //ALM(0)(0).info();
    //printf("\n");

    

    return 0;
	
  }//end Symbolic()
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Factor(Int option)
  {
    #ifdef BASKER_KOKKOS_TIME
    Kokkos::Impl::Timer timer;
    #endif
    
    factor_notoken(option);
    
    #ifdef BASKER_KOKKOS_TIME
    stats.time_nfactor += timer.seconds();
    #endif
    return 0;
  }//end Factor()

  //This is the interface for Amesos2
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::Factor(Int nrow, Int ncol,
	        Int nnz, Int *col_ptr, Int *row_idx, Entry *val) 
  {

    int err = 0;

    if (Options.verbose == BASKER_TRUE)
      {
	printf("Basker Factor Called\n");
	printf("Matrix: %d %d %d \n",
	       nrow, ncol, nnz);
      }

    
    /*
    int err = A.copy_values(nrow, ncol, nnz, col_ptr, 
			    row_idx, val);
    */

    if((Options.same_pattern == BASKER_TRUE) &&
       (Options.no_pivot == BASKER_FALSE))
      {
	printf("Warning: Same Pattern will not allow pivoting\n");
	Options.no_pivot = BASKER_TRUE;
      }

 			    
    if(Options.transpose == BASKER_FALSE)
      {
	//printf("=======NO TRANS=====\n");
	//A.init_matrix("Original Matrix",
	//	      nrow, ncol, nnz, col_ptr, row_idx, val);
	//A.scol = 0;
	//A.srow = 0;
	A.copy_values(nrow, ncol, nnz, col_ptr,
		      row_idx, val);
	//printf("Copy done\n");
	//printMTX("A_LOAD.mtx", A);
      }
    else
      {
	//printf("======TRANS=====\n");
	//Will transpose and put in A using little extra
	matrix_transpose(0, nrow,
			 0, ncol,
			 nnz,
			 col_ptr,
			 row_idx,
			 val,
			 A);
      }
    sort_matrix(A);
    if(Options.verbose_matrix_out == BASKER_TRUE)
      {
	printMTX("A_Factor.mtx", A);
      }
    matrix_flag = BASKER_TRUE;
	

    if(err == BASKER_ERROR)
      {
	return BASKER_ERROR;
      }
    //err = sfactor_copy();
    err = sfactor_copy2();
    if (Options.verbose == BASKER_TRUE)
      {
	printf("Basker Copy Structure Done \n");
      }

    //printf("Done with sfactor_copy: %d \n", err);
    if(err == BASKER_ERROR)
      {
	return BASKER_ERROR;
      }
    //printf("before notoken\n");

    //Kokkos::Impl::Timer timer;
   

    if(Options.incomplete == BASKER_FALSE)    
      {
	err = factor_notoken(0);
	//printf("Notoken called\n");
      }
    else
      {
	err = factor_inc_lvl(0);
      }
    if(err == BASKER_ERROR)
      {
	return BASKER_ERROR;
      }

    if(Options.verbose == BASKER_TRUE)
      {
	printf("Basker Factor Done \n");
      }

    /*
    std::cout << "Raw Factor Time: "
	      << timer.seconds()
	      << std::endl;
    */
    
    //DEBUG_PRINT();

    return 0;

  }//end Factor()

  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::Factor_Inc(Int Options)
  {
    factor_inc_lvl(Options);
    return 0;
  }


  //Interface for solve.... only doing paralllel solve righ now.
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::SolveTest()
  {
    test_solve();
    return 0;
  }//end SolveTest

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Solve(Entry *b, 
					   Entry *x)
  {
    // printf("Currrently not used \n");
    solve_interface(x,b);
    return 0;
  }//Solve(Entry *, Entry *);

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::Solve(Int nrhs,
                                         Entry *b,
                                         Entry *x)
  {
    solve_interface(nrhs,x,b);
    return 0;
  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Solve(ENTRY_1DARRAY b,
					   ENTRY_1DARRAY x)
  {
    printf("Currently not used \n");
    return -1;
  }//Solve(ENTRY_1D, ENTRY_1D);


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::Solve(Int nrhs, 
					   Entry *b, Entry *x, 
                                           Int option)
  {    
    int err = 0;
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
    if((nthreads != 1) && 
       (nthreads%2 != 0))
      {
	BASKER_ASSERT(0==1, 
		      "Number of thread error");
	//Set default 1
	num_threads = 1;
	return BASKER_ERROR;
      }

    //Next test if Kokkos has that many threads!
    //This is a common mistake in mpi-based apps
    #ifdef KOKKOS_HAVE_OPENMP
    int check_value = Kokkos::OpenMP::max_hardware_threads();
    if(nthreads > check_value)
      {
        BASKER_ASSERT(0==1,
                      "Number of thread not aval in Kokkos");
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
      return BASKER_ERROR;
    else
      return BASKER_SUCCESS;
  }//end GetLnnz();

  //Return nnz of U
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::GetUnnz(Int &Unnz)
  {
    (Unnz) = get_Unnz();
    if(Unnz == 0)
      return BASKER_ERROR;
    else
      return BASKER_SUCCESS;
  }//end GetUnnz()

  //Returns assembled L
  template<class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::GetL(Int &n, Int &nnz,
           Int **col_ptr, Int **row_idx,
           Entry **val)
  {
    get_L(n,nnz,col_ptr, row_idx, val);
    
    return BASKER_SUCCESS;
  }//end GetL()
  
  //returns assembles U
  template<class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::GetU(Int &n, Int &nnz,
           Int **col_ptr, Int **row_idx,
           Entry **val)
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


  //Debub tester function
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
	     order_match_array.dimension_0());
      }
    if(btf_flag == BASKER_TRUE)
      {
    printVec("btf.csc", order_btf_array,
	     order_btf_array.dimension_0());
    printVec("amdblk.csc", order_blk_amd_array,
	     order_blk_amd_array.dimension_0());
      }
    if(btf_tabs_offset != 0)
      {
	printVec("ND.csc", part_tree.permtab, 
		 part_tree.permtab.dimension_0());
      }
    if(amd_flag == BASKER_TRUE)
      {
	printVec("camd.csc", order_csym_array,
		 order_csym_array.dimension_0());
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
