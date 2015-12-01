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
#include "basker_nfactor.hpp"
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
    matrix_flag    = false;
    order_flag     = false;
    tree_flag      = false;
    symb_flag      = false;
    factor_flag    = false;
    workspace_flag = false;
    rhs_flag       = false;
    solve_flag     = false;
    nd_flag        = false;
    amd_flag       = false;


    //Default number of threads
    num_threads = 1;
    global_nnz = 0;

  }//end Basker()
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Basker<Int ,Entry, Exe_Space>::~Basker()
  {
    

  
  }//end ~Basker()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::InitMatrix(string filename)
  {    
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

    basker_barrier.init(num_threads, 16, tree.nlvls );

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

    //Init Matrix A.
    if(matrix_flag == BASKER_TRUE)
      {
	printf("YOU CANNOT RERUN SYMBOLIC\n");
	return BASKER_ERROR;
      }
    else
      {
	A.init_matrix("Original Matrix",
		  nrow, ncol, nnz, col_ptr, row_idx, val);
	A.scol = 0;
	A.srow = 0;
	sort_matrix(A);
	matrix_flag = BASKER_TRUE;
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
	btf_order();
	if(btf_tabs_offset != 0)
	  {
	    basker_barrier.init(num_threads, 16, tree.nlvls );
	  }
	order_flag = BASKER_TRUE;
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
	sfactor();
	symb_flag = BASKER_TRUE;
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
    
    //factor_tree(option);
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
    
    //printf("first in factor \n");
    //ALM(0)(0).info();

    int err = A.copy_values(nrow, ncol, nnz, col_ptr, row_idx, val);

    //ALM(0)(0).info();

    if(err == BASKER_ERROR)
      {
	return BASKER_ERROR;
      }
    err = sfactor_copy();
    printf("Done with sfactor_copy: %d \n", err);
    if(err == BASKER_ERROR)
      {
	return BASKER_ERROR;
      }
    printf("before notoken\n");
    err = factor_notoken(0);
    if(err == BASKER_ERROR)
      {
	return BASKER_ERROR;
      }
    return 0;

  }//end Factor()


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
    num_threads = nthreads;
    return 0;
  }//end SetThreads()

  //Return nnz of L
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::GetLnnz(Int *Lnnz)
  {
    (*Lnnz) = get_Lnnz();
    if(*Lnnz == 0)
      return -1;
    else
      return 0;
  }//end GetLnnz();

  //Return nnz of U
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::GetUnnz(Int *Unnz)
  {
    (*Unnz) = get_Unnz();
    if(*Unnz == 0)
      return -1;
    else
      return 0;
  }//end GetUnnz()

  //Returns assembled L
  template<class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::GetL(Int *n, Int *nnz,
           Int **col_ptr, Int **row_idx,
           Entry **val)
  {
    int err = get_L(n,nnz,col_ptr, row_idx, val);
    return err;
  }//end GetL()
  
  //returns assembles U
  template<class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::GetU(Int *n, Int *nnz,
           Int **col_ptr, Int **row_idx,
           Entry **val)
  {
    int err = get_U(n, nnz, col_ptr, row_idx, val);
    return err;
  }//end GetU()

  //returns global P
  template<class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::GetP(Int **p)
  {
    int err = get_p(p);
    return err;
  }//end GetP()

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
    #else
    //printL();
    #endif
    std::cout << "L printed " << std::endl;
    printU();
    std::cout << "U printed" << std::endl;
    //printRHS();
    std::cout << "RHS printed" << std::endl;
    //printSOL();
    std::cout << "SOL printed" << std::endl;
    //printTree();
    std::cout << "Tree printed" << std::endl;
    
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
