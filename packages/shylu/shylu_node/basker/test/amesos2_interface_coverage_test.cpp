// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
  Joshua Dennis Booth
  Test the Amesos2 Interface calls with refactoring
  Test should have no inputs
*/

#include "test_util.hpp"
#include "shylubasker_decl.hpp"
#include "shylubasker_def.hpp"
#include <Kokkos_Core.hpp>

int main(int argc, char* argv[])
{

  typedef long           Int;
  typedef double         Entry;
  typedef Kokkos::OpenMP Exe_Space;

  std::string mname;

  if(argc > 2)
  {
    std::cout <<"Test Input is only the coverage matrix"
      << std::endl;
    return -1;
  }
  else if (argc == 2)
  {
    mname = std::string(argv[1]);
  }
  else
  {
    mname = std::string("coverage.mtx");
  }
  
  //Load inital information
  //Matrix
  Int m, n, nnz; 
  m = n = nnz = 0;
  Int *col_ptr = nullptr;
  Int *row_idx = nullptr;
  Entry *val = nullptr;
  
  std::cout << "Matrix read" << std::endl;
  double rmatrix = myTime();
  readMatrix<Int,Entry>(mname, m, n, nnz, 
			&col_ptr, &row_idx, &val);
  std::cout << "Read Matrix, Time: " 
	    << totalTime(rmatrix,myTime()) << std::endl;
  
  //RHS
  Int vn, vm;
  vn = vm = 0;
  Entry *x, *xhat, *y;
  x = xhat = y = NULL;
  vn = n;
  x = new Entry[vn]();
  xhat = new Entry[vn]();
  //Populate Vector
  {
    vm = m;
    y = new Entry[m]();
    for(Int i = 0; i < vm; i++)
    {
      xhat[i] = (Entry) i;
    }
    multiply<Int,Entry>(m,n,col_ptr,row_idx,val, xhat, y);
    for(Int i = 0; i < vm; i++)
    {
      //std::cout  << "y " << y[i] << std::endl;
      xhat[i] = (Entry) 0.0;
    }
  }
  
  //Starting up Kokkos
  int nthreads = 4; // We will not use all 4 in all tests
  Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(nthreads));
  std::cout << "Kokkos Settings" << std::endl;
  std::cout << "hwloc aval: " 
	    << Kokkos::hwloc::available()<< std::endl;
  std::cout << "numa count: " 
	    << Kokkos::hwloc::get_available_numa_count() 
	    << std::endl;
  std::cout << "thrd numa:  " 
	    << Kokkos::hwloc::get_available_cores_per_numa() 
	    << std::endl;
 
  //-----------------------Start Basker (Test - 1, 1 thread)-----------------
  {
    std::cout << "==============Starting Test 1, 1 Threads===========" 
              << std::endl;

//    int result = 0; // NDE: warning unused
    BaskerNS::Basker<Int, Entry, Exe_Space> mybasker;
    //---Options
    mybasker.Options.same_pattern       = BASKER_FALSE;
    mybasker.Options.verbose            = BASKER_FALSE;
    mybasker.Options.verbose_matrix_out = BASKER_FALSE;
    mybasker.Options.realloc            = BASKER_TRUE;
    mybasker.Options.transpose          = BASKER_FALSE;
    mybasker.Options.symmetric          = BASKER_FALSE;
    mybasker.Options.AtA                = BASKER_TRUE;
    mybasker.Options.A_plus_At          = BASKER_TRUE;
    mybasker.Options.matching           = BASKER_TRUE;
    mybasker.Options.matching_type      = BASKER_MATCHING_BN;
    mybasker.Options.btf                = BASKER_TRUE;
    mybasker.Options.btf_max_percent    = BASKER_BTF_MAX_PERCENT;
    mybasker.Options.btf_large          = BASKER_BTF_LARGE;
    mybasker.Options.no_pivot           = BASKER_FALSE;
    //mybasker.Options.pivot_tol          = .001;
    //mybasker.Options.pivot_bias         = .001;
    //mybasker.Options.btf_prune_size      = 2;
   
    mybasker.SetThreads(1);
    std::cout << "Setting Threads:" << 1 << std::endl;
    double stime = myTime();
    mybasker.Symbolic(m,n,nnz,col_ptr,row_idx,val);
    std::cout << "Done with Symbolic, Time: " 
	      << totalTime(stime, myTime()) << std::endl;
    double ftime = myTime();
    mybasker.Factor(m,n,nnz,col_ptr,row_idx,val);
    std::cout << "Done with Factor, Time: "
	      << totalTime(ftime, myTime()) << std::endl;
    //mybasker.DEBUG_PRINT();
    double ttime = myTime();
    Int *lperm;
    Int *rperm;
    mybasker.GetPerm(&lperm,&rperm);
    
    mybasker.Solve(y,x);
    std::cout << "Done with Solve, Time: "
	      << totalTime(ttime, myTime()) << std::endl;

    multiply<Int,Entry>(m,n,col_ptr,row_idx,val, x, xhat);
    for(Int i = 0; i < m; i++)
    {
      xhat[i] = y[i] - xhat[i];
    }
    std::cout << "||X||: " << norm2<Int,Entry>(n,x)
	      << " ||Y-AX||: " << norm2<Int,Entry>(m,xhat)
	      << std::endl;
    
    //Refactor
    double rftime = myTime();
    mybasker.Factor(m,n,nnz,col_ptr,row_idx,val);
    std::cout << "Done with Refactor Factor, Time: "
	      << totalTime(rftime, myTime()) << std::endl;
    //ReSolve
    double rttime = myTime();
    mybasker.Solve(y,x);
    std::cout << "Done with Refactor Solve, Time: "
	      << totalTime(rttime, myTime()) << std::endl;

    multiply<Int,Entry>(m,n,col_ptr,row_idx,val, x, xhat);
    for(Int i = 0; i < m; i++)
    {
      xhat[i] = y[i] - xhat[i];
    }
    
    std::cout << "||X||: " << norm2<Int,Entry>(n,x)
	      << " ||Y-AX||: " << norm2<Int,Entry>(m,xhat)
	      << std::endl;

    mybasker.Finalize();
  }

  //-----------------------Start Basker (Test - 2, 2 thread)--------------------//
  {
    std::cout << "============Starting Test 2, 2 Threads================" 
              << std::endl;

//    int result = 0; // NDE: warning unused
    BaskerNS::Basker<Int, Entry, Exe_Space> mybasker;
    //---Options
    mybasker.Options.same_pattern       = BASKER_FALSE;
    mybasker.Options.verbose            = BASKER_FALSE;
    mybasker.Options.verbose_matrix_out = BASKER_FALSE;
    mybasker.Options.realloc            = BASKER_TRUE;
    mybasker.Options.transpose          = BASKER_FALSE;
    mybasker.Options.symmetric          = BASKER_FALSE;
    mybasker.Options.AtA                = BASKER_TRUE;
    mybasker.Options.A_plus_At          = BASKER_TRUE;
    mybasker.Options.matching           = BASKER_TRUE;
    mybasker.Options.matching_type      = BASKER_MATCHING_BN;
    mybasker.Options.btf                = BASKER_TRUE;
    mybasker.Options.btf_max_percent    = BASKER_BTF_MAX_PERCENT;
    mybasker.Options.btf_large          = BASKER_BTF_LARGE;
    mybasker.Options.no_pivot           = BASKER_FALSE;
    //mybasker.Options.pivot_tol          = .001;
    //mybasker.Options.pivot_bias         = .001;
    //mybasker.Options.btf_prune_size      = 2;
    
    mybasker.SetThreads(2);
    std::cout << "Setting Threads:" << 2 << std::endl;
    double stime = myTime();
    mybasker.Symbolic(m,n,nnz,col_ptr,row_idx,val);
    std::cout << "Done with Symbolic, Time: " 
	      << totalTime(stime, myTime()) << std::endl;
    double ftime = myTime();
    mybasker.Factor(m,n,nnz,col_ptr,row_idx,val);
    std::cout << "Done with Factor, Time: "
	      << totalTime(ftime, myTime()) << std::endl;
    //mybasker.DEBUG_PRINT();
    double ttime = myTime();
    Int *lperm;
    Int *rperm;
    mybasker.GetPerm(&lperm,&rperm);
    
    mybasker.Solve(y,x);
    std::cout << "Done with Solve, Time: "
	      << totalTime(ttime, myTime()) << std::endl;

    multiply<Int,Entry>(m,n,col_ptr,row_idx,val, x, xhat);
    for(Int i = 0; i < m; i++)
    {
      xhat[i] = y[i] - xhat[i];
    }
    std::cout << "||X||: " << norm2<Int,Entry>(n,x)
	      << " ||Y-AX||: " << norm2<Int,Entry>(m,xhat)
	      << std::endl;
    
    //Refactor
    double rftime = myTime();
    mybasker.Factor(m,n,nnz,col_ptr,row_idx,val);
    std::cout << "Done with Refactor Factor, Time: "
	      << totalTime(rftime, myTime()) << std::endl;
    //ReSolve
    double rttime = myTime();
    mybasker.Solve(y,x);
    std::cout << "Done with Refactor Solve, Time: "
	      << totalTime(rttime, myTime()) << std::endl;

    multiply<Int,Entry>(m,n,col_ptr,row_idx,val, x, xhat);
    for(Int i = 0; i < m; i++)
    {
      xhat[i] = y[i] - xhat[i];
    }

    std::cout << "||X||: " << norm2<Int,Entry>(n,x)
	      << " ||Y-AX||: " << norm2<Int,Entry>(m,xhat)
	      << std::endl;

    mybasker.Finalize();
  }


  
  //-----------------------Start Basker (Test - 4, 4 thread)--------------------//
/*
  {
    std::cout << "============Starting Test 4, 4 Threads=============" 
              << std::endl;

//    int result = 0; // NDE: warning unused
    BaskerNS::Basker<Int, Entry, Exe_Space> mybasker;
    //---Options
    mybasker.Options.same_pattern       = BASKER_FALSE;
    mybasker.Options.verbose            = BASKER_FALSE;
    mybasker.Options.verbose_matrix_out = BASKER_FALSE;
    mybasker.Options.realloc            = BASKER_TRUE;
    mybasker.Options.transpose          = BASKER_FALSE;
    mybasker.Options.symmetric          = BASKER_FALSE;
    mybasker.Options.AtA                = BASKER_TRUE;
    mybasker.Options.A_plus_At          = BASKER_TRUE;
    mybasker.Options.matching           = BASKER_TRUE;
    mybasker.Options.matching_type      = BASKER_MATCHING_BN;
    mybasker.Options.btf                = BASKER_TRUE;
    mybasker.Options.btf_max_percent    = BASKER_BTF_MAX_PERCENT;
    mybasker.Options.btf_large          = BASKER_BTF_LARGE;
    mybasker.Options.no_pivot           = BASKER_FALSE;
    //mybasker.Options.pivot_tol          = .001;
    //mybasker.Options.pivot_bias         = .001;
    //mybasker.Options.btf_prune_size      = 2;
    
    mybasker.SetThreads(4);
    std::cout << "Setting Threads:" << 4 << std::endl;
    double stime = myTime();
    mybasker.Symbolic(m,n,nnz,col_ptr,row_idx,val);
    std::cout << "Done with Symbolic, Time: " 
	      << totalTime(stime, myTime()) << std::endl;
    double ftime = myTime();
    mybasker.Factor(m,n,nnz,col_ptr,row_idx,val);
    std::cout << "Done with Factor, Time: "
	      << totalTime(ftime, myTime()) << std::endl;
    //mybasker.DEBUG_PRINT();
    double ttime = myTime();
    Int *lperm;
    Int *rperm;
    mybasker.GetPerm(&lperm,&rperm);
    
    mybasker.Solve(y,x);
    std::cout << "Done with Solve, Time: "
	      << totalTime(ttime, myTime()) << std::endl;

    multiply<Int,Entry>(m,n,col_ptr,row_idx,val, x, xhat);
    for(Int i = 0; i < m; i++)
    {
     	xhat[i] = y[i] - xhat[i];
    }
    std::cout << "||X||: " << norm2<Int,Entry>(n,x)
	      << " ||Y-AX||: " << norm2<Int,Entry>(m,xhat)
	      << std::endl;
    
    //Refactor
    double rftime = myTime();
    mybasker.Factor(m,n,nnz,col_ptr,row_idx,val);
    std::cout << "Done with Refactor Factor, Time: "
	      << totalTime(rftime, myTime()) << std::endl;
    //ReSolve
    double rttime = myTime();
    mybasker.Solve(y,x);
    std::cout << "Done with Refactor Solve, Time: "
	      << totalTime(rttime, myTime()) << std::endl;

    multiply<Int,Entry>(m,n,col_ptr,row_idx,val, x, xhat);
    for(Int i = 0; i < m; i++)
    {
     	xhat[i] = y[i] - xhat[i];
    }
    
    std::cout << "||X||: " << norm2<Int,Entry>(n,x)
	      << " ||Y-AX||: " << norm2<Int,Entry>(m,xhat)
	      << std::endl;

    mybasker.Finalize();
  }
*/
  
  Kokkos::finalize();

}//end main
