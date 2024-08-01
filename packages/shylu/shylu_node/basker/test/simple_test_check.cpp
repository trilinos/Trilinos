// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Kokkos_Core.hpp"
#include "test_util.hpp"
#include "shylubasker_decl.hpp"
#include "shylubasker_def.hpp"

int main(int argc, char* argv[])
{

  using Int = int;
  using Entry = double;
  using Exe_Space = Kokkos::OpenMP;
  using Entry_ViewType = Kokkos::View<Entry*, Kokkos::HostSpace>;
  using Int_ViewType = Kokkos::View<Int*, Kokkos::HostSpace>;

  Int m = 3;
  Int n = 3;
  Int nnz = 5;

  const Int nthreads = 1;
  Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(nthreads));
  {

    std::cout << "Create Matrix Views" << std::endl;
    std::cout << "Create col_ptrv" << std::endl;
    Int_ViewType  col_ptrv("col_ptr", m+1);
    std::cout << "Create row_idxv" << std::endl;
    Int_ViewType  row_idxv("row_idx", nnz);
    std::cout << "Create valv" << std::endl;
    Entry_ViewType valv("val", nnz);

    std::cout << "  Pointer of Matrix Views" << std::endl;
    auto col_ptr = col_ptrv.data();
    auto row_idx = row_idxv.data();
    auto val = valv.data();


    // Create 3x3 matrix:
    //
    //   1     0     0
    //   1     2     0
    //   1     0     3

    std::cout << "Init Matrix Views" << std::endl;
    col_ptrv(0) = 0;
    col_ptrv(1) = 3;
    col_ptrv(2) = 4;
    col_ptrv(3) = 5;

    row_idxv(0) = 0;
    row_idxv(1) = 1;
    row_idxv(2) = 2;
    row_idxv(3) = 1;
    row_idxv(4) = 2;

    valv(0) = 1;
    valv(1) = 1;
    valv(2) = 1;
    valv(3) = 2;
    valv(4) = 3;

    std::cout << "Create RHS Views" << std::endl;
    // non-transpose b
    Entry_ViewType nyv("ny", m);
    nyv(0) = 1;
    nyv(1) = 3;
    nyv(2) = 4;

    auto ny = nyv.data();

    // transpose b
    Entry_ViewType tyv("ty", m);
    tyv(0) = 3;
    tyv(1) = 2;
    tyv(2) = 3;

    auto ty = tyv.data();

    // Solution will be "ones" to
    // A  \ ny
    // A' \ ty

    std::cout << "Create known soln View" << std::endl;
    Entry_ViewType soln("soln", m);
    soln(0) = 1;
    soln(1) = 1;
    soln(2) = 1;

    // x - initialized to 0's
    std::cout << "Create computed soln View" << std::endl;
    Entry_ViewType xv("x", m);
    auto x = xv.data();


    //Start Basker
    // Non-transpose test
    {
      std::cout << "Non-Transpose test" << std::endl;
      BaskerNS::Basker<Int, Entry, Exe_Space> mybasker;
      //---Options
      mybasker.Options.transpose          = BASKER_FALSE;
      mybasker.Options.symmetric          = BASKER_FALSE;
      mybasker.Options.matching           = BASKER_TRUE;
      mybasker.Options.matching_type      = BASKER_MATCHING_BN;

      mybasker.SetThreads(nthreads);
      std::cout << "Setting Threads:" << nthreads << std::endl;
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
    
      mybasker.Solve(ny,x);
      std::cout << "Done with Solve, Time: "
        	      << totalTime(ttime, myTime()) << std::endl;

      //std::cout << "Non-transpose solution" << std::endl;
      //for (Int i = 0; i < m; ++i) {
      //  std::cout << "  xv(" << i << ") = " << xv(i) << std::endl;
      //}
      bool is_correct = true;
      for (Int i = 0; i < m; ++i) {
        if (xv(i) != 1)
          is_correct = false;
      }
      std::cout << "  Non-transpose: is_correct = " << is_correct << std::endl;

      mybasker.Finalize();
    }


    //Start Basker
    // Transpose test
    {
      deep_copy(xv,0);
      std::cout << "Transpose test" << std::endl;
      BaskerNS::Basker<Int, Entry, Exe_Space> mybasker;
      //---Options
      mybasker.Options.transpose          = BASKER_TRUE;
      mybasker.Options.symmetric          = BASKER_FALSE;
      mybasker.Options.matching           = BASKER_TRUE;
      mybasker.Options.matching_type      = BASKER_MATCHING_BN;

      mybasker.SetThreads(nthreads);
      std::cout << "Setting Threads:" << nthreads << std::endl;
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
    
      mybasker.Solve(ty,x);
      std::cout << "Done with Solve, Time: "
        	      << totalTime(ttime, myTime()) << std::endl;

      //std::cout << "Transpose solution" << std::endl;
      //for (Int i = 0; i < m; ++i) {
      //  std::cout << "  xv(" << i << ") = " << xv(i) << std::endl;
      //}

      bool is_correct = true;
      for (Int i = 0; i < m; ++i) {
        if (xv(i) != 1)
          is_correct = false;
      }
      std::cout << "  Transpose: is_correct = " << is_correct << std::endl;

      mybasker.Finalize();
    }

  }
  Kokkos::finalize();


} //end main
