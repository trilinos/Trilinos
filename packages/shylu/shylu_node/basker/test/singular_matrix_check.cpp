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

  Int m = 2;
  Int n = 2;
  Int nnz = 4;
  int error = 0;

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


    // Create 2x2 singular matrix:
    //
    //   2    -2
    //  -2     2
    //

    std::cout << "Init Matrix Views" << std::endl;
    col_ptrv(0) = 0;
    col_ptrv(1) = 2;
    col_ptrv(2) = 4;

    row_idxv(0) = 0;
    row_idxv(1) = 1;
    row_idxv(2) = 0;
    row_idxv(3) = 1;

    valv(0) = 2.;
    valv(1) = -2.;
    valv(2) = -2.;
    valv(3) = 2.;


    //Start Basker
    // Non-transpose test, supply as non-symmetric
    {
      std::cout << "Non-Transpose test" << std::endl;
      BaskerNS::Basker<Int, Entry, Exe_Space> mybasker;
      //---Options
      mybasker.Options.transpose          = BASKER_FALSE;
      mybasker.Options.symmetric          = BASKER_FALSE;
      mybasker.Options.matching           = BASKER_TRUE;
      mybasker.Options.matching_type      = BASKER_MATCHING_BN;
      //mybasker.Options.verbose            = BASKER_TRUE;

      mybasker.SetThreads(nthreads);
      std::cout << "Setting Threads:" << nthreads << std::endl;

      double stime = myTime();
      error = mybasker.Symbolic(m,n,nnz,col_ptr,row_idx,val);
      std::cout << "Done with Symbolic"
                << "\nError code: " << error
                << "\nTime: " 
        	      << totalTime(stime, myTime()) << std::endl;

      double ftime = myTime();
      try
      {
        error = mybasker.Factor(m,n,nnz,col_ptr,row_idx,val);
      }
      catch (std::runtime_error& e)
      {
        std::cout << " ** Factor threw exception **" << std::endl;
        error = 1;
      }
      std::cout << "Done with Factor"
                << "\nError code: " << error
                << "\nTime: " 
	              << totalTime(ftime, myTime()) << std::endl;
      //mybasker.DEBUG_PRINT();
    
      if (error == 0) {
        double ttime = myTime();
        Int *lperm;
        Int *rperm;
        mybasker.GetPerm(&lperm,&rperm);

        std::cout << "Create RHS Views" << std::endl;
        // Non-transpose b = A*xsoln
        // [-2; 2];
        Entry_ViewType nyv("ny", m);
        nyv(0) = -2.;
        nyv(1) = 2.;

        auto ny = nyv.data();

        // x - solution, initialized to 0's
        std::cout << "Create computed soln View" << std::endl;
        Entry_ViewType xv("x", m);
        auto x = xv.data();

        mybasker.Solve(ny,x);
        std::cout << "Done with Solve"
                  << "\nTime: "
                  << totalTime(ttime, myTime()) << std::endl;

        std::cout << "Non-transpose solution" << std::endl;
        for (Int i = 0; i < m; ++i) {
          std::cout << "  nyv(" << i << ") = " << nyv(i) << std::endl;
        }
        for (Int i = 0; i < n; ++i) {
          std::cout << "  xv(" << i << ") = " << xv(i) << std::endl;
        }
      }
      else {
        std::cout << "Error in factorization step - no solve attempted" << std::endl;
      }
      //std::cout << "  Non-transpose: is_correct = " << is_correct << std::endl;

      mybasker.Finalize();
    }

  }
  Kokkos::finalize();

  return (error == 0 ? 1 : 0);
} //end main
