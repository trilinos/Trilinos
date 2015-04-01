#pragma once
#ifndef __ICHOL_RIGHT_BLOCKED_HPP__
#define __ICHOL_RIGHT_BLOCKED_HPP__

/// \file ichol_right_blocked.hpp
/// \brief Blocked incomplete Chloesky factorization.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///
/// Unlike the dense matrix algebra, this blocked version of sparse
/// factorization does not lead to an efficient computation.
/// This algorithm is only for the testing and debugging purpose.

namespace Example {

  using namespace std;

  // use Upper Triangular part only
  template<>
  template<typename ParallelForType,
           typename CrsExecViewType>
  KOKKOS_INLINE_FUNCTION
  int
  IChol<Uplo::Upper,AlgoIChol::RightBlocked>
  ::invoke(const typename CrsExecViewType::policy_type::member_type &member,
           const CrsExecViewType &A) {

    typedef typename CrsExecViewType::value_type   value_type;
    typedef typename CrsExecViewType::ordinal_type ordinal_type;

    const ordinal_type mb = blocksize;

    CrsExecViewType ATL, ATR,      A00, A01, A02,
      /**/          ABL, ABR,      A10, A11, A12,
      /**/                         A20, A21, A22;

    Part_2x2(A,  ATL, ATR,
             /**/ABL, ABR,
             0, 0, Partition::TopLeft);

    while (ATL.NumRows() < A.NumRows()) {
      Part_2x2_to_3x3(ATL, ATR, /**/  A00, A01, A02,
                      /*******/ /**/  A10, A11, A12,
                      ABL, ABR, /**/  A20, A21, A22,
                      mb, mb, Partition::BottomRight);
      // -----------------------------------------------------
      A11.fillRowViewArray();
      A12.fillRowViewArray();      
      A22.fillRowViewArray();      

      int r_val = IChol<Uplo::Upper,AlgoIChol::RightUnblockedOpt1>
        ::invoke<ParallelForType,CrsExecViewType>(member, A11);

      if (r_val)
        return A00.NumRows() + r_val;

      Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,AlgoTrsm::ForRightBlocked>
        ::invoke<ParallelForType,value_type,CrsExecViewType>(member, Diag::NonUnit, 1.0, A11, A12);

      Herk<Uplo::Upper,Trans::ConjTranspose,AlgoHerk::ForRightBlocked>
        ::invoke<ParallelForType,value_type,CrsExecViewType>(member, -1.0, A12, 1.0, A22);

      // -----------------------------------------------------
      Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                       A10, A11, A12, /**/ /******/
                       A20, A21, A22, /**/ ABL, ABR,
                       Partition::TopLeft);
    }

    return 0;
  }

}

#endif
