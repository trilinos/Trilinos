#ifndef __TACHOEXP_TRSM_LEFT_UPPER_CONJTRANS_EXTERNAL_HPP__
#define __TACHOEXP_TRSM_LEFT_UPPER_CONJTRANS_EXTERNAL_HPP__


/// \file  Tacho_Trsm_Left_Upper_ConjTrans_External.hpp
/// \brief BLAS triangular solve matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#ifdef HAVE_SHYLUTACHO_TEUCHOS
#include "Teuchos_BLAS.hpp"
#endif

namespace Tacho {

  namespace Experimental {

    /// BLAS Trsm
    /// =========
    template<>
    template<typename PolicyType,
             typename MemberType,
             typename DiagType,
             typename ScalarType,
             typename ViewTypeA,
             typename ViewTypeB>
    KOKKOS_INLINE_FUNCTION
    int
    Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
    ::invoke(const PolicyType &policy,
             const MemberType &member,
             const DiagType diagA,
             const ScalarType alpha,
             const ViewTypeA &A,
             const ViewTypeB &B) {
      typedef typename ViewTypeA::non_const_value_type value_type;
      typedef typename ViewTypeB::non_const_value_type value_type_b;
      typedef typename ViewTypeA::array_layout array_layout_a;
      typedef typename ViewTypeB::array_layout array_layout_b;

      static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
      static_assert(ViewTypeB::rank == 2,"B is not rank 2 view.");
      
      static_assert(std::is_same<array_layout_a,Kokkos::LayoutLeft>::value,
                    "A does not have Kokkos::LayoutLeft.");
      static_assert(std::is_same<array_layout_b,Kokkos::LayoutLeft>::value,
                    "B does not have Kokkos::LayoutLeft.");

      static_assert(std::is_same<value_type,value_type_b>::value,
                    "A and B do not have the same value type.");
      
      const ordinal_type m = A.dimension_0(), n = B.dimension_1();

      if (m > 0 && n > 0) {
        if (get_team_rank(member) == 0) {
#if                                                     \
  defined( HAVE_SHYLUTACHO_TEUCHOS ) &&                 \
  defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
          Teuchos::BLAS<ordinal_type,value_type> blas;
          blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::CONJ_TRANS,
                    (diagA.tag == Diag::Unit::tag ? Teuchos::UNIT_DIAG : Teuchos::NON_UNIT_DIAG),
                    m, n,
                    alpha,
                    A.data(), A.stride_1(),
                    B.data(), B.stride_1());
#else
          TACHO_TEST_FOR_ABORT( true, MSG_NOT_HAVE_PACKAGE("Teuchos BLAS") );
          
#endif
        }
      }
      return 0;
    }
  }
}
#endif
