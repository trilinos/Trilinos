#ifndef __TACHOEXP_CHOL_EXTERNAL_HPP__
#define __TACHOEXP_CHOL_EXTERNAL_HPP__

/// \file  TachoExp_Chol_Upper_ExternalLapack.hpp
/// \brief LAPACK upper Cholesky factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Teuchos_LAPACK.hpp"

namespace Tacho {

  namespace Experimental {

    /// LAPACK Chol
    /// ===========
    template<typename ArgUplo>
    struct Chol<ArgUplo,Algo::External> {
      template<typename SchedType,
               typename MemberType,
               typename ViewTypeA>
      inline
      static int
      invoke(const SchedType &sched,
             const MemberType &member,
             const ViewTypeA &A) {
        typedef typename ViewTypeA::non_const_value_type value_type;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");

        int r_val = 0;      
        
        const ordinal_type m = A.dimension_0();
        if (m > 0) {
          if (get_team_rank(member) == 0) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
#if defined( HAVE_SHYLUTACHO_MKL )
            const char lapack_uplo = ArgUplo::param;
            if      (std::is_same<value_type,float>::value) 
              r_val = LAPACKE_spotrf (LAPACK_COL_MAJOR, lapack_uplo, 
                                      m, 
                                      (float *)A.data(), (lapack_int)A.stride_1());
            else if (std::is_same<value_type,double>::value) 
              r_val = LAPACKE_dpotrf (LAPACK_COL_MAJOR, lapack_uplo, 
                                      m, 
                                      (double *)A.data(), (lapack_int)A.stride_1());
            else if (std::is_same<value_type,Kokkos::complex<float> >::value ||
                     std::is_same<value_type,   std::complex<float> >::value)
              r_val = LAPACKE_cpotrf (LAPACK_COL_MAJOR, lapack_uplo, 
                                      m, 
                                      (lapack_complex_float *)A.data(), (lapack_int)A.stride_1());
            else if (std::is_same<value_type,Kokkos::complex<double> >::value ||
                     std::is_same<value_type,   std::complex<double> >::value)
              r_val = LAPACKE_zpotrf (LAPACK_COL_MAJOR, lapack_uplo, 
                                      m, 
                                      (lapack_complex_double *)A.data(), (lapack_int)A.stride_1());
            else {
              TACHO_TEST_FOR_ABORT( true, ">> Datatype is not supported.");                           
            }
#else
            typedef typename TypeTraits<value_type>::std_value_type std_value_type;
            Teuchos::LAPACK<ordinal_type,std_value_type> lapack;
            const char lapack_uplo = ArgUplo::param;
            lapack.POTRF(lapack_uplo,
                         m, 
                         (std_value_type*)A.data(), A.stride_1(),
                         &r_val);

            TACHO_TEST_FOR_EXCEPTION(r_val, std::runtime_error, 
                                     "LAPACK (potrf) returns non-zero error code.");
#endif
#else
            TACHO_TEST_FOR_ABORT( true, ">> This function is only allowed in host space." );
#endif
          }
        }
        return r_val;
      }
    };
  }
}

#endif
