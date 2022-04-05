#ifndef __TACHO_TEST_DENSE_LINEAR_ALGEBRA_HPP__
#define __TACHO_TEST_DENSE_LINEAR_ALGEBRA_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"

#include "Tacho_Lapack_External.hpp"
#include "Tacho_Lapack_Team.hpp"

#include "Tacho_Chol.hpp"
#include "Tacho_Chol_External.hpp"
#include "Tacho_Chol_Internal.hpp"

#include "Tacho_LU.hpp"
#include "Tacho_LU_External.hpp"
#include "Tacho_LU_Internal.hpp"

namespace Test {

  TEST( DenseLinearAlgebra, chol ) {
    const int m = 20;
    const char uplo = 'U';

    value_type_matrix_type_host A("A", m, m);
    fill_spd_tridiag_matrix(A);

    int r_val(0);
#if defined(KOKKOS_ENABLE_SERIAL)
    value_type_matrix_type_host A_team("A_team", m, m);
    { /// cholesky A team

      Kokkos::deep_copy(A_team, A);
      
      const auto member = host_serial_member();
      Tacho::LapackTeam<value_type>::potrf(member,
                                    uplo, m,
                                    (value_type*)A_team.data(), (int)A_team.stride(1),
                                    &r_val);
      EXPECT_TRUE(r_val == 0);
    }
#endif
    value_type_matrix_type_host A_ext("A_ext", m, m);
    { /// cholesky A external

      Kokkos::deep_copy(A_ext, A);
      
      Lapack<value_type>::potrf(uplo, m,
                                (value_type*)A_ext.data(), (int)A_ext.stride(1),
                                &r_val);
      EXPECT_TRUE(r_val == 0);
    }
    {
      auto check = [](value_type_matrix_type_host A, value_type_matrix_type_host B) {
        const magnitude_type eps = atsv::epsilon()*100; 
        for (int i=0;i<m;++i) 
          for (int j=0;j<m;++j) {
            EXPECT_NEAR(atsv::real(A(i,j)), atsv::real(B(i,j)), eps);
            EXPECT_NEAR(atsv::imag(A(i,j)), atsv::imag(B(i,j)), eps);
          }
      };
      /// TODO:: check A = U'U
#if defined(KOKKOS_ENABLE_SERIAL)
      check(A_ext, A_team);
#endif
    }
  }

}
#endif
