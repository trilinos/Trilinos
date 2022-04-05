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
    const int m = 10;
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

    auto check_this = [=](const value_type_matrix_type_host AA) {
      value_type_matrix_type_host L("L", m, m);
      {
        const bool transL(true);
        copy_lower_triangular(L, transL, value_type(0), AA);
      }
      value_type_matrix_type_host U("U", m, m);
      {
        const bool transU(false);
        copy_upper_triangular(U, transU, value_type(0), AA);
      }
      value_type_matrix_type_host A_lu("A_lu", m, m);      
      compute_A(A_lu, L, U);
      check_same_matrix(A, A_lu);      
    };

    {
      check_this(A_ext);
#if defined(KOKKOS_ENABLE_SERIAL)
      check_this(A_team);
#endif
    }
  }


  TEST( DenseLinearAlgebra, lu ) {
    const int m = 5, n = 8;

    value_type_matrix_type_host A("A", m, n);
    //fill_spd_tridiag_matrix(A);
    fill_random_matrix(A);

    int r_val(0);
#if defined(KOKKOS_ENABLE_SERIAL)
    value_type_matrix_type_host A_team("A_team", m, n);
    ordinal_type_array_type_host P_team("P_team", 4*m);
    { /// LU A team
      Kokkos::deep_copy(A_team, A);
      
      const auto member = host_serial_member();
      Tacho::LapackTeam<value_type>::getrf(member,
                                           m, n,
                                           (value_type*)A_team.data(), (int)A_team.stride(1),
                                           (int*)P_team.data(),
                                           &r_val);
      EXPECT_TRUE(r_val == 0);
    }
#endif
    value_type_matrix_type_host A_ext("A_ext", m, n);
    ordinal_type_array_type_host P_ext("P_ext", 4*m);
    { /// LU A external
      Kokkos::deep_copy(A_ext, A);
      
      Lapack<value_type>::getrf(m, n,
                                (value_type*)A_ext.data(), (int)A_ext.stride(1),
                                (int*)P_ext.data(),
                                &r_val);
      EXPECT_TRUE(r_val == 0);
    }

    auto check_this = [=](const ordinal_type_array_type_host P, 
                          const value_type_matrix_type_host AA) {
      value_type_matrix_type_host L("L", m, m);
      {
        const bool transL(false);
        auto range = Kokkos::pair<int,int>(0,m);
        auto AL = Kokkos::subview(AA, range, range);
        copy_lower_triangular(L, transL, value_type(1), AL);
      }
      value_type_matrix_type_host U("U", m, n);
      {
        const bool transU(false);
        copy_upper_triangular(U, transU, value_type(0), AA);
      }

      value_type_matrix_type_host A_lu("A_lu", m, n);
      compute_A(A_lu, L, U);
      
      value_type_matrix_type_host A_perm("A_perm", m, n);
      
      copy_matrix(A_perm, A);
      apply_lapack_pivots_left_no_trans(A_perm, P);

      check_same_matrix(A_perm, A_lu);   
    };

    {
      check_this(P_ext, A_ext);
#if defined(KOKKOS_ENABLE_SERIAL)
      check_this(P_team, A_team);
#endif
    }
  }

}
#endif
