#ifndef __TACHO_TEST_HPP__
#define __TACHO_TEST_HPP__


#include "Tacho_CrsMatrixBase.hpp"

namespace Test {
  using namespace Tacho;

  using atsv = ArithTraits<value_type>;
  using atsm = ArithTraits<typename atsv::mag_type>;

  using crs_matrix_base_type_host = CrsMatrixBase<value_type,host_device_type>;
  using crs_matrix_base_type = CrsMatrixBase<value_type,device_type>;

  using ordinal_type_array_type_host = Kokkos::View<ordinal_type*,host_device_type>; 
  using size_type_array_type_host = Kokkos::View<size_type*,host_device_type>; 

  using value_type_array_type_host = Kokkos::View<value_type*,host_device_type>;
  using value_type_matrix_type_host = Kokkos::View<value_type**,Kokkos::LayoutLeft,host_device_type>;
  
  void fill_spd_tridiag_matrix(const value_type_matrix_type_host & A) {
    const int m = A.extent(0), n = A.extent(1);
    EXPECT_TRUE(m == n);
    {
      for (int i=0;i<m;++i)
        A(i,i) = 4.0;

      for (int i=0;i<(m-1);++i) {
        A(i,i+1) = -1.0;
        A(i+1,i) = -1.0;
      }
    }
  }
  void symmetrize_with_upper(const value_type_matrix_type_host & A) {
    const int m = A.extent(0), n = A.extent(1);
    EXPECT_TRUE(m == n);
    for (int i=0;i<m;++i) 
      for (int j=i;i<m;++j)
        A(j,i) = A(i,j);
  }
  void fill_random_matrix(const value_type_matrix_type_host & A) {
    const int m = A.extent(0), n = A.extent(1);
    value_type one;
    atsv::set_real(one, 1);
    atsv::set_imag(one, 1);

    Kokkos::Random_XorShift64_Pool<typename host_device_type::execution_space> random(13718);
    Kokkos::fill_random(A, random, one);
  }
  void fill_random_symmetric_matrix(const value_type_matrix_type_host & A) {
    fill_random_matrix(A);
    symmetrize_with_upper(A);
  }
  void copy_lower_triangular(const value_type_matrix_type_host A, 
                             const value_type diag, const value_type_matrix_type_host L) {
    const int m = L.extent(0), n = L.extent(1);
    {
      const int mA = A.extent(0), nA = A.extent(1);
      EXPECT_TRUE(m == mA);
      EXPECT_TRUE(n == nA);
    }
    const bool replace_diag = diag != value_type(0);
    for (int j=0;j<n;++j)
      for (int i=j;i<m;++i) {
        if (i == j && replace_diag) A(i,j) = diag;
        else A(i,j) = L(i,j);
      }
  } 
}


#include "Tacho_TestCrsMatrixBase.hpp"
#include "Tacho_TestGraphTools.hpp"
#include "Tacho_TestSymbolicTools.hpp"
#include "Tacho_TestDenseLinearAlgebra.hpp"
// //#include "Tacho_TestNumeric.hpp"
// //#include "Tacho_TestTaskFunctor.hpp"

// #include "Tacho_TestDenseMatrixView.hpp"
// #include "Tacho_TestDenseByBlocks.hpp"

// #include "Tacho_TestDenseLinearAlgebra.hpp"

#endif
