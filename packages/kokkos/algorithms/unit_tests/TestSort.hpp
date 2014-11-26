/*
 * TestSort.hpp
 *
 *  Created on: Nov 25, 2014
 *      Author: crtrott
 */

#ifndef TESTSORT_HPP_
#define TESTSORT_HPP_

#include <gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<Kokkos_Sort.hpp>

namespace Test {

namespace Impl{

template<class ExecutionSpace, class Scalar>
struct is_sorted_struct {
  typedef unsigned int value_type;
  typedef ExecutionSpace execution_space;

  Kokkos::View<Scalar*,ExecutionSpace> keys;

  is_sorted_struct(Kokkos::View<Scalar*,ExecutionSpace> keys_):keys(keys_) {}
  KOKKOS_INLINE_FUNCTION
  void operator() (int i, unsigned int& count) const {
    if(keys(i)>keys(i+1)) count++;
  }
};

template<class ExecutionSpace, class Scalar>
struct sum {
  typedef double value_type;
  typedef ExecutionSpace execution_space;

  Kokkos::View<Scalar*,ExecutionSpace> keys;

  sum(Kokkos::View<Scalar*,ExecutionSpace> keys_):keys(keys_) {}
  KOKKOS_INLINE_FUNCTION
  void operator() (int i, double& count) const {
    count+=keys(i);
  }
};

template<class ExecutionSpace, class Scalar>
struct bin3d_is_sorted_struct {
  typedef unsigned int value_type;
  typedef ExecutionSpace execution_space;

  Kokkos::View<Scalar*[3],ExecutionSpace> keys;

  int max_bins;
  Scalar min;
  Scalar max;

  bin3d_is_sorted_struct(Kokkos::View<Scalar*[3],ExecutionSpace> keys_,int max_bins_,Scalar min_,Scalar max_):
    keys(keys_),max_bins(max_bins_),min(min_),max(max_) {
  }
  KOKKOS_INLINE_FUNCTION
  void operator() (int i, unsigned int& count) const {
    int ix1 = int ((keys(i,0)-min)/max * max_bins);
    int iy1 = int ((keys(i,1)-min)/max * max_bins);
    int iz1 = int ((keys(i,2)-min)/max * max_bins);
    int ix2 = int ((keys(i+1,0)-min)/max * max_bins);
    int iy2 = int ((keys(i+1,1)-min)/max * max_bins);
    int iz2 = int ((keys(i+1,2)-min)/max * max_bins);

    if (ix1>ix2)  count++;
    else if(ix1==ix2) {
      if (iy1>iy2)  count++;
      else if ((iy1==iy2) && (iz1>iz2))  count++;
    }
  }
};

template<class ExecutionSpace, class Scalar>
struct sum3D {
  typedef double value_type;
  typedef ExecutionSpace execution_space;

  Kokkos::View<Scalar*[3],ExecutionSpace> keys;

  sum3D(Kokkos::View<Scalar*[3],ExecutionSpace> keys_):keys(keys_) {}
  KOKKOS_INLINE_FUNCTION
  void operator() (int i, double& count) const {
    count+=keys(i,0);
    count+=keys(i,1);
    count+=keys(i,2);
  }
};

template<class ExecutionSpace, typename KeyType>
void test_1D_sort(unsigned int n,bool force_kokkos) {
  typedef Kokkos::View<KeyType*,ExecutionSpace> KeyViewType;
  typedef typename KeyViewType::memory_space::size_type size_type;
  KeyViewType keys("Keys",n);

  Kokkos::Random_XorShift64_Pool<ExecutionSpace> g(1931);
  Kokkos::fill_random(keys,g,Kokkos::Random_XorShift64_Pool<ExecutionSpace>::generator_type::MAX_URAND);

  double sum_before = 0.0;
  double sum_after = 0.0;
  unsigned int sort_fails = 0;

  Kokkos::parallel_reduce(n,sum<ExecutionSpace, KeyType>(keys),sum_before);

  Kokkos::sort(keys,force_kokkos);

  Kokkos::parallel_reduce(n,sum<ExecutionSpace, KeyType>(keys),sum_after);
  Kokkos::parallel_reduce(n-1,is_sorted_struct<ExecutionSpace, KeyType>(keys),sort_fails);

  double ratio = sum_before/sum_after;
  double epsilon = 1e-10;
  unsigned int equal_sum = (ratio > (1.0-epsilon)) && (ratio < (1.0+epsilon)) ? 1 : 0;

  ASSERT_EQ(sort_fails,0);
  ASSERT_EQ(equal_sum,1);
}

template<class ExecutionSpace, typename KeyType>
void test_3D_sort(unsigned int n) {
  typedef Kokkos::View<KeyType*[3],ExecutionSpace > KeyViewType;
  typedef typename KeyViewType::memory_space::size_type size_type;

  KeyViewType keys("Keys",n*n*n);

  Kokkos::Random_XorShift64_Pool<ExecutionSpace> g(1931);
  Kokkos::fill_random(keys,g,100.0);

  double sum_before = 0.0;
  double sum_after = 0.0;
  unsigned int sort_fails = 0;

  Kokkos::parallel_reduce(keys.dimension_0(),sum3D<ExecutionSpace, KeyType>(keys),sum_before);

  int bin_1d = 1;
  while(bin_1d*bin_1d*bin_1d*4<keys.dimension_0()) bin_1d*=2;
  int bin_max[3] = {bin_1d,bin_1d,bin_1d};
  typename KeyViewType::value_type min[3] = {0,0,0};
  typename KeyViewType::value_type max[3] = {100,100,100};

  typedef Kokkos::SortImpl::DefaultBinOp3D< KeyViewType > BinOp;
  BinOp bin_op(bin_max,min,max);
  Kokkos::BinSort< KeyViewType , BinOp >
    Sorter(keys,bin_op,false);
  Sorter.create_permute_vector();
  Sorter.template sort< KeyViewType >(keys);

  Kokkos::parallel_reduce(keys.dimension_0(),sum3D<ExecutionSpace, KeyType>(keys),sum_after);
  Kokkos::parallel_reduce(keys.dimension_0()-1,bin3d_is_sorted_struct<ExecutionSpace, KeyType>(keys,bin_1d,min[0],max[0]),sort_fails);

  double ratio = sum_before/sum_after;
  double epsilon = 1e-10;
  unsigned int equal_sum = (ratio > (1.0-epsilon)) && (ratio < (1.0+epsilon)) ? 1 : 0;

  printf("3D Sort Sum: %lf %lf Fails: %u\n",sum_before,sum_after,sort_fails);
  ASSERT_EQ(sort_fails,0);
  ASSERT_EQ(equal_sum,1);
}

template<class ExecutionSpace, typename KeyType>
void test_sort(unsigned int N)
{
  test_1D_sort<ExecutionSpace,KeyType>(N*N*N, true);
  test_1D_sort<ExecutionSpace,KeyType>(N*N*N, false);
  test_3D_sort<ExecutionSpace,KeyType>(N);
}

}
}
#endif /* TESTSORT_HPP_ */
