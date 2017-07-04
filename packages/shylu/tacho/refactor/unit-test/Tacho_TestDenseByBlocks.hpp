#ifndef __TACHO_TEST_DENSE_BYBLOCKS_HPP__
#define __TACHO_TEST_DENSE_BYBLOCKS_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_DenseMatrixView.hpp"

#include "TachoExp_Chol_ByBlocks.hpp"

using namespace Tacho::Experimental;

typedef Kokkos::View<ValueType*,HostSpaceType> value_type_array_host;
typedef Kokkos::View<ValueType*,DeviceSpaceType> value_type_array_device;

typedef DenseMatrixView<ValueType,HostSpaceType> DenseMatrixViewHostType;
typedef DenseMatrixView<ValueType,DeviceSpaceType> DenseMatrixViewDeviceType;

typedef DenseMatrixView<DenseMatrixViewHostType,HostSpaceType> DenseMatrixOfBlocksHostType;
typedef DenseMatrixView<DenseMatrixViewHostType,DeviceSpaceType> DenseMatrixOfBlocksDeviceType;

TEST( DenseByBlocks, chol ) {
  const ordinal_type m = 100, mb = 32;

  Kokkos::View<ValueType*,HostSpaceType> a("a", m*m), a1("a1", m*m), a2("a2", m*m);
  DenseMatrixViewHostType A;

  // use tridiagonal matrix for testing
  {
    A.set_view(0, m,
               0, m);

    // make tri diag for testing
    A.attach_buffer(1, m, a.data());
    for (ordinal_type i=0;i<m;++i) {
      A(i,i) = 4;
      const ordinal_type ip = i+1;
      if (ip < m) {
        A(ip,i ) = 1;
        A(i ,ip) = 1;      
      }
    }
    Kokkos::deep_copy(a1, a);
    Kokkos::deep_copy(a2, a);
  }
  
  // referece: lapack chol
  {    
    int dummy;
    Chol<Uplo::Upper,Algo::External>
      ::invoke(dummy, dummy, A);
  }

  // test: chol by blocks with attached base buffer
  const ordinal_type bm = (m/mb) + (m%mb>0);
  Kokkos::View<DenseMatrixViewHostType*,HostSpaceType> h("h", bm*bm);
  DenseMatrixOfBlocksHostType H;

  typedef Kokkos::TaskScheduler<HostSpaceType> sched_type_host;
  sched_type_host sched;
  
  typedef TaskFunctor_Chol<sched_type_host,DenseMatrixOfBlocksHostType,
    Uplo::Upper,Algo::ByBlocks> task_functor_chol;

  const ordinal_type max_functor_size = 4*sizeof(task_functor_chol);
  
  {
    const ordinal_type
      task_queue_capacity = 1024*max_functor_size,
      min_block_size  = 16,
      max_block_size  = 4*max_functor_size,
      num_superblock  = 4,
      superblock_size = task_queue_capacity/num_superblock;
    
    sched = sched_type_host(typename HostSpaceType::memory_space(),
                            task_queue_capacity,
                            min_block_size,
                            max_block_size,
                            superblock_size);
  }

  // compute chol with byblocks - attached buffer
  A.attach_buffer(1, m, a1.data());

  H.set_view(0, bm,
             0, bm);
  
  H.attach_buffer(1, bm, h.data());
  {
    setMatrixOfBlocks(H, m, m, mb);
    attachBaseBuffer(H, A.data(), A.stride_0(), A.stride_1());
    
    Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
                       task_functor_chol(sched, H));
    
    Kokkos::wait(sched);

    clearFutureOfBlocks(H);
  }

  {
    double diff = 0.0, norm = 0.0;
    for (ordinal_type k=0;k<(m*m);++k) {
      norm += a(k)*a(k);
      diff += (a(k) - a1(k))*(a(k) - a1(k));
    }
    
    const double eps = std::numeric_limits<double>::epsilon()*100;
    EXPECT_TRUE(sqrt(diff/norm) < eps);
  }

  // test: chol by blocks with memory pool
  A.attach_buffer(1, m, a2.data());
  {
    const ordinal_type
      capacity = mb*mb*bm*bm*sizeof(ValueType)+1024,
      min_block_size  = mb*mb*sizeof(ValueType),
      max_block_size  = mb*mb*sizeof(ValueType),
      num_superblock  = 1,
      superblock_size = capacity/num_superblock;
    
    Kokkos::MemoryPool<HostSpaceType> pool(typename HostSpaceType::memory_space(),
                                           capacity,
                                           min_block_size,
                                           max_block_size,
                                           superblock_size);
  
    allocateStorageByBlocks(H, pool);
    copyElementwise(H, A);

    Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
                       task_functor_chol(sched, H));
    
    Kokkos::wait(sched);

    clearFutureOfBlocks(H);
    copyElementwise(A, H);
  }

  {
    double diff = 0.0, norm = 0.0;
    for (ordinal_type k=0;k<(m*m);++k) {
      norm += a(k)*a(k);
      diff += (a(k) - a2(k))*(a(k) - a2(k));
    }
    
    const double eps = std::numeric_limits<double>::epsilon()*100;
    EXPECT_TRUE(sqrt(diff/norm) < eps);
  }

}

#endif
