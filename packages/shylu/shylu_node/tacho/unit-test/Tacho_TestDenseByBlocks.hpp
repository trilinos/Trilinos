#ifndef __TACHO_TEST_DENSE_BYBLOCKS_HPP__
#define __TACHO_TEST_DENSE_BYBLOCKS_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_DenseMatrixView.hpp"

#include "TachoExp_Chol_ByBlocks.hpp"
#include "TachoExp_Gemm_ByBlocks.hpp"
#include "TachoExp_Herk_ByBlocks.hpp"
#include "TachoExp_Trsm_ByBlocks.hpp"

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
    A.set_view(m, m);

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

  // temporary testing ldl as dummy
  {
    double a[10][10], work[10][10];
    int ipiv[10], info;
    Lapack<double>::sytrf('U', 10, &a[0][0], 10, &ipiv[0], &work[0][0], 100, &info);
    printf("ldl tested\n");
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

  H.set_view(bm, bm);
  
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
    for (ordinal_type p=0;p<(m*m);++p) {
      norm += real(a(p)*conj(a(p)));
      diff += real((a(p) - a1(p))*conj(a(p) - a1(p)));
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
    for (ordinal_type p=0;p<(m*m);++p) {
      norm += real(a(p)*conj(a(p)));
      diff += real((a(p) - a2(p))*conj(a(p) - a2(p)));
    }
    
    const double eps = std::numeric_limits<double>::epsilon()*100;
    EXPECT_TRUE(sqrt(diff/norm) < eps);
  }

}

TEST( DenseByBlocks, gemm ) {
  double alpha = 2.0, beta = 0.5;
  const ordinal_type m = 100, n = 100, k = 100, mb = 32;

  Kokkos::View<ValueType*,HostSpaceType> a("a", m*k), b("b", k*n), c("c", m*n), c1("c1", m*n);
  DenseMatrixViewHostType A, B, C;

  // use random matrix for testing
  {
    A.set_view(m, k);
    A.attach_buffer(1, m, a.data());

    B.set_view(k, n);
    B.attach_buffer(1, k, b.data());

    C.set_view(m, n);
    C.attach_buffer(1, m, c.data());

    Random<ValueType> random;
    auto randomize = [&](const DenseMatrixViewHostType &mat) {
      const ordinal_type m = mat.dimension_0(), n = mat.dimension_1();
      for (ordinal_type j=0;j<n;++j)
        for (ordinal_type i=0;i<m;++i)
          mat(i,j) = random.value();
    };
    randomize(A);
    randomize(B);
    randomize(C);

    Kokkos::deep_copy(c1, c);
  }
  
  // referece: blas gemm
  {    
    int dummy;
    Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::External>
      ::invoke(dummy, dummy, alpha, A, B, beta, C);

    Gemm<Trans::ConjTranspose,Trans::NoTranspose,Algo::External>
      ::invoke(dummy, dummy, alpha, A, B, beta, C);
  }


  // test: gemm by blocks with attached base buffer
  const ordinal_type 
    bm = (m/mb) + (m%mb>0),
    bn = (n/mb) + (n%mb>0),
    bk = (k/mb) + (k%mb>0);

  Kokkos::View<DenseMatrixViewHostType*,HostSpaceType> ha("ha", bm*bk), hb("hb", bk*bn), hc("hc", bm*bn);
  DenseMatrixOfBlocksHostType HA, HB, HC;

  typedef Kokkos::TaskScheduler<HostSpaceType> sched_type_host;
  sched_type_host sched;
  
  typedef TaskFunctor_Gemm<sched_type_host,double,DenseMatrixOfBlocksHostType,
    Trans::NoTranspose,Trans::NoTranspose,Algo::ByBlocks> task_functor_gemm_nt_nt;

  typedef TaskFunctor_Gemm<sched_type_host,double,DenseMatrixOfBlocksHostType,
    Trans::ConjTranspose,Trans::NoTranspose,Algo::ByBlocks> task_functor_gemm_ct_nt;

  const ordinal_type max_functor_size = 4*sizeof(task_functor_gemm_nt_nt);
  
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

  // compute gemm with byblocks - attached buffer
  C.attach_buffer(1, m, c1.data());

  HA.set_view(bm, bk);
  HB.set_view(bk, bn);
  HC.set_view(bm, bn);
  
  HA.attach_buffer(1, bm, ha.data());
  HB.attach_buffer(1, bk, hb.data());
  HC.attach_buffer(1, bm, hc.data());
  {
    setMatrixOfBlocks(HA, m, k, mb);
    setMatrixOfBlocks(HB, k, n, mb);
    setMatrixOfBlocks(HC, m, n, mb);

    attachBaseBuffer(HA, A.data(), A.stride_0(), A.stride_1());
    attachBaseBuffer(HB, B.data(), B.stride_0(), B.stride_1());
    attachBaseBuffer(HC, C.data(), C.stride_0(), C.stride_1());
    
    Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
                       task_functor_gemm_nt_nt(sched, alpha, HA, HB, beta, HC));
    Kokkos::wait(sched);

    Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
                       task_functor_gemm_ct_nt(sched, alpha, HA, HB, beta, HC));
    Kokkos::wait(sched);

    clearFutureOfBlocks(HC);
  }


  {
    double diff = 0.0, norm = 0.0;
    for (ordinal_type p=0;p<(m*m);++p) {
      norm += real(c(p)*conj(c(p)));
      diff += real((c(p) - c1(p))*conj(c(p) - c1(p)));
    }
    
    const double eps = std::numeric_limits<double>::epsilon()*100;
    EXPECT_TRUE(sqrt(diff/norm) < eps);
  }
}

TEST( DenseByBlocks, herk ) {
  double alpha = 2.0, beta = 0.5;
  const ordinal_type n = 100, k = 100, mb = 32; 

  Kokkos::View<ValueType*,HostSpaceType> a("a", k*n), c("c", n*n), c1("c1", n*n);
  DenseMatrixViewHostType A, C;

  // use random matrix for testing
  {
    A.set_view(k, n);
    A.attach_buffer(1, k, a.data());

    C.set_view(n, n);
    C.attach_buffer(1, n, c.data());

    Random<ValueType> random;
    auto randomize = [&](const DenseMatrixViewHostType &mat) {
      const ordinal_type m = mat.dimension_0(), n = mat.dimension_1();
      for (ordinal_type j=0;j<n;++j)
        for (ordinal_type i=0;i<m;++i)
          mat(i,j) = random.value();
    };
    randomize(A);
    randomize(C);

    Kokkos::deep_copy(c1, c);
  }
  
  // referece: blas gemm
  {    
    int dummy;
    Herk<Uplo::Upper,Trans::ConjTranspose,Algo::External>
      ::invoke(dummy, dummy, alpha, A, beta, C);
  }


  // test: herk by blocks with attached base buffer
  const ordinal_type 
    bn = (n/mb) + (n%mb>0),
    bk = (k/mb) + (k%mb>0);

  Kokkos::View<DenseMatrixViewHostType*,HostSpaceType> ha("ha", bk*bn), hc("hc", bn*bn);
  DenseMatrixOfBlocksHostType HA, HC;

  typedef Kokkos::TaskScheduler<HostSpaceType> sched_type_host;
  sched_type_host sched;
  
  typedef TaskFunctor_Herk<sched_type_host,double,DenseMatrixOfBlocksHostType,
    Uplo::Upper,Trans::ConjTranspose,Algo::ByBlocks> task_functor_herk_u_ct;

  const ordinal_type max_functor_size = 4*sizeof(task_functor_herk_u_ct);
  
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

  // compute gemm with byblocks - attached buffer
  C.attach_buffer(1, n, c1.data());

  HA.set_view(bk, bn);
  HC.set_view(bn, bn);
  
  HA.attach_buffer(1, bk, ha.data());
  HC.attach_buffer(1, bn, hc.data());
  {
    setMatrixOfBlocks(HA, k, n, mb);
    setMatrixOfBlocks(HC, n, n, mb);

    attachBaseBuffer(HA, A.data(), A.stride_0(), A.stride_1());
    attachBaseBuffer(HC, C.data(), C.stride_0(), C.stride_1());
    
    Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
                       task_functor_herk_u_ct(sched, alpha, HA, beta, HC));
    Kokkos::wait(sched);

    clearFutureOfBlocks(HC);
  }


  {
    double diff = 0.0, norm = 0.0;
    for (ordinal_type p=0;p<(n*n);++p) {
      norm += real(c(p)*conj(c(p)));
      diff += real((c(p) - c1(p))*conj(c(p) - c1(p)));
    }
    
    const double eps = std::numeric_limits<double>::epsilon()*100;
    EXPECT_TRUE(sqrt(diff/norm) < eps);
  }
}


TEST( DenseByBlocks, trsm ) {
  double alpha = 2.0;
  const ordinal_type m = 4, n = 4, mb = 4; 

  Kokkos::View<ValueType*,HostSpaceType> a("a", m*m), b("c", m*n), b1("c1", m*n);
  DenseMatrixViewHostType A, B;

  // use random matrix for testing
  {
    A.set_view(m, m);
    A.attach_buffer(1, m, a.data());

    B.set_view(m, n);
    B.attach_buffer(1, m, b.data());

    Random<ValueType> random;
    auto randomize = [&](const DenseMatrixViewHostType &mat) {
      const ordinal_type m = mat.dimension_0(), n = mat.dimension_1();
      for (ordinal_type j=0;j<n;++j)
        for (ordinal_type i=0;i<m;++i)
          mat(i,j) = random.value();
    };
    randomize(A);
    randomize(B);

    Kokkos::deep_copy(b1, b);
  }

  // referece: blas trsm
  {    
    int dummy;
    Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
      ::invoke(dummy, dummy, Diag::NonUnit(), alpha, A, B);
    Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Algo::External>
      ::invoke(dummy, dummy, Diag::NonUnit(), alpha, A, B);
  }

  // test: trsm by blocks with attached base buffer
  const ordinal_type 
    bm = (m/mb) + (m%mb>0),
    bn = (n/mb) + (n%mb>0);

  Kokkos::View<DenseMatrixViewHostType*,HostSpaceType> ha("ha", bm*bm), hb("hb", bm*bn);
  DenseMatrixOfBlocksHostType HA, HB;

  typedef Kokkos::TaskScheduler<HostSpaceType> sched_type_host;
  sched_type_host sched;
  
  typedef TaskFunctor_Trsm<sched_type_host,double,DenseMatrixOfBlocksHostType,
    Side::Left,Uplo::Upper,Trans::ConjTranspose,Diag::NonUnit,Algo::ByBlocks> 
    task_functor_trsm_l_u_ct_nd;
  typedef TaskFunctor_Trsm<sched_type_host,double,DenseMatrixOfBlocksHostType,
    Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Algo::ByBlocks> 
    task_functor_trsm_l_u_nt_nd;

  const ordinal_type max_functor_size = 4*sizeof(task_functor_trsm_l_u_ct_nd);
  
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

  // compute gemm with byblocks - attached buffer
  B.attach_buffer(1, m, b1.data());
  
  HA.set_view(bm, bm);
  HB.set_view(bm, bn);
  
  HA.attach_buffer(1, bm, ha.data());
  HB.attach_buffer(1, bm, hb.data());
  {
    setMatrixOfBlocks(HA, m, m, mb);
    setMatrixOfBlocks(HB, m, n, mb);

    attachBaseBuffer(HA, A.data(), A.stride_0(), A.stride_1());
    attachBaseBuffer(HB, B.data(), B.stride_0(), B.stride_1());
    
    Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
                       task_functor_trsm_l_u_ct_nd(sched, alpha, HA, HB));
    Kokkos::wait(sched);
    
    Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
                       task_functor_trsm_l_u_nt_nd(sched, alpha, HA, HB));
    Kokkos::wait(sched);

    clearFutureOfBlocks(HB);
  }

  {
    double diff = 0.0, norm = 0.0;
    for (ordinal_type p=0;p<(m*n);++p) {
      norm += real(b(p)*conj(b(p)));
      diff += real((b(p) - b1(p))*conj(b(p) - b1(p)));
    }
    
    const double eps = std::numeric_limits<double>::epsilon()*100;
    EXPECT_TRUE(sqrt(diff/norm) < eps);
  }
}


#endif
