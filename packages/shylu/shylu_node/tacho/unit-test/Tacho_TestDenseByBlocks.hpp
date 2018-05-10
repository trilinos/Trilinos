#ifndef __TACHO_TEST_DENSE_BYBLOCKS_HPP__
#define __TACHO_TEST_DENSE_BYBLOCKS_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_DenseMatrixView.hpp"

#include "Tacho_Chol_ByBlocks.hpp"
#include "Tacho_Gemm_ByBlocks.hpp"
#include "Tacho_Herk_ByBlocks.hpp"
#include "Tacho_Trsm_ByBlocks.hpp"

using namespace Tacho;

typedef Kokkos::View<ValueType*,HostSpaceType>   value_type_array_host;
typedef Kokkos::View<ValueType*,DeviceSpaceType> value_type_array;
typedef Kokkos::TaskScheduler<DeviceSpaceType>   sched_type;

typedef ArithTraits<ValueType> ats;

typedef DenseMatrixView<ValueType,HostSpaceType>   DenseMatrixViewTypeHost;
typedef DenseMatrixView<ValueType,DeviceSpaceType> DenseMatrixViewType;

typedef DenseMatrixView<DenseMatrixViewType,HostSpaceType>   DenseMatrixOfBlocksTypeHost;
typedef DenseMatrixView<DenseMatrixViewType,DeviceSpaceType> DenseMatrixOfBlocksType;


TEST( DenseByBlocks, ldl ) {
  TEST_BEGIN;
  // dummy for compiler test
  double aa[10][10], work[10][10];
  int ipiv[10], info;
  Lapack<double>::sytrf('U', 10, &aa[0][0], 10, &ipiv[0], &work[0][0], 100, &info);
  TEST_END;
}
#if 1
TEST( DenseByBlocks, chol ) {
  TEST_BEGIN;
  const ordinal_type m = 100, mb = 32;
  
  // a  : referece with lapack
  // a1 : byblocks with partitioned matrices
  Kokkos::DualView<ValueType*,DeviceSpaceType> a("a", m*m), a1("a1", m*m);
  
  // reference lapack 
  {
    DenseMatrixViewTypeHost A;
    a.modify<HostSpaceType>();

    A.set_view(m, m);
    A.attach_buffer(1, m, a.h_view.data());

    // make tri diag for testing
    for (ordinal_type i=0;i<m;++i) {
      A(i,i) = 4;
      const ordinal_type ip = i+1;
      if (ip < m) {
        A(ip,i ) = 1;
        A(i ,ip) = 1;      
      }
    }

    a1.modify<DeviceSpaceType>();
    
    Kokkos::deep_copy(a1.d_view, a.h_view);
  
    Chol<Uplo::Upper,Algo::External>::invoke(A);
  }

  // test: chol by blocks with attached base buffer
  {
    const ordinal_type bm = (m/mb) + (m%mb>0);

    DenseMatrixViewType A;
    A.set_view(m, m);
    A.attach_buffer(1, m, a1.d_view.data());

    a1.sync<DeviceSpaceType>();

    Kokkos::DualView<DenseMatrixViewType*,DeviceSpaceType> h("h", bm*bm);
    h.modify<HostSpaceType>();

    DenseMatrixOfBlocksTypeHost H;
    H.set_view(bm, bm);  
    H.attach_buffer(1, bm, h.h_view.data());

    DenseMatrixOfBlocksType D;
    D.set_view(bm, bm);  
    D.attach_buffer(1, bm, h.d_view.data());

    typedef TaskFunctor_Chol<sched_type,DenseMatrixOfBlocksType,
                             Uplo::Upper,
                             Algo::ByBlocks> TaskFunctorChol;
    
    const ordinal_type max_functor_size = 4*sizeof(TaskFunctorChol);
    const ordinal_type
      task_queue_capacity = 1024*max_functor_size,
      min_block_size  = 16,
      max_block_size  = 4*max_functor_size,
      num_superblock  = 4,
      superblock_size = task_queue_capacity/num_superblock;
    
    sched_type sched(typename sched_type::memory_space(),
                     task_queue_capacity,
                     min_block_size,
                     max_block_size,
                     superblock_size);

    setMatrixOfBlocks(H, m, m, mb);
    attachBaseBuffer(H, A.data(), A.stride_0(), A.stride_1());

    h.sync<DeviceSpaceType>();

    
    Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
                       TaskFunctorChol(sched, D));
    Kokkos::wait(sched);

    a1.sync<HostSpaceType>();
    clearFutureOfBlocks(H);
  }

  // check
  {
    double diff = 0.0, norm = 0.0;
    for (ordinal_type p=0;p<(m*m);++p) {
      norm += ats::abs(a.h_view(p)*ats::conj(a.h_view(p)));
      diff += ats::abs((a.h_view(p) - a1.h_view(p))*ats::conj(a.h_view(p) - a1.h_view(p)));
    }
    const double eps = std::numeric_limits<double>::epsilon()*100;
    EXPECT_TRUE(sqrt(diff/norm) < eps);
  }
  TEST_END;
}
#endif
TEST( DenseByBlocks, gemm ) {
  TEST_BEGIN;

  double alpha = 2.0, beta = 0.5;
  const ordinal_type m = 100, n = 100, k = 100, mb = 32;

  // c  : result from reference blas
  // c1 : result from partitioned matrices
  Kokkos::DualView<ValueType*,DeviceSpaceType> a("a", m*k), b("b", k*n), c("c", m*n), c1("c1", m*n);


  // reference blas
  {
    DenseMatrixViewTypeHost A, B, C;

    a.modify<HostSpaceType>();
    b.modify<HostSpaceType>();
    c.modify<HostSpaceType>();

    A.set_view(m, k);
    A.attach_buffer(1, m, a.h_view.data());

    B.set_view(k, n);
    B.attach_buffer(1, k, b.h_view.data());

    C.set_view(m, n);
    C.attach_buffer(1, m, c.h_view.data());

    Random<ValueType> random;
    auto randomize = [&](const DenseMatrixViewTypeHost &mat) {
      const ordinal_type mm = mat.dimension_0(), nn = mat.dimension_1();
      for (ordinal_type j=0;j<nn;++j)
        for (ordinal_type i=0;i<mm;++i)
          mat(i,j) = random.value();
    };

    randomize(A);
    randomize(B);
    randomize(C);

    c1.modify<DeviceSpaceType>();

    Kokkos::deep_copy(c1.d_view, c.h_view);

    Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::External>
      ::invoke(alpha, A, B, beta, C);

    Gemm<Trans::ConjTranspose,Trans::NoTranspose,Algo::External>
      ::invoke(alpha, A, B, beta, C);
  }
  

  // test: gemm by blocks with attached base buffer
  {
    // compute gemm with byblocks - attached buffer
    DenseMatrixViewType A, B, C;
    
    A.set_view(m, k);
    A.attach_buffer(1, m, a.d_view.data());

    B.set_view(k, n);
    B.attach_buffer(1, k, b.d_view.data());

    C.set_view(m, n);
    C.attach_buffer(1, m, c1.d_view.data());

    a.sync<DeviceSpaceType>();
    b.sync<DeviceSpaceType>();
    c1.sync<DeviceSpaceType>();

    const ordinal_type 
      bm = (m/mb) + (m%mb>0),
      bn = (n/mb) + (n%mb>0),
      bk = (k/mb) + (k%mb>0);
    
    Kokkos::DualView<DenseMatrixViewType*,DeviceSpaceType> ha("ha", bm*bk), hb("hb", bk*bn), hc("hc", bm*bn);

    DenseMatrixOfBlocksTypeHost HA, HB, HC;

    HA.set_view(bm, bk);
    HB.set_view(bk, bn);
    HC.set_view(bm, bn);
    
    HA.attach_buffer(1, bm, ha.h_view.data());
    HB.attach_buffer(1, bk, hb.h_view.data());
    HC.attach_buffer(1, bm, hc.h_view.data());
    
    DenseMatrixOfBlocksType DA, DB, DC;    

    DA.set_view(bm, bk);
    DB.set_view(bk, bn);
    DC.set_view(bm, bn);
    
    DA.attach_buffer(1, bm, ha.d_view.data());
    DB.attach_buffer(1, bk, hb.d_view.data());
    DC.attach_buffer(1, bm, hc.d_view.data());
    
    typedef TaskFunctor_Gemm<sched_type,double,DenseMatrixOfBlocksType,
                             Trans::NoTranspose,Trans::NoTranspose,
                             Algo::ByBlocks> TaskFunctorGemm_NT_NT;
    
    typedef TaskFunctor_Gemm<sched_type,double,DenseMatrixOfBlocksType,
                             Trans::ConjTranspose,Trans::NoTranspose,
                             Algo::ByBlocks> TaskFunctorGemm_CT_NT;

    const ordinal_type max_functor_size = 4*sizeof(TaskFunctorGemm_NT_NT);
    const ordinal_type
      task_queue_capacity = 1024*max_functor_size,
      min_block_size  = 16,
      max_block_size  = 4*max_functor_size,
      num_superblock  = 4,
      superblock_size = task_queue_capacity/num_superblock;
    
    sched_type sched(typename sched_type::memory_space(),
                     task_queue_capacity,
                     min_block_size,
                     max_block_size,
                     superblock_size);

    ha.modify<HostSpaceType>();
    hb.modify<HostSpaceType>();
    hc.modify<HostSpaceType>();

    setMatrixOfBlocks(HA, m, k, mb);
    setMatrixOfBlocks(HB, k, n, mb);
    setMatrixOfBlocks(HC, m, n, mb);

    attachBaseBuffer(HA, A.data(), A.stride_0(), A.stride_1());
    attachBaseBuffer(HB, B.data(), B.stride_0(), B.stride_1());
    attachBaseBuffer(HC, C.data(), C.stride_0(), C.stride_1());

    ha.sync<DeviceSpaceType>();
    hb.sync<DeviceSpaceType>();
    hc.sync<DeviceSpaceType>();
    
    Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
                       TaskFunctorGemm_NT_NT(sched, alpha, DA, DB, beta, DC));
    Kokkos::wait(sched);

    Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
                       TaskFunctorGemm_CT_NT(sched, alpha, DA, DB, beta, DC));
    Kokkos::wait(sched);

    c1.sync<HostSpaceType>();
    clearFutureOfBlocks(HC);
  }

  // check
  {
    double diff = 0.0, norm = 0.0;
    for (ordinal_type p=0;p<(m*m);++p) {
      norm += ats::abs(c.h_view(p)*ats::conj(c.h_view(p)));
      diff += ats::abs((c.h_view(p) - c1.h_view(p))*ats::conj(c.h_view(p) - c1.h_view(p)));
    }
    
    const double eps = std::numeric_limits<double>::epsilon()*100;
    EXPECT_TRUE(sqrt(diff/norm) < eps);
  }
  TEST_END;
}

TEST( DenseByBlocks, herk ) {
  TEST_BEGIN;

  double alpha = 2.0, beta = 0.5;
  const ordinal_type n = 100, k = 50, mb = 32;

  // c  : result from reference blas
  // c1 : result from byblocks
  Kokkos::DualView<ValueType*,DeviceSpaceType> a("a", k*n), c("c", n*n), c1("c1", n*n);

  // referece: blas herk
  {
    DenseMatrixViewTypeHost A, C;

    a.modify<HostSpaceType>();
    c.modify<HostSpaceType>();

    A.set_view(k, n);
    A.attach_buffer(1, k, a.h_view.data());

    C.set_view(n, n);
    C.attach_buffer(1, n, c.h_view.data());
    
    Random<ValueType> random;
    auto randomize = [&](const DenseMatrixViewTypeHost &mat) {
      const ordinal_type mm = mat.dimension_0(), nn = mat.dimension_1();
      for (ordinal_type j=0;j<nn;++j)
        for (ordinal_type i=0;i<mm;++i)
          mat(i,j) = random.value();
    };

    randomize(A);
    randomize(C);

    c1.modify<DeviceSpaceType>();

    Kokkos::deep_copy(c1.d_view, c.h_view);

    Herk<Uplo::Upper,Trans::ConjTranspose,Algo::External>
      ::invoke(alpha, A, beta, C);
  }


  // test: herk by blocks with attached base buffer
  {
    DenseMatrixViewType A, C;
    
    A.set_view(k, n);
    A.attach_buffer(1, k, a.d_view.data());

    C.set_view(n, n);
    C.attach_buffer(1, n, c1.d_view.data());

    a.sync<DeviceSpaceType>();
    c1.sync<DeviceSpaceType>();

    const ordinal_type 
      bn = (n/mb) + (n%mb>0),
      bk = (k/mb) + (k%mb>0);
    
    Kokkos::DualView<DenseMatrixViewType*,DeviceSpaceType> ha("ha", bk*bn), hc("hc", bn*bn);

    DenseMatrixOfBlocksTypeHost HA, HC;

    HA.set_view(bk, bn);
    HC.set_view(bn, bn);
    
    HA.attach_buffer(1, bk, ha.h_view.data());
    HC.attach_buffer(1, bn, hc.h_view.data());
    
    DenseMatrixOfBlocksType DA, DC;    
    
    DA.set_view(bk, bn);
    DC.set_view(bn, bn);
    
    DA.attach_buffer(1, bk, ha.d_view.data());
    DC.attach_buffer(1, bn, hc.d_view.data());
    
    typedef TaskFunctor_Herk<sched_type,double,DenseMatrixOfBlocksType,
                             Uplo::Upper,Trans::ConjTranspose,
                             Algo::ByBlocks> TaskFunctorHerk_U_CT;

    const ordinal_type max_functor_size = 4*sizeof(TaskFunctorHerk_U_CT);
    const ordinal_type
      task_queue_capacity = 1024*max_functor_size,
      min_block_size  = 16,
      max_block_size  = 4*max_functor_size,
      num_superblock  = 4,
      superblock_size = task_queue_capacity/num_superblock;
    
    sched_type sched(typename sched_type::memory_space(),
                     task_queue_capacity,
                     min_block_size,
                     max_block_size,
                     superblock_size);
    
    ha.modify<HostSpaceType>();
    hc.modify<HostSpaceType>();

    setMatrixOfBlocks(HA, k, n, mb);
    setMatrixOfBlocks(HC, n, n, mb);

    attachBaseBuffer(HA, A.data(), A.stride_0(), A.stride_1());
    attachBaseBuffer(HC, C.data(), C.stride_0(), C.stride_1());
    
    ha.sync<DeviceSpaceType>();
    hc.sync<DeviceSpaceType>();
    c1.modify<DeviceSpaceType>();

    Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
                       TaskFunctorHerk_U_CT(sched, alpha, DA, beta, DC));
    Kokkos::wait(sched);

    c1.sync<HostSpaceType>();
    clearFutureOfBlocks(HC);
  }

  // check: C is hermitian matrix should be forced to be zero imaginary 
  {
    double diff = 0.0, norm = 0.0;
    for (ordinal_type i=0;i<n;++i) {
      for (ordinal_type j=0;j<n;++j) {
        const ordinal_type p = i*n+j;
        if (i == j) {
          norm += ats::abs(ats::real(c.h_view(p))*ats::real(c.h_view(p)));
          diff += ats::abs(ats::real(c.h_view(p) - c1.h_view(p))*ats::real(c.h_view(p) - c1.h_view(p)));
        } else {
          norm += ats::abs(c.h_view(p)*ats::conj(c.h_view(p)));
          diff += ats::abs((c.h_view(p) - c1.h_view(p))*ats::conj(c.h_view(p) - c1.h_view(p)));
        }
      }
    }

    const double eps = std::numeric_limits<double>::epsilon()*100;
    EXPECT_TRUE(sqrt(diff/norm) < eps);
  }
  TEST_END;
}

TEST( DenseByBlocks, trsm ) {
  TEST_BEGIN;

  double alpha = 2.0;
  const ordinal_type m = 4, n = 4, mb = 4; 

  // b  : result from reference blas
  // b1 : result from byblocks
  Kokkos::DualView<ValueType*,DeviceSpaceType> a("a", m*m), b("c", m*n), b1("c1", m*n);

  // reference blas
  {
    DenseMatrixViewTypeHost A, B;

    a.modify<HostSpaceType>();
    b.modify<HostSpaceType>();

    A.set_view(m, m);
    A.attach_buffer(1, m, a.h_view.data());

    B.set_view(m, n);
    B.attach_buffer(1, m, b.h_view.data());

    Random<ValueType> random;
    auto randomize = [&](const DenseMatrixViewTypeHost &mat) {
      const ordinal_type mm = mat.dimension_0(), nn = mat.dimension_1();
      for (ordinal_type j=0;j<nn;++j)
        for (ordinal_type i=0;i<mm;++i)
          mat(i,j) = random.value();
    };

    randomize(A);
    randomize(B);

    b1.modify<DeviceSpaceType>();

    Kokkos::deep_copy(b1.d_view, b.h_view);
    
    Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
      ::invoke(Diag::NonUnit(), alpha, A, B);
    Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Algo::External>
      ::invoke(Diag::NonUnit(), alpha, A, B);
  }

  // test: trsm by blocks with attached base buffer
  {
    DenseMatrixViewType A, B;    

    A.set_view(m, m);
    A.attach_buffer(1, m, a.d_view.data());

    B.set_view(m, n);
    B.attach_buffer(1, m, b1.d_view.data());

    a.sync<DeviceSpaceType>();
    b1.sync<DeviceSpaceType>();
    
    const ordinal_type 
      bm = (m/mb) + (m%mb>0),
      bn = (n/mb) + (n%mb>0);

    Kokkos::DualView<DenseMatrixViewType*,DeviceSpaceType> ha("ha", bm*bm), hb("hb", bm*bn);

    DenseMatrixOfBlocksTypeHost HA, HB;

    HA.set_view(bm, bm);
    HB.set_view(bm, bn);
    
    HA.attach_buffer(1, bm, ha.h_view.data());
    HB.attach_buffer(1, bm, hb.h_view.data());

    DenseMatrixOfBlocksType DA, DB;

    DA.set_view(bm, bm);
    DB.set_view(bm, bn);
    
    DA.attach_buffer(1, bm, ha.d_view.data());
    DB.attach_buffer(1, bm, hb.d_view.data());
    
    typedef TaskFunctor_Trsm<sched_type,double,DenseMatrixOfBlocksType,
                             Side::Left,Uplo::Upper,Trans::ConjTranspose,Diag::NonUnit,
                             Algo::ByBlocks> TaskFunctorTrsm_L_U_CT_ND;
    typedef TaskFunctor_Trsm<sched_type,double,DenseMatrixOfBlocksType,
                             Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,
                             Algo::ByBlocks> TaskFunctorTrsm_L_U_NT_ND;
    
    const ordinal_type max_functor_size = 4*sizeof(TaskFunctorTrsm_L_U_CT_ND);
    const ordinal_type
      task_queue_capacity = 1024*max_functor_size,
      min_block_size  = 16,
      max_block_size  = 4*max_functor_size,
      num_superblock  = 4,
      superblock_size = task_queue_capacity/num_superblock;
    
    sched_type sched(typename sched_type::memory_space(),
                     task_queue_capacity,
                     min_block_size,
                     max_block_size,
                     superblock_size);

    ha.modify<HostSpaceType>();
    hb.modify<HostSpaceType>();

    setMatrixOfBlocks(HA, m, m, mb);
    setMatrixOfBlocks(HB, m, n, mb);

    attachBaseBuffer(HA, A.data(), A.stride_0(), A.stride_1());
    attachBaseBuffer(HB, B.data(), B.stride_0(), B.stride_1());

    ha.sync<DeviceSpaceType>();
    hb.sync<DeviceSpaceType>();
    
    Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
                       TaskFunctorTrsm_L_U_CT_ND(sched, alpha, DA, DB));
    Kokkos::wait(sched);
    
    Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
                       TaskFunctorTrsm_L_U_NT_ND(sched, alpha, DA, DB));
    Kokkos::wait(sched);

    b1.sync<HostSpaceType>();
    clearFutureOfBlocks(HB);
  }
  
  // check
  {
    double diff = 0.0, norm = 0.0;
    for (ordinal_type p=0;p<(m*n);++p) {
      norm += ats::abs(b.h_view(p)*ats::conj(b.h_view(p)));
      diff += ats::abs((b.h_view(p) - b1.h_view(p))*ats::conj(b.h_view(p) - b1.h_view(p)));
    }
    
    const double eps = std::numeric_limits<double>::epsilon()*100;
    EXPECT_TRUE(sqrt(diff/norm) < eps);
  }
  TEST_END;
}


#endif








//
// the following test is experimental code that I want to test with storage by blocks 
// aided by a memory pool
//

// // test: chol by blocks with memory pool
// A.attach_buffer(1, m, a2.data());
// {
//   const ordinal_type bm = (m/mb) + (m%mb>0);
//   const ordinal_type
//     capacity = mb*mb*bm*bm*sizeof(ValueType)+1024,
//     min_block_size  = mb*mb*sizeof(ValueType),
//     max_block_size  = mb*mb*sizeof(ValueType),
//     num_superblock  = 1,
//     superblock_size = capacity/num_superblock;
    
//   Kokkos::MemoryPool<HostSpaceType> pool(typename HostSpaceType::memory_space(),
//                                          capacity,
//                                          min_block_size,
//                                          max_block_size,
//                                          superblock_size);
  
//   allocateStorageByBlocks(H, pool);
//   copyElementwise(H, A);

//   Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
//                      task_functor_chol(sched, H));
//   Kokkos::wait(sched);

//   clearFutureOfBlocks(H);
//   copyElementwise(A, H);
// }

// {
//   double diff = 0.0, norm = 0.0;
//   for (ordinal_type p=0;p<(m*m);++p) {
//     norm += ats::abs(a(p)*ats::conj(a(p)));
//     diff += ats::abs((a(p) - a2(p))*ats::conj(a(p) - a2(p)));
//   }
    
//   const double eps = std::numeric_limits<double>::epsilon()*100;
//   EXPECT_TRUE(sqrt(diff/norm) < eps);
// }
