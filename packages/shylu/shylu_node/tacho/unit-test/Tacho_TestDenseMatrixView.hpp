#ifndef __TACHO_TEST_DENSE_MATRIX_VIEW_HPP__
#define __TACHO_TEST_DENSE_MATRIX_VIEW_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_DenseMatrixView.hpp"

using namespace Tacho;

typedef Kokkos::View<ValueType*,HostDeviceType> value_type_array_host;
//typedef Kokkos::View<ValueType*,DeviceType> value_type_array_device;

typedef TaskSchedulerType<typename HostDeviceType::execution_space> host_scheduler_type;
typedef TaskSchedulerType<typename DeviceType::execution_space> scheduler_type;

typedef DenseMatrixView<ValueType,host_scheduler_type> DenseMatrixViewHostType;
//typedef DenseMatrixView<ValueType,scheduler_type> DenseMatrixViewDeviceType;

typedef DenseMatrixView<DenseMatrixViewHostType,host_scheduler_type> DenseMatrixOfBlocksHostType;
//typedef DenseMatrixView<DenseMatrixViewHostType,scheduler_type> DenseMatrixOfBlocksDeviceType;

TEST( DenseMatrixView, flat ) {
  TEST_BEGIN;
  const ordinal_type m = 10, n = 10;

  Kokkos::View<ValueType*,HostDeviceType> a("a", m*n);
  DenseMatrixViewHostType A;

  A.set_view(0, m,
             0, n);
  
  A.attach_buffer(1, m, a.data());

  {
    ordinal_type cnt = 0;
    for (ordinal_type j=0;j<n;++j)  
      for (ordinal_type i=0;i<m;++i)
        A(i,j) = cnt++;
  }

  for (ordinal_type k=0;k<(m*n);++k) {
    EXPECT_TRUE(a[k] == ValueType(k));
  }  
  TEST_END;
}

TEST( DenseMatrixView, hier ) {
  TEST_BEGIN;
  const ordinal_type m = 5, n = 5, mb = 3;

  Kokkos::View<ValueType*,HostDeviceType> a("a", m*n), a1("a1", m*n);
  DenseMatrixViewHostType A;

  A.set_view(0, m,
             0, n);
  
  A.attach_buffer(1, m, a.data());

  {
    ordinal_type cnt = 0;
    for (ordinal_type j=0;j<n;++j)  
      for (ordinal_type i=0;i<m;++i)
        A(i,j) = cnt++;
  }

  const ordinal_type bm = (m/mb) + (m%mb>0), bn = (n/mb) + (n%mb>0);

  Kokkos::View<DenseMatrixViewHostType*,HostDeviceType> hbuf("hbuf", bm*bn);  
  DenseMatrixOfBlocksHostType H;

  H.set_view(0, bm,
             0, bn);
  
  H.attach_buffer(1, bm, hbuf.data());

  setMatrixOfBlocks(H, m, n, mb);
  attachBaseBuffer(H, A.data(), A.stride_0(), A.stride_1());

  DenseMatrixViewHostType A1;
  A1.set_view(0, m,
             0, n);
  
  A1.attach_buffer(1, m, a1.data());

  copyElementwise(A1, H);

  for (ordinal_type k=0;k<(m*n);++k) 
    EXPECT_EQ(a(k), a1(k));
  TEST_END;
}

TEST( DenseMatrixView, memorypool ) {
  TEST_BEGIN;
  const ordinal_type m = 5, n = 5, mb = 3;

  Kokkos::View<ValueType*,HostDeviceType> a("a", m*n), a1("a1", m*n);
  DenseMatrixViewHostType A;

  A.set_view(0, m,
             0, n);
  
  A.attach_buffer(1, m, a.data());

  {
    ordinal_type cnt = 0;
    for (ordinal_type j=0;j<n;++j)  
      for (ordinal_type i=0;i<m;++i)
        A(i,j) = cnt++;
  }

  const ordinal_type bm = (m/mb) + (m%mb>0), bn = (n/mb) + (n%mb>0);

  Kokkos::View<DenseMatrixViewHostType*,HostDeviceType> hbuf("hbuf", bm*bn);  
  DenseMatrixOfBlocksHostType H;

  H.set_view(0, bm,
             0, bn);
  
  H.attach_buffer(1, bm, hbuf.data());

  setMatrixOfBlocks(H, m, n, mb);

  Kokkos::MemoryPool<HostDeviceType> pool(typename HostDeviceType::memory_space(),
                                         1024*sizeof(ValueType),
                                         1*sizeof(ValueType),
                                         1024*sizeof(ValueType),
                                         1024*sizeof(ValueType));
  allocateStorageByBlocks(H, pool);

  DenseMatrixViewHostType A1;
  A1.set_view(0, m,
              0, n);
  
  A1.attach_buffer(1, m, a1.data());

  copyElementwise(H, A);
  copyElementwise(A1, H);

  deallocateStorageByBlocks(H, pool);

  for (ordinal_type k=0;k<(m*n);++k) 
    EXPECT_EQ(a(k), a1(k));  
  TEST_END;
}


#endif
