#ifndef __TACHO_TEST_DENSE_MATRIX_VIEW_HPP__
#define __TACHO_TEST_DENSE_MATRIX_VIEW_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_DenseMatrixView.hpp"

using namespace Tacho::Experimental;

typedef Kokkos::View<ValueType*,HostSpaceType> value_type_array_host;
typedef Kokkos::View<ValueType*,DeviceSpaceType> value_type_array_device;

typedef DenseMatrixView<ValueType,HostSpaceType> DenseMatrixViewHostType;
typedef DenseMatrixView<ValueType,DeviceSpaceType> DenseMatrixViewDeviceType;

typedef DenseMatrixView<DenseMatrixViewHostType,HostSpaceType> DenseMatrixOfBlocksHostType;
typedef DenseMatrixView<DenseMatrixViewHostType,DeviceSpaceType> DenseMatrixOfBlocksDeviceType;

TEST( DenseMatrixView, flat ) {
  const ordinal_type m = 10, n = 10;

  Kokkos::View<ValueType*,HostSpaceType> buf("buf", m*n);
  DenseMatrixViewHostType A;

  A.set_view(0, m,
             0, n);
  
  A.attach_buffer(1, m, buf.data());

  {
    ordinal_type cnt = 0;
    for (ordinal_type j=0;j<n;++j)  
      for (ordinal_type i=0;i<m;++i)
        A(i,j) = cnt++;
  }

  for (ordinal_type k=0;k<(m*n);++k) {
    EXPECT_TRUE(buf[k] == k);
  }  
}

TEST( DenseMatrixView, hier ) {
  const ordinal_type m = 5, n = 5, mb = 3;

  Kokkos::View<ValueType*,HostSpaceType> abuf("abuf", m*n);
  DenseMatrixViewHostType A;

  A.set_view(0, m,
             0, n);
  
  A.attach_buffer(1, m, abuf.data());

  {
    ordinal_type cnt = 0;
    for (ordinal_type j=0;j<n;++j)  
      for (ordinal_type i=0;i<m;++i)
        A(i,j) = cnt++;
  }

  const ordinal_type bm = (m/mb) + (m%mb>0), bn = (n/mb) + (n%mb>0);

  Kokkos::View<DenseMatrixViewHostType*,HostSpaceType> hbuf("hbuf", bm*bn);  
  DenseMatrixOfBlocksHostType H;

  H.set_view(0, bm,
             0, bn);
  
  H.attach_buffer(1, bm, hbuf.data());

  setMatrixOfPartitionedBlocks(H, m, n, mb);
  attachBaseBuffer(H, A.data(), A.stride_0(), A.stride_1());

  for (ordinal_type k=0;k<(bm*bn);++k) {
    auto &blk = hbuf[k];
    printf("block k = %d, offset (%d,%d), dim (%d,%d), stride (%d,%d)\n",
           k, 
           blk.offset_0(),blk.offset_1(),
           blk.dimension_0(),blk.dimension_1(),
           blk.stride_0(),blk.stride_1());
    
    if (blk.dimension_0() > 0 && blk.dimension_1()) {
      for (ordinal_type i=0;i<blk.dimension_0();++i) {
        for (ordinal_type j=0;j<blk.dimension_1();++j)
          printf(" %4.2f ", blk(i,j));
        printf("\n");
      }
    } else {
      printf("this block is empty\n");
    }
  }  
}

TEST( DenseMatrixView, memorypool ) {
  const ordinal_type m = 5, n = 5, mb = 3;

  Kokkos::View<ValueType*,HostSpaceType> abuf("abuf", m*n);
  DenseMatrixViewHostType A;

  A.set_view(0, m,
             0, n);
  
  A.attach_buffer(1, m, abuf.data());

  {
    ordinal_type cnt = 0;
    for (ordinal_type j=0;j<n;++j)  
      for (ordinal_type i=0;i<m;++i)
        A(i,j) = cnt++;
  }

  const ordinal_type bm = (m/mb) + (m%mb>0), bn = (n/mb) + (n%mb>0);

  Kokkos::View<DenseMatrixViewHostType*,HostSpaceType> hbuf("hbuf", bm*bn);  
  DenseMatrixOfBlocksHostType H;

  H.set_view(0, bm,
             0, bn);
  
  H.attach_buffer(1, bm, hbuf.data());

  setMatrixOfPartitionedBlocks(H, m, n, mb);

  Kokkos::MemoryPool<HostSpaceType> pool(typename HostSpaceType::memory_space(),
                                         1024*sizeof(ValueType),
                                         16*sizeof(ValueType),
                                         1024*sizeof(ValueType),
                                         512*sizeof(ValueType));
  allocateStorageByBlocks(H, pool);


  ordinal_type cnt = 0;
  for (ordinal_type j=0;j<bn;++j) 
    for (ordinal_type i=0;i<bm;++i) {
      auto &blk = H(i,j);
      for (ordinal_type ii=0;ii<blk.dimension_0();++ii)
        for (ordinal_type jj=0;jj<blk.dimension_1();++jj)
          blk(ii,jj) = cnt++;
    }

  for (ordinal_type k=0;k<(bm*bn);++k) {
    auto &blk = hbuf[k];
    printf("block k = %d, offset (%d,%d), dim (%d,%d), stride (%d,%d)\n",
           k, 
           blk.offset_0(),blk.offset_1(),
           blk.dimension_0(),blk.dimension_1(),
           blk.stride_0(),blk.stride_1());
    
    if (blk.dimension_0() > 0 && blk.dimension_1()) {
      for (ordinal_type i=0;i<blk.dimension_0();++i) {
        for (ordinal_type j=0;j<blk.dimension_1();++j)
          printf(" %4.2f ", blk(i,j));
        printf("\n");
      }
    } else {
      printf("this block is empty\n");
    }
  }  

  deallocateStorageByBlocks(H, pool);

  
}


#endif
