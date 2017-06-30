#ifndef __TACHOEXP_DENSE_BLOCK_HPP__
#define __TACHOEXP_DENSE_BLOCK_HPP__

#include "TachoExp_Util.hpp"

/// \file TachoExp_DenseBlock.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  namespace Experimental {

    template<typename ValueType, typename ExecSpace>
    struct DenseBlock {
      typedef ValueType value_type;
      typedef ExecSpace exec_space;

      typedef Kokkos::Future<int,exec_space> future_type;

      ordinal_type m, n, cs, rs;
      value_type *buf;      
      future_type future;

      KOKKOS_INLINE_FUNCTION
      value_type& operator()(const ordinal_type i,
                             const ordinal_type j) const {
        return buf[i*rs + j*cs];
      }      
    };

    KOKKOS_INLINE_FUNCTION
    template<typename MatrixOfBlocksViewType>
    setMatrixOfBlocks(const MatrixOfBlocksViewType &H,
                      const ordinal_type m,
                      const ordinal_type n,
                      const ordinal_type mb) {     
      const ordinal_type 
        bm = H.dimension_0(),
        bn = H.dimension_1();

      DenseBlock blk;
      for (ordinal_type j=0;j<bn;++j) {
        const ordinal_type
          jbeg = j*mb, jtmp = jbeg + mb,
          jend = jtmp > n ? n : jtmp,
          jdiff = (jend > jbeg)*(jend - jbeg);

        for (ordinal_type i=0;i<bm;++i) {
          const ordinal_type
            ibeg = i*mb, itmp = ibeg + mb,
            iend = itmp > m ? m : itmp,
            idiff = (iend > ibeg)*(iend - ibeg);

          H(i,j).m = idiff;
          H(i,j).n = jdiff;
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    template<typename MatrixOfBlocksViewType>
    attachBaseBuffer(const MatrixOfBlocksViewType &H,
                     const void *ptr,
                     const ordinal_type cs,
                     const ordinal_type rs) {
      const ordinal_type m = H.dimension_0(), n = H.dimension_1();
      for (ordinal_type j=0;j<n;++j)
        for (ordinal_type i=0;i<m;++i) {
          auto &blk = H(i,j);
          blk.buf = ptr;
          blk.cs  = cs;
          blk.rs  = rs;
        }
    }

    KOKKOS_INLINE_FUNCTION
    template<typename MatrixOfBlocksViewType,
             typename MemoryPoolType>
    allocateStorageByBlocks(const MatrixOfBlocksViewType &H,
                            const MemoryPoolType &pool) {
      typedef typename MatrixOfBlocksViewType::non_const_value_type dense_block_type;
      typedef typename dense_block_type::value_type value_type;

      const ordinal_type m = H.dimension_0(), n = H.dimension_1();
      for (ordinal_type j=0;j<n;++j)
        for (ordinal_type i=0;i<m;++i) {
          auto &blk = H(i,j);
          if (blk.m > 0 && blk.n > 0) {
            blk.buf = (value_type*)pool.allocate(blk.m*blk.n*sizeof(value_type));
            TACHO_TEST_FOR_ABORT(blk.buf == NULL, "memory pool allocation fails");          

            blk.cs  = m;
            blk.rs  = 1;
          }
        }
    }

    KOKKOS_INLINE_FUNCTION
    template<typename MatrixOfBlocksViewType,
             typename MemoryPoolType>
    deallocateStorageByBlocks(const MatrixOfBlocksViewType &H,
                              const MemoryPoolType &pool) {
      typedef typename MatrixOfBlocksViewType::non_const_value_type dense_block_type;
      typedef typename dense_block_type::value_type value_type;

      dense_block_type empty;
      const ordinal_type m = H.dimension_0(), n = H.dimension_1();
      for (ordinal_type j=0;j<n;++j)
        for (ordinal_type i=0;i<m;++i) {
          auto &blk = H(i,j);
          if (blk.m > 0 && blk.n > 0) 
            pool.deallocate(blk.buf, blk.m*blk.n*sizeof(value_type));
          blk = empty;
        }
    }
    
  }
}

#endif
