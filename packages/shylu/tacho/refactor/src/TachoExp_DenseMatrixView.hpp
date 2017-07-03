#ifndef __TACHOEXP_DENSE_MATRIX_VIEW_HPP__
#define __TACHOEXP_DENSE_MATRIX_VIEW_HPP__

#include "TachoExp_Util.hpp"

/// \file TachoExp_DenseMatrixView.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  namespace Experimental {

    template<typename ValueType, typename ExecSpace>
    struct DenseMatrixView {
    public:
      typedef ValueType value_type;
      typedef ExecSpace exec_space;

      typedef Kokkos::Future<int,exec_space> future_type;

    private:
      ordinal_type _offm, _offn, _m, _n, _rs, _cs;
      value_type *_buf;      
      future_type _future;

    public:
      KOKKOS_INLINE_FUNCTION
      DenseMatrixView()
        : _offm(0), _offn(0), 
          _m(0), _n(0), 
          _rs(0), _cs(0),
          _buf(NULL), _future() {}

      KOKKOS_INLINE_FUNCTION
      DenseMatrixView(const DenseMatrixView &b) 
        : _offm(b._offm), _offn(b._offn), 
          _m(b._m), _n(b._n), 
          _rs(b._rs), _cs(b._cs),
          _buf(b._buf), _future() {}

      KOKKOS_INLINE_FUNCTION
      value_type& operator()(const ordinal_type i,
                             const ordinal_type j) const {
        return _buf[(i+_offm)*_rs + (j+_offn)*_cs];
      }

      KOKKOS_INLINE_FUNCTION
      void set_view(const ordinal_type offm, const ordinal_type m, 
                    const ordinal_type offn, const ordinal_type n) {
        _offm = offm; _m = m; 
        _offn = offn; _n = n;
      }

      KOKKOS_INLINE_FUNCTION
      void attach_buffer(const ordinal_type rs, 
                         const ordinal_type cs, 
                         const value_type *buf) { 
        _rs = rs; _cs = cs; _buf = const_cast<value_type*>(buf);
      }

      KOKKOS_INLINE_FUNCTION
      void set_future(const future_type &f) { _future = f; }

      /// get methods

      KOKKOS_INLINE_FUNCTION
      ordinal_type offset_0() const { return _offm; }

      KOKKOS_INLINE_FUNCTION
      ordinal_type offset_1() const { return _offn; } 

      KOKKOS_INLINE_FUNCTION
      ordinal_type dimension_0() const { return _m; }

      KOKKOS_INLINE_FUNCTION
      ordinal_type dimension_1() const { return _n; } 

      KOKKOS_INLINE_FUNCTION
      ordinal_type stride_0() const { return _rs; }

      KOKKOS_INLINE_FUNCTION
      ordinal_type stride_1() const { return _cs; }

      KOKKOS_INLINE_FUNCTION
      value_type* data() const { return _buf; }

      KOKKOS_INLINE_FUNCTION
      ordinal_type future() const { return _future; } 
    };

    template<typename MatrixOfBlocksViewType>
    KOKKOS_INLINE_FUNCTION
    void 
    setMatrixOfPartitionedBlocks(const MatrixOfBlocksViewType &H,
                                 const ordinal_type m,
                                 const ordinal_type n,
                                 const ordinal_type mb) {     
      const ordinal_type 
        bm = H.dimension_0(),
        bn = H.dimension_1();
      
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

          H(i,j).set_view(ibeg, idiff, 
                          jbeg, jdiff);
        }
      }
    }

    template<typename MatrixOfBlocksViewType,
             typename BaseBufferPtrType>
    KOKKOS_INLINE_FUNCTION
    void 
    attachBaseBuffer(const MatrixOfBlocksViewType &H,
                     const BaseBufferPtrType ptr,
                     const ordinal_type rs,
                     const ordinal_type cs) {
      const ordinal_type m = H.dimension_0(), n = H.dimension_1();
      for (ordinal_type j=0;j<n;++j)
        for (ordinal_type i=0;i<m;++i) 
          H(i,j).attach_buffer(rs, cs, ptr);
    }
    
    template<typename MatrixOfBlocksViewType,
             typename MemoryPoolType>
    KOKKOS_INLINE_FUNCTION
    void
    allocateStorageByBlocks(const MatrixOfBlocksViewType &H,
                            const MemoryPoolType &pool) {
      typedef typename MatrixOfBlocksViewType::value_type dense_block_type;
      typedef typename dense_block_type::value_type value_type;

      const ordinal_type m = H.dimension_0(), n = H.dimension_1();
      for (ordinal_type j=0;j<n;++j)
        for (ordinal_type i=0;i<m;++i) {
          const ordinal_type mm = H(i,j).dimension_0(), nn = H(i,j).dimension_1();
          if (mm > 0 && nn > 0) {
            auto ptr = (value_type*)pool.allocate(mm*nn*sizeof(value_type));
            TACHO_TEST_FOR_ABORT(ptr == NULL, "memory pool allocation fails");          

            H(i,j).attach_buffer(1, mm, ptr);
          }
        }
    }

    template<typename MatrixOfBlocksViewType,
             typename MemoryPoolType>
    KOKKOS_INLINE_FUNCTION
    void
    deallocateStorageByBlocks(const MatrixOfBlocksViewType &H,
                              const MemoryPoolType &pool) {
      typedef typename MatrixOfBlocksViewType::value_type dense_block_type;
      typedef typename dense_block_type::value_type value_type;

      const ordinal_type m = H.dimension_0(), n = H.dimension_1();
      for (ordinal_type j=0;j<n;++j)
        for (ordinal_type i=0;i<m;++i) {
          const ordinal_type mm = H(i,j).dimension_0(), nn = H(i,j).dimension_1();
          if (mm > 0 && nn > 0) 
            pool.deallocate(H(i,j).data(), mm*nn*sizeof(value_type));
        }
    }
    
  }
}

#endif
