#ifndef __TACHO_DENSE_MATRIX_VIEW_HPP__
#define __TACHO_DENSE_MATRIX_VIEW_HPP__

#include "Tacho_Util.hpp"

/// \file Tacho_DenseMatrixView.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

    template<typename ValueType, typename ExecSpace>
    struct DenseMatrixView {
    public:
      enum : ordinal_type { rank = 2 };

      typedef ValueType value_type;
      typedef value_type non_const_value_type;

      typedef ExecSpace execution_space;

      typedef Kokkos::Future<int,execution_space> future_type;

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
      DenseMatrixView(      value_type *buf, 
                      const ordinal_type m,
                      const ordinal_type n)
        : _offm(0), _offn(0), 
          _m(m), _n(n), 
          _rs(1), _cs(m),
          _buf(buf), _future() {}

      KOKKOS_INLINE_FUNCTION
      DenseMatrixView(const DenseMatrixView &b) 
        : _offm(b._offm), _offn(b._offn), 
          _m(b._m), _n(b._n), 
          _rs(b._rs), _cs(b._cs),
          _buf(b._buf), _future() {}

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](const ordinal_type k) const {
        return _buf[k];
      }

      KOKKOS_INLINE_FUNCTION
      value_type& operator()(const ordinal_type i,
                             const ordinal_type j) const {
        return _buf[(i+_offm)*_rs + (j+_offn)*_cs];
      }

      KOKKOS_INLINE_FUNCTION
      void set_view(const DenseMatrixView &base,
                    const ordinal_type offm, const ordinal_type m, 
                    const ordinal_type offn, const ordinal_type n) {
        _rs = base._rs; _cs = base._cs; _buf = base._buf;

        _offm = offm; _m = m; 
        _offn = offn; _n = n;
      }

      KOKKOS_INLINE_FUNCTION
      void set_view(const ordinal_type offm, const ordinal_type m, 
                    const ordinal_type offn, const ordinal_type n) {
        _offm = offm; _m = m; 
        _offn = offn; _n = n;
      }

      KOKKOS_INLINE_FUNCTION
      void set_view(const ordinal_type m, 
                    const ordinal_type n) {
        _offm = 0; _m = m; 
        _offn = 0; _n = n;
      }

      KOKKOS_INLINE_FUNCTION
      void attach_buffer(const ordinal_type rs, 
                         const ordinal_type cs, 
                         const value_type *buf) { 
        _rs = rs; _cs = cs; _buf = const_cast<value_type*>(buf);
      }

      KOKKOS_INLINE_FUNCTION
      void set_future(const future_type &f) { _future = f; }

      KOKKOS_INLINE_FUNCTION
      void set_future() { _future.~future_type(); }
      
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
      value_type* data() const { return _buf+_offm*_rs+_offn*_cs; }

      KOKKOS_INLINE_FUNCTION
      future_type future() const { return _future; } 
    };
 
    template<typename MatrixOfBlocksViewType>
    KOKKOS_INLINE_FUNCTION
    void 
    clearFutureOfBlocks(const MatrixOfBlocksViewType &H) {
      const ordinal_type m = H.dimension_0();
      const ordinal_type n = H.dimension_1();
      for (ordinal_type j=0;j<n;++j)
        for (ordinal_type i=0;i<m;++i)
          H(i,j).set_future();
    }

    template<typename MemberType, 
             typename MatrixOfBlocksViewType>
    KOKKOS_INLINE_FUNCTION
    void 
    clearFutureOfBlocks(/* */ MemberType &member, 
                        const MatrixOfBlocksViewType &H) {
      const ordinal_type m = H.dimension_0();
      const ordinal_type n = H.dimension_1();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,m),[&](const int &i) {
              H(i,j).set_future();
            });
        });
    }
    
    template<typename MatrixOfBlocksViewType>
    KOKKOS_INLINE_FUNCTION
    void 
    setMatrixOfBlocks(const MatrixOfBlocksViewType &H,
                      const ordinal_type m,
                      const ordinal_type n,
                      const ordinal_type mb,
                      const ordinal_type nb) {     
      const ordinal_type bm = H.dimension_0();
      const ordinal_type bn = H.dimension_1();
      
      for (ordinal_type j=0;j<bn;++j) {
        const ordinal_type
          jbeg = j*nb, jtmp = jbeg + nb,
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

    template<typename MatrixOfBlocksViewType>
    KOKKOS_INLINE_FUNCTION
    void 
    setMatrixOfBlocks(const MatrixOfBlocksViewType &H,
                      const ordinal_type m,
                      const ordinal_type n,
                      const ordinal_type mb) {
      setMatrixOfBlocks(H, m, n, mb, mb);
    }
    
    
    template<typename MemberType,
             typename MatrixOfBlocksViewType>
    KOKKOS_INLINE_FUNCTION
    void 
    setMatrixOfBlocks(/* */ MemberType &member, 
                      const MatrixOfBlocksViewType &H,
                      const ordinal_type m,
                      const ordinal_type n,
                      const ordinal_type mb,
                      const ordinal_type nb) {     
      const ordinal_type bm = H.dimension_0();
      const ordinal_type bn = H.dimension_1();
      
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,bn),[&](const int &j) {
          const ordinal_type
            jbeg = j*nb, jtmp = jbeg + nb,
            jend = jtmp > n ? n : jtmp,
            jdiff = (jend > jbeg)*(jend - jbeg);

          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,bm),[&](const int &i) {
              const ordinal_type
                ibeg = i*mb, itmp = ibeg + mb,
                iend = itmp > m ? m : itmp,
                idiff = (iend > ibeg)*(iend - ibeg);
              
              H(i,j).set_view(ibeg, idiff, 
                              jbeg, jdiff);
            });
        });
    }
    
    template<typename MemberType, 
             typename MatrixOfBlocksViewType>
    KOKKOS_INLINE_FUNCTION
    void 
    setMatrixOfBlocks(/* */ MemberType &member, 
                      const MatrixOfBlocksViewType &H,
                      const ordinal_type m,
                      const ordinal_type n,
                      const ordinal_type mb) {
      setMatrixOfBlocks(member, H, m, n, mb, mb);
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

    template<typename MemberType, 
             typename MatrixOfBlocksViewType,
             typename BaseBufferPtrType>
    KOKKOS_INLINE_FUNCTION
    void 
    attachBaseBuffer(/* */ MemberType &member, 
                     const MatrixOfBlocksViewType &H,
                     const BaseBufferPtrType ptr,
                     const ordinal_type rs,
                     const ordinal_type cs) {
      const ordinal_type m = H.dimension_0();
      const ordinal_type n = H.dimension_1();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,m),[&](const int &i) {
              H(i,j).attach_buffer(rs, cs, ptr);
            });
        });
    }
    
    template<typename MatrixOfBlocksViewType,
             typename MemoryPoolType>
    KOKKOS_INLINE_FUNCTION
    void
    allocateStorageByBlocks(const MatrixOfBlocksViewType &H,
                            const MemoryPoolType &pool) {
      typedef typename MatrixOfBlocksViewType::value_type dense_block_type;
      typedef typename dense_block_type::value_type value_type;

      const ordinal_type m = H.dimension_0();
      const ordinal_type n = H.dimension_1();
      for (ordinal_type j=0;j<n;++j)
        for (ordinal_type i=0;i<m;++i) {
          const ordinal_type mm = H(i,j).dimension_0(), nn = H(i,j).dimension_1();
          if (mm > 0 && nn > 0) {
            auto ptr = (value_type*)pool.allocate(mm*nn*sizeof(value_type));
            TACHO_TEST_FOR_ABORT(ptr == NULL, "memory pool allocation fails");          

            H(i,j).set_view(mm, nn); // whatever offsets are defined here, they are gone.
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
          auto &blk = H(i,j);
          const ordinal_type mm = blk.dimension_0(), nn = blk.dimension_1();
          if (mm > 0 && nn > 0) 
            pool.deallocate(blk.data(), mm*nn*sizeof(value_type));
        }
    }

    template<typename ValueType, typename ExecSpace>  
    KOKKOS_INLINE_FUNCTION
    void
    copyElementwise(const DenseMatrixView<ValueType,ExecSpace> &F,
                    const DenseMatrixView<DenseMatrixView<ValueType,ExecSpace>,ExecSpace> &H) {
      const ordinal_type 
        hm = H.dimension_0(), hn = H.dimension_0(),
        fm = F.dimension_0(), fn = F.dimension_0();
        
      if (hm > 0 && hn > 0) {
        ordinal_type offj = 0;
        for (ordinal_type j=0;j<hn;++j) {
          ordinal_type offi = 0;
          for (ordinal_type i=0;i<hm;++i) {
            const auto &blk = H(i,j);
            const ordinal_type 
              mm = blk.dimension_0(), nn = blk.dimension_1();
            for(ordinal_type jj=0;jj<nn;++jj) {
              const ordinal_type jjj = offj+jj;
              for(ordinal_type ii=0;ii<mm;++ii) {
                const ordinal_type iii = offi+ii;
                if (iii < fm && jjj < fn) 
                  F(iii, jjj) = blk(ii,jj);
              }
            }
            offi += mm;
          }
          offj += H(0,j).dimension_1();
        }
      }
    }

    template<typename ValueType, typename ExecSpace>  
    KOKKOS_INLINE_FUNCTION
    void
    copyElementwise(const DenseMatrixView<DenseMatrixView<ValueType,ExecSpace>,ExecSpace> &H,
                    const DenseMatrixView<ValueType,ExecSpace> &F) {
      const ordinal_type 
        hm = H.dimension_0(), hn = H.dimension_0(),
        fm = F.dimension_0(), fn = F.dimension_0();
        
      if (hm > 0 && hn > 0) {
        ordinal_type offj = 0;
        for (ordinal_type j=0;j<hn;++j) {
          ordinal_type offi = 0;
          for (ordinal_type i=0;i<hm;++i) {
            const auto &blk = H(i,j);
            const ordinal_type 
              mm = blk.dimension_0(), nn = blk.dimension_1();
            for(ordinal_type jj=0;jj<nn;++jj) {
              const ordinal_type jjj = offj+jj;
              for(ordinal_type ii=0;ii<mm;++ii) {
                const ordinal_type iii = offi+ii;
                if (iii < fm && jjj < fn) 
                  blk(ii,jj) = F(iii, jjj);
              }
            }
            offi += mm;
          }
          offj += H(0,j).dimension_1();
        }
      }
    }

    template<typename DenseMatrixViewType,
             typename OrdinalTypeArray>
    inline
    void
    applyRowPermutation(const DenseMatrixViewType &A, 
                        const DenseMatrixViewType &B,
                        const OrdinalTypeArray &p) {
      const ordinal_type m = A.dimension_0(), n = A.dimension_1();
      typedef typename DenseMatrixViewType::execution_space execution_space;

      if (true) { //std::is_same<typename execution_space::memory_space,Kokkos::HostSpace>::value) {
        // serial copy on host
        Kokkos::RangePolicy<execution_space,Kokkos::Schedule<Kokkos::Static> > policy(0, m);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ordinal_type &i) {
            for (ordinal_type j=0;j<n;++j)
              A(p(i), j) = B(i, j);
          });
      } else {      
        // gcc has compiler errors
        // Kokkos::TeamPolicy<execution_space,Kokkos::Schedule<Kokkos::Static> > policy(m, 1);
        // Kokkos::parallel_for
        //   (policy, KOKKOS_LAMBDA (const typename Kokkos::TeamPolicy<execution_space>::member_type &member) {
        //     const ordinal_type i = member.league_rank();
        //     Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,n),[&](const int &j) {
        //         Kokkos::single(Kokkos::PerThread(member), [&]() {
        //             A(p(i), j) = B(i, j);
        //           });
        //       });
        //   });
      }
      
    }

}

#endif
