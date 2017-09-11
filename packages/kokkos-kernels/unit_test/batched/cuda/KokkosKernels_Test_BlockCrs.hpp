#include <cassert>
#include <vector>
#include <algorithm>

#if defined(__KOKKOSKERNELS_INTEL_MKL__)
#include "mkl.h"
#endif

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosKernels_Util.hpp"
#include "KokkosKernels_Vector.hpp"

#include "KokkosKernels_Gemv_Decl.hpp"
#include "KokkosKernels_Gemv_Serial_Impl.hpp"

#include "KokkosKernels_Trsv_Decl.hpp"
#include "KokkosKernels_Trsv_Serial_Impl.hpp"

#include "KokkosKernels_Gemm_Decl.hpp"
#include "KokkosKernels_Gemm_Serial_Impl.hpp"
//#include "KokkosKernels_Gemm_Team_Impl.hpp"

#include "KokkosKernels_Trsm_Decl.hpp"
#include "KokkosKernels_Trsm_Serial_Impl.hpp"

#include "KokkosKernels_LU_Decl.hpp"
#include "KokkosKernels_LU_Serial_Impl.hpp"

#include "KokkosKernels_Test_BlockCrs_Util.hpp"

namespace KokkosKernels {
  
  namespace Test {

    template <typename ExecSpace, typename ArrayLayout>
    class BlockCrsMatrixVectorProductByRow {
    public:
      typedef BlockCrsMatrix<ExecSpace,ArrayLayout> block_crs_matrix_type;
      typedef typename block_crs_matrix_type::crs_graph_type crs_graph_type;
      typedef BlockMultiVector<ExecSpace,ArrayLayout> block_multi_vector_type;
      
    private:
      ConstUnmanagedViewType<typename crs_graph_type::row_ptr_type> _rowptr;
      ConstUnmanagedViewType<typename crs_graph_type::col_idx_type> _colidx;

      ConstUnmanagedViewType<typename block_crs_matrix_type::value_array_type> _A;
      ConstUnmanagedViewType<typename block_multi_vector_type::value_array_type> _x;
      /**/ UnmanagedViewType<typename block_multi_vector_type::value_array_type> _y;

      ordinal_type _blocksize;

    public:
      // A thread maps to a point row of the matrix.
      // loop = blksize*m
      KOKKOS_INLINE_FUNCTION 
      void operator()(const ordinal_type idx) const {
        // index of blockrow and row in a block
        const ordinal_type i  = idx/_blocksize;
        const ordinal_type ii = idx%_blocksize;
        
        // loop over multivectors
        const ordinal_type jend = _y.dimension_0();
        for (ordinal_type j=0;j<jend;++j) {
          scalar_type tmp = 0;
          
          // block row 
          const ordinal_type 
            cbegin = _rowptr(i), cend = _rowptr(i+1);

          for (ordinal_type c=cbegin;c<cend;++c) {
            const ordinal_type col = _colidx(c);
            for (ordinal_type jj=0;jj<_blocksize;++jj) 
              tmp += _A(col,ii,jj)*_x(j, col, jj);
          }
          _y(j, i, ii) = tmp;
        }
      }
      
      void run(const block_crs_matrix_type A,
               const block_multi_vector_type x, 
               const block_multi_vector_type y) {
        _rowptr = A.CrsGraph().rowptr;
        _colidx = A.CrsGraph().colidx;

        _blocksize = A.BlockSize();

        _A = A.Values();
        _x = x.Values();
        _y = y.Values();

        Kokkos::parallel_for(_x.dimension_1()*_blocksize, *this);
      }
    };

    template <typename ExecSpace, typename ArrayLayout>
    class BlockCrsMatrixVectorProductByBlockRow {
    public:
      typedef BlockCrsMatrix<ExecSpace,ArrayLayout> block_crs_matrix_type;
      typedef typename block_crs_matrix_type::crs_graph_type crs_graph_type;
      typedef BlockMultiVector<ExecSpace,ArrayLayout> block_multi_vector_type;

    private:
      ConstUnmanagedViewType<typename crs_graph_type::row_ptr_type> _rowptr;
      ConstUnmanagedViewType<typename crs_graph_type::col_idx_type> _colidx;

      ConstUnmanagedViewType<typename block_crs_matrix_type::value_array_type> _A;
      ConstUnmanagedViewType<typename block_multi_vector_type::value_array_type> _x;
      /**/ UnmanagedViewType<typename block_multi_vector_type::value_array_type> _y;

      ordinal_type _blocksize;

    public:
      
      // A thread maps to a row block of the matrix.
      // loop = m
      KOKKOS_INLINE_FUNCTION 
      void operator()(const ordinal_type i) const {
        // loop over multivector colums
        const ordinal_type jend = _y.dimension_0();
        for (ordinal_type j=0;j<jend;++j) {
          // set zero
          for (ordinal_type ii=0;ii<_blocksize;++ii) 
            _y(j, i, ii) = 0;
          
          // block row 
          const ordinal_type 
            cbegin = _rowptr(i), cend = _rowptr(i+1);
          
          for (ordinal_type c=cbegin;c<cend;++c) {
            const ordinal_type col = _colidx(c);
            for (ordinal_type ii=0;ii<_blocksize;++ii) {
              scalar_type tmp = 0;
              for (ordinal_type jj=0;jj<_blocksize;++jj) 
                tmp += _A(col,ii,jj)*_x(j, col, jj);
              _y(j, i, ii) += tmp;
            }
          }
        }
      }
      
      void run(const block_crs_matrix_type A,
               const block_multi_vector_type x, 
               const block_multi_vector_type y) {
        _rowptr = A.CrsGraph().rowptr;
        _colidx = A.CrsGraph().colidx;

        _blocksize = A.BlockSize();

        _A = A.Values();
        _x = x.Values();
        _y = y.Values();

        Kokkos::parallel_for(_x.dimension_1(), *this);
      }
    };

    template <typename ExecSpace, typename ValueType, typename ArrayLayout>
    class ExtractBlockTridiagMatrices {
    public:
      typedef ExecSpace exec_space;
      typedef ValueType value_type;
      typedef ArrayLayout array_layout;

      typedef StructuredBlock structured_block_mesh_type;
      typedef BlockCrsMatrix<exec_space,array_layout> block_crs_matrix_type;
      typedef typename block_crs_matrix_type::crs_graph_type crs_graph_type;
      typedef BlockTridiagMatrices<exec_space,value_type,array_layout> block_tridiag_matrices_type;

    private:
      structured_block_mesh_type _mesh;
      ordinal_type _blocksize;
      
      ConstUnmanagedViewType<typename crs_graph_type::row_ptr_type> _rowptr;
      ConstUnmanagedViewType<typename crs_graph_type::row_idx_type> _rowidx;
      ConstUnmanagedViewType<typename crs_graph_type::col_idx_type> _colidx;
      
      ConstUnmanagedViewType<typename block_crs_matrix_type::value_array_type> _A;
      /**/ UnmanagedViewType<typename block_tridiag_matrices_type::value_array_type> _TA, _TB, _TC;
      
    public:
      ExtractBlockTridiagMatrices(const structured_block_mesh_type mesh) 
        : _mesh(mesh) {}

      template<typename TViewType,
               typename AViewType>
      KOKKOS_INLINE_FUNCTION
      void 
      elementwise_copy(const TViewType &T,
                       const AViewType &A,
                       const ordinal_type ij, 
                       const ordinal_type k,
                       const ordinal_type c,
                       const ordinal_type blocksize) const {
        for (ordinal_type ii=0;ii<blocksize;++ii)
          for (ordinal_type jj=0;jj<blocksize;++jj) 
            tdiag_val(T, ij, k, ii, jj) = A(c, ii, jj);
      }
      
      // A thread maps nonzero blocks
      KOKKOS_INLINE_FUNCTION 
      void operator()(const ordinal_type c) const {
        const ordinal_type row = _rowidx[c], col = _colidx[c];

        ordinal_type ri, rj, rk, ci, cj, ck;
        _mesh.id2ijk(row, ri, rj, rk);
        _mesh.id2ijk(col, ci, cj, ck);
  
        if (ri == ci && rj == cj) {
          const ordinal_type ij = _mesh.ij2id(ri, rj);
          // consider connectivity to k-direction
          switch (rk - ck) {
          case  1: elementwise_copy(_TC, _A, ij, ck, c, _blocksize); break;
          case  0: elementwise_copy(_TA, _A, ij, rk, c, _blocksize); break;
          case -1: elementwise_copy(_TB, _A, ij, rk, c, _blocksize); break;
          }
        }
      }
      
      void run(const block_crs_matrix_type A,
               const block_tridiag_matrices_type T) { 
        _rowptr = A.CrsGraph().rowptr;
        _rowidx = A.CrsGraph().rowidx;
        _colidx = A.CrsGraph().colidx;

        _A = A.Values();

        _TA = T.A(); 
        _TB = T.B(); 
        _TC = T.C();

        _blocksize = A.BlockSize();
        Kokkos::parallel_for(_A.dimension_0(), *this);
      }

      template<typename TViewType,
               typename AViewType>
      bool elementwise_check(const TViewType &T,
                             const AViewType &A,
                             const ordinal_type ij, 
                             const ordinal_type k,
                             const ordinal_type c,
                             const ordinal_type blocksize) const {
        const auto eps = 1e2*std::numeric_limits<scalar_type>::epsilon();
        for (ordinal_type ii=0;ii<blocksize;++ii)
          for (ordinal_type jj=0;jj<blocksize;++jj) 
            if ( std::abs(tdiag_val(T, ij, k, ii, jj) - A(c, ii, jj)) >= eps ) return false; 
        return true;
      }
      
      bool check() const {        
        auto rowptr = Kokkos::create_mirror_view(_rowptr); Kokkos::deep_copy(rowptr, _rowptr);
        auto colidx = Kokkos::create_mirror_view(_colidx); Kokkos::deep_copy(colidx, _colidx);
        auto TA     = Kokkos::create_mirror_view(_TA);     Kokkos::deep_copy(TA, _TA);
        auto TB     = Kokkos::create_mirror_view(_TB);     Kokkos::deep_copy(TB, _TB);
        auto TC     = Kokkos::create_mirror_view(_TC);     Kokkos::deep_copy(TC, _TC);
        auto A      = Kokkos::create_mirror_view(_A);      Kokkos::deep_copy(A, _A);

        const ordinal_type 
          ijend = adjustDimension<value_type>(_mesh.ni*_mesh.nj),
          kend = _mesh.nk;

        assert(ijend == TA.dimension_0()); assert((kend - 0) == TA.dimension_1()); 
        assert(ijend == TB.dimension_0()); assert((kend - 1) == TB.dimension_1());
        assert(ijend == TC.dimension_0()); assert((kend - 1) == TC.dimension_1());

        for (ordinal_type ij=0;ij<ijend;++ij) {
          ordinal_type i, j;
          _mesh.id2ij(ij, i, j);

          for (ordinal_type k=0;k<kend;++k) {
            const ordinal_type row = _mesh.ijk2id(i, j, k),
              idx_begin = rowptr[row], 
              idx_end = rowptr[row+1];

            // check
            bool found[3] = {}, same[3] = {};
            for (ordinal_type idx=idx_begin;idx<idx_end;++idx) {
              switch (row - colidx[idx]) {
              case  1: same[2] = elementwise_check(TC, A, ij, k-1, idx, _blocksize); found[2] = true; break;
              case  0: same[0] = elementwise_check(TA, A, ij, k,   idx, _blocksize); found[0] = true; break;
              case -1: same[1] = elementwise_check(TB, A, ij, k,   idx, _blocksize); found[1] = true; break;
              }
            }
            if      (k == 0)         assert(found[0] & same[0] && found[1] & same[1]); 
            else if (k == (kend-1))  assert(found[0] & same[0] && found[2] & same[2]); 
            else                     assert(found[0] & same[0] && found[1] & same[1] && found[2] & same[2]); 
          }
        }            
        return true;
      }
    };

    template <typename ExecSpace, typename ValueType, typename ArrayLayout,
              typename LU_AlgoTagType,
              typename Trsm_AlgoTagType,
              typename Gemm_AlgoTagType>
    class FactorizeBlockTridiagMatrices {
    public:
      typedef ExecSpace exec_space;
      typedef ValueType value_type;
      typedef ArrayLayout array_layout;

      typedef BlockTridiagMatrices<exec_space,value_type,array_layout> block_tridiag_matrices_type;

    private:
      ordinal_type _ntridiag, _m, _blocksize;

      UnmanagedViewType<typename block_tridiag_matrices_type::value_array_type> _TA, _TB, _TC;

    public:
      FactorizeBlockTridiagMatrices() {}

      // A thread maps nonzero blocks
      KOKKOS_INLINE_FUNCTION 
      void operator()(const ordinal_type ij) const {
        auto A = Kokkos::subview(_TA, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto B = Kokkos::subview(_TB, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto C = Kokkos::subview(_TC, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

        const ordinal_type kend = _m - 1;
        for (ordinal_type k=0;k<kend;++k) {
          auto AA = Kokkos::subview(A, k,   Kokkos::ALL(), Kokkos::ALL());
          auto BB = Kokkos::subview(B, k,   Kokkos::ALL(), Kokkos::ALL());
          auto CC = Kokkos::subview(C, k,   Kokkos::ALL(), Kokkos::ALL());
          auto DD = Kokkos::subview(A, k+1, Kokkos::ALL(), Kokkos::ALL());

          Serial::LU<LU_AlgoTagType>
            ::invoke(AA);
          Serial::Trsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsm_AlgoTagType>
            ::invoke(1.0, AA, BB);
          Serial::Trsm<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsm_AlgoTagType>
            ::invoke(1.0, AA, CC);
          Serial::Gemm<Trans::NoTranspose,Trans::NoTranspose,Gemm_AlgoTagType>
            ::invoke(-1.0, CC, BB, 1.0, DD);
        }
        auto AA = Kokkos::subview(A, kend, Kokkos::ALL(), Kokkos::ALL());
        Serial::LU<LU_AlgoTagType>
          ::invoke(AA);
      }

      double FlopCount(const block_tridiag_matrices_type T) {
        const int 
          ntridiag = T.NumTridiagMatrices(),
          m = T.NumRows(),
          blocksize = T.BlockSize();
        
        return ntridiag*( (m-1)*(LU_FlopCount(blocksize, blocksize) + 
                                 Trsm_Lower_FlopCountLower(blocksize, blocksize) +
                                 Trsm_Upper_FlopCountUpper(blocksize, blocksize) +
                                 Gemm_FlopCount(blocksize, blocksize, blocksize)) +
                          LU_FlopCount(blocksize, blocksize) );
      }

      // for batched blas check
      void run(const block_tridiag_matrices_type T, const bool fake = false) { 
        _ntridiag = T.NumTridiagMatrices();
        _m = T.NumRows(); 
        _blocksize = T.BlockSize();

        _TA = T.A(); 
        _TB = T.B(); 
        _TC = T.C();

        // parallel over the instances of tridiagonal matrices
        if (!fake)
          Kokkos::parallel_for(_ntridiag, *this);
      }

      template<typename AViewType, typename LViewType, typename UViewType>
      void a_subtract_mult_l_and_u(const ordinal_type tl, const ordinal_type il, AViewType A,
                                   const ordinal_type tr, const ordinal_type ir, LViewType L, UViewType U) {
        for (ordinal_type ii=0;ii<_blocksize;++ii)
          for (ordinal_type jj=0;jj<_blocksize;++jj)
            for (ordinal_type kk=0;kk<_blocksize;++kk) {
              const auto l = ( ii == kk ? 1 :
                               ii >  kk ? tdiag_val(L, tr, ir, ii, kk) : 0 );
              const auto u = ( kk <= jj ? tdiag_val(U, tr, ir, kk, jj) : 0 );
              tdiag_val(A, tl, il, ii, jj) -= l*u;
            }
      }

      template<typename AViewType, typename BViewType, typename CViewType>
      void a_subtract_mult_b_and_c(const ordinal_type tl, const ordinal_type il, AViewType A,
                                   const ordinal_type tr, const ordinal_type ir, BViewType B, CViewType C) {
        for (ordinal_type ii=0;ii<_blocksize;++ii)
          for (ordinal_type jj=0;jj<_blocksize;++jj)
            for (ordinal_type kk=0;kk<_blocksize;++kk) 
              tdiag_val(A, tl, il, ii, jj) -= ( tdiag_val(B, tr, ir, ii, kk)*
                                                tdiag_val(C, tr, ir, kk, jj) );
      }

      template<typename AViewType, typename LViewType, typename BViewType>
      void a_subtract_mult_l_and_b(const ordinal_type tl, const ordinal_type il, AViewType A,
                                   const ordinal_type tr, const ordinal_type ir, LViewType L, BViewType B) {
        for (ordinal_type ii=0;ii<_blocksize;++ii)
          for (ordinal_type jj=0;jj<_blocksize;++jj)
            for (ordinal_type kk=0;kk<_blocksize;++kk) {
              const auto l = ( ii == kk ? 1.0 :
                               ii >  kk ? tdiag_val(L, tr, ir, ii, kk) : 0 );
              tdiag_val(A, tl, il, ii, jj) -= l*tdiag_val(B, tr, ir, kk, jj);
            }
      }
      
      template<typename AViewType, typename BViewType, typename UViewType>
      void a_subtract_mult_b_and_u(const ordinal_type tl, const ordinal_type il, AViewType A,
                                   const ordinal_type tr, const ordinal_type ir, BViewType B, UViewType U) {
        for (ordinal_type ii=0;ii<_blocksize;++ii)
          for (ordinal_type jj=0;jj<_blocksize;++jj)
            for (ordinal_type kk=0;kk<_blocksize;++kk) {
              const auto u = ( kk <= jj ? tdiag_val(U, tr, ir, kk, jj) : 0 );
              tdiag_val(A, tl, il, ii, jj) -= tdiag_val(B, tr, ir, ii, kk)*u;
            }
      }

      bool check(const block_tridiag_matrices_type T) {
        // factors
        auto DD = Kokkos::create_mirror_view(_TA);    Kokkos::deep_copy(DD, _TA);
        auto UU = Kokkos::create_mirror_view(_TB);    Kokkos::deep_copy(UU, _TB);
        auto LL = Kokkos::create_mirror_view(_TC);    Kokkos::deep_copy(LL, _TC);

        // input A
        auto A = Kokkos::create_mirror_view(T.A());   Kokkos::deep_copy(A, T.A());
        auto B = Kokkos::create_mirror_view(T.B());   Kokkos::deep_copy(B, T.B());
        auto C = Kokkos::create_mirror_view(T.C());   Kokkos::deep_copy(C, T.C());
        
        // diffs
        Kokkos::View<value_type****,Kokkos::DefaultHostExecutionSpace> 
          AA("AA", _ntridiag, _m,   _blocksize, _blocksize),
          BB("BB", _ntridiag, _m-1, _blocksize, _blocksize),
          CC("CC", _ntridiag, _m-1, _blocksize, _blocksize);
        
        Kokkos::deep_copy(AA, A);
        Kokkos::deep_copy(BB, B);
        Kokkos::deep_copy(CC, C);

        // Check | A - L U | / | A |
        for (ordinal_type t=0;t<_ntridiag;++t) {
          a_subtract_mult_l_and_u(t, 0, AA, 
                                  t, 0, DD, DD);
          for (ordinal_type i=1;i<_m;++i) {
            a_subtract_mult_l_and_u(t, i,   AA, 
                                    t, i,   DD, DD);
            a_subtract_mult_b_and_c(t, i,   AA, 
                                    t, i-1, LL, UU);
            a_subtract_mult_l_and_b(t, i-1, BB, 
                                    t, i-1, DD, UU);
            a_subtract_mult_b_and_u(t, i-1, CC, 
                                    t, i-1, LL, DD);
          }
        }

        double norm = 0, diff = 0;
        for (ordinal_type t=0;t<_ntridiag;++t) {
          for (ordinal_type ii=0;ii<_blocksize;++ii)
            for (ordinal_type jj=0;jj<_blocksize;++jj) {
              norm += std::abs(tdiag_val(A ,t, 0, ii, jj));
              diff += std::abs(tdiag_val(AA,t, 0, ii, jj));
            }
          for (ordinal_type i=1;i<_m;++i) 
            for (ordinal_type ii=0;ii<_blocksize;++ii)
              for (ordinal_type jj=0;jj<_blocksize;++jj) {
                norm += std::abs(tdiag_val(A ,t, i,   ii, jj));
                diff += std::abs(tdiag_val(AA,t, i,   ii, jj));
                norm += std::abs(tdiag_val(B ,t, i-1, ii, jj));
                diff += std::abs(tdiag_val(BB,t, i-1, ii, jj));
                norm += std::abs(tdiag_val(C ,t, i-1, ii, jj));
                diff += std::abs(tdiag_val(CC,t, i-1, ii, jj));
              }
        }
        //std::cout << "tridiag factor check  norm = " << norm << "  diff = " << diff << std::endl;
        const bool r_val = diff/norm < 1e2*std::numeric_limits<scalar_type>::epsilon();
        return r_val;
      }
    };

    template <typename ExecSpace, typename ValueType, typename ArrayLayout,
              typename Trsv_AlgoTagType,
              typename Gemv_AlgoTagType>
    class SolveBlockTridiagMatrices {
    public:
      typedef ExecSpace exec_space;
      typedef ValueType value_type;
      typedef ArrayLayout array_layout;

      typedef BlockTridiagMatrices<exec_space,value_type,array_layout> block_tridiag_matrices_type;
      typedef PartitionedBlockMultiVector<exec_space,value_type,array_layout> partitioned_block_multi_vector_type;
      
    private:
      ordinal_type _ntridiag, _m, _blocksize, _nvectors;

      ConstUnmanagedViewType<typename block_tridiag_matrices_type::value_array_type> _TA, _TB, _TC;
      ConstUnmanagedViewType<typename partitioned_block_multi_vector_type::value_array_type> _b;
      /**/ UnmanagedViewType<typename partitioned_block_multi_vector_type::value_array_type> _x;

    public:
      SolveBlockTridiagMatrices() {}

      KOKKOS_INLINE_FUNCTION 
      void operator()(const ordinal_type ij) const {
        auto A = Kokkos::subview(_TA, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto B = Kokkos::subview(_TB, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto C = Kokkos::subview(_TC, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

        ///
        /// loop over multivectors
        ///
        for (int jvec=0;jvec<_nvectors;++jvec) {
          auto x = Kokkos::subview(_x, ij, jvec, Kokkos::ALL(), Kokkos::ALL());
          auto b = Kokkos::subview(_b, ij, jvec, Kokkos::ALL(), Kokkos::ALL());

          ///
          /// forward substitution
          ///
          {
            const bool is_same_x_and_b = (x.data() == b.data());
            {
              auto x0 = Kokkos::subview(x, 0, Kokkos::ALL());          
              auto b0 = Kokkos::subview(b, 0, Kokkos::ALL());
              if (!is_same_x_and_b)
                for (ordinal_type ii=0;ii<_blocksize;++ii)
                  x0(ii) = b0(ii);
            }
            const ordinal_type kend = _m - 1;
            for (ordinal_type k=0;k<kend;++k) {
              auto LT = Kokkos::subview(A, k,   Kokkos::ALL(), Kokkos::ALL());
              auto LB = Kokkos::subview(C, k,   Kokkos::ALL(), Kokkos::ALL());

              auto xt = Kokkos::subview(x, k,   Kokkos::ALL());
              auto xb = Kokkos::subview(x, k+1, Kokkos::ALL());

              auto bb = Kokkos::subview(b, k+1, Kokkos::ALL());

              if (!is_same_x_and_b)
                for (ordinal_type ii=0;ii<_blocksize;++ii)
                  xb(ii) = bb(ii);

              Serial::Trsv<Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsv_AlgoTagType>
                ::invoke(1.0, LT, xt);
              Serial::Gemv<Trans::NoTranspose,Gemv_AlgoTagType>
                ::invoke(-1.0, LB, xt, 1.0, xb);
            }
            auto LL = Kokkos::subview(A, kend, Kokkos::ALL(), Kokkos::ALL());
            auto xx = Kokkos::subview(x, kend, Kokkos::ALL());
            Serial::Trsv<Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsv_AlgoTagType>
              ::invoke(1.0, LL, xx);
          }

          ///
          /// backward substitution
          ///
          {
            const ordinal_type kbegin = _m - 1;
            for (ordinal_type k=kbegin;k>0;--k) {
              auto UT = Kokkos::subview(B, k-1, Kokkos::ALL(), Kokkos::ALL());
              auto UB = Kokkos::subview(A, k,   Kokkos::ALL(), Kokkos::ALL());
            
              auto xt = Kokkos::subview(x, k-1, Kokkos::ALL());
              auto xb = Kokkos::subview(x, k,   Kokkos::ALL());

              Serial::Trsv<Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsv_AlgoTagType>
                ::invoke(1.0, UB, xb);
              Serial::Gemv<Trans::NoTranspose,Gemv_AlgoTagType>
                ::invoke(-1.0, UT, xb, 1.0, xt);
            }
            auto UU = Kokkos::subview(A, 0, Kokkos::ALL(), Kokkos::ALL());
            auto xx = Kokkos::subview(x, 0, Kokkos::ALL());
            Serial::Trsv<Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsv_AlgoTagType>
              ::invoke(1.0, UU, xx);
          }
        }
      }

      void run(const block_tridiag_matrices_type T,
               const partitioned_block_multi_vector_type x,
               const partitioned_block_multi_vector_type b) {
        assert(T.NumTridiagMatrices() == x.NumPartitions());
        assert(T.NumRows() == x.NumRows());
        assert(T.BlockSize() == x.BlockSize());

        _ntridiag = T.NumTridiagMatrices();
        _m = T.NumRows(); 
        _blocksize = T.BlockSize();
        _nvectors = x.NumVectors();
        
        _TA = T.A(); 
        _TB = T.B(); 
        _TC = T.C();
        
        _x = x.Values();
        _b = b.Values();

        // parallel over the instances of tridiagonal matrices
        Kokkos::parallel_for(_ntridiag, *this);
      }

      template<typename RViewType, typename AViewType, typename XViewType>
      void r_subtract_mult_a_and_x(const ordinal_type tr, const ordinal_type ir, RViewType R,
                                   const ordinal_type ta, const ordinal_type ia, AViewType A,
                                   const ordinal_type tx, const ordinal_type ix, XViewType X) {
        for (ordinal_type kk=0;kk<_nvectors;++kk) 
          for (ordinal_type ii=0;ii<_blocksize;++ii)
            for (ordinal_type jj=0;jj<_blocksize;++jj)
              tdiag_val(R, tr, kk, ir, ii) -= tdiag_val(A, ta, ia, ii, jj) * tdiag_val(X, tx, kk, ix, jj);
      }

      bool check(const block_tridiag_matrices_type T, 
                 const partitioned_block_multi_vector_type b) {
        // input A
        auto AA = Kokkos::create_mirror_view(T.A());      Kokkos::deep_copy(AA, T.A());
        auto BB = Kokkos::create_mirror_view(T.B());      Kokkos::deep_copy(BB, T.B());
        auto CC = Kokkos::create_mirror_view(T.C());      Kokkos::deep_copy(CC, T.C());

        auto bb = Kokkos::create_mirror_view(b.Values()); Kokkos::deep_copy(bb, b.Values());
        auto xx = Kokkos::create_mirror_view(_x);         Kokkos::deep_copy(xx, _x);
        
        // diffs
        Kokkos::View<value_type****,Kokkos::DefaultHostExecutionSpace> 
          rr("rr", bb.dimension_0(), bb.dimension_1(), bb.dimension_2(), bb.dimension_3()); 
        
        Kokkos::deep_copy(rr, bb);

        // Check | Ax - b | / | b |
        for (ordinal_type t=0;t<_ntridiag;++t) {
          r_subtract_mult_a_and_x(t, 0, rr, 
                                  t, 0, AA, 
                                  t, 0, xx);
          r_subtract_mult_a_and_x(t, 0, rr, 
                                  t, 0, BB, 
                                  t, 1, xx);
          
          for (ordinal_type i=1;i<(_m-1);++i) {
            r_subtract_mult_a_and_x(t, i,   rr, 
                                    t, i-1, CC, 
                                    t, i-1, xx);
            r_subtract_mult_a_and_x(t, i,   rr, 
                                    t, i,   AA, 
                                    t, i,   xx);
            r_subtract_mult_a_and_x(t, i,   rr, 
                                    t, i,   BB, 
                                    t, i+1, xx);
          }
          r_subtract_mult_a_and_x(t, _m-1, rr, 
                                  t, _m-2, CC, 
                                  t, _m-2, xx);
          r_subtract_mult_a_and_x(t, _m-1, rr, 
                                  t, _m-1, AA, 
                                  t, _m-1, xx);
        }

        double norm = 0, diff = 0;
        for (ordinal_type t=0;t<_ntridiag;++t) 
          for (ordinal_type jvec=0;jvec<_nvectors;++jvec) 
            for (ordinal_type i=0;i<_m;++i) 
              for (ordinal_type ii=0;ii<_blocksize;++ii) {
                norm += std::abs(tdiag_val(bb, t, jvec, i, ii));
                diff += std::abs(tdiag_val(rr, t, jvec, i, ii));
              }

        //std::cout << "tridiag solve check  norm = " << norm << "  diff = " << diff << std::endl;
        const bool r_val = diff/norm < 1e2*std::numeric_limits<scalar_type>::epsilon();
        return r_val;
      }
    };


    // unit tests
    template<typename DeviceSpace, typename ValueType = scalar_type>
    void run(const ordinal_type ni, const ordinal_type nj, const ordinal_type nk, 
             const ordinal_type blocksize, 
             const ordinal_type nrhs,
             const bool test_cublas = false) {
      typedef typename DeviceSpace::array_layout DeviceArrayLayout;
      typedef Kokkos::DefaultHostExecutionSpace HostSpace;

      bool success = true;
      StructuredBlock mesh(ni, nj, nk);

      // Test StructuredBlock.
      for (ordinal_type c=0;c<mesh.size();++c) {
        ordinal_type i, j, k;
        mesh.id2ijk(c, i, j, k);
        TEST_ASSERT(i >= 0 && i < mesh.ni, success);
        TEST_ASSERT(j >= 0 && j < mesh.nj, success);
        TEST_ASSERT(k >= 0 && k < mesh.nk, success);
        TEST_ASSERT(mesh.ijk2id(i, j, k) == c, success);
      }

      // Graph construction
      CrsGraph<HostSpace,DeviceArrayLayout> graph_host 
        = create_graph_host_for_structured_block<DeviceArrayLayout>(mesh, StencilShape::cross);

      // Crs matrix and multi vector construction
      BlockCrsMatrix<HostSpace,DeviceArrayLayout> A_host(graph_host, blocksize);
      fill_block_crs_matrix_host(A_host);      

      // Device mirroring
      auto A_device = create_mirror<DeviceSpace>(A_host);
      deep_copy(A_device, A_host);

      // Test Matrix Vector product
      {
        const ordinal_type m = graph_host.NumRows();
        
        BlockMultiVector<HostSpace,DeviceArrayLayout> x_host(nrhs, m, blocksize);
        fill_block_multi_vector_host(x_host);
        
        auto x_device = create_mirror<DeviceSpace>(x_host);
        deep_copy(x_device, x_host);
        
        BlockMultiVector<DeviceSpace,DeviceArrayLayout> 
          y1_device(nrhs, m, blocksize), 
          y2_device(nrhs, m, blocksize);

        {
          BlockCrsMatrixVectorProductByRow<DeviceSpace,DeviceArrayLayout> matvec;
          matvec.run(A_device, x_device, y1_device);
        }
        {
          BlockCrsMatrixVectorProductByBlockRow<DeviceSpace,DeviceArrayLayout> matvec;
          matvec.run(A_device, x_device, y2_device);
        }

        const double rdiff = compute_relative_diff(y1_device.Values(), y2_device.Values());
        TEST_ASSERT(rdiff <= 1e2*std::numeric_limits<scalar_type>::epsilon(), success);
      }

      // Test Block TriDiag Extraction
      BlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> T_device
        = create_block_tridiag_matrices
        <DeviceSpace,ValueType,DeviceArrayLayout>(mesh.ni*mesh.nj, 
                                                  mesh.nk,
                                                  blocksize);
      {
        ExtractBlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> extblk(mesh);
        extblk.run(A_device, T_device);
        TEST_ASSERT(extblk.check(), success);
      }

      BlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> T_org_device
        = create_block_tridiag_matrices
        <DeviceSpace,ValueType,DeviceArrayLayout>(mesh.ni*mesh.nj, 
                                                  mesh.nk,
                                                  blocksize);
      
      deep_copy(T_org_device, T_device);
      
      // Test Block TriDiag Factorization
      if (test_cublas) {
        //
      } else {
        FactorizeBlockTridiagMatrices<DeviceSpace,
                                      ValueType,
                                      DeviceArrayLayout,
                                      Algo::LU::Blocked,
                                      Algo::Trsm::Blocked,
                                      Algo::Gemm::Blocked> factorblk;
        factorblk.run(T_device);
        TEST_ASSERT(factorblk.check(T_org_device), success);
      }

      // Test Block TriDiag Solve
      {
        PartitionedBlockMultiVector<HostSpace,ValueType,DeviceArrayLayout> b_host
          = create_partitioned_block_multi_vector
          <HostSpace,ValueType,DeviceArrayLayout>(mesh.ni*mesh.nj, 
                                                  nrhs,
                                                  mesh.nk, 
                                                  blocksize);
        fill_partitioned_block_multi_vector_host(b_host, mesh.ni*mesh.nj);

        auto b_device = create_mirror<DeviceSpace>(b_host);
        deep_copy(b_device, b_host);

        PartitionedBlockMultiVector<DeviceSpace,ValueType,DeviceArrayLayout> x_device
          = create_partitioned_block_multi_vector
          <DeviceSpace,ValueType,DeviceArrayLayout>(mesh.ni*mesh.nj, 
                                                    nrhs,
                                                    mesh.nk, 
                                                    blocksize);
        if (test_cublas) {
          //
        } else {
          SolveBlockTridiagMatrices<DeviceSpace,
                                    ValueType,
                                    DeviceArrayLayout,
                                    Algo::Trsv::Blocked,
                                    Algo::Gemv::Blocked> solveblk;
          
          solveblk.run(T_device, x_device, b_device);
          TEST_ASSERT(solveblk.check(T_org_device, b_device), success);
        }
      }

      if (!success)
        std::cout << "Unit Tests:: Failed:: "
                  << " ni = " << ni << " nj = " << nj << " nk = " << nk 
                  << " blocksize = " << blocksize << " nrhs = " << nrhs << " \n"; 
    }
    
    // performance tests
    template<typename DeviceSpace, typename ValueType = scalar_type>
    int run(const Input &input, const bool test_cublas = false) { 
      typedef typename DeviceSpace::array_layout DeviceArrayLayout;
      typedef Kokkos::DefaultHostExecutionSpace HostSpace;

      const ordinal_type niter = 50;      
      int dontopt = 0;
      bool success = true;
      
      /// 
      /// construct a discrete system of equations
      ///
      const ordinal_type 
        ni = input.ni, 
        nj = input.nj, 
        nk = input.nk,
        blocksize = input.bs, 
        nrhs = input.nrhs;

      StructuredBlock mesh(ni, nj, nk);

      // something is not copyable ... don't know why yet...
      BlockCrsMatrix<DeviceSpace,DeviceArrayLayout> A_device;
      double t_fill_block_crs_matrix = 0.0, t_fill_graph = 0.0;
      {
        const StencilShape::Enum stencil_shape = input.stencil_shape;
        CrsGraph<HostSpace,DeviceArrayLayout> graph_host;
        {
          Timer timer("Fill Graph _______________");
          timer.reset();
          graph_host = create_graph_host_for_structured_block<DeviceArrayLayout>(mesh, stencil_shape);
          t_fill_graph = timer.seconds();
        }
        BlockCrsMatrix<HostSpace,DeviceArrayLayout> A_host(graph_host, blocksize);
        {
          Timer timer("Fill Block CRS Matrix_______________");
          timer.reset();
          fill_block_crs_matrix_host(A_host);       
          t_fill_block_crs_matrix = timer.seconds();       
        }
        A_device = create_mirror<DeviceSpace>(A_host);
        deep_copy(A_device, A_host);
      }

      // memory size
      const double memsize_A = A_device.Values().dimension_0()*blocksize*blocksize*8;
      
      ///
      /// matrix vector multiplication test
      ///
      double t_matvec = 0.0;
      double t_fill_block_multi_vector = 0.0;
      {
        const ordinal_type m = mesh.size();

        BlockMultiVector<HostSpace,DeviceArrayLayout> x_host(nrhs, m, blocksize);
        {
          Timer timer("Fill Block Multi Vector______________");
          timer.reset();
          fill_block_multi_vector_host(x_host);
          t_fill_block_multi_vector = timer.seconds();
        }
        auto x_device = create_mirror<DeviceSpace>(x_host);
        deep_copy(x_device, x_host);
        
        BlockMultiVector<DeviceSpace,DeviceArrayLayout> y_device(nrhs, m, blocksize);
        {
          //BlockCrsMatrixVectorProductByRow<DeviceSpace> matvec;
          BlockCrsMatrixVectorProductByBlockRow<DeviceSpace,DeviceArrayLayout> matvec;
          {
            Timer timer("50 BlockCrsMatrixVectorProduct");
            timer.reset();
            for (ordinal_type i=0;i<niter;++i) {
              matvec.run(A_device, x_device, y_device);
              dontopt += i;
            }
            t_matvec = timer.seconds();
          }
        }
      }

      ///
      /// block tridiag extraction test
      ///
      const double memsize_T = ni*nj*(3*(nk-1)*blocksize*blocksize + blocksize*blocksize)*8;

      double t_extract = 0.0;
      BlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> T_device
        = create_block_tridiag_matrices<DeviceSpace,ValueType,DeviceArrayLayout>(ni*nj, nk, blocksize);
      {
        ExtractBlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> extblk(mesh);
        {
          Timer timer("ExtractBlockTridiagMatrices");
          timer.reset();
          extblk.run(A_device, T_device);
          t_extract = timer.seconds();
        }
        if (input.check) TEST_ASSERT(extblk.check(), success);
      }

      // keep original matrix for check
      BlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> T_org_device
        = create_block_tridiag_matrices<DeviceSpace,ValueType,DeviceArrayLayout>(ni*nj, nk, blocksize);
      
      deep_copy(T_org_device, T_device);
      
      ///
      /// block tridiag factorization test
      ///
      double t_factorize = 0.0, f_factorize = 0.0;
      if (test_cublas) {
        //
      } else {
        FactorizeBlockTridiagMatrices<DeviceSpace,
                                      ValueType,
                                      DeviceArrayLayout,
                                      Algo::LU::Blocked,
                                      Algo::Trsm::Blocked,
                                      Algo::Gemm::Blocked> factorblk;
        
        f_factorize = factorblk.FlopCount(T_device)*(sizeof(ValueType)/sizeof(double));
        {
          Timer timer("FactorizeBlockTridiagMatrices");
          timer.reset();
          factorblk.run(T_device);
          t_factorize = timer.seconds();
        }
        if (input.check) TEST_ASSERT(factorblk.check(T_org_device), success);
      }
      
      ///
      /// block tridiag solve test
      ///
      double t_solve = 0.0;
      {
        PartitionedBlockMultiVector<HostSpace,ValueType,DeviceArrayLayout> b_host
          = create_partitioned_block_multi_vector
          <HostSpace,ValueType,DeviceArrayLayout>(ni*nj, 
                                                  nrhs,
                                                  nk, 
                                                  blocksize);
        fill_partitioned_block_multi_vector_host(b_host, ni*nj);

        auto b_device = create_mirror<DeviceSpace>(b_host);
        deep_copy(b_device, b_host);
        
        PartitionedBlockMultiVector<DeviceSpace,ValueType,DeviceArrayLayout> x_device
          = create_partitioned_block_multi_vector
          <DeviceSpace,ValueType,DeviceArrayLayout>(ni*nj, 
                                                    nrhs,
                                                    nk, 
                                                    blocksize);
        if (test_cublas) {
          //
        } else {
          SolveBlockTridiagMatrices<DeviceSpace,
                                    ValueType,
                                    DeviceArrayLayout,
                                    Algo::Trsv::Blocked,
                                    Algo::Gemv::Blocked> solveblk;
          {
            Timer timer("50 SolveBlockTridiagMatrices");
            timer.reset();
            for (ordinal_type i=0;i<niter;++i) {
              solveblk.run(T_device, x_device, b_device);
              dontopt += i;
            }          
            t_solve = timer.seconds();
          }
          if (input.check) TEST_ASSERT(solveblk.check(T_org_device, b_device), success);
        }
      }
      
      const double t_matvec_per_iter = t_matvec/double(niter), t_solve_per_iter = t_solve/double(niter);
      std::cout << " matvec     = " << t_matvec_per_iter << std::endl; 
      std::cout << " extract    = " << t_extract        << " extract/matvec = " << (t_extract/t_matvec_per_iter) << std::endl; 
      //std::cout << " factor     = " << t_factorize      << " factor/matvec  = " << (t_factorize/t_matvec_per_iter) << std::endl; 
      std::cout << " factor     = " << t_factorize      << " factor/matvec  = " << (t_factorize/t_matvec_per_iter) << " flop = " << f_factorize << " flop/s = " << (f_factorize/t_factorize) << std::endl; 
      std::cout << " solve      = " << t_solve_per_iter << "  solve/matvec  = " << (t_solve_per_iter/t_matvec_per_iter) << std::endl; 
      std::cout << " memory used     = " << (memsize_A + memsize_T) << std::endl; 

      return dontopt + success;
    }

  }
}
