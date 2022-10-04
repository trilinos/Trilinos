// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_LOCALCRSMATRIXOPERATOR_DECL_HPP
#define TPETRA_LOCALCRSMATRIXOPERATOR_DECL_HPP

#include "Tpetra_LocalCrsMatrixOperator_fwd.hpp"
#include "Tpetra_LocalOperator.hpp"
#include "Tpetra_Spaces.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"

#include <memory> // std::shared_ptr

namespace Tpetra {

  /// \class LocalCrsMatrixOperator
  /// \brief Abstract interface for local operators (e.g., matrices
  ///   and preconditioners).
  ///
  /// \tparam MultiVectorScalar The type of the entries of the input
  ///   and output (multi)vectors.
  /// \tparam MatrixScalar The type of the entries of the sparse matrix.
  /// \tparam Device The Kokkos Device type; must be a specialization
  ///   of Kokkos::Device.
  template<class MultiVectorScalar, class MatrixScalar, class Device>
  class LocalCrsMatrixOperator :
    public LocalOperator<MultiVectorScalar, Device> {
  private:
    using mv_scalar_type =
      typename LocalOperator<MultiVectorScalar, Device>::scalar_type;
    using matrix_scalar_type =
      typename LocalOperator<MatrixScalar, Device>::scalar_type;
    using array_layout =
      typename LocalOperator<MultiVectorScalar, Device>::array_layout;
    using device_type =
      typename LocalOperator<MultiVectorScalar, Device>::device_type;
    using local_ordinal_type =
      ::Tpetra::Details::DefaultTypes::local_ordinal_type;
    using execution_space = typename Device::execution_space;
  public:
    using local_matrix_device_type =
      KokkosSparse::CrsMatrix<matrix_scalar_type,
                              local_ordinal_type,
                              device_type,
                              void,
                              size_t>;
  private:
    //The type of a matrix with offset=ordinal, but otherwise the same as local_matrix_device_type
    using local_cusparse_matrix_type =
      KokkosSparse::CrsMatrix<matrix_scalar_type,
                              local_ordinal_type,
                              device_type,
                              void,
                              local_ordinal_type>;
    using local_graph_device_type = typename local_matrix_device_type::StaticCrsGraphType;

  public:

    // construct a view of the on-rank entries in a row
    template<typename OffsetDeviceViewType>
    struct OnRankRowViewer {
      KOKKOS_INLINE_FUNCTION static KokkosSparse::SparseRowViewConst<local_matrix_device_type> view(
          const local_matrix_device_type &A,
          const OffsetDeviceViewType &aOff,
          const local_ordinal_type i
      ) {
        const typename local_matrix_device_type::size_type start = A.graph.row_map(i);
        const local_ordinal_type stride = 1;
        const local_ordinal_type length = aOff[i] - start;
        return KokkosSparse::SparseRowViewConst<local_matrix_device_type> (
          &A.values(start),
          &A.graph.entries(start),
          stride,
          length
        );
      }
    };

    // construct a view of the off-rank entries in a row
    template<typename OffsetDeviceViewType>
    struct OffRankRowViewer {
      KOKKOS_INLINE_FUNCTION static KokkosSparse::SparseRowViewConst<local_matrix_device_type> view(
          const local_matrix_device_type &A,
          const OffsetDeviceViewType &aOff,
          const local_ordinal_type i
      ) {
        const typename local_matrix_device_type::size_type start = aOff[i];
        const local_ordinal_type stride = 1;
        const local_ordinal_type length = A.graph.row_map(i+1) - start;
        return KokkosSparse::SparseRowViewConst<local_matrix_device_type> (
          &A.values(start),
          &A.graph.entries(start),
          stride,
          length
        );
      }
    };

#if 0
    // cwp 21 Sep 2022
    // A functor that does the on-rank part of a local SpMV
    // KokkosKernels does not currently have a 4-array CSR
    template<
    typename OffsetDeviceViewType, 
    typename RowViewer // OnRankRowView or OffRankRowView
    > class SpmvFunctor {

    public:

      // different SpMV execution modes
      struct TagNonTrans{};
      struct TagTrans{};
      struct TagConjTrans{};

      typedef LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device> local_crs_matrix_operator_type;
      typedef typename local_crs_matrix_operator_type::local_matrix_device_type local_matrix_device_type;
      typedef typename local_crs_matrix_operator_type::array_layout array_layout;
      
      typedef Kokkos::View<const MultiVectorScalar**, array_layout,
            Device, Kokkos::MemoryTraits<Kokkos::Unmanaged> > x_type;
      typedef Kokkos::View<MultiVectorScalar**, array_layout,
            Device, Kokkos::MemoryTraits<Kokkos::Unmanaged> > y_type;
    
    private:
      MultiVectorScalar alpha_;
      local_matrix_device_type A_;
      x_type X_;
      MultiVectorScalar beta_;
      y_type Y_;
      OffsetDeviceViewType offRankOffsets_;

      typedef typename local_matrix_device_type::non_const_value_type value_type;
      typedef typename local_matrix_device_type::non_const_ordinal_type ordinal_type; 
      typedef typename local_matrix_device_type::non_const_size_type size_type; 

      typedef typename local_matrix_device_type::execution_space execution_space;
      typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
      typedef typename team_policy::member_type team_member;
      typedef Kokkos::Details::ArithTraits<value_type> ATV;

    public:
      SpmvFunctor(const MultiVectorScalar &alpha, 
      const local_matrix_device_type &A, 
      x_type &X, 
      const MultiVectorScalar &beta, 
      y_type &Y,
      const OffsetDeviceViewType &offRankOffsets) 
        : alpha_(alpha), A_(A), X_(X), beta_(beta), Y_(Y), offRankOffsets_(offRankOffsets) {}

      /*! \brief contribution of row i

          This is only the local offsets, so the offsets we want are between
          rowPtr(i) and offRankOffsets_(i)
      */
      KOKKOS_INLINE_FUNCTION void operator()(TagNonTrans, const ordinal_type i) const {

        if (i >= A_.numRows()) {
          return;
        }

        const KokkosSparse::SparseRowViewConst<local_matrix_device_type> row = 
            RowViewer::view(A_, offRankOffsets_, i);

        // this implementation may be best for a single vector (not multivector)
        for (ordinal_type k = 0; k < Y_.extent(1); ++k) {
          MultiVectorScalar sum = 0;

          for (ordinal_type ri = 0; ri < row.length; ++ri) {
            value_type A_ij = row.value(ri);
            ordinal_type j = row.colidx(ri);
            sum += A_ij * X_(j, k); 
          }
          sum *= alpha_;

          if (0 == beta_) {
            Y_(i,k) = sum;
          } else {
            Y_(i,k) = beta_ * Y_(i,k) + sum;
          } 
        }
      }

      KOKKOS_INLINE_FUNCTION
      void operator() (TagNonTrans, const team_member& dev) const
      {
        using y_value_type = typename y_type::non_const_value_type;

        const int rowsPerTeam = dev.team_size();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(dev,0,rowsPerTeam), [&] (const ordinal_type& loop) {

          const ordinal_type iRow = static_cast<ordinal_type> ( dev.league_rank() ) * rowsPerTeam + loop;
          if (iRow >= A_.numRows ()) {
            return;
          }
          const KokkosSparse::SparseRowViewConst<local_matrix_device_type> row = 
            RowViewer::view(A_, offRankOffsets_, iRow);
          const ordinal_type row_length = static_cast<ordinal_type> (row.length);
          y_value_type sum = 0;

          Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(dev,row_length), [&] (const ordinal_type& iEntry, y_value_type& lsum) {
            const value_type val = row.value(iEntry);
            lsum += val * X_(row.colidx(iEntry), 0);
          },sum);

          Kokkos::single(Kokkos::PerThread(dev), [&] () {
            sum *= alpha_;

            #warning "assuming rank-1 vector"
            Y_(iRow, 0) = beta_ * Y_(iRow, 0) + sum;
          });
        });
      }


      KOKKOS_INLINE_FUNCTION void operator()(TagTrans, const size_t i) const {
        #warning trans unimplemented
      }

      KOKKOS_INLINE_FUNCTION void operator()(TagConjTrans, const size_t i) const {
        #warning conj trans unimplemented
      }

      /// \brief Kokkos dispatch of non-transpose
      // void launch(TagNonTrans, const Kokkos::DefaultExecutionSpace &space) {
      void launch(TagNonTrans, const execution_space &space) {

#if 1 // fancy one
        const int estNnz = A_.nnz() * 0.95; // entries are mostly on-rank.
        if (!KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()) {
          if (estNnz > 10000000) {
            Kokkos::parallel_for(Kokkos::RangePolicy<TagNonTrans, execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(space, 0, A_.numRows()), *this);
          } else {
            Kokkos::parallel_for(Kokkos::RangePolicy<TagNonTrans, execution_space, Kokkos::Schedule<Kokkos::Static>>(space, 0, A_.numRows()), *this);
          }
        } else {
          const int nnzPerRow = (estNnz + A_.numRows() - 1) / A_.numRows();
          int vectorLength = 1;
          while (vectorLength < 32 && vectorLength * 6 < nnzPerRow) {
            vectorLength *= 2;
          }
          const int teamSize = 256 / vectorLength;
          const int rowsPerTeam = teamSize;
          const int64_t worksets = (Y_.extent(0)+rowsPerTeam-1)/rowsPerTeam;

          if (estNnz > 10000000) {
            Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>, TagNonTrans>policy(space, worksets, teamSize, vectorLength);
            Kokkos::parallel_for("on-rank", policy, *this);
          } else {
            Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>, TagNonTrans>policy(space, worksets,teamSize,vectorLength);
            Kokkos::parallel_for("on-rank", policy, *this);
          }
        }
#else // boring one
        Kokkos::RangePolicy<execution_space, TagNonTrans> policy(space, 0, A_.numRows());
        Kokkos::parallel_for("on-rank", policy, *this);
#endif
      }

      // \brief Kokkos dispatch of non-transpose in default space
      void launch(TagNonTrans t) {
        launch(t, execution_space());
      }

      /// \brief Kokkos dispatch of transpose
      void launch(TagTrans, execution_space &space) {
        Kokkos::parallel_for(Kokkos::RangePolicy<TagTrans, execution_space>(space, 0, A_.numRows()), *this);
      }

      /// \brief Kokkos dispatch of conjugate transpose
      void launch(TagConjTrans, execution_space &space) {
        Kokkos::parallel_for(Kokkos::RangePolicy<TagConjTrans, execution_space>(space, 0, A_.numRows()), *this);
      }
    };
#endif

    // cwp 21 Sep 2022
    // A functor that does on-rank or off-rank part of a local Sparse-matrix multivector product
    // KokkosKernels does not currently have a 4-array CSR
    // SPMV_MV_LayoutLeft_Functor
    template<
    typename OffsetDeviceViewType, 
    typename RowViewer, // OnRankRowView or OffRankRowView
    bool CONJ
    > class SpmvMvFunctor {
    public:
      typedef LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device> local_crs_matrix_operator_type;
      typedef typename local_crs_matrix_operator_type::local_matrix_device_type local_matrix_device_type;
      typedef typename local_crs_matrix_operator_type::array_layout array_layout;
      
      typedef Kokkos::View<const MultiVectorScalar**, array_layout,
            Device, Kokkos::MemoryTraits<Kokkos::Unmanaged> > x_type;
      typedef Kokkos::View<MultiVectorScalar**, array_layout,
            Device, Kokkos::MemoryTraits<Kokkos::Unmanaged> > y_type;
    
    private:
      typedef typename local_matrix_device_type::non_const_value_type value_type;
      typedef typename local_matrix_device_type::non_const_ordinal_type ordinal_type; 
      typedef typename local_matrix_device_type::non_const_size_type size_type; 

      typedef typename local_matrix_device_type::execution_space execution_space;
      typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
      typedef typename team_policy::member_type team_member;
      typedef Kokkos::Details::ArithTraits<value_type> ATV;

      MultiVectorScalar alpha_;
      local_matrix_device_type A_;
      x_type X_;
      MultiVectorScalar beta_;
      y_type Y_;
      OffsetDeviceViewType offRankOffsets_;
      ordinal_type rowsPerThread_;
      int vectorLength_;

    public:
      SpmvMvFunctor(const MultiVectorScalar &alpha, 
      const local_matrix_device_type &A, 
      x_type &X, 
      const MultiVectorScalar &beta, 
      y_type &Y,
      const OffsetDeviceViewType &offRankOffsets,
      ordinal_type rowsPerThread,
      int vectorLength) 
        : alpha_(alpha), A_(A), X_(X), beta_(beta), Y_(Y),
         offRankOffsets_(offRankOffsets),
         rowsPerThread_(rowsPerThread), vectorLength_(vectorLength) {}

      
      /* simplified version brought from Kokkos Kernels
      */
      template <int UNROLL>
      KOKKOS_INLINE_FUNCTION void strip_mine(const team_member& dev,
                                             const ordinal_type& iRow,
                                             const ordinal_type& kk) const {
        MultiVectorScalar sum[UNROLL];

        for (int k = 0; k < UNROLL; ++k) {
          sum[k] = Kokkos::Details::ArithTraits<MultiVectorScalar>::zero();
        }

        const auto row = RowViewer::view(A_, offRankOffsets_, iRow);

        // split loop indices 0..row.length over the vector lanes of the calling thread
        // each lane will have a subset of the partial products
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(dev, row.length), [&](ordinal_type iEntry) {
              const value_type val =
                  CONJ ? Kokkos::Details::ArithTraits<value_type>::conj(
                                  row.value(iEntry))
                            : row.value(iEntry);
              const ordinal_type ind = row.colidx(iEntry);
              for (int k = 0; k < UNROLL; ++k) {
                sum[k] += val * X_(ind, kk + k);
              }
            }
        );

#if 0
        for (int ii = 0; ii < UNROLL; ++ii) {
          int tr = dev.team_rank();
          int lr = dev.league_rank();
          printf("pre:         <%d,%d> sum(%d, %d) = %d\n",
              lr, tr, int(iRow), int(ii), int(sum[ii]));
        }
#endif

        for (int ii = 0; ii < UNROLL; ++ii) {

            // in CUDA, each lane is actually a sub-warp (thread) so each thread's sum is not correct
            // in serial, each lane is a true vector lane so sum[k] is already correct
            if (Spaces::is_gpu_exec_space<execution_space>()) {
              // try to sum up the sum[ii] in each vector lane
              MultiVectorScalar sumt;
              Kokkos::parallel_reduce(
                  Kokkos::ThreadVectorRange(dev, vectorLength_),
                  [&](ordinal_type, MultiVectorScalar& lsum) {
                    lsum += sum[ii];
                  },
                  sumt);
              // sumt is broadcast to every lane of the thread

              // now every lane should have the complete product
              // of the row this thread was working on
              sum[ii] = sumt * alpha_;
            } else {
              sum[ii] *= alpha_;
            }
        }

#if 0
        for (int ii = 0; ii < UNROLL; ++ii) {
          int tr = dev.team_rank();
          int lr = dev.league_rank();
          printf("all-reduced: <%d,%d> sum(%d, %d) = %d\n",
              lr, tr, int(iRow), int(ii), int(sum[ii]));
        }
#endif

        // split 0..UNROLL over the vector lanes of the calling thread
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(dev, UNROLL), [&](ordinal_type k) {
              Y_(iRow, kk + k) = beta_ * Y_(iRow, kk + k) + sum[k];
            }
        );
        
      }

      /* simplified version brought from Kokkos Kernels
      */
      KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {
        for (ordinal_type loop = 0; loop < rowsPerThread_; ++loop) {
          const ordinal_type iRow =
              (dev.league_rank() * dev.team_size() + dev.team_rank()) *
                  rowsPerThread_ + loop;
          if (iRow >= A_.numRows()) {
            return;
          }

          // mfh 20 Mar 2015, 07 Jun 2016: This is ordinal_type because it
          // needs to have the same type as n.
          ordinal_type kk = 0;

          const ordinal_type n = X_.extent(1);

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
          if ((n > 8) && (n % 8 == 1)) {
            strip_mine<9>(dev, iRow, kk);
            kk += 9;
          }
          for (; kk + 8 <= n; kk += 8) strip_mine<8>(dev, iRow, kk);
          if (kk < n) {
            switch (n - kk) {
#else   // NOT a GPU
          if ((n > 16) && (n % 16 == 1)) {
            strip_mine<17>(dev, iRow, kk);
            kk += 17;
          }

          for (; kk + 16 <= n; kk += 16) {
            strip_mine<16>(dev, iRow, kk);
          }
          if (kk < n) {
            switch (n - kk) {
              case 15: strip_mine<15>(dev, iRow, kk); break;

              case 14: strip_mine<14>(dev, iRow, kk); break;

              case 13: strip_mine<13>(dev, iRow, kk); break;

              case 12: strip_mine<12>(dev, iRow, kk); break;

              case 11: strip_mine<11>(dev, iRow, kk); break;

              case 10: strip_mine<10>(dev, iRow, kk); break;

              case 9: strip_mine<9>(dev, iRow, kk); break;

              case 8: strip_mine<8>(dev, iRow, kk); break;
#endif  // if/else: __CUDA_ARCH__ or __HIP_DEVICE_COMPILE__
              case 7: strip_mine<7>(dev, iRow, kk); break;

              case 6: strip_mine<6>(dev, iRow, kk); break;

              case 5: strip_mine<5>(dev, iRow, kk); break;

              case 4: strip_mine<4>(dev, iRow, kk); break;

              case 3: strip_mine<3>(dev, iRow, kk); break;

              case 2: strip_mine<2>(dev, iRow, kk); break;

              case 1: strip_mine<1>(dev, iRow, kk); break; // was strip_mine_1
            }
          }
        }
      }

      static void launch(const MultiVectorScalar &alpha, 
        const local_matrix_device_type &A, 
        x_type &X, 
        const MultiVectorScalar &beta, 
        y_type &Y,
        const OffsetDeviceViewType &offRankOffsets,
        const execution_space &space) {

        /* nothing to do */
        if (0 == A.numRows()) {
          return;
        }

        const ordinal_type NNZPerRow = A.nnz() / A.numRows();

        ordinal_type vectorLength = 1;
        while ((vectorLength * 2 * 3 <= NNZPerRow) && (vectorLength < 8)) {
          vectorLength *= 2;
        }
        // std::cerr << __FILE__<<":"<<__LINE__<<": vectorLength=" << vectorLength << "\n";

        const ordinal_type rowsPerThread =
            KokkosSparse::RowsPerThread<execution_space>(NNZPerRow);

        SpmvMvFunctor op(alpha, A, X, beta, Y, offRankOffsets,
                  rowsPerThread, vectorLength);

        const ordinal_type teamSize =
            Kokkos::TeamPolicy<execution_space>(space, 
                rowsPerThread, Kokkos::AUTO, vectorLength)
                .team_size_recommended(op, Kokkos::ParallelForTag());
        const ordinal_type rowsPerTeam = rowsPerThread * teamSize;
        const size_type nteams = (A.numRows() + rowsPerTeam - 1) / rowsPerTeam;
        Kokkos::parallel_for("SpmvMvFunctor",
                            Kokkos::TeamPolicy<execution_space>(
                                space, nteams, teamSize, vectorLength),
                            op);
      }

    };


    /* Simplified SPMV_MV_Transpose_Functor from Kokkos Kernels
    */
    template<
    typename OffsetDeviceViewType, 
    typename RowViewer, // OnRankRowView or OffRankRowView
    bool CONJ
    > struct SpmvMvTransFunctor {
    public:
      typedef LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device> local_crs_matrix_operator_type;
      typedef typename local_crs_matrix_operator_type::local_matrix_device_type local_matrix_device_type;
      typedef typename local_crs_matrix_operator_type::array_layout array_layout;
      
      typedef Kokkos::View<const MultiVectorScalar**, array_layout,
            Device, Kokkos::MemoryTraits<Kokkos::Unmanaged> > x_type;
      typedef Kokkos::View<MultiVectorScalar**, array_layout,
            Device, Kokkos::MemoryTraits<Kokkos::Unmanaged> > y_type;
    
    private:
      typedef typename local_matrix_device_type::non_const_value_type value_type;
      typedef typename local_matrix_device_type::non_const_ordinal_type ordinal_type; 
      typedef typename local_matrix_device_type::non_const_size_type size_type; 

      typedef typename local_matrix_device_type::execution_space execution_space;
      typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
      typedef typename team_policy::member_type team_member;
      typedef Kokkos::Details::ArithTraits<value_type> ATV;

      MultiVectorScalar alpha_;
      local_matrix_device_type A_;
      x_type X_;
      MultiVectorScalar beta_;
      y_type Y_;
      OffsetDeviceViewType offRankOffsets_;
    public:
      ordinal_type rowsPerTeam_ = 0;
    private:

      SpmvMvTransFunctor(const MultiVectorScalar& alpha, const local_matrix_device_type& A,
                                const x_type& X, const MultiVectorScalar& beta,
                                const y_type& Y, const OffsetDeviceViewType &offRankOffsets)
          : alpha_(alpha),
            A_(A),
            X_(X),
            beta_(beta),
            Y_(Y),
            offRankOffsets_(offRankOffsets) {}

    public:
      KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {

        const ordinal_type n = X_.extent(1);

        const ordinal_type teamWork = dev.league_rank() * rowsPerTeam_;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(dev, rowsPerTeam_), [&](ordinal_type loop) {
              // iRow represents a row of the matrix, so its correct type is
              // ordinal_type.
              const ordinal_type iRow = teamWork + loop;
              if (iRow >= A_.numRows()) {
                return;
              }

              const auto row = RowViewer::view(A_, offRankOffsets_, iRow);

              Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange(dev, row.length),
                  [&](ordinal_type iEntry) {
                    const value_type val =
                        CONJ
                            ? Kokkos::Details::ArithTraits<value_type>::conj(
                                  row.value(iEntry))
                            : row.value(iEntry);
                    const ordinal_type ind = row.colidx(iEntry);

    #ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
    #pragma unroll
    #endif
                    for (ordinal_type k = 0; k < n; ++k) {
                      Kokkos::atomic_add(
                          &Y_(ind, k),
                          static_cast<MultiVectorScalar>(alpha_ * val * X_(iRow, k)));
                    }
                  });
            });
      }

      // a kind of crummy scal implementation that we can put in an execution space instance
      struct ScalTag{};
      KOKKOS_INLINE_FUNCTION void operator()(const ScalTag&, const ordinal_type i) const {
        if (i >= A_.numRows()) {return;}
        for (ordinal_type k = 0; k < Y_.extent(1); ++k) {
          Y_(i, k) *= beta_;
        }
      }

      static void launch(const MultiVectorScalar &alpha, 
        const local_matrix_device_type &A, 
        x_type &X, 
        const MultiVectorScalar &beta, 
        y_type &Y,
        const OffsetDeviceViewType &offRankOffsets,
        const execution_space &space) {


        SpmvMvTransFunctor op(alpha, A, X, beta, Y, offRankOffsets);
        Kokkos::parallel_for("SpmvMvTransFunctor scal",
                            Kokkos::RangePolicy<ScalTag, execution_space>(
                                space, 0, A.numRows()),
                            op);

        const ordinal_type NNZPerRow = A.nnz() / A.numRows();

        ordinal_type vectorLength = 1;
        while ((vectorLength * 2 * 3 <= NNZPerRow) && (vectorLength < 8)) {
          vectorLength *= 2;
        }

        ordinal_type nrow = A.numRows();

        const ordinal_type rowsPerThread =
            KokkosSparse::RowsPerThread<execution_space>(NNZPerRow);
        const ordinal_type teamSize =
            Kokkos::TeamPolicy<execution_space>(space, 
                rowsPerThread, Kokkos::AUTO, vectorLength)
                .team_size_recommended(op, Kokkos::ParallelForTag());
        const ordinal_type rowsPerTeam = rowsPerThread * teamSize;
        op.rowsPerTeam_ = rowsPerTeam;
        const size_type nteams = (nrow + rowsPerTeam - 1) / rowsPerTeam;
        Kokkos::parallel_for("SpmvMvTransFunctor",
                            Kokkos::TeamPolicy<execution_space>(
                                space, nteams, teamSize, vectorLength),
                            op);
      }

    };

    
    using ordinal_view_type = typename local_graph_device_type::entries_type::non_const_type;

    LocalCrsMatrixOperator (const std::shared_ptr<local_matrix_device_type>& A);
    LocalCrsMatrixOperator (const std::shared_ptr<local_matrix_device_type>& A, const ordinal_view_type& A_ordinal_rowptrs);
    ~LocalCrsMatrixOperator () override = default;

    void
    apply (Kokkos::View<const mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > X,
           Kokkos::View<mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Y,
           const Teuchos::ETransp mode,
           const mv_scalar_type alpha,
           const mv_scalar_type beta) const override;

    void
    applyImbalancedRows (
           Kokkos::View<const mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > X,
           Kokkos::View<mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Y,
           const Teuchos::ETransp mode,
           const mv_scalar_type alpha,
           const mv_scalar_type beta) const;

    /*! \brief
        \c apply() but only contribute entries specified in offRankOffsets,
        which should be populated by \c CrsGraph::getLocalOffRankOffsets.
        Complement of \c applyLocalColumns

        \tparam OffsetDeviceViewType should be the CrsMatrix::crs_graph_type::offset_device_view_type of the \c Tpetra::CrsMatrix that owns this LocalCrsMatrixOperator
    
        cwp 05 Apr 2022
        applyRemoteColumns() with applyLocalColumns() shall have the same effect as apply()
    */ 
    template<typename OffsetDeviceViewType>
    void
    applyRemoteColumns (execution_space &execSpace,
           Kokkos::View<const mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > X,
           Kokkos::View<mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Y,
           const Teuchos::ETransp mode,
           const mv_scalar_type alpha,
           const mv_scalar_type beta,
           const OffsetDeviceViewType &offRankOffsets) const {

      typedef OffRankRowViewer<OffsetDeviceViewType> RowViewer;
      switch(mode) {
        case Teuchos::ETransp::TRANS: {
          typedef SpmvMvTransFunctor<OffsetDeviceViewType, RowViewer, false> Op;
          Op::launch(alpha, *A_, X, beta, Y, offRankOffsets, execSpace);
          return;
        }
        case Teuchos::ETransp::NO_TRANS: {
          typedef SpmvMvFunctor<OffsetDeviceViewType, RowViewer, false> Op;
          Op::launch(alpha, *A_, X, beta, Y, offRankOffsets, execSpace);
          return;
        }
        case Teuchos::ETransp::CONJ_TRANS: {
          typedef SpmvMvTransFunctor<OffsetDeviceViewType, RowViewer, true> Op;
          Op::launch(alpha, *A_, X, beta, Y, offRankOffsets, execSpace);
          return;
        }

        default:
          throw std::runtime_error("unexpected Teuchos::ETransp mode in off-rank SpMV");
      }
    }

    /*! \brief
        Complement of \c applyRemoteColumns(). Only contribute matrix entries NOT specified in offRankOffsets,
        which should be populated by \c CrsGraph::getLocalOffRankOffsets

        \tparam OffsetDeviceViewType should be the CrsMatrix::crs_graph_type::offset_device_view_type of the \c Tpetra::CrsMatrix that owns this LocalCrsMatrixOperator
    
        cwp 05 Apr 2022
        applyRemoteColumns() with applyLocalColumns() shall have the same effect as apply()
    */ 
    template<typename OffsetDeviceViewType>
    void
    applyLocalColumns (execution_space &execSpace,
           Kokkos::View<const mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > X,
           Kokkos::View<mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Y,
           const Teuchos::ETransp mode,
           const mv_scalar_type alpha,
           const mv_scalar_type beta,
           const OffsetDeviceViewType &offRankOffsets) const {
      typedef OnRankRowViewer<OffsetDeviceViewType> RowViewer;
      switch(mode) {
        case Teuchos::ETransp::TRANS: {
          typedef SpmvMvTransFunctor<OffsetDeviceViewType, RowViewer, false> Op;
          Op::launch(alpha, *A_, X, beta, Y, offRankOffsets, execSpace);
          return;
        }
        case Teuchos::ETransp::NO_TRANS: {
          typedef SpmvMvFunctor<OffsetDeviceViewType, RowViewer, false> Op;
          Op::launch(alpha, *A_, X, beta, Y, offRankOffsets, execSpace);
          return;
        }
        case Teuchos::ETransp::CONJ_TRANS: {
          typedef SpmvMvTransFunctor<OffsetDeviceViewType, RowViewer, true> Op;
          Op::launch(alpha, *A_, X, beta, Y, offRankOffsets, execSpace);
          return;
        }
        default:
          throw std::runtime_error("unexpected Teuchos::ETransp mode in on-rank SpMV");
      }
    }


           

    bool hasTransposeApply () const override;

    const local_matrix_device_type& getLocalMatrixDevice () const;

  private:
    std::shared_ptr<local_matrix_device_type> A_;
    local_cusparse_matrix_type A_cusparse;
    const bool have_A_cusparse;
  };

} // namespace Tpetra

#endif // TPETRA_LOCALCRSMATRIXOPERATOR_DECL_HPP
