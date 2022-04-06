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
#include "KokkosSparse_CrsMatrix.hpp"
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

    // cwp 05 Apr 2022
    // A functor that does the on-rank part of a local SpMV
    // KokkosKernels does not currently have a 4-array CSR
    template<typename OffsetDeviceViewType>
    class OnRankSpmvFunctor {

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

    public:
      OnRankSpmvFunctor(const MultiVectorScalar &alpha, 
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
      KOKKOS_INLINE_FUNCTION void operator()(TagNonTrans, const size_t i) const {

        // beta * Y
        for (size_t k = 0; k < Y_.extent(1); ++k) {
          Y_(i,k) = beta_ * Y_(i,k); 
        }       

        // + alpha A X
        for (size_type ji = A_.graph.row_map(i); ji < offRankOffsets_(i); ++ji) {
          value_type A_ij = A_.values(ji);
          ordinal_type j = A_.graph.entries(ji);
          for (size_t k = 0; k < Y_.extent(1); ++k) {
            Y_(i,k) += alpha_ * A_ij * X_(j, k); 
          }
        }
      }

      KOKKOS_INLINE_FUNCTION void operator()(TagTrans, const size_t i) const {
        #warning trans unimplemented
      }

      KOKKOS_INLINE_FUNCTION void operator()(TagConjTrans, const size_t i) const {
        // typedef Kokkos::Details::ArithTraits<value_type> KAT; 
        #warning conj trans unimplemented
      }

      /// \brief Kokkos dispatch of non-transpose
      void launch(TagNonTrans) {
        Kokkos::parallel_for(Kokkos::RangePolicy<TagNonTrans>(0, A_.numRows()), *this);
      }

      /// \brief Kokkos dispatch of transpose
      void launch(TagTrans) {
        Kokkos::parallel_for(Kokkos::RangePolicy<TagTrans>(0, A_.numRows()), *this);
      }

      /// \brief Kokkos dispatch of conjugate transpose
      void launch(TagConjTrans) {
        Kokkos::parallel_for(Kokkos::RangePolicy<TagConjTrans>(0, A_.numRows()), *this);
      }
    };

    // cwp 06 Apr 2022
    // A functor that does the off-rank part of a local SpMV
    // KokkosKernels does not currently have a 4-array CSR
    template<typename OffsetDeviceViewType>
    class OffRankSpmvFunctor {

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

    public:
      OffRankSpmvFunctor(const MultiVectorScalar &alpha, 
      const local_matrix_device_type &A, 
      x_type &X, 
      const MultiVectorScalar &beta, 
      y_type &Y,
      const OffsetDeviceViewType &offRankOffsets) 
        : alpha_(alpha), A_(A), X_(X), beta_(beta), Y_(Y), offRankOffsets_(offRankOffsets) {}

      /*! \brief contribution of row i

          This is only the remote offsets, so the offsets we want are between
          offRankOffsets_(i) and entries(i+1)
      */
      KOKKOS_INLINE_FUNCTION void operator()(TagNonTrans, const size_t i) const {

        // beta * Y
        for (size_t k = 0; k < Y_.extent(1); ++k) {
          Y_(i,k) = beta_ * Y_(i,k); 
        }       

        // + alpha A x
        for (size_type ji = offRankOffsets_(i); ji < A_.graph.row_map(i+1); ++ji) {
          value_type A_ij = A_.values(ji);
          ordinal_type j = A_.graph.entries(ji);
          for (size_t k = 0; k < Y_.extent(1); ++k) {
            Y_(i,k) += alpha_ * A_ij * X_(j, k); 
          }
        }
      }

      KOKKOS_INLINE_FUNCTION void operator()(TagTrans, const size_t i) const {
        #warning trans unimplemented
      }

      KOKKOS_INLINE_FUNCTION void operator()(TagConjTrans, const size_t i) const {
        // typedef Kokkos::Details::ArithTraits<value_type> KAT; 
        #warning conj trans unimplemented
      }

      /// \brief Kokkos dispatch of non-transpose
      void launch(TagNonTrans) {
        Kokkos::parallel_for(Kokkos::RangePolicy<TagNonTrans>(0, A_.numRows()), *this);
      }

      /// \brief Kokkos dispatch of transpose
      void launch(TagTrans) {
        Kokkos::parallel_for(Kokkos::RangePolicy<TagTrans>(0, A_.numRows()), *this);
      }

      /// \brief Kokkos dispatch of conjugate transpose
      void launch(TagConjTrans) {
        Kokkos::parallel_for(Kokkos::RangePolicy<TagConjTrans>(0, A_.numRows()), *this);
      }
    };



  public:
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
    applyRemoteColumns (Kokkos::View<const mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > X,
           Kokkos::View<mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Y,
           const Teuchos::ETransp mode,
           const mv_scalar_type alpha,
           const mv_scalar_type beta,
           const OffsetDeviceViewType &offRankOffsets) const {

      typedef OffRankSpmvFunctor<OffsetDeviceViewType> ORSF;
      ORSF orsf(alpha, *A_, X, beta, Y, offRankOffsets);
      switch(mode) {
        case Teuchos::ETransp::TRANS: {
          orsf.launch(typename ORSF::TagTrans{});
          return;
        }
        case Teuchos::ETransp::NO_TRANS: {
          orsf.launch(typename ORSF::TagNonTrans{});
          return;
        }
        case Teuchos::ETransp::CONJ_TRANS: {
          orsf.launch(typename ORSF::TagConjTrans{});
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
    applyLocalColumns (Kokkos::View<const mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > X,
           Kokkos::View<mv_scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Y,
           const Teuchos::ETransp mode,
           const mv_scalar_type alpha,
           const mv_scalar_type beta,
           const OffsetDeviceViewType &offRankOffsets) const {

      typedef OnRankSpmvFunctor<OffsetDeviceViewType> ORSF;
      ORSF orsf(alpha, *A_, X, beta, Y, offRankOffsets);
      switch(mode) {
        case Teuchos::ETransp::TRANS: {
          orsf.launch(typename ORSF::TagTrans{});
          return;
        }
        case Teuchos::ETransp::NO_TRANS: {
          orsf.launch(typename ORSF::TagNonTrans{});
          return;
        }
        case Teuchos::ETransp::CONJ_TRANS: {
          orsf.launch(typename ORSF::TagConjTrans{});
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
