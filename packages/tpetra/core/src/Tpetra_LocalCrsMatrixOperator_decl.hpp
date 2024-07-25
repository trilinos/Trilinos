// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

    bool hasTransposeApply () const override;

    const local_matrix_device_type& getLocalMatrixDevice () const;

  private:
    std::shared_ptr<local_matrix_device_type> A_;
    local_cusparse_matrix_type A_cusparse;
    const bool have_A_cusparse;
  };

} // namespace Tpetra

#endif // TPETRA_LOCALCRSMATRIXOPERATOR_DECL_HPP
