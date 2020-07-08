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
    using local_matrix_type =
      KokkosSparse::CrsMatrix<matrix_scalar_type,
                              local_ordinal_type,
                              device_type>;
  private:
    using local_graph_type = typename local_matrix_type::StaticCrsGraphType;

  public:
    LocalCrsMatrixOperator (const std::shared_ptr<local_matrix_type>& A);
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

    const local_matrix_type& getLocalMatrix () const;

  private:
    std::shared_ptr<local_matrix_type> A_;
  };

} // namespace Tpetra

#endif // TPETRA_LOCALCRSMATRIXOPERATOR_DECL_HPP
