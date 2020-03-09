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

#ifndef TPETRA_DETAILS_LOCALCRSMATRIXOPERATORWITHSETUP_HPP
#define TPETRA_DETAILS_LOCALCRSMATRIXOPERATORWITHSETUP_HPP

#include "Tpetra_Details_LocalCrsMatrixOperatorWithSetup_fwd.hpp"
#include "Tpetra_LocalCrsMatrixOperator.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace Tpetra {
namespace Details {

  /// \class LocalCrsMatrixOperatorWithSetup
  /// \brief Subclass of LocalCrsMatrixOperator that has a set-up
  ///   phase ("resume fill") and a ready-for-apply phase ("fill
  ///   complete"), that users may freely toggle.
  ///
  /// This object has two states: "is fill resumed" and "is fill
  /// completed."  It is in the "is fill resumed" state after calling
  /// resumeFill, but before the next call to fillComplete.  It is in
  /// the "is fill completed" state after calling fillComplete, but
  /// before the next call to resumeFill.  After creating this object,
  /// but before the next call to fillComplete, this object is in the
  /// "is fill resumed" state.  Calls to resumeFill or fillComplete
  /// are idempotent.
  ///
  /// <ol>
  /// <li> After creating this object, you may not change the matrix's
  ///      allocations or modify its graph structure.
  /// </li>
  /// <li> After calling \c resumeFill, you may not call \c apply
  ///      unless you first call \c fillComplete.
  /// </li>
  /// <li> After calling \c fillComplete, you may not modify the
  ///      matrix's values unless you first call \c resumeFill.
  /// </li>
  /// </ol>
  ///
  /// \tparam MultiVectorScalar The type of the entries of the input
  ///   and output (multi)vectors.
  /// \tparam MatrixScalar The type of the entries of the sparse matrix.
  /// \tparam Device The Kokkos Device type; must be a specialization
  ///   of Kokkos::Device.
  template<class MultiVectorScalar, class MatrixScalar, class Device>
  class LocalCrsMatrixOperatorWithSetup :
    public ::Tpetra::LocalCrsMatrixOperator<MultiVectorScalar,
                                            MatrixScalar, Device>
  {
  private:
    using base_type =
      ::Tpetra::LocalCrsMatrixOperator<MultiVectorScalar,
                                       MatrixScalar, Device>;
    using LO = ::Tpetra::Details::DefaultTypes::local_ordinal_type;

  public:
    using local_matrix_type = typename base_type::local_matrix_type;

    LocalCrsMatrixOperatorWithSetup(
      const std::shared_ptr<local_matrix_type>& A) : base_type(A)
    {}
    virtual ~LocalCrsMatrixOperatorWithSetup() override = default;

    /// \brief Tell this object that you intend to modify the values
    ///   in the matrix.
    ///
    /// After creating this object, you are not allowed to change the
    /// matrix's allocations or modify its graph structure.
    virtual void resumeFill() {}

    /// \brief Tell this object that you are done modifying the values
    ///   in the matrix for now, and are ready to use the base class
    ///   method apply().
    ///
    /// After creating this object, you are not allowed to change the
    /// matrix's allocations or modify its graph structure.
    virtual void fillComplete() {}

    /// \brief Tell this object the minimum and maximum number of
    ///   matrix entries per row in the sparse matrix given to its
    ///   constructor.
    ///
    /// This call will take effect on the next call to fillComplete.
    ///
    /// This is useful for telling the computational kernel how to
    /// load-balance thread parallelism over rows of the matrix.
    virtual void
    setMinMaxNumberOfEntriesPerRow(const LO /* minNumEnt */,
                                   const LO /* maxNumEnt */)
    {}
  };

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_LOCALCRSMATRIXOPERATORWITHSETUP_HPP
