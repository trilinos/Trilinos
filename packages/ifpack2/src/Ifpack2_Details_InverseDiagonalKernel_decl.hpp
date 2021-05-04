/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

/// \file Ifpack2_Details_DiagonalKernel_decl.hpp
/// \brief Declaration of Ifpack2::Details::DiagonalKernel
///
/// \note This file and its contents are implementation details of
///   Ifpack2.

#ifndef IFPACK2_DETAILS_INVERSEDIAGONALKERNEL_DECL_HPP
#define IFPACK2_DETAILS_INVERSEDIAGONALKERNEL_DECL_HPP

#include "Ifpack2_config.h"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_MultiVector_fwd.hpp"
#include "Tpetra_Operator_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"
#include "Tpetra_Export_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Teuchos_RCP.hpp"
#include <memory>

namespace Ifpack2 {
namespace Details {

/// \brief Compute scaled damped residual for Chebyshev.
///
/// This is an implementation detail of Ifpack2::Chebyshev.  Given a
/// linear system A*X=B and an "inverse diagonal" matrix (stored as a
/// vector) D_inv, it computes a "scaled damped" residual vector W :=
/// alpha*D_inv*(B-A*X) + beta*W.
///
/// \tparam TpetraOperatorType Specialization of Tpetra::Operator.
///
/// \note To Ifpack2 developers: We can't fuse this with X := X + W,
///   because data dependencies in the input X are not elementwise
///   (unless A is diagonal).
template<class TpetraOperatorType>
class InverseDiagonalKernel {
private:
  using SC = typename TpetraOperatorType::scalar_type;
  using LO = typename TpetraOperatorType::local_ordinal_type;
  using GO = typename TpetraOperatorType::global_ordinal_type;
  using NT = typename TpetraOperatorType::node_type;

  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using multivector_type = Tpetra::MultiVector<SC, LO, GO, NT>;
  using operator_type = Tpetra::Operator<SC, LO, GO, NT>;
  using vector_type = Tpetra::Vector<SC, LO, GO, NT>;
  using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;
  using device_type = typename NT::device_type;
  using offset_type = Kokkos::View<size_t*, device_type>;

public:
  InverseDiagonalKernel (const Teuchos::RCP<const operator_type>& A);

  void
  setMatrix (const Teuchos::RCP<const operator_type>& A);

  void
  compute (vector_type& D_inv,
           bool do_l1, magnitude_type L1Eta,
           bool fixTinyDiagEntries, magnitude_type MinDiagonalValue);

private:
  Teuchos::RCP<const operator_type> A_op_;
  Teuchos::RCP<const crs_matrix_type> A_crs_;
  offset_type offsets_;

};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_INVERSEINVERSEDIAGONALKERNEL_DECL_HPP
