// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REPLACEDIAGONALCRSMATRIX_DECL_HPP
#define TPETRA_REPLACEDIAGONALCRSMATRIX_DECL_HPP

/// \file Tpetra_replaceDiagonalCrsMatrix_decl.hpp
/// \brief Declaration of Tpetra::replaceDiagonalCrsMatrix

#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"

namespace Tpetra {

/// \brief Replace diagonal entries of an input Tpetra::CrsMatrix \c matrix
///   with values given in \c newDiag
///
/// \tparam Scalar The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LocalOrdinal The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GlobalOrdinal The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam Node The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// \param[in/out] matrix Tpetra::CrsMatrix to be modified
/// \param[in] newDiag Tpetra::Vector with new values for the diagonal; must have same Tpetra::Map as matrix's rowMap
///
///
/// \return Local number of successfully replaced diagonal entries
template<class SC, class LO, class GO, class NT>
LO
replaceDiagonalCrsMatrix(::Tpetra::CrsMatrix<SC, LO, GO, NT>& matrix,
    const ::Tpetra::Vector<SC, LO, GO, NT>& newDiag);

} // namespace Tpetra

#endif // #ifndef TPETRA_REPLACEDIAGONALCRSMATRIX_DECL_HPP
