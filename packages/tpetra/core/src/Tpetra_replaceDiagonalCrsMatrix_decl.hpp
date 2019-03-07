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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_REPLACEDIAGONALCRSMATRIX_DECL_HPP
#define TPETRA_REPLACEDIAGONALCRSMATRIX_DECL_HPP

/// \file Tpetra_replaceDiagonalCrsMatrix_decl.hpp
/// \brief Declaration of Tpetra::repalceDiagonalCrsMatrix

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
/// \paran[in] newDiag Tpetra::Vector with new values for the diagonal
///
/// \return Local number of successfully replaced diagonal entries
template<class SC, class LO, class GO, class NT>
LO
replaceDiagonalCrsMatrix(::Tpetra::CrsMatrix<SC, LO, GO, NT>& matrix,
    const ::Tpetra::Vector<SC, LO, GO, NT>& newDiag);

} // namespace Tpetra

#endif // #ifndef TPETRA_REPLACEDIAGONALCRSMATRIX_DECL_HPP
