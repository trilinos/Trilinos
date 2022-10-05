/*@HEADER
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_DETAILS_GETCRSMATRIX_HPP
#define IFPACK2_DETAILS_GETCRSMATRIX_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Ifpack2_Details_AdditiveSchwarzFilter.hpp"

namespace Ifpack2 {
namespace Details {

//Helper to get A as a Tpetra::CrsMatrix, if that is possible cheaply and without copying.
//In the simplest case (when A is already a CrsMatrix), this is just a dynamic cast.
template<typename SC, typename LO, typename GO, typename NO>
Teuchos::RCP<const Tpetra::CrsMatrix<SC, LO, GO, NO>> getCrsMatrix(const Teuchos::RCP<const Tpetra::RowMatrix<SC, LO, GO, NO>>& A)
{
  using row_matrix_type = Tpetra::RowMatrix<SC, LO, GO, NO>;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NO>;
  auto Acrs = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(A);
  if(!Acrs.is_null())
    return Acrs;
  auto Aasf = Teuchos::rcp_dynamic_cast<const AdditiveSchwarzFilter<row_matrix_type>>(A);
  if(!Aasf.is_null())
    return Aasf->getFilteredMatrix();
  return Teuchos::null;
}

} // namespace Details
} // namespace Ifpack2

#endif
