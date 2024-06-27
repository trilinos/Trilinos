// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_GETCRSMATRIX_HPP
#define IFPACK2_DETAILS_GETCRSMATRIX_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
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

//Helper to get A as a Tpetra::BlockCrsMatrix, if that is possible cheaply and without copying.
//This is just a dynamic cast.
template<typename SC, typename LO, typename GO, typename NO>
Teuchos::RCP<const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>> getBcrsMatrix(const Teuchos::RCP<const Tpetra::RowMatrix<SC, LO, GO, NO>>& A)
{
  using bcrs_matrix_type = Tpetra::BlockCrsMatrix<SC, LO, GO, NO>;
  auto Abcrs = Teuchos::rcp_dynamic_cast<const bcrs_matrix_type>(A);
  if(!Abcrs.is_null())
    return Abcrs;

  return Teuchos::null;
}

} // namespace Details
} // namespace Ifpack2

#endif
