// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_ROWMATRIX_HPP
#define IFPACK2_DETAILS_ROWMATRIX_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"

namespace Ifpack2 {
namespace Details {

/// \class RowMatrix
/// \brief All Ifpack2 implementations of Tpetra::RowMatrix must
///   inherit from this class.
/// \tparam MatrixType Tpetra::RowMatrix specialization.
///
/// \warning This class is an implementation detail of Ifpack2.  Users
///   should not rely on its interface.
///
/// This class exists to facilitate Tpetra interface changes.  See
/// e.g., GitHub Issue #2630.
template<class MatrixType>
class RowMatrix :
    public Tpetra::RowMatrix<typename MatrixType::scalar_type,
                             typename MatrixType::local_ordinal_type,
                             typename MatrixType::global_ordinal_type,
                             typename MatrixType::node_type> {
public:
  //! \name Typedefs
  //@{
  using scalar_type = typename MatrixType::scalar_type;
  using local_ordinal_type = typename MatrixType::local_ordinal_type;
  using global_ordinal_type = typename MatrixType::global_ordinal_type;
  using node_type = typename MatrixType::node_type;

  //@}
  //! \name Destructor
  //@{

  //! Destructor (virtual for memory safety of derived classes)
  virtual ~RowMatrix () = default;

  //@}

};

} // namespace Details
} // namespace Ifpack2

#endif /* IFPACK2_DETAILS_ROWMATRIX_HPP */
