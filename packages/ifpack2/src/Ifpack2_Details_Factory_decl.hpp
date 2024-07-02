// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_FACTORY_DECL_HPP
#define IFPACK2_DETAILS_FACTORY_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"

namespace Ifpack2 {
namespace Details {


template<class SC, class LO, class GO, class NT>
class Factory {
public:
  typedef Tpetra::RowMatrix<SC, LO, GO, NT> row_matrix_type;
  typedef ::Ifpack2::Preconditioner<SC, LO, GO, NT> prec_type;

  /// \brief Create an instance of Ifpack2::Preconditioner given the
  ///   string name of the preconditioner type.
  ///
  /// \param precType [in] Name of preconditioner type to be created.
  /// \param matrix [in] Matrix used to define the preconditioner
  ///
  /// Throw an exception if the preconditioner with that input name
  /// does not exist.  Otherwise, return a newly created
  /// preconditioner object.
  Teuchos::RCP<prec_type>
  create (const std::string& precType,
          const Teuchos::RCP<const row_matrix_type>& matrix);

  /** \brief Create an instance of Ifpack2::Preconditioner given the
   *   string name of the preconditioner type.
   *
   * \warning This version of the constructor is DEPRECATED, because
   *   the single-argument version suffices; users may specify the
   *   overlap level via the "schwarz: overlap level" parameter.
   *
   * \param precType [in] Name of preconditioner type to be created.
   * \param matrix [in] Matrix used to define the preconditioner
   * \param overlap (in) Specified overlap; defaults to 0.
   *
   * Throw an exception if the preconditioner with that input name
   * does not exist.  Otherwise, return a newly created preconditioner
   * object.
   */
  Teuchos::RCP<prec_type>
  create (const std::string& precType,
          const Teuchos::RCP<const row_matrix_type>& matrix,
          const int overlap);

  bool
  isSupported (const std::string& precType);
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_FACTORY_DECL_HPP
