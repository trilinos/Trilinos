// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_FACTORY_DECL_HPP
#define IFPACK2_FACTORY_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_Factory.hpp"

#include "Ifpack2_Chebyshev.hpp"
#include "Ifpack2_RILUK.hpp"
#include "Ifpack2_Experimental_RBILUK.hpp"

#include <type_traits>

namespace Ifpack2 {

//! \C true if the specified preconditioner type supports nonsymmetric matrices, else false.
bool supportsUnsymmetric (const std::string& prec_type);

/*!
\class Factory
\brief "Factory" for creating Ifpack2 preconditioners.

This class' create() method lets users create an instance of any Ifpack2 preconditioner.

The create() method has three arguments:
  - a string, indicating the type of preconditioner to compute;
  - a pointer to a Tpetra::RowMatrix, representing the matrix to be
    used to define the preconditioner;
  - an optional integer (defaults to 0), that specifies the amount of
    overlap between processes (if the input matrix is distributed over
    multiple processes).

The first argument can assume the following values:
  - "DIAGONAL": returns an instance of Ifpack2::Diagonal.
  - "RELAXATION": returns an instance of Ifpack2::Relaxation.
  - "CHEBYSHEV": returns an instance of Ifpack2::Chebyshev (overlap is ignored).
  - "ILUT": returns an instance of Ifpack2::ILUT.
  - "MDF": returns an instance of Ifpack2::MDF.
  - "RILUK": returns an instance of Ifpack2::RILUK.
  - "RBILUK": returns an instance of Ifpack2::Experimental::RBILUK.

The following fragment of code shows the basic usage of this class.
\code
#include "Ifpack2_Factory.hpp"
// ...

using Teuchos::ParameterList;
using Teuchos::RCP;
typedef double Scalar;
typedef Tpetra::CrsMatrix<Scalar> crs_matrix_type;
typedef Ifpack2::Preconditioner<Scalar> prec_type;
// ...

Ifpack2::Factory factory;

RCP<crs_matrix_type> A;
// ... fill A and call fillComplete() ...

// Use ILUT (incomplete LU with thresholding) on each process.
const std::string precType = "ILUT";
RCP<prec_type> prec = factory.create (precType, A);

ParameterList params;
params.set ("fact: ilut level-of-fill", 5.0); // ILUT(fill=5, drop=0)

prec->setParameters (params);
prec->initialize ();
prec->compute ();

// now prec can be used as a preconditioner
\endcode
*/
class Factory {
public:
  /** \brief Create an instance of Ifpack2_Preconditioner given the string
   * name of the preconditioner type.
   *
   * \param precType [in] Name of preconditioner type to be created.
   * \param matrix [in] Matrix used to define the preconditioner
   *
   * Throw an exception if the preconditioner with that input name
   * does not exist.  Otherwise, return a newly created preconditioner
   * object.
   */
  template<class MatrixType>
  static
  Teuchos::RCP<Preconditioner<typename MatrixType::scalar_type,
                              typename MatrixType::local_ordinal_type,
                              typename MatrixType::global_ordinal_type,
                              typename MatrixType::node_type> >
  create (const std::string& precType,
          const Teuchos::RCP<const MatrixType>& matrix)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_implicit_cast;
    typedef typename MatrixType::scalar_type SC;
    typedef typename MatrixType::local_ordinal_type LO;
    typedef typename MatrixType::global_ordinal_type GO;
    typedef typename MatrixType::node_type NT;
    typedef Tpetra::RowMatrix<SC, LO, GO, NT> row_matrix_type;

    RCP<const row_matrix_type> A;
    if (! matrix.is_null ()) {
      A = rcp_implicit_cast<const row_matrix_type> (matrix);
    }
    Ifpack2::Details::Factory<SC, LO, GO, NT> factory;
    return factory.create (precType, A);
  }

  /** \brief Create an instance of Ifpack2_Preconditioner given the string
   * name of the preconditioner type.
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
  template<class MatrixType>
  static
  Teuchos::RCP<Preconditioner<typename MatrixType::scalar_type,
                              typename MatrixType::local_ordinal_type,
                              typename MatrixType::global_ordinal_type,
                              typename MatrixType::node_type> >
  create (const std::string& precType,
          const Teuchos::RCP<const MatrixType>& matrix,
          const int overlap)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_implicit_cast;
    typedef typename MatrixType::scalar_type SC;
    typedef typename MatrixType::local_ordinal_type LO;
    typedef typename MatrixType::global_ordinal_type GO;
    typedef typename MatrixType::node_type NT;
    typedef Tpetra::RowMatrix<SC, LO, GO, NT> row_matrix_type;

    RCP<const row_matrix_type> A;
    if (! matrix.is_null ()) {
      A = rcp_implicit_cast<const row_matrix_type> (matrix);
    }
    Ifpack2::Details::Factory<SC, LO, GO, NT> factory;
    return factory.create (precType, A, overlap);
  }

  template<class MatrixType>
  static
  bool
  isSupported (const std::string& precType)
  {
    typedef typename MatrixType::scalar_type SC;
    typedef typename MatrixType::local_ordinal_type LO;
    typedef typename MatrixType::global_ordinal_type GO;
    typedef typename MatrixType::node_type NT;

    Ifpack2::Details::Factory<SC, LO, GO, NT> factory;
    return factory.isSupported (precType);
  }

};

} // namespace Ifpack2

#endif // IFPACK2_FACTORY_DECL_HPP
