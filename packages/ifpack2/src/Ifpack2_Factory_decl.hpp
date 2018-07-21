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

#ifndef IFPACK2_FACTORY_DECL_HPP
#define IFPACK2_FACTORY_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_Factory.hpp"

// FIXME (mfh 28 Jul 2015) I would very much not like to include ANY
// specific preconditioner header files here.  However, these are
// unfortunately necessary for the largely useless clone() method.
#include "Ifpack2_Chebyshev.hpp"
#include "Ifpack2_RILUK.hpp"
#include "Ifpack2_Experimental_RBILUK.hpp"

#include <type_traits>

namespace Ifpack2 {

//! \c true if the specified preconditioner type supports nonsymmetric matrices, else false.
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

  /// \brief Clones a preconditioner for a different node type from an
  ///   Ifpack2 RILUK or Chebyshev preconditioner
  template<class InputMatrixType, class OutputMatrixType>
  static
  Teuchos::RCP<Preconditioner<typename OutputMatrixType::scalar_type,
                              typename OutputMatrixType::local_ordinal_type,
                              typename OutputMatrixType::global_ordinal_type,
                              typename OutputMatrixType::node_type> >
  clone (const Teuchos::RCP<Preconditioner<typename InputMatrixType::scalar_type,
                                           typename InputMatrixType::local_ordinal_type,
                                           typename InputMatrixType::global_ordinal_type,
                                           typename InputMatrixType::node_type> >& prec,
         const Teuchos::RCP<const OutputMatrixType>& matrix,
         const Teuchos::ParameterList& params = Teuchos::ParameterList ())
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;

    // FIXME (mfh 09 Nov 2013) The code below assumes that the old and
    // new scalar, local ordinal, and global ordinal types are the same.

    typedef typename InputMatrixType::scalar_type scalar_type;
    typedef typename InputMatrixType::local_ordinal_type local_ordinal_type;
    typedef typename InputMatrixType::global_ordinal_type global_ordinal_type;
    typedef typename InputMatrixType::node_type old_node_type;
    typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type,
      global_ordinal_type, old_node_type> input_row_matrix_type;

    static_assert (std::is_same<typename OutputMatrixType::scalar_type, scalar_type>::value,
                   "Input and output scalar_type must be the same.");
    static_assert (std::is_same<typename OutputMatrixType::local_ordinal_type, local_ordinal_type>::value,
                   "Input and output local_ordinal_type must be the same.");
    static_assert (std::is_same<typename OutputMatrixType::global_ordinal_type, global_ordinal_type>::value,
                   "Input and output global_ordinal_type must be the same.");
    typedef typename OutputMatrixType::node_type new_node_type;
    typedef Preconditioner<scalar_type, local_ordinal_type,
      global_ordinal_type, new_node_type> output_prec_type;

    // FIXME (mfh 09 Nov 2013) The code below only knows how to clone
    // three different kinds of preconditioners.  This is probably because
    // only two subclasses of Preconditioner implement a clone() method.

    RCP<output_prec_type> new_prec;
    RCP<Chebyshev<input_row_matrix_type> > chebyPrec =
      rcp_dynamic_cast<Chebyshev<input_row_matrix_type> > (prec);
    if (! chebyPrec.is_null ()) {
      new_prec = chebyPrec->clone (matrix, params);
      return new_prec;
    }
    RCP<RILUK<input_row_matrix_type> > luPrec;
    luPrec = rcp_dynamic_cast<RILUK<input_row_matrix_type> > (prec);
    if (luPrec != null) {
      new_prec = luPrec->clone (matrix);
      return new_prec;
    }
    RCP<Experimental::RBILUK<input_row_matrix_type> > rbilukPrec;
    rbilukPrec = rcp_dynamic_cast<Experimental::RBILUK<input_row_matrix_type> > (prec);
    if (rbilukPrec != null) {
      new_prec = rbilukPrec->clone (matrix);
      return new_prec;
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::logic_error, "Ifpack2::Factory::clone: Not implemented for the "
       "current preconditioner type.  The only supported types thus far are "
       "Chebyshev, RILUK, and RBILUK.");
  }
};

} // namespace Ifpack2

#endif // IFPACK2_FACTORY_DECL_HPP
