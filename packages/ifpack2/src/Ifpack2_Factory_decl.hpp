/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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
  - "KRYLOV": returns an instance of Ifpack2::Krylov.

The following fragment of code shows the basic usage of this class.
\code
#include "Ifpack2_Factory.hpp"
// ...

using Teuchos::ParameterList;
using Teuchos::RCP;
typedef double Scalar;
typedef int    LocalOrdinal;
typedef int    GlobalOrdinal;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> prec_type;
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
          const Teuchos::RCP<const MatrixType>& matrix);

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
          const int overlap);

  //! Clones a preconditioner for a different node type from an Ifpack2 RILUK or Chebyshev preconditioner
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
         const Teuchos::ParameterList& params = Teuchos::ParameterList ());
};

} // namespace Ifpack2

#endif // IFPACK2_FACTORY_DECL_HPP
