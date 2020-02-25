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
