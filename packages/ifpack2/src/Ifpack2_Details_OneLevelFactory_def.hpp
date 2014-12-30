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

#ifndef IFPACK2_DETAILS_ONELEVELFACTORY_DEF_HPP
#define IFPACK2_DETAILS_ONELEVELFACTORY_DEF_HPP

#include "Ifpack2_Chebyshev.hpp"
#include "Ifpack2_Details_DenseSolver.hpp"
#include "Ifpack2_Diagonal.hpp"
#include "Ifpack2_IdentitySolver.hpp"
#include "Ifpack2_ILUT.hpp"
#include "Ifpack2_Relaxation.hpp"
#include "Ifpack2_RILUK.hpp"
#include "Ifpack2_Krylov.hpp"
#include "Ifpack2_BlockRelaxation.hpp"
#include "Ifpack2_DenseContainer.hpp"
#include "Ifpack2_Details_OneLevelFactory.hpp"

#ifdef HAVE_IFPACK2_AMESOS2
#  include "Ifpack2_Details_Amesos2Wrapper.hpp"
#endif // HAVE_IFPACK2_AMESOS2

namespace Ifpack2 {
namespace Details {

template<class MatrixType>
Teuchos::RCP<typename OneLevelFactory<MatrixType>::prec_type>
OneLevelFactory<MatrixType>::create (const std::string& precType,
                                     const Teuchos::RCP<const row_matrix_type>& matrix) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<prec_type> prec;

  // precTypeUpper is the upper-case version of precType.
  std::string precTypeUpper (precType);
  if (precTypeUpper.size () > 0) {
    std::locale locale;
    for (size_t k = 0; k < precTypeUpper.size (); ++k) {
      precTypeUpper[k] = std::toupper<char> (precTypeUpper[k], locale);
    }
  }

  if (precTypeUpper == "CHEBYSHEV") {
    // We have to distinguish Ifpack2::Chebyshev from its
    // implementation class Ifpack2::Details::Chebyshev.
    prec = rcp (new ::Ifpack2::Chebyshev<MatrixType> (matrix));
  }
  else if (precTypeUpper == "DENSE" || precTypeUpper == "LAPACK") {
    prec = rcp (new Details::DenseSolver<MatrixType> (matrix));
  }
  else if (precTypeUpper == "AMESOS2") {
#ifdef HAVE_IFPACK2_AMESOS2
    prec = rcp (new Details::Amesos2Wrapper<MatrixType> (matrix));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "Ifpack2::Details::OneLevelFactory: "
      "You may not ask for the preconditioner \"AMESOS2\" unless "
      "you have built Trilinos with the Amesos2 package enabled.");
#endif // HAVE_IFPACK2_AMESOS2
  }
  else if (precTypeUpper == "DIAGONAL") {
    prec = rcp (new Diagonal<MatrixType> (matrix));
  }
  else if (precTypeUpper == "ILUT") {
    prec = rcp (new ILUT<MatrixType> (matrix));
  }
  else if (precTypeUpper == "RELAXATION") {
    prec = rcp (new Relaxation<MatrixType> (matrix));
  }
  else if (precTypeUpper == "RILUK") {
    prec = rcp (new RILUK<MatrixType> (matrix));
  }
  else if (precTypeUpper == "KRYLOV") {
    prec = rcp (new Krylov<MatrixType> (matrix));
  }
  else if (precTypeUpper == "BLOCK_RELAXATION" ||
           precTypeUpper == "BLOCK RELAXATION" ||
           precTypeUpper == "BLOCKRELAXATION" ) {
    // FIXME (mfh 22 May 2014) We would prefer to have the choice of
    // dense or sparse blocks (the "container type") be a run-time
    // decision.  This will require refactoring BlockRelaxation so
    // that the "container type" is not a template parameter.  For
    // now, we default to use dense blocks.
    typedef DenseContainer<MatrixType, scalar_type> container_type;
    prec = rcp (new BlockRelaxation<MatrixType, container_type> (matrix));
  }
  else if (precTypeUpper == "IDENTITY" || precTypeUpper == "IDENTITY_SOLVER") {
    prec = rcp (new IdentitySolver<MatrixType> (matrix));
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "Ifpack2::Details::OneLevelFactory::create: "
      "Invalid preconditioner type \"" << precType << "\".");
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    prec.is_null (), std::logic_error, "Ifpack2::Details::OneLevelFactory::"
    "create: Return value is null right before return.  This should never "
    "happen.  Please report this bug to the Ifpack2 developers.");
  return prec;
}

} // namespace Details
} // namespace Ifpack2

#define IFPACK2_DETAILS_ONELEVELFACTORY_INSTANT(S,LO,GO,N)              \
  template class Ifpack2::Details::OneLevelFactory< Tpetra::CrsMatrix<S, LO, GO, N> >; \
  template class Ifpack2::Details::OneLevelFactory< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_DETAILS_ONELEVELFACTORY_DEF_HPP
