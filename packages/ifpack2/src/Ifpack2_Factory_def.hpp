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

#ifndef IFPACK2_FACTORY_DEF_HPP
#define IFPACK2_FACTORY_DEF_HPP

#include "Ifpack2_Relaxation.hpp"
#include "Ifpack2_Diagonal.hpp"
#include "Ifpack2_Chebyshev.hpp"
#include "Ifpack2_RILUK.hpp"
#include "Ifpack2_ILUT.hpp"
#include "Ifpack2_Krylov.hpp"
#include "Ifpack2_AdditiveSchwarz.hpp"
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
#  include "Ifpack2_SupportGraph.hpp"
#endif


namespace Ifpack2 {

template<class MatrixType>
Teuchos::RCP<Preconditioner<typename MatrixType::scalar_type,
                            typename MatrixType::local_ordinal_type,
                            typename MatrixType::global_ordinal_type,
                            typename MatrixType::node_type> >
Factory::create (const std::string& precType,
                 const Teuchos::RCP<const MatrixType>& matrix,
                 const int overlap)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;
  typedef Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type, node_type> prec_base_type;

  RCP<prec_base_type> prec;

  // precTypeUpper is the upper-case version of precType.
  std::string precTypeUpper (precType);
  if (precTypeUpper.size () > 0) {
    std::locale locale;
    for (size_t k = 0; k < precTypeUpper.size (); ++k) {
      precTypeUpper[k] = std::toupper<char> (precTypeUpper[k], locale);
    }
  }

  const bool one_mpi_rank = (matrix->getComm ()->getSize () == 1);

  if (precTypeUpper == "ILUT") {
    // Note: ILUT doesn't work for multiple MPI ranks... you have to use AdditiveSchwarz.
    if (one_mpi_rank) {
      prec = rcp (new ILUT<MatrixType> (matrix));
    }
    else {
      prec = rcp (new AdditiveSchwarz<MatrixType, ILUT<MatrixType> > (matrix, overlap));
    }
  }
  else if (precTypeUpper == "RILUK") {
    prec = rcp (new RILUK<MatrixType> (matrix));
  }
  else if (precTypeUpper == "RELAXATION") {
    prec = rcp (new Relaxation<MatrixType> (matrix));
  }
  else if (precTypeUpper == "CHEBYSHEV") {
    prec = rcp (new Chebyshev<MatrixType> (matrix));
  }
  else if (precTypeUpper == "DIAGONAL") {
    prec = rcp (new Diagonal<MatrixType> (matrix));
  }
  else if (precTypeUpper == "SCHWARZ") {
    prec = rcp (new AdditiveSchwarz<MatrixType,ILUT<MatrixType> > (matrix, overlap));
  }
  else if (precTypeUpper == "KRYLOV") {
    prec = rcp (new Krylov<MatrixType, prec_base_type> (matrix));
  }
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
  else if (precTypeUpper == "SUPPORTGRAPH") {
    if (one_mpi_rank) {
      prec = rcp (new SupportGraph<MatrixType> (matrix));
    }
    else {
      prec = rcp (new AdditiveSchwarz<MatrixType, SupportGraph<MatrixType> > (matrix, overlap));
    }
  }
#endif
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "Ifpack2::Factory::create: "
      "Invalid preconditioner type \"" << precType << "\".");
  }
  return prec;
}


template<class MatrixType>
Teuchos::RCP<Preconditioner<typename MatrixType::scalar_type,
                                     typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
Factory::create (const std::string& precType,
                 const Teuchos::RCP<const MatrixType>& matrix)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;
  typedef Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type, node_type> prec_base_type;

  RCP<prec_base_type> prec;

  // precTypeUpper is the upper-case version of precType.
  std::string precTypeUpper (precType);
  if (precTypeUpper.size () > 0) {
    std::locale locale;
    for (size_t k = 0; k < precTypeUpper.size (); ++k) {
      precTypeUpper[k] = std::toupper<char> (precTypeUpper[k], locale);
    }
  }

  const bool one_mpi_rank = (matrix->getComm ()->getSize () == 1);

  if (precTypeUpper == "ILUT") {
    // Note: ILUT doesn't work for multiple MPI ranks... you have to use AdditiveSchwarz.
    if (one_mpi_rank) {
      prec = rcp (new ILUT<MatrixType> (matrix));
    }
    else {
      prec = rcp (new AdditiveSchwarz<MatrixType, ILUT<MatrixType> > (matrix));
    }
  }
  else if (precTypeUpper == "RILUK") {
    prec = rcp (new RILUK<MatrixType> (matrix));
  }
  else if (precTypeUpper == "RELAXATION") {
    prec = rcp (new Relaxation<MatrixType> (matrix));
  }
  else if (precTypeUpper == "CHEBYSHEV") {
    prec = rcp (new Chebyshev<MatrixType> (matrix));
  }
  else if (precTypeUpper == "DIAGONAL") {
    prec = rcp (new Diagonal<MatrixType> (matrix));
  }
  else if (precTypeUpper == "SCHWARZ") {
    prec = rcp (new AdditiveSchwarz<MatrixType, ILUT<MatrixType> > (matrix));
  }
  else if (precTypeUpper == "KRYLOV") {
    prec = rcp (new Krylov<MatrixType, prec_base_type> (matrix));
  }
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
  else if (precTypeUpper == "SUPPORTGRAPH") {
    if (one_mpi_rank) {
      prec = rcp (new SupportGraph<MatrixType> (matrix));
    }
    else {
      prec = rcp (new AdditiveSchwarz<MatrixType, SupportGraph<MatrixType> > (matrix));
    }
  }
#endif
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "Ifpack2::Factory::create: "
      "Invalid preconditioner type \"" << precType << "\".");
  }
  return prec;
}

template<class MatrixType, class M2>
Teuchos::RCP<Preconditioner<typename M2::scalar_type,
                                     typename M2::local_ordinal_type,
                                     typename M2::global_ordinal_type,
                                     typename M2::node_type> >
Factory::clone (const Teuchos::RCP<Preconditioner<typename MatrixType::scalar_type,
                                                           typename MatrixType::local_ordinal_type,
                                                           typename MatrixType::global_ordinal_type,
                                                           typename MatrixType::node_type> >& prec,
                const Teuchos::RCP<const M2>& matrix, const Teuchos::ParameterList& params)
{
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // FIXME (mfh 09 Nov 2013) The code below assumes that the old and
  // new scalar, local ordinal, and global ordinal types are the same.

  typedef typename M2::scalar_type scalar_type;
  typedef typename M2::local_ordinal_type local_ordinal_type;
  typedef typename M2::global_ordinal_type global_ordinal_type;
  typedef typename M2::node_type new_node_type;

  // FIXME (mfh 09 Nov 2013) The code below only knows how to clone
  // two different kinds of preconditioners!

  RCP<Preconditioner<scalar_type, local_ordinal_type,global_ordinal_type, new_node_type> > new_prec;
  RCP<Chebyshev<MatrixType> > chebyPrec;
  chebyPrec = rcp_dynamic_cast<Chebyshev<MatrixType> >(prec);
  if (chebyPrec != null){
        new_prec = chebyPrec->clone(matrix, params);
        return new_prec;
  }
  RCP<RILUK<MatrixType> > luPrec;
  luPrec = rcp_dynamic_cast<RILUK<MatrixType> >(prec);
  if (luPrec != null){
        new_prec = luPrec->clone(matrix);
        return new_prec;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error, "Ifpack2::Factory::clone: Not implemented for the "
    "current preconditioner type.  The only supported types thus far are "
    "Chebyshev and RILUK.");
}

} // namespace Ifpack2

#endif // IFPACK2_FACTORY_DEF_HPP
