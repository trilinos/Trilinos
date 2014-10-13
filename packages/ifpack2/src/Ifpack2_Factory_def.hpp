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

#include "Ifpack2_Details_OneLevelFactory.hpp"
#include "Ifpack2_AdditiveSchwarz.hpp"
#include "Ifpack2_Krylov.hpp"
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
#  include "Ifpack2_SupportGraph.hpp"
#endif // defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
#include "Ifpack2_RILUK.hpp"
#include "Ifpack2_Factory.hpp"

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
  typedef Preconditioner<scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> prec_base_type;

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
  // Forestall "unused variable" warnings.
  (void) one_mpi_rank;

  if (precTypeUpper == "SCHWARZ") {
    // Discussion related to Bug 5987: The line of code below will
    // give AdditiveSchwarz a default subdomain solver by default.
    // However, you can change it later via setParameters() or
    // setInnerPreconditioner().  In the former case, AdditiveSchwarz
    // will not create the subdomain solver until you call
    // initialize(), so there is no performance loss in waiting after
    // calling AdditiveSchwarz's constructor before specifying the
    // subdomain solver's type.
    //
    // FIXME (mfh 14 Jan 2014) Use of "CUSTOM" in AdditiveSchwarz may
    // destroy information needed for fixing Bug 5987.  In particular,
    // the input ParameterList needs to keep its subdomain solver
    // info.  setInnerPreconditioner must _not_ destroy that info _if_
    // the Factory creates the AdditiveSchwarz instance.
    prec = rcp (new AdditiveSchwarz<MatrixType> (matrix, overlap));
  }
  else if (precTypeUpper == "KRYLOV") {
    prec = rcp (new Krylov<MatrixType> (matrix));
  }
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
  else if (precTypeUpper == "SUPPORTGRAPH") {
    prec = rcp (new SupportGraph<MatrixType> (matrix));
  }
#endif
  else {
    try {
      Details::OneLevelFactory<MatrixType> factory;
      prec = factory.create (precType, matrix);
    } catch (std::invalid_argument&) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Ifpack2::Factory::create: "
        "Invalid preconditioner type \"" << precType << "\".");
    }
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
  typedef Preconditioner<scalar_type,
                         local_ordinal_type,
                         global_ordinal_type,
                         node_type> prec_base_type;

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
  // Forestall "unused variable" warnings.
  (void) one_mpi_rank;

  if (precTypeUpper == "SCHWARZ") {
    // Discussion related to Bug 5987: The line of code below will
    // give AdditiveSchwarz a default subdomain solver by default.
    // However, you can change it later via setParameters() or
    // setInnerPreconditioner().  In the former case, AdditiveSchwarz
    // will not create the subdomain solver until you call
    // initialize(), so there is no performance loss in waiting after
    // calling AdditiveSchwarz's constructor before specifying the
    // subdomain solver's type.
    //
    // FIXME (mfh 14 Jan 2014) Use of "CUSTOM" in AdditiveSchwarz may
    // destroy information needed for fixing Bug 5987.  In particular,
    // the input ParameterList needs to keep its subdomain solver
    // info.  setInnerPreconditioner must _not_ destroy that info _if_
    // the Factory creates the AdditiveSchwarz instance.
    //
    // "CUSTOM" isn't necessary.  If Inverse_ is not null, then
    // AdditiveSchwarz's initialize() should just use the inner
    // preconditioner as it is.  If Inverse_ is null, then we assume
    // you want the default inner preconditioner.  You shouldn't have
    // called setInnerPreconditioner() with a null argument if that's
    // not what you meant!
    prec = rcp (new AdditiveSchwarz<MatrixType> (matrix));
  }
  else if (precTypeUpper == "KRYLOV") {
    prec = rcp (new Krylov<MatrixType> (matrix));
  }
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
  else if (precTypeUpper == "SUPPORTGRAPH") {
    prec = rcp (new SupportGraph<MatrixType> (matrix));
  }
#endif
  else {
    bool success = false;
    std::ostringstream err;
    try {
      Details::OneLevelFactory<MatrixType> factory;
      prec = factory.create (precType, matrix);
      success = true;
    } catch (std::invalid_argument& e) {
      err << "Ifpack2::Factory::create: Invalid preconditioner type \""
          << precType << "\".  More information for Ifpack2 developers: "
          << e.what ();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(! success, std::invalid_argument, err.str ());
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    prec.is_null (), std::logic_error, "Ifpack2::Factory::create: "
    "Return value is null right before return.  This should never happen.  "
    "Please report this bug to the Ifpack2 developers.");
  return prec;
}

template<class InputMatrixType, class OutputMatrixType>
Teuchos::RCP<Preconditioner<typename OutputMatrixType::scalar_type,
                            typename OutputMatrixType::local_ordinal_type,
                            typename OutputMatrixType::global_ordinal_type,
                            typename OutputMatrixType::node_type> >
Factory::clone (const Teuchos::RCP<Preconditioner<typename InputMatrixType::scalar_type,
                                                  typename InputMatrixType::local_ordinal_type,
                                                  typename InputMatrixType::global_ordinal_type,
                                                  typename InputMatrixType::node_type> >& prec,
                const Teuchos::RCP<const OutputMatrixType>& matrix,
                const Teuchos::ParameterList& params)
{
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // FIXME (mfh 09 Nov 2013) The code below assumes that the old and
  // new scalar, local ordinal, and global ordinal types are the same.

  typedef typename OutputMatrixType::scalar_type scalar_type;
  typedef typename OutputMatrixType::local_ordinal_type local_ordinal_type;
  typedef typename OutputMatrixType::global_ordinal_type global_ordinal_type;
  typedef typename OutputMatrixType::node_type new_node_type;
  typedef Preconditioner<scalar_type, local_ordinal_type,
                         global_ordinal_type, new_node_type> output_prec_type;

  // FIXME (mfh 09 Nov 2013) The code below only knows how to clone
  // two different kinds of preconditioners.  This is probably because
  // only two subclasses of Preconditioner implement a clone() method.

  RCP<output_prec_type> new_prec;
  RCP<Chebyshev<InputMatrixType> > chebyPrec;
  chebyPrec = rcp_dynamic_cast<Chebyshev<InputMatrixType> > (prec);
  if (chebyPrec != null) {
    new_prec = chebyPrec->clone (matrix, params);
    return new_prec;
  }
  RCP<RILUK<InputMatrixType> > luPrec;
  luPrec = rcp_dynamic_cast<RILUK<InputMatrixType> > (prec);
  if (luPrec != null) {
    new_prec = luPrec->clone (matrix);
    return new_prec;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error, "Ifpack2::Factory::clone: Not implemented for the "
    "current preconditioner type.  The only supported types thus far are "
    "Chebyshev and RILUK.");
}

} // namespace Ifpack2

#define IFPACK2_FACTORY_INSTANT(S,LO,GO,N)                              \
  template Teuchos::RCP<Ifpack2::Preconditioner<S, LO, GO, N> >         \
  Ifpack2::Factory::create<Tpetra::CrsMatrix< S, LO, GO, N> > (         \
    const std::string&,                                                 \
    const Teuchos::RCP<const Tpetra::CrsMatrix<S, LO, GO, N> >&);       \
  template Teuchos::RCP<Ifpack2::Preconditioner<S, LO, GO, N> >         \
  Ifpack2::Factory::create<Tpetra::CrsMatrix< S, LO, GO, N> > (         \
    const std::string&,                                                 \
    const Teuchos::RCP<const Tpetra::CrsMatrix<S, LO, GO, N> >&,        \
    const int);

#endif // IFPACK2_FACTORY_DEF_HPP
