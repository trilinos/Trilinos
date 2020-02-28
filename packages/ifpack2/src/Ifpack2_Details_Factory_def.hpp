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

#ifndef IFPACK2_DETAILS_FACTORY_DEF_HPP
#define IFPACK2_DETAILS_FACTORY_DEF_HPP

#include "Ifpack2_Factory.hpp"
#include "Ifpack2_Details_OneLevelFactory.hpp"
#include "Ifpack2_AdditiveSchwarz.hpp"
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
#  include "Ifpack2_SupportGraph.hpp"
#endif // defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)

namespace Ifpack2 {
namespace Details {

template<class SC, class LO, class GO, class NT>
Teuchos::RCP<typename Factory<SC, LO, GO, NT>::prec_type>
Factory<SC, LO, GO, NT>::
create (const std::string& precType,
        const Teuchos::RCP<const row_matrix_type>& matrix,
        const int overlap)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  RCP<prec_type> prec;

  // precTypeUpper is the upper-case version of precType.
  std::string precTypeUpper = canonicalize(precType);

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
    prec = rcp (new AdditiveSchwarz<row_matrix_type> (matrix, overlap));
  }
  else if (precTypeUpper == "KRYLOV") {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::invalid_argument, "The \"KRYLOV\" preconditioner option has "
       "been deprecated and removed.  If you want a Krylov solver, use the "
       "Belos package.");
  }
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
  else if (precTypeUpper == "SUPPORTGRAPH") {
    prec = rcp (new SupportGraph<row_matrix_type> (matrix));
  }
#endif
  else {
    try {
      Details::OneLevelFactory<row_matrix_type> factory;
      prec = factory.create (precType, matrix);
    } catch (std::invalid_argument&) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Ifpack2::Factory::create: "
        "Invalid preconditioner type \"" << precType << "\".");
    }
  }
  return prec;
}

template<class SC, class LO, class GO, class NT>
Teuchos::RCP<typename Factory<SC, LO, GO, NT>::prec_type>
Factory<SC, LO, GO, NT>::
create (const std::string& precType,
        const Teuchos::RCP<const row_matrix_type>& matrix)
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
    prec = rcp (new AdditiveSchwarz<row_matrix_type> (matrix));
  }
  else if (precTypeUpper == "KRYLOV") {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::invalid_argument, "The \"KRYLOV\" preconditioner option has "
       "been deprecated and removed.  If you want a Krylov solver, use the "
       "Belos package.");
  }
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
  else if (precTypeUpper == "SUPPORTGRAPH") {
    prec = rcp (new SupportGraph<row_matrix_type> (matrix));
  }
#endif
  else {
    bool success = false;
    std::ostringstream err;
    try {
      Details::OneLevelFactory<row_matrix_type> factory;
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


template<class SC, class LO, class GO, class NT>
bool
Factory<SC, LO, GO, NT>::
isSupported (const std::string& precType)
{
  // precTypeUpper is the upper-case version of precType.
  std::string precTypeUpper = canonicalize(precType);

  std::vector<std::string> supportedNames = {
    "SCHWARZ",
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
    "SUPPORTGRAPH"
#endif
  };
  auto it = std::find(std::begin(supportedNames), std::end(supportedNames), precTypeUpper);

  if (it != std::end(supportedNames)) {
    return true;
  } else {
    Details::OneLevelFactory<row_matrix_type> factory;
    return factory.isSupported (precType);
  }
}

} // namespace Details
} // namespace Ifpack2

#define IFPACK2_DETAILS_FACTORY_INSTANT(S, LO, GO, N) \
  template class Ifpack2::Details::Factory<S, LO, GO, N>;

#endif // IFPACK2_DETAILS_FACTORY_DEF_HPP
