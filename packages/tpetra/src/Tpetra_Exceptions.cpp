/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER
*/

#include "Tpetra_Exceptions.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

// FIXME (mfh 16 Dec 2012) See note below.
//#include "Tpetra_ETIHelperMacros.h"
#include "Tpetra_Exceptions_def.hpp"
#include "Teuchos_ConfigDefs.hpp"

namespace Tpetra {
namespace Details {

  // FIXME (mfh 16 Dec 2012) The TPETRA_INSTANTIATE_G macro (for
  // explicit instantiations over all desired GlobalOrdinal (== G)
  // types) does not exist yet.  Thus, rather than relying on the
  // Tpetra explicit instantation (ETI) macro system, I just roll my
  // own explicit instantiations for reasonable types.  I don't use
  // the explicit instantiation macros defined in
  // Tpetra_Exceptions_def.hpp because that would require using a
  // typedef for "long long" (since the type has a space in its name,
  // and the explicit instantiation macros don't like that).  I don't
  // want to clobber any existing typedef for "long long" for this
  // purpose.

  ////////////////////////////////////////////////////////////////////
  // Explicit instantiations for InvalidGlobalIndex
  ////////////////////////////////////////////////////////////////////

  //template class InvalidGlobalIndex<int>;
  template InvalidGlobalIndex<int>::InvalidGlobalIndex (const std::string& msg, const int globalIndex);
  template int InvalidGlobalIndex<int>::offendingIndex () const;

  //template class InvalidGlobalIndex<long>;
  template InvalidGlobalIndex<long>::InvalidGlobalIndex (const std::string& msg, const long globalIndex);
  template long InvalidGlobalIndex<long>::offendingIndex () const;

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  //template class InvalidGlobalIndex<long long>;
  template InvalidGlobalIndex<long long>::InvalidGlobalIndex (const std::string& msg, const long long globalIndex);
  template long long InvalidGlobalIndex<long long>::offendingIndex () const;
#endif // HAVE_TEUCHOS_LONG_LONG_INT

  ////////////////////////////////////////////////////////////////////
  // Explicit instantiations for InvalidGlobalRowIndex
  ////////////////////////////////////////////////////////////////////

  //template class InvalidGlobalRowIndex<int>;
  template InvalidGlobalRowIndex<int>::InvalidGlobalRowIndex (const std::string& msg, const int globalIndex);

  //template class InvalidGlobalRowIndex<long>;
  template InvalidGlobalRowIndex<long>::InvalidGlobalRowIndex (const std::string& msg, const long globalIndex);

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  //template class InvalidGlobalRowIndex<long long>;
  template InvalidGlobalRowIndex<long long>::InvalidGlobalRowIndex (const std::string& msg, const long long globalIndex);
#endif // HAVE_TEUCHOS_LONG_LONG_INT

} // namespace Details
} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
