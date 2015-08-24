//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
//@HEADER

#include "Belos_Details_registerLinearSolverFactory.hpp"
#include "Belos_Details_LinearSolverFactory.hpp"
#include "Teuchos_ConfigDefs.hpp" // for __attribute__((weak)) check

// FIXME (mfh 23 Aug 2015) Belos' main library is upstream from
// libraries where Belos' Epetra and Tpetra specializations live.
// That's why we need to use weak symbols here.  Otherwise, link
// errors result.  Fixing this requires the same run-time registration
// solution that Bug 6392 is all about.

#if defined(HAVE_BELOS_EPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
namespace Belos {
namespace Details {
namespace Epetra {
  extern void __attribute__((weak)) registerLinearSolverFactory ();
} // namespace Epetra
} // namespace Details
} // namespace Belos
#endif // defined(HAVE_BELOS_EPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)

#if defined(HAVE_BELOS_TPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
namespace Belos {
namespace Details {
namespace Tpetra {
  extern void __attribute__((weak)) registerLinearSolverFactory ();
} // namespace Tpetra
} // namespace Details
} // namespace Belos
#endif // defined(HAVE_BELOS_TPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)

//
// FIXME (mfh 23 Aug 2015) We should add Thyra as well.
//

namespace Belos {
namespace Details {

void
registerLinearSolverFactory ()
{
#if defined(HAVE_BELOS_EPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
  // It's a weak symbol, so it might be NULL.
  if (::Belos::Details::Tpetra::registerLinearSolverFactory != NULL) {
    ::Belos::Details::Tpetra::registerLinearSolverFactory ();
  }
#endif // defined(HAVE_BELOS_EPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)

#if defined(HAVE_BELOS_EPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
  // It's a weak symbol, so it might be NULL.
  if (::Belos::Details::Epetra::registerLinearSolverFactory != NULL) {
    ::Belos::Details::Epetra::registerLinearSolverFactory ();
  }
#endif // defined(HAVE_BELOS_EPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
}

} // namespace Details
} // namespace Belos


