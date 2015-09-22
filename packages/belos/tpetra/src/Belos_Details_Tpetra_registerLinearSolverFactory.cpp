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

#include "BelosConfigDefs.hpp"

#ifdef HAVE_BELOS_TPETRA
#  include "BelosTpetraAdapter.hpp" // must go first
#  include "Belos_Details_registerLinearSolverFactory.hpp"
#  include "Belos_Details_LinearSolverFactory.hpp"
#  include "TpetraCore_ETIHelperMacros.h"

// Define Tpetra instantiation macros and typedefs that make the
// macros work.  The fix for Bug 6380 makes this work whether or not
// ETI is ON.  We use the Tpetra macros because Belos doesn't have
// its own macros.
TPETRA_ETI_MANGLING_TYPEDEFS()

// Macro that registers Belos's LinearSolverFactory for Tpetra
// objects, for the given four template parameters (Scalar = SC,
// LocalOrdinal = LO, GlobalOrdinal = GO, Node = NT).  The macro is
// local to this file.
//
// NOTE: This macro does NOT do explicit instantiation!  That's why I
// call it LCL_CALL and not LCL_INST.  We are just using the macros to
// invoke this class method over the set of enabled template
// parameters.
#define LCL_CALL( SC, LO, GO, NT ) \
  ::Belos::Details::LinearSolverFactory< ::Tpetra::MultiVector<SC, LO, GO, NT>, \
                                         ::Tpetra::Operator<SC, LO, GO, NT>, \
                                         SC, \
                                         typename ::Tpetra::MultiVector<SC, LO, GO, NT>::mag_type>::registerLinearSolverFactory ();

namespace Belos {
namespace Details {
namespace Tpetra {

void
registerLinearSolverFactory ()
{
  // Fill in the body of the function with all the type-specific
  // run-time registration functions, for registering Belos's
  // LinearSolverFactory with Tpetra objects.
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( LCL_CALL )
}

} // namespace Tpetra
} // namespace Details
} // namespace Belos

#endif // HAVE_BELOS_TPETRA

