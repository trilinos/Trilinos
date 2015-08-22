// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#include "Trilinos_Details_LinearSolverFactory.hpp"
#include <set>

namespace Trilinos {
namespace Details {
namespace Impl {

namespace { // (anonymous)

// All package names that registerLinearSolverFactory has seen,
// for any combination of template parameters MV and OP.
static std::set<std::string>* packageNames_ = NULL;

// atexit() hook for freeing packageNames_.
void freePackageNames ()
{
  if (packageNames_ != NULL) {
    delete packageNames_;
    packageNames_ = NULL;
  }
}

void createPackageNames ()
{
  if (packageNames_ == NULL) {
    packageNames_ = new std::set<std::string> ();
    // It _is_ possible for atexit() to fail (e.g., because it ran
    // out of memory for storing callbacks).  We could throw an
    // exception here in that case, but I think it's better just
    // to let the minor memory leak happen.
    (void) atexit (freePackageNames);
  }
}

} // namespace (anonymous)

bool rememberRegisteredSomeLinearSolverFactory (const std::string& packageName)
{
  createPackageNames ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (packageNames_ == NULL, std::logic_error, "Trilinos::Details::"
     "Impl::rememberRegisteredSomeLinearSolverFactory: "
     "Should never get here!  packageNames_ is NULL.");

  std::pair<std::set<std::string>::iterator, bool> ret =
    packageNames_->insert (packageName);
  // ret.second is true if the item was NOT in the set before.
  return ! ret.second;
}

bool registeredSomeLinearSolverFactory (const std::string& packageName)
{
  createPackageNames ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (packageNames_ == NULL, std::logic_error, "Trilinos::Details::"
     "Impl::rememberRegisteredSomeLinearSolverFactory: "
     "Should never get here!  packageNames_ is NULL.");

  std::set<std::string>::const_iterator it = packageNames_->find (packageName);
  return it != packageNames_->end ();
}

bool haveLinearSolverFactoryRunTimeRegistration ()
{
#if defined(TRILINOS_HAVE_LINEAR_SOLVER_FACTORY_REGISTRATION)
  return true;
#else // NOT defined(TRILINOS_HAVE_LINEAR_SOLVER_FACTORY_REGISTRATION)
  return false;
#endif // defined(TRILINOS_HAVE_LINEAR_SOLVER_FACTORY_REGISTRATION)
}

} // namespace Impl
} // namespace Details
} // namespace Trilinos


