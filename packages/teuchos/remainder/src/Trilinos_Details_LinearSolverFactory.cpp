// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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


