// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestMain.cpp

\brief Ifpack2 Unit testing main program.

This file is the main for the unit test executable.

NOTE: This file should *not* be built and included as part of the Ifpack2
library.  It is instead to be directly included in the build files for
specific unit test suites.

*/


#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestRepository.hpp>
#include <Teuchos_GlobalMPISession.hpp>

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
#include <qd/dd_real.h>
#endif

int main( int argc, char* argv[] )
{
#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
  //If we're using extended-precision types, we have to run a QD-specific
  //function at the beginning and end of our main:
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
#endif

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
  fpu_fix_end(&old_cw);
#endif

  return ret;
}
