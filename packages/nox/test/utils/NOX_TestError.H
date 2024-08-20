// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* The purpose of this file is to allow for
informative error reporting that will not abort a test program after a single
error the way that an assert command does. */

/* The macro takes two arguments.  The first is the error code to be examined.  The second is an int that can be viewed either as an error count, or as a bool that simply indicates if any errors have occurred as of the time that the macro was run (zero -> no prior errors, non-zero -> prior errors). */

/* This is a rip off of the EPETRA_TEST_ERR function written by Mike Heroux for Epetra regression testing */

#ifndef _NOX_TEST_ERR_H_
#define _NOX_TEST_ERR_H_
#include "NOX_Common.H"

// This function is to be used when first identifying an error.
#define NOX_TEST_ERR(a,b) { { int nox_testing_err = a; \
  if (nox_testing_err != 0) {\
    std::cerr << "Non zero error code " << nox_testing_err << \
       ", file: " << __FILE__ << ", line: " << __LINE__ << std::endl;\
 b+=1;\
  }\
  }\
}

#endif /*_NOX_TEST_ERR_H_ */
