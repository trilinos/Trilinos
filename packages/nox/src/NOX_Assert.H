// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_ASSERT_H
#define NOX_ASSERT_H

#include "NOX_Common.H"
#include "Teuchos_Assert.hpp"

// Uses TEUCHOS_ASSERT but only in debug mode.  Check is disabled in
// released mode TEUCHOS_ASSERT is normally not disabled in optimized
// builds.
#ifdef HAVE_NOX_DEBUG
#define NOX_ASSERT(test) TEUCHOS_ASSERT(test)
#else
#define NOX_ASSERT(test)
#endif

#endif
