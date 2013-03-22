// @HEADER
// ***********************************************************************
// 
//      Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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
