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
