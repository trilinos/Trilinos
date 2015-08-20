/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

#include "Ifpack2_Details_registerLinearSolverFactory.hpp"
#include "Ifpack2_Details_LinearSolverFactory.hpp"
#include "Ifpack2_ETIHelperMacros.h"

#if defined(HAVE_IFPACK2_EXPLICIT_INSTANTIATION)

// Define the Ifpack2 ETI macros and Tpetra typedefs that go along
// with them.
IFPACK2_ETI_MANGLING_TYPEDEFS()

#else

// Unfortunate hack for working around Bug 6380.  The Ifpack2 which
// defines mangled (for use in macros) typedefs for certain Scalar,
// GlobalOrdinal, and Node types, also has an empty definition if ETI
// is OFF.

#ifdef HAVE_TPETRA_COMPLEX_DOUBLE
typedef std::complex<double> std_complex0double0;
#endif // HAVE_TPETRA_COMPLEX_DOUBLE

#ifdef HAVE_TPETRA_COMPLEX_FLOAT
typedef std::complex<float> std_complex0float0;
#endif // HAVE_TPETRA_COMPLEX_FLOAT

#ifdef HAVE_TPETRA_INT_LONG_LONG
typedef long long longlong;
#endif // HAVE_TPETRA_INT_LONG_LONG

#ifdef HAVE_TPETRA_INT_UNSIGNED_LONG
typedef unsigned long unsignedlong;
#endif // HAVE_TPETRA_INT_UNSIGNED_LONG

#ifdef HAVE_TPETRA_@NT_MACRO_NAME@
typedef @NT@ @NT_MANGLED@;
#endif // HAVE_TPETRA_@NT_MACRO_NAME@

#endif // defined(HAVE_IFPACK2_EXPLICIT_INSTANTIATION)


// Local (to this file) macro that registers Ifpack2's
// LinearSolverFactory for the given four template parameters (Scalar
// = SC, LocalOrdinal = LO, GlobalOrdinal = GO, Node = NT).
#define LCL_CALL( SC, LO, GO, NT ) \
  ::Ifpack2::Details::LinearSolverFactory<SC, LO, GO, NT>::registerLinearSolverFactory ();

namespace Ifpack2 {
namespace Details {

void
registerLinearSolverFactory ()
{
#ifdef HAVE_IFPACK2_DEBUG
  std::cerr << " *** Ifpack2::Details::registerLinearSolverFactory *** " << std::endl;
#endif // HAVE_IFPACK2_DEBUG

  // Fill in the body of the function with all the type-specific
  // run-time registration functions.
  IFPACK2_INSTANTIATE_SLGN( LCL_CALL )
}

} // namespace Details
} // namespace Ifpack2


namespace { // (anonymous)

// \class RegisterLinearSolverFactory
// \brief Register Ifpack2's solver factory/ies with the central registry.
//
// \warning FOR EXPERT USE ONLY.
//
// Invoke this class' constructor to register Ifpack2's solver
// factory/ies with the central registry, for all template parameter
// combinations that Ifpack2 enabled.  You need not keep the instance
// of the class around; the constructor has a side effect if it
// returns.  (This is the C++ way of doing
// <tt>__attribute__((constructor))</tt>, without actually requiring
// the syntax extension.)
class RegisterLinearSolverFactory {
public:
  RegisterLinearSolverFactory () {
    Ifpack2::Details::registerLinearSolverFactory ();
  }
};

// Creating an instance of RegisterLinearSolverFactory invokes its
// constructor, which has the side effect of calling
// Ifpack2::Details::registerLinearSolverFactory().
RegisterLinearSolverFactory registerIt;

} // namespace (anonymous)



// static void
// #ifdef HAVE_TEUCHOS_CXX_ATTRIBUTE_CONSTRUCTOR
// __attribute__((constructor))
// #endif // HAVE_TEUCHOS_CXX_ATTRIBUTE_CONSTRUCTOR


