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

#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_Assert.hpp"

// Define this if you want to force name demangling if supported
//#define HAVE_TEUCHOS_DEMANGLE

#if defined(HAVE_GCC_ABI_DEMANGLE) && defined(HAVE_TEUCHOS_DEMANGLE)
#  include <cxxabi.h>
#endif


std::string Teuchos::demangleName( const std::string &mangledName )
{
#if defined(HAVE_GCC_ABI_DEMANGLE) && defined(HAVE_TEUCHOS_DEMANGLE)
  int status;
  char* _demangledName = abi::__cxa_demangle (mangledName.c_str (), 0, 0, &status);
  if (status != 0 || 0 == _demangledName) {
#ifdef TEUCHOS_DEBUG
    // In a debug build, we check if demangling succeeded.
    std::string nullstr ("NULL");
    const char* demangle_output = _demangledName ? _demangledName : nullstr.c_str ();
    using std::endl;
#endif // TEUCHOS_DEBUG
    if (_demangledName != NULL) {
      // The C library standard requires that free() do the right
      // thing with NULL input, but it doesn't hurt to check.
      free (_demangledName);
    }
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error, name demangling with g++ has been enabled but the function "
      "abi::__cxa_demangle(" << mangledName << ") returned returnVal = " <<
      demangle_output <<" and status = " << status << "." << endl << "This "
      "suggests that name demangling may be broken on this platform." << endl
      << "You may prevent this exception from being thrown in the future by "
      "turning off name demangling for this build at configure time." << endl
      << "Do this by setting the CMake configuration option "
      "Teuchos_ENABLE_GCC_DEMANGLE to OFF." << endl << "Add the following to "
      "your list of CMake options:" << endl << endl <<
      "  -D Teuchos_ENABLE_GCC_DEMANGLE:BOOL=OFF" << endl);
#else // NOT TEUCHOS_DEBUG
    return (mangledName + "<demangle-failed>");
#endif // TEUCHOS_DEBUG
  }
  const std::string demangledName (_demangledName);
  free (_demangledName); // We have to free this before we return!
  return demangledName;
#else
  return mangledName;
#endif
}
