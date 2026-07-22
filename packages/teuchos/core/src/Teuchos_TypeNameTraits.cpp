// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
