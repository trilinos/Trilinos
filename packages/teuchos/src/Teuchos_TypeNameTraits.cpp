// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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

#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_TestForException.hpp"

// Define this if you want to fource name demangling if supported
//#define HAVE_TEUCHOS_DEMANGLE

#if defined(HAVE_GCC_ABI_DEMANGLE) && defined(HAVE_TEUCHOS_DEMANGLE)
#  include <cxxabi.h>
#endif


std::string Teuchos::demangleName( const std::string &mangledName )
{
#if defined(HAVE_GCC_ABI_DEMANGLE) && defined(HAVE_TEUCHOS_DEMANGLE)
  int status;
  char *_demangledName = abi::__cxa_demangle(mangledName.c_str(), 0, 0, &status);
  if (status != 0 || 0==_demangledName) {
#ifdef TEUCHOS_DEBUG
    std::string nullstr("NULL");
    const char* demangle_output = _demangledName ? _demangledName : nullstr.c_str();
    TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error, name demangling with g++ has been enabled but the function "
      "abi::__cxa_demangle("<<mangledName<<") returned returnVal = "<<demangle_output
      <<" and status = "<<status<<".  Name demangling for this build"
      "can be turned off by using --disable-teuchos-demangle at configure time." );
#endif
    if (_demangledName != NULL)
      free(_demangledName);
    return ( mangledName + "<demangle-failed>" );
  }
  const std::string demangledName(_demangledName);
  free(_demangledName); // We have to free this before we return!
  return demangledName;
#else
  return mangledName;
#endif
}
