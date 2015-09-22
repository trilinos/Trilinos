INCLUDE(CheckCXXSourceCompiles)

# This checks two things.  First, it checks whether #pragma weak
# $FUNCTION compiles in C++, when dealing with an extern void()
# function $FUNCTION in a namespace.  Second, it checks whether it is
# legit syntax to test the resulting function pointer and call it if
# it's not NULL.
#
# The first thing might be a bit more standard than the
# "__attribute__((weak))" syntax of GCC and compatible compilers.  The
# second might not work on Mac, but should work on Linux.  That's OK,
# because we only need to use this technique with static libraries
# (see Bug 6392), and Macs don't like static libraries.
FUNCTION(CHECK_CXX_PRAGMA_WEAK_COMPILES VARNAME)
  SET(SOURCE
  "
#include <iostream>

namespace A {
// theFunction never gets defined, because we
// don't link with a library that defines it.
// That's OK, because it's weak linkage.
#pragma weak theFunction
extern void theFunction ();
}

int main() {
  std::cout << \"Hi!  I am main.\" << std::endl;
  if (A::theFunction != NULL) {
    // Should never be called, since we don't link
    // with a library that defines A::theFunction.
    A::theFunction ();
  }
  return 0;
}
  "
  )
  
  # This is a local variable name.  ${VARNAME} is the output variable.
  CHECK_CXX_SOURCE_COMPILES("${SOURCE}" HAVE_CXX_PRAGMA_WEAK)

  IF(HAVE_CXX_PRAGMA_WEAK)
    GLOBAL_SET(${VARNAME} TRUE)
    MESSAGE(STATUS "C++ compiler supports #pragma weak syntax and testing weak functions")
  ELSE()
    GLOBAL_SET(${VARNAME} FALSE)
    MESSAGE(STATUS "C++ compiler does NOT support #pragma weak syntax and testing weak functions")
  ENDIF()
ENDFUNCTION()
