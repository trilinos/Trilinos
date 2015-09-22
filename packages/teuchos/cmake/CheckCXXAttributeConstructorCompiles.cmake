INCLUDE(CheckCXXSourceCompiles)

# Check whether __attribute__((constructor)) compiles in C++.  This
# should work with GCC and compilers that claim GCC compatibility
# (Clang, Intel, recent PGI).  We don't try to test whether it
# actually works; this is purely a syntax check.
#
# TODO (mfh 19 Aug 2015) It would make sense actually to check that
# theConstructor() runs and produces the expected output.
FUNCTION(CHECK_CXX_ATTRIBUTE_CONSTRUCTOR_COMPILES VARNAME)
  SET(SOURCE
  "
#include <iostream>

__attribute__((constructor)) static void theConstructor () {
  std::cout << \"Hi!  I am the static constructor.\" << std::endl;
}

int main() {
  std::cout << \"Hi!  I am main.\" << std::endl;
  return 0;
}
  "
  )
  
  # This is a local variable name.  ${VARNAME} is the output variable.
  CHECK_CXX_SOURCE_COMPILES("${SOURCE}" HAVE_CXX_ATTRIBUTE_CONSTRUCTOR)

  IF(HAVE_CXX_ATTRIBUTE_CONSTRUCTOR)
    GLOBAL_SET(${VARNAME} TRUE)
    MESSAGE(STATUS "C++ compiler supports __attribute__((constructor)) syntax")
  ELSE()
    GLOBAL_SET(${VARNAME} FALSE)
    MESSAGE(STATUS "C++ compiler does NOT support __attribute__((constructor)) syntax")
  ENDIF()
ENDFUNCTION()
