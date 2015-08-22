INCLUDE(CheckCXXSourceCompiles)

# This checks two things.  First, it checks whether
# __attribute__((weak)) compiles in C++, when dealing with an extern
# void() function in a namespace.  Second, it checks whether it is
# legit syntax to test the resulting function pointer and call it if
# it's not NULL.
#
# The first thing should work with GCC and compilers that claim GCC
# compatibility (Clang, Intel, recent PGI).  The second doesn't work
# on Mac, but does work on Linux.  That's OK, because we only need to
# use this technique with static libraries (see Bug 6392), and Macs
# don't like static libraries.
#
# #pragma weak $FUNCTION (see CheckCXXPragmaWeakCompiles.cmake in this
# #directory) seems to be implemented by more compilers.  GCC claims
# #to implement it "[f]or compatibility with SVR4":
#
# https://gcc.gnu.org/onlinedocs/gcc/Weak-Pragmas.html
# 
# Other compilers implement this as well, probably for the same
# reason.  For example (all links tested 20 Aug 2015):
#
#   - Intel's C++ compiler:
#     https://software.intel.com/en-us/node/524560
#   - Sun Studio 12: C++ User's Guide:
#     http://docs.oracle.com/cd/E19205-01/819-5267/bkbkr/index.html
#   - PGI: I can't find anything in the User's Guide or Reference, but
#     the following bug fix report suggests that #pragma weak works
#     (otherwise, how did PGI 2014 fix a problem with it?):
#     https://www.pgroup.com/support/release_tprs_2014.htm
#   - "IBM XL C/C++ for Linux, V11.1 Compiler Reference, Version 11.1":
#     the compiler does support #pragma weak, BUT unfortunately requires
#     the mangled C++ name:
#     http://www-01.ibm.com/support/docview.wss?uid=swg27018970&aid=1
#     That means it won't pass this test, though.  However, IBM XL
#     C/C++ for AIX 13.1.2 supports __attribute__((weak)):
#     http://www-01.ibm.com/support/knowledgecenter/SSGH3R_13.1.2/com.ibm.xlcpp131.aix.doc/language_ref/fn_attrib_weak.html
FUNCTION(CHECK_CXX_ATTRIBUTE_WEAK_COMPILES VARNAME)
  SET(SOURCE
  "
#include <iostream>

namespace A {
// theFunction never gets defined, because we
// don't link with a library that defines it.
// That's OK, because it's weak linkage.
extern void __attribute__((weak)) theFunction ();
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
  CHECK_CXX_SOURCE_COMPILES("${SOURCE}" HAVE_CXX_ATTRIBUTE_WEAK)

  IF(HAVE_CXX_ATTRIBUTE_WEAK)
    GLOBAL_SET(${VARNAME} TRUE)
    MESSAGE(STATUS "C++ compiler supports __attribute__((weak)) syntax and testing weak functions")
  ELSE()
    GLOBAL_SET(${VARNAME} FALSE)
    MESSAGE(STATUS "C++ compiler does NOT support __attribute__((weak)) syntax and testing weak functions")
  ENDIF()
ENDFUNCTION()
