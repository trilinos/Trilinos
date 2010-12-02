#ifndef CTHULHU_EPETRA_EXCEPTIONS_HPP
#define CTHULHU_EPETRA_EXCEPTIONS_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_EPETRA
#error This file should be included only if HAVE_CTHULHU_EPETRA is defined.
#endif

#include "Cthulhu_Exceptions.hpp"

#define CATCH_EPETRA_EXCEPTION_AND_THROW_INVALID_ARG(sourceCode)        \
  try {                                                                 \
    sourceCode                                                          \
      } catch (int epetraErrCode) {                                     \
    TEST_FOR_EXCEPTION(true,std::invalid_argument, "Epetra threw exception: " << epetraErrCode); \
  }

#endif

