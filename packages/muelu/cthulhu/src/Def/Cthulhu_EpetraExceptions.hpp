#ifndef CTHULHUEPETRA_EXCEPTIONS_HPP
#define CTHULHUEPETRA_EXCEPTIONS_HPP

#include <Cthulhu_Exceptions.hpp>

#define CATCH_EPETRA_EXCEPTION_AND_THROW_INVALID_ARG(sourceCode)        \
  try {                                                                 \
    sourceCode                                                          \
      } catch (int epetraErrCode) {                                     \
    TEST_FOR_EXCEPTION(true,std::invalid_argument, "Epetra threw exception: " << epetraErrCode); \
  }


#endif
