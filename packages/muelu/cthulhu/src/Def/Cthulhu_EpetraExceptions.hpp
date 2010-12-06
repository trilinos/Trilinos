#ifndef CTHULHU_EPETRA_EXCEPTIONS_HPP
#define CTHULHU_EPETRA_EXCEPTIONS_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_EPETRA
#error This file should be included only if HAVE_CTHULHU_EPETRA is defined.
#endif

#include "Cthulhu_Exceptions.hpp"

// This macro takes in argument a source code line.
// It catchs exceptions that could be throwed by 'sourceCode'
// If an exception is throw in any node, then all the node throws
// an std::invalid_argument exceptions.
#define IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG(sourceCode)   \
  {                                                                     \
    int localSuccess = 1; /* 1 == true */                               \
    try {                                                               \
      sourceCode;                                                       \
    }                                                                   \
    catch (int epetraErrCode) {                                         \
      localSuccess = 0; /* 0 == false */                                \
    }                                                                   \
                                                                        \
    {                                                                   \
      int globalSuccess = 0; /* 0 == success */                         \
      Teuchos::reduceAll<int>(*comm, Teuchos::REDUCE_SUM, localSuccess, Teuchos::outArg(globalSuccess)); \
      TEST_FOR_EXCEPTION(globalSuccess==0, std::invalid_argument, "Epetra threw exception"); \
    }                                                                   \
  }
#endif

