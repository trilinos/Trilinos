// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EPETRAEXCEPTIONS_HPP
#define XPETRA_EPETRAEXCEPTIONS_HPP

#include "Xpetra_ConfigDefs.hpp"

#ifndef HAVE_XPETRA_EPETRA
#error This file should be included only if HAVE_XPETRA_EPETRA is defined.
#endif

#include "Xpetra_Exceptions.hpp"

// This macro takes in argument a source code line.
// It catchs exceptions that could be throwed by 'sourceCode'
// If an exception is throw in any node, then all the node throws
// an std::invalid_argument exceptions.
#define IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG(sourceCode)                                    \
  {                                                                                                      \
    int localFailure = 0; /* 0 == success */                                                             \
    try {                                                                                                \
      sourceCode;                                                                                        \
    } catch (int /*epetraErrCode*/) {                                                                    \
      localFailure = 1; /* 1 == failure */                                                               \
    }                                                                                                    \
                                                                                                         \
    {                                                                                                    \
      int globalFailure = 0; /* 0 == success */                                                          \
      Teuchos::reduceAll<int>(*comm, Teuchos::REDUCE_SUM, localFailure, Teuchos::outArg(globalFailure)); \
      TEUCHOS_TEST_FOR_EXCEPTION(globalFailure != 0, std::invalid_argument, "Epetra threw exception");   \
    }                                                                                                    \
  }
#endif
