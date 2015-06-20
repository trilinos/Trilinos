// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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
#define IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG(sourceCode)   \
  {                                                                     \
    int localFailure = 0; /* 0 == success */                            \
    try {                                                               \
      sourceCode;                                                       \
    }                                                                   \
    catch (int /*epetraErrCode*/) {                                     \
      localFailure = 1; /* 1 == failure */                              \
    }                                                                   \
                                                                        \
    {                                                                   \
      int globalFailure = 0; /* 0 == success */                         \
      Teuchos::reduceAll<int>(*comm, Teuchos::REDUCE_SUM, localFailure, Teuchos::outArg(globalFailure)); \
      TEUCHOS_TEST_FOR_EXCEPTION(globalFailure != 0, std::invalid_argument, "Epetra threw exception"); \
    }                                                                   \
  }
#endif

