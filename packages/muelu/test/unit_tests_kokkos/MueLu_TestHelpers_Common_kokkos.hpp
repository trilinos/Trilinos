// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
#ifndef MUELU_TEST_HELPERS_COMMON_KOKKOS_HPP
#define MUELU_TEST_HELPERS_COMMON_KOKKOS_HPP

#if 1
#define MUELU_DEBUGGER_MACRO                                         \
  out << "debug: scalar = " << typeid(GO).name() << std::endl;       \
  out << "debug: node   = " << typeid(NO).name() << std::endl;       \
  if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra) \
    out << "debug: linAlgebra = Epetra" << std::endl;                \
  else                                                               \
    out << "debug: linAlgebra = Tpetra" << std::endl;
#else
#define MUELU_DEBUGGER_MACRO
#endif

#define TYPE_EQUAL(TYPE1, TYPE2) \
  (typeid(TYPE1).name() == typeid(TYPE2).name())

#ifdef HAVE_MUELU_EPETRA

#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
// <double,int,int,EpetraNode>: run Epetra, do not run Tpetra
#define MUELU_TESTING_LIMIT_SCOPE(SC, GO, NO)                                                  \
  MUELU_DEBUGGER_MACRO                                                                         \
  if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra) {                         \
    if (!TYPE_EQUAL(SC, double)) {                                                             \
      out << "Skipping Epetra for SC != double" << std::endl;                                  \
      return;                                                                                  \
    }                                                                                          \
    if (!TYPE_EQUAL(GO, int)) {                                                                \
      out << "Skipping Epetra for GO != int" << std::endl;                                     \
      return;                                                                                  \
    }                                                                                          \
    if (!TYPE_EQUAL(NO, Xpetra::EpetraNode)) {                                                 \
      out << "Skipping Epetra for NO != EpetraNode" << std::endl;                              \
      return;                                                                                  \
    }                                                                                          \
  } else if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseTpetra) {                  \
    if (TYPE_EQUAL(SC, double) && TYPE_EQUAL(GO, int) && TYPE_EQUAL(NO, Xpetra::EpetraNode)) { \
      out << "Skipping Tpetra for <double, int, EpetraNode>" << std::endl;                     \
      return;                                                                                  \
    }                                                                                          \
  }

#define MUELU_TESTING_LIMIT_EPETRA_SCOPE_TPETRA_IS_DEFAULT(SC, GO, NO)                                           \
  {                                                                                                              \
    out << "Skipping as we cannot run Epetra and Tpetra for the same combination of template args" << std::endl; \
    return;                                                                                                      \
  }

#else

#define MUELU_TESTING_LIMIT_SCOPE(SC, GO, NO)                          \
  MUELU_DEBUGGER_MACRO                                                 \
  if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra) { \
    if (!TYPE_EQUAL(SC, double)) {                                     \
      out << "Skipping Epetra for SC != double" << std::endl;          \
      return;                                                          \
    }                                                                  \
    if (!TYPE_EQUAL(GO, int)) {                                        \
      out << "Skipping Epetra for GO != int" << std::endl;             \
      return;                                                          \
    }                                                                  \
    if (!TYPE_EQUAL(NO, Xpetra::EpetraNode)) {                         \
      out << "Skipping Epetra for NO != EpetraNode" << std::endl;      \
      return;                                                          \
    }                                                                  \
  }

// If linAlgebra==Tpetra, but the test also requires Epetra, this macro will cause the test
// to return early if SC!=double, GO!={int}, or NO!=Serial.
#define MUELU_TESTING_LIMIT_EPETRA_SCOPE_TPETRA_IS_DEFAULT(SC, GO, NO) \
  if (!TYPE_EQUAL(SC, double)) {                                       \
    out << "Skipping Epetra for SC != double" << std::endl;            \
    return;                                                            \
  }                                                                    \
  if (!TYPE_EQUAL(GO, int)) {                                          \
    out << "Skipping Epetra for GO != int" << std::endl;               \
    return;                                                            \
  }                                                                    \
  if (!TYPE_EQUAL(NO, Xpetra::EpetraNode)) {                           \
    out << "Skipping Epetra for NO != EpetraNode" << std::endl;        \
    return;                                                            \
  }
#endif

#else  // HAVE_MUELU_EPETRA

#define MUELU_TESTING_LIMIT_SCOPE(SC, GO, NO)
#define MUELU_TESTING_LIMIT_EPETRA_SCOPE_TPETRA_IS_DEFAULT(SC, GO, NO)

#endif  // HAVE_MUELU_EPETRA

// Macro to set MueLu's internal oh-so FancyOStream to be the same as the one used by Teuchos' unit testing framework.
// This prevents MueLu's output from intermingling with with the unit test pass/fail summary lines.
#define MUELU_TESTING_SET_OSTREAM \
  MueLu::VerboseObject::SetMueLuOStream(Teuchos::fancyOStream(out.getOStream()));

#define MUELU_TESTING_DO_NOT_TEST(lib, packagesNotEnabled)                                                                                                                                                                  \
  if (TestHelpers_kokkos::Parameters::getLib() == lib) {                                                                                                                                                                    \
    out << "Skipping test for " << (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra ? "Epetra" : "Tpetra") << " because some required packages are not enabled (" << packagesNotEnabled << ")." << std::endl; \
    return;                                                                                                                                                                                                                 \
  }

#endif  // ifndef MUELU_TEST_HELPERS_COMMON_KOKKOS_HPP
