// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TEST_HELPERS_COMMON_KOKKOS_HPP
#define MUELU_TEST_HELPERS_COMMON_KOKKOS_HPP

#if 1
#define MUELU_DEBUGGER_MACRO                                   \
  out << "debug: scalar = " << typeid(GO).name() << std::endl; \
  out << "debug: node   = " << typeid(NO).name() << std::endl; \
  out << "debug: linAlgebra = Tpetra" << std::endl;
#else
#define MUELU_DEBUGGER_MACRO
#endif

#define TYPE_EQUAL(TYPE1, TYPE2) \
  (typeid(TYPE1).name() == typeid(TYPE2).name())

#define MUELU_TESTING_LIMIT_SCOPE(SC, GO, NO)
#define MUELU_TESTING_LIMIT_EPETRA_SCOPE_TPETRA_IS_DEFAULT(SC, GO, NO)

// Macro to set MueLu's internal oh-so FancyOStream to be the same as the one used by Teuchos' unit testing framework.
// This prevents MueLu's output from intermingling with with the unit test pass/fail summary lines.
#define MUELU_TESTING_SET_OSTREAM \
  MueLu::VerboseObject::SetMueLuOStream(Teuchos::fancyOStream(out.getOStream()));

#define MUELU_TESTING_DO_NOT_TEST(lib, packagesNotEnabled)                                                 \
  if (TestHelpers_kokkos::Parameters::getLib() == lib) {                                                   \
    out << "Skipping test for "                                                                            \
        << "Tpetra"                                                                                        \
        << " because some required packages are not enabled (" << packagesNotEnabled << ")." << std::endl; \
    return;                                                                                                \
  }

#endif  // ifndef MUELU_TEST_HELPERS_COMMON_KOKKOS_HPP
