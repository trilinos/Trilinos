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
#ifndef MUELU_TEST_ETI_H
#define MUELU_TEST_ETI_H

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include "MueLu_TestHelpers_Common.hpp"

#ifndef MUELU_AUTOMATIC_TEST_ETI_NAME
#error "The macro MUELU_AUTOMATIC_TEST_ETI_NAME was not defined"
#endif


// This function serves to manage the ETI (or not) for tests that are *not* unit tests.  This
// allows us to isolate this particular bit of functionality, rather than have it cut-and-paste
// duplicated over half of the testing tree.  See documentation at
// https://github.com/muelu/Developer-Guide/wiki/Writing-a-non-unit-test.
// TAW: 7/24/17 Please be aware of the fact that this routine only tests one enabled configuration
//              either Epetra or Tpetra (depending on the user choice with --linAlgebra)
//              and only one instantiation (depending on the choice of --instantiation).
//              It does NOT test all enabled configurations as the unit tests would do.
bool Automatic_Test_ETI(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using MueLu::Exceptions::RuntimeError;

  // MPI initialization using Teuchos
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  // Tpetra nodes call Kokkos::execution_space::initialize if the execution
  // space is not initialized, but they don't call Kokkos::initialize.
  // Teuchos::GlobalMPISession captures its command-line arguments for later
  // use that Tpetra takes advantage of.
  //
  // We call Kokkos::initialize() after MPI so that MPI has the chance to bind
  // processes correctly before Kokkos touches things.
#ifdef HAVE_MUELU_KOKKOSCORE
  Kokkos::initialize(argc, argv);
#endif

  bool success = true;
  bool verbose = true;
  try {
    // Parameters
    Teuchos::CommandLineProcessor clp(false);
    std::string node   = "";    clp.setOption("node",               &node,   "node type (serial | openmp | cuda)");
    bool        config = false; clp.setOption("config", "noconfig", &config, "display kokkos configuration");
    Xpetra::Parameters xpetraParameters(clp);

    clp.recogniseAllOptions(false);
    switch (clp.parse(argc, argv, NULL)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR:                return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:         break;
    }

    auto lib = xpetraParameters.GetLib();
    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      // TAW: we might want to simplify the following logic block.
      //      In fact, there are examples/tests which only run with Tpetra
      //      We might need a feature that allows to run Epetra/Tpetra only
      //      We still need to make sure that the test compiles (i.e., we
      //      need some preprocessor flags/macros RUN_WITH_EPETRA and RUN_WITH_TPETRA
#  ifdef HAVE_MUELU_TPETRA
#    if defined(HAVE_MUELU_INST_DOUBLE_INT_INT) || defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
      // Both Epetra and Tpetra (with double, int, int) enabled
      return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,int,Xpetra::EpetraNode>(clp, lib, argc, argv);
#    else
      std::cout << "Skip running with Epetra since both Epetra and Tpetra are enabled but Tpetra is not instantiated on double, int, int." << std::endl;
#    endif // end Tpetra instantiated on double, int, int
#  else
      // only Epetra enabled. No Tpetra instantiation possible
      return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,int,Xpetra::EpetraNode>(clp, lib, argc, argv);
#  endif // HAVE_MUELU_TPETRA
#else
      throw RuntimeError("Epetra is not available");
#endif
    }

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      auto inst = xpetraParameters.GetInstantiation();
      if (node == "") {
        typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long,Node>(clp, lib, argc, argv);
#else
#  if defined(HAVE_MUELU_INST_DOUBLE_INT_INT) || defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
        if (inst == Xpetra::DOUBLE_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,int,Node> (clp,  lib, argc, argv);
#  endif
#  if defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT) || defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG)
        if (inst == Xpetra::DOUBLE_INT_LONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long,Node>(clp,  lib, argc, argv);
#  endif
#  if defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT) || defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        if (inst == Xpetra::DOUBLE_INT_LONGLONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long long,Node>(clp,  lib, argc, argv);
#  endif
#  if defined(HAVE_MUELU_INST_COMPLEX_INT_INT) || defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
        if (inst == Xpetra::COMPLEX_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<std::complex<double>,int,int,Node>(clp,  lib, argc, argv);
#  endif
#  if defined(HAVE_MUELU_INST_FLOAT_INT_INT) || defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_TPETRA_INST_INT_INT)
        if (inst == Xpetra::FLOAT_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<float,int,int,Node>(clp,  lib, argc, argv);
#  endif
        throw RuntimeError("Found no suitable instantiation");
#endif
      } else if (node == "serial") {
#ifdef KOKKOS_ENABLE_SERIAL
        typedef Kokkos::Compat::KokkosSerialWrapperNode Node;

        if (config)
          Kokkos::Serial::print_configuration(std::cout, true/*details*/);

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long,Node>(clp,  lib, argc, argv);
#  else
#    if defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        if (inst == Xpetra::DOUBLE_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,int,Node> (clp,  lib, argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long,Node>(clp,  lib, argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGLONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long long,Node>(clp,  lib, argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_COMPLEX_INT_INT)
        if (inst == Xpetra::COMPLEX_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<std::complex<double>,int,int,Node>(clp,  lib, argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_FLOAT_INT_INT)
        if (inst == Xpetra::FLOAT_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<float,int,int,Node>(clp,  lib, argc, argv);
#    endif
        throw RuntimeError("Found no suitable instantiation");
#  endif
#else
        throw RuntimeError("Serial node type is disabled");
#endif
      } else if (node == "openmp") {
#ifdef KOKKOS_ENABLE_OPENMP
        typedef Kokkos::Compat::KokkosOpenMPWrapperNode Node;

        if (config) {
          Kokkos::OpenMP::print_configuration(std::cout, true/*details*/);
          std::cout << "OpenMP Max Threads = " << omp_get_max_threads() << std::endl;
        }

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long,Node>(clp,  lib, argc, argv);
#  else
#    if defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        if (inst == Xpetra::DOUBLE_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,int,Node> (clp, lib,  argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long,Node>(clp,  lib, argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGLONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long long,Node>(clp,  lib, argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_COMPLEX_INT_INT)
        if (inst == Xpetra::COMPLEX_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<std::complex<double>,int,int,Node>(clp,  lib, argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_FLOAT_INT_INT)
        if (inst == Xpetra::FLOAT_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<float,int,int,Node>(clp,  lib, argc, argv);
#    endif
        throw RuntimeError("Found no suitable instantiation");
#  endif
#else
        throw RuntimeError("OpenMP node type is disabled");
#endif
      } else if (node == "cuda") {
#ifdef KOKKOS_ENABLE_CUDA
        typedef Kokkos::Compat::KokkosCudaWrapperNode Node;

        if (config)
          Kokkos::Cuda::print_configuration(std::cout, true/*details*/);

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long,Node>(clp, lib, argc, argv);
#  else
#    if defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        if (inst == Xpetra::DOUBLE_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,int,Node> (clp, lib, argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long,Node>(clp, lib, argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGLONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double,int,long long,Node>(clp, lib, argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_COMPLEX_INT_INT)
        if (inst == Xpetra::COMPLEX_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<std::complex<double>,int,int,Node>(clp,  lib, argc, argv);
#    endif
#    if defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_FLOAT_INT_INT)
        if (inst == Xpetra::FLOAT_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<float,int,int,Node>(clp,  lib, argc, argv);
#    endif
        throw RuntimeError("Found no suitable instantiation");
#  endif
#else
        throw RuntimeError("CUDA node type is disabled");
#endif
      } else {
        throw RuntimeError("Unrecognized node type");
      }
#else
      throw RuntimeError("Tpetra is not available");
#endif
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef HAVE_MUELU_KOKKOSCORE
  Kokkos::finalize();
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

// Undefine the function macro
#undef MUELU_AUTOMATIC_TEST_ETI_NAME

#endif
