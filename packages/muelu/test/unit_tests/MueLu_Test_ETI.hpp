// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TEST_ETI_H
#define MUELU_TEST_ETI_H

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include "MueLu_TestHelpers_Common.hpp"

// need this to have the ETI defined macros
#if defined(HAVE_MUELU_EXPLICIT_INSTANTIATION)
#include <MueLu_ExplicitInstantiation.hpp>
#endif

#include <TpetraCore_config.h>
#include <Tpetra_Details_KokkosTeuchosTimerInjection.hpp>

#include <KokkosKernels_config.h>
#include <KokkosKernels_Controls.hpp>

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
  using MueLu::Exceptions::RuntimeError;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // MPI initialization using Teuchos
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)));
#ifdef HAVE_MPI
  Teuchos::RCP<const Teuchos::MpiComm<int> > comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());
  if (comm->getSize() > 1) {
    out->setOutputToRootOnly(0);
  }
#endif

  // Tpetra nodes call Kokkos::execution_space::initialize if the execution
  // space is not initialized, but they don't call Kokkos::initialize.
  // Teuchos::GlobalMPISession captures its command-line arguments for later
  // use that Tpetra takes advantage of.
  //
  // We call Kokkos::initialize() after MPI so that MPI has the chance to bind
  // processes correctly before Kokkos touches things.
  Kokkos::initialize(argc, argv);

  // Create handles for cuBLAS and cuSPARSE. Otherwise they get
  // created on the first call to these libraries, and that can mess
  // up timings.
  KokkosKernels::Experimental::Controls controls;
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
  controls.getCublasHandle();
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  controls.getCusparseHandle();
#endif
  Kokkos::fence();
  bool success = true;
  bool verbose = true;
  try {
    // Parameters
    Teuchos::CommandLineProcessor clp(false);
    std::string node = "";
    clp.setOption("node", &node, "node type (serial | openmp | cuda | hip)");
    bool config = false;
    clp.setOption("config", "noconfig", &config, "display kokkos configuration");
#ifdef HAVE_TEUCHOS_STACKTRACE
    bool stacktrace = true;
    clp.setOption("stacktrace", "nostacktrace", &stacktrace, "display stacktrace");
#endif

    bool timedeepcopy = false;
    clp.setOption("timedeepcopy", "notimedeepcopy", &timedeepcopy, "instrument Kokkos::deep_copy() with Teuchos timers.  This can also be done with by setting the environment variable TPETRA_TIME_KOKKOS_DEEP_COPY=ON");
    bool timefence = false;
    clp.setOption("timefence", "notimefence", &timefence, "instrument Kokkos::fence() with Teuchos timers.  This can also be done with by setting the environment variable TPETRA_TIME_KOKKOS_FENCE=ON");
    Xpetra::Parameters xpetraParameters(clp);

    clp.recogniseAllOptions(false);
    switch (clp.parse(argc, argv, NULL)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: break;
    }

    if (timedeepcopy)
      Tpetra::Details::AddKokkosDeepCopyToTimeMonitor(true);
    if (timefence)
      Tpetra::Details::AddKokkosFenceToTimeMonitor(true);

#ifdef HAVE_TEUCHOS_STACKTRACE
    if (stacktrace)
      Teuchos::print_stack_on_segfault();
#endif

    auto lib = xpetraParameters.GetLib();
    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      // TAW: we might want to simplify the following logic block.
      //      In fact, there are examples/tests which only run with Tpetra
      //      We might need a feature that allows to run Epetra/Tpetra only
      //      We still need to make sure that the test compiles (i.e., we
      //      need some preprocessor flags/macros RUN_WITH_EPETRA and RUN_WITH_TPETRA
#if defined(HAVE_MUELU_INST_DOUBLE_INT_INT) || defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
      // Both Epetra and Tpetra (with double, int, int) enabled
      return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, int, Xpetra::EpetraNode>(clp, lib, argc, argv);
#else
      *out << "Skip running with Epetra since both Epetra and Tpetra are enabled but Tpetra is not instantiated on double, int, int." << std::endl;
#endif  // end Tpetra instantiated on double, int, int
#else
      throw RuntimeError("Epetra is not available");
#endif
    }

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
      auto inst = xpetraParameters.GetInstantiation();
#endif
      if (node == "") {
        typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType Node;

        if (config) {
          *out << "Node type: " << Node::execution_space::name() << std::endl;
          Node::execution_space().print_configuration(*out, true /*details*/);
        }

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_MUELU_INST_DOUBLE_INT_INT) || defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
        if (inst == Xpetra::DOUBLE_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT) || defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG)
        if (inst == Xpetra::DOUBLE_INT_LONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT) || defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        if (inst == Xpetra::DOUBLE_INT_LONGLONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_MUELU_INST_COMPLEX_INT_INT) || defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
        if (inst == Xpetra::COMPLEX_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<std::complex<double>, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_MUELU_INST_FLOAT_INT_INT) || defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_TPETRA_INST_INT_INT)
        if (inst == Xpetra::FLOAT_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<float, int, int, Node>(clp, lib, argc, argv);
#endif
        throw RuntimeError("Found no suitable Default instantiation");
#endif
      } else if (node == "serial") {
#ifdef KOKKOS_ENABLE_SERIAL
        typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode Node;

        if (config) {
          *out << "Node type: " << Node::execution_space::name() << std::endl;
          Kokkos::Serial().print_configuration(*out, true /*details*/);
        }

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        if (inst == Xpetra::DOUBLE_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGLONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_COMPLEX_INT_INT)
        if (inst == Xpetra::COMPLEX_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<std::complex<double>, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_FLOAT_INT_INT)
        if (inst == Xpetra::FLOAT_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<float, int, int, Node>(clp, lib, argc, argv);
#endif
        throw RuntimeError("Found no suitable Serial instantiation");
#endif
#else
        throw RuntimeError("Serial node type is disabled");
#endif
      } else if (node == "openmp") {
#ifdef KOKKOS_ENABLE_OPENMP
        typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode Node;

        if (config) {
          *out << "Node type: " << Node::execution_space::name() << std::endl;
          Kokkos::OpenMP().print_configuration(*out, true /*details*/);
          *out << "OpenMP Max Threads = " << omp_get_max_threads() << std::endl;
        }

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        if (inst == Xpetra::DOUBLE_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGLONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_COMPLEX_INT_INT)
        if (inst == Xpetra::COMPLEX_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<std::complex<double>, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_FLOAT_INT_INT)
        if (inst == Xpetra::FLOAT_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<float, int, int, Node>(clp, lib, argc, argv);
#endif
        throw RuntimeError("Found no suitable OpenMP instantiation");
#endif
#else
        throw RuntimeError("OpenMP node type is disabled");
#endif
      } else if (node == "cuda") {
#ifdef KOKKOS_ENABLE_CUDA
        typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode Node;

        if (config) {
          *out << "Node type: " << Node::execution_space::name() << std::endl;
          Kokkos::Cuda().print_configuration(*out, true /*details*/);
        }

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        if (inst == Xpetra::DOUBLE_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGLONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_COMPLEX_INT_INT)
        if (inst == Xpetra::COMPLEX_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<std::complex<double>, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_FLOAT_INT_INT)
        if (inst == Xpetra::FLOAT_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<float, int, int, Node>(clp, lib, argc, argv);
#endif
        throw RuntimeError("Found no suitable Cuda instantiation");
#endif
#else
        throw RuntimeError("CUDA node type is disabled");
#endif
      } else if (node == "hip") {
#ifdef KOKKOS_ENABLE_HIP
        typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode Node;

        if (config) {
          *out << "Node type: " << Node::execution_space::name() << std::endl;
          Kokkos::HIP().print_configuration(*out, true /*details*/);
        }

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_HIP) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        if (inst == Xpetra::DOUBLE_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_HIP) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_HIP) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGLONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_HIP) && defined(HAVE_MUELU_INST_COMPLEX_INT_INT)
        if (inst == Xpetra::COMPLEX_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<std::complex<double>, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_HIP) && defined(HAVE_MUELU_INST_FLOAT_INT_INT)
        if (inst == Xpetra::FLOAT_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<float, int, int, Node>(clp, lib, argc, argv);
#endif
        throw RuntimeError("Found no suitable HIP instantiation");
#endif
#else
        throw RuntimeError("HIP node type is disabled");
#endif
      } else if (node == "sycl") {
#ifdef KOKKOS_ENABLE_SYCL
        typedef Tpetra::KokkosCompat::KokkosSYCLWrapperNode Node;

        if (config) {
          *out << "Node type: " << Node::execution_space::name() << std::endl;
          Kokkos::Experimental::SYCL().print_configuration(*out, true /*details*/);
        }

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_SYCL) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        if (inst == Xpetra::DOUBLE_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_SYCL) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_SYCL) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        if (inst == Xpetra::DOUBLE_INT_LONGLONGINT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<double, int, long long, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_SYCL) && defined(HAVE_MUELU_INST_COMPLEX_INT_INT)
        if (inst == Xpetra::COMPLEX_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<std::complex<double>, int, int, Node>(clp, lib, argc, argv);
#endif
#if defined(HAVE_TPETRA_INST_SYCL) && defined(HAVE_MUELU_INST_FLOAT_INT_INT)
        if (inst == Xpetra::FLOAT_INT_INT)
          return MUELU_AUTOMATIC_TEST_ETI_NAME<float, int, int, Node>(clp, lib, argc, argv);
#endif
        throw RuntimeError("Found no suitable SYCL instantiation");
#endif
#else
        throw RuntimeError("SYCL node type is disabled");
#endif
      } else {
        throw RuntimeError("Unrecognized node type");
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  Kokkos::finalize();

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

// Undefine the function macro
#undef MUELU_AUTOMATIC_TEST_ETI_NAME

#endif
