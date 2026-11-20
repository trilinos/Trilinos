// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER

// Tests of KokkosParser initializing or joining Kokkos (see the KokkosParser documentation). Kokkos can be
// initialized only once per process, so each mode, given as the first argument, runs in its own process:
//   args               KokkosParser(argc, argv) initializes Kokkos; later constructors with arguments throw
//   no_args            KokkosParser(), with Kokkos not initialized, warns and initializes from the environment
//   external           the program initializes Kokkos: KokkosParser() joins, KokkosParser(argc, argv) throws
//   embedded           (Python) pycompadre.KokkosParser() in an embedded interpreter, Kokkos not initialized
//   embedded_external  (Python) pycompadre.KokkosParser() in an embedded interpreter, after the program
//                      initialized Kokkos

#include <Compadre_KokkosParser.hpp>

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#ifdef COMPADRE_TEST_EMBEDDED_PYTHON
#include <Python.h>
#endif

using namespace Compadre;

namespace {

int failures = 0;

void check(const bool condition, const std::string& what) {
    if (!condition) {
        std::cout << "FAILED: " << what << std::endl;
        failures++;
    }
}

template <typename F>
bool throwsRuntimeError(F&& f) {
    try {
        f();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

// captures what is written to std::cerr while it is alive
struct CaptureCerr {
    std::stringstream captured;
    std::streambuf* original;
    CaptureCerr() : original(std::cerr.rdbuf(captured.rdbuf())) {}
    ~CaptureCerr() { std::cerr.rdbuf(original); }
    bool contains(const std::string& text) const { return captured.str().find(text) != std::string::npos; }
};

// part of the warning KokkosParser() gives when it initializes Kokkos without the command line
const std::string ignored_arguments_warning = "command-line arguments are ignored";

void checkThreads(const int expected) {
#ifdef KOKKOS_ENABLE_OPENMP
    if (std::is_same_v<Kokkos::DefaultHostExecutionSpace, Kokkos::OpenMP>) {
        check(Kokkos::num_threads() == expected, "Kokkos uses the thread count given on the command line");
    }
#endif
}

void testArgs(int argc, char* argv[]) {
    const std::vector<std::string> before(argv, argv + argc);
    {
        KokkosParser kp(argc, argv);
        check(Kokkos::is_initialized(), "KokkosParser(argc, argv) initializes Kokkos");
        checkThreads(2);
        check(std::vector<std::string>(argv, argv + argc) == before,
                "KokkosParser(argc, argv) does not modify the caller's argv");
        check(throwsRuntimeError([&]() { KokkosParser again(argc, argv); }),
                "KokkosParser(argc, argv) throws once Kokkos is initialized");
        check(throwsRuntimeError([&]() { KokkosParser again(before); }),
                "KokkosParser(args) throws once Kokkos is initialized");
        check(throwsRuntimeError([&]() { KokkosParser again{KokkosInitArguments()}; }),
                "KokkosParser(settings) throws once Kokkos is initialized");
        CaptureCerr captured;
        {
            KokkosParser joined;
            check(Kokkos::is_initialized(), "KokkosParser() joins the running Kokkos");
        }
        check(!captured.contains(ignored_arguments_warning), "KokkosParser() joins without a warning");
        check(Kokkos::is_initialized(), "a KokkosParser that joined does not finalize Kokkos");
    }
    check(Kokkos::is_finalized(), "Kokkos is finalized with the last KokkosParser");
    check(throwsRuntimeError([]() { KokkosParser after; }), "KokkosParser() throws once Kokkos is finalized");
    check(throwsRuntimeError([&]() { KokkosParser after(before); }),
            "KokkosParser(args) throws once Kokkos is finalized");
}

void testNoArgs() {
    {
        CaptureCerr captured;
        KokkosParser kp;
        check(Kokkos::is_initialized(), "KokkosParser() initializes Kokkos when it is not initialized");
        check(captured.contains(ignored_arguments_warning),
                "KokkosParser() warns that command-line arguments are ignored");
    }
    check(Kokkos::is_finalized(), "Kokkos is finalized with the last KokkosParser");
}

void testExternal(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    check(throwsRuntimeError([&]() { KokkosParser kp(argc, argv); }),
            "KokkosParser(argc, argv) throws when the program initialized Kokkos");
    {
        CaptureCerr captured;
        {
            KokkosParser kp;
        }
        check(!captured.contains(ignored_arguments_warning),
                "KokkosParser() joins the program's Kokkos without a warning");
    }
    check(Kokkos::is_initialized(), "KokkosParser() does not finalize the program's Kokkos");
    Kokkos::finalize();
}

#ifdef COMPADRE_TEST_EMBEDDED_PYTHON
void testEmbedded() {
    CaptureCerr captured;
    Py_Initialize();
    check(PyRun_SimpleString(
            "import sys, pycompadre\n"
            "assert sys.orig_argv == [], sys.orig_argv\n"
            "kp = pycompadre.KokkosParser()\n") == 0, "pycompadre.KokkosParser() runs embedded");
    check(Kokkos::is_initialized(), "embedded, pycompadre.KokkosParser() initializes Kokkos when it is not");
    check(captured.contains(ignored_arguments_warning),
            "embedded, pycompadre.KokkosParser() warns that command-line arguments are ignored");
    check(PyRun_SimpleString("del kp\n") == 0, "pycompadre.KokkosParser is released");
    check(Kokkos::is_finalized(), "Kokkos is finalized with the last pycompadre.KokkosParser");
    Py_Finalize();
}

void testEmbeddedExternal(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        CaptureCerr captured;
        Py_Initialize();
        check(PyRun_SimpleString(
                "import sys, pycompadre\n"
                "kp = pycompadre.KokkosParser()\n"
                "kp_argv = pycompadre.KokkosParser(sys.argv)\n"
                "try:\n"
                "    pycompadre.KokkosParser(['program', '--kokkos-disable-warnings'])\n"
                "    raise AssertionError('a list other than sys.argv or sys.orig_argv must raise')\n"
                "except RuntimeError:\n"
                "    pass\n"
                "del kp, kp_argv\n") == 0, "embedded, pycompadre.KokkosParser joins the program's Kokkos");
        Py_Finalize();
        check(!captured.contains(ignored_arguments_warning),
                "embedded, pycompadre.KokkosParser() joins the program's Kokkos without a warning");
    }
    check(Kokkos::is_initialized(), "pycompadre.KokkosParser does not finalize the program's Kokkos");
    Kokkos::finalize();
}
#endif

} // anonymous namespace

int main(int argc, char* argv[]) {
    const std::string mode = (argc > 1) ? argv[1] : "";
    if (mode == "args") {
        testArgs(argc, argv);
    } else if (mode == "no_args") {
        testNoArgs();
    } else if (mode == "external") {
        testExternal(argc, argv);
#ifdef COMPADRE_TEST_EMBEDDED_PYTHON
    } else if (mode == "embedded") {
        testEmbedded();
    } else if (mode == "embedded_external") {
        testEmbeddedExternal(argc, argv);
#endif
    } else {
        std::cout << "FAILED: unknown mode '" << mode << "'" << std::endl;
        failures++;
    }
    if (failures == 0) {
        std::cout << "Passed test" << std::endl;
        return 0;
    }
    std::cout << "Failed test" << std::endl;
    return 1;
}
