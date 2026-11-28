// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_KOKKOSPARSER_HPP_
#define _COMPADRE_KOKKOSPARSER_HPP_

#include "Compadre_Config.h"
#include "Compadre_Typedefs.hpp"
#include <memory>
#include <sstream>

namespace Compadre {

/*! \class KokkosParser
    \brief Initializes Kokkos from the command line, or joins the Kokkos already running.

    Kokkos is initialized from the program's command-line arguments whenever they are available, together
    with the KOKKOS_* environment variables (which Kokkos always reads), and never from arguments that were
    not given on the command line.

    In C++:

    | Constructor                                    | Kokkos not initialized                 | Kokkos already initialized |
    |------------------------------------------------|----------------------------------------|----------------------------|
    | KokkosParser(argc, argv), KokkosParser(args)   | initializes with these arguments       | throws std::runtime_error  |
    | KokkosParser(KokkosInitArguments)              | initializes with these settings        | throws std::runtime_error  |
    | KokkosParser()                                 | warns that command-line arguments are  | joins the running Kokkos   |
    |                                                | ignored, and initializes from the      |                            |
    |                                                | KOKKOS_* environment variables only    |                            |

    - A program that initializes Kokkos itself needs no KokkosParser: GMLS and the other classes use the
      running Kokkos. Constructing KokkosParser(argc, argv) after that throws, since Kokkos would be initialized
      twice with arguments in one program. Code without access to the command line uses KokkosParser().
    - The arguments are copied before they are given to Kokkos, so the caller's argc and argv are not modified
      (Kokkos::initialize removes the arguments it recognizes).
    - Kokkos cannot be initialized again once it is finalized, so any KokkosParser constructed after that throws.

    From Python (pycompadre, Python 3.10+), the arguments come from sys.orig_argv, the command line given to the
    Python executable, which cannot change during a run:
    - run by the Python executable (sys.orig_argv is not empty): joins the running Kokkos if it is initialized,
      and otherwise initializes Kokkos with sys.orig_argv;
    - embedded in a C++ program (sys.orig_argv is empty): behaves as KokkosParser() does in C++.
    An argument list passed from Python must be sys.argv or sys.orig_argv (the same list or an equal one);
    any other list throws, since Kokkos is always initialized from sys.orig_argv and the list would be ignored.

    Lifetime: Kokkos started by a KokkosParser is finalized when the last KokkosParser and the last object
    holding runtime() (a GMLS object, for example) are destroyed, in whatever order. Kokkos initialized outside
    Compadre is never finalized by Compadre.
*/
class KokkosParser {

private:

  // shared ownership of the Kokkos runtime started by a KokkosParser (empty if Kokkos was initialized elsewhere)
  std::shared_ptr<void> _runtime;

public:

  //! Initializes Kokkos with these settings (combined with the KOKKOS_* environment variables, as Kokkos does).
  //! Throws if Kokkos is already initialized, or was finalized.
  KokkosParser(KokkosInitArguments args);

  //! Initializes Kokkos with a copy of the command-line arguments (argc and argv are not modified).
  //! Throws if Kokkos is already initialized, or was finalized.
  KokkosParser(int argc, char* args[]);

  //! Initializes Kokkos with the command-line arguments given as strings (argument 0 is the program name).
  //! Throws if Kokkos is already initialized, or was finalized.
  KokkosParser(std::vector<std::string> args);

  //! For code without access to the command line: joins the running Kokkos if it is initialized; otherwise
  //! warns that command-line arguments are ignored and initializes Kokkos from the KOKKOS_* environment
  //! variables only. Throws if Kokkos was finalized.
  KokkosParser();

  //! Shares ownership of the Kokkos runtime started by a KokkosParser, so that an object owning Kokkos
  //! memory can keep Kokkos initialized until it is destroyed: Kokkos is finalized when the last
  //! KokkosParser and the last holder of runtime() are gone, in whatever order. Empty if Kokkos was not
  //! started by a KokkosParser (or is not running); then whoever initialized Kokkos finalizes it.
  static std::shared_ptr<void> runtime();

  // prints Kokkos configuration
  static std::string status();

  // KokkosParser is not copyable. Deleting the copy constructor explicitly (not only the assignment)
  // keeps it from being generated implicitly, e.g. by nanobind for pycompadre, which clang reports as
  // deprecated (-Wdeprecated-copy) when a copy assignment operator is user-declared.
  KokkosParser(const KokkosParser&) = delete;
  KokkosParser& operator=( const KokkosParser& ) = delete;

};

} // Compadre

#endif
