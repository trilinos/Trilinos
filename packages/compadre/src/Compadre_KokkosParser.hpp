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
#include <sstream>

namespace Compadre {

/*! \class KokkosParser
    \brief Class handling Kokkos command line arguments and returning parameters.
*/
class KokkosParser {

private:

  // prevent default constructor
  KokkosParser();

  Kokkos::ScopeGuard* ksg;

public:

  // call with command line arguments
  KokkosParser(KokkosInitArguments args, bool print_status = false);

  // call with command line arguments
  KokkosParser(int argc, char* args[], bool print_status = false);

  // call with std::vector of std::string's
  KokkosParser(std::vector<std::string> args, bool print_status = false);

  // call for default arguments
  KokkosParser(bool print_status = false);

  ~KokkosParser() {
      delete ksg;
  }

  // prints Kokkos configuration
  static std::string status();

  // prohibit using the assignment constructor
  KokkosParser& operator=( const KokkosParser& ) = delete;

};

} // Compadre

#endif
