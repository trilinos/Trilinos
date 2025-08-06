// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Amesos2_Util.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultSerialComm.hpp>

#ifdef HAVE_MPI
#  include <Teuchos_DefaultMpiComm.hpp>
#endif // HAVE_MPI

/// Prints a line of 80 "-"s on out.
void Amesos2::Util::printLine( Teuchos::FancyOStream& out )
{
  out << "----------------------------------------"
      << "----------------------------------------"
      << std::endl;
}


