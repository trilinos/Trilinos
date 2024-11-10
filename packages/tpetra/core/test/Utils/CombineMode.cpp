// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CombineMode.hpp"

namespace { // (anonymous)

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( TpetraUtils, CombineModeToString )
  {
    {
      const std::string cms = Tpetra::combineModeToString (Tpetra::ADD);
      TEST_EQUALITY( cms, "ADD" );
    }
    {
      const std::string cms = Tpetra::combineModeToString (Tpetra::INSERT);
      TEST_EQUALITY( cms, "INSERT" );
    }
    {
      const std::string cms = Tpetra::combineModeToString (Tpetra::REPLACE);
      TEST_EQUALITY( cms, "REPLACE" );
    }
    {
      const std::string cms = Tpetra::combineModeToString (Tpetra::ABSMAX);
      TEST_EQUALITY( cms, "ABSMAX" );
    }
    {
      const std::string cms = Tpetra::combineModeToString (Tpetra::ZERO);
      TEST_EQUALITY( cms, "ZERO" );
    }
  }

} // namespace (anonymous)


