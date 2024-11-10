// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Util.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include <vector>

namespace { // (anonymous)

  TEUCHOS_UNIT_TEST( Utils, VerbosePrintArray_threshold )
  {
    using Tpetra::Details::verbosePrintArray;

    std::vector<int> x {{3, 5, 7, 9, 11}};
    std::ostringstream os;

    verbosePrintArray(os, x, "x", 10);
    os << ", ";
    verbosePrintArray(os, x, "x2", 3);
    os << ", ";
    verbosePrintArray(os, x, "x3", 5);

    const std::string expected
      ("x: [3, 5, 7, 9, 11], x2: [3, 5, 7, ...], x3: [3, 5, 7, 9, 11]");
    TEST_EQUALITY( os.str(), expected );
  }

  TEUCHOS_UNIT_TEST( Utils, VerbosePrintArray_empty )
  {
    using Tpetra::Details::verbosePrintArray;

    std::vector<int> x;
    std::ostringstream os;

    verbosePrintArray(os, x, "x", 10);
    os << ", ";
    verbosePrintArray(os, x, "x2", 3);
    os << ", ";
    verbosePrintArray(os, x, "x3", 5);

    const std::string expected("x: [], x2: [], x3: []");
    TEST_EQUALITY( os.str(), expected );
  }

  TEUCHOS_UNIT_TEST( Utils, VerbosePrintArray_zero_threshold )
  {
    using Tpetra::Details::verbosePrintArray;

    std::vector<int> x {{3, 5, 7, 9, 11}};
    std::vector<double> y;
    std::ostringstream os;

    verbosePrintArray(os, x, "x", 0);
    os << ", ";
    verbosePrintArray(os, x, "x2", 0);
    os << ", ";
    verbosePrintArray(os, x, "x3", 0);
    os << ", ";
    verbosePrintArray(os, y, "y", 0);

    const std::string expected
      ("x: [...], x2: [...], "
       "x3: [...], y: []");
    TEST_EQUALITY( os.str(), expected );
  }

} // namespace (anonymous)
