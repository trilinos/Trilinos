// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tTpetraOperatorWrapper_hpp__
#define __tTpetraOperatorWrapper_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"

#include <string>

#include "Test_Utils.hpp"
#include "Teko_ConfigDefs.hpp"

namespace Teko {
namespace Test {

class tTpetraOperatorWrapper : public UnitTest {
 public:
  virtual ~tTpetraOperatorWrapper() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_functionality(int verbosity, std::ostream& os);
};

}  // namespace Test
}  // end namespace Teko

#endif
