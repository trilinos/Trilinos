// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tLSCHIntegrationTest_tpetra_hpp__
#define __tLSCHIntegrationTest_tpetra_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tLSCHIntegrationTest_tpetra : public UnitTest {
 public:
  virtual ~tLSCHIntegrationTest_tpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  bool test_hScaling(int verbosity, std::ostream& os);

 protected:
  ST tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
