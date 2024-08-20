// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tNeumannSeries_tpetra_hpp__
#define __tNeumannSeries_tpetra_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Teko includes
#include "Teko_Utilities.hpp"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tNeumannSeries_tpetra : public UnitTest {
 public:
  virtual ~tNeumannSeries_tpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  bool test_simpleOp(int verbosity, std::ostream& os);
  bool test_scaledOp(int verbosity, std::ostream& os);

 protected:
  ST tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
