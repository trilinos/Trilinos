// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tLSCStabilized_tpetra_hpp__
#define __tLSCStabilized_tpetra_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include <string>

#include "Test_Utils.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {
namespace Test {

class tLSCStabilized_tpetra : public UnitTest {
 public:
  virtual ~tLSCStabilized_tpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  // non-member tests
  bool test_diagonal(int verbosity, std::ostream& os);
  bool test_diagonalNotSym(int verbosity, std::ostream& os);
  bool test_strategy(int verbosity, std::ostream& os);

 protected:
  // some simple matrix subblocks
  Teuchos::RCP<const Teuchos::Comm<int> > comm;

  ST tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
