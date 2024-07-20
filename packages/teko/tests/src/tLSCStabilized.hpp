// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tLSCStabilized_hpp__
#define __tLSCStabilized_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include "Epetra_SerialComm.h"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tLSCStabilized : public UnitTest {
 public:
  virtual ~tLSCStabilized() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  // non-member tests
  bool test_diagonal(int verbosity, std::ostream& os);
  bool test_diagonalNotSym(int verbosity, std::ostream& os);
  bool test_strategy(int verbosity, std::ostream& os);

 protected:
  // some simple matrix subblocks
  Teuchos::RCP<Epetra_SerialComm> comm;

  double tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
