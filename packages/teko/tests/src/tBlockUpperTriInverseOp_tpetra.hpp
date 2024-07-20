// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tBlockUpperTriInverseOp_tpetra_hpp__
#define __tBlockUpperTriInverseOp_tpetra_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#include <string>
#include <vector>

#include "Test_Utils.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {
namespace Test {

class tBlockUpperTriInverseOp_tpetra : public UnitTest {
 public:
  virtual ~tBlockUpperTriInverseOp_tpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  bool test_apply(int verbosity, std::ostream& os);
  bool test_alphabeta(int verbosity, std::ostream& os);

 protected:
  ST tolerance_;
  Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<ST> > A_;
  Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<ST> > invA_;

  std::vector<Teuchos::RCP<const Thyra::LinearOpBase<ST> > > invDiag_;
};

}  // namespace Test
}  // end namespace Teko

#endif
