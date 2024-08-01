// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tLU2x2PreconditionerFactory_tpetra_hpp__
#define __tLU2x2PreconditionerFactory_tpetra_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include <string>

#include "Test_Utils.hpp"

#include "Teko_Utilities.hpp"

namespace Teko {
namespace Test {

class tLU2x2PreconditionerFactory_tpetra : public UnitTest {
 public:
  virtual ~tLU2x2PreconditionerFactory_tpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  bool test_createPrec(int verbosity, std::ostream& os);
  bool test_initializePrec(int verbosity, std::ostream& os);
  bool test_uninitializePrec(int verbosity, std::ostream& os);
  bool test_isCompatable(int verbosity, std::ostream& os);

  // non-member tests
  bool test_alphabeta(int verbosity, std::ostream& os);
  bool test_result(int verbosity, std::ostream& os);
  bool test_identity(int verbosity, std::ostream& os);
  bool test_diagonal(int verbosity, std::ostream& os);

 protected:
  // some simple matrix subblocks
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > A_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > F_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > B_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > Bt_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > invF_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > invS_;
  Teuchos::RCP<const Teuchos::Comm<int> > comm;

  double tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
