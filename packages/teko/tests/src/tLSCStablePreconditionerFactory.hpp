// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tLSCStablePreconditionerFactory_hpp__
#define __tLSCStablePreconditionerFactory_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include "Epetra_SerialComm.h"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tLSCStablePreconditionerFactory : public UnitTest {
 public:
  virtual ~tLSCStablePreconditionerFactory() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  bool test_createPrec(int verbosity, std::ostream& os);
  bool test_initializePrec(int verbosity, std::ostream& os);
  bool test_uninitializePrec(int verbosity, std::ostream& os);
  bool test_isCompatable(int verbosity, std::ostream& os);

  // non-member tests
  bool test_result(int verbosity, std::ostream& os);
  bool test_identity(int verbosity, std::ostream& os);
  bool test_diagonal(int verbosity, std::ostream& os);

 protected:
  // some simple matrix subblocks
  Teuchos::RCP<const Thyra::LinearOpBase<double> > A_;
  Teuchos::RCP<const Thyra::LinearOpBase<double> > F_;
  Teuchos::RCP<const Thyra::LinearOpBase<double> > B_;
  Teuchos::RCP<const Thyra::LinearOpBase<double> > Bt_;
  Teuchos::RCP<const Thyra::LinearOpBase<double> > invF_;
  Teuchos::RCP<const Thyra::LinearOpBase<double> > invBQBt_;
  Teuchos::RCP<const Thyra::LinearOpBase<double> > invMass_;
  Teuchos::RCP<Epetra_SerialComm> comm;

  double tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
