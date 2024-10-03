// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tBlockJacobiPreconditionerFactory_tpetra_hpp__
#define __tBlockJacobiPreconditionerFactory_tpetra_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include "Epetra_SerialComm.h"

#include <string>

#include "Test_Utils.hpp"

#include "Teko_Utilities.hpp"

namespace Teko {
namespace Test {

class tBlockJacobiPreconditionerFactory_tpetra : public UnitTest {
 public:
  virtual ~tBlockJacobiPreconditionerFactory_tpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_createPrec(int verbosity, std::ostream& os);
  bool test_initializePrec(int verbosity, std::ostream& os);
  bool test_uninitializePrec(int verbosity, std::ostream& os);
  bool test_isCompatible(int verbosity, std::ostream& os);
  bool test_iterativeSolves(int verbosity, std::ostream& os);

 protected:
  double tolerance_;

  Teuchos::RCP<const Thyra::LinearOpBase<ST> > F_, C_, B_, Bt_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > invF_, invC_;
};

}  // namespace Test
}  // end namespace Teko

#endif
