// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tJacobi2x2PreconditionerFactory_tpetra_hpp__
#define __tJacobi2x2PreconditionerFactory_tpetra_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tJacobi2x2PreconditionerFactory_tpetra : public UnitTest {
 public:
  virtual ~tJacobi2x2PreconditionerFactory_tpetra() {}

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
  bool test_initializeFromParameterList(int verbosity, std::ostream& os);

 protected:
  // some simple matrix subblocks
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > A_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > F_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > D_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > G_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > C_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > invF_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > invC_;
  Teuchos::RCP<const Teuchos::Comm<int> > comm;

  ST tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
